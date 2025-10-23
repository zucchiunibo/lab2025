#include "TGraph.h"
#include "TGraphErrors.h"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TRandom.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <TFitResult.h>

class macro {
  double k_;
  double phi_;
  double b_;

 public:
  macro(double k = 5.2, double phi = 1.8, double b = 0.2)
      : k_{k}
      , phi_{phi}
      , b_{b} {
  }

  TF1* cos_function(double norm = 1.) {
    TF1* cos = new TF1("Funzione coseno", "[3]*((cos([0]*x + [1]))^2 + [2])", 0., 0.6);
    cos->SetParameters(k_, phi_, b_, norm);
    return cos;
  }

  TGraph* random_generation_graph(int n) {
    std::vector<double> vx, vy;
    for (int i{0}; i < n; ++i) {
      double x           = gRandom->Uniform(0., 0.6);
      double upper_bound = cos_function()->Eval(x);
      double y           = gRandom->Uniform(0., 1.2);
      if (y <= upper_bound) {
        vx.push_back(x);
        vy.push_back(y);
      }
    }
    TGraph* graph = new TGraph(vx.size(), &vx[0], &vy[0]);
    return graph;
  }

  TH1F* random_generation_hist(int n, int b) {
    std::vector<double> vx;
    int entries = 0;
    for (int i = 0; i < n; ++i) {
      double x           = gRandom->Uniform(0., 0.6);
      double y           = gRandom->Uniform(0., 1.2);
      double upper_bound = cos_function()->Eval(x);
      if (y <= upper_bound) {
        vx.push_back(x);
        entries += 1;
      }
    }
    static int histCount = 0;
    TH1F* hist           = new TH1F(Form("hist_%d", histCount++), "Istogramma Occorrenze", b, 0, 0.6);
    for (double val : vx)
      hist->Fill(val);

    return hist;
  }

  void accordo(int n = 10000, int b = 50) { // normalizza HIST e TF1
    TH1F* hist1 = (TH1F*)(random_generation_hist(n, b)->Clone("hist1"));
    hist1->Scale(1. / hist1->Integral()), "width";

    TF1* cos1     = (TF1*)(cos_function()->Clone("cos1"));
    double cosInt = cos1->Integral(0., 0.6);

    TF1* cosScaled = (TF1*)(cos_function(1 / cosInt)->Clone("cosScaled"));
    std::cout << "Hist integral: " << hist1->Integral() << "\n";
    std::cout << "Cos integral: " << cosScaled->Integral(0., 0.6) << "\n";

    std::vector<double> diff, sigma;
    auto binWidth = 0.6 / b;
    for (int i = 0; i < b; ++i) {
      double xlow        = hist1->GetBinLowEdge(i + 1);
      double xup         = hist1->GetBinLowEdge(i + 2);
      double cosIntegral = cosScaled->Integral(xlow, xup);
      diff.push_back(cosIntegral - hist1->GetBinContent(i + 1));
      sigma.push_back(cosIntegral);
    }

    double chiSquared;

    for (int i = 0; i < b; ++i) {
      chiSquared += std::pow(diff[i] / std::sqrt(sigma[i]), 2);
    }

    std::cout << "Chi quadro: " << chiSquared << "\n";

    TCanvas* c4 = new TCanvas("c4", "Hist scalato", 800, 600);
    hist1->Draw("HIST");
    c4->SaveAs("istrogramma_norm.png");
    TCanvas* c5 = new TCanvas("c5", "Coseno scalato", 800, 600);
    cosScaled->Draw();
    c5->SaveAs("cos_scalato.png");
  }

  struct bin_mean_sigma{
    std::vector<double> media;
    std::vector<double> sigma;
  };

  bin_mean_sigma rigenerazione_incertezze(int nGenerazioni = 100, int nEventi = 10000, int nBin = 50) {
    std::vector<TH1F*> histSet;

    for (int i = 0; i < nGenerazioni; ++i) {
      histSet.push_back(random_generation_hist(nEventi, nBin));
    }

    std::vector<double> media(nBin, 0.0);
    std::vector<double> sigma(nBin, 0.0);

    for (int i = 0; i < nBin; ++i) {
      for (auto& h : histSet)
        media[i] += h->GetBinContent(i + 1);
      media[i] /= nGenerazioni;
    }

    for (int i = 0; i < nBin; ++i) {
      for (auto& h : histSet)
        sigma[i] += pow(h->GetBinContent(i + 1) - media[i], 2);
      sigma[i] = sqrt(sigma[i] / (nGenerazioni - 1));
    }

    TGraphErrors* gSigma = new TGraphErrors(nBin);
    for (int i = 0; i < nBin; ++i) {
      gSigma->SetPoint(i, i, media[i]);
      gSigma->SetPointError(i, 0., sigma[i]); // Argomenti: pos in lista, x, y
    }
    
    TCanvas* c = new TCanvas("c_sigma", "Incertezze per bin", 800, 600);
    gSigma->SetTitle("Fluttuazioni bin; Bin; Deviazione standard");
    gSigma->SetMarkerStyle(20);
    gSigma->Draw("AP");
    c->SaveAs("sigma.png");

    return {media, sigma};
  }

  void binSmearing(int b = 50, int gauss=30) {
    bin_mean_sigma bin = rigenerazione_incertezze();
    std::vector<double> g_media(b, 0.0);
    std::vector<double> g_media_2(b, 0.0);
    std::vector<double> g_sigma(b, 0.0);
 
    for (int i{0}; i < b; ++i) {
      for (int j{0}; j < gauss; ++j) {
        double bincontent = gRandom->Gaus(bin.media[i], bin.sigma[i]);
        g_media[i]+=bincontent;
        g_media_2[i]+= bincontent * bincontent;
      }
      g_media[i] /= gauss;
      g_media_2[i] /= gauss;
      g_sigma[i] = std::sqrt(g_media_2[i] - std::pow(g_media[i], 2));
      std::cout << "sigma_dist: " << std::abs(bin.sigma[i] - g_sigma[i]) << '\n';
    }
  }

  void fit() {
    TF1* cos = new TF1("Funzione coseno", "[3]*((cos([0]*x + [1]))^2 + [2])", 0., 0.6);
    cos->SetParameters(k_, phi_, b_);
    TH1F *hist = random_generation_hist(10000, 50);
    auto freepar = hist->Fit(cos, "RSQ");
    TF1* cos1 = new TF1("Funzione coseno", "[3]*((cos([0]*x + [1]))^2 + [2])", 0., 0.6);
    cos1->FixParameter(0, k_);
    cos1->FixParameter(1, phi_);
    cos1->FixParameter(2, b_);
    auto fixpar = hist->Fit(cos1, "RSQ");

    int statusfree = freepar->Status();
    int statusfix = fixpar->Status();
    std::ofstream ofs("Fit.txt");
      if(!ofs.is_open()){
        std::cout << "Error" << '\n';
      }
    ofs << "# Fit results\n";
    ofs << "# Status: " << statusfree << "\n";
    ofs << "# Function: " << cos->GetName() << "\n";

    double chi2 = cos->GetChisquare();
    int ndf     = cos->GetNDF();
    ofs << "# Chi2 NDF\n" << chi2 << " " << ndf << "\n";

    int npar = cos->GetNpar();
    ofs << "# NParameters: " << npar << "\n";
    ofs << "Index,Name,Value,Error\n";
    for (int i = 0; i < npar; ++i) {
      const char* pname = cos->GetParName(i) ? cos->GetParName(i) : "";
      ofs << i << "," << pname << "," << cos->GetParameter(i) << "," << cos->GetParError(i) << "\n";
    }
    ofs.close();
  }

  void draw() {
    TCanvas* c1 = new TCanvas("c1", "Funzione coseno", 800, 600);
    cos_function()->Draw();
    c1->SaveAs("grafico.png");

    TCanvas* c2 = new TCanvas("c2", "Estrazione punti", 800, 600);
    random_generation_graph(10000)->Draw("AP");
    c2->SaveAs("punti.png");

    TCanvas* c3 = new TCanvas("c3", "Istogramma", 800, 600);
    random_generation_hist(10000, 50)->Draw();
    // cos_function()->Draw();
    c3->SaveAs("istogramma.png");
  }
};

