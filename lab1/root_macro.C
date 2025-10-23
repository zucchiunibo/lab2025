#include "TGraph.h"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TRandom.h>
#include <iostream>
#include <vector>

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
    // for (float i{0.f}; i < 0.6f; i += 0.6 / n) {
    //   double upper_bound = cos_function()->Eval(i);
    //   double y           = gRandom->Uniform(0., 1.2);
    //   if (y <= upper_bound) {
    //     vx.push_back(i);
    //     vy.push_back(y);
    //   }
    // }
    // OPPURE
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

  // TH1F* random_generation_hist(int n, int b) {
  //   std::vector<double> vx, vy;
  //   for (float i{0.f}; i < 0.6f; i += 0.6 / n) {
  //     double upper_bound = cos_function()->Eval(i);
  //     double y           = gRandom->Uniform(0., 1.2);
  //     if (y <= upper_bound) {
  //       vx.push_back(i);
  //     }
  //   }
  //   TH1F* hist = new TH1F("hist", "Istogramma Occorrenze", b, 0, 0.6);
  //   for (unsigned int i{0}; i < vx.size(); ++i) {
  //     hist->Fill(vx[i], 1);
  //   }
  //   return hist;
  // }

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
    std::cout << "entries: " << entries << '\n';
    static int histCount = 0;
    TH1F* hist = new TH1F(Form("hist_%d", histCount++), "Istogramma Occorrenze", b, 0, 0.6);
    for (double val : vx)
      hist->Fill(val);

    return hist;
  }

  void accordo(int n = 1000, int b = 50) { // normalizza HIST e TF1
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
      diff.push_back(cosIntegral - hist1->GetBinContent(i));
      sigma.push_back(cosIntegral);
      // std::cout << "Differnenza bin " << i << " " << cosIntegral - hist1->GetBinContent(i) << "\n";
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

  void rigenerazione_incertezze(int nGenerazioni = 100, int nEventi = 10000, int nBin = 50) {
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

    TGraph* gSigma = new TGraph(nBin);
    for (int i = 0; i < nBin; ++i)
      gSigma->SetPoint(i, i, sigma[i]); //Argomenti: pos in lista, x, y

    TCanvas* c = new TCanvas("c_sigma", "Incertezze per bin", 800, 600);
    gSigma->SetTitle("Fluttuazioni bin; Bin; Deviazione standard");
    gSigma->SetMarkerStyle(20);
    gSigma->Draw("AP");
    c->SaveAs("sigma.png");
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