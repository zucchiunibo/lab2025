#include "TGraph.h"
#include "TGraphErrors.h"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1F.h>
#include <TRandom.h>
#include <fstream>
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

  TH1F* random_generation_hist(int n, int b, TF1* f = nullptr) {
    if (!f)
      f = cos_function();
    std::vector<double> vx;
    int entries = 0;
    for (int i = 0; i < n; ++i) {
      double x           = gRandom->Uniform(0., 0.6);
      double y           = gRandom->Uniform(0., 1.2);
      double upper_bound = f->Eval(x);
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
    cosScaled->SetTitle("Funzione normalizzata");
    cosScaled->Draw();
    c5->SaveAs("cos_scalato.png");
  }

  struct bin_mean_sigma {
    std::vector<double> media;
    std::vector<double> sigma;
  };

  bin_mean_sigma rigenerazione_incertezze(int nGenerazioni = 100, int nEventi = 10000, int nBin = 50, TF1* f = nullptr) {
    std::vector<TH1F*> histSet;

    for (int i = 0; i < nGenerazioni; ++i) {
      histSet.push_back(random_generation_hist(nEventi, nBin, f));
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

  void binSmearing(int b = 50, int gauss = 30) {
    bin_mean_sigma bin = rigenerazione_incertezze();
    std::vector<double> g_media(b, 0.0);
    std::vector<double> g_media_2(b, 0.0);
    std::vector<double> g_sigma(b, 0.0);

    for (int i{0}; i < b; ++i) {
      for (int j{0}; j < gauss; ++j) {
        double bincontent = gRandom->Gaus(bin.media[i], bin.sigma[i]);
        g_media[i] += bincontent;
        g_media_2[i] += bincontent * bincontent;
      }
      g_media[i] /= gauss;
      g_media_2[i] /= gauss;
      g_sigma[i] = std::sqrt(g_media_2[i] - std::pow(g_media[i], 2));
      std::cout << "sigma_dist: " << std::abs(bin.sigma[i] - g_sigma[i]) << '\n';
    }
  }

  void gaussian_parameters(int n_generazioni = 100, int n_eventi = 10000, int bins = 50) {
    std::vector<TH1F*> hist_list;
    // TF1* cos_g = new TF1("Funzione coseno", "[3]*((cos([0]*x + [1]))^2 + [2])", 0., 0.6);
    // double cosInt = cos_g->Integral(0., 0.6);
    // TF1* cosScaled_g = (TF1*)(cos_g->Clone("cosScaled_g"));
    // cosScaled_g->SetParameter(3, 1 / cosInt);

    for (int i{0}; i < n_generazioni; ++i) {
      double k   = gRandom->Gaus(k_, k_ * 0.01);
      double phi = gRandom->Gaus(phi_, phi_ * 0.05);
      double b   = gRandom->Gaus(b_, b_ * 0.01);

      TF1* cos_g = new TF1("Funzione coseno", "[3]*((cos([0]*x + [1]))^2 + [2])", 0., 0.6);
      cos_g->SetParameters(k, phi, b);
      double cosInt    = cos_g->Integral(0., 0.6);
      TF1* cosScaled_g = (TF1*)(cos_g->Clone("cosScaled_g"));
      cosScaled_g->SetParameter(3, 1 / cosInt);
      hist_list.push_back(random_generation_hist(n_eventi, bins, cosScaled_g));
    }

    std::vector<double> media(bins, 0.0);
    std::vector<double> sigma(bins, 0.0);

    for (int i = 0; i < bins; ++i) {
      for (auto& h : hist_list)
        media[i] += h->GetBinContent(i + 1);
      media[i] /= n_generazioni;
    }

    for (int i = 0; i < bins; ++i) {
      for (auto& h : hist_list)
        sigma[i] += pow(h->GetBinContent(i + 1) - media[i], 2);
      sigma[i] = sqrt(sigma[i] / (n_generazioni - 1));
    }

    TGraphErrors* gSigma = new TGraphErrors(bins);
    for (int i = 0; i < bins; ++i) {
      gSigma->SetPoint(i, i, media[i]);
      gSigma->SetPointError(i, 0., sigma[i]); // Argomenti: pos in lista, x, y
    }

    TCanvas* c6 = new TCanvas("c_sigma_g", "Incertezze per bin con parametri aleatori", 800, 600);
    gSigma->SetTitle("Fluttuazioni bin; Bin; Deviazione standard");
    gSigma->SetMarkerStyle(20);
    gSigma->Draw("AP");
    c6->SaveAs("sigma_g.png");
  }

<<<<<<< HEAD:lab1/root_macro.C
  void fit() {
    TF1* cos = new TF1("Funzione coseno", "[3]*((cos([0]*x + [1]))^2 + [2])", 0., 0.6);
    cos->SetParameters(k_, phi_, b_);
    TH1F* hist = random_generation_hist(10000, 50);
=======
  void gaussian_smearing(int n_generazioni = 100, int n_eventi = 10000, int bins = 50, int gauss = 30) {
    double k   = gRandom->Gaus(k_, k_ * 0.01);
    double phi = gRandom->Gaus(phi_, phi_ * 0.05);
    double b   = gRandom->Gaus(b_, b_ * 0.01);
>>>>>>> 91172c6f64e5df3838ba19927d6b47f2115c83c9:modulo1/root_macro.C

    TF1* cos_g = new TF1("Funzione coseno", "[3]*((cos([0]*x + [1]))^2 + [2])", 0., 0.6);
    cos_g->SetParameters(k, phi, b);
    double cosInt    = cos_g->Integral(0., 0.6);
    TF1* cosScaled_g = (TF1*)(cos_g->Clone("cosScaled_g"));
    cosScaled_g->SetParameter(3, 1 / cosInt);
    auto bin = rigenerazione_incertezze(n_generazioni, n_eventi, bins, cosScaled_g);
    std::vector<double> g_media(bins, 0.0);
    std::vector<double> g_media_2(bins, 0.0);
    std::vector<double> g_sigma(bins, 0.0);

<<<<<<< HEAD:lab1/root_macro.C
    // Fit con parametri fissati
    TF1* cos1 = new TF1("Funzione coseno", "[3]*((cos([0]*x + [1]))^2 + [2])", 0., 0.6);
    cos1->FixParameter(0, k_);
    cos1->FixParameter(1, phi_);
    cos1->FixParameter(2, b_);
    auto fixpar = hist->Fit(cos1, "RSQ");

    int statusfree = freepar->Status();
    int statusfix  = fixpar->Status();

    std::cout << "residuo: " << residuo(cos1, hist) << '\n';


    std::ofstream ofs("Fit.md");
    if (!ofs.is_open()) {
      std::cout << "Errore apertura file\n";
      return;
    }

    ofs << "# Risultati del Fit\n\n";

    ofs << "## Fit con parametri liberi\n\n";
    ofs << "- **Status:** " << statusfree << "\n";
    ofs << "- **Funzione:** `" << cos->GetName() << "`\n";

    double chi2_free = cos->GetChisquare();
    int ndf_free     = cos->GetNDF();
    ofs << "- **Chi² / NDF:** " << chi2_free << " / " << ndf_free << "\n\n";

    int npar_free = cos->GetNpar();
    ofs << "| Index | Name | Value | Error |\n";
    ofs << "|:------:|:------|:------:|:------:|\n";
    for (int i = 0; i < npar_free; ++i) {
      const char* pname = cos->GetParName(i) ? cos->GetParName(i) : "";
      ofs << "| " << i << " | " << pname << " | " << cos->GetParameter(i) << " | " << cos->GetParError(i) << " |\n";
    }

    ofs << "\n---\n\n";
    ofs << "## Fit con parametri fissati\n\n";
    ofs << "- **Status:** " << statusfix << "\n";
    ofs << "- **Funzione:** `" << cos1->GetName() << "`\n";

    double chi2_fix = cos1->GetChisquare();
    int ndf_fix     = cos1->GetNDF();
    ofs << "- **Chi² / NDF:** " << chi2_fix << " / " << ndf_fix << "\n\n";

    int npar_fix = cos1->GetNpar();
    ofs << "| Index | Name | Value | Error |\n";
    ofs << "|:------:|:------|:------:|:------:|\n";
    for (int i = 0; i < npar_fix; ++i) {
      const char* pname = cos1->GetParName(i) ? cos1->GetParName(i) : "";
      ofs << "| " << i << " | " << pname << " | " << cos1->GetParameter(i) << " | " << cos1->GetParError(i) << " |\n";
    }

    ofs.close();
  }

  double residuo(TF1* cos, TH1F* hist, int b = 50) {
    double diff;
    for (int i = 0; i < b; ++i) {
      double xlow        = hist->GetBinLowEdge(i + 1);
      double xup         = hist->GetBinLowEdge(i + 2);
      double cosIntegral = cos->Integral(xlow, xup);
      double mean        = cosIntegral / (xup - xlow);

      diff += std::pow(mean - hist->GetBinContent(i + 1), 2);
    }
    return std::sqrt(diff);
  }

  void draw() {
    TCanvas* c1 = new TCanvas("c1", "Funzione coseno", 800, 600);
    cos_function()->SetTitle("Funzione");
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
=======
    for (int i{0}; i < bins; ++i) {
      for (int j{0}; j < gauss; ++j) {
        double bincontent = gRandom->Gaus(bin.media[i], bin.sigma[i]);
        g_media[i] += bincontent;
        g_media_2[i] += bincontent * bincontent;
      }
      g_media[i] /= gauss;
      g_media_2[i] /= gauss;
      g_sigma[i] = std::sqrt(g_media_2[i] - std::pow(g_media[i], 2));
    }

  TH1F* hist_smear = new TH1F("hist_smear", "Istogramma da bin smearing; x (0-0.6); Conteggio medio", bins, 0., 0.6);
  for (int i = 0; i < bins; ++i) {
    hist_smear->SetBinContent(i + 1, g_media[i]);
    hist_smear->SetBinError(i + 1,g_sigma[i]);
  }

  TCanvas* c7 = new TCanvas("smearing_sigma_g", "Istogramma bin smearing", 800, 600);
  hist_smear->SetMarkerStyle(20);
  hist_smear->Draw("HIST E"); // mostra istogramma con barre di errore
  c7->SaveAs("hist_smearing.png");
}

  void fit() {
  TF1* cos = new TF1("Funzione coseno", "[3]*((cos([0]*x + [1]))^2 + [2])", 0., 0.6);
  cos->SetParameters(k_, phi_, b_);
  TH1F* hist = random_generation_hist(10000, 50);

  // Fit con parametri liberi
  auto freepar = hist->Fit(cos, "RSQ");

  // Fit con parametri fissati
  TF1* cos1 = new TF1("Funzione coseno", "[3]*((cos([0]*x + [1]))^2 + [2])", 0., 0.6);
  cos1->FixParameter(0, k_);
  cos1->FixParameter(1, phi_);
  cos1->FixParameter(2, b_);
  auto fixpar = hist->Fit(cos1, "RSQ");

  int statusfree = freepar->Status();
  int statusfix  = fixpar->Status();

  std::ofstream ofs("Fit.md");
  if (!ofs.is_open()) {
    std::cout << "Errore apertura file\n";
    return;
  }

  ofs << "# Risultati del Fit\n\n";

  ofs << "## Fit con parametri liberi\n\n";
  ofs << "- **Status:** " << statusfree << "\n";
  ofs << "- **Funzione:** `" << cos->GetName() << "`\n";

  double chi2_free = cos->GetChisquare();
  int ndf_free     = cos->GetNDF();
  ofs << "- **Chi² / NDF:** " << chi2_free << " / " << ndf_free << "\n\n";

  int npar_free = cos->GetNpar();
  ofs << "| Index | Name | Value | Error |\n";
  ofs << "|:------:|:------|:------:|:------:|\n";
  for (int i = 0; i < npar_free; ++i) {
    const char* pname = cos->GetParName(i) ? cos->GetParName(i) : "";
    ofs << "| " << i << " | " << pname << " | " << cos->GetParameter(i) << " | " << cos->GetParError(i) << " |\n";
  }

  ofs << "\n---\n\n";
  ofs << "## Fit con parametri fissati\n\n";
  ofs << "- **Status:** " << statusfix << "\n";
  ofs << "- **Funzione:** `" << cos1->GetName() << "`\n";

  double chi2_fix = cos1->GetChisquare();
  int ndf_fix     = cos1->GetNDF();
  ofs << "- **Chi² / NDF:** " << chi2_fix << " / " << ndf_fix << "\n\n";

  int npar_fix = cos1->GetNpar();
  ofs << "| Index | Name | Value | Error |\n";
  ofs << "|:------:|:------|:------:|:------:|\n";
  for (int i = 0; i < npar_fix; ++i) {
    const char* pname = cos1->GetParName(i) ? cos1->GetParName(i) : "";
    ofs << "| " << i << " | " << pname << " | " << cos1->GetParameter(i) << " | " << cos1->GetParError(i) << " |\n";
  }

  ofs.close();
}

void draw() {
  TCanvas* c1 = new TCanvas("c1", "Funzione coseno", 800, 600);
  cos_function()->SetTitle("Funzione");
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
}
;
>>>>>>> 91172c6f64e5df3838ba19927d6b47f2115c83c9:modulo1/root_macro.C
