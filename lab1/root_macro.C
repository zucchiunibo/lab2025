#include <TCanvas.h>
#include <TRandom.h>
#include <iostream>
#include "TGraph.h"
#include <vector>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>


class macro{
    double k_;
    double phi_;
    double b_;
public:

macro(double k = 5.2, double phi = 1.8, double b = 0.2) : k_{k} , phi_{phi} , b_{b} {}
 
TF1* cos_function() {
    TF1 *cos = new TF1("Funzione coseno", "(cos([0]*x + [1]))^2 + [2]", 0., 0.6);
    cos->SetParameters(k_, phi_, b_);
    return cos;
}

TGraph* random_generation_graph(int n) {
    std::vector<double> vx , vy;
    for (float i{0.f}; i<0.6f; i += 0.6/n) {
        double upper_bound = cos_function()->Eval(i);
        double y = gRandom->Uniform(0., 1.2);
        if (y <= upper_bound){
            vx.push_back(i);
            vy.push_back(y);
        } 
    }
    TGraph* graph = new TGraph(vx.size(), &vx[0], &vy[0]);
    return graph;
}

TH1F* random_generation_hist(int n, int b) {
    std::vector<double> vx , vy;
    for (float i{0.f}; i<0.6f; i += 0.6/n) {
        double upper_bound = cos_function()->Eval(i);
        double y = gRandom->Uniform(0., 1.2);
        if (y <= upper_bound){
            vx.push_back(i);
        } 
    }
    TH1F* hist = new TH1F("hist", "Istogramma Occorrenze", b, 0, 0.6);
    for(unsigned int i{0}; i < vx.size(); ++i) {
        hist->Fill(vx[i], 1);
    }
    return hist;
}

TH1F* accordo(int n, int b) {
    auto hist_a = random_generation_hist(n, b);
    TH1F* hist1 = (TH1F*)(hist_a->Clone("hist1"));
    hist1->Scale(1./hist1->Integral()), "width";
    return hist1;
}

void draw() {
    TCanvas* c1 = new TCanvas("c1", "Funzione coseno", 800, 600);
    cos_function()->Draw();
    c1->SaveAs("grafico.png");

    TCanvas* c2 = new TCanvas("c2", "Estrazione punti", 800, 600);
    random_generation_graph(10000)->Draw();
    c2->SaveAs("punti.png");

    TCanvas* c3 = new TCanvas("c3", "Istogramma", 800, 600);
    random_generation_hist(10000, 50)->Draw();
    // cos_function()->Draw();
    c3->SaveAs("istogramma.png");

    TCanvas* c4 = new TCanvas("c4", "Istogramma norm", 800, 600);
    accordo(10000, 50)->Draw("HIST");
    c4->SaveAs("istogramma_norm.png");
}
};