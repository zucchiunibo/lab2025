#include <TCanvas.h>
#include <TRandom.h>
#include <iostream>
#include "TGraph.h"
#include <vector>
#include <TF1.h>
#include <TFile.h>

class macro{
    double k_{5.2};
    double phi_{1.8};
    double b_{0.2};
public:

//macro(double k, double phi, double b) : k_{k} , phi_{phi} , b_{b} {}
 
TF1* cos_function() {
    TF1 *cos = new TF1("Funzione coseno", "(cos([0]*x + [1]))^2 + [2]", -0.5, 5.);
    cos->SetParameters(k_, phi_, b_);
    return cos;
}

void draw() {
    TCanvas* c = new TCanvas("c", "Funzione coseno", 800, 600);
    cos_function()->Draw();
    c->SaveAs("test.png");
}



};