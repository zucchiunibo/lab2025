#include "TGraph.h"
#include "TGraphErrors.h"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TRandom.h>
#include <fstream>
#include <iostream>
#include <vector>

Double_t Diffrazione(double *x, double *par);

void myfunc(double bkg, double I0, double lambda, double d, double x0, double L);

void mydata(TString fname);

void myfit(TString fname = "diff1(step100)", double bkg = 0., double I0 = 1.41,
           double lambda = 678.E-9, double d = 5.8E-5, double x0 = 0.0376,
           double L = 0.935);