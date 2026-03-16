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

void myfunc(double bkg, double I0, double lambda,
            double d, double x0, double L);

void mydata(TString fname);

void myfit(TString fname, double bkg, double I0, double lambda, double d,
           double x0, double L);