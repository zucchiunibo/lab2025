#include "Fit_Diffrazione.h"
int main() {

  FitParams p1 = {
      "diff1(step100).txt", 0., 1.4, 677.e-9, 5.8e-5, 0.0376, 0.935};
  FitParams p2 = {
      "diff2(step100).txt", 0., 1.42, 695.e-9, 5.8e-5, 0.0376, 0.935};
  FitParams p3 = {
      "diff3(step100)fenditura1", 0., 3.7012, 702.e-9, 1.01e-4, 0.0377, 0.935};

  FitParams p4 = {"diff2.txt", 0., 1.33, 683.e-9, 5.8e-5, 0.037, 0.935};

  FitParams p5 = {"diff3.txt", 0., 1.39, 714.e-9, 5.8e-5, 0.0365, 0.935};

  FitParams p6 = {"diffrazione_1mm", 0., 2.76, 683.e-9, 5.8e-5, 0.0375, 0.935};

  mydata(p6);
  myfunc(p6);
  myfit(p6);

  std::cout << "Fine!" << '\n';
}
