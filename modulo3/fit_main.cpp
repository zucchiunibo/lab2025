#include "Fit_Diffrazione.h"
int main() {
  std::tuple<TString, double, double, double, double, double, double> parametri{
      "diff1(step100)", 0., 1.41, 678.E-9, 5.8E-5, 0.0376, 0.935};

  mydata("diff1(step100)");
  myfunc(0., 1.41, 678.e-9, 5.8e-5, 0.0376, 0.935);
  myfit("diff1(step100)", 0., 1.41, 678.e-9, 5.8e-5, 0.0376, 0.935);
  std::cout << "Fine!" << '\n';
}
