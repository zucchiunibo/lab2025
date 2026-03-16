#include "Fit_Diffrazione.h"
int main(){
    mydata("diff1(step100)");
    myfunc(0., 1.41, 678.e9, 5.8e-5, 0.0376, 0.935);
    myfit("diff1(step100)", 0., 1.41, 678.e9, 5.8e-5, 0.0376, 0.935);
    std::cout << "Fine!" << '\n';
}
