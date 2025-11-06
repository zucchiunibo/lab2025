#include "root_macro.C"
#include "double_slit.C"

int main(int argc, char** argv) {
    gRandom->SetSeed(0);
    macro test;
    // test.draw();
    // test.accordo(10000, 100);
    // test.rigenerazione_incertezze();
    // test.binSmearing();
    // test.gaussian_parameters();
    // test.fit();
    // test.gaussian_smearing();
    // test.chi2_different_method();

    double_slit DB_Slit;
    DB_Slit.fit();
    
    return 0;
}