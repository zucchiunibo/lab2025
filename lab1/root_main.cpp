#include "root_macro.C"

int main(int arcg, char** argv) {
    gRandom->SetSeed(0);
    macro test;
    test.draw(); // ciao
    test.accordo(10000, 100);
    test.rigenerazione_incertezze();
    test.binSmearing();
    return 0;
}