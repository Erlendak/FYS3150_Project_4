#ifndef UNITTESTS_H
#define UNITTESTS_H
#include <armadillo>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
using namespace std;

void test_temp_1(){
    int n_spins = 4; //La stå for senere bruk!
    double temp= 1;
    double j = 1;
    double kb = 1;//1.38064852*pow(10,-23)
    double beta = 1/(kb*temp);
    double Cv = ( (4*64*j*j)/(kb*temp*temp) )*( ( cosh(8*j*beta)/(12+4*sinh(8*j*beta) ) ) - ( cosh(2*8*j*beta)/pow((12+4*sinh(8*j*beta) ) ,2) ) ) ;
    cout <<"Analytisk beregnet varme ved temperatur 1. ; ";
    cout<<Cv <<endl;
    cout<<"\n\n" <<endl;
}


#endif // UNITTESTS_H
