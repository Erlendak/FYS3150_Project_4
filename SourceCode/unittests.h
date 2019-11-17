#ifndef UNITTESTS_H
#define UNITTESTS_H
#include <armadillo>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <monte_carlo.h>

using namespace std;
using namespace arma;

void test_temp_1(){
    int n_spins = 2; //La stå for senere bruk!
    double temp= 1;
    double j = 1;
    double kb = 1;//1.38064852*pow(10,-23)
    double beta = 1/(kb*temp);

    double z = 4*cosh(beta*j*8)+12;

    double E_base = (32*j*sinh(8*j*beta))/z;
    double M_base = (16+8*exp(-8*j*beta))/z;

    double E = (4*64*j*j*cosh(8*j*beta))/z;
    double E2 = (1024*j*j*(cosh(8*j*beta)*cosh(8*j*beta)))/(z*z);

    double Cv = (1/(kb*temp*temp))*(E-E2);

    double Chi = (1/(kb*temp))*(32*(1+exp(-8*j*beta)))/z;

    double Chi2 = (1/(kb*temp))*(32*(1+exp(-8*j*beta))/z)-((16+8*exp(-8*j*beta))*(16+8*exp(-8*j*beta)))/(z*z);

    cout <<"Analytisk beregnet gjennomsnittsenergi ved 1: ";
    cout<<E_base/4 <<endl;

    cout <<"Analytisk beregnet absolutt gjennomsnittsmagnetisasjon 1: ";
    cout<<M_base/4 <<endl;

    cout <<"Analytisk beregnet varme ved temperatur 1: ";
    cout<<Cv/4 <<endl;

    cout <<"Analytisk beregnet susceptibilitet ved temperatur 1:";
    cout<<Chi/4 <<endl;

    cout <<"Analytisk beregnet absolutt susceptibilitet ved temperatur 1. ; ";
    cout<<Chi2/4 <<endl;

    cout<<"\n\n" <<endl;
    int mcs = 1000000;
    vec average =isingmodel_cold_start(n_spins, mcs, temp);
    cout<<"Antall elektroner som spin; ";
    cout<<n_spins*n_spins<<endl;

    cout<<"Antall Monte Carlo integrasjonspoeng ; ";
    cout<<mcs<<endl;

    cout<<"Temperatur i simuleringen ; ";
    cout<<temp<<endl;

    cout<<"<e> ; ";
    cout<<average(0)/(mcs*n_spins*n_spins)<<endl;

    cout<<"Cv ; ";
    cout<<(( average(1)/(mcs) ) - ((average(0)/(mcs))*(average(0)/(mcs))))/(n_spins*n_spins)<<endl;

    cout<<"<m> ; ";
    cout<<average(2)/(mcs*n_spins*n_spins)<<endl;

    cout<<"<|m|> ; ";
    cout<<average(4)/(mcs*n_spins*n_spins)<<endl;

    cout<<"Susceptibilitet ; ";
    cout<<(( average(3)/(mcs) ) - ((average(4)/(mcs))*(average(4)/(mcs))))/(n_spins*n_spins)<<endl;
    cout<<"\n\n\n";


}


#endif // UNITTESTS_H
