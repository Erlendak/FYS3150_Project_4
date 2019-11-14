#ifndef UNITTESTS_H
#define UNITTESTS_H
#include <armadillo>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <monte_carlo.h>
using namespace std;
/*
void test_temp_1(){
    int n_spins = 2; //La stå for senere bruk!
    double temp= 1;
    double j = 1;
    double kb = 1;//1.38064852*pow(10,-23)
    double beta = 1/(kb*temp);
    double z = 4*cosh(beta*j*8)+12;

    //double Cv = ( ( (64*4* cosh(8*j*beta)/z) -
    //(pow(2,10)*j*j*pow((cosh(8*j*beta)), 2) ) )   ) ;//((1)/(kb*temp*temp))
    double Cv = ( (1)/(kb*temp*temp) )*( (64*4* cosh(8*j*beta)/(12+4*sinh(8*j*beta) ) ) - (pow(2,10)* cosh(2*8*j*beta)/pow((12+4*sinh(8*j*beta) ) ,2) ) ) ;
    //Cv = 0;
    cout <<"Analytisk beregnet varme ved temperatur 1. ; ";
    cout<<Cv/4 <<endl;

  cout<<z <<endl;
    cout<<"\n\n" <<endl;
    int mcs = 1000000;
    vec average = temperature_integration(n_spins, mcs, temp);
    cout<<"Antall elektroner som spinner i positiv rettning ; ";
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

*/
#endif // UNITTESTS_H
