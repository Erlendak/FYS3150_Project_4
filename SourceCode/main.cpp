#include <iostream>
#include <fstream>
#include <iomanip>
//#include "lib.h"
#include <string>
#include <unittests.h>
#include <monte_caro.h>
#include <armadillo>
  //ofile.close();  // close output file

using namespace  std;
using namespace  arma;

int main()
{
    test_temp_1();

    double temp= 1;
    int n_spins = 2;
    int mcs;

    int M = 8;
    int _x = 10;
    ofstream afile;
    string afilename = "Temperature_1.dat";
    afile.open(afilename);
    afile << setiosflags(ios::showpoint | ios::uppercase);
    for (int i=2;i<M;i++){
            mcs = pow(_x,i);
            cout <<  mcs <<endl;
            afile << setw(15) << setprecision(8) <<  mcs;
            vec average = temperature_integration(n_spins, mcs, temp);
            afile << setw(20) << setprecision(8) <<  average(0)/(mcs*n_spins*n_spins);
            afile << setw(20) << setprecision(8) <<  (( average(1)/(mcs) ) - ((average(0)/(mcs))*(average(0)/(mcs))))/(n_spins*n_spins);
            afile << setw(20) << setprecision(8) << average(2)/(mcs*n_spins*n_spins);
            afile << setw(20) << setprecision(8) << average(4)/(mcs*n_spins*n_spins);
            afile << setw(20) << setprecision(8) <<  (( average(3)/(mcs) ) - ((average(4)/(mcs))*(average(4)/(mcs))))/(n_spins*n_spins)<<endl;
       }
         afile.close();



    /*
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
    cout<<"\n\n\n"; */
    return 0;
};
