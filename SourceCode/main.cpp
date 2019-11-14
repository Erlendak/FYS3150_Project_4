#include <iostream>
#include <fstream>
#include <iomanip>
//#include "lib.h"
#include <string>
#include <unittests.h>
#include <monte_carlo.h>
#include <armadillo>
  //ofile.close();  // close output file

using namespace  std;
using namespace  arma;

int main()
{
    //test_temp_1();

    double temp= 1;
    int n_spins = 20;
    int mcs;
    vec average;
    int _E = 60;
    double _x = 1.3;
    int _start = 7;
    ofstream afile;
    string afilename = "20x20_temp_1.dat";
    afile.open(afilename);
    afile << setiosflags(ios::showpoint | ios::uppercase);
    for (int i=_start;i<_E;i++){
            mcs = pow(_x,i);
            cout <<  mcs <<endl;
            afile << setw(15) << setprecision(8) <<  mcs;
            average = isingmodel_cold_start(n_spins, mcs, temp);
            afile << setw(20) << setprecision(8) <<  average(0)/(mcs*n_spins*n_spins);
            afile << setw(20) << setprecision(8) <<  (( average(1)/(mcs) ) - ((average(0)/(mcs))*(average(0)/(mcs))))/(n_spins*n_spins);
            afile << setw(20) << setprecision(8) << average(2)/(mcs*n_spins*n_spins);
            afile << setw(20) << setprecision(8) << average(4)/(mcs*n_spins*n_spins);
            afile << setw(20) << setprecision(8) <<  (( average(3)/(mcs) ) - ((average(4)/(mcs))*(average(4)/(mcs))))/(n_spins*n_spins);
            afile << setw(20) << setprecision(8) << average(5);
            afile << setw(20) << setprecision(8) << average(6)/(average(5))*100;
            afile << setw(20) << setprecision(8) << average(7)/(average(5))*100;
            afile << setw(20) << setprecision(8) << average(8)/(average(5))*100;
            afile << setw(20) << setprecision(8) << average(9)/(average(5))*100;
            afile << setw(20) << setprecision(8) << average(10)/(average(5))*100<<endl;
       }
      afile.close();

      temp = 2.4;

      ofstream bfile;
      string bfilename = "20x20_temp_2_4.dat";
      bfile.open(bfilename);
      bfile << setiosflags(ios::showpoint | ios::uppercase);
      for (int j=_start;j<_E;j++){
         mcs = pow(_x,j);
         cout <<  mcs <<endl;
         bfile << setw(15) << setprecision(8) <<  mcs;
         average = isingmodel_cold_start(n_spins, mcs, temp);
         //cout <<  average[0] <<endl;

         bfile << setw(20) << setprecision(8) <<  average(0)/(mcs*n_spins*n_spins);
         bfile << setw(20) << setprecision(8) <<  (( average(1)/(mcs) ) - ((average(0)/(mcs))*(average(0)/(mcs))))/(n_spins*n_spins);
         bfile << setw(20) << setprecision(8) << average(2)/(mcs*n_spins*n_spins);
         bfile << setw(20) << setprecision(8) << average(4)/(mcs*n_spins*n_spins);
         bfile << setw(20) << setprecision(8) <<  (( average(3)/(mcs) ) - ((average(4)/(mcs))*(average(4)/(mcs))))/(n_spins*n_spins);
         bfile << setw(20) << setprecision(8) << average(5);
         bfile << setw(20) << setprecision(8) << average(6)/(average(5))*100;
         bfile << setw(20) << setprecision(8) << average(7)/(average(5))*100;
         bfile << setw(20) << setprecision(8) << average(8)/(average(5))*100;
         bfile << setw(20) << setprecision(8) << average(9)/(average(5))*100;
         bfile << setw(20) << setprecision(8) << average(10)/(average(5))*100<<endl;
      }
         bfile.close();



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
