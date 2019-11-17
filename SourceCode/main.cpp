#include <iostream>
#include <fstream>
#include <iomanip>
//#include "lib.h"
#include <string>
#include <unittests.h>
#include <monte_carlo.h>
#include <armadillo>
#include <omp.h>
  //ofile.close();  // close output file

using namespace  std;
using namespace  arma;

int main()
{
    int nthreads = omp_get_num_threads();
    printf("Using %d threads outside parallel loop\n",nthreads);
    int mthreads = omp_get_max_threads();
    printf("There are %d threads available at a time\n",mthreads);
    //test_temp_1();

    double temp= 1;
    int n_spins = 20;
    int mcs;

    vec average;

    int _E = 70;
    double _x = 1.2;
    int _start = 7;
    ofstream afile;
    string afilename = "20x20_cold_start_temp_1.dat";
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
      string bfilename = "20x20_cold_start_temp_2_4.dat";
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

      temp= 1;
      ofstream cfile;
      string cfilename = "20x20_random_start_temp_1.dat";
      cfile.open(cfilename);
      cfile << setiosflags(ios::showpoint | ios::uppercase);
      for (int i=_start;i<_E;i++){

          mcs = pow(_x,i);
          cout <<  mcs <<endl;
          cfile << setw(15) << setprecision(8) <<  mcs;
          average = isingmodel_random_start(n_spins, mcs, temp);
          cfile << setw(20) << setprecision(8) <<  average(0)/(mcs*n_spins*n_spins);
          cfile << setw(20) << setprecision(8) <<  (( average(1)/(mcs) ) - ((average(0)/(mcs))*(average(0)/(mcs))))/(n_spins*n_spins);
          cfile << setw(20) << setprecision(8) << average(2)/(mcs*n_spins*n_spins);
          cfile << setw(20) << setprecision(8) << average(4)/(mcs*n_spins*n_spins);
          cfile << setw(20) << setprecision(8) <<  (( average(3)/(mcs) ) - ((average(4)/(mcs))*(average(4)/(mcs))))/(n_spins*n_spins);
          cfile << setw(20) << setprecision(8) << average(5);
          cfile << setw(20) << setprecision(8) << average(6)/(average(5))*100;
          cfile << setw(20) << setprecision(8) << average(7)/(average(5))*100;
          cfile << setw(20) << setprecision(8) << average(8)/(average(5))*100;
          cfile << setw(20) << setprecision(8) << average(9)/(average(5))*100;
          cfile << setw(20) << setprecision(8) << average(10)/(average(5))*100<<endl;
       }
       cfile.close();

       temp = 2.4;
       ofstream dfile;
       string dfilename = "20x20_random_start_temp_2_4.dat";
       dfile.open(dfilename);
       dfile << setiosflags(ios::showpoint | ios::uppercase);
       for (int j=_start;j<_E;j++){

           mcs = pow(_x,j);
           cout <<  mcs <<endl;
           dfile << setw(15) << setprecision(8) <<  mcs;
           average = isingmodel_random_start(n_spins, mcs, temp);
           //cout <<  average[0] <<endl;

           dfile << setw(20) << setprecision(8) <<  average(0)/(mcs*n_spins*n_spins);
           dfile << setw(20) << setprecision(8) <<  (( average(1)/(mcs) ) - ((average(0)/(mcs))*(average(0)/(mcs))))/(n_spins*n_spins);
           dfile << setw(20) << setprecision(8) << average(2)/(mcs*n_spins*n_spins);
           dfile << setw(20) << setprecision(8) << average(4)/(mcs*n_spins*n_spins);
           dfile << setw(20) << setprecision(8) <<  (( average(3)/(mcs) ) - ((average(4)/(mcs))*(average(4)/(mcs))))/(n_spins*n_spins);
           dfile << setw(20) << setprecision(8) << average(5);
           dfile << setw(20) << setprecision(8) << average(6)/(average(5))*100;
           dfile << setw(20) << setprecision(8) << average(7)/(average(5))*100;
           dfile << setw(20) << setprecision(8) << average(8)/(average(5))*100;
           dfile << setw(20) << setprecision(8) << average(9)/(average(5))*100;
           dfile << setw(20) << setprecision(8) << average(10)/(average(5))*100<<endl;
        }
        dfile.close();


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

    ising_model_simulation();

    return 0;
};
