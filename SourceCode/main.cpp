#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <unittests.h>
#include <monte_carlo.h>
#include <armadillo>
#include <omp.h>


using namespace  std;
using namespace  arma;



int main()
   {
    //We start with checking if Open MP works, and to tell us how many threads are available.
    int mthreads = omp_get_max_threads();
    printf("There are %d threads available at a time\n",mthreads);

    // Run unit tests.
    test_temp_1();

    /*
    Then we begin our prosses with studying how the Monte Carlo simulations changes with the amount of
    Monte Carlo integrationpoints. Also the way we initialize our spin orintaion wil impact the burn
    in time of the simulation. so we also run to see if a cold start or a random start wil be beneficial.

    For this study we find it natural to just simulate for a cold temperature 1. and a warmer temperature
    2.4, so we are able to see how this would evolve. In order to study the evolution of the simulation,
    we are simulating on a logaritmic basis.
    */

    int n_spins = 20; // This is the L, in LxL latice.
    int samples = 70;  // Samplepoints in each run.
    double logarotmic_base = 1.2; // This is the base we use when calculation how many, Monte Carlo integrationpoint we shall use.
    int shift = 7; // We shift the first samples to start our simulation for a little bigger numbers.
    int mcs;
    vec average; // Result vector to temporarly keep our results.



//  Temperature ;  1.0,      Cold start.

    double temp= 1;
    ofstream afile;
    string afilename = "20x20_cold_start_temp_1.dat";
    afile.open(afilename);
    afile << setiosflags(ios::showpoint | ios::uppercase);

    for (int i=shift;i<samples;i++){
      mcs = pow(logarotmic_base,i); // Get the appropriate number of Monte Carlo integration points
      cout <<  mcs <<endl;
      average = isingmodel_cold_start(n_spins, mcs, temp); // Run simulation with that amount of integrationpoints

      // Writes to file
      afile << setw(15) << setprecision(8) <<  mcs; // Monte Carlo integrationpoints
      afile << setw(20) << setprecision(8) <<  average(0)/(mcs*n_spins*n_spins); // <E>
      afile << setw(20) << setprecision(8) <<  (( average(1)/(mcs) ) - ((average(0)/(mcs))*(average(0)/(mcs))))/(n_spins*n_spins); // Cv
      afile << setw(20) << setprecision(8) << average(2)/(mcs*n_spins*n_spins); // <M>
      afile << setw(20) << setprecision(8) << average(4)/(mcs*n_spins*n_spins); // <|M|>
      afile << setw(20) << setprecision(8) <<  (( average(3)/(mcs) ) - ((average(4)/(mcs))*(average(4)/(mcs))))/(n_spins*n_spins); // Chi
      afile << setw(20) << setprecision(8) << average(5); // Approved Metropolispoints
      afile << setw(20) << setprecision(8) << average(6)/(average(5))*100; // Approved Delta energy = -j8
      afile << setw(20) << setprecision(8) << average(7)/(average(5))*100; // Approved Delta energy = -j4
      afile << setw(20) << setprecision(8) << average(8)/(average(5))*100; // Approved Delta energy = j0
      afile << setw(20) << setprecision(8) << average(9)/(average(5))*100; // Approved Delta energy = j4
      afile << setw(20) << setprecision(8) << average(10)/(average(5))*100<<endl; // Approved Delta energy = j8
      }
    afile.close();




//  Temperature ; 2.4,       Cold start.

    temp = 2.4;
    ofstream bfile;
    string bfilename = "20x20_cold_start_temp_2_4.dat";
    bfile.open(bfilename);
    bfile << setiosflags(ios::showpoint | ios::uppercase);
    for (int j=shift;j<samples;j++){
      mcs = pow(logarotmic_base,j); // Get the appropriate number of Monte Carlo integration points
      cout <<  mcs <<endl;
      bfile << setw(15) << setprecision(8) <<  mcs;
      average = isingmodel_cold_start(n_spins, mcs, temp); // Run simulation with that amount of integrationpoints

      // Writes to file
      bfile << setw(20) << setprecision(8) <<  average(0)/(mcs*n_spins*n_spins); // <E>
      bfile << setw(20) << setprecision(8) <<  (( average(1)/(mcs) ) - ((average(0)/(mcs))*(average(0)/(mcs))))/(n_spins*n_spins); // Cv
      bfile << setw(20) << setprecision(8) << average(2)/(mcs*n_spins*n_spins); // <M>
      bfile << setw(20) << setprecision(8) << average(4)/(mcs*n_spins*n_spins); // <|M|>
      bfile << setw(20) << setprecision(8) <<  (( average(3)/(mcs) ) - ((average(4)/(mcs))*(average(4)/(mcs))))/(n_spins*n_spins); // Chi
      bfile << setw(20) << setprecision(8) << average(5); // Approved Metropolispoints
      bfile << setw(20) << setprecision(8) << average(6)/(average(5))*100; // Approved Delta energy = -j8
      bfile << setw(20) << setprecision(8) << average(7)/(average(5))*100; // Approved Delta energy = -j4
      bfile << setw(20) << setprecision(8) << average(8)/(average(5))*100; // Approved Delta energy = j0
      bfile << setw(20) << setprecision(8) << average(9)/(average(5))*100; // Approved Delta energy = j4
      bfile << setw(20) << setprecision(8) << average(10)/(average(5))*100<<endl; // Approved Delta energy = j8
    }
    bfile.close();




//  Temperature ; 1.0,       Rabdom start.

    temp= 1;
    ofstream cfile;
    string cfilename = "20x20_random_start_temp_1.dat";
    cfile.open(cfilename);
    cfile << setiosflags(ios::showpoint | ios::uppercase);
    for (int i=shift;i<samples;i++){
      mcs = pow(logarotmic_base,i); // Get the appropriate number of Monte Carlo integration points
      cout <<  mcs <<endl;
      cfile << setw(15) << setprecision(8) <<  mcs;
      average = isingmodel_random_start(n_spins, mcs, temp); // Run simulation with that amount of integrationpoints

      // Writes to file
      cfile << setw(20) << setprecision(8) <<  average(0)/(mcs*n_spins*n_spins); // <E>
      cfile << setw(20) << setprecision(8) <<  (( average(1)/(mcs) ) - ((average(0)/(mcs))*(average(0)/(mcs))))/(n_spins*n_spins); // Cv
      cfile << setw(20) << setprecision(8) << average(2)/(mcs*n_spins*n_spins); // <M>
      cfile << setw(20) << setprecision(8) << average(4)/(mcs*n_spins*n_spins); // <|M|>
      cfile << setw(20) << setprecision(8) <<  (( average(3)/(mcs) ) - ((average(4)/(mcs))*(average(4)/(mcs))))/(n_spins*n_spins); // Chi
      cfile << setw(20) << setprecision(8) << average(5); // Approved Metropolispoints
      cfile << setw(20) << setprecision(8) << average(6)/(average(5))*100; // Approved Delta energy = -j8
      cfile << setw(20) << setprecision(8) << average(7)/(average(5))*100; // Approved Delta energy = -j4
      cfile << setw(20) << setprecision(8) << average(8)/(average(5))*100; // Approved Delta energy = j0
      cfile << setw(20) << setprecision(8) << average(9)/(average(5))*100; // Approved Delta energy = j4
      cfile << setw(20) << setprecision(8) << average(10)/(average(5))*100<<endl; // Approved Delta energy = j8
    }
    cfile.close();




//  Temperature ; 2.4,       Random start.

    temp = 2.4;
    ofstream dfile;
    string dfilename = "20x20_random_start_temp_2_4.dat";
    dfile.open(dfilename);
    dfile << setiosflags(ios::showpoint | ios::uppercase);
    for (int j=shift;j<samples;j++){
      mcs = pow(logarotmic_base,j); // Get the appropriate number of Monte Carlo integration points
      cout <<  mcs <<endl;
      dfile << setw(15) << setprecision(8) <<  mcs;
      average = isingmodel_random_start(n_spins, mcs, temp); // Run simulation with that amount of integrationpoints

      // Writes to file
      dfile << setw(20) << setprecision(8) <<  average(0)/(mcs*n_spins*n_spins); // <E>
      dfile << setw(20) << setprecision(8) <<  (( average(1)/(mcs) ) - ((average(0)/(mcs))*(average(0)/(mcs))))/(n_spins*n_spins); // Cv
      dfile << setw(20) << setprecision(8) << average(2)/(mcs*n_spins*n_spins); // <M>
      dfile << setw(20) << setprecision(8) << average(4)/(mcs*n_spins*n_spins); // <|M|>
      dfile << setw(20) << setprecision(8) <<  (( average(3)/(mcs) ) - ((average(4)/(mcs))*(average(4)/(mcs))))/(n_spins*n_spins); // Chi
      dfile << setw(20) << setprecision(8) << average(5); // Approved Metropolispoints
      dfile << setw(20) << setprecision(8) << average(6)/(average(5))*100; // Approved Delta energy = -j8
      dfile << setw(20) << setprecision(8) << average(7)/(average(5))*100; // Approved Delta energy = -j4
      dfile << setw(20) << setprecision(8) << average(8)/(average(5))*100; // Approved Delta energy = j0
      dfile << setw(20) << setprecision(8) << average(9)/(average(5))*100; // Approved Delta energy =  j4
      dfile << setw(20) << setprecision(8) << average(10)/(average(5))*100<<endl; // Approved Delta energy = j8
    }
    dfile.close();




    /*
    Now we run our real simulation,

    */

    ising_model_simulation();

    return 0;
};
