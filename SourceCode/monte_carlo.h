#ifndef MONTE_CARO_H
#define MONTE_CARO_H
#include <cmath>
#include <iostream>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <lib.h>
#include <string>
#include <armadillo>

using namespace arma;
using namespace std;

inline int periodic(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}

void Counting_Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w, double &count)
{
  // loop over all spins
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      // Find random position
      int ix = (int) (ran1(&idum)*(double)n_spins);
      int iy = (int) (ran1(&idum)*(double)n_spins);
      int deltaE =  2*spin_matrix[iy][ix]*
    (spin_matrix[iy][periodic(ix,n_spins,-1)]+
     spin_matrix[periodic(iy,n_spins,-1)][ix] +
     spin_matrix[iy][periodic(ix,n_spins,1)] +
     spin_matrix[periodic(iy,n_spins,1)][ix]);
      // Here we perform the Metropolis test
      if ( ran1(&idum) <= w[deltaE+8] ) {
    spin_matrix[iy][ix] *= -1;  // flip one spin and accept new spin config
        // update energy and magnetization
        M += (double) 2*spin_matrix[iy][ix];
        E += (double) deltaE;
        count +=1;
      }
    }
  }
} // end of Metropolis sampling over spins

// function to initialise energy, spin matrix and magnetization
void cold_start_initialize(int n_spins, double temp, int **spin_matrix,
        double& E, double& M)
{
  // setup spin matrix and intial magnetization
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      spin_matrix[y][x] = 1; // spin orientation for the ground state
      M +=  (double) spin_matrix[y][x];
    }
  }
  // setup initial energy
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      E -=  (double) spin_matrix[y][x]*
    (spin_matrix[periodic(y,n_spins,-1)][x] +
     spin_matrix[y][periodic(x,n_spins,-1)]);
    }
  }
}// end function initialise

// function to initialise energy, spin matrix and magnetization
void dynamic_initialize(int n_spins, double temp, int **spin_matrix,
        double& E, double& M, double initial_temp)
{
  // setup spin matrix and intial magnetization
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      if (temp <= initial_temp) spin_matrix[y][x] = 1; // spin orientation for the ground state
      M +=  (double) spin_matrix[y][x];
    }
  }
  // setup initial energy
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      E -=  (double) spin_matrix[y][x]*
    (spin_matrix[periodic(y,n_spins,-1)][x] +
     spin_matrix[y][periodic(x,n_spins,-1)]);
    }
  }
}// end function initialise

// function to initialise energy, spin matrix and magnetization
void random_initialize(int n_spins, double temp, int **spin_matrix,
        double& E, double& M)
{
  random_device ran;
  mt19937_64 gen(ran());
  uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
  // setup spin matrix and intial magnetization
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      if ( RandomNumberGenerator(gen)>0.5){
          spin_matrix[y][x] = 1; // spin orientation.
      }
      else {
      spin_matrix[y][x] = -1; // spin orientation.
      }
      M +=  (double) spin_matrix[y][x];
    }
  }
  // setup initial energy
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      E -=  (double) spin_matrix[y][x]*
    (spin_matrix[periodic(y,n_spins,-1)][x] +
     spin_matrix[y][periodic(x,n_spins,-1)]);
    }
  }
}// end function initialise

ofstream ofile;
// inline function for periodic boundary conditions
void read_input(int &n_spins, int &mcs, double &initial_temp,double &final_temp,double &temp_step){
   /* cout << "Size of latice (LxL) L ; ";
    cin >> n_spins;
    cout << "Numbers off Montecarlo integration points ; ";
    cin >> mcs;
    cout << "Numbers off  Initial temperature ; ";
    cin >> initial_temp;
    cout << "Numbers off final temperature ; ";
    cin >> final_temp;
    cout << "Numbers off temperature step ; ";
    cin >> temp_step;
*/
    n_spins = 20;
    mcs = 1000000;
    initial_temp = 1;
    final_temp = 2.4;
    temp_step = 1.4;

};

void output(int n_spins, int mcs, double temp, double *average, double count){
    cout<<"Antall Monte Carlo integrasjonspoeng ; ";
    cout<<mcs<<endl;
    /*cout<<"Antall elektroner som spinner i positiv rettning ; ";
    cout<<n_spins*n_spins<<endl;

    cout<<"Temperatur i simuleringen ; ";
    cout<<temp<<endl;

    cout<<"<e> ; ";
    cout<<average[0]/(mcs*n_spins*n_spins)<<endl;

    cout<<"Cv ; ";
    cout<<(( average[1]/(mcs) ) - ((average[0]/(mcs))*(average[0]/(mcs))))/(n_spins*n_spins)<<endl;


    cout<<"<m> ; ";
    cout<<average[2]/(mcs*n_spins*n_spins)<<endl;

    cout<<"<|m|> ; ";
    cout<<average[4]/(mcs*n_spins*n_spins)<<endl;

    cout<<"Susceptibilitet ; ";
    cout<<(( average[3]/(mcs) ) - ((average[4]/(mcs))*(average[4]/(mcs))))/(n_spins*n_spins)<<endl;
    */

    ofile <<temp;
    ofile << setw(20) << setprecision(8) <<  mcs;
    ofile << setw(20) << setprecision(8) <<  average[0]/(mcs*n_spins*n_spins);
    ofile << setw(20) << setprecision(8) <<  (( average[1]/(mcs) ) - ((average[0]/(mcs))*(average[0]/(mcs))))/(n_spins*n_spins);
    ofile << setw(20) << setprecision(8) << average[2]/(mcs*n_spins*n_spins);
    ofile << setw(20) << setprecision(8) << average[4]/(mcs*n_spins*n_spins);
    ofile << setw(20) << setprecision(8) <<  (( average[3]/(mcs) ) - ((average[4]/(mcs))*(average[4]/(mcs))))/(n_spins*n_spins)<<endl;
    ofile << setw(20) << setprecision(8) <<  count;
};

// Function to read in data from screen
void read_input(int&, int&, double&, double&, double&);
// Function to initialise energy and magnetization
void initialize(int, double, int **, double&, double&);
// The metropolis algorithm
void Metropolis(int, long&, int **, double&, double&, double *);
// prints to file the results of the calculations
void output(int, int, double, double *);

//  main program
void test()
{
  char *outfilename;
  long idum;
  int **spin_matrix, n_spins, ipoint;
  double w[17], average[5], initial_temp, final_temp, E, M, temp_step, count;


  outfilename="test.dat";

  ofile.open(outfilename);
  //    Read in initial values such as size of lattice, temp and cycles
  read_input(n_spins, ipoint, initial_temp, final_temp, temp_step);
  spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
  idum = -1; // random starting point
  for (int mcs = 100; mcs<= ipoint; mcs *= 5){
  for ( double temp = initial_temp; temp <= final_temp; temp+=temp_step){
    //    initialise energy and magnetization
    E = count = M = 0.;
    // setup array for possible energy changes
    for( int de =-8; de <= 8; de++) w[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temp);
    // initialise array for expectation values
    for( int i = 0; i < 5; i++) average[i] = 0.;
    random_initialize(n_spins, (double) temp, spin_matrix, E, M);
    // start Monte Carlo computation
    for (int cycles = 1; cycles <= mcs; cycles++){
      Counting_Metropolis(n_spins, idum, spin_matrix, E, M, w, count);
      // update expectation values
      average[0] += E;    average[1] += E*E;
      average[2] += M;    average[3] += M*M; average[4] += fabs(M);
    }
    cout<< count<<endl;
    output(n_spins, mcs, temp, average, count);
  }}
  free_matrix((void **) spin_matrix); // free memory
  ofile.close();  // close output file

}
#endif // MONTE_CARO_H
