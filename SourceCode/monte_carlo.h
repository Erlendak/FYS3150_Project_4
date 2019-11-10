#ifndef MONTE_CARO_H
#define MONTE_CARO_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
#include <string>
#include <armadillo>

using namespace arma;

inline int periodic(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}

void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w)
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
      }
    }
  }
}

void initialize(int n_spins, double temp, int **spin_matrix,
        double& E, double& M)
{
  // setup spin matrix and intial magnetization
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      if (temp < 1.5) spin_matrix[y][x] = 1; // spin orientation for the ground state
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
}





//void read_input(int&, int&, double&, double&, double&);
// Function to initialise energy and magnetization
//void initialize(int, double, int **, double&, double&);
// The metropolis algorithm
//void Metropolis(int, long&, int **, double&, double&, double *);
// prints to file the results of the calculations
//void output(int, int, double, double *);


/*
   Program to solve the two-dimensional Ising model
   The coupling constant J = 1
   Boltzmann's constant = 1, temperature has thus dimension energy
   Metropolis sampling is used. Periodic boundary conditions.
*/




vec temperature_integration(int n_spins, int mcs, double temp)
{
  long idum;
  int **spin_matrix;
  double w[17],  E, M;
  vec average(5);


  spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
  idum = -1; // random starting point

    //    initialise energy and magnetization
    E = M = 0.;
    // setup array for possible energy changes
    for( int de =-8; de <= 8; de++) w[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temp);
    // initialise array for expectation values
    for( int i = 0; i < 5; i++) average(i) = 0.;
    initialize(n_spins, temp, spin_matrix, E, M);
    // start Monte Carlo computation
    for (int cycles = 1; cycles <= mcs; cycles++){
      Metropolis(n_spins, idum, spin_matrix, E, M, w);
      // update expectation values
      average(0) += E;    average(1) += E*E;
      average(2) += M;    average(3) += M*M; average(4) += fabs(M);
    }


  return (average);
  free_matrix((void **) spin_matrix); // free memory
};
#endif // MONTE_CARO_H
