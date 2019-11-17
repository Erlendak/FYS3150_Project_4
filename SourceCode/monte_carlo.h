#ifndef MONTE_CARO_H
#define MONTE_CARO_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
#include <string>
#include <armadillo>
#include <omp.h>

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

void Counting_Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w, double *count)
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
        count[0] +=1;
        if((double) deltaE == -8){
            count[1] += 1;
        }
        if((double) deltaE == -4){
            count[2] += 1;
        }
        if((double) deltaE == 0){
            count[3] += 1;
        }
        if((double) deltaE == 8){
            count[4] += 1;
        }
        if((double) deltaE == 4){
            count[5] += 1;
        }
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




vec isingmodel_cold_start(int n_spins, int mcs, double temp)
{
  long idum;
  int **spin_matrix;
  double w[17],  E, M;
  double average[5];
  double count[6];
  ofstream cfile;
  string strtemp = to_string(temp);
  string cfilename = "Pe_temp"+strtemp +"_cold.dat";
  cfile.open(cfilename);
  cfile << setiosflags(ios::showpoint | ios::uppercase);
  spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
  idum = -1; // random starting point

    //    initialise energy and magnetization
    E = M = 0.;
    // setup array for possible energy changes
    for( int de =-8; de <= 8; de++) w[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temp);
    // initialise array for expectation values
    for( int i = 0; i < 5; i++) average[i] = 0.;
    for( int i = 0; i < 6; i++) count[i] = 0.;
    cold_start_initialize(n_spins, temp, spin_matrix, E, M);
    // start Monte Carlo computation
    for (int cycles = 1; cycles <= mcs; cycles++){
      Counting_Metropolis(n_spins, idum, spin_matrix, E, M, w, count);
      // update expectation values
      average[0] += E;    average[1] += E*E;
      average[2] += M;    average[3] += M*M; average[4] += fabs(M);

      cfile << setw(15) << setprecision(8) <<  E<<endl;
    }
  //cout<< count<<endl;
  vec ans(11);
  for (int i=0; i<5; i++){
  ans(i) = average[i];
  }
  for (int i=5; i<11; i++){
  ans(i) = count[i-5];
  }
  //ans(5) = count;
  free_matrix((void **) spin_matrix); // free memory
  return (ans);

};

vec isingmodel_random_start(int n_spins, int mcs, double temp)
{
  long idum;
  int **spin_matrix;
  double w[17],  E, M;
  double average[5];
  double count[6];
  ofstream cfile;
  string strtemp = to_string(temp);
  string cfilename = "Pe_temp"+strtemp +"_random.dat";
  cfile.open(cfilename);
  cfile << setiosflags(ios::showpoint | ios::uppercase);
  spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
  idum = -1; // random starting point

    //    initialise energy and magnetization
    E = M = 0.;
    // setup array for possible energy changes
    for( int de =-8; de <= 8; de++) w[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temp);
    // initialise array for expectation values
    for( int i = 0; i < 5; i++) average[i] = 0.;
    for( int i = 0; i < 6; i++) count[i] = 0.;
    random_initialize(n_spins, temp, spin_matrix, E, M);
    // start Monte Carlo computation
    for (int cycles = 1; cycles <= mcs; cycles++){
      Counting_Metropolis(n_spins, idum, spin_matrix, E, M, w, count);
      // update expectation values
      average[0] += E;    average[1] += E*E;
      average[2] += M;    average[3] += M*M; average[4] += fabs(M);

      cfile << setw(15) << setprecision(8) <<  E<<endl;
    }
  //cout<< count<<endl;
  vec ans(11);
  for (int i=0; i<5; i++){
  ans(i) = average[i];
  }
  for (int i=5; i<11; i++){
  ans(i) = count[i-5];
  }
  //ans(5) = count;
  free_matrix((void **) spin_matrix); // free memory
  return (ans);

};




//
ofstream ofile;
// inline function for periodic boundary conditions
void read_input(int &initial_spins, int &final_spins, int &mcs, double &initial_temp,double &final_temp,double &temp_step){
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
    initial_spins = 40;
    final_spins = 100;
    mcs = 100000;
    initial_temp = 2;
    final_temp = 2.6;
    temp_step = 0.01;

};

void output(int n_spins, int mcs, double temp, double *average){

    cout<<"Antall elektroner som spinner ; ";
    cout<<n_spins*n_spins<<endl;

    cout<<"Temperatur i simuleringen ; ";
    cout<<temp<<endl;
    /*
    cout<<"Antall Monte Carlo integrasjonspoeng ; ";
    cout<<mcs<<endl;
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
    ofile << setw(20) << setprecision(8) << n_spins ;
    ofile << setw(20) << setprecision(8) <<  average[0]/(mcs*n_spins*n_spins);
    ofile << setw(20) << setprecision(8) <<  (( average[1]/(mcs) ) - ((average[0]/(mcs))*(average[0]/(mcs))))/(n_spins*n_spins);
    ofile << setw(20) << setprecision(8) << average[2]/(mcs*n_spins*n_spins);
    ofile << setw(20) << setprecision(8) << average[4]/(mcs*n_spins*n_spins);
    ofile << setw(20) << setprecision(8) <<  (( average[3]/(mcs) ) - ((average[4]/(mcs))*(average[4]/(mcs))))/(n_spins*n_spins)<<endl;
};



void ising_model_simulation()
{
  char *outfilename;
  long idum;
  int **spin_matrix, initial_spins, final_spins, mcs;
  double w[17], average[5], initial_temp, final_temp, E, M, temp_step,temp;

 //ofstream ofile;
  outfilename="simulation_4e.dat";

  ofile.open(outfilename);
  //    Read in initial values such as size of lattice, temp and cycles
  read_input(initial_spins, final_spins,mcs, initial_temp, final_temp, temp_step);
  int _n = int((final_temp-initial_temp)/temp_step)+1;




  for ( int n_spins = initial_spins; n_spins <= final_spins; n_spins+=20){

      // random starting point

  #pragma omp parallel for  private(w,average,spin_matrix,E,M,temp,idum) shared(n_spins) //ordered
  for ( int k = 0; k <= _n; k+=1){
  //for ( double temp = initial_temp; temp <= final_temp; temp+=temp_step){
    //    initialise energy and magnetization
    spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
    E = M = 0.;
    idum = -1;
    //cout<<omp_get_thread_num()<<endl;
    temp = initial_temp+(temp_step*k);
    //int nthreads = omp_get_num_threads();
    //printf("Using %d threads outside parallel loop\n",nthreads);

    // setup array for possible energy changes
    for( int de =-8; de <= 8; de++) w[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temp);
    // initialise array for expectation values
    for( int i = 0; i < 5; i++) average[i] = 0.;
    random_initialize(n_spins, (double) temp, spin_matrix, E, M);
    // start Monte Carlo computation
    for (int cycles = 1; cycles <= mcs; cycles++){
      Metropolis(n_spins, idum, spin_matrix, E, M, w);
      // update expectation values
      average[0] += E;    average[1] += E*E;
      average[2] += M;    average[3] += M*M; average[4] += fabs(M);
    }
    #pragma omp critical //ordered
    output(n_spins, mcs, temp, average);
}

   //free_matrix((void **) spin_matrix); // free memory
  }

  ofile.close();  // close output file

}


#endif // MONTE_CARO_H
