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

/*
The implementation of the periodic boundary conditions, for the spin latice.
*/


    return (i+limit+add) % (limit);
}





void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w)
{

/*
Implementation of the metropolis algorythm.
*/


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

/*
Alternativ implementation of the metropolis algorythm, we modified this one to count,
both how many times the metropolis algorithm approves the new configuration and what
the change in energy is when approved.
*/


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

        count[0] +=1; // Counts number of approved configurations.


        if((double) deltaE == -8){ // Counts number of approved configurations with energy -8j.
            count[1] += 1;
        }
        if((double) deltaE == -4){ // Counts number of approved configurations with energy -4j.
            count[2] += 1;
        }
        if((double) deltaE == 0){ // Counts number of approved configurations with energy 0j.
            count[3] += 1;
        }
        if((double) deltaE == 4){ // Counts number of approved configurations with energy 4j.
            count[4] += 1;
        }
        if((double) deltaE == 8){ // Counts number of approved configurations with energy 8j.
            count[5] += 1;
        }
      }
    }
  }
} // end of Metropolis sampling over spins





void cold_start_initialize(int n_spins, double temp, int **spin_matrix,
        double& E, double& M)
{

/*
Function to initialise energy, spin matrix and magnetization,
for spesificly a could start configuration.
*/

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





void random_initialize(int n_spins, double temp, int **spin_matrix,
        double& E, double& M)
{

/*
Function to initialise energy, spin matrix and magnetization,
for spesificly a random start configuration.
*/


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





vec isingmodel_cold_start(int n_spins, int mcs, double temp)
{

/*
    Function to solve the two-dimensional Ising model, specificly for the
    cold start configuration. Also for a specific number of Monte Carlo
    points, temperature and Latice size. Periodic boundary conditions is in
    use.

    Assumes that the eaquations has been scaled such as the coupling constant
    J = 1, Boltzmann's constant = 1, temperature has thus dimension energy.
    We are using metropolis sampling, specificly our modified version that captures
    metropolis algorithms approval rate.
*/

  long idum; // Declear memory for the seed  in our random generator function.
  int **spin_matrix; // Declear memory for spinnmatrix
  double w[17],  E, M, average[5], count[6]; // Initialize variables for loop

  ofstream efile; // Initialize file so that we could save the diffrent energy stages to it
  string strtemp = to_string(temp);
  string efilename = "Pe_temp"+strtemp +"_cold.dat";
  efile.open(efilename);
  efile << setiosflags(ios::showpoint | ios::uppercase);
  spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int)); //initialize the spin matrix
  idum = -1; // random starting point

    //  initialise energy and magnetization
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
      average[0] += E;    average[1] += E*E; // Setup <E> and <E^2>
      average[2] += M;    average[3] += M*M; average[4] += fabs(M);  // Setup <M>, <M^2> and <|M|>

      efile << setw(15) << setprecision(8) <<  E<<endl; // Save the energy configurations
    }

  vec ans(11);  //initialize result vector
  for (int i=0; i<5; i++){
  ans(i) = average[i]; // Store <E>, <E^2>, <M>, <M^2> and <|M|>
  }
  for (int i=5; i<11; i++){
  ans(i) = count[i-5];  // Store approved configurations
  }

  free_matrix((void **) spin_matrix); // free memory
  return (ans);

};





vec isingmodel_random_start(int n_spins, int mcs, double temp)
{

/*
    Function to solve the two-dimensional Ising model, specificly for the
    random start configuration. Also for a specific number of Monte Carlo
    points, temperature and Latice size. Periodic boundary conditions is in
    use.

    Assumes that the eaquations has been scaled such as the coupling constant
    J = 1, Boltzmann's constant = 1, temperature has thus dimension energy.
    We are using metropolis sampling, specificly our modified version that captures
    metropolis algorithms approval rate.
*/

  long idum; // Declear memory for the seed  in our random generator function.
  int **spin_matrix; // Declear memory for spinnmatrix
  double w[17],  E, M, average[5], count[6]; // Initialize variables for loop

  ofstream efile; // Initialize file so that we could save the diffrent energy stages to it
  string strtemp = to_string(temp);
  string efilename = "Pe_temp"+strtemp +"_random.dat";
  efile.open(efilename);
  efile << setiosflags(ios::showpoint | ios::uppercase);
  spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int)); //initialize the spin matrix
  idum = -1; // random starting point

  // initialise energy and magnetization
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
    average[0] += E;    average[1] += E*E; // Setup <E> and <E^2>
    average[2] += M;    average[3] += M*M; average[4] += fabs(M); // Setup <M>, <M^2> and <|M|>

    efile << setw(15) << setprecision(8) <<  E<<endl; // Save the energy configurations
  }
  efile.close();

  vec ans(11); //initialize result vector
  for (int i=0; i<5; i++){
  ans(i) = average[i]; // Store <E>, <E^2>, <M>, <M^2> and <|M|>
  }
  for (int i=5; i<11; i++){
  ans(i) = count[i-5]; // Store approved configurations
  }
  free_matrix((void **) spin_matrix); // free memory
  return (ans);

};





ofstream ofile;
void read_input(int &initial_spins, int &final_spins, int &mcs, double &initial_temp,double &final_temp,double &temp_step){

/*
Function to deside what we are simulation in our big main simulation. Tells us how far we want to simulate both in latice
size domane, temperature domane and Monte Carlo integrationpoints domane.
*/

  initial_spins = 40;
  final_spins = 100;
  mcs =  1000000;
  initial_temp = 2;
  final_temp = 2.6;
  temp_step = 0.025;

/*

 // Old implementation of this function.


cout << "Size of latice (LxL) L ; ";
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
};





void output(int n_spins, int mcs, double temp, double *average){

/*
  Function to save our results in our big main simulation. Also this function tells us approximatly
  how far the simulation has come at a given point
*/

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

    ofile <<temp; // Saving temperature
    ofile << setw(20) << setprecision(8) << n_spins ; // Saving size L in the LxL latice
    ofile << setw(20) << setprecision(8) <<  average[0]/(mcs*n_spins*n_spins); // Saving <E>
    ofile << setw(20) << setprecision(8) <<  (( average[1]/(mcs) ) - ((average[0]/(mcs))*(average[0]/(mcs))))/(n_spins*n_spins); // Saving variance <E^2> + <E>^2
    ofile << setw(20) << setprecision(8) << average[2]/(mcs*n_spins*n_spins); // Saving <M>
    ofile << setw(20) << setprecision(8) << average[4]/(mcs*n_spins*n_spins); // Saving <|M|>
    ofile << setw(20) << setprecision(8) <<  (( average[3]/(mcs) ) - ((average[4]/(mcs))*(average[4]/(mcs))))/(n_spins*n_spins)<<endl; // Saving <M^2> + <|M|>^2
};





void ising_model_simulation()
{

/*
  The main simulation from temperature 2.0 to 2.4, with a step size smaller than 0.05, where we are simulating in parallel
  because the simulation itself is very large and therefore taxing on the computer hardware.
*/

  char *outfilename;
  long idum; // Declear memory for the seed  in our random generator function.
  int **spin_matrix, initial_spins, final_spins, mcs; // Declear memory for the spin matrices, initial L, final L and the amount of Monte Carlo integration points.
  double w[17], average[5], initial_temp, final_temp, E, M, temp_step,temp; // Initialize variables for loop

  outfilename="simulation_4e.dat"; // Initialize the file to save our results from the simulation
  ofile.open(outfilename);

  //    Read in initial values such as size of lattice, temp and cycles
  read_input(initial_spins, final_spins,mcs, initial_temp, final_temp, temp_step);


  int samples = int((final_temp-initial_temp)/temp_step)+1; // Total amount of temperature samples per latice.



  for ( int n_spins = initial_spins; n_spins <= final_spins; n_spins+=20){ // loop for the different latices

  // Setup loop for parallel calculation, need to declare wich variables that needs to be cloned to each threads so that do not overwrite each other.
  //#pragma omp parallel for private(w,average,spin_matrix,E,M,temp,idum) schedule(auto)// nowait//shared(n_spins) schedule(static) //nowait //ordered
  for ( int k = 0; k < samples; k+=1){
  //for ( double temp = initial_temp; temp <= final_temp; temp+=temp_step){

    spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int)); //initialize the spin matrix

    //  initialise energy and magnetization
    E = M = 0.;
    idum = (-1- omp_get_thread_num()); // Random starting point
    temp = initial_temp+(temp_step*k); // Calculating the temperature because loop need to be a int to be parallelized.

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
      average[0] += E;    average[1] += E*E;// Setup <E> and <E^2>
      average[2] += M;    average[3] += M*M; average[4] += fabs(M);  // Setup <M>, <M^2> and <|M|>
    }
    //#pragma omp critical // Make sure that the threads do not overwrite each other when saving results.
    output(n_spins, mcs, temp, average); // Save results.
//}

   //free_matrix((void **) spin_matrix); // free memory
  }
}
  ofile.close();  // close output file

}





#endif // MONTE_CARO_H
