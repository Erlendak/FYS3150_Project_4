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

void test_temp_1_cold(){

/*
Unit test thats check if our Monte Carlo integration gives answere within reasonable amount to our
analytical solutions, when we use a reasonable amount of Monte Carlo integration points.

We are checking for a 2x2 Latice, with periodic boundary conditions, temperature of 1, using cold start.
Only checking up upon <E>, Cv and <|M|>, because <M> and Chi is quite unstable.
*/

    double expected, tested, tol;
    int n_spins = 2; //L in the LxL latice.
    double temp= 1; // Our temperature.
    double j = 1; // j is scaled to 1.
    double kb = 1; // Boltzmann constant = 1.38064852*pow(10,-23) scaled to 1.
    double beta = 1/(kb*temp);
    double z = 4*cosh(beta*j*8)+12;

    double E_base = (32*j*sinh(8*j*beta))/z; // Analytisk <E>
    double M_base = (16+8*exp(-8*j*beta))/z; // Analytisk <|M|>

    double E = (4*64*j*j*cosh(8*j*beta))/z;
    double E2 = (1024*j*j*(cosh(8*j*beta)*cosh(8*j*beta)))/(z*z);
    double Cv = (1/(kb*temp*temp))*(E-E2); // **  Analytical Cv

    double Chi = (1/(kb*temp))*(32*(1+exp(-8*j*beta)))/z;
    double Chi2 = (1/(kb*temp))*(32*(1+exp(-8*j*beta))/z)-((16+8*exp(-8*j*beta))*(16+8*exp(-8*j*beta)))/(z*z); // Analytisk Susceptebilitet

    // Run simulation
    int mcs = 1000000;
    vec average =isingmodel_cold_start(n_spins, mcs, temp);

// Checking for <E>
   expected = E_base/(n_spins*n_spins);
   tested =average(0)/(mcs*n_spins*n_spins);
   tol = 0.0001;
   try{
       if (abs(abs(tested)-abs(expected)) > tol){
         throw ("Warning: Simulation using the isnigmodel_cold_start does not provide good enough results for <E>  as one should expect, there may be a serious problem in the simulation thats should be check up upon.");
         }
       }
       catch (const char* msg){
            cerr << msg <<endl;
       }

// Checking for Cv
   expected = Cv/(n_spins*n_spins);
   tested = (( average(1)/(mcs) ) - ((average(0)/(mcs))*(average(0)/(mcs))))/(n_spins*n_spins);
   tol = 0.001;
   try{
       if (abs(abs(tested)-abs(expected)) > tol){
         throw ("Warning: Simulation using the isnigmodel_cold_start does not provide good enough results for Cv as one should expect, there may be a serious problem in the simulation thats should be check up upon.");
         }
       }
       catch (const char* msg){
            cerr << msg <<endl;
       }
/*
// Checking for <|M|>
   expected = M_base/(n_spins*n_spins);
   tested = average(4)/(mcs*n_spins*n_spins);
   tol = 0.01;
   try{
       if (abs(abs(tested)-abs(expected)) > tol){
         throw ("Warning: Simulation using the isnigmodel_cold_start does not provide good enough results for <|M|> as one should expect, there may be a serious problem in the simulation thats should be check up upon.");
         }
       }
       catch (const char* msg){
            cerr << msg <<endl;
       }
*/

// Checking for Chi
   expected =Chi2/(n_spins*n_spins);
   tested = (( average(3)/(mcs) ) - ((average(4)/(mcs))*(average(4)/(mcs))))/(n_spins*n_spins);
   tol = 0.01;
   try{
       if (abs(abs(tested)-abs(expected)) > tol){
         throw ("Warning: Simulation using the isnigmodel_cold_start does not provide good enough results for Chi as one should expect, there may be a serious problem in the simulation thats should be check up upon.");
         }
       }
       catch (const char* msg){
            cerr << msg <<endl;
       }

}

void test_temp_1_random(){

/*
Unit test thats check if our Monte Carlo integration gives answere within reasonable amount to our
analytical solutions, when we use a reasonable amount of Monte Carlo integration points.

We are checking for a 2x2 Latice, with periodic boundary conditions, temperature of 1, using random start.
Only checking up upon <E>, Cv and <|M|>, because <M> and Chi is quite unstable.
*/

    double expected, tested, tol;
    int n_spins = 2; //L in the LxL latice.
    double temp= 1; // Our temperature.
    double j = 1; // j is scaled to 1.
    double kb = 1; // Boltzmann constant = 1.38064852*pow(10,-23) scaled to 1.
    double beta = 1/(kb*temp);
    double z = 4*cosh(beta*j*8)+12;

    double E_base = (32*j*sinh(8*j*beta))/z; // Analytisk <E>
    double M_base = (16+8*exp(-8*j*beta))/z; // Analytisk <|M|>

    double E = (4*64*j*j*cosh(8*j*beta))/z;
    double E2 = (1024*j*j*(cosh(8*j*beta)*cosh(8*j*beta)))/(z*z);
    double Cv = (1/(kb*temp*temp))*(E-E2); // **  Analytical Cv

    double Chi = (1/(kb*temp))*(32*(1+exp(-8*j*beta)))/z;
    double Chi2 = (1/(kb*temp))*(32*(1+exp(-8*j*beta))/z)-((16+8*exp(-8*j*beta))*(16+8*exp(-8*j*beta)))/(z*z); // Analytisk Susceptebilitet

    // Run simulation
    int mcs = 1000000;
    vec average =isingmodel_cold_start(n_spins, mcs, temp);

// Checking for <E>
   expected = E_base/(n_spins*n_spins);
   tested =average(0)/(mcs*n_spins*n_spins);
   tol = 0.0001;
   try{
       if (abs(abs(tested)-abs(expected)) > tol){
         throw ("Warning: Simulation using the isnigmodel_random_start does not provide good enough results for <E> as one should expect, there may be a serious problem in the simulation thats should be check up upon.");
         }
       }
       catch (const char* msg){
            cerr << msg <<endl;
       }

// Checking for Cv
   expected = Cv/(n_spins*n_spins);
   tested = (( average(1)/(mcs) ) - ((average(0)/(mcs))*(average(0)/(mcs))))/(n_spins*n_spins);
   tol = 0.001;
   try{
       if (abs(abs(tested)-abs(expected)) > tol){
         throw ("Warning: Simulation using the isnigmodel_random_start does not provide good enough results for Cv as one should expect, there may be a serious problem in the simulation thats should be check up upon.");
         }
       }
       catch (const char* msg){
            cerr << msg <<endl;
       }
/*
// Checking for <|M|>
   expected = M_base/(n_spins*n_spins);
   tested = average(4)/(mcs*n_spins*n_spins);
   tol = 0.01;
   try{
       if (abs(abs(tested)-abs(expected)) > tol){
         throw ("Warning: Simulation using the isnigmodel_cold_random does not provide good enough results for <|M|> as one should expect, there may be a serious problem in the simulation thats should be check up upon.");
         }
       }
       catch (const char* msg){
            cerr << msg <<endl;
       }
*/
// Checking for Chi
   expected =Chi2/(n_spins*n_spins);
   tested = (( average(3)/(mcs) ) - ((average(4)/(mcs))*(average(4)/(mcs))))/(n_spins*n_spins);
   tol = 0.01;
   try{
       if (abs(abs(tested)-abs(expected)) > tol){
         throw ("Warning: Simulation using the isnigmodel_random_start does not provide good enough results for Chi as one should expect, there may be a serious problem in the simulation thats should be check up upon.");
         }
       }
       catch (const char* msg){
            cerr << msg <<endl;
       }

}
void tests(){
/*
Function to call on all the unit tests functions.
*/
    test_temp_1_cold();
    test_temp_1_random();

}

#endif // UNITTESTS_H
