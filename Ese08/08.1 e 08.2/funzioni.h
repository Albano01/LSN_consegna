/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __fluid__
#define __fluid__
#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const double xmax=3.;                            //Per istogramma funzione d'onda quadra
const int nbins=100;
double bin_size = 2*xmax/nbins;                 //Campiono da -3 a 3, quindi la lunghezza di ciascuno dei 100 nbins è 0.06    
double Psi2[nbins];           

double walker;                          //In measure e accumulate

// averages
int ACCETTATO, CONTATO;                                 //In SA
double blk_av, blk_norm, accepted, attempted;
double glob_av, glob_av2;
double stima_ham, en_ham, err_ham;                                 //In averages

//configuration (Metodo "Move")
double x;                                       //posizione corrente
double x_old, x_new;                            //posizione vecchia e nuova
double p, psi2_old, psi2_new;                   //probabilità

// thermodynamical state
double sigma, mu;
double temp;                                                            //temperatura fittizia 
vector<double> Energia(2);                                              //Energia ed incertezza
double Delta_sigma, Delta_mu;                                           //Parametri in input per passo su sigma e mu
//Per simulated annealing
double probability, randd;  
double sigma_old, mu_old, sigma_new, mu_new;
vector<double> Energia_Old(2), Energia_New(2);

// simulation
int nstep, nblk, restart;
double Delta_x;
int last = 0;                                               //Mi dice quando posso campionare istogramma
const int wd=12;

//funzioni
void Input(void);
void Simulated_annealing (void);
vector<double> Energia_per_simulated_annealing(void);
void Istogramma_Funz_onda_quadra(void);
void Reset(int);
void Move(void);
double Prob_density(double);
void Measure(void);
void Accumulate(void);
void Averages(int);
double Error(double,double,int);
void ConfFinal(void);
void ConfX(int);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
