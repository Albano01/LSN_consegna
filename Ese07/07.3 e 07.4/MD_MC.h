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

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props, iv, ik, it, ie, iw, ip, ig;
const int nbins=100;
double vtail, ptail, larghezza_bin, rI, rII;
double walker[m_props];

//averages
double blk_av[m_props], blk_norm, accepted, attempted;
double glob_av[m_props], glob_av2[m_props];
double stima_pot, stima_kin, stima_etot, stima_temp, stima_press;
double err_pot, err_kin, err_etot, err_temp, err_gdir, err_press;

//g(r): funzione di distribuzione radiale
double g[nbins];                        //vettore per l'istogramma di g(r) (la campiono in 100 bin)
double glob_av_g[nbins], glob_av2_g[nbins];
double stima_gi;
double err_gi;

//configuration
const int m_part=108;
double x[m_part],    y[m_part],    z[m_part];
double xold[m_part], yold[m_part], zold[m_part];
double vx[m_part],  vy[m_part],   vz[m_part];

// thermodynamical state
int npart;
double beta,temp,energy,vol,rho,box,rcut;

// simulation
int iNVET, nstep, nblk, restart, eq_time;
double delta;

//pigreco
const double pi=3.1415927;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Boltzmann(double, double, double, int);
double Pbc(double);
double Error(double,double,int);
double Force(int, int);

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