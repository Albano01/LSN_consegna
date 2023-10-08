#ifndef __Funzioni__
#define __Funzioni__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;

//Funzioni
double Massimo (double a, double b){
   if (a>=b)
      return a;
   else
      return b;
}

double Calcolo_S_T (Random &rnd, double S_0, double mu, double vol, double T){
   double W = rnd.Gauss (0, T);
   return S_0*exp( (mu-(1./2.)*pow(vol,2))*T + vol*W );
}
double Calcolo_S_t2 (Random &rnd, double S_t1, double mu, double vol, double delta_t){
   double Z = rnd.Gauss (0, 1);
   return S_t1*exp( (mu-(1./2.)*pow(vol,2))*delta_t + vol*Z*sqrt(delta_t) );
}
double Calcolo_C (double r, double T, double S_T, double K){
   return exp(-r*T)*Massimo(0., S_T-K);
}
double Calcolo_P (double r, double T, double S_T, double K){
   return exp(-r*T)*Massimo(0., K-S_T);
}


void VettoriSuFile (vector<double> V, vector<double> T, string nome){         //Funzione per scrivere 2 vector su file output
   ofstream File;
   File.open(nome);
   if (File.is_open()){
      for (int i=0; i<V.size(); i++){
         File <<V[i]<<" "<<T[i]<< endl;
      }
      cout<<"Il file "<<nome<<" è stato stampato"<<endl;
   } else cerr << "PROBLEMA: Impossibile aprire "<<nome<< endl;
  File.close();
  return;
}

double CalcolaDevStMedia (double m_quadr, double m_media, int n){       //Funzione per calcolare dev st della media   
   if(n==0){
      return 0;
   }
   else{
      return sqrt((m_quadr - pow(m_media,2))/n);
   }
}

void Inizializzazione_random (Random &rnd){
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
}

#endif // __Funzioni__

