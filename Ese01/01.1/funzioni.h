#ifndef __funzioni_h__
#define __funzioni_h__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include "random.h"

using namespace std;

//Funzioni usate nel main

double CalcolaDevStMedia (double m_quadr, double m_media, int n){       //Funzione per calcolare dev st della media
    
    if(n==0){
        return 0;
    }
    else{
        return sqrt((m_quadr - pow(m_media,2))/n);
    }
}

void VettoriSuFile (vector<double> V, vector<double> T, string nome){      //Funzione per scrivere 2 vector su file output
   
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

void VettoreSuFile (vector<double> V, string nome){               //Funzione per scrivere 1 vector su file output
    
   ofstream File;
   File.open(nome);
   if (File.is_open()){
    for (int i=0; i<V.size(); i++){
      File <<V[i]<<endl;
    }
    cout<<"Il file "<<nome<<" è stato stampato"<<endl;
   } else cerr << "PROBLEMA: Impossibile aprire "<<nome<< endl;
  File.close();
  return;
}

double ChiQuadroCella (int conti, int Numeri_tot_estratti, int nro_celle){
    double A = double (Numeri_tot_estratti)/double(nro_celle);
    return (pow((conti- A),2)/(A));  
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

#endif // __funzioni_h__