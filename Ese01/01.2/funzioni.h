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

void UniformiSuFile (Random &rnd, const int quanti, string nome, int SOMME){           //Funzione per scrivere numeri casuali uniformi su file output
   ofstream File;
   File.open(nome);
   if (File.is_open()){
      for (int i=0; i<quanti; i++){
         double S_N = 0;
         for (int j=0; j<SOMME; j++){
            S_N += rnd.Rannyu()/SOMME;
         }
         File <<S_N<< endl;
      }
      cout<<"Il file "<<nome<<" è stato stampato"<<endl;
   } else cerr << "PROBLEMA: Impossibile aprire "<<nome<< endl;
  File.close();
  return;
}

void ExpSuFile (Random &rnd, const double lambda, const int quanti, string nome, int SOMME){                //Funzione per scrivere numeri casuali distribuiti esponenzialmente su file output
   ofstream File;
   File.open(nome);
   if (File.is_open()){
      for (int i=0; i<quanti; i++){
         double S_N = 0;
         for (int j=0; j<SOMME; j++){
            S_N += rnd.Exp(lambda)/SOMME;
         }
         File <<S_N<< endl;
      }
      cout<<"Il file "<<nome<<" è stato stampato"<<endl;
   } else cerr << "PROBLEMA: Impossibile aprire "<<nome<< endl;
  File.close();
  return;
}

void LorentzSuFile (Random &rnd, const double mu, const double Gamma, const int quanti, string nome, int SOMME){            //Funzione per scrivere numeri casuali uniformi su file output 
   ofstream File;
   File.open(nome);
   if (File.is_open()){
      for (int i=0; i<quanti; i++){
         double S_N = 0;
         for (int j=0; j<SOMME; j++){
            S_N += rnd.Lorentz(mu, Gamma)/SOMME;
         }
         File <<S_N<< endl;
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
   return;
}

#endif // __funzioni_h__