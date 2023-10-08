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

double error ( double ave , double ave2, int n){            //Funzione per la deviazione standard della media di n+1 blocchi
    if (n==0)
        return 0;
    else
        return sqrt((ave2 - ave*ave)/n);
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