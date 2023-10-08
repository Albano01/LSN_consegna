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

struct Posizione{
  double x ;
  double y ;
  double z ;
};

double DistanzaOrigineQuadra (Posizione r){
   return (pow (r.x, 2) + pow (r.y, 2) + pow (r.z, 2));
}

//Nuovo passo discreto
Posizione PassoDiscreto (Random& rnd, Posizione r, double lunghezza){

   int direzione = int ( rnd.Rannyu(0, 3));      //Estraggo direzione ([0]=x, [1]=y, [2]=z)
   double verso = rnd.Rannyu(0,2);             //Estraggo verso ([0]=-, [1]=+)
   
   if(verso<1){
      verso = -1;
   }
   else{verso = 1;}
   
   double passo[3] = {0., 0., 0.}; 
   passo[direzione] = verso*lunghezza;

   r.x += passo[0];
   r.y += passo[1];
   r.z += passo[2];
   
   return r;
}

//Nuovo passo continuo
Posizione PassoContinuo (Random& rnd, Posizione r, double lunghezza){

   double phi = rnd.Rannyu(0, 2* M_PI);            //Estraggo angolo azimutale
   double theta = rnd.Rannyu(0, M_PI);           //Estraggo angolo polare

   r.x += lunghezza * sin(theta) * cos(phi);
   r.y += lunghezza * sin(theta) * sin(phi);
   r.z += lunghezza * cos(theta);
   
   return r;
}

double CalcolaDevStMedia (double m_quadr, double m_media, int n){
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
   return;
}

#endif // __funzioni_h__