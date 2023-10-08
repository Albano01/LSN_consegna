#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"
#include "funzioni.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   Inizializzazione_random (rnd);

   const double lambda = 1.;
   const double mu = 0.;
   const double Gamma = 1.;
   

   const int TOT = pow(10,4);

   const int somme1 = 1;
   const int somme2 = 2;
   const int somme3 = 10;
   const int somme4 = 100;

   string a = "../2.Uniforme_1.out";
   UniformiSuFile (rnd, TOT, a, somme1);

   string b = "../2.Uniforme_2.out";
   UniformiSuFile (rnd, TOT, b, somme2);

   string c = "../2.Uniforme_10.out";
   UniformiSuFile (rnd, TOT, c, somme3);

   string d = "../2.Uniforme_100.out";
   UniformiSuFile (rnd, TOT, d, somme4);

   string e = "../2.Exponential_1.out";
   ExpSuFile (rnd, lambda, TOT, e, somme1);

   string f = "../2.Exponential_2.out";
   ExpSuFile (rnd, lambda, TOT, f, somme2);

   string g = "../2.Exponential_10.out";
   ExpSuFile (rnd, lambda, TOT, g, somme3);

   string h = "../2.Exponential_100.out";
   ExpSuFile (rnd, lambda, TOT, h, somme4);

   string i = "../2.Lorentziana_1.out";
   LorentzSuFile (rnd, mu, Gamma, TOT, i, somme1);

   string j = "../2.Lorentziana_2.out";
   LorentzSuFile (rnd, mu, Gamma, TOT, j, somme2);

   string k = "../2.Lorentziana_10.out";
   LorentzSuFile (rnd, mu, Gamma, TOT, k, somme3);

   string l = "../2.Lorentziana_100.out";
   LorentzSuFile (rnd, mu, Gamma, TOT, l, somme4);

   rnd.SaveSeed();
   return 0;
}

