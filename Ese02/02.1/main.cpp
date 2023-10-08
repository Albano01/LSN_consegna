#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"
#include "funzioni.h"

using namespace std;

int main (){

   Random rnd;
   Inizializzazione_random (rnd);

   int M = 10000;              //numero tiri per blocco
   int N = 100;                //numero di blocchi

   //CASO DENSITà di PROBABILITà = 1
   vector<double> I (N);      //vettore con stime integrale per ogni blocco
   vector<double> I_sum;      //vettore con somme progressive di I (ciò che vado a graficare)
   vector<double> I_sigma;    //vettore con dev st della media dell'I_sum

   //CASO DENSITà di PROBABILITà = 3/2(1-x^2)
   vector<double> Y (N);      //vettore con stime integrale per ogni blocco
   vector<double> Y_sum;      //vettore con somme progressive di I (ciò che vado a graficare)
   vector<double> Y_sigma;    //vettore con dev st della media dell'I_sum



   //Integro con MonteCarlo con densità di probabilità Uniforme e NON
   for(int i=0; i<N; i++){
      for(int j=0; j<M; j++){
         I[i] += (1./M)*funzione1(rnd.Rannyu());
         Y[i] += (1./M)*funzione2(rnd.Funz2_1());           //Funz2_1 in random.cpp è l'inversa della cumulativa della densità di probabilità NON uniforme
      }
   }

   //Metodo della somma a blocchi
   double sum_I =0;
   double sum_I_quadro =0;
   double media, media_quadr;        //Variabili d'appoggio

   double sum_Y =0;
   double sum_Y_quadro =0;
   double mediaY, media_quadrY;        //Variabili d'appoggio

   for(int i=0; i<N; i++){
      sum_I += I[i];
      sum_I_quadro += pow(I[i],2);
      media = sum_I/(i+1);
      media_quadr = sum_I_quadro/(i+1);
      I_sum.push_back(media);
      I_sigma.push_back(CalcolaDevStMedia(media_quadr, media, i));

      sum_Y += Y[i];
      sum_Y_quadro += pow(Y[i],2);
      mediaY = sum_Y/(i+1);
      media_quadrY = sum_Y_quadro/(i+1);
      Y_sum.push_back(mediaY);
      Y_sigma.push_back(CalcolaDevStMedia(media_quadrY, mediaY, i));
   }

   //Stampo i vettori su file di output
   string c = "../1.I_uniforme.out";
   VettoriSuFile (I_sum, I_sigma, c);


   string t = "../1.I_NONuniforme.out";
   VettoriSuFile (Y_sum, Y_sigma, t);

   cout<<"Programma eseguito!"<<endl;

   rnd.SaveSeed();
   return 0;
}