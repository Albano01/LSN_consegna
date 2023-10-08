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

   const double S_0 = 100.;         //Prezzo del bene al tempo iniziale
   const double T=1.;               //Tempo totale
   const double K=100.;             //Prezzo d'esercizio (strike price): a quanto si può effettuare il diritto
   const double r=0.1;              //Interesse senza rischio
   const double vol=0.25;           //Volatilità (sigma)

   const int B=100;                    //Numero blocchi
   const int N_ST = pow(10,6);         //Numero valori finali di S_T
   const int L = N_ST/B;               //Numero di estrazioni per blocco

   //Variabili di appoggio
   double C_sum=0;
   double P_sum=0;
   

   vector<double> C (B);         //vettore con medie per ogni blocco
   vector<double> C_sommati (B);      //vettore con somme progressive dei blocchi (ciò che vado a graficare)
   vector<double> C_sigma (B);       //vettore con dev st della media del C_sommati
   
   vector<double> P (B);         //vettore con medie per ogni blocco
   vector<double> P_sommati (B);      //vettore con somme progressive dei blocchi (ciò che vado a graficare)
   vector<double> P_sigma (B);       //vettore con dev st della media del C_sommati

   //-------------------------------------------------------------------------------------------------------------------------------------
   //MODO "UNICO" S_T
   double S_T;

   for (int j=0; j<B; j++){
      C_sum=0;
      P_sum=0;
      for(int i=0; i<L; i++){
         S_T=Calcolo_S_T(rnd, S_0, r, vol, T);

         C_sum += Calcolo_C(r,T,S_T,K);
         P_sum += Calcolo_P(r,T,S_T,K);

      }
      C[j] = C_sum/L;
      P[j] = P_sum/L;
   }

   double C_i =0;
   double P_i=0;
   double C_quadro_sum=0;
   double P_quadro_sum=0;
   double media_C;
   double media_P;
   double media_quadrC;
   double media_quadrP;

   for (int i=0; i<B; i++){
      C_i += C[i];
      C_quadro_sum += pow(C[i],2);
      media_C = C_i/(i+1);
      media_quadrC = C_quadro_sum/(i+1);
      C_sommati[i] = media_C;
      C_sigma[i] = CalcolaDevStMedia(media_quadrC, media_C, i);

      P_i += P[i];
      P_quadro_sum += pow(P[i],2);
      media_P = P_i/(i+1);
      media_quadrP = P_quadro_sum/(i+1);
      P_sommati[i]= media_P;
      P_sigma[i] = CalcolaDevStMedia(media_quadrP, media_P, i);
   }
 
   //Stampo i vettori su file di output
   string c = "../1.C__ST_unico.out";
   VettoriSuFile (C_sommati, C_sigma, c);

   //Stampo i vettori su file di output
   string p = "../1.P__ST_unico.out";
   VettoriSuFile (P_sommati, P_sigma, p);

   //-------------------------------------------------------------------------------------------------------------------------------------
   //MODO DISCRETO S_ti (Riuso gli stessi vettori e variabili d'appoggio)
   double S_ti;
   const int passi = 100;


   for (int j=0; j<B; j++){
      C_sum=0;
      P_sum=0;
      for(int i=0; i<L; i++){
         S_ti=100;
         for (int k=0; k<passi; k++){        //Calcolo gli S_ti per arrivare al S_T in modo discreto
            S_ti = Calcolo_S_t2(rnd, S_ti, r, vol, T/passi);
         }

         C_sum += Calcolo_C(r,T,S_ti,K);
         P_sum += Calcolo_P(r,T,S_ti,K);

      }
      C[j] = C_sum/L;
      P[j] = P_sum/L;
   }

   //Rimetto a zero
   C_i =0;
   P_i=0;
   C_quadro_sum=0;
   P_quadro_sum=0;

   for (int i=0; i<B; i++){
      C_i += C[i];
      C_quadro_sum += pow(C[i],2);
      media_C = C_i/(i+1);
      media_quadrC = C_quadro_sum/(i+1);
      C_sommati[i] = media_C;
      C_sigma[i] = CalcolaDevStMedia(media_quadrC, media_C, i);

      P_i += P[i];
      P_quadro_sum += pow(P[i],2);
      media_P = P_i/(i+1);
      media_quadrP = P_quadro_sum/(i+1);
      P_sommati[i] = media_P;
      P_sigma[i] = CalcolaDevStMedia(media_quadrP, media_P, i);
   }


   //Stampo i vettori su file di output
   string q = "../1.C__ST_discreto.out";
   VettoriSuFile (C_sommati, C_sigma, q);

   //Stampo i vettori su file di output
   string qq = "../1.P__ST_discreto.out";
   VettoriSuFile (P_sommati, P_sigma, qq);
         

   rnd.SaveSeed();
   return 0;
}