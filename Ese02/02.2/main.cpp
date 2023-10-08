#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "funzioni.h"

using namespace std;

int main (){

   Random rnd;
   Inizializzazione_random (rnd);


   const int step = 100;                    //Numero passi per random walk
   const int M_RW = 10000;                 //Numero simulazioni random walks
   const int B = 100;                     //Numero di blocchi in cui divido i random walks
   int L = int(M_RW/B);                  //Numero di random walks per blocco
   const double a = 1.;                 //Lunghezza del passo

   //CASO DISCRETO
   vector <double> D_SommaRquad(step, 0.0);                          //Vectors per metodo a blocchi, Inizializzati con zeri 
   vector <double> D_MediaBlocchi(step, 0.0);
   vector <double> D_MediaQuadrBlocchi(step, 0.0); 
   double D_DevStMedia;

   Posizione R_D {0. , 0. , 0.};

   ofstream D_out;                                                   //File di output
   string nome1 = "../2.RW_Discreto.out";
   D_out.open(nome1);

   //CASO CONTINUO
   vector<double> C_SommaRquad(step, 0.0);                               //Vectors per metodo a blocchi, Inizializzati con zeri 
   vector <double> C_MediaBlocchi(step, 0.0);
   vector <double> C_MediaQuadrBlocchi(step, 0.0); 
   double C_DevStMedia;

   Posizione R_C {0. , 0. , 0.};

   ofstream C_out;                                                   //File di output
   string nome2 = "../2.RW_Continuo.out";
   C_out.open(nome2);

   for (int k = 0; k < B; k++){                 //Ciclo sui blocchi
      for (int i = 0; i < L; i++){              //Ciclo sulle RWs

         R_D = {0., 0., 0.};                    //Origine come punto di partenza per ogni RW
         R_C = {0., 0., 0.};

         for (int j = 0; j < step; j++){        //Ciclo sui passi

            R_D = PassoDiscreto( rnd, R_D, a);      //Nuovo passo
            R_C = PassoContinuo( rnd, R_C, a);

            D_SommaRquad[j] += DistanzaOrigineQuadra( R_D );             //Distanza dall'origine al quadrato ad ogni passo: sommo le corrispondenti per tutte le RWs di un blocco
            C_SommaRquad[j] += DistanzaOrigineQuadra( R_C );             
         }
      }
      
      for (int j = 0; j < step; j++){

         D_MediaBlocchi[j] += sqrt( D_SommaRquad[j]/L );    //Media
         C_MediaBlocchi[j] += sqrt( C_SommaRquad[j]/L );    

         D_MediaQuadrBlocchi[j] += D_SommaRquad[j]/L;     //Quadrato della media
         C_MediaQuadrBlocchi[j] += C_SommaRquad[j]/L;     
         
         D_SommaRquad[j] = 0.;                  //Reset prima di cambiare blocco
         C_SommaRquad[j] = 0.;
      }
   }

   for (int j=0; j<step; j++){

      D_MediaBlocchi[j] /= B;                //Media sui blocchi
      D_MediaQuadrBlocchi[j] /= B;               //Media al quadrato sui blocchi
               
      D_DevStMedia = CalcolaDevStMedia (D_MediaQuadrBlocchi[j], D_MediaBlocchi[j], B-1);

      D_out << D_MediaBlocchi[j] << "\t" << D_DevStMedia << endl;


      C_MediaBlocchi[j] /= B;                
      C_MediaQuadrBlocchi[j] /= B;               
               
      C_DevStMedia = CalcolaDevStMedia (C_MediaQuadrBlocchi[j], C_MediaBlocchi[j], B-1);

      C_out << C_MediaBlocchi[j] << "\t" << C_DevStMedia << endl;
   }
   
   if (D_out.is_open()){
      cout<<"Il file "<<nome1<<" è stato stampato"<<endl;
   } else cerr << "PROBLEMA: Impossibile aprire "<<nome1<< endl;
   if (C_out.is_open()){
      cout<<"Il file "<<nome2<<" è stato stampato"<<endl;
   } else cerr << "PROBLEMA: Impossibile aprire "<<nome2<< endl;

   D_out.close();
   C_out.close();

   rnd.SaveSeed();
   return 0;
}