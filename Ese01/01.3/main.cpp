#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "random.h"
#include "funzioni.h"

using namespace std;

struct needle {
   double y_center;              //y del centro dell'ago
   double y_extr;                //Delta y massimo in modulo dal centro
};

int main (int argc, char *argv[]){

   Random rnd;
   Inizializzazione_random (rnd);

   //Parametri
   double d = 4;             //distanza tra le linee del tavolo
   double L = 1.;              //lunghezza dell'ago
   int Nhit = 0;              //numero di intersezioni tra ago e righe del tavolo

   int Nblocks = 100;              //numero di blocchi
   int Npi = 100;                  //numero di stime di pigreco in ogni blocco
   int Nthrows = 10000;             //numero di tiri per ogni stima di pigreco

   vector<double> ave;       //vettore delle medie di ogni blocco (<Pi>)
   vector<double> ave2;      //vettore delle medie al quadrato di ogni blocco (<Pi>^2)

   vector<double> sum_prog(Nblocks);         //vettore delle somme progressive
   vector<double> sum2_prog(Nblocks);        //vettore delle somme dei quadrati progressive
   vector<double> err_prog(Nblocks);         //vettore delle somme degli errori progressive

   double x, y, r2, Pi, piSum=0, piSum2 = 0;

   //Tiri dell'ago
   needle ndl;

   for (int j=0; j < Nblocks; j++){             //cicli sui blocchi
      piSum =0;
      piSum2 =0;

      for (int i =0; i < Npi; i++){               //stime di pigreco
         Nhit =0;
            
         for (int k=0; k < Nthrows; k++){           //tiri in ogni stima di pigreco
                        
            ndl.y_center = rnd.Rannyu( -0.5*d, 0.5*d );        //the fundamental block of the plane is the spacing d with 0.5d on each side

            if ( abs(ndl.y_center) <= L/2. ){         //se è impossibile che l'ago intersechi, non estraggo la direzione
               do{                     
                  x = rnd.Rannyu();                   //estraggo punto a caso nel primo quarto del cerchio unitario 
                  y = rnd.Rannyu();                   //sin(theta) = y/rad(x^2 + y^2)

                  r2 = x*x + y*y;
               } while (r2 > 1 || r2 == 0);

               ndl.y_extr = L/2. * y/sqrt(r2);           //0 = orizzontale, 0.5*L = verticale
   
               if (ndl.y_center <= 0 && ndl.y_center + ndl.y_extr >= 0)    //centro dell'ago sotto la linea ma estremità superiore è sopra: c'è intersezione!
                  Nhit += 1;
               if (ndl.y_center > 0 && ndl.y_center - ndl.y_extr <= 0)     //centro dell'ago sopra la linea ma estremità inferiore è sotto: c'è intersezione!
                  Nhit += 1;
            }
         }

         if (Nhit > 0){
            Pi = 2*L*Nthrows/(Nhit*d);         //stima di pigreco
            
            piSum += Pi;
            piSum2 += Pi*Pi;

         } else
            cout << "Non c'è stata intersezione nel blocco " << j+1 << endl;
      }
      
      ave.push_back( piSum/Npi );
      ave2.push_back( ave[j]*ave[j] );
   }

   ofstream out;
   string a = "../3.Pigreco.out";
   out.open(a);

   //metodo a blocchi per valori ed errori finali
   for (int i=0; i<Nblocks; i++){

      if (i == 0){
         sum_prog[0] = ave[0];
         sum2_prog[0] = ave2[0];
      } else {
         sum_prog[i] = sum_prog[i-1]*i + ave[i];      
         sum2_prog[i] = sum2_prog[i-1]*i + ave2[i];
      }

      sum_prog[i] = sum_prog[i]/(i+1);           //normalizzazione
      sum2_prog[i] = sum2_prog[i]/(i+1);
      err_prog[i] = error(sum_prog[i], sum2_prog[i], i);

      out << sum_prog[i] << " " << err_prog[i] << endl;
   }

   if (out.is_open()){
    cout<<"Il file "<<a<<" è stato stampato"<<endl;
   } else cerr << "PROBLEMA: Impossibile aprire "<<a<< endl;

   out.close();

   rnd.SaveSeed();
   return 0;
}