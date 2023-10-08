#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "funzioni.h"

using namespace std;

int main(){

    //Funzionamento generatore numeri casuali
    Random rnd;
    Inizializzazione_random (rnd);


    //Estrazione di 100.000 numeri in 100 blocchi (L = numeri per blocco = 1000)
    const int M=100000;              //Numero totale di tiri
    const int N=100;                 //Numero di blocchi
    const int L=int(M/N);            //Numero di tiri per ciascun blocco

    vector<double> A;               //Vettore che conterrà le medie
    vector<double> medie_medie;     //Vettore che conterrà le medie delle medie, al crescere di blocchi considerati
    vector<double> sigma_m;         //Vettore che conterrà le deviazioni standard della media del vettore medie_medie

    vector<double> ave_sigma;       //Vettore che conterrà le varianze di ogni blocco
    vector<double> ave_sigma_medie;      //Vettore che conterrà le medie di ave_sigma, al crescere dei blocchi considerati
    vector<double> err_ave_sigma_m;   //Vettore che conterrà gli errori sulle varianze di ogni blocco (la loro dev stand della media)
   
    double r, sum, sum_sigma;       //Variabili d'appoggio

    for (int i=0; i<N; i++){
        sum = 0;
        sum_sigma = 0;
            for (int j=0; j<L; j++){
                r = rnd.Rannyu();
                sum += r;  
                sum_sigma += pow((r-0.5),2);
            }
        A.push_back(sum/L); 
        ave_sigma.push_back(sum_sigma/L);
    }

    double sum_medie = 0;
    double sum_medie_quadro = 0;

    double sum_ave_sigma = 0;
    double sum_ave_sigma_quadro = 0;


    //variabili di appoggio 
    double media_m, media_quadr, MMs, MQs;

    for (int i = 0; i < N; i++){
        sum_medie += A[i];
        sum_medie_quadro += (A[i])*(A[i]);
        media_m = sum_medie/(i+1);
        media_quadr = sum_medie_quadro/(i+1);
        medie_medie.push_back(media_m);
        sigma_m.push_back(CalcolaDevStMedia(media_quadr, media_m, i));


        sum_ave_sigma += ave_sigma[i];
        sum_ave_sigma_quadro += ave_sigma[i]*ave_sigma[i];
        MMs = sum_ave_sigma/(i+1);
        MQs = sum_ave_sigma_quadro/(i+1);
        ave_sigma_medie.push_back(MMs);
        err_ave_sigma_m.push_back(CalcolaDevStMedia(MQs, MMs, i));
    }


    //Stampo i vettori su due file di output
    string c = "../1.medie_medie.out";
    VettoriSuFile (medie_medie, sigma_m, c);

    string k = "../1.ave_sigma_medie.out";
    VettoriSuFile (ave_sigma_medie, err_ave_sigma_m, k);

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //Calcolo il chi quadro da cento blocchi da 10.000 numeri ciascuno, dividendo l'intervallo [0,1) in 100 sottointervalli
    const int tot = 10000;            //Numeri estratti per ogni blocco
    const int blocchi = 100;           //Numero di blocchi
    const int celle = 100;              //Numero di celle in cui divido ogni blocco

    vector<int> conti (celle);            //Vettore di 100 celle in cui tengo conto delle estrazioni in ogni sottointervallo
    vector<double> chi_quadro (blocchi);     //Vettore che contiene i chi quadro da graficare
    
    int y;                              //variabile d'appoggio

    for(int q=0;q<blocchi; q++){
        fill(conti.begin(), conti.end(), 0);
        for(int i=0; i<tot; i++){
                r = rnd.Rannyu();
                y = floor(r*100);       //Prendo la parte intera dei numeri estratti, così da sapere a quale sottointervallo appartengono
                conti[y] ++;
        }
        for(int p=0; p<celle; p++){
           chi_quadro[q] +=  ChiQuadroCella (conti[p], tot, celle);           //Sommo i chi-quadro per ogni cella, per trovare il chi-quadro del blocco q
        }
    }

    //Stampo il vettore chi quadro su file di output
    string t = "../1.chi_quadro.out";
    VettoreSuFile (chi_quadro, t);

    double media_chiquadro =0.;
    for(int i=0; i<blocchi; i++){
        media_chiquadro += chi_quadro[i]/blocchi;
    }
    cout<<"Valor medio del Chi quadro (l'atteso è 100):  "<<media_chiquadro<<endl;

    rnd.SaveSeed();
    return 0;
}