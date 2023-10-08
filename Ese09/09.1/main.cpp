#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include "file_classi.h"
#include "random.h"

using namespace std;
void Input();
void Stampa_coordinate_citta( const Territorio&, const Popolazione&);

int main(){
  Input();

  Territorio TERR;
  Popolazione POP;
 
  for (int i = 0; i < N_Generazioni; i++){     //Per ogni generazione
    POP.Nuova_Generazione();                   //Proponi crossover e/o mutazioni ai 3/4 peggiori della popolazione
    POP.Distanza(TERR);                           //Calcola la funzione costo per tutti i cromosomi della popolazione
    POP.Riordino();                           //Riordina i cromosomi dal migliore (percorso più corto) al peggiore (percorso più lungo)
    POP.Stampa_Lungh();                       //Stampa la lunghezza del migliore e la media della metà migliore per ogni generazione
    cout<<"------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"GENERAZIONE "<< i+1 <<endl;
    cout<<"------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
  }

  Stampa_coordinate_citta(TERR, POP);              //Stampa le coordinate delle città nell'ordine in cui appaiono nel miglior cromosoma

  return 0;
}

void Input(){
  //GENERATORE NUMERI CASUALI
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
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];         //sets up the RNG
        rnd.SetRandom(seed,p1,p2);
      }
    }
    input.close();
  } else cerr << "PROBLEM: Unable to open seed.in" << endl;

  //LETTURA FILE DI INPUT
  ifstream ReadInput;
  cout << "Traveling Salesman Problem:" << endl;
  cout << "Vincoli:" << endl;
  cout <<"1) La città da cui parto deve essere anche quella in cui arrivo"<<endl;
  cout <<"2) Devo visitare ogni città solo una volta"<<endl<<endl;

  //Read input informations
  ReadInput.open("input.in");

  ReadInput >> tipo_figura;
  ReadInput >> N_Citta;
  ReadInput >> N_Cromosomi;
  ReadInput >> N_Generazioni;
  ReadInput >> pot_selezione;
  ReadInput >> p_crossover;
  ReadInput >> p_mutazioni;

  ReadInput.close();

    //CONTROLLO
  if (N_Citta <= 0){
    cerr << "ERRORE: Hai inserito un numero di città negativo o uguale a zero!" <<endl;
    exit(1);
  }
  if (N_Cromosomi <= 0){
    cerr << "ERRORE: Hai inserito un numero di cromosomi negativo o uguale a zero!" <<endl;
    exit(1);
  }
  if (N_Generazioni <= 0){
    cerr << "ERRORE: Hai inserito un numero di generazioni negativo o uguale a zero!" <<endl;
    exit(1);
  }

  //MESSAGGI DI VERIFICA
  if (tipo_figura == 0){
    cout << endl << "Le tue " << N_Citta << " città sono posizionate su una circonferenza" << endl;
  }
  else{
    cout << "Le tue " << N_Citta << " città sono posizionate all'interno di un quadrato" << endl;
  }
  
  cout << "Hai " << N_Cromosomi << " possibili percorsi (cromosomi) che evolveranno alla ricerca del tragitto migliore"<< endl;
  
  cout << "I cromosomi evolvono per "<< N_Generazioni << " generazioni, poi scelgo il migliore!"<<endl<<endl;

  cout << "La potenza dell'operatore di selezione è posta pari a "<< pot_selezione <<endl;
  cout << "La probabilità di effettuare un crossover è posta pari a "<< p_crossover <<endl;
  cout << "La probabilità di effettuare una mutazione è posta pari a "<< p_mutazioni <<endl<<endl;

  //Calcolo quanto deve fare la somma e la somma dei quadrati del vettore 1D contenente N_Citta tutte 1 sola volta (riutilizzo in check)
  somma_giusta = (N_Citta * (N_Citta - 1)) / 2;                                     //Perché le città vanno da 0 a N_Citta-1
  somma_quadrati_giusta = ((N_Citta-1) * (N_Citta) * (2*(N_Citta - 1)+1)) / 6;          //Perché le città vanno da 0 a N_Citta-1
}
void Stampa_coordinate_citta( const Territorio& terr, const Popolazione& pop ){

  ofstream PrintOutput;

  string RR;
  
  if (tipo_figura == 0){                          
    RR = "1.Circ.coordinate.out"; 
  }
  else{
    RR = "1.Quad.coordinate.out";
  }

  PrintOutput.open(RR);

  cout << "Percorso migliore:" << endl;  

  Cromosoma best_cromosoma_finale = pop.GetCromo(0);     //copia il percorso migliore
  Citta c;

  for (int i = 0; i < N_Citta; i++){            //Stampa le città nell'ordine del cromosoma migliore
    cout << best_cromosoma_finale.GetOrdine(i) << ", " ;    //Stampa a schermo l'ordine del percorso migliore

    c = terr.GetCitta(best_cromosoma_finale.GetOrdine(i) ); //Coordinate della città
    PrintOutput << best_cromosoma_finale.GetOrdine(i) <<"\t" << c.X << "\t" << c.Y << endl;
  }

  PrintOutput.close();
  cout << endl << "Stampa del file: "<<RR<<endl;
}