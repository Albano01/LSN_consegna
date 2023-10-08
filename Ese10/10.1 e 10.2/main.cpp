#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "file_classi.h"
#include <mpi.h>                  //Libreria per calcolo parallelo


using namespace std;

//Funzioni
void Input( int );
void Stampa_coordinate_citta( const Territorio&, const Popolazione& );
vector<int>& newSchedule ( int );


int main(int argc, char* argv[]){

  //Parallelizzazione
  int size, rank;
  MPI_Init(&argc, &argv);                   //Inizializzo l'ambiente MPI

  MPI_Comm_size(MPI_COMM_WORLD, &size);     //Dice ai processi quanti sono
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);     //Dice ai processi che size sono [0,...,size-1]

  Input( rank );                      //Inizializza la simulazione
  
  
  Territorio terr ( "American_capitals.dat" );    
  Popolazione pop;

 
  Cromosoma bestchromo;                                                    //variabili e vettori di supporto
  vector<int> send_bestpath(N_Citta), receive_bestpath(N_Citta);
  vector<int> all_bestpaths(size * N_Citta);
  int ricevitore;

  //Simulazione in parallelo con le migrazioni
  for (int imigr = 0; imigr < N_Migrazioni; imigr ++){         //ciclo sulle migrazioni

    for (int igen = 0; igen < N_Gen_statiche; igen ++){   //ciclo sulle generazioni statiche
      pop.Nuova_Generazione();                                     //Proponi crossover e/o mutazioni ai 3/4 peggiori della popolazione
      pop.Distanza( terr );                           //Calcola la funzione costo per tutti i cromosomi della popolazione
      pop.Riordino();                               //Riordina i cromosomi dal migliore (percorso più corto) al peggiore (percorso più lungo)
    }
    
    bestchromo = pop.GetCromo(0);              //ogni continente ha il suo cromosoma migliore
    send_bestpath = bestchromo.GetPercorso();       //estrazione del percorso migliore
    
    schedule = newSchedule( size );             //Nuova schedule creata
    ricevitore = schedule[rank];                  //ogni continente vede chi riceve il suo messaggio

    
    MPI_Sendrecv(send_bestpath.data(), N_Citta, MPI_INT, ricevitore, 0,                                                   //Invio e ricezione
                 receive_bestpath.data(), N_Citta, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    bestchromo.SetPercorso(receive_bestpath);       //Set del percorso migliore ricevuto
    pop.SetCromo( 0, bestchromo );             //ogni continente mette il nuovo cromosoma come migliore
    
  }

  MPI_Barrier(MPI_COMM_WORLD);           //attesa che ogni continente abbia inviato/ricevuto

  pop.Distanza( terr );                //Calcola la funzione costo per tutti i cromosomi della popolazione
  pop.Riordino();                    //Riordina i cromosomi dal migliore (percorso più corto) al peggiore (percorso più lungo)

  bestchromo = pop.GetCromo(0);         //ogni continente ha il suo cromosoma migliore
  send_bestpath = bestchromo.GetPercorso();  //estrazione del percorso migliore


  MPI_Gather(send_bestpath.data(), N_Citta, MPI_INT,   //Invio e ricezione
             all_bestpaths.data(), N_Citta, MPI_INT,
             0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);          //attesa che ogni continente abbia inviato/ricevuto

  
  if (rank == 0){
    N_Cromosomi = size;
    Popolazione bestpop;                           //Continente zero crea una popolazione formata dagli individui migliori

    for(int i = 0; i < size; i++){
      bestchromo = bestpop.GetCromo(i);

      for (int j = 0; j < N_Citta; j++){
        bestchromo.SetOrdine(j, all_bestpaths[i * N_Citta + j]);     //copia i percorsi migliori
      }

      bestpop.SetCromo(i, bestchromo);
      bestpop.GetCromo(i).Check();
    }


    bestpop.Distanza( terr );                //Calcola la funzione costo per tutti i cromosomi della popolazione
    bestpop.Riordino();                     //Riordina i cromosomi dal migliore (percorso più corto) al peggiore (percorso più lungo)

    Stampa_coordinate_citta( terr, bestpop );                   //Stampa le coordinate del migliore!
  }

  MPI_Finalize();

  return 0;
}


//Inizializzo RNG e leggo file di input
void Input( int rank ){

  //reads two primes to set up the random numbers generator (RNG)
  ifstream Primes("Primes");
  if (Primes.is_open()){
    Primes >> p1 >> p2 ;
    for (int i =0; i < rank; i++)
      Primes >> pskip >> pskip;
    Primes >> p3 >> p4;

  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("seed.in");
  string property;
  if (input.is_open()){
    while ( !input.eof() ){
      input >> property;
      if( property == "RANDOMSEED" ){
        input >> world_seed[0] >> world_seed[1] >> world_seed[2] >> world_seed[3];         //sets up the RNG
        world_rnd.SetRandom(world_seed,p1,p2);

        for (int j =0; j < 4; j++){
          seed[j] = world_seed[j] + int(world_rnd.Rannyu(0, 100)) *rank;        //cambia casualmente il seme in funzione del rank
        }
          
        rnd.SetRandom(seed,p3,p4);
      }
    }
    input.close();
  } else cerr << "PROBLEM: Unable to open seed.in" << endl;
  

  ifstream ReadInput;
  ReadInput.open("input.in");

  ReadInput >> N_Cromosomi;        
  ReadInput >> N_Generazioni;            
  ReadInput >> N_Migrazioni;          

  N_Gen_statiche = (int) N_Generazioni/N_Migrazioni;  //Numero di generazioni prima della migrazione

  ReadInput >> pot_selezione;       
  ReadInput >> p_crossover;    
  ReadInput >> p_mutazioni;    

  //CONTROLLO
  if (N_Cromosomi <= 0){
    cerr << "ERRORE: Hai inserito un numero di cromosomi negativo o uguale a zero!" <<endl;
    exit(1);
  }
  if (N_Generazioni <= 0){
    cerr << "ERRORE: Hai inserito un numero di generazioni negativo o uguale a zero!" <<endl;
    exit(1);
  }
  if (N_Migrazioni <= 0){
    cerr << "ERRORE: Hai inserito un numero di migrazioni negativo o uguale a zero!" <<endl;
    exit(1);
  }

  //MESSAGGI DI VERIFICA
  
  cout << "Hai " << N_Cromosomi << " possibili percorsi (cromosomi) che evolveranno alla ricerca del tragitto migliore"<< endl;
  
  cout << "I cromosomi evolvono per "<< N_Generazioni << " generazioni, poi scelgo i migliori!"<<endl;
  cout << "Hai "<< N_Migrazioni << " migrazioni tra i continenti!"<<endl<<endl;

  cout << "La potenza dell'operatore di selezione è posta pari a "<< pot_selezione <<endl;
  cout << "La probabilità di effettuare un crossover è posta pari a "<< p_crossover <<endl;
  cout << "La probabilità di effettuare una mutazione è posta pari a "<< p_mutazioni <<endl<<endl;

  ReadInput.close();
}

//Metodo per creare nuova schedule
vector<int>& newSchedule ( int N_continenti ){
  schedule.clear();                       //pulisci la precedente

  for (int i = 0; i < N_continenti; i++){     
    schedule.push_back(i);                 //crea una nuova schedule ordinata [0,1,2,3,...,size-1] -> ogni continente migra in se stesso
  }

  int pos, appoggio_i;
  for (int i=1; i < N_continenti; i++){                  //disordina casualmente
    
    appoggio_i = schedule[i];                              

    do{ pos = int( world_rnd.Rannyu(0, N_continenti) );     
    } while ( i == schedule[pos] );                    

    schedule[i] = schedule[pos];    
    schedule[pos] = appoggio_i;
  }

  return schedule;
}


void Stampa_coordinate_citta( const Territorio& terr, const Popolazione& pop ){

  ofstream PrintOutput; 
  PrintOutput.open("2.USACapitals.coordinate.out");

  cout << endl << "Miglior percorso: " << endl;  

  Cromosoma bestchromo = pop.GetCromo(0);     //copia il percorso migliore
  Citta c;

  for (int i = 0; i < N_Citta; i++){            //stampa le città nell'ordine detto dal miglior cromosoma
    cout << bestchromo.GetOrdine(i);    
    if (i < N_Citta -1){
      cout << ", " ;
    }
    c = terr.GetCitta( bestchromo.GetOrdine(i) ); 
    PrintOutput << c.X << "\t" << c.Y << endl;
  }

  cout<<endl;

  PrintOutput.close();
}
