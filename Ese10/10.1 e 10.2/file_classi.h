#ifndef __fluid__
#define __fluid__

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <ostream>
#include <iomanip>
#include "random.h"

using namespace std;

//Variabili globali ed inizializzazione funzioni
int seed[4];
int world_seed[4];
int p1, p2, p3, p4, pskip;;
Random rnd, world_rnd;

int tipo_figura;                     //se è zero posiziono le città su una circonferenza, se no dentro un quadrato
int N_Citta;                          //quante città posiziono (ovvero devo visitare)
int N_Cromosomi;                        //numero di cromosomi (ovvero percorsi) che evolvo. Costituiscono la popolazione
int N_Generazioni;                         //numero di generazioni dopo le quali scelgo il cromosoma migliore
int N_Migrazioni;                            //numero di migrazioni tra continenti
int N_Gen_statiche;                             //numero di generazioni prima della migrazione

double pot_selezione, p_crossover, p_mutazioni;       //probabilità

double check_somma, check_somma_quadrati, somma_giusta, somma_quadrati_giusta;           //Per funzione di check
double somma_progressiva;                                                                //Per funzione calcolo distanza

vector<int> Sequenza_padre;                                         //Per crossover
vector<int> Sequenza_madre;
vector<int> schedule;

struct Citta {                                                  //struct con coordinate di ciascuna città
    double X,Y;
};

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
//CLASSI
class Territorio {

 public:

    //Costruttore
    Territorio(){                  //Costruttore che posiziona le città: se tipo_figura == 0 sul cerchio, else dentro quadrato di lato 1 con vertice in (0,0)
        m_figura = tipo_figura;
        for (int i=0; i<N_Citta; i++){
            if (m_figura == 0){
               m_angolo = rnd.Rannyu(0, 2*M_PI);
               m_Posizioni.push_back({cos(m_angolo), sin(m_angolo)});
            }
            else{
                m_Posizioni.push_back({rnd.Rannyu(), rnd.Rannyu()});
            }
        }
    }

    //Costruttore che prende coordinate città da file esterno
    Territorio(const string& filename) {
        double lon, lat;                    //longitudine e latitudine

        ifstream File_Citta;
        File_Citta.open(filename);

        if (!File_Citta) {
            cerr << "Error: File not found." << endl;
            exit(-11);
        }

        while (File_Citta >> lon >> lat) {                  //Legge le coordinate dal file finché non finiscono...
            m_Posizioni.push_back({lon, lat});
            N_Citta++;                                      //...aumentando il contatore del numero di città
        }

        cout<<"Il tuo territorio ha "<<N_Citta<<" città!"<<endl;
        //Calcolo quanto deve fare la somma e la somma dei quadrati del vettore 1D contenente N_Citta tutte 1 sola volta (riutilizzo in check)
        somma_giusta = (N_Citta * (N_Citta - 1)) / 2;                                     //Perché le città vanno da 0 a N_Citta-1
        somma_quadrati_giusta = ((N_Citta-1) * (N_Citta) * (2*(N_Citta - 1)+1)) / 6;          //Perché le città vanno da 0 a N_Citta-1

        File_Citta.close();
    }

    //Metodi
    Citta GetCitta( int pos ) const {                                 //Metodo che restituisce città, dato il suo intero di riferimento (con controllo)
        if (pos < N_Citta )
            return m_Posizioni[pos];
        else{
            cerr << "Errore: chiamata città NON esistente!" << endl;
            exit(-3);
        }
    }

 protected:
    int m_figura;
    double m_angolo;
    vector<Citta> m_Posizioni;
};
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
class Cromosoma{

 public:

    Cromosoma(){                               //Costruttore (crea la generazione zero: dapprima creo cromosoma (percorso) ordinato...
        for (int i = 0; i < N_Citta; i++) {
            m_Ordine.push_back(i);
        }

        int appoggio1, appoggio2;                   //...poi lo mischio casualmente)
        for (int j = 1; j<N_Citta; j++){
            appoggio1 = (int)(rnd.Rannyu(1, N_Citta));
            appoggio2 = m_Ordine[j];
            m_Ordine[j] = m_Ordine[appoggio1]; 
            m_Ordine[appoggio1] = appoggio2;
        }

        m_distanza_percorso = 0.;                          //Inizializza la lunghezza del percorso a 0 (dopo la calcolo!)
    } 

    //Metodi 
    void Check(){                                   //Verifico che la città di partenza sia sempre la stessa e che non ci sono ripetizioni
        if (m_Ordine[0] != 0){
            cerr<<"ERRORE: UN CROMOSOMA HA LA CITTÀ INIZIALE CAMBIATA!"<<endl;
            exit(-1);
        }
        check_somma = 0;
        check_somma_quadrati = 0;                               //Azzero le somme

        for (int i=0; i<N_Citta; i++){
            check_somma += m_Ordine[i];
            check_somma_quadrati += pow(m_Ordine[i], 2);
        }
        
        if(check_somma != somma_giusta || check_somma_quadrati != somma_quadrati_giusta){           //Verifico che corrispondano a quelle attese (calcolate in Input())
            cerr<<"ERRORE: UN CROMOSOMA HA DELLE CITTÀ RIPETUTE!"<<endl;
            cerr<<"Infatti, somma ottenuta: "<< check_somma << "\t" << "somma corretta: "<<somma_giusta<<endl;
            cerr<<"Somma quadrati ottenuta: "<< check_somma_quadrati << "\t" << "somma quadrati corretta: "<<somma_quadrati_giusta<<endl;
            exit(-1);
        }
    }

    
    void SetDistanza(double L){
      m_distanza_percorso = L;
    }

    
    double GetDistanza() const{
      return m_distanza_percorso;
    }

    
    int GetOrdine(int pos) const{               
      return m_Ordine[pos];
    }

    vector<int> GetPercorso() const{
      return m_Ordine;
    }

    void SetPercorso(vector<int> newpath) {
      m_Ordine = newpath;
      Check();
    }


    void SetOrdine(int pos, int new_ord){
      m_Ordine[pos] = new_ord;
    }

    
    Cromosoma& operator= (const Cromosoma& other) {                                 //overriding dell'operatore "=" affinché possa copiare cromosomi
      if (this != &other) {
        m_distanza_percorso = other.GetDistanza();      
      }

      for(int i = 1; i<N_Citta; i++){
        m_Ordine[i] = other.GetOrdine(i);       
      }
      return *this;
    }

    //GENETIC-MUTATION OPERATORS
    void Shuffle(){                                 //Scambio due città (a e b) casualmente (ECCETTO la città di partenza)
        int a = (int)(rnd.Rannyu(1, N_Citta));
        int b = (int)(rnd.Rannyu(1, N_Citta));
        int appoggio = m_Ordine[a];
        m_Ordine[a] = m_Ordine[b]; 
        m_Ordine[b] = appoggio;
    }
    void Trasla_a_blocchi() {
      vector<int> Old_Cromo;                    //support cromosome

      for (int i=0; i<N_Citta; i++){
        Old_Cromo.push_back(m_Ordine[i]);        //saves the old chromosome
      }

      int m = int (rnd.Rannyu(1, N_Citta-1));   //extracts the number of cities to be shifted, from 1 (only the city in position 1) to Ncities-2 (every city but the first and last)
      int n = int (rnd.Rannyu(1, N_Citta-m));   //extracts the shift length, from 1 (only the next city goes before) to Ncities-m (all the next cities go before)

      for (int i=1; i<=m+n; i++){               //only the cities in positions [1,m+n] are affected from a shift
        if (i <= m)                             
          m_Ordine[i+n] = Old_Cromo[i];          //the city labels in positions [1,m] are copied in the forwards shifted positions [1+n,m+n]
        else
          m_Ordine[i-m] = Old_Cromo[i];          //the city labels in positions [m+1, Ncities-1] are copied in the backwards shifted positions [1, Ncities-1 -m]
      }
    }



    void Shuffle_a_segmenti(){
      int dim = int (rnd.Rannyu(2, N_Citta/2));         //Dimensione del blocco da scambiare
      int j = int (rnd.Rannyu(dim+1, N_Citta-dim+1));     //Posizione del secondo blocco

      int appoggio;
      for (int i = 1; i <= dim; i++){
        appoggio = m_Ordine[i];                  

        m_Ordine[i] = m_Ordine[j];               
        m_Ordine[j] = appoggio;
        j++;                                   
      }    
    }
    void Inversione(){                              //Inverto l'ordine di m città contigue (ESCLUSA la città di partenza)
        //Indici del segmento da invertire
        int indiceInizio = (int)rnd.Rannyu(1, N_Citta-1);
        int indiceFine = (int)rnd.Rannyu(indiceInizio+1, N_Citta);

        //Inversione del segmento
        reverse(m_Ordine.begin() + indiceInizio, m_Ordine.begin() + indiceFine);       
    }


 protected:
    vector<int> m_Ordine;            //Vettore contenente l'ordine in cui il viaggiatore visita le città, il primo (zero) è fisso
    double m_distanza_percorso;
};
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------

class Popolazione : public Cromosoma, public Territorio{

 public:

    
    Popolazione(){                                              //Costruttore (crea la popolazione: insieme di cromosomi)
        for (int i = 0; i < N_Cromosomi; i++) {
            Cromosoma CROMO;                                //Creo un individuo (un cromosoma della classe cromosoma)
            m_Popolazione.push_back(CROMO);
        }
    }

    Cromosoma GetCromo(int pos) const{                      //Ritorna il cromosoma nella posizione pos nel vettore popolazione
        return m_Popolazione[pos];
    }

    void SetCromo(int pos, Cromosoma& newchromo ){
        m_Popolazione[pos] = newchromo;
    }

    //Metodi
    void Distanza( const Territorio& terr){
        for(int j=0; j < N_Cromosomi; j++){          //ciclo sui cromosomi nella popolazione
            somma_progressiva = 0.;                   

            for (int i = 0; i < N_Citta-1; i++){    //ciclo sulle città nel cromosoma
                somma_progressiva += Distanza_tra_due_citta( terr.GetCitta( m_Popolazione[j].GetOrdine(i) ), terr.GetCitta( m_Popolazione[j].GetOrdine(i+1) ) );        
            }

            somma_progressiva += Distanza_tra_due_citta( terr.GetCitta( m_Popolazione[j].GetOrdine(N_Citta-1) ), terr.GetCitta( m_Popolazione[j].GetOrdine(0) ) );    //Sommo separatamente la distanza dell'ultima con la prima
            m_Popolazione[j].SetDistanza(somma_progressiva);       //imposta la m_ditanza_percorso al cromosoma j-esimo
        }
    }

    double Distanza_tra_due_citta( Citta A, Citta B ){                  //Calcola la distanza tra due città
        return sqrt( pow((A.X - B.X), 2) + pow((A.Y - B.Y), 2) );
    }

    static bool Comparatore(const Cromosoma& a, const Cromosoma& b) {              //Restituisce true se la distanza del primo cromosoma è minore della distanza del secondo
        return a.GetDistanza() < b.GetDistanza();
    }   

    //Ordina dal migliore (distanza più piccola) al peggiore (distanza più grande)
    void Riordino(){
        sort(m_Popolazione.begin(), m_Popolazione.end(), Comparatore);                 //Ordina il vettore popolazione in base al data membro 'm_distanza_percorso'
    }

    int Operatore_di_Selezione(int pot){                                                       //Seleziona un cromosoma nella popolazione ORDINATA, ritornando la sua posizione
        //IL VETTORE DEVE ESSERE ORDINATO!
        //Se pot<1 favorisco le prime posizioni (quindi i cromosomi migliori), se pot>1 il contrario
        return (int)(N_Cromosomi * pow(rnd.Rannyu(), pot));
    }


    //ALGORITMO DI CROSSOVER
    void Crossover(int bad_cromosoma1, int bad_cromosoma2) {
        int madre = Operatore_di_Selezione(pot_selezione);                                  //seleziona madre e padre
        int padre = Operatore_di_Selezione(pot_selezione);                                

        Cromosoma MaCromo = m_Popolazione[madre];                      //Variabili di appoggio per la madre ed il padre, che vengono copiati
        Cromosoma PaCromo = m_Popolazione[padre];                 

        int cut = int(rnd.Rannyu(1, N_Citta - 1));              //taglio dell'algoritmo (decide dove tagliare i cromosomi madre e padre)

        Sequenza_madre.clear();                                      //azzera vettori globali
        Sequenza_padre.clear();                                     

        //Metto in sequenza_madre e sequenza_padre le posizioni delle città nel nuovo ordine reciproco per le posizioni dopo il taglio
        for (int i = cut; i < N_Citta; i++) {                  //Ciclo sulle città tagliate

            for (int j = 1; j < N_Citta; j++) {                  
                if (Sequenza_madre.size() < (unsigned int)(i-cut +1)){                   //la ricerca termina quando l'i-esima città della madre è stata trovata
                    if ( MaCromo.GetOrdine(i) == PaCromo.GetOrdine(j) ) {
                        Sequenza_madre.push_back(j);                         
                    }
                }
            
                if (Sequenza_padre.size() < (unsigned int)(i-cut +1)){                   //la ricerca termina quando l'i-esima città del padre è stata trovata
                    if ( PaCromo.GetOrdine(i) == MaCromo.GetOrdine(j) ) { 
                        Sequenza_padre.push_back(j);                       
                    }
                 }
            }
        }

        //CONTROLLI
        if (Sequenza_madre.size() != Sequenza_padre.size()){             
            cerr << "ERRORE: differenti dimensioni nel crossover" << endl;
            exit(-10);
        }
        if (Sequenza_madre.size() != (unsigned int)(N_Citta - cut)){              
            cerr << "ERRORE: dimensione non corretta nel crossover" << endl;
            exit(-11);
        }

        //Riordino 
        sort(Sequenza_madre.begin(), Sequenza_madre.end());   
        sort(Sequenza_padre.begin(), Sequenza_padre.end());

        //Sostituisco nell'ordine trovato nel partner
        for (unsigned int i = 0; i<Sequenza_madre.size(); i++){
            Sequenza_madre[i] = PaCromo.GetOrdine(Sequenza_madre[i]);
            Sequenza_padre[i] = MaCromo.GetOrdine(Sequenza_padre[i]);
        }

   

        //Metto i figli al posto dei cromosomi "cattivi"
        int k = 0;
        for (int i = 1; i < N_Citta; i++) {

            if (i < cut) {                                        //Prima del taglio copio dai genitori
                m_Popolazione[bad_cromosoma1].SetOrdine(i, MaCromo.GetOrdine(i));
                m_Popolazione[bad_cromosoma2].SetOrdine(i, PaCromo.GetOrdine(i));
            } else {                                              //Dopo il taglio copio dai figli
                m_Popolazione[bad_cromosoma1].SetOrdine(i, (Sequenza_madre[k]));
                m_Popolazione[bad_cromosoma2].SetOrdine(i, (Sequenza_padre[k]));
                k++;
            }
        }
    }


    
    void Scegli_Mutazione( int m, double prob_rand ){             //Muto il cromosoma numerato m confrontando una probabilità estratta tra [0,prob_muta) in Nuova_Generazione
        if (prob_rand < p_mutazioni/2.){
            if (prob_rand < p_mutazioni/4.){                                        //Se prob_rand in [0, p_mutazioni/4) -> SHUFFLE
                m_Popolazione[m].Shuffle();
            } else {                                                                //Se prob_rand in [p_mutazioni/4, p_mutazioni/2) -> TRASLAZIONE
                m_Popolazione[m].Trasla_a_blocchi();
            }
        } else{
            if (prob_rand < 3.*p_mutazioni/4.){                                     //Se prob_rand in [p_mutazioni/2, 3*p_mutazioni/4) -> SHUFFLE_A_SEGMENTI
                m_Popolazione[m].Shuffle_a_segmenti();
            } else {                                                                //Se prob_rand in [3*p_mutazioni/4, p_mutazioni) -> INVERSIONE
                m_Popolazione[m].Inversione();
            }
        }
    }


  void Nuova_Generazione(){                                   //Metodo che crea nuova generazione proponendo crossover e/o mutazione ai 3/4 della popolazione esistente
    double r_mutazione1, r_mutazione2, r_crossover;
    

    for (int i = int(N_Cromosomi/4); i<N_Cromosomi-1; i=i+2){     //ciclo sui pari "brutti" (tolgo il primo quarto, parte bella della popolazione)
      r_crossover = rnd.Rannyu();

      if (r_crossover < p_crossover){           //Se estratto è minore della probzbilità di Crossover...
        Crossover(i, i+1);                  //...Crossover accettato!
      }

      r_mutazione1 = rnd.Rannyu();  
      r_mutazione2 = rnd.Rannyu();

      if (r_mutazione1 < p_mutazioni){
        Scegli_Mutazione( i, r_mutazione1 );       //Mutazione per cromosoma i
      }
      if (r_mutazione2 < p_mutazioni){
        Scegli_Mutazione( i+1, r_mutazione2 );     //Mutazione per cromosoma i+1
      }

      m_Popolazione[i].Check();             //Funzioni di controllo sui nuovi cromosomi
      m_Popolazione[i+1].Check();
    }
  }

  void Stampa_Lungh(){
    ofstream LunghezzaOutput;

    string EE;
  
    if (tipo_figura == 0){                          
        EE = "1.Circ.lunghezza.out"; 
    }
    else{
        EE = "1.Quad.lunghezza.out";
    }

    LunghezzaOutput.open(EE, ios::app);

    double Lunghezza_migliore_parziale;
    double Media_meta_lunghezze_migliori_parziale = 0.;
    int k = 0;

    Lunghezza_migliore_parziale = m_Popolazione[0].GetDistanza();
    for (int i = 0; i<N_Citta/2.; i++){
        Media_meta_lunghezze_migliori_parziale += m_Popolazione[i].GetDistanza();
        k++;
    }
    Media_meta_lunghezze_migliori_parziale = Media_meta_lunghezze_migliori_parziale/k;

    LunghezzaOutput << Lunghezza_migliore_parziale <<"\t" << Media_meta_lunghezze_migliori_parziale << endl;
    
    LunghezzaOutput.close();
  }


 protected:
    vector<Cromosoma> m_Popolazione;
};



#endif