#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include "funzioni.h"

using namespace std;

int main(){ 
  Input(); 							                              //Inizialization

  Simulated_annealing();                                //Algoritmo di ottimizzazione per trovare mu e sigma del minimo (stampo i due file Hamiltoniana e mu,sigma in funzione dei passi del SA)

  last++;                                                 //Mi dice che ho finito e posso campionare istogramma
  Energia=Energia_per_simulated_annealing();              //Stampo file con valori finali energia in funzione del numero di blocchi
  cout<<"STAMPA DEL FILE: 2.Energie_finali_Hamiltoniana(numero_blocchi).out"<<endl;

  Istogramma_Funz_onda_quadra();                          //Stampo file funzione d'onda quadra in funzione di x

  ConfFinal();                                              //Write final configuration

  return 0;
}



//FUNZIONI------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



void Input(){
  ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed;             //File da leggere
  
  cout << "Principio variazionale di Ritz in Meccanica Quantistica" << endl;
  cout << "Simulazione Montecarlo: particella 1D in potenziale definito" << endl;
  cout << "Potenziale: V(x) = x^4 - 5/2 x^2" << endl;
  cout << "Peso: modulo quadro della funzione d'onda" << endl;
  cout << "Il programma usa unità naturali" << endl;


//Read seed for random numbers
  int p1, p2;
  Primes.open("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

//Read input informations
  ReadInput.open("input.in");

  ReadInput >> restart;                                        //0=NON RIPARTO, 1=RIPARTO DA CONFIG PRECEDENTE

  if(restart) Seed.open("seed.out");
  else Seed.open("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  Seed.close();
    
  ReadInput >> mu;

  ReadInput >> sigma;
  
  ReadInput >> Delta_x;
  ReadInput >> Delta_mu;
  ReadInput >> Delta_sigma;
  
  ReadInput >> temp;

  ReadInput >> nblk;

  ReadInput >> nstep;

  cout << "Il programma esegue mosse con l'algoritmo di Metropolis con traslazioni uniformi" << endl << endl;
  cout <<"Elenco dei parametri dati in input:"<<endl;
  cout << "Mu iniziale = " << mu << endl;
  cout << "Sigma iniziale = " << sigma << endl; 
  cout << "Mossa (Delta_x) = " << Delta_x << endl;
  cout << "Mossa per mu (Delta_mu) = " << Delta_mu << endl;
  cout << "Mossa per sigma (Delta_sigma) = " << Delta_sigma << endl;
  cout << "Temperatura fittizia iniziale per Simulated Annealing = " << temp << endl;
  cout << "Numero di blocchi = " << nblk << endl;
  cout << "Numero di passi in ciascun blocco = " << nstep << endl << endl;
  ReadInput.close();                

  //Read initial configuration
  cout << "Il programma sta leggendo la configurazione iniziale" << endl << endl;
  if(restart){
    ReadConf.open("config.out");
  }
  else{
    ReadConf.open("config.in");
  }

  ReadConf >> x;                    //Leggo le configurazioni (ho messo particella che parte in  x=0)

  ReadConf.close();
  x_old = x;            

  return;
}


void Simulated_annealing (){                            //algoritmo di ottimizzazione per trovare mu e sigma del minimo, da usare poi in MonteCarlo

  ACCETTATO = 0;
  CONTATO = 0;           

  cout << "Uso il Simulated Annealing per trovare i valori di mu e sigma che minimizzano <H>" << endl;
  cout << "-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-" <<endl<<endl;
  ofstream Ham, Param;                                //Apro file di output per registrare passi del SA
  string aaa = "2.Hamiltoniana(passi).out";
  string bbb = "2.Parametri(passi).out";
  Ham.open(aaa, ios::app);
  Param.open(bbb, ios::app);

  while(temp>pow(10, -10)){

    CONTATO ++;

    mu_old = mu;
    sigma_old = sigma;
    Energia_Old = Energia_per_simulated_annealing();                  //Calcolo energia dei parametri mu e sigma correnti

    mu_new = mu + (rnd.Rannyu()-0.5)*Delta_mu;                              //Aggiorno i parametri casualmente
    sigma_new = sigma + (rnd.Rannyu()-0.5)*Delta_sigma;

    mu = mu_new;
    sigma = sigma_new;
    Energia_New = Energia_per_simulated_annealing();                    //VETTORE: Calcolo energia dei nuovi parametri
                                             
    probability = exp(-(1./temp)*(Energia_New[0] - Energia_Old[0]));                  //Nota: Se l'energia nuova è minore della vecchia probability > 1 (va sempre bene!)

    randd = rnd.Rannyu();
    if(randd > probability){                                                //Testo la probabilità e se NON va bene torno a mu e sigma vecchi
      mu = mu_old;
      sigma = sigma_old;
      Energia = Energia_Old;
    }
    else{
      ACCETTATO ++;
      Energia = Energia_New;
    }

    //Stampa temperatura, step del SA, Energia a quello step, incertezza statistica sull'energia a quello step
    Ham <<  setw(wd) << temp <<  setw(wd) << CONTATO  <<  setw(wd) << Energia[0] <<  setw(wd) << Energia[1] << endl;
    //Stampa temperatura, step del SA, parametri mu e sigma
    Param <<  setw(wd) << temp <<  setw(wd) << CONTATO <<  setw(wd) << mu <<  setw(wd) << sigma << endl;
    
    temp = temp*0.95;                                                         //Aggiorno temperatura
  }

  cout << endl <<"Valori ottimali di mu e sigma: " << mu << "     " << sigma << endl;
  cout << "Acceptance rate del Simulated Annealing: " << (double)ACCETTATO/(double)CONTATO << endl;
  cout << "Numero di temperature provate (ovvero passi del SA compiuti): " << CONTATO << endl;

  Ham.close();
  Param.close();

  cout<<"STAMPA DEL FILE: "<<aaa<<endl;
  cout<<"STAMPA DEL FILE: "<<bbb<<endl<<endl;

  return;
}

vector<double> Energia_per_simulated_annealing (){
  vector<double> Energia_Incertezza(2);

  for(int iblk=1; iblk <= nblk; iblk++){                            //Simulation
    Reset(iblk);                                                    //Reset block averages
        
    for(int istep=1; istep <= nstep; istep++){
      Move();
      Measure();
      Accumulate();                                                 //Update block averages
    }

    Averages(iblk);                                                 //Print results for current block
    Energia_Incertezza[0] = en_ham;
    Energia_Incertezza[1] = err_ham;
  }

  return Energia_Incertezza;
}


void Istogramma_Funz_onda_quadra(){
  int norm = 0;                                             //Calcolo della normalizzazione

  for(int istep=1; istep <= nstep; istep++){
    Move();
    if(x < abs(xmax)){
      Psi2[ int((x + xmax) / bin_size) ] += 1;                //Riempio l'istogramma
      norm ++;
    } 
  }
  norm = norm/nbins;
  
  ofstream Isto;                                           //Stampo su file di output
  string ccc = "2.Istogramma_funz_onda_quadra.out";
  Isto.open(ccc);
  for(int j=0; j< nbins; j++){
    Isto << -xmax + (bin_size * j) << " \t" << Psi2[j]/(2*xmax*norm) << endl;
  }
  Isto.close();
  cout<<"STAMPA DEL FILE: "<<ccc<<endl;

  return;
}

void Reset(int iblk){                                    //Reset block averages
  if(iblk == 1){
    glob_av = 0;
    glob_av2 = 0;
  }
  blk_av = 0;
  blk_norm = 0;
  attempted = 0;
  accepted = 0;
}

void Move(){
    //Mossa con Monte Carlo 
    //Old
    psi2_old = Prob_density(x);                           //Calcola la vecchia densità di probabilità

    //New
    x_new = x + Delta_x*(rnd.Rannyu() - 0.5);             //Propone nuova posizione
    psi2_new = Prob_density(x_new);                       //Calcola la nuova densità di probabilità

    //Metropolis test
    p = psi2_new/psi2_old;                                //Probabilità con cui fare confronto

    if(p >= rnd.Rannyu()){                                //Se trovo probabilità minore evolvo il sistema
      //Update
      x = x_new;
      accepted = accepted + 1.0;
    } 
    attempted = attempted + 1.0;
  return;
}

double Prob_density(double x){                           //Calcolo del modulo quadro della funzione d'onda
  //Non serve normalizzazione perchè Metropolis fa rapporto!
  return pow( exp( - pow( x - mu, 2) / (2.*sigma*sigma) ) + exp( - pow( (x + mu), 2) / (2.*sigma*sigma) ), 2 );     
}


void Measure(void) {                                    //FORMULA TROVATA ANALITICAMENTE FACENDO DERIVATA SECONDA DELLA PSI FRATTO LA PSI + somma potenziale)
    double En_potenziale = pow(x, 4) - 2.5 * pow(x, 2);
    double En_cinetica = exp( -(x-mu)*(x-mu)/(2.0* sigma*sigma) ) * ( (x-mu)*(x-mu)/(sigma * sigma) - 1 )/(sigma*sigma) + exp( -(x+mu)*(x+mu)/(2.0* sigma*sigma) ) * ( (x+mu)*(x+mu)/(sigma * sigma) - 1 )/(sigma*sigma);
    En_cinetica = - 0.5 * En_cinetica /( exp( - (x-mu)*(x-mu)/(2.0 * sigma*sigma) ) + exp( - (x+mu)*(x+mu)/(2.0 * sigma*sigma) ) );    //h_tagliato = 1, m = 1

    walker = (En_potenziale + En_cinetica);
    return;
}

void Accumulate(void){                                  //Update block averages

  blk_av = blk_av + walker;
  blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) {

  stima_ham = blk_av/blk_norm;                    //energia totale media del blocco
  glob_av += stima_ham;                           //di tutti i blocchi
  glob_av2 += stima_ham*stima_ham;

if (iblk == nblk){                             //for all temperatures, only the final block results are printed
            cout << "Temperature " << temp << endl;
            cout << "Acceptance rate " << (double) accepted/attempted << endl << endl;
            cout << "----------------------------" << endl << endl;
        }

  en_ham = glob_av/(double)iblk;                    //Energia
  err_ham = Error(glob_av,glob_av2,iblk);           //Incertezza 
    
  if(last){
    ofstream BestConf;
    BestConf.open("2.Energie_finali_Hamiltoniana(numero_blocchi).out", ios::app);
    BestConf << iblk << "  " << stima_ham << "  " << en_ham << "  " << err_ham << endl;
    BestConf.close();
  }
}

double Error(double sum, double sum2, int iblk){
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}


void ConfFinal(void){
  ofstream WriteConf, WriteSeed;

  cout << "Stampa la configurazione finale sul file config.out" << endl << endl;
  WriteConf.open("config.out");
  WriteConf << x << endl;
  WriteConf.close();

  rnd.SaveSeed();
}