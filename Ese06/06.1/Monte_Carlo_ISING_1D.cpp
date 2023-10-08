/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  Input(); //Inizializzazione
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro;                                           //Se metro=1 usa algoritmo di Metropolis, altrimenti di Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration (MODIFICATO IO)
  for (int i=0; i<nspin; ++i)
  {
    s[i] = 1;
  }
  //CONVIENE PARTIRE TUTTI UP (O TUTTI DOWN) COSì EQUILIBRAZIONE ALLE BASSE T è PIù VELOCE!!
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, Boltz, randd;
  //double energy_up, energy_down; (NON USATO, potevo usarli in Gibbs richiamando Boltzmann e facendo delta, ma ho messo formula compatta calcolata a mano)

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1){ //Metropolis (IMPLEMENTATO DA ME)
      energy_old = Boltzmann(s[o], o);              //Salvo l'energia precedente
      energy_new = Boltzmann(-s[o], o);         //Salvo la nuova energia

      Boltz = pow(M_E, (-beta*(energy_new - energy_old)));        //Calcolo l'esponenziale della probabilità secondo Boltzmann
      if (Boltz > 1){
        p=1;
      }
      else{
        p=Boltz;
      }

      randd=rnd.Rannyu();         //Estraggo numero casuale tra 0 ed 1 da confrontare con la probabilità
      if (p>randd){               //Se la probabilità è maggiore evolvo il sistema
        s[o]=-s[o];
        accepted ++;                //Conto che è stato accettato  (accepted è posto giustamente a zero nel Reset)
      }
      attempted ++;                //Conto che è stato tentato   (attempted è posto giustamente a zero nel Reset)
    }
    else{ //Gibbs sampling (IMPLEMENTATO DA ME)

      p = 1./(1. + exp(-2*beta*J*(s[Pbc(o-1)]+s[Pbc(o+1)])-2*beta*h));     //Calcolo probabilità per porre spin nella posizione o in +1
    
      randd=rnd.Rannyu();         //Estraggo numero casuale tra 0 ed 1 da confrontare con la probabilità
      if (p>randd){
        s[o]=1;
      }
      else{
        s[o]=-1;
      }
      accepted ++;                //Conto che è stato accettato  (accepted è posto giustamente a zero nel Reset)
      attempted ++;                //Conto che è stato tentato   (attempted è posto giustamente a zero nel Reset)
      //In Gibbs l'accettazione è sempre 1 (quindi inutile, ma è già scritta per Metropolis)
    } 
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure(){
  //Aggiungo io c,m,x. u già implementata. bin?? h è il campo esterno, beta è 1/T
  //int bin;  (NON USATO)
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i){
    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
    m += s[i];
  }
  walker[iu] = u;
  walker[ic] = pow(u, 2);   //Vedi fine implementazione in Averages
  walker[im] = m;
  walker[ix] = beta*pow (m,2);
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    if(h==0){                                             //Se campo esterno è spento calcolo energia interna, capacità termica e suscettività magnetica
      //Già implementato ENERGIA INTERNA
      if (metro == 1){                                    //Se Metropolis metto la M nel nome del file
        Ene.open("../M.output.ene.0",ios::app);
      }
      else{
        Ene.open("../G.output.ene.0",ios::app);              //Se Gibbs metto la G nel nome del file
      }
      stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
      glob_av[iu]  += stima_u;
      glob_av2[iu] += stima_u*stima_u;
      err_u=Error(glob_av[iu],glob_av2[iu],iblk);
      Ene << setw(wd) << iblk <<  setw(wd) << temp << setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
      Ene.close();

      //Implemento io CAPACITà TERMICA (vedi formula nella consegna, iniziata in measures e finita qui)
      if (metro == 1){                                      //Se Metropolis metto la M nel nome del file
        Heat.open("../M.output.heat.0",ios::app);
      }
      else{
        Heat.open("../G.output.heat.0",ios::app);              //Se Gibbs metto la G nel nome del file
      }
      stima_c = pow(beta, 2)*(blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm, 2))/(double)nspin; 
      glob_av[ic]  += stima_c;
      glob_av2[ic] += stima_c*stima_c;
      err_c=Error(glob_av[ic],glob_av2[ic],iblk);
      Heat << setw(wd) << iblk <<  setw(wd) << temp << setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
      Heat.close();


      //Implemento io SUSCETTIVITà MAGNETICA
      if (metro == 1){                                    //Se Metropolis metto la M nel nome del file
        Chi.open("../M.output.chi.0",ios::app);
      }
      else{
        Chi.open("../G.output.chi.0",ios::app);              //Se Gibbs metto la G nel nome del file
      }
      stima_x = blk_av[ix]/blk_norm/(double)nspin;        
      glob_av[ix]  += stima_x;
      glob_av2[ix] += stima_x*stima_x;
      err_x=Error(glob_av[ix],glob_av2[ix],iblk);
      Chi << setw(wd) << iblk << setw(wd) << temp << setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
      Chi.close();
      cout << "-----------------------------------------------------------------------------------------------------" << endl << endl;
    } 

    else{                                                     //Se il campo esterno è acceso calcolo magnetizzazione
      //Implemento io MAGNETIZZAZIONE CON CAMPO ACCESO
      if (metro == 1){                                    //Se Metropolis metto la M nel nome del file
        Mag.open("M.output.mag.0,02",ios::app);
      }
      else{
        Mag.open("G.output.mag.0,02",ios::app);              //Se Gibbs metto la G nel nome del file
      }
      stima_m = blk_av[im]/blk_norm/(double)nspin;
      glob_av[im]  += stima_m;
      glob_av2[im] += stima_m*stima_m;
      err_m=Error(glob_av[im],glob_av2[im],iblk);
      Mag << setw(wd) << iblk << setw(wd) << temp << setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
      Mag.close();
      cout << "-----------------------------------------------------------------------------------------------------" << endl << endl;
    }

      
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
