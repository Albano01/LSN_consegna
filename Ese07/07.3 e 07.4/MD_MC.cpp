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
#include "MD_MC.h"

using namespace std;

int main(){

  Input();                                        //Inizializzazione
  int nconf = 1;

  for(int iblk=1; iblk <= nblk; iblk++){        
    Reset(iblk);                                  
  
    for(int istep=1; istep <= nstep; istep++){    //cycles over the steps in each block
      Move();                                     //Muovere molecole con algoritmo Verlet o Metropolis 
      Measure();                                  //Misurare proprietà
      Accumulate();                               //Updates the block averages
      if(istep%10 == 0){
        //ConfXYZ(nconf);                         //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
      }
    }
  
    Averages(iblk);       //Stampa risultati per blocco corrente
  }

  ConfFinal();              //Scrive le configurazioni finali in output

  return 0;
}

//Inizializza prendendo configurazioni iniziali e stampa i valori iniziali delle osservabili
void Input(void){
  ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed;

  //Stampa informazioni sul codice
  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "MD(NVE) / MC(NVT) simulation       " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  //RNG
  int p1, p2;
  Primes.open("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  //Legge informazioni in input
  ReadInput.open("input.in");

  ReadInput >> iNVET;       //0=Molecular Dynamics (NVE), 1=MonteCarlo (NVT)
  ReadInput >> restart;      //0=NON RIPARTO, 1=RIPARTO DA CONFIG PRECEDENTE

  if(restart)                                             //Caso riparto da config precedente
    Seed.open("seed.out");                                
  else 
    Seed.open("seed.in");                                 //Caso parto da nuova config
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);                              
  Seed.close();

  ReadInput >> temp;                                      //Legge temperatura iniziale
  beta = 1.0/temp;
  cout << "Input temperature = " << temp << endl;

  ReadInput >> eq_time;                                   //Legge tempo di equilibrazione (MESSO IO)
  cout << "Number of steps to reach equilibrium = " << eq_time << endl;

  ReadInput >> npart;                                     //Legge numero particelle
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;                                       //Legge la densità numerica
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;                                //Calcola il volume
  box = pow(vol,1.0/3.0);
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;

  larghezza_bin = (box/2.) / nbins;                            //calcola la lunghezza del bin per istogramma di g(r) (MESSO IO)

  ReadInput >> rcut;                                      //legge raggio di cutoff
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
  ReadInput >> delta;                                    

  ReadInput >> nblk;                                    

  ReadInput >> nstep;                                   

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

  //prepare arrays for measurements
  iv = 0;      //potential energy
  it = 1;      //temperature
  ik = 2;      //kinetic energy
  ie = 3;      //total energy
  ip = 4;      //pressure
  ig = 5;      //g(r)
  n_props = 6; //total number of observables

  //Read initial configuration and velocities
  cout << "Read initial configuration" << endl << endl;
  if(restart) {
    ReadConf.open("config.out");                
    ReadVelocity.open("velocity.out");         

    for (int i=0; i<npart; ++i)
      ReadVelocity >> vx[i] >> vy[i] >> vz[i];    
  
  } else {
    ReadConf.open("config.in");                 

    cout << "Prepare velocities with center of mass velocity equal to zero" << endl;
    double sumv[3] = {0.0, 0.0, 0.0};            

    for (int i=0; i<npart; ++i) {
      vx[i] = rnd.Gauss(0., sqrt(temp));          //Distribuzione della velocità secondo Maxwell-Boltzmann
      vy[i] = rnd.Gauss(0., sqrt(temp));
      vz[i] = rnd.Gauss(0., sqrt(temp));

      sumv[0] += vx[i];                  
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }

    for (int idim=0; idim<3; ++idim)
      sumv[idim] /= (double)npart;        
    
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){
      vx[i] = vx[i] - sumv[0];            
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];

      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];       
    }

    sumv2 /= (double)npart;             
    fs = sqrt(3 * temp / sumv2);          //fs = velocity scale factor
    cout << "velocity scale factor: " << fs << endl << endl;

    for (int i=0; i<npart; ++i){
      vx[i] *= fs;                        //Scales velocities 
      vy[i] *= fs;
      vz[i] *= fs;
    }
  }

  for (int i=0; i<npart; ++i) {
    ReadConf >> x[i] >> y[i] >> z[i];     //Leggo le configurazioni (unità di misura del lato)
    x[i] = Pbc( x[i] * box );              //Riscalo in unità di Lennard-Jones, CON ANCHE PERIODIC BOUNDARY CONDITION PER PRUDENZA
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  }

  ReadConf.close();


  for (int i=0; i<npart; ++i) {
    if(iNVET) {         //MC
      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];
    }
    else                //MD
    {
      xold[i] = Pbc(x[i] - vx[i] * delta);    
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }

  //TAILS CORRECTION: FATTO IO
  vtail = 8. * M_PI * rho * (1./(9.*pow(rcut,9)) - 1./(3.*pow(rcut,3)));
  ptail = 32.*M_PI*rho*rho*( 1./(9.* pow(rcut, 9)) - 1./(6.* pow(rcut, 3)) );

  //Equilibrazione
  cout << "Equilibrazione in corso!" << endl;
  for(int t=0; t<eq_time; t++){
    Move();
  }
  cout << "Il sistema ha raggiunto l'equilibrio." << endl<<endl;
  
  //Evaluate the properties of the initial configuration
  Measure();

  //Prints initial values for measured properties
  cout << "Initial potential energy = " << walker[iv]/(double)npart << endl;
  cout << "Initial temperature      = " << walker[it] << endl;
  //cout << "Initial kinetic energy   = " << walker[ik]/(double)npart << endl;
  //cout << "Initial total energy     = " << walker[ie]/(double)npart << endl;
  cout << "Initial pressure         = " << walker[ip] << endl << endl;

  return;
}


//Evoluzione verso nuova configurazione
void Move(){
  int o;
  double p, energy_old, energy_new;
  double xnew, ynew, znew;

  if(iNVET) {   //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Monte Carlo (NVT) move!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for(int i=0; i<npart; ++i){
      
      o = (int)(rnd.Rannyu() * npart);            //selects randomly a particle (for C++ syntax, 0 <= o <= npart-1)

      energy_old = Boltzmann(x[o],y[o],z[o],o);   //Old

      x[o] = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );    //new
      y[o] = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
      z[o] = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

      energy_new = Boltzmann(x[o],y[o],z[o],o);

      //Metropolis test
      p = exp(beta * (energy_old-energy_new));      //Peso di Boltzmann come probabilità per decidere se accettare la mossa
      if(p >= rnd.Rannyu()) {

        xold[o] = x[o];                             //Nuova nella vecchia se accetto
        yold[o] = y[o];
        zold[o] = z[o];

        accepted = accepted + 1.0;                  //rate accettanza +1
      } else {

        x[o] = xold[o];                             //Mantengo vecchia se rifiuto
        y[o] = yold[o];
        z[o] = zold[o];
      }
      attempted = attempted + 1.0;
    }
  } else    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Molecular Dynamics (NVE) move!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  {
    double fx[m_part], fy[m_part], fz[m_part];     

    for(int i=0; i<npart; ++i){     
      fx[i] = Force(i,0);
      fy[i] = Force(i,1);
      fz[i] = Force(i,2);
    }

    //Verlet
    for(int i=0; i<npart; ++i){ 

      xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
      ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
      znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

      vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
      vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
      vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

      xold[i] = x[i];       //vecchia con presente
      yold[i] = y[i];
      zold[i] = z[i];

      x[i] = xnew;          //presente con nuova
      y[i] = ynew;
      z[i] = znew;

      accepted = accepted + 1.0;          //accettanza DEVE essere 1
      attempted = attempted + 1.0;
    }
  }
  return;
}


//Calcola energia della particella ip nel potenziale di LJ con raggio di cutoff
double Boltzmann(double xx, double yy, double zz, int ip) {

  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i){            
    if(i != ip){

      dx = Pbc(xx - x[i]);                //distance ip-i in PBC
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut) {                             //entro il raggio di cutoff??
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);   
      }
    }
  }

  return 4.0*ene;
}


//Calcola forza in direzione idir (x=0, y=1, z=2) agente sulla ip-esima particella come -Grad_ip V(r)
double Force(int ip, int idir){

  double f=0.0;
  double dvec[3], dr;                     

  for (int i=0; i<npart; ++i){            
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );      //distance ip-i in PBCs
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){                                                //è dentro il raggio di cutoff??
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8));       
      }
    }
  }
 
  return f;
}


//Misura delle proprietà termodinamiche in una certa config
void Measure() {

  double v = 0.0, kin=0.0, p=0.0;     
  double dx, dy, dz, dr;

  for (int i=0; i<npart-1; ++i){     
    for (int j=i+1; j<npart; ++j){ 

      dx = Pbc(x[i] - x[j]);              
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = sqrt(dx*dx + dy*dy + dz*dz);

      //update istogramma di g(r)
      if ( dr < box/2. ){                      
        g[ int (dr/larghezza_bin) ] += 2;            
      }

      if(dr < rcut){                            //è nel raggio di cutoff??
        v += 1.0/pow(dr,12) - 1.0/pow(dr,6);    
        p += 1.0/pow(dr,12) - 0.5/pow(dr,6);   
      }
    }          
  }

  for (int i=0; i<npart; ++i)
    kin += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);                   //measures kinetic energy

  walker[iv] = 4.0 * v;                                                       //current value of potential energy
  //walker[ik] = kin;                                                         //current value of kinetic energy
  walker[it] = (2.0 / 3.0) * kin/(double)npart;                               //current value of temperature
  //walker[ie] = 4.0 * v + kin;                                               //current value of total energy
  walker[ip] = rho * walker[it] + (16./vol) * p;                              //current value of pressure

  return;
}


//Reset block averages
void Reset(int iblk) {
   
  if(iblk == 1) {                      
    for(int i=0; i<n_props; ++i) {

      glob_av[i] = 0;
      glob_av2[i] = 0;
    }

    for(int i=0; i<nbins; ++i){
      g[i] =0;
      glob_av_g[i]=0;
      glob_av2_g[i]=0;
    }
  }

  for(int i=0; i<n_props; ++i) {
    blk_av[i] = 0;                   
  }
  
  for(int i=0; i<nbins; i++){
    g[i] =0;                     
  }
   
  blk_norm = 0;
  attempted = 0;
  accepted = 0;
}


//Update the block averages
void Accumulate(void) {

  for(int i=0; i<n_props; ++i){
    blk_av[i] = blk_av[i] + walker[i];     
  }

  blk_norm = blk_norm + 1.0;          
}


//Medie e incertezza del blocco corrente!
void Averages(int iblk) {
    
  ofstream Epot, Ekin, Etot, Temp, Press, RadDis;
  const int wd=12; 

  cout << "Block number " << iblk << endl;
  cout << "Acceptance rate " << accepted/attempted << endl << endl;

  if(iNVET){
    Epot.open("../4.MC.epot.out",ios::app);
    Press.open("../4.MC.press.out",ios::app);
    RadDis.open("../4.MC.dis_rad.out",ios::app);
  } else{
    Epot.open("../4.MD.epot.out",ios::app);
    Press.open("../4.MD.press.out",ios::app);
    RadDis.open("../4.MD.dis_rad.out",ios::app);
  }
    
  /*
  NON SERVONO PER QUESTA ESERCITAZIONE
  stima_kin = blk_av[ik]/blk_norm/(double)npart;    //kinetic energy
  glob_av[ik] += stima_kin;                       
  glob_av2[ik] += stima_kin*stima_kin;
  err_kin=Error(glob_av[ik],glob_av2[ik],iblk);    

  stima_etot = blk_av[ie]/blk_norm/(double)npart;   //total energy
  glob_av[ie] += stima_etot;                    
  glob_av2[ie] += stima_etot*stima_etot;
  err_etot=Error(glob_av[ie],glob_av2[ie],iblk);    
  */

  stima_temp = blk_av[it]/blk_norm;                 //temperature
  glob_av[it] += stima_temp;                        
  glob_av2[it] += stima_temp*stima_temp;
  err_temp=Error(glob_av[it],glob_av2[it],iblk);    

  stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail;  //Potential energy + tail correction
  glob_av[iv] += stima_pot;                              
  glob_av2[iv] += stima_pot*stima_pot;
  err_pot=Error(glob_av[iv],glob_av2[iv],iblk);     

  stima_press = blk_av[ip]/blk_norm + ptail;        //Pressure + tail correction
  glob_av[ip] += stima_press;                      
  glob_av2[ip] += stima_press*stima_press;
  err_press=Error(glob_av[ip],glob_av2[ip],iblk);   
  

  //potential energy
  Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
  //pressure
  Press << setw(wd) << iblk <<  setw(wd) << stima_press << setw(wd) << glob_av[ip]/(double)iblk << setw(wd) << err_press << endl;

  //g(r): ho nbins per blocco, stampo su file di output BLOCCO, limite inferiore del RAGGIO, VALORE a quella distanza per quel blocco, MEDIA a quella distanza per i blocchi fatti, INCERTEZZA
  for(int i=0; i<nbins; i++){
    rI = i * larghezza_bin;              //limite inferiore dell'i-esimo bin
    rII = rI + larghezza_bin;            //limite superiore dell'i-esimo bin

    
    stima_gi = (g[i]/blk_norm) *3./(4.*M_PI*rho*npart * (pow(rII, 3) - pow(rI, 3)));

    glob_av_g[i] += stima_gi;                                  
    glob_av2_g[i] += stima_gi*stima_gi;
    err_gi = Error(glob_av_g[i],glob_av2_g[i],iblk);     //incertezza statistica
  
    //stampo quanto scritto sopra
    RadDis << iblk <<  setw(wd) << rI << setw(wd) << stima_gi << setw(wd) << glob_av_g[i]/(double)iblk << setw(wd) << err_gi << endl;
  }
  
  cout << "----------------------------" << endl << endl;

  Epot.close();
  Press.close();
  RadDis.close();
}


//congif finali
void ConfFinal(void){
  ofstream WriteConf, WriteVelocity, WriteSeed;

  cout << "Print final configuration to file config.out" << endl << endl;
  WriteConf.open("config.out");
  WriteVelocity.open("velocity.out");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;  
    WriteVelocity << vx[i] << "   " <<  vy[i] << "   " << vz[i] << endl;       
  }

  WriteConf.close();
  WriteVelocity.close();

  rnd.SaveSeed();          
}


//writes the nconf configuration in .xyz format 
void ConfXYZ(int nconf){ 
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

//algorithm for periodic boundary conditions with side L=box
double Pbc(double r){
  return r - box * rint(r/box);   
}

//Dev st della media di iblk blocchi
double Error(double sum, double sum2, int iblk){
  return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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