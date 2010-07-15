/*
 
Copyright (C) 2007 Jindrich Kolorenc

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 
*/


#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>  // formating output stream
#include <cmath>    // pow, abs
#include <limits>   // to figure out precision of double
#include <cstdlib>
#include <cstring>
#include "indexx.h" // Quicksort from Numerical Recipes

//#define PRINT_VECS 1        // to print k-vectors of basis functions
//#define PRINT_SHELL_DETAILS // to print extra summary for each shell 
//#define PBC_ONLY 1          // to compile Gamma-point only version


using namespace std;

void write_basis(string oname, int twist, double k0,
		 int* kx, int* ky, int* kz, int* idx,
		 int basissize, int outprec);
void write_slater(string oname, int twist, int* Ne);
void write_bcs(string oname, int twist, int* Ne, int nmo);
void write_orb(string oname, int twist, int basissize, int nmo,
	       int first_k2);
void write_sys(string oname, int twist, int tw[4][3], int* Ne, double L,
	       int outprec);
void write_jast2_simple(string oname, double L);
void write_jast2(string oname, int twist, double L, double k0,
		 int* kx, int* ky, int* kz, int* idx, int* k2,
		 int* shells, int nshells, int basissize, int outprec);


int main(int argc, char ** argv) {
  // {{{ .

  const int ikmax=8;    // upper limit for planewave search
  const double pi=3.1415926535897932385;

  // four twists that lead to a real wavefunction
  char pbc[2]={ 'A', 'P' };
  int tw[4][3]={ {1,1,1}, {1,1,0}, {1,0,0}, {0,0,0} };
  int tw_weight[4]={ 1, 3, 3, 1 };
  int ikmax3dat[4]={ (2*ikmax+1)*(2*ikmax+1)*ikmax+(2*ikmax+1)*ikmax+ikmax+1,
		     (2*ikmax+1)*(2*ikmax+1)*ikmax,
		     (2*ikmax+1)*2*ikmax*ikmax,
		     2*ikmax*2*ikmax*ikmax };

  // find out how many decimal places has double and set precision
  // for printing k-vectors and cell size
  typedef numeric_limits< double > dl;
  int outprec=dl::digits10;  
  //cout << "double has " << outprec << " digits" << endl;
  
  // base for output files is given on the command line
  string oname;
  bool oflag=false;
  bool bcs_switch=false;
  for(int i=1; i< argc; i++) {
    if(!strcmp(argv[i], "-o") && argc > i+1) {
      oname=argv[++i];
      oflag=true;
    }
    else if(!strcmp(argv[i], "-bcs")) {
      bcs_switch=true;
    }
    else if(argv[i][0]=='-') {
      cerr << "Unknown option: " << argv[i] << endl;
    }
  }

  if ( !oflag ) {
    cout << "Usage: heg2qmc -o <output base> [-bcs]" << endl;
    cout << endl;
    cout << "Takes N_up, N_down and rs from standard input and ";
    cout << "writes" << endl;
    cout << "   <output base>.i.sys" << endl;
    cout << "   <output base>.i.slater" << endl;
    cout << "   <output base>.i.basis" << endl;
    cout << "   <output base>.i.orb" << endl;
    cout << "   <output base>.i.jast2   (Jastrow with long-range part)"
	 << endl;  
    cout << "   <output base>.jast2     (short-range only Jastrow)" 
	 << endl;
    cout << "   <output base>.centers" << endl;
    cout << "where i indexes k-points (twists of boundary conditions) and ";
    cout << "runs from 0 to 3;" << endl;
    cout << "0 is the Gamma-point." << endl;
    cout << endl;
    cout << "If -bcs option is given, extra input is taken (number ";
    cout << "of virtual orbitals), " << endl;
    cout << "and extra files are written: <output base>.i.bcs). These "
	 << "are in Pfaffian format, which is now obsolete for this "
	 << "application." << endl;
    cout << endl;
    return 1;
  }   

  // figure out the parameters we need (read from standard input)
  int Ne[2];
  int Nup, Ndown;
  cout << "Number of spin-up electrons" << endl;
  cin >> Nup;
  Ne[0]=Nup;
  cout << "Number of spin-down electrons" << endl;
  cin >> Ndown;
  Ne[1]=Ndown;
  int nmo;           // number of molecular orbitals
  if ( Ne[0] > Ne[1] ) { nmo=Ne[0]; } else { nmo=Ne[1]; }
  if ( bcs_switch ) {
    int nvirt;
    cout << "Number of virtual orbitals" << endl;
    cin >> nvirt;
    nmo+=nvirt;
  }

  cout << "rs (in atomic units, i.e., in bohrs)" << endl;
  double rs;
  cin >> rs;
  cout << endl;
  double L=rs*pow(4.0/3*(Ne[0]+Ne[1])*pi,1.0/3);
#ifdef PBC_ONLY
  double k0=2*pi/L;  // for the PBC only version
#else
  double k0=pi/L;    // for PBC and ABC mixtures
#endif
  double Ekin=0;     // Kinetic energy for constructed slater wavefunction

  // now we construct all k-vectors in a cube and sort them with
  // respect to their magnitude; then we occuppy as many of them as
  // needed
  //
  // there are two branches of the code: periodic boundary conditions
  // (PBC) only and more general case with either PBC or antiperiodic
  // boundary conditions (ABC) along each direction; the PBC part could
  // be dropped, but it is not bad to have a check

#ifdef PBC_ONLY
  // PBC in all directions
  //
  // I need only k from k,-k pair...
  // BTW, at last I understand what is that strange index game
  // Lucas plays when setting up Ewald interaction
  //
  //const int ikmax3=(2*ikmax+1)*(2*ikmax+1)*ikmax+(2*ikmax+1)*ikmax+ikmax+1;
  const int ikmax3=ikmax3dat[0];
  int kx[ikmax3], ky[ikmax3], kz[ikmax3], k2[ikmax3];
  int idx[ikmax3], shells[ikmax3];

  int jkmin,kkmin;
  int seq=0;
  for( int i=0; i <= ikmax; i++ ) {    // half the cube
    jkmin=-ikmax;
    if ( i==0 ) jkmin=0;              // half the kx=0 plane
    for( int j=jkmin; j <= ikmax; j++ ) {
      kkmin=-ikmax;
      if ( i==0 && j==0 ) kkmin=0;    // half the kx=0,ky=0 line 
      for( int k=kkmin; k <= ikmax; k++ ) {
	kx[seq]=i;
	ky[seq]=j;
	kz[seq]=k;
	k2[seq]=i*i+j*j+k*k;
	seq++;
      }
    }
  }
  // twist is not used in PBC only version, but needs to be declared for
  // write_basis, write_slater etc.
  int twist;
  //
#else // PBC_ONLY
  // choice of PBC or ABC in each direction

  // loop over all defined twists
  for (int twist=0; twist<4; twist++) {
    cout << "Twist ";
    for (int i=0; i < 3; i++) {
      cout << pbc[tw[twist][i]];
    }
    cout << " (weight " << tw_weight[twist] << "), ";
#ifdef PRINT_SHELL_DETAILS
    cout << endl;
#endif // PRINT_SHELL_DETAILS
  
  double Ekin_tw=0;

  int ikmax3=ikmax3dat[twist];
  int kx[ikmax3], ky[ikmax3], kz[ikmax3], k2[ikmax3];
  int idx[ikmax3], shells[ikmax3];
  
  int ikmin, jkmin, kkmin;
  int seq=0;
  if ( tw[twist][0] ) { ikmin=0; } else { ikmin=1; }
  for( int i=ikmin; i <= 2*ikmax; i+=2 ) {
    if ( tw[twist][1] ) { jkmin=-2*ikmax; } else { jkmin=-2*ikmax+1; }
    if ( i==0 ) { jkmin+=2*ikmax; }
    for( int j=jkmin; j <= 2*ikmax; j+=2 ) {
      if ( tw[twist][2] ) { kkmin=-2*ikmax; } else { kkmin=-2*ikmax+1; }
      if ( i==0 && j==0 ) { kkmin+=2*ikmax; }
      for( int k=kkmin; k <= 2*ikmax; k+=2 ) {
	kx[seq]=i;
	ky[seq]=j;
	kz[seq]=k;
	k2[seq]=i*i+j*j+k*k;
	seq++;
      }
    }
  }  

#endif // PBC_ONLY

  if ( seq != ikmax3 ) {
    cerr << "Something's wrong: got incorrect number of k-vectors." << endl;
    cerr << seq << " instead of " << ikmax3 << endl;
    exit(1);
  }

  // order according to length of vectors (i.e., kinetic energy)
  // LKW note..this could be done better using objects and the STL sort,
  //but I'll leave it for now.
  // JK note..be my guest
  indexx(k2,idx,ikmax3);

  // find closed shells
  int shell=0;     // counter of closed shells
  int inshell=0;   // counter of states in between consecutive shells
  for ( seq=0; seq < ikmax3; seq++ ) {

    // we do not have complete shells for vectors from corners of the
    // "test" cube, throw these vectors away
#ifdef PBC_ONLY
    if ( k2[idx[seq]] > ikmax*ikmax ) break;
#else
    if ( k2[idx[seq]] > (2*ikmax-1)*(2*ikmax-1) ) break;
#endif
    
#ifdef PRINT_VECS
    cout << setw(4) << seq 
	 << setw(4) << kx[idx[seq]] 
	 << setw(4) << ky[idx[seq]] 
	 << setw(4) << kz[idx[seq]] 
	 << setw(5) << k2[idx[seq]] << endl;
#endif

    inshell++;

    if ( seq+1 < ikmax3 && k2[idx[seq+1]] > k2[idx[seq]] ) {
      // every k-vector we found has opposite brother, except (0,0,0),
      // which is present only in the case of pure PBC
      if ( shell > 0 ) { 
	shells[shell]=2*inshell+shells[shell-1];
      } else {
	if ( k2[idx[0]] == 0 ) { shells[0]=1; } else { shells[0]=2*inshell; }
      }
#ifdef PRINT_VECS
      cout << "----- " << setw(2) << shell+1 << " ("
	   << setw(4) << shells[shell] << ") -----" << endl;
#endif // PRINT_VECS
#ifdef PRINT_SHELL_DETAILS
      cout << "shell: states ";
      if ( k2[idx[seq]] == 0 ) { cout << 1; } else { cout << 2*inshell; }
      cout << "; |k|=" << setprecision(outprec) << k0*sqrt(1.0*k2[idx[seq]])
	   << endl;
#endif // PRINT_SHELL_DETAILS
      inshell=0;
      shell++;
    }  
  }  // for ( seq=0; seq < ikmax3; seq++ )
  int nshells=shell;

  // print out closed shells
  cout << "first " << nshells << " closed shells:" << endl;
  const int lcut=12;
  if ( nshells < lcut ) {
    for (int shell=0; shell < nshells; shell++) {
      cout << setw(6) << shells[shell];
    }
    cout << endl;
  } else {
    shell=0;
    for (int i=0; i < nshells/lcut+1; i++) {
      int jmax=nshells-shell;
      if ( jmax > lcut ) jmax=lcut;
      for(int j=0; j < jmax; j++) {
	cout << setw(6) << shells[shell];
	shell++;
      }
      cout << endl;
    }
  }
  cout << endl;
  
  // test whether user requested more states than we have found
  if ( Ne[0] > shells[nshells-1] || Ne[1] > shells[nshells-1] ) {
    cerr << "Maximum number of electrons of one kind is currently "
	 << shells[nshells-1] << endl;
    cerr << "To go higher, change hardcoded limit of basis size (ikmax)."
	 << endl;
    exit(1);
  }

  int basissize=(nmo+1)/2;
  // if we have (0,0,0) k-vector in the game, then the second
  // basis function, sin(0)=0, is not usable and hence we might
  // need to add one more function at the top
  if ( k2[idx[0]]==0 && 2*basissize==nmo ) basissize+=1;

  // write basis set file
  write_basis(oname, twist, k0, kx, ky, kz, idx, basissize,  outprec);

  // write slater determinant
  write_slater(oname, twist, Ne);
  if ( bcs_switch ) write_bcs(oname, twist, Ne, nmo);
  
  // write orb file
  write_orb(oname, twist, basissize, nmo, k2[idx[0]]);

  // write sys file
  write_sys(oname, twist, tw, Ne, L, outprec);

  // write simple two-body Jastrow
  write_jast2(oname, twist, L, k0, kx, ky, kz, idx, k2,
	      shells, nshells, basissize, outprec);  

  // calculate the (kinetic) energy of non-interacting particles
  // to have an idea of this part of finite-size errors
  for (int i=0; i<2; i++) {
    int shift=0;
    if ( k2[0]==0 ) shift=1;
    for (int j=0; j<Ne[i]-shift; j++) {
      double dEkin=k2[idx[j/2+shift]]*k0*k0/2;
#ifndef PBC_ONLY
      dEkin*=tw_weight[twist];
      Ekin_tw+=dEkin;
#endif
      Ekin+=dEkin;
    }
  }
  
#ifndef PBC_ONLY
  // inetic energy for a given twist
  cout << "Approximate kinetic energy (hartree): " 
       << Ekin_tw/tw_weight[twist] << " ("
       << Ekin_tw/(Ne[0]+Ne[1])/tw_weight[twist] << " per particle)" << endl;
  cout << endl << "---" << endl << endl;
  Ekin_tw=0;

  } // for (int twist=1; twist<4; twist++)
#endif

  // write the file with definition of a single center
  string centername=oname+".centers";
  ofstream centerout(centername.c_str());
  centerout << "1" << endl; 
  centerout << "origin 0.0 0.0 0.0" << endl;
  centerout.close();

  // write simple two-body Jastrow, this one is k-point independent
  write_jast2_simple(oname, L);

  // compare our approximate kinetic energy with exact one
#ifndef PBC_ONLY
  int tot_weight=0;
  for( int i=0; i < 4; i++) tot_weight+=tw_weight[i];
  Ekin/=tot_weight;
#endif
  double Ekin_ex=0;
  for (int i=0; i<2; i++) {
    double kF=pow(6*pi*pi*Ne[i],1.0/3)/L;
    Ekin_ex+=L*L*L*pow(kF,5.0)/(20*pi*pi);
  }
  cout << "Kinetic energy (hartree):" << endl;
  cout << "  approximate: " << Ekin    << "  ("
       << Ekin/(Ne[0]+Ne[1]) << " per particle)" << endl;
  cout << "  exact      : " << Ekin_ex << "  ("
       << Ekin_ex/(Ne[0]+Ne[1]) << " per particle)" << endl;
  cout << "  error      : " << 100*abs((Ekin-Ekin_ex)/Ekin) << "%"
       << endl << endl;

  cout << "---" << endl << endl;
  cout << "BTW, magnification factor in *.slater might need to be lowered for "
       << "large" << endl << "number of particles (cca 200+ per spin channel)."
       << endl << endl;

  return 0;

  // }}}
}


void write_basis(string oname, int twist, double k0,
		 int* kx, int* ky, int* kz, int* idx,
		 int basissize, int outprec) {
  // {{{ .

#ifdef PBC_ONLY
  string basisname=oname+".basis";
#else
  stringstream str_twist;
  str_twist << twist;
  string basisname=oname+"."+str_twist.str()+".basis";
#endif
  ofstream basisout(basisname.c_str());
  basisout << "  BASIS {" << endl;
  basisout << "   origin" << endl;
  basisout << "   PLANEWAVE" << endl; 
  basisout << "   GVECTOR {" << endl;
  
  for (int i=0; i < basissize; i++) {
    basisout << "   ";
    basisout << setprecision(outprec) << setw(outprec+6) << k0*kx[idx[i]]
             << setprecision(outprec) << setw(outprec+6) << k0*ky[idx[i]]
             << setprecision(outprec) << setw(outprec+6) << k0*kz[idx[i]]
	     << endl;
  }

  basisout << "   }" << endl;
  basisout << "  }" << endl;  
  basisout.close();

  // }}}
}

void write_slater(string oname, int twist, int* Ne) {
  // {{{ .

  int nmo_slater;
  if ( Ne[0] > Ne[1] ) { nmo_slater=Ne[0]; } else { nmo_slater=Ne[1]; }

#ifdef PBC_ONLY
  string slatername=oname+".slater";
#else
  stringstream str_twist;
  str_twist << twist;
  string slatername=oname+"."+str_twist.str()+".slater";
#endif
  ofstream slaterout(slatername.c_str());
  slaterout << " SLATER" << endl;
  slaterout << " ORBITALS {" << endl;
  slaterout << "  BASFUNC_MO" << endl;
  slaterout << "  MAGNIFY 0.5" << endl;
  slaterout << "  INCLUDE " << oname;
#ifndef PBC_ONLY
  slaterout << "." << twist;
#endif
  slaterout << ".basis" << endl;
  slaterout << "  CENTERS { READ " << oname << ".centers }" << endl;
  slaterout << "  NMO " << nmo_slater << endl;
  slaterout << "  ORBFILE " << oname;
#ifndef PBC_ONLY
  slaterout << "." << twist;
#endif
  slaterout << ".orb" << endl;
  slaterout << " }" << endl;
  slaterout << " DETWT { 1.0 }" << endl;
  slaterout << " STATES {" << endl;
  
  for(int k=0; k < 2; k++) {
    if ( k==0 ) {
      slaterout << "  # spin-ups" << endl;
    } else {
      slaterout << "  # spin-downs" << endl;
    }
    if ( Ne[k] < 10 ) {
      slaterout << "  " ;
      for (int orb=0; orb < Ne[k]; orb++) {
	slaterout << setw(6) << orb+1;
      }
      slaterout << endl;
    } else {
      int orb=0;
      for (int i=0; i < Ne[k]/10+1; i++) {
	int jmax=Ne[k]-orb;
	if ( jmax > 10 ) jmax=10;
	if ( jmax > 0 ) slaterout << "  " ;
	for(int j=0; j < jmax; j++) {
	  orb++;
	  slaterout << setw(6) << orb;
	}
	if ( jmax > 0 ) slaterout << endl;
      }
    }
  }

  slaterout << " }" << endl;
  slaterout.close();

  // }}}
}

void write_bcs(string oname, int twist, int* Ne, int nmo) {
  // {{{ .

  int nmo_slater;
  if ( Ne[0] > Ne[1] ) { nmo_slater=Ne[0]; } else { nmo_slater=Ne[1]; }
  
#ifdef PBC_ONLY
  string bcsname=oname+".bcs";
#else
  stringstream str_twist;
  str_twist << twist;
  string bcsname=oname+"."+str_twist.str()+".bcs";
#endif
  ofstream bcsout(bcsname.c_str());
  
  bcsout << " PFAFFIAN" << endl;
  bcsout << " ORBITALS {" << endl;
  bcsout << "  BASFUNC_MO" << endl;
  bcsout << "  MAGNIFY 0.5" << endl;
  bcsout << "  INCLUDE " << oname;
#ifndef PBC_ONLY
  bcsout << "." << twist;
#endif
  bcsout << ".basis" << endl;
  bcsout << "  CENTERS { READ " << oname << ".centers }" << endl;
  bcsout << "  NMO " << nmo << endl;
  bcsout << "  ORBFILE " << oname;
#ifndef PBC_ONLY
  bcsout << "." << twist;
#endif
  bcsout << ".orb" << endl;
  bcsout << " }" << endl;
  // following NPAIRS is for unpolarized case only, will investigate later
  bcsout << " NPAIRS { " << Ne[0] << " " << Ne[1] << " 0 }" << endl;
  bcsout << " PFWT { 1 }" << endl;
  bcsout << " PAIRING_ORBITAL {" << endl;
  bcsout << "  ORBITALS_IN_PAIRING {" << endl;
  if ( nmo < 10 ) {
    bcsout << "  " ;
    for (int orb=0; orb < nmo; orb++) {
      bcsout << setw(6) << orb+1;
    }
    bcsout << endl;
  } else {
    int orb=0;
    for (int i=0; i < nmo/10+1; i++) {
      int jmax=nmo-orb;
      if ( jmax > 10 ) jmax=10;
      if ( jmax > 0 ) bcsout << "  " ;
      for(int j=0; j < jmax; j++) {
	orb++;
	bcsout << setw(6) << orb;
      }
      if ( jmax > 0 ) bcsout << endl;
    }
  }
  bcsout << "  }" << endl;
  bcsout << "  OPTIMIZE_PF {" << endl;
  bcsout << "   SINGLET_DIAG" << endl;
  bcsout << "  }" << endl;
  bcsout << "  SINGLET_COEF {" << endl;
  for ( int i=nmo; i > 0; i--) {
    if ( nmo-i < Ne[0] ) { bcsout << "  1 "; } else { bcsout << "  0 "; }
    for ( int j=1; j < i; j++) { bcsout << "0 "; }
    bcsout << endl;
  }
  bcsout << "  }" << endl;
  bcsout << " }" << endl;

  // }}}
}

void write_orb(string oname, int twist, int basissize, int nmo,
	       int first_k2) {
  // {{{ .
#ifdef PBC_ONLY

  string orbname=oname+".orb";
  ofstream orbout(orbname.c_str());
  orbout << 1 << "  " << 1 << endl;
  // second basis function is invalid in this case, sin(0)=0
  for( int i=3; i < nmo+2; i++) {
    orbout << i << "  " << 1 << endl;
  }
  orbout.close();
  
  // let's write also orb file for STANDARD_MO for benchmarking purposes
  orbname=oname+".STANDARD_MO.orb";
  orbout.open(orbname.c_str());
  int coefno=0;
  for ( int i=1; i < nmo+1; i++ ) {
    for ( int j=1; j <= 2*basissize; j++ ) {
      coefno++;
      orbout << i << " " << j << " " << 1 << " " << coefno << endl;
    }
  }
  orbout << "COEFFICIENTS" << endl;
  for ( int i=1; i < nmo+2; i++ ) {
    for ( int j=1; j <= 2*basissize; j++ ) {
      if ( i == 2 ) i++;
      if ( i==j ) {
	orbout << "1.0 ";
      } else {
	orbout << "0.0 ";
      }
    }
    orbout << endl;
  }
  orbout.close();

#else // PBC_ONLY

  stringstream str_twist;
  str_twist << twist;
  string orbname=oname+"."+str_twist.str()+".orb";
  ofstream orbout(orbname.c_str());
  int istart=1;
  int istop=nmo+1;
  // if we have (0,0,0) k-vector in the game, then the second
  // basis function, sin(0)=0, is not a valid orbital and has
  // to be skipped
  //if ( k2[idx[0]]==0 ) {
  if ( first_k2==0 ) {
    orbout << 1 << "  " << 1 << endl;
    istart=3;
    istop=nmo+2;
  }
  for( int i=istart; i < istop; i++) {
    orbout << i << "  " << 1 << endl;
  }
  orbout.close();
  
#endif // PBC_ONLY

  // }}}
}

void write_sys(string oname, int twist, int tw[4][3], int* Ne, double L,
	       int outprec) {
  // {{{ .
#ifdef PBC_ONLY
  string sysname=oname+".sys";
#else
  stringstream str_twist;
  str_twist << twist;
  string sysname=oname+"."+str_twist.str()+".sys";
#endif
  ofstream sysout(sysname.c_str());
  sysout << "SYSTEM { HEG" << endl;
  sysout << " NSPIN { " << Ne[0] << " " << Ne[1] << " }" << endl;
  /*
  sysout << " LATTICEVEC {" << endl;
  sysout << setprecision(outprec) << setw(outprec+6) << L
	 << setprecision(outprec) << setw(outprec+6) << 0.0
	 << setprecision(outprec) << setw(outprec+6) << 0.0 << endl;
  sysout << setprecision(outprec) << setw(outprec+6) << 0.0
	 << setprecision(outprec) << setw(outprec+6) << L
	 << setprecision(outprec) << setw(outprec+6) << 0.0 << endl;
  sysout << setprecision(outprec) << setw(outprec+6) << 0.0
	 << setprecision(outprec) << setw(outprec+6) << 0.0
	 << setprecision(outprec) << setw(outprec+6) << L << endl;
  sysout << " }" << endl;
  */
  // There are a few places in HEG system (notably calculation of e-e
  // distances) that are optimized for orthogonal simulation cell
  // and general lattice vectors would lead to incorrect answers.
  sysout << " boxsize {" << endl;
  sysout << setprecision(outprec) << setw(outprec+6) << L
	 << setprecision(outprec) << setw(outprec+6) << L
	 << setprecision(outprec) << setw(outprec+6) << L << endl;
  sysout << " }" << endl;
  sysout << " origin { 0.0 0.0 0.0 }" << endl;
#ifdef PBC_ONLY
  sysout << " kpoint { 0.0 0.0 0.0 }" << endl;
#else
  {
    int inverse[2]={1,0};
    sysout << " kpoint { "
	   << inverse[tw[twist][0]] << " "
	   << inverse[tw[twist][1]] << " "
	   << inverse[tw[twist][2]] << " }" << endl;
  }
#endif
  sysout << " interaction { truncCoul }" << endl;
  sysout << "}" << endl;
  sysout.close();
  // }}}
}

void write_jast2_short_range_part(ofstream& jastout, double cutoff) {
  // {{{ .
  jastout << "  # e-e cusp conditions" << endl;
  jastout << "  GROUP {" << endl;
  jastout << "    OPTIMIZEBASIS" << endl;
  jastout << "    TWOBODY_SPIN {" << endl;
  jastout << "      FREEZE" << endl;
  jastout << "       LIKE_COEFFICIENTS   { 0.25  0.00 }" << endl;
  jastout << "       UNLIKE_COEFFICIENTS { 0.00  0.50 }" << endl;
  jastout << "    }" << endl;
  jastout << "    EEBASIS {" << endl;
  jastout << "      EE" << endl;
  jastout << "      CUTOFF_CUSP" << endl;
  jastout << "      GAMMA 0.1" << endl;   // maybe zero here
  jastout << "      CUSP 1" << endl;
  jastout << "      CUTOFF " << cutoff << endl;
  jastout << "    }" << endl;
  jastout << "    EEBASIS {" << endl;
  jastout << "      EE" << endl;
  jastout << "      CUTOFF_CUSP" << endl;
  jastout << "      GAMMA 0.5" << endl;
  jastout << "      CUSP 1" << endl;
  jastout << "      CUTOFF " << cutoff << endl;
  jastout << "    }" << endl;
  jastout << "  }" << endl << endl;
  jastout << "  # isotropic (short-range) e-e term" << endl; 
  jastout << "  GROUP {" << endl;
  jastout << "    OPTIMIZEBASIS" << endl;
  jastout << "    TWOBODY {" << endl;
  jastout << "      COEFFICIENTS { 0.0 }" << endl;
  jastout << "    }" << endl;
  jastout << "    EEBASIS {" << endl;
  jastout << "      EE" << endl;
  jastout << "      POLYPADE" << endl;
  jastout << "      RCUT " << cutoff << endl;
  jastout << "      BETA0 -0.5" << endl;
  jastout << "      NFUNC 1" << endl;
  jastout << "    }" << endl;
  jastout << "  }" << endl << endl;
  // }}}
}

void write_jast2_simple(string oname, double L) {
  // {{{ .
  string jastname=oname+".jast2";
  ofstream jastout(jastname.c_str());
  double cutoff=L/2.001;  // will work as long as we use the simple cubic cell
  jastout << "JASTROW2" << endl << endl;
  write_jast2_short_range_part(jastout, cutoff);
  // }}}
}

void write_jast2(string oname, int twist, double L, double k0,
		 int* kx, int* ky, int* kz, int* idx, int* k2,
		 int* shells, int nshells, int basissize, int outprec) {
  // {{{ .
#ifdef PBC_ONLY  
  string jastname=oname+".jast2";
#else
  stringstream str_twist;
  str_twist << twist;
  string jastname=oname+"."+str_twist.str()+".jast2";
#endif
  ofstream jastout(jastname.c_str());

  double cutoff=L/2.001;  // will work as long as we use the simple cubic cell
  jastout << "JASTROW2" << endl << endl;
  write_jast2_short_range_part(jastout, cutoff);

  // long-range part uses one parameter per shell/star
  // N.B.: this is probably the right thing to do only in the case of
  // closed-shell occupation (of course, closed-shell with respect to given
  // boundary condition)
  int first_shell=0;
  int first_k_vec=0;
  int extra_vec=0;
  // cos(0) would be just normalization so it is not included
  if ( k2[idx[0]]==0 ) {
    first_shell=1;
    first_k_vec=shells[0];
    extra_vec=1;
  }
  // a decision is needed how many shells/stars will be included in the Jastrow
  // factor; using occupied shells seems reasonable
  int last_shell;
  for ( int i=0; i < nshells; i++ ) {
    last_shell=i;
    if ( shells[i]/2+extra_vec >= basissize ) break;
  }  

  jastout << "  # long-range e-e term (plane-wave expansion, k-point dependent)"
	  << endl;
  jastout << "  GROUP {" << endl;
  jastout << "    TWOBODY {" << endl;
  jastout << "      COEFFICIENTS { ";
  for ( int i=first_shell; i < last_shell+1; i++ ) jastout << "0.0 ";
  jastout << "}" << endl;
  jastout << "    }" << endl;
  jastout << "    EEBASIS {" << endl;
  jastout << "      EE" << endl;
  jastout << "      BASIS_GROUPS" << endl;
  for ( int i=first_shell; i < last_shell+1; i++) {
    jastout << "      BASIS_GROUP {" << endl;
    jastout << "        EE" << endl;
    jastout << "        COSINE" << endl;
    jastout << "        GVECTOR {" << endl;
    for ( int j=first_k_vec; j < shells[i]/2+extra_vec; j++ ) {
      jastout << "     ";
      jastout << setprecision(outprec) << setw(outprec+6) << k0*kx[idx[j]]
	      << setprecision(outprec) << setw(outprec+6) << k0*ky[idx[j]]
	      << setprecision(outprec) << setw(outprec+6) << k0*kz[idx[j]]
	      << endl;
      first_k_vec=j+1;
    }
    //first_k_vec=shells[i]/2+extra_vec;
    jastout << "        }" << endl;
    jastout << "      }" << endl;
  }
  jastout << "    }" << endl;
  jastout << "  }" << endl; 
  jastout.close();
  // }}}
}


// Local variables:
// folded-file: t
