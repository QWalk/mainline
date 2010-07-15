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

//#define PRINT_VECS   1  // to print k-vectors of basis functions
//#define PRINT_SHELLS 1  // to print (closed) shells
#define TETRAHEDRON  1  // to use tetrahedron method to generate k-points

using namespace std;

void write_basis(string oname, string twist,
                 double* kx, double* ky, double* kz, int* idx,
		 int basissize, int outprec);
void write_slater(string oname, string twist, int* Ne);
void write_orb(string oname, string twist, int nmo);
void write_sys(string oname, string twist, int* Ne, double L,
	       double jKx, double jKy, double jKz, double weight_tw,
               int outprec);
void write_jast2(string oname, double L);

extern "C" {
  // soubroutine from TETPACK library
  void setk01_(int* na, double* a, double* ptk, int* nptk, int*idef,
	       int* ntet, int* nkmax, int* ntmax);
}

int main(int argc, char ** argv) {
  // {{{ .

  const int ikmax=8;    // upper limit for planewave search
  const double pi=3.1415926535897932385;

  // find out how many decimal places has double and set precision
  // for printing k-vectors and cell size
  typedef numeric_limits< double > dl;
  int outprec=dl::digits10;  
  //cout << "double has " << outprec << " digits" << endl;
  
  // base for output files is given on the command line
  string oname;
  bool oflag=false;
  for(int i=1; i< argc; i++) {
    if(!strcmp(argv[i], "-o") && argc > i+1) {
      oname=argv[++i];
      oflag=true;
    }
    else if(argv[i][0]=='-') {
      cerr << "Unknown option: " << argv[i] << endl;
    }
  }

  if ( !oflag ) {
    cout << "Usage: heg2qmcC -o <output base> [-bcs]" << endl;
    cout << endl;
    cout << "Takes N_up, N_down and rs from standard input and ";
    cout << "writes" << endl;
    cout << "   <output base>.i.sys" << endl;
    cout << "   <output base>.i.slater" << endl;
    cout << "   <output base>.i.basis" << endl;
    cout << "   <output base>.i.orb" << endl;
    cout << "   <output base>.centers" << endl;
    cout << "   <output base>.jast2" << endl;
    cout << "where i indexes k-points (twists of boundary conditions) and ";
    cout << "runs from 0 to 3;" << endl;
    cout << "0 is the Gamma-point." << endl;
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

  cout << "rs (in atomic units, i.e., in bohrs)" << endl;
  double rs;
  cin >> rs;
#ifdef TETRAHEDRON
  cout << "K-mesh discretization parameter (tetrahedron method)" << endl;
#else
  cout << "Number of K-points in each direction" << endl;
#endif
  int NK;
  cin >> NK;
  cout << endl;
  double L=rs*pow(4.0/3*(Ne[0]+Ne[1])*pi,1.0/3);
  double k0=2*pi/L;

  //  for (NK=2; NK<40; NK++) {
  //    cout << "---------------------------------------------------------"
  //	 << "-----------------------" << endl;
  //    cout << "NK=" << NK << endl;
 
  double Ekin=0;     // Kinetic energy for constructed slater wavefunction
  double weight=0;

  int NKpt;

#ifdef TETRAHEDRON

  int ntet;
  int nkmax=(NK+1)*(NK+2)*(NK+3)/6;
  int ntmax=NK*NK*NK;
  int idef[ntmax*5];
  double Kpt[nkmax*4];
  setk01_(&NK,&L,Kpt,&NKpt,idef,&ntet,&nkmax,&ntmax);
  cout << endl;

#else // TETRAHEDRON
  
  NKpt=NK*NK*NK;
  double Kpt[NKpt*4];
  int Kcount=0;
  for (int iKx=0; iKx<NK; iKx++) {
    for (int iKy=0; iKy<NK; iKy++) {
      for (int iKz=0; iKz<NK; iKz++) {
	Kpt[Kcount*4+0]=iKx*k0/NK;
	Kpt[Kcount*4+1]=iKy*k0/NK;
	Kpt[Kcount*4+2]=iKz*k0/NK;
	Kpt[Kcount*4+3]=1.0;
	Kcount++;
      }
    }
  }

#endif // TETRAHEDRON

  cout << "NKpt=" << NKpt << endl;

  for (int iK=0; iK<NKpt; iK++) {
    stringstream twist;
    twist << iK;
    double Kx=Kpt[iK*4+0];
    double Ky=Kpt[iK*4+1];
    double Kz=Kpt[iK*4+2];
    double weight_tw=Kpt[iK*4+3];

    double Ekin_tw=0;

  // now we construct all k-vectors in a cube and sort them with
  // respect to their magnitude; then we occuppy as many of them as
  // needed

  int ikmax3=(2*ikmax-1)*(2*ikmax-1)*(2*ikmax-1);
  double kx[ikmax3], ky[ikmax3], kz[ikmax3], k2[ikmax3];
  int idx[ikmax3], shells[ikmax3];

  int seq=0;

  for ( int i=-ikmax+1; i<ikmax; i++) {
    for ( int j=-ikmax+1; j<ikmax; j++) {
      for ( int k=-ikmax+1; k<ikmax; k++) {
	kx[seq]=Kx+i*k0;
	ky[seq]=Ky+j*k0;
	kz[seq]=Kz+k*k0;
	k2[seq]=kx[seq]*kx[seq]+ky[seq]*ky[seq]+kz[seq]*kz[seq];
	seq++;
      }
    }
  }  

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
    // "test" cube, throw these vectors away (the -2 might be too conservative
    // but there is the K-point shift we have to consider
    if ( k2[idx[seq]] > (ikmax-2)*(ikmax-2)*k0*k0 ) break;
    
#ifdef PRINT_VECS
    cout << setw(4) << seq 
	 << setw(4) << kx[idx[seq]] 
	 << setw(4) << ky[idx[seq]] 
	 << setw(4) << kz[idx[seq]] 
	 << setw(5) << k2[idx[seq]] << endl;
#endif

    inshell++;

    if ( seq+1 < ikmax3 && k2[idx[seq+1]] > k2[idx[seq]] ) {
      if ( shell > 0 ) { 
	shells[shell]=inshell+shells[shell-1];
      } else {
	shells[0]=inshell;
      }
#ifdef PRINT_VECS
      cout << "----- " << setw(2) << shell+1 << " ("
	   << setw(4) << shells[shell] << ") -----" << endl;
#endif
      inshell=0;
      shell++;
    }  
  }  // for ( seq=0; seq < ikmax3; seq++ )
  int nshells=shell;

#ifdef PRINT_SHELLS
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
#endif // PRINT_SHELLS

  
  // test whether user requested more states than we have found
  if ( Ne[0] > shells[nshells-1] || Ne[1] > shells[nshells-1] ) {
    cerr << "Maximum number of electrons of one kind is currently "
	 << shells[nshells-1] << endl;
    cerr << "To go higher, change hardcoded limit of basis size (ikmax)."
	 << endl;
    exit(1);
  }

  int basissize=nmo;

  // write basis set file
  write_basis(oname, twist.str(), kx, ky, kz, idx, basissize, outprec);

  // write slater determinant
  write_slater(oname, twist.str(), Ne);
  
  // write orb file
  write_orb(oname, twist.str(), nmo);

  // write sys file
  write_sys(oname, twist.str(), Ne, L, 2*Kx/k0, 2*Ky/k0, 2*Kz/k0,
	    weight_tw, outprec); 

  // calculate the (kinetic) energy of non-interacting particles
  // to have an idea of this part of finite-size errors
  for (int i=0; i<2; i++) {
    for (int j=0; j<Ne[i]; j++) {
      double dEkin=k2[idx[j]]/2;
      //cout << dEkin << endl;
      Ekin_tw+=dEkin;
    }
  }
  
  // kinetic energy for a given twist
  /*
  cout << Kx*NK*L/pi << " " << Ky*NK*L/pi << " " << Kz*NK*L/pi
       << " (weight "
       << setiosflags(ios::fixed) << weight_tw << resetiosflags(ios::fixed)
       << "),  ";
  cout << "kinetic energy: " 
       << Ekin_tw << " ("
       << Ekin_tw/(Ne[0]+Ne[1]) << " per particle)" << endl;
  */
  //cout << endl << "---" << endl << endl;
  Ekin+=Ekin_tw*weight_tw;
  weight+=weight_tw;
  //cout << weight_tw << " " << weight <<endl;;
    
  Ekin_tw=0;

  }

  // write the file with definition of a single center
  string centername=oname+".centers";
  ofstream centerout(centername.c_str());
  centerout << "1" << endl; 
  centerout << "origin 0.0 0.0 0.0" << endl;
  centerout.close();

  // write simple two-body Jastrow
  write_jast2(oname, L);

  // compare our approximate kinetic energy with exact one
  double Ekin_ex=0;
  for (int i=0; i<2; i++) {
    double kF=pow(6*pi*pi*Ne[i],1.0/3)/L;
    Ekin_ex+=L*L*L*pow(kF,5.0)/(20*pi*pi);
  }
  Ekin=Ekin/weight;
  cout << setprecision(8);
  cout << endl << "Total kinetic energy (hartree):" << endl;
  cout << "  approximate: "
       << setw(12) << Ekin << "  ("
       << setw(12) << Ekin/(Ne[0]+Ne[1])
       << " per particle)" << endl;
  cout << "  exact      : " 
       << setw(12) << Ekin_ex << "  ("
       << setw(12) << Ekin_ex/(Ne[0]+Ne[1])
       << " per particle)" << endl;
  cout << "  error      : " 
       << setw(12) << 100*abs((Ekin-Ekin_ex)/Ekin) << "%"
       << endl << endl;

  Ekin=0.0;

  //  } for (NK=2; NK<40; NK++)

  cout << "---" << endl << endl;
  cout << "BTW, magnification factor in *.slater might need to be lowered for "
       << "large" << endl << "number of particles (cca 200+ per spin channel)."
       << endl << endl;

  return 0;

  // }}}
}


void write_basis(string oname, string twist,
                 double* kx, double* ky, double* kz, int* idx,
		 int basissize, int outprec) {
  // {{{ .

  string basisname=oname+"."+twist+".basis";
  ofstream basisout(basisname.c_str());
  basisout << "  CBASIS {" << endl;
  basisout << "   origin" << endl;
  basisout << "   CPLANEWAVE" << endl; 
  basisout << "   GVECTOR {" << endl;
  
  for (int i=0; i < basissize; i++) {
    basisout << "   ";
    basisout << setprecision(outprec) << setw(outprec+6) << kx[idx[i]]
             << setprecision(outprec) << setw(outprec+6) << ky[idx[i]]
             << setprecision(outprec) << setw(outprec+6) << kz[idx[i]]
	     << endl;
  }

  basisout << "   }" << endl;
  basisout << "  }" << endl;  
  basisout.close();

  // }}}
}

void write_slater(string oname, string twist, int* Ne) {
  // {{{ .

  int nmo_slater;
  if ( Ne[0] > Ne[1] ) { nmo_slater=Ne[0]; } else { nmo_slater=Ne[1]; }

  string slatername=oname+"."+twist+".slater";
  ofstream slaterout(slatername.c_str());
  slaterout << " SLATER" << endl;
  slaterout << " CORBITALS {" << endl;
  slaterout << "  CBASFUNC_MO" << endl;
  slaterout << "  MAGNIFY 0.5" << endl;
  slaterout << "  INCLUDE " << oname;
  slaterout << "." << twist;
  slaterout << ".basis" << endl;
  slaterout << "  CENTERS { READ " << oname << ".centers }" << endl;
  slaterout << "  NMO " << nmo_slater << endl;
  slaterout << "  ORBFILE " << oname;
  slaterout << "." << twist;
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


void write_orb(string oname, string twist, int nmo) {
  // {{{ .
  string orbname=oname+"."+twist+".orb";
  ofstream orbout(orbname.c_str());
  int istart=1;
  int istop=nmo+1;
  for( int i=istart; i < istop; i++) {
    orbout << i << "  " << 1 << endl;
  }
  orbout.close();
  
  // }}}
}

void write_sys(string oname, string twist, int* Ne, double L,
	       double jKx, double jKy, double jKz, double weight_tw,
               int outprec) {
  // {{{ .
  string sysname=oname+"."+twist+".sys";
  ofstream sysout(sysname.c_str());
  sysout << "# weight of this k-point (twist): " << weight_tw << endl;
  sysout << "SYSTEM { HEG" << endl;
  sysout << " NSPIN { " << Ne[0] << " " << Ne[1] << " }" << endl;
  // There are a few places in HEG system (notably calculation of e-e
  // distances) that are optimized for orthogonal simulation cell
  // and general lattice vectors would lead to incorrect answers.
  sysout << " boxsize {" << endl;
  sysout << setprecision(outprec) << setw(outprec+6) << L
	 << setprecision(outprec) << setw(outprec+6) << L
	 << setprecision(outprec) << setw(outprec+6) << L << endl;
  sysout << " }" << endl;
  sysout << " origin { 0.0 0.0 0.0 }" << endl;
  // all kpoint related information should be in the basis itself
  //  <-- this is probably NOT TRUE - phase change across boundary
  sysout << " kpoint { " << jKx << " " << jKy << " " << jKz << " }" << endl;
  sysout << " interaction { truncCoul }" << endl;
  sysout << "# interaction { coulomb } " << endl;
  sysout << "}" << endl;
  sysout.close();
  // }}}
}

void write_jast2(string oname, double L) {
  // {{{ we could certainly use RPA to obtain better initial guess,
  // but right now I just need something
  string jastname=oname+".jast2";
  ofstream jastout(jastname.c_str());
  double cutoff=L/2.001;  // will work as long as we use the simple cubic cell
  jastout << "JASTROW2" << endl;
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
  jastout << "  }" << endl;
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
  jastout << "  }" << endl;
  jastout.close();
  // }}}
}


// Local variables:
// folded-file: t
