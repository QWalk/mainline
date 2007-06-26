/*
 
Copyright (C) 2007 Lucas K. Wagner

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
#include <cstdio>
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cmath>
#include "JeepBasis.h"
using namespace std;



/*
  (format comment copied from T. Ogitsu and B. Militzer)
Jeep wf format

spin up data

int = 0
int nst
double[nst]               = occ[n], n=0, .. , nst-1   --number of electrons in each orbital
int ngw
complex<double>[nst*ngw]  = wf(t)
complex<double>[nst*ngw]  = wf(t-dt)
 
next lines only if nspin == 2:
 
spin down data

int = 0
int nst
double[nst]               = occ[n], n=0, .. , nst-1
int ngw
complex<double>[nst*ngw]  = wf(t)
complex<double>[nst*ngw]  = wf(t-dt)

*/


//---------------------------------------------
//read/write stuff
int read_int(FILE * file) {
  int i;
  fread(&i, sizeof(int), 1, file);
  return i;
}

double read_double(FILE * file) {
  double f;
  fread(&f, sizeof(double), 1, file);
  return f;
}

void write_int(FILE *outfile, int i) {
  fwrite(&i,sizeof(int),1,outfile);
}

void write_double(FILE *outfile, double d) {
  fwrite(&d,sizeof(double),1,outfile);
}

//---------------------------------------------------

int main(int argc, char ** argv) {

  if(argc < 4) {
    cout << "usage:  " << argv[0] 
	 << "  -ecut <ecut in Hartrees>  <base name> " << endl;
    cout << "other options: " << endl;
    cout << "--no_derivatives         Don't print the derivatives" << endl;
    cout << "--no_laplacian           Don't print the laplacian" << endl;
    exit(1);
  }

  bool do_dx=true;
  bool do_dy=true;
  bool do_dz=true;
  bool do_lap=true;

  double ecut=0;
  for(int i=1; i < argc-1; i++) {
    if(!strcmp(argv[i], "-ecut")) {
      if( (++i) < argc) {
	ecut=atof(argv[i]);
      }
      else {
	cout << "Need number after -ecut" << endl;
	exit(1);
      }
    }
    else if(!strcmp(argv[i], "--no_derivatives")) {
      do_dx=do_dy=do_dz=false;
    }
    else if(!strcmp(argv[i], "--no_laplacian") ) {
      do_lap=false;
    }
    else {
      cout << "don't understand argument " << argv[i] << endl;
      exit(1);
    }
  }
  string basename=argv[argc-1];
  string wf_name=basename+".wf";
  string sys_name=basename+".sys";




  int nspin=1;

  //Get the cell from the sys file
  D3vector cell;
  string dummy;
  ifstream sysfile(sys_name.c_str());
  while(sysfile >> dummy) {
    if(dummy=="set_cell") {
      sysfile >> cell.x >> cell.y >> cell.z;
      break;
    }
  }
  sysfile.close();
  
  
  //make the basis set

  RealBasis basis;
  if(!basis.resize(cell, cell, ecut)) {
    cout << "Basis resize failed " << endl;
    exit(1);
  }

  cout << basis.size()  << " plane waves " << endl;



  

  FILE * wfin = fopen(wf_name.c_str(), "r");
  FILE * wfoutdx, * wfoutdy, * wfoutdz,* wfoutlap;
  string lapfile=wf_name+"lap";
  string dxfile=wf_name+"dx";
  string dyfile=wf_name+"dy";
  string dzfile=wf_name+"dz";

  if(do_lap) {
    wfoutlap=fopen(lapfile.c_str(), "w");
  }
  if(do_dx) {
    wfoutdx=fopen(dxfile.c_str(), "w");
  }
  if(do_dy) {
    wfoutdy=fopen(dyfile.c_str(), "w");
  }
  if(do_dz) {
    wfoutdz=fopen(dzfile.c_str(), "w");
  }
    
  int nst;
  for(int s=0; s < nspin; s++) {
    int test=read_int(wfin);
    cout << "should be zero " << test << endl;
    nst=read_int(wfin);
    cout << "number of states " << nst << endl;
    double * occupation=new double[nst];
    for(int i=0; i< nst; i++) {
      occupation[i]=read_double(wfin);
      cout << "occupation " << i << "   " <<  occupation[i] << endl;
    }
    int ngw=read_int(wfin);
    assert(ngw==basis.size());
    cout << "Number of plane waves " << ngw << endl;


    //-------------------------------
    //Write header
    if(do_lap) {
      write_int(wfoutlap, test);
      write_int(wfoutlap, nst);
      for(int i=0; i< nst; i++) {
	write_double(wfoutlap, occupation[i]);
      }
      write_int(wfoutlap, ngw);
    }

    if(do_dx) {
      write_int(wfoutdx, test);
      write_int(wfoutdx, nst);
      for(int i=0; i< nst; i++) {
	write_double(wfoutdx, occupation[i]);
      }
      write_int(wfoutdx, ngw);
    }

    if(do_dy) {
      write_int(wfoutdy, test);
      write_int(wfoutdy, nst);
      for(int i=0; i< nst; i++) {
	write_double(wfoutdy, occupation[i]);
      }
      write_int(wfoutdy, ngw);
    }

    if(do_dz) {
      write_int(wfoutdz, test);
      write_int(wfoutdz, nst);
      for(int i=0; i< nst; i++) {
	write_double(wfoutdz, occupation[i]);
      }
      write_int(wfoutdz, ngw);
    }
    
    //------------------------------
    //read/write states times the appropriate g
    //note that we're just taking the current time wf.  
    for(int i=0; i < nst; i++) {

      cout << "State number " << i << endl;
      for(int j=0; j < ngw; j++) {
	double real, imag;
	real=read_double(wfin);
	imag=read_double(wfin);
	
	double gx=basis.gx[j];
	double gy=basis.gx[ngw+j];
	double gz=basis.gx[ngw+ngw+j];

	if(do_lap) {
	  double gsquared=gx*gx+gy*gy+gz*gz;
	  write_double(wfoutlap, -gsquared*real);
	  write_double(wfoutlap, -gsquared*imag);
	}
	
	if(do_dx) {
	  write_double(wfoutdx, -gx*imag);
	  write_double(wfoutdx, gx*real);
	}
	if(do_dy) {
	  write_double(wfoutdy, -gy*imag);
	  write_double(wfoutdy, gy*real);
	}
	if(do_dz) {
	  write_double(wfoutdz, -gz*imag);
	  write_double(wfoutdz, gz*real);
	}

      }
    }

    delete [] occupation;
    
  }

  if(do_dx) {
    fclose(wfoutdx);
  }
  if(do_dy) {
    fclose(wfoutdy);
  }
  if(do_dz) {
    fclose(wfoutdz);
  }
  if(do_lap) {
    fclose(wfoutlap);
  }
  fclose(wfin);


  ofstream jeeplot("plot.in");
  jeeplot << sys_name << endl;
  jeeplot << "set ecut " << ecut*2 << endl;

  
  jeeplot << "\nload " << wf_name << endl;
  for(int i=0; i< nst; i++) {
    jeeplot << "plot orbital " << i+1 << " orb" << i+1 << ".plt\n";
  }
  jeeplot << endl;

  if(do_dx) {
    jeeplot << "\nload " << dxfile << endl;
    for(int i=0; i< nst; i++) {
      jeeplot << "plot orbital " << i+1 << " orb" << i+1 << "dx.plt\n";
    }
  }
  jeeplot << endl;

  if(do_dy) {
    jeeplot << "\nload " << dyfile << endl;
    for(int i=0; i< nst; i++) {
      jeeplot << "plot orbital " << i+1 << " orb" << i+1 << "dy.plt\n";
    }
  }
  jeeplot << endl;

  if(do_dz) {
    jeeplot << "\nload " << dzfile << endl;
    for(int i=0; i< nst; i++) {
      jeeplot << "plot orbital " << i+1 << " orb" << i+1 << "dz.plt\n";
    }
  }
  jeeplot << endl;

  if(do_lap) {
    jeeplot << "\nload " << lapfile << endl;
    for(int i=0; i< nst; i++) {
      jeeplot << "plot orbital " << i+1 << " orb" << i+1 << "lap.plt\n";
    }
  }
  jeeplot << endl;


  jeeplot.close();

}


