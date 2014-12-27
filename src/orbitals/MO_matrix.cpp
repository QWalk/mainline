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

#include "Qmc_std.h"
#include "MO_matrix.h"
#include "MO_matrix_cutoff.h"
#include "MO_matrix_standard.h"
#include "Sample_point.h"
#include "qmc_io.h"
#include "MO_1d.h"
#include "MO_matrix_blas.h"
#include "MO_matrix_basfunc.h"
#include "MO_matrix_Cbasfunc.h"
#include "MO_matrix_einspline.h"
#include <algorithm>

int allocate(vector <string> & words, System * sys, MO_matrix *& moptr) {
  assert(moptr==NULL);

  if(caseless_eq(words[0],"CUTOFF_MO"))
    moptr=new MO_matrix_cutoff<doublevar>;
  else if(caseless_eq(words[0],"STANDARD_MO"))
    moptr=new MO_matrix_standard;
  else if(caseless_eq(words[0],"BLAS_MO"))
    moptr=new MO_matrix_blas;
  else if(caseless_eq(words[0],"BASFUNC_MO"))
    moptr=new MO_matrix_basfunc;
  else if(caseless_eq(words[0],"EINSPLINE_MO"))
    moptr=new MO_matrix_einspline<doublevar>;
  
  else {
    error("Didn't  understand ",words[0]);
  }

  unsigned int pos=0;
  moptr->read(words,pos, sys);
  return 1;
}

int allocate(vector <string> & words, System * sys, 
             Complex_MO_matrix *& moptr) {
  if(caseless_eq(words[0],"MO_1D"))
    moptr=new MO_1d;
  else if(caseless_eq(words[0],"CBASFUNC_MO"))
    moptr=new MO_matrix_Cbasfunc;
  else if(caseless_eq(words[0],"CUTOFF_MO"))
    moptr=new MO_matrix_cutoff<dcomplex>;
  else if(caseless_eq(words[0],"EINSPLINE_MO"))
    moptr=new MO_matrix_einspline<dcomplex>;
  else 
    error("Unknown complex MO: ", words[0]);

  unsigned int pos=0;
  moptr->read(words, pos, sys);
  return 1;
}


//----------------------------------------------------------------------
/*!

 */
void rotate_orb(istream & orbin, ostream & orbout,
                Array2 <doublevar> & rotation,
                Array1 <int>  & moList, int nfunctions) {
  int nmo_write=moList.GetDim(0);
  assert(nmo_write==rotation.GetDim(1));
  assert(nmo_write==rotation.GetDim(0));
  if(!orbin) error("couldn't open orb input file");
  if(!orbout) error("couldn't open orb output file");



  Array2 <doublevar> rotated_mo(nmo_write, nfunctions);

  string dummy;
  while (orbin >> dummy) {
    if(dummy=="COEFFICIENTS") break;
  }
  if(!orbin) error("rotate_orb::Didn't find COEFFICIENTS in orbin");
  int nmo_read=0;
  for(int i=0; i< nmo_write; i++) {
    if(moList(i) > nmo_read) nmo_read=moList(i);
  }
  nmo_read++;

  cout << "nmo_write " << nmo_write << "  nmo_read " << nmo_read << endl;

  Array2 <doublevar> moCoeff(nmo_read, nfunctions);

  //bug here.
  for(int mo=0; mo < nmo_read; mo++) {
    for(int f=0; f< nfunctions; f++) {
      if(!(orbin >> moCoeff(mo,f))) error("rotate_orb::orb file ended before I expected.");
    }
  }

  debug_write(cout,"rotating mo's \n");
  rotated_mo=0;
  for(int m=0; m< nmo_write; m++) {
    
      for(int f=0; f< nfunctions; f++) {
        for(int m2=0; m2 < nmo_write; m2++) {
          int mo=moList(m2);
          //          cout << "m " << m << " f " << f  << " m2 " << m2 << endl;
          rotated_mo(m, f)+=rotation(m,m2)*moCoeff(mo,f);
        }
      }
  }

  debug_write(cout,"outputing mo's \n");
  for(int m=0; m < nmo_write; m++) {
    int counter2=1;
    for(int f=0; f< nfunctions; f++) {
      orbout << rotated_mo(m, f) << "   ";
      if(counter2 % 5 ==0) orbout << endl;
      counter2++;
    }
    orbout << endl;
  }


}

//----------------------------------------------------------------------


/*!
Implementation of 'rotate_orb' for complex orbital coefficients.
 */
void rotate_Corb(istream & orbin, ostream & orbout,
		 Array2 <doublevar> & rotation,
		 Array1 <int>  & moList, int nfunctions) {
  int nmo_write=moList.GetDim(0);
  assert(nmo_write==rotation.GetDim(1));
  assert(nmo_write==rotation.GetDim(0));
  if(!orbin) error("couldn't open orb input file");
  if(!orbout) error("couldn't open orb output file");

  Array2 <dcomplex> rotated_mo(nmo_write, nfunctions);

  string dummy;
  while (orbin >> dummy) {
    if(dummy=="COEFFICIENTS") break;
  }
  if(!orbin) error("rotate_orb::Didn't find COEFFICIENTS in orbin");
  int nmo_read=0;
  for(int i=0; i< nmo_write; i++) {
    if(moList(i) > nmo_read) nmo_read=moList(i);
  }
  nmo_read++;

  Array2 <dcomplex> moCoeff(nmo_read, nfunctions);

  //bug here. ??? 
  for(int mo=0; mo < nmo_read; mo++) {
    for(int f=0; f< nfunctions; f++) {
      if(!(orbin >> moCoeff(mo,f))) error("rotate_orb::orb file ended before I expected.");
    }
  }

  debug_write(cout,"rotating mo's \n");
  rotated_mo=0;
  for(int m=0; m< nmo_write; m++) {
    
      for(int f=0; f< nfunctions; f++) {
        for(int m2=0; m2 < nmo_write; m2++) {
          int mo=moList(m2);
          //          cout << "m " << m << " f " << f  << " m2 " << m2 << endl;
          rotated_mo(m, f)+=rotation(m,m2)*moCoeff(mo,f);
        }
      }
  }

  debug_write(cout,"outputing mo's \n");
  for(int m=0; m < nmo_write; m++) {
    int counter2=1;
    for(int f=0; f< nfunctions; f++) {
      orbout << rotated_mo(m, f) << "   ";
      if(counter2 % 5 ==0) orbout << endl;
      counter2++;
    }
    orbout << endl;
  }


}

//----------------------------------------------------------------------


/*!

*/
//------------------------------------------------------------------------------------------
#ifdef USE_MPI
void overloaded_broadcast(Array1 <doublevar> & v) { 
  MPI_Bcast(v.v,v.GetDim(0), MPI_DOUBLE,0,MPI_Comm_grp);
}
void overloaded_broadcast(Array1 <dcomplex> & v) { 
  MPI_Bcast(v.v,v.GetDim(0), MPI_DOUBLE_COMPLEX,0,MPI_Comm_grp);
}
#endif

template <class T> int readorb(istream & input, Center_set & centers, 
                                  int nmo, int maxbasis, Array1 <doublevar> & kpoint,
                                  Array3 <int> & coeffmat, Array1 <T> & coeff) {
  int nmo_read=0;
  int maxlabel=0; 
  coeffmat.clear(); //important to do this so that we know exactly how big the array v will be
                    //This enables us to use relatively fast Bcast() operations.
  coeff.clear();
  if(mpi_info.node==0) { 
    string dummy;
    vector <int> mo,center,basis,label;
    while(true) { 
      input >> dummy;
      if(dummy=="COEFFICIENTS") break;
      int currmo=atoi(dummy.c_str())-1;
      if(currmo > nmo) { 
        input.ignore(300,'\n');
      }
      else { 
        mo.push_back(currmo);
        input >> dummy;
        basis.push_back(atoi(dummy.c_str())-1);
        if(basis.back() >= maxbasis) 
          error("Basis function greater than maxbasis requested:",basis.back()+1);
        else if(basis.back() < 0) 
          error("Basis function cannot be less than 1:",basis.back()+1);
        input >> dummy;
        center.push_back(atoi(dummy.c_str())-1);
        if(center.back() > centers.equiv_centers.GetDim(0) )
          error("Center number in orb file greater than maximum number ", 
                 centers.equiv_centers.GetDim(0));
        
        input >> dummy;
        label.push_back(atoi(dummy.c_str())-1);
      }
      if(!input) 
        error("Unexpected end of file; did not find COEFFICIENTS while reading orbitals");
    }
    nmo_read=*std::max_element(mo.begin(),mo.end())+1;
    coeffmat.Resize(nmo_read, centers.size(), maxbasis);
    coeffmat=-1;
    {
      vector<int>::iterator m=mo.begin(),
        c=center.begin(),
        b=basis.begin(),
        l=label.begin();

      for(  ; m!=mo.end() && c!=center.end() && b!=basis.end() && l!=label.end();
          m++,c++,b++,l++) { 
//        coeffmat(*m,*c,*b)=*l;
        for(int c_eq=0; c_eq < centers.ncenters_atom(*c); c_eq++) {
          int cen2=centers.equiv_centers(*c, c_eq);
          coeffmat(*m, cen2,*b)=*l;
        }
        
      }
    }

    maxlabel=*std::max_element(label.begin(),label.end())+1;
    coeff.Resize(maxlabel);
    for(int i=0; i< maxlabel; i++) { 
      if(! (input >> coeff(i) ) )
        error("unexpected end of file when reading orbital coefficients");
    }
  }
#ifdef USE_MPI
  MPI_Bcast(&nmo_read,1,MPI_INT,0,MPI_Comm_grp);
  MPI_Bcast(&maxlabel,1,MPI_INT,0,MPI_Comm_grp);
  int coeffmatsize;
  
  if(mpi_info.node!=0) { 
    coeffmat.Resize(nmo_read,centers.size(),maxbasis);
    coeff.Resize(maxlabel);
  }
  MPI_Bcast(coeffmat.v,coeffmat.size,MPI_INT,0,MPI_Comm_grp);
  overloaded_broadcast(coeff);
#endif
  return nmo_read;
}



//----------------------------------------------------------------------
