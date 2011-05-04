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
#include "MO_matrix_bspline.h"
#include "MO_matrix_Cbspline.h"
#include "MO_matrix_Cbasfunc.h"
#include "MO_matrix_Ccutoff.h"
#include "MO_matrix_einspline.h"
#define COMPLEX_WF
#include "MO_matrix_einspline.h"
#undef COMPLEX_WF

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
  else if(caseless_eq(words[0],"BSPLINE_MO"))
    moptr=new MO_matrix_bspline;
  else if(caseless_eq(words[0],"EINSPLINE_MO"))
    moptr=new MO_matrix_einspline;
  
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
  else if(caseless_eq(words[0],"CBSPLINE_MO"))
    moptr=new MO_matrix_Cbspline;
  else if(caseless_eq(words[0],"EINSPLINE_MO"))
    moptr=new MO_matrix_Ceinspline;
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

  //cout << "nmo_write " << nmo_write << "  nmo_read " << nmo_read << endl;

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


/*!
foo.readorb(input stream)
 Gets the coefficients and sets everything up from a formatted
 input file(ORB).  The file should look like: <br>
 MO#  AO#(for center) Center# Coeff# <br>
 for all MO's, then a listing of the coefficients in sequential order

 All the listings for a given MO must be in one block.

\todo
Check the sums of functions versus how many we think should
be there, and make sure everything adds up correctly.

*/
int readorb(istream & input, Center_set & centers, 
            int nmo, int maxbasis, Array1 <doublevar> & kpoint, Array3 <int > & coeffmat, 
            Array1 <doublevar> & coeff
           )
{
  
  int center;
  int max=0;
  int mo=0;
  int molast=0;
  int fn;

  //cout << "readorb" << endl;
  int ncenters=centers.ncenters_atom.GetDim(0);
  
  coeffmat.Resize(nmo, centers.size(), maxbasis);
  //cout << "nmo " << nmo << "  ncenters "
  //    << centers.size() << "  maxbasis " << maxbasis << endl;
  coeffmat=-1;
  
  string dummy;
  input >> mo;

 
  for(int i=0; i< nmo; i++)
  {
    int firstDone=0;
    while(1)
    {
      if(firstDone)
      {
        if(!(input >> dummy)) error("Unexpected end of orbital file\n");
        mo=atoi(dummy.c_str());
      }


      if(mo != molast+1 &&  firstDone)
      {
        break;
      }
      mo-=1;

      
      if(mo != i)
      {
        cout << "mo " << mo << " i " << i << "   dummy " << dummy << endl;
        error("Bad formatting in the orb file, or not enough mo's there.");
      }

      if(!input) cout << " here" << endl;

      input >> fn;
      if(fn > maxbasis) error("Function in orb file greater than maximum number"
                              "of basis functions", fn);
      fn-=1;
      input >> center;
      if(center > centers.equiv_centers.GetDim(0) )
        error("Center number in orb file greater than maximum number ", 
               centers.equiv_centers.GetDim(0));
      center-=1;
       
      int coeffnum;

            
      input >> coeffnum;

      if(max < coeffnum) max=coeffnum;
      coeffnum-=1;

      if(coeffnum < 0)
         error("Coefficient pointer less than one in orb file");

      //loop through equivalent centers..

      if(center > ncenters) error("center number too high in orb file: ", center+1);
      //cout << " mo " << mo << " cen " << center << " fn " << fn << endl;
      
      for(int c=0; c < centers.ncenters_atom(center); c++) {
        int cen2=centers.equiv_centers(center, c);
        coeffmat(mo, cen2, fn)=coeffnum;
      }

      firstDone=1;
      molast=mo;
      
    }
  }
    
  //Find the coefficients section.
  while(dummy!="COEFFICIENTS")
  {
    if(!(input >> dummy))
      error("Couldn't find COEFFICIENTS section in the orb file.\n");
  }

  //cout << "half " << endl;

  coeff.Resize(max);
  //coeff=1e99;
  for(int i=0; i<max; i++)
  {
    //input >> dummy;
    //cout << dummy << endl;
    if(!(input >> coeff(i)))
      error("Didn't find all the MO coefficients I should've");
  }

  //cout << "done readorb" << endl;

  return max;
}

//----------------------------------------------------------------------------------


/*!
Version of 'readorb' that reads complex orbital expansion coefficients. Just a
plain copy of above real version, the only change is in declaration of parameters.
 */
int readorb(istream & input, Center_set & centers, 
            int nmo, int maxbasis, Array1 <doublevar> & kpoint,
	    Array3 <int > & coeffmat, Array1 <dcomplex> & coeff
           )
{
  
  int center;
  int max=0;
  int mo=0;
  int molast=0;
  int fn;

  //cout << "readorb" << endl;
  int ncenters=centers.ncenters_atom.GetDim(0);
  
  coeffmat.Resize(nmo, centers.size(), maxbasis);
  //cout << "nmo " << nmo << "  ncenters "
  //    << centers.size() << "  maxbasis " << maxbasis << endl;
  coeffmat=-1;
  
  string dummy;
  input >> mo;

 
  for(int i=0; i< nmo; i++)
  {
    int firstDone=0;
    while(1)
    {
      if(firstDone)
      {
        if(!(input >> dummy)) error("Unexpected end of orbital file\n");
        mo=atoi(dummy.c_str());
      }


      if(mo != molast+1 &&  firstDone)
      {
        break;
      }
      mo-=1;

      
      if(mo != i)
      {
        cout << "mo " << mo << " i " << i << "   dummy " << dummy << endl;
        error("Bad formatting in the orb file, or not enough mo's there.");
      }

      if(!input) cout << " here" << endl;

      input >> fn;
      if(fn > maxbasis) error("Function in orb file greater than maximum number"
                              "of basis functions", fn);
      fn-=1;
      input >> center;
      if(center > centers.equiv_centers.GetDim(0) )
        error("Center number in orb file greater than maximum number ", 
               centers.equiv_centers.GetDim(0));
      center-=1;
       
      int coeffnum;

            
      input >> coeffnum;

      if(max < coeffnum) max=coeffnum;
      coeffnum-=1;

      if(coeffnum < 0)
         error("Coefficient pointer less than one in orb file");

      //loop through equivalent centers..

      if(center > ncenters) error("center number too high in orb file: ", center+1);
      //cout << " mo " << mo << " cen " << center << " fn " << fn << endl;
      
      for(int c=0; c < centers.ncenters_atom(center); c++) {
        int cen2=centers.equiv_centers(center, c);
        coeffmat(mo, cen2, fn)=coeffnum;
      }

      firstDone=1;
      molast=mo;
      
    }
  }
    
  //Find the coefficients section.
  while(dummy!="COEFFICIENTS")
  {
    if(!(input >> dummy))
      error("Couldn't find COEFFICIENTS section in the orb file.\n");
  }

  coeff.Resize(max);
  //coeff=1e99;
  for(int i=0; i<max; i++)
  {
    if(!(input >> coeff(i)))
      error("Didn't find all the MO coefficients I should've");
  }

  //cout << "done readorb" << endl;

  return max;
}


//----------------------------------------------------------------------
