/*
 
Copyright (C) 2007 Jindrich Kolorenc, Michal Bajdich

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


//--------------------------------------------------------------------------
//
// The reason for this exercise is to efficiently evaluate molecular orbitals,
// when a molecular orbital is directly a basis function. Good for homogeneous
// electron gas and related stuff. Real version of this MO evaluator was
// constructed as a simplification of MO_matrix_standard. Complex version
// is a modification of the real version.
//
// Obviously, orb-file format is different from standard, cutoff and blas MOs.
// Here, each MO is defined by a pair AO# Center# -> two columns of integers
// in the orb-file.
//
//--------------------------------------------------------------------------

#ifndef MO_MATRIX_CBASFUNC_H_INCLUDED
#define MO_MATRIX_CBASFUNC_H_INCLUDED

#include "Array.h"
#include "Qmc_std.h"
#include "CBasis_function.h"
#include "Center_set.h"
#include "MO_matrix.h"

class System;
class Sample_point;
//----------------------------------------------------------------------------

class MO_matrix_Cbasfunc: public Complex_MO_matrix
{
protected:
  void init();
private:
  
  Center_set centers;
  Array1 <CBasis_function *> basis;
  int nmo;
  doublevar magnification_factor;
  string orbfile;
  int totbasis;
  int maxbasis;

  // basisMO replaces moCoeff, links given MO to its basis function
  Array2 <int> basisMO;
  // moLists(spin,.) lists MO's corresponding to the given spin,
  // the same thing as totoccupation in Slat_wf_data
  Array1 <int> eq_centers;
  Array1 < Array1 <int> > moLists;
  Array1 <dcomplex>  kptfac;
  
  Array1 <doublevar> obj_cutoff; //!< cutoff for each basis object
  Array1 <doublevar> cutoff;  //!< Cutoff for individual basis functions
    
  Array1 <doublevar> kpoint; //!< the k-point of the orbitals in fractional units(1 0 0) is X, etc..

public:

  /*!
    Build several sets of MO's to be evaluated in updateVal and updateLap.
    Each element in occupations should be a list of the MO's that should
    be evaluated.  For example, one can create a list of spin up and spin
    down MO's, and only evaluate up when an up electron is moved.
   */
  virtual void buildLists(Array1 <Array1 <int> > & occupations);


  virtual int showinfo(ostream & os);

  virtual int writeinput(string &, ostream &);

  /*!
    get the number of molecular orbitals
  */
  int getNmo() {
    return nmo;
  }

  virtual void read(vector <string> & words, unsigned int & startpos, 
                    System * sys);

  // finally, the two key functions of the class that actually evaluate the
  // molecular orbitals (and their derivatives)
  virtual void updateVal(Sample_point * sample,
			 int e,
			 //!< electron number
			 int listnum,
			 //!< Choose the list that was built in buildLists
			 Array2 <dcomplex> & newvals
			 //!< The return: in form (MO)
			 );
  
  virtual void updateLap(Sample_point * sample,
			 int e,
			 //!< electron number
			 int listnum,
			 //!< Choose the list that was built in buildLists
			 Array2 <dcomplex> & newvals
			 //!< The return: in form ([value gradient lap], MO)
			 );

  virtual void updateHessian(Sample_point * sample,
			     int e,
			     //!< electron number
			     int listnum,
			     //!< Choose the list that was built in buildLists
			     Array2 <dcomplex>& newvals
			     //!< in form ([value gradient, dxx,dyy,dzz,dxy,dxz,dyz], MO)
			     ) {
    error("CBASFUNC_MO: updateHessian not implemented/adopted.");
  }

  MO_matrix_Cbasfunc()
  {}

};


#endif // MO_MATRIX_CBASFUNC_H_INCLUDED

//--------------------------------------------------------------------------
