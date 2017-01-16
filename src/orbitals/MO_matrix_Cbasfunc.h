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

/*!
The reason for this exercise is to efficiently evaluate molecular orbitals,
when a molecular orbital is directly a basis function. Good for homogeneous
electron gas and related stuff (HEG_system). This MO evaluator is constructed
as a retype of MO_matrix_basfunc, which in turn is a simplification of
MO_matrix_standard.

There are some extra declarations in MO_matrix_Cbasfunc compared to its real
valued version MO_matrix_basfunc, since Complex_MO_matrix lacks some
functionality of MO_matrix. Maybe it would have been better to enhance
Complex_MO_matrix than to put it here. After all, we want to have more complex
functionality in the future.
*/

//--------------------------------------------------------------------------

class MO_matrix_Cbasfunc: public Complex_MO_matrix
{
protected:
  void init();
private:
  
  // transferred from MO_matrix
   Array1 <CBasis_function *> basis;
  // -----

  Array2 <int> basisMO;
  //!< basisMO replaces moCoeff, links given MO to its basis function
  Array1 <int> eq_centers;
  Array1 < Array1 <int> > moLists;
  /*!< \brief moLists(spin,.) lists MO's corresponding to the given spin,
    the same thing as totoccupation in Slat_wf_data
  */
  Array1 <dcomplex>  kptfac;
  /*!< \brief
    phase factor multiplying basis functions associated with
    a given center, N.B. equivalent centers differ by integer
    multiple of a lattice vector
  */
  
  Array1 <doublevar> obj_cutoff; //!< cutoff for each basis object
  Array1 <doublevar> cutoff;     //!< cutoff for individual basis functions
    
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
