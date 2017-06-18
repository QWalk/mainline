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

#ifndef MO_MATRIX_BASFUNC_H_INCLUDED
#define MO_MATRIX_BASFUNC_H_INCLUDED

#include "Array.h"
#include "Qmc_std.h"
#include "Basis_function.h"
#include "Center_set.h"
#include "MO_matrix.h"

class System;
class Sample_point;
//----------------------------------------------------------------------------

/*!
The reason for this exercise is to efficiently evaluate molecular orbitals,
when a molecular orbital is directly a basis function. Good for homogeneous
electron gas and related stuff (HEG_system). This MO evaluator is constructed
as a simplification of MO_matrix_standard. Equivalent complex-valued
functionality is provided by MO_matrix_Cbasfunc.

Obviously, orb-file format is different from MO_matrix_standard,
MO_matrix_cutoff and MO_matrix_blas. Here, each MO is defined by a pair
AO# center# -> two columns of integers in the orb-file.

JK: I am not sure if this class should be MO_matrix type, since it represents
different data (no coefficients, just links MO -> basis function), more
appropriate might be General_MO_matrix, if it is possible (would need
one more allocate in MO_matrix.h and likely other changes as well).
*/

//--------------------------------------------------------------------------

class MO_matrix_basfunc: public MO_matrix
{
protected:
  void init();
private:
  Array2 <int> basisMO;
  //!< basisMO replaces moCoeff, links given MO to its basis function
  Array1 <int> eq_centers;
  Array1 < Array1 <int> > moLists;
  /*!< \brief moLists(spin,.) lists MO's corresponding to the given spin,
       the same thing as totoccupation in Slat_wf_data
  */
  Array1 <doublevar> kptfac;
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

  virtual void writeorb(ostream &, Array2 <doublevar> & rotation, 
                        Array1 <int> &) {
    error("BASFUNC_MO: writeorb not implemented");
  }

  // finally, the two key functions of the class that actually evaluate the
  // molecular orbitals (and their derivatives)
  virtual void updateVal(Sample_point * sample,
			 int e,
			 //!< electron number
			 int listnum,
			 //!< Choose the list that was built in buildLists
			 Array2 <doublevar> & newvals
			 //!< The return: in form (MO)
			 );
  
  virtual void updateLap(Sample_point * sample,
			 int e,
			 //!< electron number
			 int listnum,
			 //!< Choose the list that was built in buildLists
			 Array2 <doublevar> & newvals
			 //!< The return: in form ([value gradient lap], MO)
			 );

  virtual void updateHessian(Sample_point * sample,
			     int e,
			     //!< electron number
			     int listnum,
			     //!< Choose the list that was built in buildLists
			     Array2 <doublevar>& newvals
			     //!< in form ([value gradient, dxx,dyy,dzz,dxy,dxz,dyz], MO)
			     );


  MO_matrix_basfunc()
  {}

};


#endif // MO_MATRIX_BASFUNC_H_INCLUDED

//--------------------------------------------------------------------------
