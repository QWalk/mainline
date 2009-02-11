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


#ifndef MO_MATRIX_CBSPLINE_H_INCLUDED
#define MO_MATRIX_CBSPLINE_H_INCLUDED

#include "Array.h"
#include "Qmc_std.h"
#include "Basis_function.h"
#include "Center_set.h"
#include "MO_matrix.h"
#include "jeep_utils.h"
#include "qmc_io.h"
#include "MatrixAlgebra.h"
#include "MO_matrix_bspline.h"

class System;
class Sample_point;
class MO_matrix_Cbspline: public Complex_MO_matrix
{
protected:
  void init();
private:
  Array1 < Array1 <int> > moLists;
  MultiEinsplineOrb* MultiEinspline;
  MultiEinsplineOrbComplex* MultiEinsplineComplex;

  vector <string> bandfiles;
  Array1  <doublevar> kpoint;
  doublevar kpoint_square;
  Array1 < Array1 < doublevar> > kvectors_linear;
  Array1 < Array1 < doublevar> > kvectors;
  vector < vector < vector <string> > > bands_per_kvectors;
  Array1  <doublevar> origin;
  Array2  <doublevar> RecipLatVec;
  Array2  <doublevar> LatVec;
  Array2  <doublevar> PrimRecipLatVec;
  Array2  <doublevar> PrimLatVec;
  int nsplines; //number of splines

public:
  virtual void buildLists(Array1 <Array1 <int> > & occupations);
  virtual int showinfo(ostream & os);
  virtual int writeinput(string &, ostream &);
  virtual void read(vector <string> & words, unsigned int & startpos, System * sys);
  virtual void writeorb(ostream &, Array2 <doublevar> & rotation, 
                        Array1 <int> &) {
    error("CBSPLINE_MO: writeorb not implemented");
  }

  // I guess the following three are ment for direct orbital optimization,
  // there is nothing to optimize here, only basis functions
  virtual void getMoCoeff(Array2 <doublevar> & coeff) {
    error("CBSPLINE_MO: getMoCoeff not implemented");
  }
  virtual void setMoCoeff(Array2 <doublevar> & coeff) {
    error("CBSPLINE_MO: setMoCoeff not implemented");
  }
  virtual int nMoCoeff() {
    error("CBSPLINE_MO: nMoCoeff not implemented");
    return 0;
  }

  // No problem to implement/copy, but what for is it
  // anyway?
  virtual void getBasisVal(Sample_point * sample, 
			   int e,
			   Array1 <dcomplex> & newvals
			   ) {
    error("CBSPLINE_MO: getBasisVal not implemented");
  }

  // finally, the three key functions of the class that actually evaluate the
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
			     );


  MO_matrix_Cbspline()
  {}

};





#endif // MO_MATRIX_CBSPLINE_H_INCLUDED

//--------------------------------------------------------------------------
