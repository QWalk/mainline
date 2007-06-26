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

#ifndef MO_MATRIX_CUTOFF_H_INCLUDED
#define MO_MATRIX_CUTOFF_H_INCLUDED

#include "Array.h"
#include "Qmc_std.h"
#include "Basis_function.h"
#include "Center_set.h"
#include "MO_matrix.h"

class System;
class Sample_point;
//----------------------------------------------------------------------------

class MO_matrix_cutoff: public MO_matrix
{
protected:
  void init();
private:
  //Center_set centers;
  //Array1 <Basis_function *> basis;
  //int nmo;
  //int totbasis;
  //int maxbasis;
 // doublevar magnification_factor;
  //string orbfile;
  Array2 <int> mofill;
  Array2 <doublevar> moCoeff2;
  Array1 <int> nbasis;

  Array1 <doublevar> obj_cutoff; //!< cutoff for each basis object
  Array1 <doublevar> cutoff;  //!< Cutoff for individual basis functions
  Array1 <int> nfunctions; //!< number of functions in each basis
  //Array1 <int> basismo;
  //Array2 <doublevar> moCoeff;
  //Array2 <int> basisfill;

  Array1 < Array2 <int> > basisfill_list;
  Array1 < Array2 <doublevar> > moCoeff_list;
  Array1 < Array1 <int> > basismo_list;



public:

  /*!
    Build several sets of MO's to be evaluated in updateVal and updateLap.
    Each element in occupations should be a list of the MO's that should
    be evaluated.  For example, one can create a list of spin up and spin
    down MO's, and only evaluate up when an up electron is moved.
   */
  virtual void buildLists(Array1 <Array1 <int> > & occupations);

  /*!
    get the number of molecular orbitals
   */
  //virtual int getNmo()
  //{
  //  return nmo;
  //}

  virtual int showinfo(ostream & os);

  virtual int writeinput(string &, ostream &);


  //virtual void read(vector <string> & words, unsigned int & startpos, System * sys);

  //! Takes an ORB file and inserts all the coefficients.
  //virtual int readorb(istream &);


  virtual void writeorb(ostream &, Array2 <doublevar> & rotation, Array1 <int> &);

  virtual void getMoCoeff(Array2 <doublevar> & coeff) {
    error("Cutoff_Mo doesn't support optimization yet");
  }
  
  /*!
    
   */
  virtual void setMoCoeff(Array2 <doublevar> & coeff) {
    error("Cutoff MO doesn't support optimization yet");
  }

  virtual int nMoCoeff() {
    error("Need to implement MO_matrix_cutoff::nMoCoeff()");
    return 0;
  }


  virtual void updateVal(
    Sample_point * sample,
    int e,
    //!< electron number
    int listnum,
    Array2 <doublevar> & newvals
    //!< The return: in form (MO)
  );
  
  virtual void getBasisVal(
    Sample_point * sample,
    int e,
    Array1 <doublevar> & newvals
    ){
    error("Need to implement MO_matrix_cutoff::getBasisVal()");
  }

  virtual void updateLap(
    Sample_point * sample,
    int e,
    int listnum,
    Array2 <doublevar> & newvals
  );
  virtual void updateHessian(Sample_point * sample,
			     int e,
			     int listnum,
			     Array2 <doublevar>& newvals
			     //!< in form ([value gradient, dxx,dyy,dzz,dxy,dxz,dyz], MO)
			     );
  MO_matrix_cutoff()
  {}

};



#endif // MO_MATRIX_CUTOFF_H_INCLUDED

//--------------------------------------------------------------------------
