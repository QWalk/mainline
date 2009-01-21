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

#ifndef MO_MATRIX_H_INCLUDED
#define MO_MATRIX_H_INCLUDED

#include "Array.h"
#include "Qmc_std.h"
#include "Basis_function.h"
#include "Center_set.h"

class System;
class Sample_point;

class General_MO_matrix {
 public:
  /*!
    Build several sets of MO's to be evaluated in updateVal and updateLap.
    Each element in occupations should be a list of the MO's that should
    be evaluated.  For example, one can create a list of spin up and spin
    down MO's, and only evaluate up when an up electron is moved.
   */
  virtual void buildLists(Array1 <Array1 <int> > & occupations)=0;

  /*!
    get the number of molecular orbitals
   */
  virtual int getNmo()=0;
  virtual int showinfo(ostream & os)=0;
  virtual int writeinput(string &, ostream &)=0;
  virtual void read(vector <string> & words, unsigned int & startpos, 
                    System * sys)=0;

  virtual ~General_MO_matrix()  {  }

};

//----------------------------------------------------------------------

class Complex_MO_matrix : public General_MO_matrix {

 protected:
  // transferred from MO_matrix
  Center_set centers;
  Array1 <CBasis_function *> basis;
  int nmo;
  doublevar magnification_factor;
  string orbfile;
  int totbasis;
  int maxbasis;
  string oldsofile;
  Array1 <doublevar> kpoint;
  //!< the k-point of the orbitals in fractional units(1 0 0) is X, etc..
  // -----

 public:
  /*!
    Build several sets of MO's to be evaluated in updateVal and updateLap.
    Each element in occupations should be a list of the MO's that should
    be evaluated.  For example, one can create a list of spin up and spin
    down MO's, and only evaluate up when an up electron is moved.
   */
  virtual void buildLists(Array1 <Array1 <int> > & occupations)=0;

  /*!
    get the number of molecular orbitals
  */
  int getNmo() { return nmo; }
  virtual int showinfo(ostream & os)=0;
  virtual int writeinput(string &, ostream &)=0;
  virtual void read(vector <string> & words, unsigned int & startpos, 
                    System * sys)=0;
  virtual void updateVal(Sample_point * sample,
			 int e,
			 //!< electron number
			 int listnum,
			 Array2 <dcomplex> & newvals
			 //!< The return: in form (MO)
			 )=0;
  
  virtual void updateLap(Sample_point * sample,
			 int e,
			 //!< electron number
			 int listnum,
			 //!< Choose the list that was built in buildLists
			 Array2 <dcomplex> & newvals
			 //!< The return: in form ([value gradient lap], MO)
			 )=0;

  virtual ~Complex_MO_matrix() {  }

};

//----------------------------------------------------------------------

/*!
\brief
MO_matrix is a class to hold everything necessary to calculate the
molecular orbitals for a situation.

*/
class MO_matrix: public General_MO_matrix
{
protected:
  Center_set centers;
  Array1 <Basis_function *> basis;
  int nmo;
  doublevar magnification_factor;
  string orbfile;
  int totbasis;
  int maxbasis;
  virtual void init()=0;

  string oldsofile;

  Array1 <doublevar> kpoint; //!< the k-point of the orbitals in fractional units(1 0 0) is X, etc..
  //Array1 <doublevar> mo_counter;//!< Count how many basis functions are evaluated per MO
  //unsigned int n_calls; //!< number of calls
public:

  /*!
    Build several sets of MO's to be evaluated in updateVal and updateLap.
    Each element in occupations should be a list of the MO's that should
    be evaluated.  For example, one can create a list of spin up and spin
    down MO's, and only evaluate up when an up electron is moved.
   */
  virtual void buildLists(Array1 <Array1 <int> > & occupations)=0;


  virtual void setOrbfile(string & x) {
    orbfile=x;
  }
  /*!
    get the number of molecular orbitals
   */
  int getNmo() {
    return nmo;
  }

  virtual int showinfo(ostream & os)=0;

  virtual int writeinput(string &, ostream &)=0;


  virtual void read(vector <string> & words, unsigned int & startpos, 
                    System * sys);

  //! Takes an ORB file and inserts all the coefficients.
  // virtual int readorb(istream &, Array3 <int> &, Array1 <doublevar> & );


  virtual void writeorb(ostream &, Array2 <doublevar> & rotation, Array1 <int> &)=0;

  /*!
    Get the molecular orbital coefficients
   */
  virtual void getMoCoeff(Array2 <doublevar> & coeff)=0;
  
  /*!
    
   */
  virtual void setMoCoeff(Array2 <doublevar> & coeff)=0;
  virtual int nMoCoeff()=0;


  virtual void updateVal(
    Sample_point * sample,
    int e,
    //!< electron number
    int listnum,
    Array2 <doublevar> & newvals
    //!< The return: in form (MO)
  )=0;

  virtual void getBasisVal(
    Sample_point * sample,
    int e,
    //!< electron number
    Array1 <doublevar> & newvals
    )=0;

  virtual void updateLap(
    Sample_point * sample,
    int e,
    //!< electron number
    int listnum,
    //!< Choose the list that was built in buildLists
    Array2 <doublevar> & newvals
    //!< The return: in form (MO,[value gradient lap])
  )=0;

  virtual void updateHessian(Sample_point * sample,
			     int e,
			     int listnum,
			     Array2 <doublevar>& newvals
			     //!< in form (MO, [value gradient, dxx,dyy,dzz,dxy,dxz,dyz])
			     ) { 
    error("this MO_matrix doesn't support Hessians");
  }

  MO_matrix()
  {}

  virtual ~MO_matrix()
  {
    //doublevar totcalls=0;
    //for(int i=0; i< nmo; i++) {
    //  cout << "mo_counter " << mo_counter(i)/n_calls << endl;
    //  totcalls+=mo_counter(i)/n_calls;
    //}
    //cout << " average # of basis functions evaluated " << totcalls/nmo << endl;
    for(int i=0; i< basis.GetDim(0); i++)
      deallocate(basis(i));
  }

};

//----------------------------------------------------------------------------


int allocate(vector <string> & words, System * sys, MO_matrix *& moptr);
int allocate(vector <string> & words, System * sys, 
             Complex_MO_matrix *& moptr);

void rotate_orb(istream & orbin, ostream & orbout,
                Array2 <doublevar> & rotation,
                Array1 <int>  & moList, int nfunctions);
void rotate_Corb(istream & orbin, ostream & orbout,
		 Array2 <doublevar> & rotation,
		 Array1 <int>  & moList, int nfunctions);


int readorb(istream & input, Center_set & centers, 
            int nmo, int maxbasis, Array1 <doublevar> & kpoint, 
	    Array3 <int > & coeffmat, Array1 <doublevar> & coeff);
int readorb(istream & input, Center_set & centers, 
            int nmo, int maxbasis, Array1 <doublevar> & kpoint, 
	    Array3 <int > & coeffmat, Array1 <dcomplex> & coeff);

#endif // MO_MATRIX_H_INCLUDED

//----------------------------------------------------------------------------
