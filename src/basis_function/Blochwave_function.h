/*
 
Copyright (C) 2008 Jindrich Kolorenc

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

#ifndef BLOCHWAVE_FUNCTION_H_INCLUDED
#define BLOCHWAVE_FUNCTION_H_INCLUDED

#include "CBasis_function.h"

/*!

Provides exp(igr) for a number of g's.
 
*/
class Blochwave_function: public CBasis_function
{
public:

  Blochwave_function()
  {}
  ~Blochwave_function()
  {}

  //-----------------------------------------------------


  virtual int read(
    vector <string> & words,
    //!< The words from the basis section that will create this basis function
    unsigned int & pos
    //!< The current position in the words(important if one basis section makes several functions); will be incremented as the Basis_function reads the words.
  );




  int nfunc();
  virtual string label()
  {
    return centername;
  }
  doublevar cutoff(int )
  {
    return 1e99;
  }

  int showinfo(string & indent, ostream & os);
  int writeinput(string &, ostream &);

  void calcVal(const Array1 <doublevar> & r,
               Array1 <dcomplex> & symvals,
               const int startfill=0);

  void calcLap(
    const Array1 <doublevar> & r,
    Array2 <dcomplex> & symvals,
    const int startfill=0
  );

  virtual void calcHessian(
			   const Array1 <doublevar> & r,
			   Array2 <dcomplex> & symvals,
			   const int startfill=0
			   ) { 
    error("This function doesn't support Hessians");
  };
  

  virtual void getVarParms(Array1 <doublevar> & parms);
  virtual void setVarParms(Array1 <doublevar> & parms);
  virtual int nparms() {
    return 0;
  }

private:

  Array2 <doublevar> g_vector;
  Array2 <doublevar> k_vector;
  Array3 <dcomplex> ckg; // complex Ck(g) coefficients
  int nmax; // number of G vectors
  int kmax; // number of k vectors
  int bn_total; // total number of Bloch waves
  vector <int>  bmax; // number of bands for all k points
  string centername;
};

#endif // BLOCHWAVE_FUNCTION_H_INCLUDED
//--------------------------------------------------------------------------
