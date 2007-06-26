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

#ifndef EXPONENT_CUSP_H_INCLUDED
#define EXPONENT_CUSP_H_INCLUDED

#include "Basis_function.h"

/*!
 
*/
class Exponent_cusp: public Basis_function
{
public:


  //-----------------------------------------------------


  virtual int read(
    vector <string> & words,
    //!< The words from the basis section that will create this basis function
    unsigned int & pos

  );




  int nfunc();
  virtual string label()
  {
    return centername;
  }
  doublevar cutoff(int )
  {
    return rcut;
  }

  int showinfo(string & indent, ostream & os);
  int writeinput(string &, ostream &);

  void raw_input(ifstream & input);

  void calcVal(const Array1 <doublevar> & r,
               Array1 <doublevar> & symvals,
               const int startfill=0);

  void calcLap(
    const Array1 <doublevar> & r,
    Array2 <doublevar> & symvals,
    const int startfill=0
  );

  virtual void getVarParms(Array1 <doublevar> & parms);
  virtual void setVarParms(Array1 <doublevar> & parms);
  virtual int nparms() {
    return 1;
  }

private:

  doublevar rcut;
  string centername;
  doublevar gamma;
  doublevar cusp;
};

#endif // EXPONENT_CUSP_H_INCLUDED
//--------------------------------------------------------------------------
