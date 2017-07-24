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

#ifndef MO_1D_H_INCLUDED
#define MO_1D_H_INCLUDED

#include "MO_matrix.h"
#include "Spline_fitter.h"

class MO_1d : public Complex_MO_matrix {
 public:
  
  virtual void buildLists(Array1 <Array1 <int> > & occupations);

  virtual int getNmo() { return nmo; }
  virtual int showinfo(ostream & os);
  virtual int writeinput(string &, ostream &);
  virtual void read(vector <string> & words, unsigned int & startpos, 
                    System * sys);

  virtual void updateVal(Sample_point * sample,int e, int listnum,
                         Array2 <dcomplex> & newvals);

  virtual void updateLap(Sample_point * sample,int e,int listnum,
                         Array2 <dcomplex> & newvals);


 private:
  int nmo;  
  Array1 < Array1 <int> > moLists;
  Array2 <Spline_fitter_nonuniform> splines;
  int npoints;
  string readfile;
};

#endif //MO_1D_H_INCLUDED
