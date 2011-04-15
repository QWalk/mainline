/*
 
Copyright (C) 2011 Lucas K. Wagner (based on work by Michal Bajdich)

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

#if (!defined(MO_MATRIX_CEINSPLINE_H_INCLUDED) && defined(COMPLEX_WF))|| (!defined(MO_MATRIX_EINSPLINE_H_INCLUDED) && !defined(COMPLEX_WF)) 


#ifdef COMPLEX_WF
#define MO_MATRIX_CEINSPLINE_H_INCLUDED
#else 
#define MO_MATRIX_EINSPLINE_H_INCLUDED
#endif 

#include "MO_matrix.h"
#include "Array45.h"
// USE_RESTRICT with -std=gnu++98 flag for g++; needed for bspline.h 
#ifdef USE_RESTRICT
#define restrict __restrict__
#else
#define restrict 
#endif

// USE_EINSPLINE to add einspline library; tested on einspline-0.8.2
#ifdef USE_EINSPLINE 
#include <bspline.h>
#include <multi_bspline.h>
#endif


/*!
Represents a periodic set of orbitals using Ken Esler's EINSPLINE library. 
 */
#ifdef COMPLEX_WF
class MO_matrix_Ceinspline:public Complex_MO_matrix { 
#else 
class MO_matrix_einspline: public MO_matrix {
#endif
protected:
  void init() { }
private:
#ifdef USE_EINSPLINE
#ifdef COMPLEX_WF
  Array1 <multi_UBspline_3d_z *> spline;
#else
  Array1 <multi_UBspline_3d_d *> spline;
#endif
#endif 
  Array2 <doublevar> latvec;
  Array2 <doublevar> latvecinv;
  Array1 <int> npoints;
  Array1 <doublevar> resolution;
  Array1 <doublevar> kpoint;
#ifdef COMPLEX_WF
  Array4 <dcomplex > modata;
#else
  Array4 <doublevar> modata; //indices: mo,x,y,z
#endif
  Array1 <int> nmo_lists;
  int ndim;
public:
  virtual void buildLists(Array1 <Array1 <int> > & occupations);
  virtual void read(vector <string> & words, unsigned int & startpos, 
                    System * sys);
  virtual int showinfo(ostream & os);
  virtual int writeinput(string &, ostream &);
  virtual void writeorb(ostream &, Array2 <doublevar> & rotation, Array1 <int> & tmp) { } 
  virtual void getMoCoeff(Array2 <doublevar> & coeff) {
    error("Einspline_Mo doesn't support optimization yet");
  }
  virtual void setMoCoeff(Array2 <doublevar> & coeff) {
    error("Einspline MO doesn't support optimization yet");
  }
  virtual int nMoCoeff() {
    error("Need to implement MO_matrix_einspline::nMoCoeff()");
    return 0;
  }
#ifdef COMPLEX_WF
  virtual void updateVal(Sample_point *,int e,int listnum,Array2<dcomplex>&);
  virtual void updateLap(Sample_point *,int e, int listnum, Array2<dcomplex>&);
  MO_matrix_Ceinspline() { } 
#else 
  virtual void updateVal(Sample_point *,int e, int listnum, Array2<doublevar>&);
  virtual void getBasisVal(Sample_point *,int e,Array1<doublevar> &) {
    error("Don't support getBasisVal() in MO_matrix_einspline");
  }
  virtual void updateLap(Sample_point *,int e,int listnum,Array2<doublevar> &);
  MO_matrix_einspline() {}
#endif
};


#endif //MO_MATRIX_EINSPLINE_H_INCLUDED

