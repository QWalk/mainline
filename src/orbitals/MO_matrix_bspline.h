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


#ifndef MO_MATRIX_BSPLINE_H_INCLUDED
#define MO_MATRIX_BSPLINE_H_INCLUDED

#include "Array.h"
#include "Qmc_std.h"
#include "Basis_function.h"
#include "Center_set.h"
#include "MO_matrix.h"
#include "jeep_utils.h"
#include "qmc_io.h"

// USE_RESTRICT with -std=gnu++98 flag for g++; needed for bspline.h 
#ifdef USE_RESTRICT
#define restrict __restrict__
#else
#define restrict 
#endif

// USE_EINSPLINE to add einspline library; tested on einspline-0.8.2
#ifdef USE_EINSPLINE 
#include <bspline.h>
#endif






class System;
class Sample_point;
class EinsplineOrb {
 private:
  #ifdef USE_EINSPLINE
  UBspline_3d_d *Spline;
  #endif
 public:
  Array1 <doublevar> origin, box_size, spacing;
  int periodic;
  void evaluate_spline (const Array1 < doublevar > & r, Array1 <doublevar> & newvals) {
    assert( newvals.GetSize()<=5);
    doublevar val;
#ifdef USE_EINSPLINE
    eval_UBspline_3d_d (Spline, r(0), r(1), r(2), &val);
#else
    error("need EINSPLINE library for eval_UBspline_3d_d");
#endif
    newvals(0)=val;
  }

  void evaluate_spline_vgl (const Array1 < doublevar > & r, Array1 <doublevar> & newvals){
    assert( newvals.GetSize()==5);
    doublevar val,lap;
    Array1 <doublevar> grad(3);
#ifdef USE_EINSPLINE
    eval_UBspline_3d_d_vgl (Spline, r(0), r(1), r(2),
			    &val, grad.v, &lap);
#else
    error("need EINSPLINE library for eval_UBspline_3d_d_vgl");
#endif
    newvals(0)=val;
    for(int d=1;d<4;d++)
       newvals(d)=grad(d-1);
    newvals(4)=lap;
  }

  void evaluate_spline_vgh (const Array1 < doublevar > & r, Array1 <doublevar> &  newvals){
    assert( newvals.GetSize()==10);
    doublevar val;
    Array1 <doublevar> grad(3);
    Array2 <doublevar> hessian(3,3);
#ifdef USE_EINSPLINE
    eval_UBspline_3d_d_vgh (Spline, r(0), r(1), r(2), &val, grad.v, hessian.v);
#else
 error("need EINSPLINE library for eval_UBspline_3d_d_vgh");
#endif
    newvals(0)=val;
     for(int d=1;d<4;d++)
       newvals(d)=grad(d-1);
     int d=4;
     for(int i=0;i<3;i++)
       for(int j=i;j<3;j++)
	 newvals(d++)=hessian(i,j);
  }

  void read (ifstream & ORB_0){
    single_write(cout, "reading in plotfile\n");
    Array3 <doublevar> grid;
    get_grid_from_plot(ORB_0, grid, box_size, origin);
#ifdef USE_EINSPLINE
    BCtype_d xBC, yBC, zBC;
    if(periodic){
      xBC.lCode = xBC.rCode = PERIODIC;
      yBC.lCode = yBC.rCode = PERIODIC;
      zBC.lCode = zBC.rCode = PERIODIC;
    }
    else{
      xBC.lCode = xBC.rCode = NATURAL;
      yBC.lCode = yBC.rCode = NATURAL;
      zBC.lCode = zBC.rCode = NATURAL;
    }

    Ugrid x_grid, y_grid, z_grid;
    x_grid.start = origin(0); x_grid.end = origin(0)+box_size(0); x_grid.num =grid.GetDim(0) ;
    y_grid.start = origin(1); y_grid.end = origin(1)+box_size(1); y_grid.num =grid.GetDim(1) ;
    z_grid.start = origin(2); z_grid.end = origin(2)+box_size(2); z_grid.num =grid.GetDim(2) ;
    
    spacing.Resize(3); 
    if(periodic){
      for(int d=0;d<3;d++)
	spacing(d)=box_size(d)/(grid.GetDim(d));
    }
    else
      for(int d=0;d<3;d++)
	spacing(d)=box_size(d)/(grid.GetDim(d)-1); 

    single_write(cout, "creating spline\n");
    Spline=create_UBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC, grid.v);
#else
 error("need EINSPLINE library for create_UBspline_3d_d");
#endif
    grid.clear() ;
  }
  
 ~EinsplineOrb() 
   {}
};

//----------------------------------------------------------------------------

class MO_matrix_bspline: public MO_matrix
{
protected:
  void init();
private:
  Array1 < Array1 <int> > moLists;
  Array1 <EinsplineOrb *> Einspline;
  vector <string> valfiles;
 
public:
  int periodic;
  virtual void buildLists(Array1 <Array1 <int> > & occupations);

  virtual int showinfo(ostream & os);

  virtual int writeinput(string &, ostream &);
  virtual void read(vector <string> & words, unsigned int & startpos, System * sys);
  virtual void writeorb(ostream &, Array2 <doublevar> & rotation, 
                        Array1 <int> &) {
    error("BSPLINE_MO: writeorb not implemented");
  }

  // I guess the following three are ment for direct orbital optimization,
  // there is nothing to optimize here, only basis functions
  virtual void getMoCoeff(Array2 <doublevar> & coeff) {
    error("BSPLINE_MO: getMoCoeff not implemented");
  }
  virtual void setMoCoeff(Array2 <doublevar> & coeff) {
    error("BSPLINE_MO: setMoCoeff not implemented");
  }
  virtual int nMoCoeff() {
    error("BSPLINE_MO: nMoCoeff not implemented");
    return 0;
  }

  // No problem to implement/copy, but what for is it
  // anyway?
  virtual void getBasisVal(Sample_point * sample, 
			   int e,
			   Array1 <doublevar> & newvals
			   ) {
    error("BSPLINE_MO: getBasisVal not implemented");
  }

  // finally, the three key functions of the class that actually evaluate the
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


  MO_matrix_bspline()
  {}

};





#endif // MO_MATRIX_BSPLINE_H_INCLUDED

//--------------------------------------------------------------------------
