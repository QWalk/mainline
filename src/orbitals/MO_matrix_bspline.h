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
#include "MatrixAlgebra.h"

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






class System;
class Sample_point;

//----------------------------------------------------------------------------
class MultiEinsplineOrb {
 protected:
  //Array2 <doublevar> RecipLatVec, LatVec;
  Array1 < Array1 <doublevar> > kpointt;
 public:
  virtual void kpoint (Array1 < Array1 < doublevar > > & kvector){
    kpointt.Resize(kvector.GetDim(0));
    for(int i=0;i<kvector.GetDim(0);i++){
      kpointt(i).Resize(3);
      kpointt(i)=kvector(i);
    }
  }
  virtual void evaluate_spline (const Array1 < doublevar > & r, Array1 <dcomplex> & val)
  { error ("evaluate_spline not defined!");};
  virtual void evaluate_spline  (const Array1 < doublevar > & r, Array1 <dcomplex> & val, 
				 Array1 < Array1 <dcomplex> > & grad, Array1 <dcomplex> & lap)
  {error ("evaluate_spline not defined!");};
  virtual void evaluate_spline  (const Array1 < doublevar > & r, Array1 <dcomplex> & val, 
				 Array1 < Array1 <dcomplex> > & grad, Array1 < Array2 <dcomplex> > & hess)
  {error ("evaluate_spline not defined!");};
  virtual void read (vector <string> & valfiles, int & nsplines)=0;
  virtual ~MultiEinsplineOrb() { };
};
//-------------------------------------------------------------------

class MultiEinsplineOrbComplex : public MultiEinsplineOrb {
   private:
#ifdef USE_EINSPLINE
  multi_UBspline_3d_z *MultiSpline; 
#endif
 public:
  virtual void evaluate_spline (const Array1 < doublevar > & r, Array1 < dcomplex > & val) {
    val=dcomplex(0.0,0.0);
    for(int d=0;d<3;d++){
      if((r(d)<0) || (r(d)>1.0)){
	cout  <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	cout  <<" is outside the box, orbital value is set to zero"<<endl;
	return;
      }
    }

#ifdef USE_EINSPLINE
    eval_multi_UBspline_3d_z (MultiSpline, r(0), r(1), r(2), val.v);
#else
    error("need EINSPLINE library for eval_UBspline_3d_d");
#endif
  }

  virtual void evaluate_spline (const Array1 < doublevar > & r, Array1 <dcomplex> & val, 
				Array1 < Array1 <dcomplex> >& grad, Array1 <dcomplex> & lap){
    val=lap=dcomplex(0.0,0.0);
    for(int m=0;m<grad.GetSize();m++){
      grad(m)=dcomplex(0.0,0.0);
    }
    for(int d=0;d<3;d++){
      if((r(d)<0) || (r(d)>1.0)){
	cout  <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	cout  <<" is outside the box, orbitals val, grad and lap are set to zero"<<endl;
	return;
	}
    }
   
     Array1 <dcomplex> linear_grad(3*grad.GetSize());
    
#ifdef USE_EINSPLINE
    eval_multi_UBspline_3d_z_vgl (MultiSpline, r(0), r(1), r(2),
			    val.v, linear_grad.v, lap.v);

    for(int m=0;m<grad.GetSize();m++){
      for(int i=0;i<3;i++){
	grad(m)(i)=linear_grad(3*m+i);
      }
    }

#else
    error("need EINSPLINE library for eval_UBspline_3d_d_vgl");
#endif
  }

  virtual void evaluate_spline  (const Array1 < doublevar > & r, Array1 <dcomplex> & val, 
				 Array1 <Array1 <dcomplex> > & grad, 
				 Array1 <Array2 <dcomplex> > & hess){
    val=dcomplex(0.0,0.0);
    for(int m=0;m<grad.GetSize();m++){
      grad(m)=dcomplex(0.0,0.0);
      hess(m)=dcomplex(0.0,0.0);
    }
    for(int d=0;d<3;d++){
      if((r(d)<0) || (r(d)>1.0)){
	cout  <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	cout  <<" is outside the box, orbitals val, grad and hess are set to zero"<<endl;
	return;
      }
    }

     Array1 <dcomplex> linear_grad(3*grad.GetSize());
     Array1 <dcomplex> linear_hess(9*hess.GetSize());
    
#ifdef USE_EINSPLINE
    eval_multi_UBspline_3d_z_vgh (MultiSpline, r(0), r(1), r(2), val.v, linear_grad.v, linear_hess.v);

    for(int m=0;m<grad.GetSize();m++){
      for(int i=0;i<3;i++){
	grad(m)(i)=linear_grad(3*m+i);
	for(int j=0;j<3;j++){
	  hess(m)(i,j)=linear_hess(9*m+3*i+j);
	}
      }
    }
#else
 error("need EINSPLINE library for eval_UBspline_3d_d_vgh");
#endif
  }

 
  virtual void read (vector <string> & valfiles, int & nsplines){
    Array1 <doublevar> resolution_array(3);
    Array1 <doublevar> box_size(3);
    Array1 <doublevar> origin(3);
#ifdef USE_EINSPLINE
    for(int m=0; m < nsplines; m++) {
      //read a real part
      string orbfile_0=valfiles[2*m];
      ifstream ORB_0(orbfile_0.c_str());
      if(!ORB_0){
	error("couldn't find orb file ", orbfile_0);
      } 
      single_write(cout,"Reading Orbfile: ",orbfile_0,"\n"); 
      Array3 <doublevar> grid_real;
      get_grid_from_plot(ORB_0, grid_real, box_size, origin);

      //read a imag part
      string orbfile_1=valfiles[2*m+1];
      ifstream ORB_1(orbfile_1.c_str());
      if(!ORB_1){
	error("couldn't find orb file ", orbfile_1);
      } 
      single_write(cout,"Reading Orbfile: ",orbfile_1,"\n"); 
      Array3 <doublevar> grid_imag;
      get_grid_from_plot(ORB_1, grid_imag, box_size, origin);

      //check the size of imag wrt. to real grid
      for (int d=0;d<3;d++)
	if(grid_real.GetDim(d)!=grid_imag.GetDim(d))
	  error("Real and Imag grids have diffrent sizes");

      if(m==0){
	BCtype_z xBC, yBC, zBC;
	Ugrid x_grid, y_grid, z_grid;
	xBC.lCode = xBC.rCode = PERIODIC;
	yBC.lCode = yBC.rCode = PERIODIC;
	zBC.lCode = zBC.rCode = PERIODIC;
	  
	//spline calculation is in fractional coordinates of lattice vectors
	x_grid.start = 0.0; x_grid.end = 1.0; x_grid.num =grid_real.GetDim(0) ;
	y_grid.start = 0.0; y_grid.end = 1.0; y_grid.num =grid_real.GetDim(1) ;
	z_grid.start = 0.0; z_grid.end = 1.0; z_grid.num =grid_real.GetDim(2) ;
	  
	for(int d=0;d<3;d++){
	  resolution_array(d)=1.0/grid_real.GetDim(d);
	}
	//make complex grid by multiplying by exp(-2*i*pi*k*r)
	//be aware that kvector in the code includes factor 2!
	MultiSpline = create_multi_UBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC, nsplines);
      }//m==0
      else{
	if(MultiSpline->x_grid.num!=grid_real.GetDim(0) || 
	   MultiSpline->y_grid.num!=grid_real.GetDim(1) ||
	   MultiSpline->z_grid.num!=grid_real.GetDim(2))
	  error("Number of grid points is different that in the first plotfile");
      }
	
      Array3 <dcomplex> cgrid(MultiSpline->x_grid.num,MultiSpline->y_grid.num,MultiSpline->z_grid.num);
      Array1 <doublevar> xyz(3);  //fractional coordinates of Lattice vectors
      for(int xx=0;xx<grid_real.GetDim(0);xx++){
	xyz(0)=xx*resolution_array(0);
	for(int yy=0;yy<grid_real.GetDim(1);yy++){
	  xyz(1)=yy*resolution_array(1);
	  for(int zz=0;zz<grid_real.GetDim(2);zz++){
	    xyz(2)=zz*resolution_array(2);
	    //cout <<xyz(0)<<"  "<<xyz(1)<<"  "<<xyz(2)<<endl;
	    dcomplex eikr=exp(-I*dot(xyz, kpointt(m)));
	    cgrid(xx,yy,zz)=(grid_real(xx,yy,zz)+I*grid_imag(xx,yy,zz))*eikr;
	    //if(xx==10 && yy==15)
	    //cout << xyz(2) <<"  "<<grid_real(xx,yy,zz)<<"  "<<grid_imag(xx,yy,zz)
	    //  <<"  "<<cgrid(xx,yy,zz).real()<<"  "<<cgrid(xx,yy,zz).imag()<<"  "<<eikr.real()<<"  "<<eikr.imag()<<endl;
	  } 
	}
      }
      grid_real.clear();
      grid_imag.clear();
      single_write(cout, "creating spline\n");
      set_multi_UBspline_3d_z (MultiSpline, m, cgrid.v);
      single_write(cout, "done creating spline\n");
      cgrid.clear();
    }//m
#else
 error("need EINSPLINE library for create_UBspline_3d_d");
#endif
  }
  virtual ~MultiEinsplineOrbComplex() { };
};

//----------------------------------------------------------------------------

class MO_matrix_bspline: public MO_matrix
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
