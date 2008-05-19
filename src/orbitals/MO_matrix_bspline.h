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
class EinsplineOrb {
 protected:
  int periodic;
  Array1 <doublevar> origin, box_size, spacing;
  Array2 <doublevar> RecipLatVec, LatVec;
  Array1 <doublevar> kpointt;
 public:
  virtual void getorigin(Array1 <doublevar> & a){
    a.Resize(3);
    a=origin;
  }
  virtual void getbox_size(Array1 <doublevar> & a){
    a.Resize(3);
    a=box_size;
  }
  virtual void getspacing(Array1 <doublevar> & a){
    a.Resize(3);
    a=spacing;
  }
  virtual void setperiodic(int & i){
    periodic=i;
  }
  virtual void kpoint (Array1 < doublevar > & kvector){
    kpointt.Resize(3);
    kpointt=kvector;
  }
  virtual void setorigin (Array1 < doublevar > & o){
    origin.Resize(3);
    origin=o;
  }
  virtual void setRecipLattice(Array2 <doublevar> & a){
    RecipLatVec.Resize(3,3);
    RecipLatVec=a;
  }
  virtual void setBounds(Array2 <doublevar> & a){
    LatVec.Resize(3,3);
    LatVec=a;
  }
  virtual void evaluate_spline (const Array1 < doublevar > & r, dcomplex & val)
  { error ("evaluate_spline not defined!");};
  virtual void evaluate_spline (const Array1 < doublevar > & r, dcomplex & val, Array1 <dcomplex> & grad, dcomplex & lap)
  {error ("evaluate_spline not defined!");};
  virtual void evaluate_spline (const Array1 < doublevar > & r, dcomplex & val, Array1 <dcomplex> & grad, Array2 <dcomplex> & hess)
  {error ("evaluate_spline not defined!");};
  virtual void evaluate_spline (const Array1 < doublevar > & r, doublevar & val)
  {error ("evaluate_spline not defined!");};
  virtual void evaluate_spline (const Array1 < doublevar > & r, doublevar & val, Array1 <doublevar> & grad, doublevar & lap)
  {error ("evaluate_spline not defined!");};
  virtual void evaluate_spline (const Array1 < doublevar > & r, doublevar & val, Array1 <doublevar> & grad, Array2 <doublevar> & hess)
  {error ("evaluate_spline not defined!");};
  virtual void read (ifstream & ORB_0)=0;
  virtual ~EinsplineOrb() { };
};


//-------------------------------------------------------------------

class EinsplineOrbReal : public EinsplineOrb {
 private:
  #ifdef USE_EINSPLINE
  UBspline_3d_d *Spline;
  #endif
 public:
  virtual void evaluate_spline (const Array1 < doublevar > & r, doublevar & val) {
    val=0;
    for(int d=0;d<3;d++){
      if(!periodic){
	if((r(d)<origin(d)) || (r(d)>origin(d)+box_size(d))){
	  cout  <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout <<" is outside the box, orbitals value is set zero"<<endl;
	  return;
	}
      }
      else{
	if((r(d)<0) || (r(d)>1.0)){
	  cout  <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout <<" is outside the box, orbitals value is set zero"<<endl;
	  return;
	}
      }
    }

#ifdef USE_EINSPLINE
    eval_UBspline_3d_d (Spline, r(0), r(1), r(2), &val);
#else
    error("need EINSPLINE library for eval_UBspline_3d_d");
#endif
  }

  virtual void evaluate_spline (const Array1 < doublevar > & r, doublevar & val, Array1 <doublevar> & grad, doublevar & lap){ 
    val=lap=0;
    grad=0;
    for(int d=0;d<3;d++){
      if(!periodic){
	if((r(d)<origin(d)) || (r(d)>origin(d)+box_size(d)) ){
	  cout <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout <<" is outside the box, orbitals val, grad and lap are set to zero"<<endl;
	  return;
	}
      }
      else{
	if((r(d)<0) || (r(d)>1.0)){
	  cout  <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout  <<" is outside the box, orbitals val, grad and lap are set to zero"<<endl;
	  return;
	}
      }
    }
   
#ifdef USE_EINSPLINE
    eval_UBspline_3d_d_vgl (Spline, r(0), r(1), r(2),
			    &val, grad.v, &lap);
#else
    error("need EINSPLINE library for eval_UBspline_3d_d_vgl");
#endif
  }

  virtual void evaluate_spline (const Array1 < doublevar > & r, doublevar & val, Array1 <doublevar> & grad, Array2 <doublevar> & hess){
    val=0;
    grad=0;
    hess=0;
    for(int d=0;d<3;d++){
      if(!periodic){
	if((r(d)<origin(d)) || (r(d)>origin(d)+box_size(d)) ){
	  cout <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout <<" is outside the box, orbitals val, grad and hess are set to zero"<<endl;
	  return;
	}
      }
      else{
	if((r(d)<0) || (r(d)>1.0)){
	  cout  <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout  <<" is outside the box, orbitals val, grad and hess are set to zero"<<endl;
	  return;
	}
      }
    }
    
#ifdef USE_EINSPLINE
    eval_UBspline_3d_d_vgh (Spline, r(0), r(1), r(2), &val, grad.v, hess.v);
#else
 error("need EINSPLINE library for eval_UBspline_3d_d_vgh");
#endif
  }

  virtual void read (ifstream & ORB_0){
    single_write(cout, "reading in plotfile\n");
    Array3 <doublevar> grid;
    get_grid_from_plot(ORB_0, grid, box_size, origin);
#ifdef USE_EINSPLINE
    BCtype_d xBC, yBC, zBC;
    Ugrid x_grid, y_grid, z_grid;
    if(periodic){
      xBC.lCode = xBC.rCode = PERIODIC;
      yBC.lCode = yBC.rCode = PERIODIC;
      zBC.lCode = zBC.rCode = PERIODIC;

      //spline calculation is in fractional coordinates of lattice vectors
      x_grid.start = 0.0; x_grid.end = 1.0; x_grid.num =grid.GetDim(0) ;
      y_grid.start = 0.0; y_grid.end = 1.0; y_grid.num =grid.GetDim(1) ;
      z_grid.start = 0.0; z_grid.end = 1.0; z_grid.num =grid.GetDim(2) ;

      //overwrite box_size from get_grid_from_plot
      for(int d=0;d<3;d++)
	box_size(d)=LatVec(0,d)+LatVec(1,d)+LatVec(2,d);
    }
    else{
      xBC.lCode = xBC.rCode = NATURAL;
      yBC.lCode = yBC.rCode = NATURAL;
      zBC.lCode = zBC.rCode = NATURAL;
      
      x_grid.start = origin(0); x_grid.end = origin(0)+box_size(0); x_grid.num =grid.GetDim(0) ;
      y_grid.start = origin(1); y_grid.end = origin(1)+box_size(1); y_grid.num =grid.GetDim(1) ;
      z_grid.start = origin(2); z_grid.end = origin(2)+box_size(2); z_grid.num =grid.GetDim(2) ;
    }
    single_write(cout, "creating spline\n");
    Spline=create_UBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC, grid.v);
    single_write(cout, "done creating spline\n");
    grid.clear() ;
   
    //getting the spacing from the Spline object
    spacing.Resize(3); 

    if(!periodic){
      spacing(0)=Spline->x_grid.delta;
      spacing(1)=Spline->y_grid.delta;
      spacing(2)=Spline->z_grid.delta;

    }
    else{
      //adjusted by lenght of each latice vector;    
      Array1 <doublevar> LatVecNorm(3);
      LatVecNorm=0;
      for(int i=0;i<3;i++){
	for(int j=0;j<3;j++)
	  LatVecNorm(i)+=LatVec(i,j)*LatVec(i,j);	
	LatVecNorm(i)=sqrt(LatVecNorm(i));
      }
      spacing(0)=LatVecNorm(0)*Spline->x_grid.delta;
      spacing(1)=LatVecNorm(1)*Spline->y_grid.delta;
      spacing(2)=LatVecNorm(2)*Spline->z_grid.delta;
    }

#else
    error("need EINSPLINE library for create_UBspline_3d_d");
#endif
  }
  virtual ~EinsplineOrbReal() { };
};


//-------------------------------------------------------------------



class EinsplineOrbComplex : public EinsplineOrb {
 private:
  #ifdef USE_EINSPLINE
  UBspline_3d_z *Spline; 
  #endif
 public:
  virtual void evaluate_spline (const Array1 < doublevar > & r, dcomplex & val) {
    val=dcomplex(0.0,0.0);
    for(int d=0;d<3;d++){
      if(!periodic){
	if((r(d)<origin(d)) || (r(d)>origin(d)+box_size(d)) ){
	  cout <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout <<" is outside the box, orbital value is set to zero"<<endl;
	  return;
	}
      }
      else{
	if((r(d)<0) || (r(d)>1.0)){
	  cout  <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout  <<" is outside the box, orbital value is set to zero"<<endl;
	  return;
	}
      }
    }

#ifdef USE_EINSPLINE
    eval_UBspline_3d_z (Spline, r(0), r(1), r(2), &val);
#else
    error("need EINSPLINE library for eval_UBspline_3d_d");
#endif
  }

  virtual void evaluate_spline (const Array1 < doublevar > & r, dcomplex & val, Array1 <dcomplex> & grad, dcomplex & lap){
    val=lap=0;
    grad=dcomplex(0.0,0.0);
    for(int d=0;d<3;d++){
      if(!periodic){
	if((r(d)<origin(d)) || (r(d)>origin(d)+box_size(d)) ){
	  cout <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout <<" is outside the box, orbitals val, grad and lap are set to zero"<<endl;
	  return;
	}
      }
      else{
	if((r(d)<0) || (r(d)>1.0)){
	  cout  <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout  <<" is outside the box, orbitals val, grad and lap are set to zero"<<endl;
	  return;
	}
      }
    }
   
#ifdef USE_EINSPLINE
    eval_UBspline_3d_z_vgl (Spline, r(0), r(1), r(2),
			    &val, grad.v, &lap);
#else
    error("need EINSPLINE library for eval_UBspline_3d_d_vgl");
#endif
  }

  virtual void evaluate_spline (const Array1 < doublevar > & r, dcomplex & val, Array1 <dcomplex> & grad, Array2 <dcomplex> & hess){
    val=0;
    grad=dcomplex(0.0,0.0);
    hess=dcomplex(0.0,0.0);;
    for(int d=0;d<3;d++){
      if(!periodic){
	if((r(d)<origin(d)) || (r(d)>origin(d)+box_size(d)) ){
	  cout <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout <<" is outside the box, orbitals val, grad and hess are set to zero"<<endl;
	  return;
	}
      }
      else{
	if((r(d)<0) || (r(d)>1.0)){
	  cout  <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout  <<" is outside the box, orbitals val, grad and hess are set to zero"<<endl;
	  return;
	}
      }
    }
    
#ifdef USE_EINSPLINE
    eval_UBspline_3d_z_vgh (Spline, r(0), r(1), r(2), &val, grad.v, hess.v);
#else
 error("need EINSPLINE library for eval_UBspline_3d_d_vgh");
#endif
  }

 
  virtual void read (ifstream & ORB_0){
    //single_write(cout, "reading in plotfile\n");
    Array3 <doublevar> grid;
    get_grid_from_plot(ORB_0, grid, box_size, origin);
#ifdef USE_EINSPLINE
    BCtype_z xBC, yBC, zBC;
    Ugrid x_grid, y_grid, z_grid;
    if(periodic){
      xBC.lCode = xBC.rCode = PERIODIC;
      yBC.lCode = yBC.rCode = PERIODIC;
      zBC.lCode = zBC.rCode = PERIODIC;

      //spline calculation is in fractional coordinates of lattice vectors
      x_grid.start = 0.0; x_grid.end = 1.0; x_grid.num =grid.GetDim(0) ;
      y_grid.start = 0.0; y_grid.end = 1.0; y_grid.num =grid.GetDim(1) ;
      z_grid.start = 0.0; z_grid.end = 1.0; z_grid.num =grid.GetDim(2) ;

      //overwrite box_size from get_grid_from_plot
      for(int d=0;d<3;d++)
	box_size(d)=LatVec(0,d)+LatVec(1,d)+LatVec(2,d);

      Array1 <doublevar> resolution_array(3);
      for(int i=0;i<3;i++){
	  resolution_array(i)=1.0/grid.GetDim(i);
      }
       
      //make complex grid by multiplying by exp(-2*i*pi*k*r)
      //be aware that kvector in the code includes factor 2!

      Array3 <dcomplex> cgrid(x_grid.num,y_grid.num,z_grid.num);
      Array1 <doublevar> xyz(3);  //fractional coordinates of Lattice vectors
      for(int xx=0;xx<grid.GetDim(0);xx++){
	xyz(0)=xx*resolution_array(0);
	for(int yy=0;yy<grid.GetDim(1);yy++){
	  xyz(1)=yy*resolution_array(1);
	  for(int zz=0;zz<grid.GetDim(2);zz++){
	    xyz(2)=zz*resolution_array(2);
	    //cout <<xyz(0)<<"  "<<xyz(1)<<"  "<<xyz(2)<<endl;
	    dcomplex eikr=exp(-I*dot(xyz, kpointt));
	    cgrid(xx,yy,zz)=grid(xx,yy,zz)*eikr;
	    //if(xx==10 && yy==15)
	    // cout << xyz(2) <<"  "<<grid(xx,yy,zz)<<"  "<<cgrid(xx,yy,zz).real()<<"  "<<cgrid(xx,yy,zz).imag()<<"  "<<eikr.real()<<"  "<<eikr.imag()<<endl;
	  } 
	}
      }
      grid.clear();
      single_write(cout, "creating spline\n");
      Spline=create_UBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC, cgrid.v);
      single_write(cout, "done creating spline\n");
      cgrid.clear();

      //getting the spacing from the Spline object
      //adjusted by lenght of each latice vector;       
      Array1 <doublevar> LatVecNorm(3);
      LatVecNorm=0;
      for(int i=0;i<3;i++){
	for(int j=0;j<3;j++)
	  LatVecNorm(i)+=LatVec(i,j)*LatVec(i,j);	
	LatVecNorm(i)=sqrt(LatVecNorm(i));
      }
      
      spacing.Resize(3);
      spacing(0)=LatVecNorm(0)*Spline->x_grid.delta;
      spacing(1)=LatVecNorm(1)*Spline->y_grid.delta;
      spacing(2)=LatVecNorm(2)*Spline->z_grid.delta;
      
      
    }
    else{ //not periodic 
      xBC.lCode = xBC.rCode = NATURAL;
      yBC.lCode = yBC.rCode = NATURAL;
      zBC.lCode = zBC.rCode = NATURAL;
      x_grid.start = origin(0); x_grid.end = origin(0)+box_size(0); x_grid.num =grid.GetDim(0) ;
      y_grid.start = origin(1); y_grid.end = origin(1)+box_size(1); y_grid.num =grid.GetDim(1) ;
      z_grid.start = origin(2); z_grid.end = origin(2)+box_size(2); z_grid.num =grid.GetDim(2) ;
      
      //in the future will read in complex data on the grid, its all real now;
      Array3 <dcomplex> cgrid(x_grid.num,y_grid.num,z_grid.num);
      for(int xx=0;xx<grid.GetDim(0);xx++){
	for(int yy=0;yy<grid.GetDim(1);yy++){
	  for(int zz=0;zz<grid.GetDim(0);zz++){
	    cgrid(xx,yy,zz)=dcomplex(grid(xx,yy,zz),0.0);
	  }
	}
      }
      grid.clear();
      single_write(cout, "creating spline\n");
      Spline=create_UBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC, cgrid.v);
      single_write(cout, "done creating spline\n");
      cgrid.clear();
      //getting the spacing from the Spline object
      
      spacing(0)=Spline->x_grid.delta;
      spacing(1)=Spline->y_grid.delta;
      spacing(2)=Spline->z_grid.delta;
    }
#else
 error("need EINSPLINE library for create_UBspline_3d_d");
#endif
  }
  virtual ~EinsplineOrbComplex() { };
};


//----------------------------------------------------------------------------
class MultiEinsplineOrb {
 protected:
  int periodic;
  Array1 <doublevar> origin, box_size, spacing;
  Array2 <doublevar> RecipLatVec, LatVec;
  Array1 <doublevar> kpointt;
 public:
  virtual void getorigin(Array1 <doublevar> & a){
    a.Resize(3);
    a=origin;
  }
  virtual void getbox_size(Array1 <doublevar> & a){
    a.Resize(3);
    a=box_size;
  }
  virtual void getspacing(Array1 <doublevar> & a){
    a.Resize(3);
    a=spacing;
  }
  virtual void setperiodic(int & i){
    periodic=i;
  }
  virtual void kpoint (Array1 < doublevar > & kvector){
    kpointt.Resize(3);
    kpointt=kvector;
  }
  virtual void setorigin (Array1 < doublevar > & o){
    origin.Resize(3);
    origin=o;
  }
  virtual void setRecipLattice(Array2 <doublevar> & a){
    RecipLatVec.Resize(3,3);
    RecipLatVec=a;
  }
  virtual void setBounds(Array2 <doublevar> & a){
    LatVec.Resize(3,3);
    LatVec=a;
  }
  virtual void evaluate_spline (const Array1 < doublevar > & r, Array1 <dcomplex> & val)
  { error ("evaluate_spline not defined!");};
  virtual void evaluate_spline  (const Array1 < doublevar > & r, Array1 <dcomplex> & val, 
				 Array1 < Array1 <dcomplex> > & grad, Array1 <dcomplex> & lap)
  {error ("evaluate_spline not defined!");};
  virtual void evaluate_spline  (const Array1 < doublevar > & r, Array1 <dcomplex> & val, 
				 Array1 < Array1 <dcomplex> > & grad, Array1 < Array2 <dcomplex> > & hess)
  {error ("evaluate_spline not defined!");};
  virtual void evaluate_spline (const Array1 < doublevar > & r, Array1 <doublevar> & val)
  {error ("evaluate_spline not defined!");};
  virtual void evaluate_spline  (const Array1 < doublevar > & r, Array1 <doublevar> & val, 
				 Array1 < Array1 <doublevar> > & grad, Array1 <doublevar> & lap)
  {error ("evaluate_spline not defined!");};
  virtual void evaluate_spline  (const Array1 < doublevar > & r, Array1 <doublevar> & val, 
				 Array1 < Array1 <doublevar> > & grad, Array1 < Array2 <doublevar> > & hess)
  {error ("evaluate_spline not defined!");};
  virtual void read (vector <string> & valfiles, int & nmo)=0;
  virtual ~MultiEinsplineOrb() { };
};
//-------------------------------------------------------------------
class MultiEinsplineOrbReal : public MultiEinsplineOrb {
 private:
  #ifdef USE_EINSPLINE
  multi_UBspline_3d_d *MultiSpline;
  #endif
 public:
  virtual void evaluate_spline (const Array1 < doublevar > & r, Array1 <doublevar> & val) {
    val=0;
    for(int d=0;d<3;d++){
      if(!periodic){
	if((r(d)<origin(d)) || (r(d)>origin(d)+box_size(d))){
	  cout  <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout <<" is outside the box, orbitals value is set zero"<<endl;
	  return;
	}
      }
      else{
	if((r(d)<0) || (r(d)>1.0)){
	  cout  <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout <<" is outside the box, orbitals value is set zero"<<endl;
	  return;
	}
      }
    }

#ifdef USE_EINSPLINE
    eval_multi_UBspline_3d_d (MultiSpline, r(0), r(1), r(2), val.v);
#else
    error("need EINSPLINE library for eval_UBspline_3d_d");
#endif
  }

  virtual void evaluate_spline  (const Array1 < doublevar > & r, Array1 <doublevar> & val, 
				    Array1 < Array1 <doublevar> > & grad, 
				    Array1 <doublevar> & lap){ 
    val=lap=0;
    for(int m=0;m<grad.GetSize();m++){
      grad(m)=0;
    }

    for(int d=0;d<3;d++){
      if(!periodic){
	if((r(d)<origin(d)) || (r(d)>origin(d)+box_size(d)) ){
	  cout <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout <<" is outside the box, orbitals val, grad and lap are set to zero"<<endl;
	  return;
	}
      }
      else{
	if((r(d)<0) || (r(d)>1.0)){
	  cout  <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout  <<" is outside the box, orbitals val, grad and lap are set to zero"<<endl;
	  return;
	}
      }
    }
   

    Array1 <doublevar> linear_grad(3*grad.GetSize());
#ifdef USE_EINSPLINE
    eval_multi_UBspline_3d_d_vgl (MultiSpline, r(0), r(1), r(2),
			    val.v, linear_grad.v, lap.v);

    for(int m=0;m<grad.GetSize();m++){
      for(int i=0;i<3;i++)
	grad(m)(i)=linear_grad(3*m+i);
    }

#else
    error("need EINSPLINE library for eval_UBspline_3d_d_vgl");
#endif
  }

  virtual void evaluate_spline  (const Array1 < doublevar > & r, Array1 <doublevar> & val, 
				 Array1 < Array1 <doublevar> > & grad, 
				 Array1 < Array2 <doublevar> > & hess){
    val=0;
    for(int m=0;m<grad.GetSize();m++){
      grad(m)=0;
      hess(m)=0;
    }
    for(int d=0;d<3;d++){
      if(!periodic){
	if((r(d)<origin(d)) || (r(d)>origin(d)+box_size(d)) ){
	  cout <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout <<" is outside the box, orbitals val, grad and hess are set to zero"<<endl;
	  return;
	}
      }
      else{
	if((r(d)<0) || (r(d)>1.0)){
	  cout  <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout  <<" is outside the box, orbitals val, grad and hess are set to zero"<<endl;
	  return;
	}
      }
    }
     Array1 <doublevar> linear_grad(3*grad.GetSize());
     Array1 <doublevar> linear_hess(9*grad.GetSize());

#ifdef USE_EINSPLINE
    eval_multi_UBspline_3d_d_vgh (MultiSpline, r(0), r(1), r(2), val.v, linear_grad.v, linear_hess.v);

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

  virtual void read (vector <string> & valfiles, int & nmo){
#ifdef USE_EINSPLINE
    for(int m=0; m < nmo; m++) {
      string orbfile_0=valfiles[m];
      ifstream ORB_0(orbfile_0.c_str());
      if(!ORB_0){
	error("couldn't find orb file ", orbfile_0);
      } 
      single_write(cout,"Reading Orbfile: ",orbfile_0,"\n"); 
      Array3 <doublevar> grid;
      get_grid_from_plot(ORB_0, grid, box_size, origin);
      BCtype_d xBC, yBC, zBC;
      Ugrid x_grid, y_grid, z_grid;
      if(m==0){
	if(periodic){
	  xBC.lCode = xBC.rCode = PERIODIC;
	  yBC.lCode = yBC.rCode = PERIODIC;
	  zBC.lCode = zBC.rCode = PERIODIC;
	  
	  //spline calculation is in fractional coordinates of lattice vectors
	  x_grid.start = 0.0; x_grid.end = 1.0; x_grid.num =grid.GetDim(0) ;
	  y_grid.start = 0.0; y_grid.end = 1.0; y_grid.num =grid.GetDim(1) ;
	  z_grid.start = 0.0; z_grid.end = 1.0; z_grid.num =grid.GetDim(2) ;
	  
	  //overwrite box_size from get_grid_from_plot
	  for(int d=0;d<3;d++)
	    box_size(d)=LatVec(0,d)+LatVec(1,d)+LatVec(2,d);
	}
	else{
	  xBC.lCode = xBC.rCode = NATURAL;
	  yBC.lCode = yBC.rCode = NATURAL;
	  zBC.lCode = zBC.rCode = NATURAL;
	  
	  x_grid.start = origin(0); x_grid.end = origin(0)+box_size(0); x_grid.num =grid.GetDim(0) ;
	  y_grid.start = origin(1); y_grid.end = origin(1)+box_size(1); y_grid.num =grid.GetDim(1) ;
	  z_grid.start = origin(2); z_grid.end = origin(2)+box_size(2); z_grid.num =grid.GetDim(2) ;
	}
	MultiSpline = create_multi_UBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC, nmo);
      }
      single_write(cout, "creating spline\n");
      set_multi_UBspline_3d_d (MultiSpline, m, grid.v);
      single_write(cout, "done creating spline\n");
      grid.clear() ;
    }//m

    //getting the spacing from the Spline object
    spacing.Resize(3); 
    if(!periodic){
      spacing(0)=MultiSpline->x_grid.delta;
      spacing(1)=MultiSpline->y_grid.delta;
      spacing(2)=MultiSpline->z_grid.delta;
      
    }
    else{
      //adjusted by lenght of each latice vector;    
      Array1 <doublevar> LatVecNorm(3);
      LatVecNorm=0;
      for(int i=0;i<3;i++){
	for(int j=0;j<3;j++)
	  LatVecNorm(i)+=LatVec(i,j)*LatVec(i,j);	
	  LatVecNorm(i)=sqrt(LatVecNorm(i));
      }
      spacing(0)=LatVecNorm(0)*MultiSpline->x_grid.delta;
      spacing(1)=LatVecNorm(1)*MultiSpline->y_grid.delta;
      spacing(2)=LatVecNorm(2)*MultiSpline->z_grid.delta;
    }
#else
	error("need EINSPLINE library for create_UBspline_3d_d");
#endif
  }
  virtual ~MultiEinsplineOrbReal() { };
};

//---------------------------------------------------------------------------------
class MultiEinsplineOrbComplex : public MultiEinsplineOrb {
   private:
  #ifdef USE_EINSPLINE
  multi_UBspline_3d_z *MultiSpline; 
  #endif
 public:
  virtual void evaluate_spline (const Array1 < doublevar > & r, Array1 < dcomplex > & val) {
    val=dcomplex(0.0,0.0);
    for(int d=0;d<3;d++){
      if(!periodic){
	if((r(d)<origin(d)) || (r(d)>origin(d)+box_size(d)) ){
	  cout <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout <<" is outside the box, orbital value is set to zero"<<endl;
	  return;
	}
      }
      else{
	if((r(d)<0) || (r(d)>1.0)){
	  cout  <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout  <<" is outside the box, orbital value is set to zero"<<endl;
	  return;
	}
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
      if(!periodic){
	if((r(d)<origin(d)) || (r(d)>origin(d)+box_size(d)) ){
	  cout <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout <<" is outside the box, orbitals val, grad and lap are set to zero"<<endl;
	  return;
	}
      }
      else{
	if((r(d)<0) || (r(d)>1.0)){
	  cout  <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout  <<" is outside the box, orbitals val, grad and lap are set to zero"<<endl;
	  return;
	}
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
      if(!periodic){
	if((r(d)<origin(d)) || (r(d)>origin(d)+box_size(d)) ){
	  cout <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout <<" is outside the box, orbitals val, grad and hess are set to zero"<<endl;
	  return;
	}
      }
      else{
	if((r(d)<0) || (r(d)>1.0)){
	  cout  <<"WARNING! electron at "<<r(0)<<"  "<<r(1)<<"  "<<r(2);
	  cout  <<" is outside the box, orbitals val, grad and hess are set to zero"<<endl;
	  return;
	}
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

 
  virtual void read (vector <string> & valfiles, int & nmo){
    Array1 <doublevar> resolution_array(3);
#ifdef USE_EINSPLINE
    for(int m=0; m < nmo; m++) {
      string orbfile_0=valfiles[m];
      ifstream ORB_0(orbfile_0.c_str());
      if(!ORB_0){
	error("couldn't find orb file ", orbfile_0);
      } 
      single_write(cout,"Reading Orbfile: ",orbfile_0,"\n"); 
      Array3 <doublevar> grid;
      get_grid_from_plot(ORB_0, grid, box_size, origin);

      if(periodic){
	if(m==0){
	   BCtype_z xBC, yBC, zBC;
	   Ugrid x_grid, y_grid, z_grid;
	  xBC.lCode = xBC.rCode = PERIODIC;
	  yBC.lCode = yBC.rCode = PERIODIC;
	  zBC.lCode = zBC.rCode = PERIODIC;
	  
	  //spline calculation is in fractional coordinates of lattice vectors
	  x_grid.start = 0.0; x_grid.end = 1.0; x_grid.num =grid.GetDim(0) ;
	  y_grid.start = 0.0; y_grid.end = 1.0; y_grid.num =grid.GetDim(1) ;
	  z_grid.start = 0.0; z_grid.end = 1.0; z_grid.num =grid.GetDim(2) ;
	  
	  //overwrite box_size from get_grid_from_plot
	  for(int d=0;d<3;d++)
	    box_size(d)=LatVec(0,d)+LatVec(1,d)+LatVec(2,d);
	  
	  for(int i=0;i<3;i++){
	    resolution_array(i)=1.0/grid.GetDim(i);
	  }
	  //make complex grid by multiplying by exp(-2*i*pi*k*r)
	  //be aware that kvector in the code includes factor 2!
	  MultiSpline = create_multi_UBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC, nmo);
	}//m==0
	else{
	  if(MultiSpline->x_grid.num!=grid.GetDim(0) || 
	     MultiSpline->y_grid.num!=grid.GetDim(1) ||
	     MultiSpline->z_grid.num!=grid.GetDim(2))
	    error("Number of grid points is different that in the first plotfile");
	}
	
	Array3 <dcomplex> cgrid(grid.GetDim(0),grid.GetDim(1),grid.GetDim(2));
	Array1 <doublevar> xyz(3);  //fractional coordinates of Lattice vectors
	for(int xx=0;xx<grid.GetDim(0);xx++){
	  xyz(0)=xx*resolution_array(0);
	  for(int yy=0;yy<grid.GetDim(1);yy++){
	    xyz(1)=yy*resolution_array(1);
	    for(int zz=0;zz<grid.GetDim(2);zz++){
	      xyz(2)=zz*resolution_array(2);
	      //cout <<xyz(0)<<"  "<<xyz(1)<<"  "<<xyz(2)<<endl;
	      dcomplex eikr=exp(-I*dot(xyz, kpointt));
	      cgrid(xx,yy,zz)=grid(xx,yy,zz)*eikr;
	      //   if(xx==10 && yy==15)
	      //	cout << xyz(2) <<"  "<<grid(xx,yy,zz)<<"  "<<cgrid(xx,yy,zz).real()<<"  "<<cgrid(xx,yy,zz).imag()<<"  "<<eikr.real()<<"  "<<eikr.imag()<<endl;
	    } 
	  }
	}
	grid.clear();
	single_write(cout, "creating spline\n");
	set_multi_UBspline_3d_z (MultiSpline, m, cgrid.v);
	single_write(cout, "done creating spline\n");
	cgrid.clear();

	if(m==0){
	  //getting the spacing from the Spline object
	  //adjusted by lenght of each latice vector;       
	  Array1 <doublevar> LatVecNorm(3);
	  LatVecNorm=0;
	  for(int i=0;i<3;i++){
	    for(int j=0;j<3;j++)
	      LatVecNorm(i)+=LatVec(i,j)*LatVec(i,j);	
	    LatVecNorm(i)=sqrt(LatVecNorm(i));
	  }
	  
	  spacing.Resize(3);
	  spacing(0)=LatVecNorm(0)*MultiSpline->x_grid.delta;
	  spacing(1)=LatVecNorm(1)*MultiSpline->y_grid.delta;
	  spacing(2)=LatVecNorm(2)*MultiSpline->z_grid.delta;
	}
      }
      else{ //not periodic 
	if(m==0){
	  BCtype_z xBC, yBC, zBC;
	  Ugrid x_grid, y_grid, z_grid;
	  xBC.lCode = xBC.rCode = NATURAL;
	  yBC.lCode = yBC.rCode = NATURAL;
	  zBC.lCode = zBC.rCode = NATURAL;
	  x_grid.start = origin(0); x_grid.end = origin(0)+box_size(0); x_grid.num =grid.GetDim(0) ;
	  y_grid.start = origin(1); y_grid.end = origin(1)+box_size(1); y_grid.num =grid.GetDim(1) ;
	  z_grid.start = origin(2); z_grid.end = origin(2)+box_size(2); z_grid.num =grid.GetDim(2) ;
	  MultiSpline=create_multi_UBspline_3d_z (x_grid, y_grid, z_grid, xBC, yBC, zBC, nmo);
	}
	else{
	  if(MultiSpline->x_grid.num!=grid.GetDim(0) || 
	     MultiSpline->y_grid.num!=grid.GetDim(1) ||
	     MultiSpline->z_grid.num!=grid.GetDim(2))
	    error("Number of grid points is different that in the first plotfile");
	}
      
	//in the future will read in complex data on the grid, its all real now;
	Array3 <dcomplex> cgrid(grid.GetDim(0),grid.GetDim(1),grid.GetDim(2));
	for(int xx=0;xx<grid.GetDim(0);xx++){
	  for(int yy=0;yy<grid.GetDim(1);yy++){
	    for(int zz=0;zz<grid.GetDim(0);zz++){
	      cgrid(xx,yy,zz)=dcomplex(grid(xx,yy,zz),0.0);
	    }
	  }
	}
	grid.clear();
	single_write(cout, "creating spline\n");
	set_multi_UBspline_3d_z (MultiSpline, m, cgrid.v);
	single_write(cout, "done creating spline\n");
	cgrid.clear();
	//getting the spacing from the Spline object
	if(m==0){
	  spacing.Resize(3);
	  spacing(0)=MultiSpline->x_grid.delta;
	  spacing(1)=MultiSpline->y_grid.delta;
	  spacing(2)=MultiSpline->z_grid.delta;
	}
      }//not periodic
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
  Array1 <EinsplineOrb *> Einspline;
  Array1 <EinsplineOrbReal *> EinsplineReal;
  Array1 <EinsplineOrbComplex *> EinsplineComplex;
  MultiEinsplineOrb* MultiEinspline;
  MultiEinsplineOrbReal* MultiEinsplineReal;
  MultiEinsplineOrbComplex* MultiEinsplineComplex;
 
  vector <string> valfiles;
  Array1  <doublevar> kpoint;
  doublevar kpoint_square;
  Array1  <doublevar> origin;
  int complexspline;
  int periodic;
  int multi_spline;
  Array2  <doublevar> RecipLatVec;
  Array2  <doublevar> LatVec;
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
