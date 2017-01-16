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
//----------------------------------------------------------------------
//src/Cubic_spline3D.cpp
#include "Array.h"
#include "Qmc_std.h"
using namespace std;


//use delta in interval 0.05 to 0.25 (100 to 25 points on 6x6x6 space)
static Array1 <doublevar> a(8);


inline doublevar spline_1D_my(doublevar x,doublevar delta, Array1 <doublevar> & y)
  //given 7 y values on the grid with equal spacing delta, routine calculates 2 splines
  //around y(3)=f(0) and returns value at point x
{


  //calculating spline
  a(0)=a(4)=y(3);
  a(1)=a(5)=((-y(0)-37.0*y(2)+37.0*y(4)+y(6))/48.0 +(y(1)-y(5))/6.0)/delta;
  a(2)=a(6)=((y(0)+37.0*y(2)+37.0*y(4)+y(6))/28.0 +(-2.0*y(1)-15.0*y(3)-2*y(5))/7.0)
    /(delta*delta);
  a(3)=((19.0*y(0)+367.0*y(2)+185.0*y(4)+5.0*y(6))/336.0 +
	(-19.0*y(1)-48.0*y(3)-5.0*y(5))/42.0)/(delta*delta*delta);

  a(7)=((-5.0*y(0)-185.0*y(2)-367.0*y(4)-19.0*y(6))/336.0 +
	 (5.0*y(1)+48.0*y(3)+19.0*y(5))/42.0)/(delta*delta*delta);

  //determining the returned value
   if (x>0.0) {
     return a(4)+ a(5)*x+ a(6)*x*x+ a(7)*x*x*x;
   }

  else {
    return a(0)+ a(1)*x+ a(2)*x*x+ a(3)*x*x*x;
  }

}


static Array2 <doublevar> yya_tmp(7,7);
static Array1 <doublevar> y1_tmp(7),y2_tmp(7), y3_tmp(7);
static Array1 <doublevar> x_center(3);
static Array1 <int> i_center(3);


void splin3_my(Array1 <doublevar> & rmin, doublevar delta, Array3 <doublevar> & yya,
		     Array1 <doublevar> & xx, doublevar & y)
  //given x_min, y_min, z_min (in vector rmin), mesh step delta and mesh values yya
  //routine gives function value y at vector xx
{

  //finding a center of 4x4x4 cube
  i_center(0)=int(0.5+(xx(0)-rmin(0))/delta);
  i_center(1)=int(0.5+(xx(1)-rmin(1))/delta);
  i_center(2)=int(0.5+(xx(2)-rmin(2))/delta);

  x_center(0)=rmin(0)+delta*i_center(0);
  x_center(1)=rmin(1)+delta*i_center(1);
  x_center(2)=rmin(2)+delta*i_center(2);

  //calculating splines in z direction (fastest running)
  for (int i=0; i<7 ;i++)
    for (int j=0; j<7 ;j++){
      for (int k=0; k<7 ;k++) {
	y3_tmp(k)=yya(i+i_center(0)-3,j+i_center(1)-3,k+i_center(2)-3);
      }
      yya_tmp(i,j)=spline_1D_my(xx(2)-x_center(2), delta,y3_tmp);
    }

  //calculating splines in y direction (second fastest running)
  for (int i=0;i<7;i++){
    for (int j=0;j<7;j++) y2_tmp(j)=yya_tmp(i,j);
    y1_tmp(i)=spline_1D_my(xx(1)-x_center(1), delta ,y2_tmp);
  }
  //final spline over x direction
  y=spline_1D_my(xx(0)-x_center(0), delta, y1_tmp);
}

void splin3_my(Array1 <doublevar> & rmin, Array1 <doublevar> & delta, Array3 <doublevar> & yya,
		     Array1 <doublevar> & xx, doublevar & y)
  //given x_min, y_min, z_min (in vector rmin), mesh step delta and mesh values yya
  //routine gives function value y at vector xx
{
  //static Array2 <doublevar> yya_tmp(7,7);
  //static Array1 <doublevar> y1_tmp(7),y2_tmp(7), y3_tmp(7);
  //static Array1 <doublevar> x_center(3);
  //static Array1 <int> i_center(3);

  //finding a center of 4x4x4 cube
  for (int i=0;i<3;i++){
    i_center(i)=int(0.5+(xx(i)-rmin(i))/delta(i));
    x_center(i)=rmin(i)+delta(i)*i_center(i);
  }

  //calculating splines in z direction (fastest running)
  for (int i=0; i<7 ;i++)
    for (int j=0; j<7 ;j++){
      for (int k=0; k<7 ;k++) {
	y3_tmp(k)=yya(i+i_center(0)-3,j+i_center(1)-3,k+i_center(2)-3);
      }
      yya_tmp(i,j)=spline_1D_my(xx(2)-x_center(2), delta(2),y3_tmp);
    }

  //calculating splines in y direction (second fastest running)
  for (int i=0;i<7;i++){
    for (int j=0;j<7;j++) y2_tmp(j)=yya_tmp(i,j);
    y1_tmp(i)=spline_1D_my(xx(1)-x_center(1), delta(1) ,y2_tmp);
  }
  //final spline over x direction
  y=spline_1D_my(xx(0)-x_center(0), delta(0), y1_tmp);
}

//----------------------------------------------------------------------

