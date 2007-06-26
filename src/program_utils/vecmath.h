/*
 A vector library for the Arrays
 Note: Some of these return arrays, so they're slow..don't use them 
 in speed-sensitive routines!

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
#ifndef VECMATH_H_INCLUDED
#define VECMATH_H_INCLUDED

#include "Array.h"
#include <cmath>
/*!
 */

inline void cross(const Array1 <double> & a, const Array1 <double> & b, Array1 <double> & c) {
  assert(a.GetDim(0)==3);
  assert(b.GetDim(0)==3);
  assert(c.GetDim(0)==3);
  c(0)=(a[1]*b[2]-a[2]*b[1]);
  c(1)=(a[2]*b[0]-a[0]*b[2]);
  c(2)=(a[0]*b[1]-a[1]*b[0]);
}

inline double dot(const Array1 <double> & a,const  Array1 <double> & b) {
  assert(a.GetDim(0)==3);
  assert(b.GetDim(0)==3);
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

inline double length_vec(const Array1 <double> & a) {
  assert(a.GetDim(0)==3);
  return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

inline Array1 <double> operator+(const Array1 <double> & a, 
                                 const Array1 <double> & b) {
  assert(a.GetDim(0)==3);
  assert(b.GetDim(0)==3);
  Array1 <double> c(3);
  for(int i=0; i< 3; i++) {
    c(i)=(a[i]+b[i]);
  }
  return c;
}

inline double distance_vec(const Array1 <double> & a,
                           const Array1 <double> & b) {
  assert(a.GetDim(0)==3);
  assert(b.GetDim(0)==3);
  return sqrt(
	      (a[0]-b[0])*(a[0]-b[0])
	      +(a[1]-b[1])*(a[1]-b[1])
	      +(a[2]-b[2])*(a[2]-b[2]));
}

inline Array1 <double> operator-(const Array1 <double> & a, 
                                 const Array1 <double> & b) {
  assert(a.GetDim(0)==3);
  assert(b.GetDim(0)==3);
  Array1 <double> c;
  for(int i=0; i< 3; i++) {
    c(i)=a[i]-b[i];
  }
  return c;
}

inline Array1 <double> operator*(const Array1 <double> & v,const  int & i) {
  assert(v.GetDim(0)==3);
  Array1 <double> c(3);
  for(int i=0; i< 3; i++) c(i)=(i*v[i]);
  return c;
}


//projection of a onto b (returns a dot b/abs(b))
inline double projection(const Array1 <double> & a,const  Array1 <double> & b) {
  assert(a.GetDim(0)==3);
  assert(b.GetDim(0)==3);
  double c;
  c=dot(a,b);
  double bsize=dot(b,b);
  return c/sqrt(bsize);
}






#endif //VECMATH_H_INCLUDED

//----------------------------------------------------------------------

