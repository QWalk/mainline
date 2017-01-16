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
//include/Cubic_spline3D.h


#ifndef CUBIC_SPLINE3D_H_INCLUDED
#define CUBIC_SPLINE3D_H_INCLUDED

#include "Array.h"
#include "Qmc_std.h"
using namespace std;

doublevar spline_1D_my(doublevar x,doublevar delta, Array1 <doublevar> & y);

void splin3_my(Array1 <doublevar> & rmin, doublevar delta, Array3 <doublevar> & yya,
	       Array1 <doublevar> & xx, doublevar & y);

void splin3_my(Array1 <doublevar> & rmin, Array1 <doublevar> & delta, Array3 <doublevar> & yya,
	       Array1 <doublevar> & xx, doublevar & y);

#endif //CUBIC_SPLINE3D_H_INCLUDED
//----------------------------------------------------------------------
