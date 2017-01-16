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

#ifndef PBC_ENFORCER_H_INCLUDED
#define PBC_ENFORCER_H_INCLUDED

#include "Array.h"


/*!
A class to enforce periodic boundary conditions, given lattice vectors
and an origin
*/
class Pbc_enforcer {
public:
  Pbc_enforcer():init_called(false), ndim(3) {}

  /*!
    Set up with lattice vectors.  A 3x3 matrix with the first lattice vector
    as [latVec(0,0), latVec(0,1), latVec(0,2)], the second as latVec(1,...) and so on.
  */
  void init(Array2 <doublevar> & latVec);
  
  /*!
    Set the origin to measure from.  Takes a 3-component array with the position of
    the origin.
  */
  void setOrigin(Array1 <doublevar> & origin);

  /*!
    pos is a 3-coordinate array.  Will change the position so that it is
    within the equivalent position in the simulation cell.  Returns 1 if
    the position was changed, 0 if not
  */
  int enforcePbc(Array1 <doublevar> & pos);


  /*!
    Returns whether or not pos is inside the simulation cell.  1 if
    it is, 0 if not.
  */
  int isInside(Array1 <doublevar>   & pos);
private:
  Array2 <doublevar> _latVec;
  Array2 <doublevar> _normVec;
  Array2 <doublevar> _corners;
  Array1 <doublevar> _origin;
  bool init_called;
  int ndim;
};

double find_centers(const Array1 <doublevar> & origin_,
                    const Array2 <doublevar> & latvec_,
                    const Array2 <doublevar> & atompos,
                    Array2 <doublevar> & centerpos,
                    Array1 <int> & equiv_atom, 
                    Array2 <int> & center_displacements,
                    doublevar cutoff_divider);

#endif //PBC_ENFORCER_H_INCLUDED

//----------------------------------------------------------------------
