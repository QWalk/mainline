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


#ifndef JEEP_UTILS_H_INCLUDED
#define JEEP_UTILS_H_INCLUDED
#include "Qmc_std.h"
using namespace std;

void get_grid_from_plot(istream & plotfile,
                        Array3 <doublevar> & grid,
                        Array1 <doublevar> & box_size,
                        Array1 <doublevar> & origin);

void get_grid_from_plot(istream & plotfile,
                        Array3 <doublevar> & grid,
                        Array1 <doublevar> & box_size,
                        Array1 <doublevar> & origin, 
                        Array1 <int > & padding, Array1 <int> & startfill);


void get_function_header(istream & plotfile,
                Array1 <int> & npoints,
                Array1 <doublevar> & box_size,
                Array1 <doublevar> & origin);

void get_global_header(istream & plotfile, int & nfunctions);

doublevar get_cutoff_radius(Array3 <doublevar> & grid, Array1 <doublevar> & box_size,
                       Array1 <doublevar> & origin, Array1 <doublevar> & center);

void get_wannier_centers(istream & is, Array2 <doublevar> & centers);

doublevar get_periodic_distance(const Array1 <doublevar> & r1,
				const Array1 <doublevar> & r2,
				const Array2 <doublevar> & latVec);


/*!
  Get the information contained in the header of the JEEP file.
 */
void get_jeep_header(FILE * wfin, int & nst, int & ngw, Array1 <doublevar> & occupation);

/*!
  Read the next molecular orbital into the coefficient array; 
  returns in the form (g-vector, [real, imaginary])
 */
void get_next_pw_mo(FILE * wfin, int ngw,  Array2 <doublevar> & coeff);


/*!
  Index the places where the molecular orbitals start in the 
  JEEP wavefunction file.
 */
void summarize_jeep_wf(string & filename, int & nst, int & ngw, 
                       Array1 <int> & mo_places //!< location in file where each mo is
                       );

#endif //JEEP_UTILS_H_INCLUDED


//----------------------------------------------------------------------
