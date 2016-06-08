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
#ifndef OVERLAPS_H_INCLUDED
#define OVERLAPS_H_INCLUDED

#include "Array.h"
#include "gaussian_set.h"

void pw_overlap(string & wfin, //JEEP wave function file
                Array2 <doublevar> & gvec, //g-vectors
                Array1 <Center> & centers, 
                Array1 <Contracted_gaussian> & basis,
                Gaussian_lookups & lookup,
                Array2 <doublevar> & pw_overlap, //return:
                Array1 <int> & mos //range of mo's to calculate on this node
                );

void calculate_overlap(Array2 <doublevar> & latvec, Array1 <doublevar> & origin,
                       Array1 <Center> & centers, Array1 <Contracted_gaussian> & basis,
                       Gaussian_lookups & lookup, Array2 <doublevar> & lcao_overlap);
                       
bool compare_mo(Array2 <doublevar> & oldMOCoeff,
                Array2 <doublevar> & newMOCoeff,
                Array2 <doublevar> & lcao_overlap,
                Array1 <int> compare_list);


void read_mo(string & orbin, int nmo, int nfunc, Array2 <doublevar> & moCoeff);
#endif //OVERLAPS_H_INCLUDED
