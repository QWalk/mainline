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
//------------------------------------------------------------------------
//include/Localize_method.h

#ifndef LOCALIZE_METHOD_H_INCLUDED
#define LOCALIZE_METHOD_H_INCLUDED

#include <vector>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include "Qmc_method.h"
#include "qmc_io.h"
#include "Array.h"
#include "MO_matrix.h"
#include "Sample_point.h"
class System;
class Wavefunction;
class Program_options;

/*!
\brief 
Generate a grid of sample points of MO functions. Keyword PLOT
in METHOD section

 the PLOT method generates a .xyz file that can be imported by gOpenMOl
 and "formatted" .plt files that need to be "unformatted"
 using gOpenMol's "Pltfile (conversion)" utility


<h3> Options </h3> 
<b> Required </b>

<b> ORBITALS </b> A list of the orbitals to be plotted

<b> RESOLUTION </b> Grid coarsness, from 0.1(fine) to 10(coarse)

<b> MINMAX </b> list of 6 numbers: xmin xmax ymin ymax zmin zmax

In the qmc input file
write "PLOT" in the method section and include 3 subsections "ORBITALS",
"RESOLUTION", and "MINMAX". Should output FORMATTED .plt files that can
be both unformatted and viewed in gOpenMol.
 */
class Localize_method : public Qmc_method
{
public:
  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  void run(Program_options & options, ostream & output);
  int showinfo(ostream & os);
  Localize_method():mymomat(NULL), mywalker(NULL), sysprop(NULL) {}
  ~Localize_method()
  {
    if(mymomat != NULL) delete mymomat;
    if(mywalker != NULL) delete mywalker;
    if(sysprop != NULL) delete sysprop;
  }

private:
  Array1 <int> orbs;	//!<orbital #'s to be plotted
  Array1 <int> all_orbs; 
  Array1 < Array1 <doublevar> > center; //!<positions of localization centers
  int nmo;
  doublevar resolution;	//!<grid coarsness: 10=coarser 0.1=finer
  doublevar radius; //!<radius of integration
  doublevar delta_norm; //!<norm of orbital left outside of radius of localization. 
  MO_matrix * mymomat; //!<will hold MO information
  Sample_point * mywalker; //!<a single configuration/walker
  Array2 <doublevar> mymovals; //!<(i,j) where i=MO#, j=0 default (for now)
  System * sysprop;
  Array2 <doublevar> latVec;
  Array1 <doublevar> origin;
  int ncenters;
  int norbs;
};

#endif //LOCALIZE_METHOD_H_INCLUDED
//------------------------------------------------------------------------
