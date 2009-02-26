/*
 
Copyright (C) 2007 Zachary Helms
 with further modifications by Lucas K. Wagner

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


#ifndef PLOT_METHOD_H_INCLUDED
#define PLOT_METHOD_H_INCLUDED

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

<b> RESOLUTION </b> Grid coarsness. for example 0.1(fine) or 10(coarse)

<b> MINMAX </b> list of 6 numbers: xmin xmax ymin ymax zmin zmax



<b> Optional </b>

<b> PRINT_DERIVATIVES </b> also print the gradient and laplacian of the molecular orbitals

<b> PLOTORBITALS </b>  A list of orbitals to print.  default: all of them


for example:<br>
<pre>
METHOD {
 PLOT

 PLOTORBITALS { 5 6 }
 RESOLUTION 0.5
 MINMAX { -10 10 -10 10 -10 10 }
 ORBITALS{
  INCLUDE n2.basis
  NMO 5
  ORBFILE n2.orb
  CENTERS { USEATOMS }
 }
}
</pre>
<br>

 */
class Plot_method : public Qmc_method
{
public:
  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  void run(Program_options & options, ostream & output);
  int showinfo(ostream & os);
  Plot_method():mymomat(NULL), mywalker(NULL) {}
  ~Plot_method()
  {
    if(mymomat!=NULL) delete mymomat;
    if(mywalker != NULL) delete mywalker;
    //deallocate(wf);
    //deallocate(wfdata);
    //for(int i=0; i< nconfig; i++) { delete electrons(i); }
  }

private:
  Array1 <int> orbs;	//orbital #'s to be plotted
  Array1 < Array1 <int> > orblist_pernode; //orbital list to be plotted on each node 
  Array1 <doublevar> minmax;	//xmin xmax ymin ymax zmin zmax
  doublevar resolution;	//grid coarsness: 10=coarser 0.1=finer
  MO_matrix * mymomat; //will hold MO information
  Complex_MO_matrix * cmymomat ;//will hold MO information
  Sample_point * mywalker; //a single configuration/walker
  Array2 <doublevar> mymovals; //(i,j) where i=MO#, j=0 default (for now)
  Array2 <dcomplex> cmymovals; //(i,j) where i=MO#, j=0 default (for now)
  System * sysprop;
  int print_derivatives;
  int jeep_like_cube_file;
  int periodic; 
  int use_complex;
  int totmo;
  Array2 <doublevar> LatticeVec;
  Array1 <doublevar> kpoint,origin;
};

#endif //PLOT_METHOD_H_INCLUDED
//------------------------------------------------------------------------
