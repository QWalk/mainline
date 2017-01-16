/*
 
Copyright (C) 2008 Jindrich Kolorenc

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


#ifndef PLOT1D_METHOD_H_INCLUDED
#define PLOT1D_METHOD_H_INCLUDED

#include "Qmc_method.h"
#include "Wavefunction.h"
#include "Wavefunction_data.h"
class System;
class Wavefunction;
class Program_options;

/*!
\brief 
Plotting of one-dimensional parts of a wave function (mainly Jastrow).
Functions are plotted from 0 to CUTOFF, number of points is NGRID.
Default output file is runid.plt1d, suffix plt1d can be overriden by
setting SUFFIX.
 */
class Plot1D_method : public Qmc_method
{
 public:

  void read(vector <string> words,
	    unsigned int & pos,
	    Program_options & options);
  int generateVariables(Program_options & options);
  void run(Program_options & options, ostream & output);
  int showinfo(ostream & os);
  Plot1D_method() {}
  ~Plot1D_method() {}
  

 private:

  System * sysprop;
  Wavefunction * wf;  
  Wavefunction_data * mywfdata;

  int ngrid;
  doublevar cutoff;
  string suffix;

};

#endif //PLOT1D_METHOD_H_INCLUDED
//------------------------------------------------------------------------
