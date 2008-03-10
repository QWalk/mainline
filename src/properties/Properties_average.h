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
#ifndef PROPERTIES_AVERAGE_INCLUDED
#define PROPERTIES_AVERAGE_INCLUDED
#include "Qmc_std.h"
#include "Properties_block.h"
#include <iomanip>
#include "Average_generator.h"

struct Properties_final_average {

  Properties_final_average():show_autocorr(0) {}
  void setSize(int nwf, int naux, int n_cvg=1);


  void blockReduce(Array1 <Properties_block> & block_avg,
                   int start_block, int end_block, 
                   int find_equil=0);

  void averageReduce(Array1 <Properties_final_average> & avg,
		     int start, int end);


  void twoPointForces(Array2 <doublevar> & forces);

  void showSummary(ostream & os, Array1 <Average_generator*> avg_gen);

  void showAutocorr(int i) {
    assert(i==0 || i==1);
    show_autocorr=i;
  }

  int threw_out;

  doublevar totweight;
  
  Array2 <doublevar> avg;
  Array2 <doublevar> err;

  Array2 <doublevar> avgvar; //!< average variances


  Array1 <Average_return> avgavg;
  Array1 <Average_return> avgerr;
  
  Array1 <doublevar> diff_energy;
  Array1 <doublevar> diff_energyerr;

  Array1 <dcomplex> z_pol;
  Array1 <dcomplex> z_pol_err;

  Array2 <doublevar> autocorr; //autocorrelation of energy
  Array2 <doublevar> autocorrerr; 

  Array2 <doublevar> aux_energy;
  Array2 <doublevar> aux_energyerr;
  Array2 <doublevar> aux_energyvar;

  Array2 <doublevar> aux_diff;
  Array2 <doublevar> aux_differr;

  Array2 <doublevar> aux_autocorr; 
  Array2 <doublevar> aux_autocorrerr; 
  
  Array2 <dcomplex> aux_z_poldiff;
  Array2 <dcomplex> aux_z_poldiff_err;

  Array1 <doublevar> aux_size;

  private:
  int show_autocorr;


};

#endif
