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
#ifndef PROPERTIES_BLOCK_INCLUDED
#define PROPERTIES_BLOCK_INCLUDED
#include "Qmc_std.h"

namespace Properties_types {
enum quantities {total_energy,
                 kinetic,
                 potential,
                 nonlocal, 
                 weight,
                 NUM_QUANTITIES };
extern const char*  names[];
//const char * names[]={"total_energy", 
//            "kinetic",
//            "potential",
//            "nonlocal",
//           "weight"};

};
                 



/*!
  The average over a block.
 */
struct Properties_block {
  


  void setSize(int nwf, int naux, int n_aux_cvg=1) {
    using namespace Properties_types;

    avg.Resize(NUM_QUANTITIES, nwf);
    var.Resize(NUM_QUANTITIES, nwf);

    totweight=0;

    aux_energy.Resize(naux, n_aux_cvg);
    aux_weight.Resize(naux, n_aux_cvg);
    aux_energyvar.Resize(naux, n_aux_cvg);
    aux_weightvar.Resize(naux, n_aux_cvg);
    aux_diff.Resize(naux, n_aux_cvg);
    aux_diffvar.Resize(naux, n_aux_cvg);
    
    aux_size.Resize(naux);
    aux_size=1;
    z_pol.Resize(3);
    aux_z_pol.Resize(naux, 3);
    aux_weight_correlation=0;
  }

  void storeToLog(string & indent, ostream & os, string & label);
  void restoreFromLog(vector <string> & words);
  void printBlockSummary(ostream & os);
  void reduceBlocks(Array1 <Properties_block> & blocks, 
                    int start, int end);


  doublevar totweight;
  Array2 <doublevar> avg;
  Array2 <doublevar> var;


  Array2 <doublevar> autocorr; //autocorrelation of energy

  Array2 <doublevar> aux_energy;
  Array2 <doublevar> aux_weight;
  Array2 <doublevar> aux_energyvar;
  Array2 <doublevar> aux_weightvar;

  Array2 <doublevar> aux_diff;
  Array2 <doublevar> aux_diffvar;
  Array1 <doublevar> aux_size;

  Array2 <doublevar> aux_autocorr;
  
  Array1 <dcomplex> z_pol;

  Array2 <dcomplex> aux_z_pol;

  doublevar aux_weight_correlation;

};




#endif //PROPERTIES_BLOCK_INCLUDED
