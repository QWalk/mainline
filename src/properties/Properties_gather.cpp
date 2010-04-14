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
#include "Properties_gather.h"
#include "ulec.h"
#include "qmc_io.h"
#include "Wavefunction_data.h"
//----------------------------------------------------------------------

Properties_gather::~Properties_gather() {
}

//----------------------------------------------------------------------

void Properties_gather::read(vector < string> & words) {

}


//----------------------------------------------------------------------
#include "Split_sample.h"

#include "Force_fitter.h"




void Properties_gather::gatherData(Properties_point & myprop,
                                   Pseudopotential * psp, 
                                   System * sys, 
                                   Wavefunction_data * wfdata,
                                   Wavefunction * wf, 
                                   Sample_point * sample, 
                                   Guiding_function * guide, int n_converge,
                                   int aux_updated) {  
  int nwf=wf->nfunc();
  
  myprop.setSize(nwf, 0,0);
  

  wf->updateLap(wfdata, sample);
  wf->getVal(wfdata, 0, myprop.wf_val);

  sys->calcKinetic(wfdata, sample, wf, myprop.kinetic);
  myprop.potential=sys->calcLoc(sample);

  for(int w=0; w< nwf; w++) {
    myprop.weight(w)=guide->getWeight(myprop.wf_val, 
                                      myprop.wf_val, w);
    myprop.sign(w)=myprop.wf_val.sign(w);      
  }

  int nrandvar=psp->nTest();
  Array1 <doublevar> rand_num(nrandvar);
  for(int i=0; i< nrandvar; i++) 
    rand_num(i)=rng.ulec();

  psp->calcNonlocWithTest(wfdata, sys,sample, wf,
                          rand_num,  myprop.nonlocal);

  myprop.count=1;


}
