/*
 
Copyright (C) 2012 Lucas K. Wagner

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

#include "Generate_sample.h"
#include "qmc_io.h"

void generate_sample(Sample_point * sample,
    Wavefunction * wf,
    Wavefunction_data * wfdata,
    Guiding_function * guidewf,
    int nconfig,
    Array1 <Config_save_point> & configs) { 

  Split_sampler sampler;
  sampler.setRecursionDepth(2);
  int nelectrons=sample->electronSize();
  
  //Warmup
  sample->randomGuess();
  bool have_moved_enough=false;
  Array1 <int> nmoves(nelectrons,0);
  int enough_moves=5;
  doublevar timestep=1.0;
  doublevar target_acceptance=0.7;
  doublevar tstep_delta=0.01;
  Dynamics_info dinfo;
  int nsteps=0;
  while(!have_moved_enough) { 
    for(int e=0; e< nelectrons; e++) { 
      int acc=sampler.sample(e,sample,wf,wfdata,guidewf,dinfo,timestep);
      if(dinfo.acceptance > target_acceptance) timestep+=tstep_delta;
      else timestep-=tstep_delta;
      if(acc>0) nmoves(e)++;
    }
    have_moved_enough=true;
    for(int e=0; e< nelectrons; e++) { 
      if(nmoves(e) < enough_moves) have_moved_enough=false; 
    }
    nsteps++;

  }

  debug_write(cout,"Took ",nsteps," steps to warm up.");

  //Now we generate the configurations
  enough_moves=3;
  configs.Resize(nconfig);
  for(int config=0; config < nconfig; config++) { 
    have_moved_enough=false;
    nmoves=0.0;
    while(!have_moved_enough) { 
      for(int e=0; e< nelectrons; e++) { 
        int acc=sampler.sample(e,sample,wf,wfdata,guidewf,dinfo,timestep);
// Adjusting the timestep here biases the walk slighly.
//        if(dinfo.acceptance > target_acceptance) timestep+=tstep_delta;
//        else timestep-=tstep_delta;
        if(acc>0) nmoves(e)++;
      }
      have_moved_enough=true;
      for(int e=0; e< nelectrons; e++) { 
        if(nmoves(e) < enough_moves) have_moved_enough=false; 
      }
    }
    configs(config).savePos(sample);
  }

}

