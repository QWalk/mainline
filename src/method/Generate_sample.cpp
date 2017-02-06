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

doublevar generate_sample(Sample_point * sample,
    Wavefunction * wf,
    Wavefunction_data * wfdata,
    Guiding_function * guidewf,
    int nconfig,
    Array1 <Config_save_point> & configs,
    int warmup_moves,
    int decorrelation_moves,
    doublevar timestep) { 

  Split_sampler sampler;
  sampler.setRecursionDepth(2);
  sampler.setDivider(2.0);
  int nelectrons=sample->electronSize();
  
  //Warmup
  sample->randomGuess();
  Array1 <int> nmoves(nelectrons,0);
  doublevar target_acceptance=0.7;
  doublevar tstep_delta=0.01;
  Dynamics_info dinfo;
  int nsteps=0;
  for(int step=0; step < warmup_moves; step++) { 
    for(int e=0; e< nelectrons; e++) { 
      int acc=sampler.sample(e,sample,wf,wfdata,guidewf,dinfo,timestep);
      if(dinfo.acceptance > target_acceptance) timestep+=tstep_delta;
      else timestep-=tstep_delta;
    }
  }

  debug_write(cout,"Timestep is ",timestep,"\n");

  //Now we generate the configurations. 
  //We don't change the time step from here out.
  configs.Resize(nconfig);
  for(int config=0; config < nconfig; config++) { 
    for(int step=0; step < decorrelation_moves; step++) { 
      for(int e=0; e< nelectrons; e++) { 
        int acc=sampler.sample(e,sample,wf,wfdata,guidewf,dinfo,timestep);
      }
    }
    configs(config).savePos(sample);
  }
  
  return timestep;

}

//######################################################################

#include "ulec.h"

doublevar objective(MO_matrix * mo,Sample_point * sample,int list, 
    Array2 <doublevar> & vals) { 
  
  vals=0.0;
  mo->updateVal(sample,0,list,vals);
  doublevar f=0.0;
  for(int m=0; m < vals.GetDim(0); m++) {
    f+=vals(m,0)*vals(m,0);
  }
  return f;
}


void generate_mo_sample(Sample_point * sample, System * sys,
    MO_matrix * mo, int list, int nconfig, Array1 <Array1 <doublevar> > & r) {
  int norb=mo->getNmo();
  Array2 <doublevar> vals(norb,1);
  r.Resize(nconfig);
  for(int i=0; i< nconfig; i++) r[i].Resize(3);

  //warmup
  int nacc=0;
  int n_must_acc=20;
  doublevar tstep=1.0,tstep_delta=0.01;
  sample->randomGuess();
  Array1 <doublevar> delta(3),oldpos(3),newpos(3);
  doublevar old_func=objective(mo,sample,list,vals);
  doublevar new_func=0;
  int e=0;
  while(nacc < n_must_acc) { 
    for(int d=0; d< 3; d++) delta[d]=rng.gasdev()*sqrt(tstep);
    sample->getElectronPos(e,oldpos);
    for(int d=0; d< 3; d++) newpos(d)=oldpos(d)+delta(d);
    sample->setElectronPos(e,newpos);
    new_func=objective(mo,sample,list,vals);
    if(rng.ulec() < new_func/old_func) { 
      old_func=new_func;
      tstep+=tstep_delta;
      nacc++;
    }
    else { 
      sample->setElectronPos(e,oldpos);
      tstep-=tstep_delta;
    }
  }

  for(int config=0; config < nconfig; config++) { 
    nacc=0;
    while(nacc < n_must_acc) { 
      for(int d=0; d< 3; d++) delta[d]=rng.gasdev()*sqrt(tstep);
      sample->getElectronPos(e,oldpos);
      for(int d=0; d< 3; d++) newpos(d)=oldpos(d)+delta(d);
      sample->setElectronPos(e,newpos);
      new_func=objective(mo,sample,list,vals);
      if(rng.ulec() < new_func/old_func) { 
        old_func=new_func;
        nacc++;
      }
      else { 
        sample->setElectronPos(e,oldpos);
      }
    }
    sample->getElectronPos(e,r(config));
  }
}

