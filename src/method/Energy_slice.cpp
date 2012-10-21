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

#include "Program_options.h"
#include "Energy_slice.h"
#include "qmc_io.h"
#include "System.h"
#include "Sample_point.h"
#include "ulec.h"
#include "Generate_sample.h"
#include "Wavefunction_data.h"
#include "Properties_gather.h"
#include <algorithm>
#include "Properties.h"
/*!
Read the "words" from the method section in the input file
via doinput() parsing, and store section information in private
variables orbs, resolution, and minmax.
Set up MO_matrix and Sample_point objects for wavefunction from input
*/
void Energy_slice::read(vector <string> words,
                       unsigned int & pos,
                       Program_options & options)
{
  sys=NULL;
  allocate(options.systemtext[0],  sys);
  sample=NULL;
  sys->generateSample(sample);

  wf=NULL;
  allocate(options.twftext[0],sys,wfdata);
  

  pseudo=NULL;
  sys->generatePseudo(options.pseudotext, pseudo);
  
  if(!readvalue(words,pos=0,nconfig,"NCONFIG"))
    nconfig=100;
  if(!readvalue(words,pos=0, nbin, "NBINS"))
    nbin=5;

  pos=0;
  vector <string> tmp_dens;
  while(readsection(words, pos, tmp_dens, "AVERAGE")) {
    avg_words.push_back(tmp_dens);
  }

  average_var.Resize(avg_words.size());
  average_var=NULL;
  for(int i=0; i< average_var.GetDim(0); i++) { 
    allocate(avg_words[i], sys, wfdata, average_var(i));
  }


}


//----------------------------------------------------------------------

void Energy_slice::run(Program_options & options, ostream & output) {
  Properties_gather mygather;

  sys->generateSample(sample);
  wfdata->generateWavefunction(wf);
  sample->attachObserver(wf);
  Array1 <Config_save_point> configs(nconfig);
  Vmc_sum_squares guide;
  generate_sample(sample,wf,wfdata,&guide,nconfig, configs);
  vector <string> tmp;
  mygather.read(tmp);
  vector <pair<double,int> > en;
  Properties_point pt;
  
  cout << "first energy calculation " << endl;
  for(int w=0; w< nconfig; w++) { 
    configs(w).restorePos(sample);
    wf->updateLap(wfdata,sample);
    mygather.gatherData(pt, pseudo, sys, wfdata, wf, 
                            sample, &guide);
    //oen(w)=pt.energy(0);
    en.push_back(pair<double,int>(pt.energy(0),w));
  }

        
  cout << "sorting " << endl;
  //Find the n-tiles:
  sort(en.begin(), en.end());
  int step=nconfig/nbin;
  if(nconfig%nbin!=0) error("Need to have nbin divide nconfig");
  cout << "step " << step << endl;
  for(int b=0; b< nbin; b++) { 
    cout << "set " << b << " bottom energy " << en[b*step].first << " top energy " << en[(b+1)*step-1].first << endl;
  }


  vector<pair<double,int> >::iterator wen=en.begin();
  for(int b=0; b< nbin; b++) { 
    Properties_manager prop;

    string log_label="post";
    append_number(log_label,b);
    string logfile=options.runid+".log";
    prop.setLog(logfile,log_label);
    int nblock=20;
    int nstep=step/nblock;
    prop.setSize(wf->nfunc(), nblock, nstep, 1, sys, wfdata);
    prop.initializeLog(average_var);
    
    int currstep=0;
    for(int s=0; s< step; s++) { 
      int w=wen->second;
      Properties_point pt;
      configs(w).restorePos(sample);
      //cout << "updatelap" << endl;
      wf->updateLap(wfdata,sample);
      //cout << "gather " << endl;
      mygather.gatherData(pt, pseudo, sys, wfdata, wf, 
          sample, &guide);

      pt.avgrets.Resize(1,average_var.GetDim(0));
      for(int i=0; i< average_var.GetDim(0); i++) { 
        average_var(i)->randomize(wfdata,wf,sys,sample);
        average_var(i)->evaluate(wfdata, wf, sys, sample,pt, pt.avgrets(0,i));
      }
      pt.parent=0;
      pt.nchildren=1; pt.children(0)=1;
      pt.weight=1.0;
      prop.insertPoint(currstep, 0, pt);
      

      wen++;

      currstep++;
      if(currstep==nstep) { 
        prop.endBlock();
        currstep=0;
      }
      
    }
  }
    
    

}


/*!

*/
int Energy_slice::showinfo(ostream & os)
{
  os<<"#############Energy_slice#################\n";
  sys->showinfo(os);

  return 1;
}
