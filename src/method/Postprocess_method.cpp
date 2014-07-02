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
#include "Postprocess_method.h"
#include "qmc_io.h"
#include "System.h"
#include "Sample_point.h"
#include "ulec.h"

/*!
*/
void Postprocess_method::read(vector <string> words,
                       unsigned int & pos,
                       Program_options & options)
{

  if(!readvalue(words, pos=0, configfile, "READCONFIG"))
    error("Need READCONFIG in Postprocess");


  vector <vector < string> > dens_words;
  vector<string> tmp_dens;
  pos=0;
  while(readsection(words, pos, tmp_dens, "DENSITY")) {
    dens_words.push_back(tmp_dens);
  }


  vector <vector < string> > avg_words;
  pos=0;
  while(readsection(words, pos, tmp_dens, "AVERAGE")) {
    avg_words.push_back(tmp_dens);
  }
  
  sys=NULL;
  allocate(options.systemtext[0],  sys);
  sys->generatePseudo(options.pseudotext, pseudo);
  debug_write(cout, "wfdata allocate\n");
  wfdata=NULL;
  if(options.twftext.size() < 1) error("Need TRIALFUNC section for POSTPROCESS");
  allocate(options.twftext[0], sys, wfdata);



   
  average_var.Resize(avg_words.size());
  average_var=NULL;
  for(int i=0; i< average_var.GetDim(0); i++) { 
    allocate(avg_words[i], sys, wfdata, average_var(i));
  }

  densplt.Resize(dens_words.size());
  for(int i=0; i< densplt.GetDim(0); i++) {
    allocate(dens_words[i], sys, options.runid,densplt(i));
  }

  
}

void weighted_update_average(doublevar weight, doublevar xnew,doublevar &  xavg, doublevar &  xvar, doublevar & totweight)
  {
  doublevar nweight=totweight+weight;
  doublevar navg=(weight*xnew+totweight*xavg)/nweight;
  doublevar nvar=(weight*xnew*xnew+totweight*(xavg*xavg+xvar))/nweight-navg*navg;
  xavg=navg;
  xvar=nvar;
}

void Postprocess_method::run(Program_options & options, ostream & output) {
  Sample_point * sample=NULL;
  Wavefunction * wf=NULL;
  sys->generateSample(sample);
  wfdata->generateWavefunction(wf);
  sample->attachObserver(wf);
  Properties_gather gather;
  Primary guide;
  
  int nelec=sample->electronSize();
  int ndim=3;
  Config_save_point tmpconfig;
  FILE * f=fopen(configfile.c_str(),"r");
  Properties_point pt;
  pt.setSize(1);
  doublevar avgen=0,varen=0,avgwt=0,varwt=0;
  Array1 <Average_return> avgavg(average_var.GetDim(0)), varavg(average_var.GetDim(0));
  int npoints=0;
  int skip=2048; 
  fseek(f,0,SEEK_END);
  long int lSize=ftell(f);
  rewind(f);
  int npoints_tot=lSize/(sizeof(doublevar)*(nelec*3+1));
  output << "Estimated number of samples in this file: " << npoints_tot << endl;
  
  while(tmpconfig.readBinary(f,nelec,ndim)) {
    doublevar weight;
    if(!fread(&weight,sizeof(doublevar),1,f)) error("Misformatting in binary file",configfile);
    tmpconfig.restorePos(sample);

    ///--------------------------------Gather the data
    gather.gatherData(pt,pseudo,sys,wfdata,wf,sample,&guide);
    pt.avgrets.Resize(1,average_var.GetDim(0));
    for(int i=0; i< average_var.GetDim(0); i++) { 
      average_var(i)->randomize(wfdata,wf,sys,sample);
      average_var(i)->evaluate(wfdata, wf, sys, sample,pt, pt.avgrets(0,i));
    }
    for(int i=0; i< densplt.GetDim(0); i++) densplt(i)->accumulate(sample,weight);
    //-------------------Done gathering data

    //-----------------------------------Update our averages and variances
    weighted_update_average(pt.weight(0),pt.energy(0),avgen,varen,avgwt);
    for(int i=0; i< average_var.GetDim(0); i++) { 
      if(avgavg(i).vals.GetDim(0)==0) { 
        avgavg(i)=pt.avgrets(0,i);
        varavg(i)=pt.avgrets(0,i);
        avgavg(i).vals=0;
        varavg(i).vals=0;
      }
      for(int j=0; j< pt.avgrets(0,i).vals.GetDim(0); j++) { 
        weighted_update_average(pt.weight(0),pt.avgrets(0,i).vals(j),avgavg(i).vals(j),varavg(i).vals(j),avgwt);
      }
    }
    avgwt+=pt.weight(0);

    //--------------------------------Done with the update
    
    npoints++;
    doublevar progress=doublevar(npoints)/doublevar(npoints_tot);
    if(fabs(progress*10-int(progress*10)) < 1e-7) { 
      cout << "progress: " << progress*100 << "% done" << endl;
    }
    
  }
  fclose(f);
  for(int i=0; i< densplt.GetDim(0); i++) densplt(i)->write();

  output << "Averages" << endl;
  output << "total_energy " << avgen << " +/- " << sqrt(varen/npoints) << " (sigma " << sqrt(varen) << ") "<< endl;
  output << "weight " << avgwt/npoints << endl;
  for(int i=0; i< average_var.GetDim(0); i++) { 
    for(int j=0; j< varavg(i).vals.GetDim(0); j++) 
      varavg(i).vals(j)=sqrt(varavg(i).vals(j)/npoints);
    average_var(i)->write_summary(avgavg(i),varavg(i),output);
  }

  delete sample;
  delete wf;
}


/*!

*/
int Postprocess_method::showinfo(ostream & os)
{
  os<<"#############Postprocess_method#################\n";
  sys->showinfo(os);
  wfdata->showinfo(os);
  return 1;
}
