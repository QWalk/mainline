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
  
  if(!readvalue(words, pos=0, nskip, "NSKIP"))
    nskip=0;
  

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
//----------------------------------------------------------------------

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
  int npoints_tot=0;
  FILE * f;
  if(mpi_info.node==0) { 
    f=fopen(configfile.c_str(),"r");
    if(ferror(f)) error("Could not open",configfile);
    cout << "here " << endl;
    fseek(f,0,SEEK_END);
    long int lSize=ftell(f);
    rewind(f);
    npoints_tot=lSize/(sizeof(doublevar)*(nelec*3+1+4));
    output << "Estimated number of samples in this file: " << npoints_tot << endl;
    output << "We are skipping the first " << nskip << " of these " << endl;
    Config_save_point tmpconfig;
    for(int i=0; i< nskip; i++) { 
      doublevar weight;
      tmpconfig.readBinary(f,nelec,ndim,weight);
//     doublevar weight;
//      if(!fread(&weight,sizeof(doublevar),1,f)) error("Misformatting in binary file",configfile, " perhaps nskip is too large?");
    }
  }
    

#ifdef USE_MPI
  if(mpi_info.nprocs<2) error("POSTPROCESS must be run with at least 2 processes if it is run in parallel.");
  if(mpi_info.node==0) { 
    master(wf,sample,f,output);
  }
  else { 
    worker(wf,sample);
  }
#else
//  fseek(f,0,SEEK_END);
//  long int lSize=ftell(f);
//  rewind(f);
//  int npoints_tot=lSize/(sizeof(doublevar)*(nelec*3+1));
//  output << "Estimated number of samples in this file: " << npoints_tot << endl;
  Config_save_point tmpconfig;
  Properties_point pt;
  pt.setSize(1);
  int npoints=0;
  Postprocess_average postavg(average_var.GetDim(0));
  while(tmpconfig.readBinary(f,nelec,ndim)) {
    doublevar weight;
    if(!fread(&weight,sizeof(doublevar),1,f)) error("Misformatting in binary file",configfile);
    tmpconfig.restorePos(sample);
    gen_point(wf,sample,tmpconfig,weight,pt);
    postavg.update_average(pt);

    npoints++;
    doublevar progress=doublevar(npoints)/doublevar(npoints_tot);
    if(fabs(progress*10-int(progress*10)) < 0.5/npoints_tot) { 
      cout << "progress: " << progress*100 << "% done" << endl;
    }
    
  }
  postavg.print(average_var,output);
#endif //USE_MPI
  for(int i=0; i< densplt.GetDim(0); i++) densplt(i)->write();

  if(mpi_info.node==0) fclose(f);

  delete sample;
  delete wf;
}

//----------------------------------------------------------------------

/*!

*/
int Postprocess_method::showinfo(ostream & os) {
  os<<"#############Postprocess_method#################\n";
  sys->showinfo(os);
  wfdata->showinfo(os);
  return 1;
}
//----------------------------------------------------------------------
int Postprocess_method::gen_point(Wavefunction * wf, Sample_point * sample,
    Config_save_point & configpos, doublevar weight, Properties_point & pt) { 
  Primary guide;
  Properties_gather gather;
  
  configpos.restorePos(sample);
  
  pseudo->randomize();
  gather.gatherData(pt,pseudo,sys,wfdata,wf,sample,&guide);
  pt.avgrets.Resize(1,average_var.GetDim(0));
  //cout << mpi_info.node << " generating a point with " << average_var.GetDim(0) << " avgrets " << endl;
  for(int i=0; i< average_var.GetDim(0); i++) { 
    average_var(i)->randomize(wfdata,wf,sys,sample);
    average_var(i)->evaluate(wfdata, wf, sys, sample,pt, pt.avgrets(0,i));
  }
  for(int i=0; i< densplt.GetDim(0); i++) densplt(i)->accumulate(sample,weight);
  pt.weight(0)=weight;
  
}

//----------------------------------------------------------------------
int Postprocess_method::worker(Wavefunction * wf, Sample_point * sample) { 
#ifdef USE_MPI
  Config_save_point tmpconfig;
  doublevar weight;
  Properties_point pt;
  pt.setSize(1);
  tmpconfig.mpiReceive(0);
  MPI_Recv(weight,0);
  MPI_Status status;
  
  while(true) { 
    gen_point(wf,sample,tmpconfig,weight,pt);
    int done=1;
    MPI_Send(done,0);
    pt.mpiSend(0);
    MPI_Recv(done,0);
    if(done==0) break;
    tmpconfig.mpiReceive(0);
    MPI_Recv(weight,0);
  }
  cout << mpi_info.node << " : done " << endl;

#endif //USE_MPI
}
//----------------------------------------------------------------------

int Postprocess_method::master(Wavefunction * wf, Sample_point * sample,FILE * f, ostream & os) { 
#ifdef USE_MPI
  Config_save_point tmpconfig;
  doublevar weight;
  Properties_point pt;
  pt.setSize(1);
  int nelec=sample->electronSize();
  int ndim=3;
  MPI_Status status;
  Postprocess_average postavg(average_var.GetDim(0));
  

  //Get everyone started with data
  cout << "master: sending initial data" << endl;
  for(int r=1; r < mpi_info.nprocs; r++) { 
     if(!tmpconfig.readBinary(f,nelec,ndim,weight)) {  
       error("Binary file may not contain enough walkers; finished after ",r);
     }
     tmpconfig.mpiSend(r);
     MPI_Send(weight,r);
  }
  
  int totcount=0;
  cout << "master : going through file " << endl;
  while(tmpconfig.readBinary(f,nelec,ndim,weight)) {
//    doublevar weight;
//    if(!fread(&weight,sizeof(doublevar),1,f)) error("Misformatting in binary file",configfile);

    //Is anyone done?
    //When done, receive completed point and send out new point
    int done;
    MPI_Recv(&done,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_Comm_grp,&status);
    //cout << "master: received " << done << " from " << status.MPI_SOURCE << endl;
    done=1;
    pt.mpiReceive(status.MPI_SOURCE);
    MPI_Send(done,status.MPI_SOURCE);
    tmpconfig.mpiSend(status.MPI_SOURCE);
    MPI_Send(weight,status.MPI_SOURCE);
    //introduce completed point into the average
    //cout << "master: updating average " << endl;
    postavg.update_average(pt);
    totcount++;
    if(totcount%1000==0) cout << "Completed " << totcount << " walkers " << endl;
  }
  
  cout << "master: collecting final averages " << endl;
  
  //Loop through all the nodes and collect their last points, adding them in
  //Write out the final averages.
  for(int r=1; r < mpi_info.nprocs; r++) { 
    cout << "collecting from " << r << endl;
    int done;
    MPI_Recv(done,r);
    done=0;
    pt.mpiReceive(r);
    MPI_Send(done,r);
    postavg.update_average(pt);
  }
  cout << "printing " << endl;
  postavg.print(average_var,os);
  
#endif //USE_MPI
}
    

//----------------------------------------------------------------------
