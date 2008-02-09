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

#include "Vmc_method.h"
#include "qmc_io.h"
#include "ulec.h"
#include <iomanip>
#include "Program_options.h"
#include "System.h"
#include "Split_sample.h"
#include "Properties.h"

//----------------------------------------------------------------------

void Vmc_method::read(vector <string> words,
                      unsigned int & pos,
                      Program_options & options)
{

  have_read_options=1;

  ndim=3;

  if(!readvalue(words, pos=0, nblock, "NBLOCK"))
    error("Need NBLOCK in VMC section");
  if(!readvalue(words, pos=0, nconfig, "NCONFIG"))
    error("Need NCONFIG in METHOD section.");
  if(!readvalue(words, pos=0, nstep, "NSTEP"))
    error("Need NSTEP in VMC section");
  if(!readvalue(words, pos=0, timestep, "TIMESTEP"))
    error("Need TIMESTEP in VMC section");

  //optional options

  if(!readvalue(words, pos=0, ndecorr, "NDECORR"))
    ndecorr=2;

  if(readvalue(words, pos=0, storeconfig, "STORECONFIG"))
    canonical_filename(storeconfig);

  if(readvalue(words, pos=0, readconfig, "READCONFIG"))
    canonical_filename(readconfig);
  

  if(!readvalue(words, pos=0, log_label, "LABEL"))
    log_label="vmc";


  vector <string> proptxt;
  if(readsection(words, pos=0, proptxt, "PROPERTIES")) 
    mygather.read(proptxt);
  //myprop.read(proptxt, options.systemtext[0], options.twftext[0]);


  vector <string> dynsec;
  if(!readsection(words, pos=0, dynsec, "DYNAMICS")) 
    dynsec.push_back("SPLIT");

  allocate(dynsec, sampler);
  
  vector<string> tmp_dens;
  pos=0;
  while(readsection(words, pos, tmp_dens, "DENSITY")) {
    dens_words.push_back(tmp_dens);
  }
  
  pos=0;
  while(readsection(words, pos, tmp_dens, "AVERAGE")) {
    avg_words.push_back(tmp_dens);
  }
    
  
  if(readvalue(words, pos=0, config_trace, "CONFIG_TRACE"))
    canonical_filename(config_trace);


  pos=0;
  if(readvalue(words, pos, guidetype, "GUIDETYPE"))
  {
    if(guidetype=="SUMSQUARES")
    {
      guidewf=new Vmc_sum_squares;
    }
    else if(guidetype=="FIRST") {
      guidewf=new Primary;
    }
    else
    {
      error("I don't understand ", guidetype, " after GUIDETYPE.");
    }
  }
  else
  {
    guidetype="SUMSQUARES";
    guidewf=new Vmc_sum_squares;
  }

  print_wf_vals=0;
  if(haskeyword(words, pos=0, "PRINT_WF_VALS")) {
    print_wf_vals=1;
    cout << " enabling wf printing " << endl;
  }

  low_io=0;
  if(haskeyword(words,pos=0,"LOW_IO")) low_io=1;
 
}


//----------------------------------------------------------------------

int Vmc_method::generateVariables(Program_options & options) {

  
  if(!have_read_options) {
    error("Need to call Vmc_method::read() before generateVariables()");
  }
  have_generated_variables=1;

  debug_write(cout, "System\n");
  allocate(options.systemtext[0], sysprop );
  debug_write(cout, "Wave function\n");
  if(options.twftext.size() < 1)
    error("Need TRIALFUNC section for VMC");
  allocate(options.twftext[0], sysprop, mywfdata);
  debug_write(cout, "Pseudopotential\n");
  sysprop->generatePseudo(options.pseudotext, pseudo);
  
  densplt.Resize(dens_words.size());
  for(int i=0; i< densplt.GetDim(0); i++) {
    allocate(dens_words[i], sysprop, options.runid,densplt(i));
  }
  


  return 1;
}


//----------------------------------------------------------------------

int Vmc_method::allocateIntermediateVariables(System * locsys, 
                                              Wavefunction_data * locwfdata) {

  debug_write(cout, "temporary variables\n");
  locsys->generateSample(sample);
  locwfdata->generateWavefunction(wf);
  sample->attachObserver(wf);

  
  average_var.Resize(avg_words.size());
  average_var=NULL;
  for(int i=0; i< average_var.GetDim(0); i++) { 
    allocate(avg_words[i], locsys, locwfdata, average_var(i));
  }
  
  
  return 1;
}

//----------------------------------------------------------------------

int Vmc_method::deallocateIntermediateVariables() {
  if(wf) delete wf;
  wf=NULL;
  if(sample) delete sample;
  sample=NULL;
  
  for(int i=0; i< average_var.GetDim(0); i++) { 
    if(average_var(i)) delete average_var(i);
    average_var(i)=NULL;
  }
  return 1;
}

//----------------------------------------------------------------------

int Vmc_method::showinfo(ostream & os)
{


  if(os)
  {
    if(have_generated_variables) {
      sysprop->showinfo(os);
      os << endl << endl;
      os << "                               Wavefunction  "
         << endl << endl;
      mywfdata->showinfo(os);
      pseudo->showinfo(os);
      os << endl;
    }
    os << "Variational Monte Carlo:\n";
    os << "Configurations per processor: " <<  nconfig   << endl;
    os << "Number of processors: "        <<  mpi_info.nprocs << endl;
    os << "Total configurations: " <<          nconfig*mpi_info.nprocs << endl;
    os << "Blocks: " <<                        nblock    << endl;
    os << "Steps per block: " <<               nstep     << endl;
    os << "Number of decorrelation steps: " << ndecorr   << endl;
    os << "Timestep: " <<                      timestep  << endl;
    os << "Guiding wavefunction: " << guidetype << endl;
    string indent="  ";
    sampler->showinfo(indent, os);

    os << endl;
    return 1;
  }
  else
    return 0;
}

//----------------------------------------------------------------------

void Vmc_method::readcheck(string & filename) {
  int configsread=0;

  config_pos.Resize(nconfig);

  if(filename != "") { 

    ifstream checkfile(filename.c_str());
    if(!checkfile) 
      error("Couldn't open ", filename);
    
    long int is1, is2;
    string dummy;
    checkfile >> dummy;
    if(dummy != "RANDNUM") error("Expected RANDNUM in checkfile");
    checkfile >> is1 >> is2;
    //rng.seed(is1, is2);  


    while(checkfile >>dummy && configsread < nconfig) {
    //if(read_config(dummy, checkfile, electrons(configsread)))
    //  configsread++;
      if(read_config(dummy, checkfile, sample)) {
        config_pos(configsread++).savePos(sample);
      }
      
    }
    checkfile.close();
  }

  for(int i=configsread; i< nconfig; i++) {
    //electrons(i)->randomGuess();
    sample->randomGuess();
    config_pos(i).savePos(sample);
  }
  debug_write(cout, configsread, " configs read ", nconfig-configsread,
              " configs randomly generated \n");
}

//----------------------------------------------------------------------

void Vmc_method::storecheck(string & filename, int append) {
  if(filename=="") return;

  ofstream checkfile;
  if(append) checkfile.open(filename.c_str(), ios::app);
  else checkfile.open(filename.c_str());
  if(!checkfile) error("Couldn't open ", filename);
  checkfile.precision(15);
  
  long int is1, is2;
  rng.getseed(is1, is2);
  checkfile << "RANDNUM " << is1 << "  " << is2 << endl;
  for(int i=0; i< nconfig; i++) {
    checkfile << "SAMPLE_POINT { \n";
    //write_config(checkfile, electrons(i));
    config_pos(i).restorePos(sample);
    write_config(checkfile, sample);
    checkfile << "}\n\n";
  }

  checkfile.close();  

}


//----------------------------------------------------------------------

void Vmc_method::run(Program_options & options, ostream & output) {
  if(!have_generated_variables) 
    error("Must generate variables to use Vmc_method::run");
  
  string logfile=options.runid+".log";

  if(mpi_info.node==0 ) {
    ofstream logout(logfile.c_str(), ios::app);
    logout << "#-------------------------------------------------\n";
    logout << "#VMC run: timestep " << timestep 
           << " nconfig " << mpi_info.nprocs*nconfig
           << " steps " << nstep << endl;
    logout << "#-------------------------------------------------\n\n\n";
    logout.close();
  }

  myprop.setLog(logfile, log_label);
  runWithVariables(myprop, sysprop, mywfdata, pseudo,output);
}


/*!

*/
void Vmc_method::runWithVariables(Properties_manager & prop, 
                                  System * sys, 
                                  Wavefunction_data * wfdata,
                                  Pseudopotential * psp,
                                  ostream & output)
{

  nelectrons=sys->nelectrons(0)+sys->nelectrons(1);

  
  allocateIntermediateVariables(sys, wfdata);
  readcheck(readconfig);

  cout.precision(10);

  if(output) output << "Entering VMC" << endl;

  prop.setSize(wf->nfunc(), nblock, nstep, nconfig, sys, wfdata,
	       mygather.nAux());
  mygather.squareWeight(1);

  output.precision(10);
  prop.initializeLog(average_var);
  //our averaging variables
  Array1 <doublevar> avglifetime(nblock);
  unsigned int maxlife=0;
  avglifetime=0;
  Array2 <unsigned int> age(nconfig, nelectrons);
  age=0;
  Array1 <doublevar> diffusion_rate(nblock,0.0);

 
  for(int block=0; block< nblock; block++) {
    int nwf_guide=wf->nfunc();
    Dynamics_info dinfo;

    for(int walker=0; walker<nconfig; walker++) {  
      
      config_pos(walker).restorePos(sample);
      wf->notify(all_electrons_move,0);
     
      wf->updateLap(wfdata, sample);
      
      if(print_wf_vals) { 
        Wf_return wfval(nwf_guide,2);
        wf->getVal(wfdata,0, wfval);
        cout << "node " << mpi_info.node << "  amp " << wfval.amp(0,0) << endl;
      }
      
      for(int step=0; step< nstep; step++) {
        Array1 <doublevar> rotx(3), roty(3), rotz(3);
        generate_random_rotation(rotx, roty, rotz);
        psp->rotateQuadrature(rotx, roty, rotz);

        //------------------------------------------
        
        Array1 <doublevar> oldpos(3);
        Array1 <doublevar> newpos(3);
        for(int decorr=0; decorr< ndecorr; decorr++) {

            for(int e=0; e<nelectrons; e++) {
              sample->getElectronPos(e,oldpos);
              
              int acc=sampler->sample(e,sample, wf, 
                                 wfdata, guidewf,dinfo, timestep);
              sample->getElectronPos(e,newpos);
              
              for(int d=0; d< 3; d++) {
                diffusion_rate(block)+=(newpos(d)-oldpos(d))                  
                  *(newpos(d)-oldpos(d));
              }


              if(acc>0) {
                age(walker, e)=0;
              }
              else {
                age(walker,e)++;
                if(age(walker,e) > maxlife) maxlife=age(walker,e);
              }
              avglifetime(block)+=age(walker, e);
           }  //electron
        }  //decorrelation
          
        
        Properties_point pt;
        mygather.gatherData(pt, psp, sys, wfdata, wf, 
                            sample, guidewf);
        
        for(int i=0; i< densplt.GetDim(0); i++)
          densplt(i)->accumulate(sample,1.0);

        pt.avgrets.Resize(average_var.GetDim(0));
        for(int i=0; i< average_var.GetDim(0); i++) { 
          average_var(i)->evaluate(wfdata, wf, sys, sample, pt.avgrets(i));
        }
        pt.parent=walker;
        pt.nchildren=1; pt.children(0)=1;
        prop.insertPoint(step, walker, pt);
        
        //This may screw up if we have >1 walker!
        if(config_trace!="" && block >0) {
          if(nconfig !=1) error("trace only works with nconfig=1");
          config_pos(walker).savePos(sample);
          storecheck(config_trace,1);
        }
        
      }   //step
      
      config_pos(walker).savePos(sample);
      
    }   //walker

    prop.endBlock();
    
    if(!low_io || block==nblock-1) { 
      storecheck(storeconfig);
      for(int i=0; i< densplt.GetDim(0); i++)
        densplt(i)->write();
    }



    //Output for the block
    if(output) {
      if(global_options::rappture ) { 
        ofstream rapout("rappture.log");
        rapout << "=RAPPTURE-PROGRESS=>" << int(100.0*doublevar(block+1)/doublevar(nblock))
               << "  Sampling wave function" << endl;
        cout << "=RAPPTURE-PROGRESS=>" << int(100.0*doublevar(block+1)/doublevar(nblock))
               << "  Sampling wave function" << endl;
        rapout.close();
      }
      output << "****Block " << block 
             << " avg life " << avglifetime(block)/(nstep*nconfig*ndecorr*nelectrons)
             << " max life " << double(maxlife)/double(ndecorr)
             << endl;

      sampler->showStats(output);
      sampler->resetStats();
      prop.printBlockSummary(output);
      output << endl;
    }
  }              //blocks done



  if(output) {
    output << "-----------------Run Ended------------\n" << endl;

    prop.printSummary(output, average_var);
    output << endl;
    output << "Maximum lifetime " << double(maxlife)/double(ndecorr)
           << " steps "<<  endl;
    doublevar avglife=0;
    for(int b=0; b< nblock; b++) {
      avglife+=avglifetime(b)/(nblock*nstep*nconfig*ndecorr*ndecorr*nelectrons);
    }

    output << "Average lifetime : " << avglife << " steps " << endl;
    output << "VMC Done. \n";
  }
 
  wfdata->clearObserver();
  deallocateIntermediateVariables();
}

//------------------------------------------------------------------------

