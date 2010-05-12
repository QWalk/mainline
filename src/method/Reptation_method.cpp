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


#include "Reptation_method.h"
#include "qmc_io.h"
#include "ulec.h"
#include <iomanip>
#include "Program_options.h"
#include "System.h"
#include "Split_sample.h"
#include "Properties.h"



//----------------------------------------------------------------------


struct Reptile_point {

  Reptile_point() { 
    age=0;
  }
  Properties_point prop;

  Array1 < Array1 <doublevar> > electronpos;

  doublevar age;
  doublevar branching; //!< branching weight

  Array2 <doublevar> drift;


  //--------------------------------------------------
  void savePos(Sample_point * sample) {
    int nelectrons=sample->electronSize();
    electronpos.Resize(nelectrons);
    for(int i=0; i< nelectrons; i++) {
      electronpos(i).Resize(3);
      sample->getElectronPos(i,electronpos(i));
    }
  }
  //--------------------------------------------------

  void restorePos(Sample_point * sample) {
    int nelectrons=sample->electronSize();
    assert(nelectrons==electronpos.GetDim(0));
    for(int i=0; i< nelectrons; i++) 
      sample->setElectronPos(i,electronpos(i));
  }

  //--------------------------------------------------
  void write(string & indent, ostream & os) {
    int nelectrons=electronpos.GetDim(0);
    os << indent <<  "age " << age << endl;
    os << indent << "branching " << branching << endl;
	

    os << indent << "numElectrons " << nelectrons <<  endl;
    for(int e=0; e< nelectrons; e++) {
      os << indent;
      for(int d=0; d< 3; d++) 
        os << electronpos(e)(d) << "   ";
      os << endl;
    }
    os << indent << "drift " << endl;
    for(int e=0; e< nelectrons; e++) {
      os << indent;
      for(int d=0; d< 3; d++) 
        os << drift(e,d) << "   ";
      os << endl;
    }

    string indent2=indent+"  ";
    os << indent << "Properties_point  {" << endl;
    prop.write(indent2,os);
    os << indent << "}" << endl;

  }
  //--------------------------------------------------
  void read(istream & is) { 
    string dummy; int nelectrons;
    const char *errmsg="misformatting in checkpoint read";
    
    is >> dummy >> age;
    if(dummy != "age") error(errmsg);
    is >> dummy >> branching;
    if(dummy != "branching") 
      error("expected branching, got ", dummy, " This probably means that your"
            " config files are too old.  Sorry!");
    is >> dummy >> nelectrons;
    if(dummy!="numElectrons") error(errmsg);
    
    electronpos.Resize(nelectrons);
    for(int e=0; e < nelectrons; e++) {
      read_array(is, 3, electronpos(e));
    }

    is >> dummy;
    if(dummy!="drift") error(errmsg);
    drift.Resize(nelectrons,3);
    for(int e=0; e< nelectrons; e++) { 
      for(int d=0; d< 3; d++) { 
	is >> drift(e,d);
      }
    }

    
    is >> dummy; 
    if(dummy!="Properties_point") error(errmsg);
    is >> dummy;
    prop.read(is);
    is >> dummy;
  }
  //--------------------------------------------------

};

//----------------------------------------------------------------------


void Reptation_method::read(vector <string> words,
                      unsigned int & pos,
                      Program_options & options)
{

  have_read_options=1;

  if(!readvalue(words, pos=0, nblock, "NBLOCK"))
    error("Need NBLOCK in RETPTATION section");

  if(!readvalue(words, pos=0, nstep, "NSTEP"))
    error("Need NSTEP in REPTATION section");

  if(!readvalue(words, pos=0, timestep, "TIMESTEP"))
    error("Need TIMESTEP in REPTATION section");

  if(!readvalue(words, pos=0, reptile_length, "LENGTH"))
    error("Need LENGTH in REPTATION");


  if(!readvalue(words, pos=0, readconfig, "READCONFIG"))
    error("Need READCONFIG in REPTATION"); 
  canonical_filename(readconfig);



  //------------------optional stuff

  vector <string> proptxt;
  if(readsection(words, pos=0, proptxt, "PROPERTIES")) 
    mygather.read(proptxt);
  
  if(!readvalue(words, pos=0, log_label, "LABEL"))
    log_label="rmc";

  if(readvalue(words, pos=0, storeconfig, "STORECONFIG"))
    canonical_filename(storeconfig);

  if(readvalue(words, pos=0, center_trace, "CENTER_TRACE"))
    canonical_filename(center_trace);

  if(!readvalue(words, pos=0, trace_wait, "TRACE_WAIT"))
    trace_wait=int(.3/timestep)+1;

    
  vector <string> dynamics_words;
  if(!readsection(words, pos=0, dynamics_words, "DYNAMICS") ) 
    dynamics_words.push_back("SPLIT");

  allocate(dynamics_words, sampler);

  
  vector<string> tmp_dens;
  pos=0;
  while(readsection(words, pos, tmp_dens, "DENSITY")) {
    dens_words.push_back(tmp_dens);
  }
  
  pos=0;
  while(readsection(words, pos, tmp_dens, "AVERAGE")) {
    avg_words.push_back(tmp_dens);
  }
  

  sampler->enforceNodes(1);

  guidewf=new Primary;
}


//----------------------------------------------------------------------

int Reptation_method::generateVariables(Program_options & options) {

  
  if(!have_read_options) {
    error("Need to call Reptation_method::read() before generateVariables()");
  }
  have_generated_variables=1;

  debug_write(cout, "properties\n");
  allocate(options.systemtext[0], sys );
  allocate(options.twftext[0], sys, mywfdata);
  debug_write(cout, "Pseudopotential\n");

  sys->generatePseudo(options.pseudotext, pseudo);

  densplt.Resize(dens_words.size());
  for(int i=0; i< densplt.GetDim(0); i++) {
    allocate(dens_words[i], sys, options.runid,densplt(i));
  }


  return 1;
}


//----------------------------------------------------------------------

int Reptation_method::allocateIntermediateVariables(System * locsys, 
                                              Wavefunction_data * locwfdata) {


  //-----------------------------------
  //Sample point initialization
  debug_write(cout, "electrons\n");

  sample=NULL;
  locsys->generateSample(sample);

  wf=NULL;
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

int Reptation_method::deallocateIntermediateVariables() {

  if(wf) delete wf;
  if(sample) delete sample;
  wf=NULL;
  sample=NULL;
  
  for(int i=0; i< average_var.GetDim(0); i++) { 
    if(average_var(i)) delete average_var(i);
    average_var(i)=NULL;
  }
  
  return 1;
}

//----------------------------------------------------------------------

int Reptation_method::showinfo(ostream & os) {

  if(os){
    if(have_generated_variables) {
      sys->showinfo(os);
      os << endl << endl;
      os << "                               Wavefunction  "
         << endl << endl;
      mywfdata->showinfo(os);
      pseudo->showinfo(os);
      os << endl;
    }
    os << "Reptation settings:\n";
    os << "Number of processors: "        <<  mpi_info.nprocs << endl;
    os << "Blocks: " <<                        nblock    << endl;
    os << "Steps per block: " <<               nstep     << endl;
    os << "Timestep: " <<                      timestep  << endl;
    os << "Reptile size: " <<                  reptile_length << endl;
    string indent="  ";
    os << "Dynamics generator:" << endl;
    sampler->showinfo(indent, os);
    os << endl;
    return 1;
  }
  else
    return 0;
}

//----------------------------------------------------------------------



int Reptation_method::readcheck(string & filename, 
                                 int & direction, 
                                 deque <Reptile_point> & reptile) {
  int configsread=0;

  if(filename == "") return 0;

  ifstream checkfile(filename.c_str());
  if(!checkfile) 
    error("Couldn't open config file ", filename);

  long int is1, is2;
  string dummy;
  checkfile >> dummy;
  if(dummy != "RANDNUM") error("Expected RANDNUM in checkfile");
  checkfile >> is1 >> is2;
  rng.seed(is1, is2);  


  while(checkfile >>dummy) {
    
    //------------------------------------------------------
    //Start from a (single) VMC configuration
    if(dummy=="SAMPLE_POINT") { 
      while(checkfile >> dummy)
        if(read_config(dummy, checkfile, sample))
          configsread++;
      if(configsread < 1) 
        error("Need at least one configuration to start reptation");
      checkfile.close();
      return 0;
    }
    
    //------------------------------------------------------
    //Read in the reptile if we're restarting from a previous RMC run
    else if(dummy=="REPTILE") {
      //cout << "reading reptile" << endl;
      checkfile >> dummy; 
      if(dummy != "{") error("expected {, got ", dummy);
      
      checkfile >> dummy;
      if(dummy != "direction") error("expected direction, got ",dummy);
      checkfile >> direction;
      
      checkfile >> dummy;
      if(dummy != "length") error("expected length, got ", dummy);
      int length;
      checkfile >> length;
      if(length != reptile_length)
        error("can't use a reptile of different length..");
      //cout << "reptile direction " << direction 
      //    << " length " << length << endl;
      Reptile_point pt;
      for(int i=0; i< length; i++) {
        
        checkfile >> dummy;
        if(dummy != "Reptile_point") error("expected Reptile_point, got ",dummy);
        checkfile >> dummy; //open brace;
        pt.read(checkfile);
        checkfile >> dummy; //close brace
        if(dummy != "}") error("expected }, got ", dummy);
        reptile.push_back(pt);
      }
      return 1;
    }
  }
  return 0;
}

//----------------------------------------------------------------------
 
void Reptation_method::storecheck(int direction, 
                                  deque <Reptile_point> & reptile,
                                  string & filename) {
  if(filename=="") return;
  ofstream checkfile(filename.c_str());
  if(!checkfile) error("Couldn't open ", filename);
  checkfile.precision(15);
  
  long int is1, is2;
  rng.getseed(is1, is2);
  checkfile << "RANDNUM " << is1 << "  " << is2 << endl;
  checkfile << "REPTILE { " << endl;
  checkfile << "direction " << direction << endl;
  checkfile << "length "  << reptile.size() << endl;
  string indent="  ";
  for(deque<Reptile_point>::iterator i=reptile.begin();
      i!=reptile.end(); i++) {
    checkfile << "Reptile_point { \n";
    i->write(indent,checkfile);
    checkfile << "} " << endl;
  }

  checkfile << "}\n";
  
  checkfile.close();  
}


//----------------------------------------------------------------------

void Reptation_method::run(Program_options & options, ostream & output) {
  if(!have_generated_variables) 
    error("Must generate variables to use Reptation_method::run");
  
  string logfile=options.runid+".log";

  if(mpi_info.node==0 ) {
    ofstream logout(logfile.c_str(), ios::app);
    logout << "#-------------------------------------------------\n";
    logout << "#Reptation run: timestep " << timestep 
           << " steps " << nstep << endl;
    logout << "#-------------------------------------------------\n\n\n";
    logout.close();
  }

  Properties_manager myprop;
  myprop.setLog(logfile, log_label);
  runWithVariables(myprop, sys, mywfdata, pseudo,output);
}


//----------------------------------------------------------------------


void Reptation_method::get_avg(deque <Reptile_point> & reptile, 
			       Properties_point & pt) {
  int size=reptile.size();
  int nwf=reptile[0].prop.kinetic.GetDim(0);
  if(nwf >1) error("nwf > 0 not supported yet");

  pt.setSize(nwf);
  
  Reptile_point & last(reptile[size-1]);

  //How to do averaging at either end.  Not doing this right
  //now because of correlated sampling..if we really want energies,
  //usually DMC is a better choice.
  pt=last.prop;
  Reptile_point & first(reptile[0]);  
  //pt.kinetic(0)=.5*(first.prop.kinetic(0)+last.prop.kinetic(0));
  //pt.nonlocal(0)=.5*(first.prop.nonlocal(0)+last.prop.nonlocal(0));
  //pt.potential(0)=.5*(first.prop.potential(0) + last.prop.potential(0));
  pt.count=1;
  pt.weight=1;
  
}


//---------------------------------------------------------------------

void Reptation_method::get_center_avg(deque <Reptile_point> & reptile, 
				      Properties_point & pt) { 
  int nwf=reptile[0].prop.kinetic.GetDim(0);
  int size=reptile.size();
  if(nwf >1) error("nwf > 0 not supported yet");

  pt.setSize(nwf);
  int num=size/2+1;
  pt=reptile[num].prop;

  pt.count=1;

  pt.weight=1;
}

//----------------------------------------------------------------------


doublevar Reptation_method::slither(int direction,
                                    deque <Reptile_point> & reptile,
                                    Properties_gather & mygather, 
                                    Reptile_point & pt, 
                                    doublevar & main_diffusion) {
  
  Dynamics_info dinfo;
  int nelectrons=sample->electronSize();
  
  for(int e=0; e<nelectrons; e++) {
    sampler->sample(e,sample, wf, 
                    mywfdata, guidewf, dinfo, timestep);

  }
  //Update all the quantities we want
  mygather.gatherData(pt.prop, pseudo, sys, mywfdata, wf, sample, guidewf);


  pt.savePos(sample);

  
  doublevar nen=pt.prop.energy(0);
  doublevar olden=0.0;
  deque<Reptile_point>::iterator last_point;
  if(reptile.size() >0) {
    if(direction==1 ) 
      last_point=(reptile.end()-1);
    else 
      last_point=reptile.begin();
      olden=last_point->prop.energy(0);
  }

  //Cutting off the energy
  doublevar fbet=max(eref-nen,eref-olden);
  doublevar cutoff=1;
  if(fbet > 1.5*energy_cutoff) 
    cutoff=0;
  else if(fbet > energy_cutoff) 
    cutoff*=(1.-(fbet-energy_cutoff)/(.5*energy_cutoff));

  pt.branching=-0.5*timestep*cutoff*(olden+nen);

  
  
  //Acceptance..
  doublevar lost_branching=0;
  if(reptile.size()>0) { 
    if(direction==1 )
      lost_branching=(reptile.begin()+1)->branching;
    else
      lost_branching=(reptile.end()-1)->branching;
  }
  doublevar accept=exp(pt.branching-lost_branching);

  if(pt.prop.wf_val.amp(0,0) < -1e12 ) {
    cout  << "rejecting based on wf value" << endl;
    accept=0;
  }
  
  pt.prop.avgrets.Resize(1,average_var.GetDim(0));
  for(int i=0; i< average_var.GetDim(0); i++) { 
    average_var(i)->evaluate(mywfdata, wf, sys, sample, pt.prop.avgrets(0,i));
  }
  
  return accept;
}


//----------------------------------------------------------------------

/*!

*/
void Reptation_method::runWithVariables(Properties_manager & prop, 
                                  System * sys, 
                                  Wavefunction_data * wfdata,
                                  Pseudopotential * psp,
                                  ostream & output)
{


  
  int nelectrons=sys->nelectrons(0)+sys->nelectrons(1);
  allocateIntermediateVariables(sys, wfdata);
  

  
  prop.setSize(wf->nfunc(), nblock, nstep, 1, 
               sys, wfdata);
  prop.initializeLog(average_var);

  Properties_manager prop_center;
  string logfile, label_temp;
  prop.getLog(logfile, label_temp);
  label_temp+="_cen";
  prop_center.setLog(logfile, label_temp);
  
  prop_center.setSize(wf->nfunc(), nblock, nstep, 1, sys, 
                      wfdata);
  prop_center.initializeLog(average_var);


  cout.precision(10);
  output.precision(10);

  Sample_point * center_samp(NULL);
  sys->generateSample(center_samp);
  
  deque <Reptile_point> reptile;
  int direction=1;
  Reptile_point pt;
  Dynamics_info dinfo;
  //Generate a new reptile if we're starting from 
  // a single VMC point.
  if(!readcheck(readconfig, direction, reptile)) {
    
    wf->notify(all_electrons_move, 0);
    wf->updateLap(wfdata, sample);
    
    for(int i=0; i< reptile_length; i++) {
      doublevar main_diffusion;
      slither(direction,reptile, mygather,pt,main_diffusion);
      reptile.push_back(pt);
    }
  }
   


  assert(reptile.size()==reptile_length);
  //Branch limiting variables
  //we start off with no limiting, and establish the parameters after the
  //first block.  This seems to be reasonably stable, since it's mostly
  //to keep the reptile from getting stuck.
  eref=0;
  energy_cutoff=1e16;

  //--------begin averaging..
  
  Array3 <doublevar> derivatives_block(nblock, sys->nIons(), 3);
  
  for(int block=0; block< nblock; block++) {

    //clock_t block_start_time=clock();
    doublevar avg_age=0;
    doublevar max_age=0;

    doublevar main_diff=0;
    double ntry=0, naccept=0;
    double nbounce=0;

    //Control variable that will be set to one when 
    //we change direction, which signals to recalculate
    //the wave function
    int recalc=1;

    for(int step=0; step< nstep; step++) {
      psp->randomize();
      

      if(recalc) { 
        if(direction==1) 
          reptile[reptile_length-1].restorePos(sample);
        else
          reptile[0].restorePos(sample);          
      }

      
      doublevar main_diffusion;
      doublevar accept=slither(direction, reptile,mygather, pt,
                               main_diffusion);

      ntry++;
      if(accept+rng.ulec() > 1.0) {
        recalc=0;
        naccept++;
        
        
        
        main_diff+=main_diffusion;
        if(direction==1) {
          reptile.pop_front();
          reptile.push_back(pt);
        }
        else {
          reptile.pop_back();
          reptile[0].branching=pt.branching;
          reptile.push_front(pt);
        }
      }
      else {
        recalc=1;
        direction*=-1;
        nbounce++;
      }

      for(deque<Reptile_point>::iterator i=reptile.begin();
          i!=reptile.end(); i++) {
        i->age++;
        avg_age+=i->age/reptile_length;
        if(i->age > max_age) max_age=i->age;
      }

      
      Properties_point avgpt;
      get_avg(reptile, avgpt);
      avgpt.parent=0; avgpt.nchildren=1; //just one walker
      avgpt.children(0)=0;
      prop.insertPoint(step, 0, avgpt);
      
      int cpt=reptile_length/2+1;      
      Properties_point centpt;
      get_center_avg(reptile, centpt);
      centpt.parent=0; centpt.nchildren=1;
      centpt.children(0)=0;
      
      prop_center.insertPoint(step, 0, centpt);
      
      reptile[cpt].restorePos(center_samp);

      for(int i=0; i< densplt.GetDim(0); i++) 
        densplt(i)->accumulate(center_samp,1.0);
      
      
      if(center_trace != "" 
         && (block*nstep+step)%trace_wait==0) {
        ofstream checkfile(center_trace.c_str(), ios::app);
        if(!checkfile)error("Couldn't open ", center_trace);
        checkfile << "SAMPLE_POINT { \n";
        write_config(checkfile, sample);
        checkfile << "}\n\n";
      }
      
      
    }   //step

    prop.endBlock();
    prop_center.endBlock();
    double ntot=parallel_sum(nstep);

    Properties_block lastblock;
    prop.getLastBlock(lastblock);
    eref=lastblock.avg(Properties_types::total_energy,0);
    energy_cutoff=10*sqrt(lastblock.var(Properties_types::total_energy,0));
    
    
    nbounce=parallel_sum(nbounce);
    naccept=parallel_sum(naccept);
    ntry=parallel_sum(ntry);
    avg_age=parallel_sum(avg_age);

    for(int i=0; i< densplt.GetDim(0); i++) 
      densplt(i)->write();

    storecheck(direction, reptile, storeconfig);
    main_diff=parallel_sum(main_diff);
    if(output) {
      output << "****Block " << block 
             << " acceptance " << naccept/ntry 
             << "  average distance before bounce " << ntot/nbounce
             << endl;
      output << "average age " << avg_age/ntot 
             << "   max age " << max_age <<  endl;
      output << "eref " << eref << " cutoff " << energy_cutoff << endl;
      output << "Green's function sampler:" << endl;
      sampler->showStats(output);
      prop.printBlockSummary(output);
      output << "Center averaging: " << endl;
      prop_center.printBlockSummary(output);
    }
    sampler->resetStats();

    //clock_t block_end_time=clock();
    
    //cout << mpi_info.node << ":CPU block time " 
    //// << double(block_end_time-block_start_time)/double(CLOCKS_PER_SEC)
    // << endl;

  }   //block


  if(output) {
    output << "############## Reptation Done ################\n";
    output << "End averages " << endl;
    prop.printSummary(output,average_var);
    output << "Center averages " << endl;
    prop_center.printSummary(output,average_var);

  }


  delete center_samp;
  wfdata->clearObserver();
  deallocateIntermediateVariables();
}

//------------------------------------------------------------------------

