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

  readvalue(words, pos=0, storeconfig, "STORECONFIG");
  readvalue(words, pos=0, readconfig, "READCONFIG");
  

  if(!readvalue(words, pos=0, log_label, "LABEL"))
    log_label="vmc";


  vector <string> proptxt;
  if(readsection(words, pos=0, proptxt, "PROPERTIES")) 
    mygather.read(proptxt);
  //myprop.read(proptxt, options.systemtext[0], options.twftext[0]);


  vector <string> dynsec;
  if(!readsection(words, pos=0, dynsec, "DYNAMICS")) 
    dynsec.push_back("SPLIT");

  int nsys=options.systemtext.size();
  sampler.Resize(nsys);
  sampler=NULL;
  for(int i=0; i< nsys; i++) 
    allocate(dynsec, sampler(i));
  
  vector<string> tmp_dens;
  pos=0;
  while(readsection(words, pos, tmp_dens, "DENSITY")) {
    dens_words.push_back(tmp_dens);
  }

  pos=0;
  while(readsection(words, pos, tmp_dens, "NONLOCAL_DENSITY")) {
    nldens_words.push_back(tmp_dens);
  }

  pos=0;
  while(readsection(words, pos, tmp_dens, "AVERAGE")) {
    avg_words.push_back(tmp_dens);
  }
    

  if(readvalue(words, pos=0, config_trace, "CONFIG_TRACE"))
    canonical_filename(config_trace);


  pos=0;
  if(readvalue(words, pos, guidetype, "GUIDETYPE")) {
    if(guidetype=="SUMSQUARES") 
      guidewf=new Vmc_sum_squares;
    else if(guidetype=="FIRST") 
      guidewf=new Primary;
    else error("I don't understand ", guidetype, " after GUIDETYPE.");
  }
  else {
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
  int nsys=options.systemtext.size();
  if(options.twftext.size()!=nsys)  
    error("Need the same number of SYSTEM and TRIALFUNC sections in VMC");
  sys.Resize(nsys);
  wfdata.Resize(nsys);
  sys=NULL;
  wfdata=NULL;


  for(int s=0; s < nsys; s++) {  
    allocate(options.systemtext[s], sys(s) );
    allocate(options.twftext[s], sys(s), wfdata(s));
  }
  debug_write(cout, "Pseudopotential\n");
  sys(0)->generatePseudo(options.pseudotext, pseudo);
  
  densplt.Resize(dens_words.size());
  for(int i=0; i< densplt.GetDim(0); i++) {
    allocate(dens_words[i], sys(0), options.runid,densplt(i));
  }
  nldensplt.Resize(nldens_words.size());
  for(int i=0; i< nldensplt.GetDim(0); i++) {
    allocate(nldens_words[i], sys(0), options.runid,nldensplt(i));
  }

  return 1;
}


//----------------------------------------------------------------------

int Vmc_method::allocateIntermediateVariables() { 

  debug_write(cout, "temporary variables\n");
  int nsys=sys.GetDim(0);
  sample.Resize(nsys);
  wf.Resize(nsys);
  sample=NULL;
  wf=NULL;
  for(int s=0; s< nsys; s++) { 
    sys(s)->generateSample(sample(s));
    wfdata(s)->generateWavefunction(wf(s));
    sample(s)->attachObserver(wf(s));
  } 
  
  average_var.Resize(nsys,avg_words.size());
  average_var=NULL;
  for(int i=0; i< average_var.GetDim(1); i++) { 
    for(int s=0;s < nsys; s++) { 
      allocate(avg_words[i], sys(s), wfdata(s), average_var(s,i));
    }
  }
  
  
  return 1;
}

//----------------------------------------------------------------------

int Vmc_method::deallocateIntermediateVariables() {
  int nsys=wf.GetDim(0);
  for(int s=0; s< nsys; s++) {
    if(wf(s)) delete wf(s);
    wf(s)=NULL;
    if(sample(s)) delete sample(s);
    sample(s)=NULL;
  }
  
  for(int i=0; i< average_var.GetDim(1); i++) { 
    for(int s=0; s< nsys; s++) { 
      if(average_var(s,i)) delete average_var(s,i);
      average_var(s,i)=NULL;
    }
  }
  return 1;
}

//----------------------------------------------------------------------

int Vmc_method::showinfo(ostream & os)
{


  if(os)
  {
    if(have_generated_variables) {
      sys(0)->showinfo(os);
      os << endl << endl;
      os << "                               Wavefunction  "
         << endl << endl;
      wfdata(0)->showinfo(os);
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
    sampler(0)->showinfo(indent, os);

    os << endl;
    return 1;
  }
  else
    return 0;
}

//----------------------------------------------------------------------

void Vmc_method::readcheck(string & filename) {
  config_pos.Resize(0);
  if(filename!="") { 
    read_configurations(filename, config_pos);
    if(config_pos.GetDim(0) != nconfig) { 
      error("Number of configurations is different in ",filename);
    /*
    Array1 <Config_save_point> tmpconfig=config_pos;
    config_pos.Resize(nconfig);
    for(int i=0; i< tmpconfig.GetDim(0); i++) { config_pos(i)=tmpconfig(i);} 
    for(int i=tmpconfig.GetDim(0); i< nconfig; i++) {
      sample->randomGuess();
      config_pos(i).savePos(sample);
    }
    */
    }   
  }
  else { 
    config_pos.Resize(nconfig);
    int nsys=sys.GetDim(0);
    if(nconfig%nsys!=0) error("nconfig needs to be a constant multiple of nsys");
    for(int i=0; i< nconfig; i++) { 
      int s=0;
      config_pos(i).system=s;
      sample(s)->randomGuess();
      config_pos(i).configs.savePos(sample(s));
    }
  }
}

//----------------------------------------------------------------------

void Vmc_method::storecheck(string & filename, int append) {
  
  if(filename=="") return;
  write_configurations(filename, config_pos);
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
  runSample(myprop, output);
}

//----------------------------------------------------------------------

/*!

*/
void Vmc_method::runWithVariables(Properties_manager & prop, 
                                  System * sys_, 
                                  Wavefunction_data * wfdata_,
                                  Pseudopotential * pseudo_,
                                  ostream & output) { 
  int nsys=1;
  sys.Resize(nsys);
  wfdata.Resize(nsys);
  sys(0)=sys_;
  wfdata(0)=wfdata_;
  pseudo=pseudo_;
  runSample(prop,output);
  sys=NULL;
  pseudo=NULL;
  wfdata=NULL;

  
}
//----------------------------------------------------------------------


void Vmc_method::runSample(Properties_manager & prop,
                                  ostream & output)
{

  nelectrons=sys(0)->nelectrons(0)+sys(0)->nelectrons(1);
  int nsys=sys.GetDim(0);
  allocateIntermediateVariables();
  readcheck(readconfig);
  cout.precision(10);
  
  if(output) output << "Entering VMC" << endl;

  prop.setSize(nsys, nblock, nstep, nconfig, sys(0), wfdata(0));
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
    int nwf_guide=wf(0)->nfunc();
    Dynamics_info dinfo;
    for(int walker=0; walker<nconfig; walker++) {  
      int s=config_pos(walker).system; 
      config_pos(walker).configs.restorePos(sample(s));
      wf(s)->notify(all_electrons_move,0);
     
      wf(s)->updateLap(wfdata(s), sample(s));
      if(print_wf_vals) { 
        Wf_return wfval(nwf_guide,2);
        wf(s)->getVal(wfdata(s),0, wfval);
        cout << "node " << mpi_info.node << "  amp " << wfval.amp(0,0) 
          << " phase " << cos(wfval.phase(0,0)) << endl;
      }
      for(int step=0; step< nstep; step++) {
        Array1 <doublevar> rotx(3), roty(3), rotz(3);
        generate_random_rotation(rotx, roty, rotz);
        pseudo->rotateQuadrature(rotx, roty, rotz);

        //------------------------------------------
        
        Array1 <doublevar> oldpos(3);
        Array1 <doublevar> newpos(3);
        for(int decorr=0; decorr< ndecorr; decorr++) {

            for(int e=0; e<nelectrons; e++) {
              sample(s)->getElectronPos(e,oldpos);
              int acc=sampler(s)->sample(e,sample(s), wf(s), 
                                 wfdata(s), guidewf,dinfo, timestep);
              sample(s)->getElectronPos(e,newpos);
              
              for(int d=0; d< 3; d++) {
                diffusion_rate(block)+=(newpos(d)-oldpos(d))                  
                  *(newpos(d)-oldpos(d));
              }
              if(print_wf_vals) { 
                Wf_return wfval(nwf_guide, 2);
                wf(s)->getVal(wfdata(s),0,wfval);
                cout << "step " << e << " amp " << wfval.amp(0,0) 
                  << " phase " << cos(wfval.phase(0,0)) << endl;
                cout << "pos " << newpos(0) << " " << newpos(1) << " " 
                  << newpos(2) << endl;
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
        pt.setSize(nsys);
        pt.avgrets.Resize(average_var.GetDim(0),average_var.GetDim(1));
        Array1 <doublevar> jacobian(nsys,1.0);
        Array1 <Wf_return> wf_value(nsys);
        for(int si=0; si< nsys; si++) { 
          if(s!=si) { 
            jacobian(si)=warper.warp_all(sample(s),sample(si));
            wf(si)->updateLap(wfdata(si), sample(si));
          }
          wf_value(si).Resize(wf(si)->nfunc(),2);
          wf(si)->getVal(wfdata(si),0,wf_value(si));
          //We'll only use the first average_var in order to preserve
          //correlated sampling.  Note that in some corner cases, this
          //may give the wrong behavior.
          for(int a=0; a < average_var.GetDim(1); a++) { 
            average_var(0,a)->randomize(wfdata(si),wf(si),sys(si),sample(si));
          }
        }
        int nrandvar=pseudo->nTest();
        Array1 <doublevar> rand_num(nrandvar);
        for(int i=0; i< nrandvar; i++) rand_num(i)=rng.ulec();
        for(int si=0; si < nsys; si++) {
          Array1 <doublevar> kinetic(wf(si)->nfunc()),nonlocal(wf(si)->nfunc());
          sys(si)->calcKinetic(wfdata(si),sample(si),wf(si),kinetic);
          pt.kinetic(si)=kinetic(0);
          pt.potential(si)=sys(si)->calcLoc(sample(si));
          pseudo->calcNonlocWithTest(wfdata(si), sys(si), sample(si), 
              wf(si),rand_num,nonlocal);
          pt.nonlocal(si)=nonlocal(0);
          //jacobian_save(i)=tot_jacobian;
          for(int a=0; a < average_var.GetDim(1); a++) { 
            average_var(0,a)->evaluate(wfdata(si),wf(si),sys(si),sample(si),pt,pt.avgrets(si,a));
          }
          
          pt.weight(si)=jacobian(si)*exp(2*(wf_value(si).amp(0,0)-wf_value(s).amp(0,0)));
        }
          

        for(int i=0; i< densplt.GetDim(0); i++)
          densplt(i)->accumulate(sample(0),1.0);
        for(int i=0; i< nldensplt.GetDim(0); i++)
          nldensplt(i)->accumulate(sample(0),1.0,wfdata(0),wf(0));
        pt.parent=walker;
        pt.nchildren=1; pt.children(0)=1;
        prop.insertPoint(step, walker, pt);
        
        //This may screw up if we have >1 walker!
        if(config_trace!="" && block >0) {
          if(nconfig !=1) error("trace only works with nconfig=1");
          config_pos(walker).configs.savePos(sample(s));
          storecheck(config_trace,1);
        }
        
      }   //step
      config_pos(walker).configs.savePos(sample(s));
    }   //walker
    prop.endBlock();
    if(!low_io || block==nblock-1) { 
      storecheck(storeconfig);
      for(int i=0; i< densplt.GetDim(0); i++)
        densplt(i)->write();
      for(int i=0; i< nldensplt.GetDim(0); i++)
        nldensplt(i)->write(log_label);
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

      for(int s=0; s< nsys; s++) { 
        sampler(s)->showStats(output);
        sampler(s)->resetStats();
      }
      prop.printBlockSummary(output);
      output << endl;
    }
  }              //blocks done


  if(output) {
    output << "-----------------Run Ended------------\n" << endl;

    //prop.printSummary(output, average_var);
    output << endl;
    output << "Maximum lifetime " << double(maxlife)/double(ndecorr)
           << " steps "<<  endl;
    doublevar avglife=0;
    for(int b=0; b< nblock; b++) {
      avglife+=avglifetime(b)/(nblock*nstep*nconfig*ndecorr*nelectrons);
    }

    output << "Average lifetime : " << avglife << " steps " << endl;
    output << "VMC Done. \n";
  }
  cout << "done output " << endl; 
  for(int s=0; s< nsys; s++) {
    wfdata(s)->clearObserver();
  }
  cout << "deallocate " << endl;
  deallocateIntermediateVariables();
}

//------------------------------------------------------------------------
//######################################################################

void Vmc_point::mpiSend(int node) { 
  configs.mpiSend(node);
  MPI_Send(system,node);
}

void Vmc_point::mpiReceive(int node) { 
  configs.mpiReceive(node);
  MPI_Recv(system,node);
}

void Vmc_point::read(istream & is) { 
  configs.read(is);
  long int filepos=is.tellg();
  string dum; is >> dum;
  if(!caseless_eq(dum,"VMC")) { 
    is.seekg(filepos);
    return;
  }
  else system=0;
  is >> dum;
  is >> dum >> system;
  //is >> dum; //clear the } 

}

void Vmc_point::write(ostream & os) { 
  string indent="";
  configs.write(os);
  os << "VMC  { \n";
     os << "system " << system << endl;
  os << " } \n";
}
//######################################################################


