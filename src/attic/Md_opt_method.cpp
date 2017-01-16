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
//------------------------------------------------------------------------
//src/Md_opt_opt_method.cpp
//Author:Lucas Wagner



#include "Md_opt_method.h"
#include "qmc_io.h"

#include <iomanip>

#include "Program_options.h"
#include "System.h"

#include "Properties.h"

//----------------------------------------------------------------------

void Md_opt_method::read(vector <string> words,
                      unsigned int & pos,
                      Program_options & options)
{
  pos=0;
  vector <string> vmcsec;
  if(!readsection(words, pos=0, vmcsec, "EMBED")) {
    error("Need EMBED section");
  }
  allocate(vmcsec,options, qmc_avg);

  if(!readvalue(words, pos=0, nstep, "NSTEP")) {
    error("Need NSTEP");
  }

  if(!readvalue(words, pos=0, tstep, "TIMESTEP")) {
    error("Need TIMESTEP in MD method");
  }

  if(readvalue(words, pos=0, damp, "DAMP")) {
    if(damp > 1 || damp < 0) error("DAMP must be in [0,1]");
  }
  else damp=0.0;

}


//----------------------------------------------------------------------

int Md_opt_method::generateVariables(Program_options & options) {


  allocate(options.systemtext[0], sys );
  allocate(options.twftext[0], sys, wfdata);

  sys->generatePseudo(options.pseudotext, pseudo);
  
  return 1;
}

//----------------------------------------------------------------------



int Md_opt_method::showinfo(ostream & os)
{

  if(os)
  {

    sys->showinfo(os);
    os << endl << endl;
    os << "                               Wavefunction  "
       << endl << endl;
    wfdata->showinfo(os);
    pseudo->showinfo(os);
    os << endl;
    
    os << "Embedded QMC" << endl;
    qmc_avg->showinfo(os);


    return 1;
  }
  else
    return 0;
}

//----------------------------------------------------------------------




/*!

*/
void Md_opt_method::run(Program_options & options, ostream & output)
{

  output.precision(10);


  
  string vmcout=options.runid+".embed";
  ofstream vmcoutput;
  if(output) 
	vmcoutput.open(vmcout.c_str());

 
  
  int nparms=wfdata->nparms();
  Array1 <int > parm_num(2*nparms);
  Array1 <doublevar> displace(2*nparms);
  for(int i=0; i< 2*nparms; i+=2) {
    parm_num(i)=i/2;
    parm_num(i+1)=i/2;
    displace(i)=.00025;
    displace(i+1)=-.00025;
  }

  prop.setParmDisplacement(parm_num, displace);

  Array2 <doublevar> parms(nstep+2, nparms);
  Array1 <doublevar> orig_parms(nparms);
  wfdata->getVarParms(orig_parms);
  for(int i=0; i< 2; i++) {
    for(int p=0; p < nparms; p++) 
      parms(i,p)=orig_parms(p);
  }
  Array1 <doublevar> weights(nparms);
  weights=100;

  Properties_final_average curravg;

  for(int step=0; step < nstep; step++) {

    int currt=step+1; //current time

    //cout << "set parms " << endl;
    Array1 <doublevar> curr_parms(nparms);
    for(int p=0; p < nparms; p++) {
      curr_parms(p)=parms(currt,p);
    }
    wfdata->setVarParms(curr_parms);

    //cout << "averaging " << endl;
    qmc_avg->runWithVariables(prop, sys, wfdata, pseudo, vmcoutput);
    //cout << "getting final " << endl;
    prop.getFinal(curravg);
    
    if(output) 
      output << "*****Step " << step << endl;
    Array1 <doublevar> force(nparms, 0.0);
    Array1 <doublevar> force_err(nparms, 0.0);


    //cout << "forces " << endl;
    for(int f=0; f< 2*nparms; f+=2) {
      int p=parm_num(f);
      doublevar fin_diff=(curravg.aux_diff(f,0)-curravg.aux_diff(f+1,0))
        /(2*curravg.aux_size(f));      
      force(p)+= -fin_diff;
          
      force_err(p)=(curravg.aux_differr(f,0)
                          +curravg.aux_differr(f+1,0))
        /(2*curravg.aux_size(f))/(2*curravg.aux_size(f));      
    }
   
    for(int p=0; p < nparms; p++) {
      force_err(p)=sqrt(force_err(p));
      if(fabs(force(p)) < force_err(p) ) {
        force(p)=0;
      }
    }

    //cout << "move " << endl;
    for(int p=0;p < nparms; p++) {

        parms(currt+1, p)=parms(currt,p)
          +(1-damp)*(parms(currt,p)-parms(currt-1, p))
          +tstep*tstep*force(p)/weights(p);
    }

    
    //cout << "write " << endl;
    string wfoutput(options.runid+".wfout");
    string indent="";
    ofstream wfout(wfoutput.c_str());
    wfdata->writeinput(indent, wfout);
    wfout.close();

    if(output) {
      for(int p=0; p < nparms; p++) {
        output << "parm" << p << "   " << parms(currt, p) << " force " 
             << force(p) << " +/- " << force_err(p) << endl;
      }
    }


  }



}

//------------------------------------------------------------------------

