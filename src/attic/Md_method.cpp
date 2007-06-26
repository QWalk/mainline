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
//src/Md_method.cpp
//Author:Lucas Wagner



#include "Md_method.h"
#include "qmc_io.h"

#include <iomanip>

#include "Program_options.h"
#include "System.h"

#include "Properties.h"

//----------------------------------------------------------------------

void Md_method::read(vector <string> words,
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

  readvalue(words, pos=0, readcheckfile, "READCHECK");
  readvalue(words, pos=0, writecheckfile, "STORECHECK");

  if(!readvalue(words, pos=0, log_label, "LABEL"))
    log_label="md";


  if(readvalue(words, pos=0, damp, "DAMP")) {
    if(damp > 1 || damp < 0) error("DAMP must be in [0,1]");
  }
  else damp=0.0;

  string resdim;
  pos=0;
  restrict_dimension.Resize(3);
  restrict_dimension=0;
  if(readvalue(words, pos, resdim, "RESTRICT")) {
    if(caseless_eq(resdim, "X")) restrict_dimension(0)=1;
    else if(caseless_eq(resdim, "Y")) restrict_dimension(1)=1;
    else if(caseless_eq(resdim, "Z")) restrict_dimension(2)=1;
    else error("Unknown argument to RESTRICT:", resdim);
  }

  vector <string> weighttxt;
  if(!readsection(words, pos=0, weighttxt, "ATOMIC_WEIGHTS") ) {
    error("Need ATOMIC_WEIGHTS in MD method");
  }
  atomic_weights.Resize(weighttxt.size());
  for(int i=0; i< atomic_weights.GetDim(0); i++) {
    atomic_weights(i)=atof(weighttxt[i].c_str());
  }

}


//----------------------------------------------------------------------

int Md_method::generateVariables(Program_options & options) {


  allocate(options.systemtext[0], sys );
  allocate(options.twftext[0], sys, wfdata);

  sys->generatePseudo(options.pseudotext, pseudo);
  
  return 1;
}

//----------------------------------------------------------------------



int Md_method::showinfo(ostream & os)
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


void recenter(Array2 <doublevar> & pos, Array1 <doublevar> & mass) {
  int natoms=pos.GetDim(0);
  assert(pos.GetDim(0)==mass.GetDim(0));
  assert(pos.GetDim(1)==3);

  Array1 <doublevar> center(3,0.0); //center of mass
  doublevar tot=0.0;
  for(int at=0; at< natoms; at++) {
    for(int d=0; d< 3; d++) {
      center(d)+=mass(at)*pos(at,d);
    }
    tot+=mass(at);
  }

  for(int d=0; d< 3; d++) {
    center(d)/=tot;
  }

  for(int at=0; at< natoms; at++) {
    for(int d=0; d< 3; d++) {
     pos(at,d)-=center(d);
     }
  }

}

/*!

*/
void Md_method::run(Program_options & options, ostream & output)
{

  output.precision(10);


  string logfile=options.runid+".log";
  prop.setLog(logfile, log_label);

  int natoms=sys->nIons();

  Array1 < Array1 <int> > atom_move;
  Array1 < Array2 <doublevar> > displace;

  //Even numbered displacements are in + direction, odd are in 
  // - direction

  if(natoms==2) {
    //For a dimer, assume they are oriented in the z 
    //direction and only move in that direction.
    atom_move.Resize(4);
    displace.Resize(4);
    int count=0;
    for(int at=0; at < natoms; at++) {
      for(int s=0; s< 2; s++) {
        atom_move(count).Resize(1);
        displace(count).Resize(1,3);
        atom_move(count)(0)=at;
        displace(count)=0;
          if(s==0) 
            displace(count)(0,2)=0.00025;
          else
            displace(count)(0,2)=-0.00025;
          count++;
      }
    }
    
  }
  else {
    int ndim=0;
    for(int d=0; d< 3; d++) {
      if(restrict_dimension(d) ==0) ndim++;
    }
    atom_move.Resize(2*ndim*natoms);
    displace.Resize(2*ndim*natoms);
    int count=0;
    for(int at=0; at< natoms; at++) {
      for(int d=0; d< 3; d++) {
        if(!restrict_dimension(d) ) {
        for(int s=0; s< 2; s++) {
          atom_move(count).Resize(1);
          displace(count).Resize(1,3);
          atom_move(count)(0)=at;
          displace(count)=0;
          if(s==0) 
            displace(count)(0,d)=0.00025;
          else
            displace(count)(0,d)=-0.00025;
          count++;
        }
      }
      }
    }
  }
  
  prop.setDisplacement(atom_move, displace);
  Properties_final_average curravg;
  
  string vmcout=options.runid+".embed";
  ofstream vmcoutput;
  if(output) 
	vmcoutput.open(vmcout.c_str());

 
  Array3 <doublevar> ionpos(nstep+2, natoms, 3, 0.0);
  Array1 <doublevar> temppos(3);


  for(int s=0; s< 2; s++) {
    for(int at=0; at< natoms; at++) {
      sys->getIonPos(at, temppos);
      for(int d=0; d< 3; d++) {
        ionpos(s,at,d)=temppos(d);
      }
    }
  }

  if(readcheckfile != "") {
    read_check(ionpos);
  }


  for(int s=0; s< 2; s++) {
    Array2 <doublevar> pos(ionpos(s));
    recenter(pos, atomic_weights);
  }

  for(int step=0; step < nstep; step++) {

    int currt=step+1; //current time

    for(int at=0; at< natoms; at++) {
      for(int d=0; d< 3; d++) {
        temppos(d)=ionpos(currt, at, d);
      }
      sys->setIonPos(at, temppos);
    }

    qmc_avg->runWithVariables(prop, sys, wfdata, pseudo, vmcoutput);
    prop.getFinal(curravg);
    
    if(output) 
      output << "*****Step " << step << endl;
    Array2 <doublevar> force(natoms, 3, 0.0);
    Array2 <doublevar> force_err(natoms, 3, 0.0);

    for(int f=0; f< atom_move.GetDim(0); f+=2 ) {
      for(int m=0; m < atom_move(f).GetDim(0); m++) {
        int at=atom_move(f)(m);
        for(int d=0; d< 3; d++) {
          //Take a finite difference between the two forces
          doublevar prop=fabs(displace(f)(m,d)/curravg.aux_size(f));
          doublevar fin_diff=(curravg.aux_diff(f,0)-curravg.aux_diff(f+1,0))
            /(2*curravg.aux_size(f));
          force(at,d)+= -prop*fin_diff;
          
          force_err(at,d)+=prop*(curravg.aux_differr(f,0)
                                 +curravg.aux_differr(f+1,0))
            /(2*curravg.aux_size(f))/(2*curravg.aux_size(f));
        }
      }
    }


   
    for(int at=0; at< natoms; at++) {
      for(int d=0; d< 3; d++) 
        force_err(at,d)=sqrt(force_err(at,d));
    }


    //Make sure that Newton's laws are actually followed for
    //two atoms; we can do this for more, as well.
    if(natoms==2) {
      doublevar average=(force(0,2)-force(1,2))/2.0;
      force(0,2)=average;
      force(1,2)=-average;
    }
    
    //Verlet algorithm..
    for(int at=0; at< natoms; at++) {
      for(int d=0; d< 3; d++) {
        ionpos(currt+1, at, d)=ionpos(currt, at,d)
          +(1-damp)*(ionpos(currt, at, d)-ionpos(currt-1, at,d))
          +tstep*tstep*force(at,d)/atomic_weights(at);
        //cout << "pos " << ionpos(currt, at,d)  
        //     << "  last " <<  ionpos(currt-1, at,d) 
        //     << "  weight " << atomic_weights(at) << endl;
      }
    }

    Array2 <doublevar> currpos(ionpos(currt+1));
    recenter(currpos, atomic_weights);
    
    doublevar kinen=0;
    int field=output.precision() + 15;
    
    for(int at=0; at < natoms; at++) {
      if(output) 

        output << "position" << at << " " 
             << setw(field)  << ionpos(currt+1, at, 0)
             << setw(field) << ionpos(currt+1, at, 1) 
             << setw(field) << ionpos(currt+1, at, 2) << endl;

      Array1 <doublevar> velocity(3);
      for(int d=0; d< 3; d++) {
        velocity(d)=(ionpos(currt+1, at, d)-ionpos(currt-1, at,d))/(2*tstep);
      }
      
      for(int d=0;d < 3; d++) {
        kinen+=.5*atomic_weights(at)*velocity(d)*velocity(d);
      }
      if(output ) {
        output << "velocity" << at << " "
               << setw(field) << velocity(0) 
               << setw(field) << velocity(1) 
               << setw(field) << velocity(2) << endl;
        
        output << "force" << at << "    "
               << setw(field) << force(at,0) 
               << setw(field) << force(at, 1) 
               << setw(field) << force(at,2) << endl;
        output << "force_err" << at 
               << setw(field) << force_err(at, 0)
               << setw(field) << force_err(at, 1)
               << setw(field) << force_err(at, 2) << endl;
      //output << "force" << at << "z  " << force(at,2) << " +/- "
      //      << force_err(at,2) << endl;
      }

    }
    if(output ) { 
    output << "kinetic_energy " << kinen << endl;
    output << "electronic_energy " << curravg.total_energy(0) << " +/- " 
           << sqrt(curravg.total_energyerr(0)) << endl;

    output << "total_energy " << curravg.total_energy(0)+kinen << " +/- "
           << sqrt(curravg.total_energyerr(0)) <<  endl;

    if(writecheckfile != "") 
      if(output) 
         write_check(ionpos, currt+1);
    }

  }
  if(output) 
    vmcoutput.close();


}

//------------------------------------------------------------------------

void Md_method::read_check(Array3 <doublevar> & last_pos) {
  ifstream is(readcheckfile.c_str());
  if(!is) error("Couldn't open ", readcheckfile);
  is.precision(15);
  string dummy;
  is >> dummy;
  assert(dummy=="natoms");
  int natoms; is >> natoms;
  assert(last_pos.GetDim(1)==natoms);
  assert(last_pos.GetDim(2)==3);
  is >> dummy; //current_pos
  for(int at=0; at< natoms; at++) {
    for(int d=0; d < 3; d++) 
      is >> last_pos(1,at,d);
  }
  is >> dummy; //last_pos
  for(int at=0; at< natoms; at++) {
    for(int d=0; d < 3; d++) 
      is >> last_pos(0,at,d);
  }
  is.close();
}

void Md_method::write_check(Array3 <doublevar> & pos, int step) {
  assert(step > 0);

  ofstream os(writecheckfile.c_str());
  os.precision(15);
  if(!os) error("Couldn't open ", writecheckfile);
  assert(pos.GetDim(2)==3);
  int natoms=pos.GetDim(1);
  os << "natoms " << natoms << endl;
  os << "current_pos" << endl;
  for(int at=0; at< natoms; at++) {
    for(int d=0; d < 3; d++) {
      os << pos(step, at, d) << "  ";
    }
    os << endl;
  }
  os << "last_pos" << endl;
  for(int at=0; at< natoms; at++) {
    for(int d=0; d < 3; d++) {
      os << pos(step-1, at, d) << "  ";
    }
    os << endl;
  }
  os.close();
}

//----------------------------------------------------------------------
