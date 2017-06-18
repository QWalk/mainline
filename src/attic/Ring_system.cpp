/*
 
Copyright (C) 2007 Lucas K. Wagner
 with modifications by Pavel Vagner

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


#include "Ring_system.h"
#include "Array.h"
#include "Wavefunction.h"
#include "Sample_point.h"
#include "Ring_sample.h"
#include "qmc_io.h"

void Ring_system::notify(change_type change, int n)
{
  switch(change)
  {
  default:
    cout << "WARNING: Ring system got a signal that it doesn't know: "
         << change << endl;
  }
}


int Ring_system::generateSample(Sample_point *& samptr) {
  assert(samptr==NULL);
  samptr=new Ring_sample;
  samptr->init(this);
  return 1;
}

int Ring_system::showinfo(ostream & os) {
  os << "Ring system " << endl;
  return 1;
}

int Ring_system::read(vector <string> & words,
                           unsigned int & pos){
  nspin.Resize(2);
  vector <string> spintxt;
  if(!readsection(words, pos=0, spintxt, "NSPIN")) {
    error("Need NSPIN in RING system");
  }
  nspin(0)=atoi(spintxt[0].c_str());
  nspin(1)=atoi(spintxt[1].c_str());


  if(!readvalue(words, pos=0, radius, "RADIUS")) {
    error("need RADIUS in RING system");
  }
  length=2*pi*radius;
  //single_write(cout, "length ", length, "\n");

  if(!readvalue(words, pos=0, ee_0, "EE_0"))
    error("Need EE_0 in RING system");

  if(!readvalue(words, pos=0, flux, "FLUX")) 
    error("need FLUX in RING system");


  if(haskeyword(words, pos=0, "NO_EE"))
    ee_interaction=0;
  else ee_interaction=1;


  if(haskeyword(words, pos=0, "EXP_INTERACTION"))
  { 
    exp_interaction=1;
    if(!readvalue(words, pos=0, ee_range, "EE_RANGE"))
      error("Need EE_RANGE in RING system");
  }
  else exp_interaction=0;

  

  //cout << "ee-interaction " << ee_interaction << endl;
  flux=flux/radius;


  vector < vector < string > > barriertxt;
  vector <string> tmp;
  pos=0;
  while(readsection(words, pos, tmp, "BARRIER"))
    barriertxt.push_back(tmp);

  nbarrier=barriertxt.size();
  barrier_strength.Resize(nbarrier);
  barrier_lower.Resize(nbarrier);
  barrier_upper.Resize(nbarrier);

  for(int i=0; i< nbarrier; i++) {
    if(!readvalue(barriertxt[i], pos=0, barrier_strength(i), "STRENGTH"))
      error("Need STRENGTH in BARRIER section");
    if(!readvalue(barriertxt[i], pos=0, barrier_lower(i), "LOWER"))
      error("Need LOWER in BARRIER section");
    if(!readvalue(barriertxt[i], pos=0, barrier_upper(i), "UPPER"))
      error("Need UPPER in BARRIER section");
    if(barrier_lower(i) > barrier_upper(i) ) 
      error("barrier lower must be less than barrier upper");
  }
    


  return 1;
}

//----------------------------------------------------------------------


doublevar Ring_system::calcLoc(Sample_point * sample)
{
  int nelectrons=sample->electronSize();
  sample->updateEEDist();
  doublevar pot=0;

  doublevar elecElec=0;
  Array1 <doublevar> R2(5);
  if(ee_interaction) {
    for(int i=0; i< nelectrons; i++) {
      for(int j=0; j<i; j++) {
        sample->getEEDist(j,i,R2);
        doublevar dist=R2(0);
        if(exp_interaction){
          if(2.0*dist>length) { cout << "wrapping " << endl; dist=length-dist;}
          elecElec+=ee_0 * exp( -dist/ee_range );
        }
        else{
          // doublevar dist=2*radius*sin(pi*R2(0)/length);
          dist=2*radius*sin(pi*dist/length);
          elecElec+= 1/(dist+ee_0);
        }
      }
    }
    pot+=elecElec;
  }

  Array1 <doublevar> pos(3);

  for(int e=0; e< nelectrons; e++) {
    sample->getElectronPos(e,pos);
    for(int b=0; b < nbarrier; b++) {
      if(pos(0) > barrier_lower(b) && pos(0) < barrier_upper(b))
        pot+=barrier_strength(b);
    }
  }

  return pot;
}


//------------------------------------------------------------------------


void Ring_system::calcKinetic(Wavefunction_data * wfdata, 
                              Sample_point * sample,
                              Wavefunction * wf,
                              Array1 <doublevar> & kinetic) {


  wf->updateLap(wfdata, sample);
  int nwf=wf->nfunc();
  Wf_return lap(nwf, 5);
  assert(kinetic.GetDim(0) >= nwf);
  kinetic=0;
  for(int e=0; e< nelectrons(0)+nelectrons(1); e++) {
    for(int w=0; w< nwf; w++) {
      wf->getLap(wfdata, e, lap);
      
      kinetic(w)-=0.5*lap.cvals(w,4).real();
      kinetic(w)+=flux*flux*.5;
      kinetic(w)+=flux*lap.cvals(w,1).imag();

    }
  }    
  
}

//----------------------------------------------------------------------
