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


#include "SHO_system.h"
#include "Array.h"
#include "Wavefunction.h"
#include "Sample_point.h"
#include "SHO_sample.h"
#include "qmc_io.h"

void SHO_system::notify(change_type change, int n)
{
  switch(change)
  {
  default:
    cout << "WARNING: SHO system got a signal that it doesn't know: "
    << change << endl;
  }
}


int SHO_system::generateSample(Sample_point *& samptr) {
  assert(samptr==NULL);
  samptr=new SHO_sample;
  samptr->init(this);
  return 1;
}

int SHO_system::showinfo(ostream & os) {
  os << "SHO system " << endl;
  return 1;
}

int SHO_system::read(vector <string> & words,
                           unsigned int & pos){

  pos=0;
  vector <string> omegatxt;
  if(!readsection(words, pos=0, omegatxt, "OMEGA")) 
    omegatxt.push_back("1.0");

  int nd=omegatxt.size();
  omega.Resize(nd);
  for(int d=0; d< nd; d++) {
    omega(d)=atof(omegatxt[d].c_str());
  }
 

  vector <string> nspintxt;
  if(!readsection(words, pos=0, nspintxt, "NSPIN"))
    error("Must have NSPIN in system section");

  if(nspintxt.size()!=2) 
    error("NSPIN must have 2 elements");

  nspin.Resize(2);
  nspin(0)=atoi(nspintxt[0].c_str());
  nspin(1)=atoi(nspintxt[1].c_str());

  
  vector <string> origintxt;
  if(readsection(words, pos=0, origintxt, "ORIGIN")) {
    if(origintxt.size()!=3) error("ORIGIN must have 3 components");
    for(int d=0; d< 3; d++) 
      origin(d)=atof(origintxt[d].c_str());
  }
  
  return 1;
}

//----------------------------------------------------------------------


doublevar SHO_system::calcLoc(Sample_point * sample)
{

  doublevar pot=0;
  int nelectrons=sample->electronSize();

  /*
  sample->updateEEDist();


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
  */
  Array1 <doublevar> pos(3);
  int nd=omega.GetDim(0);
  for(int e=0; e< nelectrons; e++) {
    sample->getElectronPos(e,pos);
    for(int d=0; d< nd; d++) {
      pot+=.5*omega(d)*omega(d)*(pos(d)-origin(d))*(pos(d)-origin(d));
    }
  }

  return pot;
}


//------------------------------------------------------------------------


