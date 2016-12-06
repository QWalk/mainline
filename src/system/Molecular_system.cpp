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


#include "Molecular_system.h"
#include "Array.h"
#include "Wavefunction.h"
#include "Sample_point.h"
#include "Molecular_sample.h"
#include "qmc_io.h"

void Molecular_system::notify(change_type change, int n)
{
  switch(change) {
  default:
    cout << "WARNING: Molecular system got a signal that it doesn't know: "
    << change << endl;
  }
}


int Molecular_system::generateSample(Sample_point *& samptr)
{
  assert(samptr==NULL);
  samptr=new Molecular_sample;
  samptr->init(this);
  return 1;
}

int Molecular_system::showinfo(ostream & os)
{
  ions.showinfo(os);
  return 1;
}

int Molecular_system::read(vector <string> & words,
                           unsigned int & pos)
{
  nspin.Resize(2);
  unsigned int startpos=pos;
  vector <string> spintxt;
  if(!readsection(words, pos, spintxt, "NSPIN")) {
    error("Need NSPIN in molecular system");
  }
  nspin(0)=atoi(spintxt[0].c_str());
  nspin(1)=atoi(spintxt[1].c_str());

  //restrcit initial walkers to a given range
  vector <string> inirangetxt;
  if(readsection(words, pos=0, inirangetxt, "INIRANGE")) {
    cout << "You are using a range of initial walkers other than default." << endl;
    if(inirangetxt.size()==0){
      error("You want to restrict the range of initial walkers, but you didn't give a range");
    }
    inirange=atof(inirangetxt[0].c_str());
  }
  else {
    inirange=3.0;
  }
  
  //use a bounding box if it's given
  vector <string> boxtxt;
  if(readsection(words, pos=0, boxtxt, "BOUNDING_BOX")) {
    const int ndim=3;
    if(boxtxt.size() != ndim*ndim) {
      error("need ", ndim*ndim," elements in BOUNDING_BOX");
    }
    Array2 <doublevar> latVec(ndim, ndim);
    for(int i=0; i< ndim; i++) {
      for(int j=0; j< ndim; j++) {
        latVec(i,j)=atof(boxtxt[i*ndim+j].c_str());
      }
    }
    Array1 <doublevar> origin(ndim);
    origin=0;
    vector <string> origintxt;
    if(readsection(words, pos=0, origintxt, "ORIGIN")) {
      for(int i=0; i < ndim; i++) origin(i)=atof(origintxt[i].c_str());
    }
    bounding_box.init(latVec);
    bounding_box.setOrigin(origin);
    use_bounding_box=1;
  }
  else {
    use_bounding_box=0;
  }
  
  
  
  ions.read(words, pos=0);
  int natoms=ions.size();
  for(int i=0; i< natoms; i++) {
    atomLabels.push_back(ions.getLabel(i));
  }


  electric_field.Resize(3);
  electric_field=0;
  vector <string> electxt;
  if(readsection(words, pos=0, electxt, "ELECTRIC_FIELD")) { 
    if(electxt.size()!=3) error("ELECTRIC_FIELD must have 3 components");
    for(int d=0; d< 3; d++) {
      electric_field(d)=atof(electxt[d].c_str());
      //cout << "efield " << electric_field(d) << endl;
    }
  }


  //---Pseudopotential
  pos=startpos;
  vector < vector <string> > pseudotext;
  vector <string> pseudotexttmp;
  if(readsection(words, pos, pseudotexttmp, "PSEUDO") != 0) {
    error("PSEUDO should go in the global space");
    pseudotext.push_back(pseudotexttmp);
  }

  //pseudo.read(pseudotext, this);
  return 1;
}

//----------------------------------------------------------------------

void Molecular_system::setIonPos(int ion, Array1 <doublevar> & r)
{
  assert(r.GetDim(0) >=3);
  Array1 <doublevar> temp(r.GetDim(0));
  temp=r;
  if(use_bounding_box) {
    bounding_box.enforcePbc(temp);
  }


  for(int i=0; i< 3; i++) {

    ions.r(i,ion)=temp(i);
  }
}



//------------------------------------------------------------------------



void Molecular_system::calcLocWithTestPos(Sample_point * sample,
                                          Array1 <doublevar> & tpos,
                                          Array1 <doublevar> & Vtest) {

  int nions=sample->ionSize();
  int nelectrons=sample->electronSize();
  Vtest.Resize(nelectrons + 1);
  Vtest = 0; 
  Array1 <doublevar> oldpos(3);
  //cout << "Calculating local energy\n";
  sample->getElectronPos(0, oldpos);
  sample->setElectronPosNoNotify(0, tpos);
  
  Array1 <doublevar> R(5);

  sample->updateEIDist();
  sample->updateEEDist();

  for(int i=0; i < nions; i++) {
    sample->getEIDist(0, i, R);
    Vtest(nelectrons)+=-sample->getIonCharge(i)/R(0);
  }
  doublevar dist = 0.0; 
  for(int d=0; d<3; d++)  {
    dist += (oldpos(d) - tpos(d))*(oldpos(d) - tpos(d)); 
  }
  dist = sqrt(dist); 
  Vtest(0) = 1/dist; 
  Array1 <doublevar> R2(5);
  for(int i=1; i< nelectrons; i++) {
    sample->getEEDist(0,i,R2);
    Vtest(i)= 1/R2(0);
  }
  sample->setElectronPosNoNotify(0, oldpos);
  sample->updateEIDist();
  sample->updateEEDist();
  //cout << "elec-elec: " << elecElec << endl;
  //cout << "pot " << pot << endl;

  for(int d=0; d< 3; d++) 
    Vtest(nelectrons) -=electric_field(d)*tpos(d);
}


//------------------------------------------------------------------------

doublevar Molecular_system::calcLoc(Sample_point * sample) {
  int nions=sample->ionSize();
  int nelectrons=sample->electronSize();

  //cout << "Calculating local energy\n";

  Array1 <doublevar> R(5);
  doublevar pot=0;
  
  doublevar elecIon=0;
  sample->updateEIDist();
  sample->updateEEDist();

  for(int e=0; e< nelectrons; e++)
  {
    for(int i=0; i < nions; i++)
    {

      sample->getEIDist(e,i, R);
      elecIon+=sample->getIonCharge(i)/R(0);
    }
  }
  elecIon*=-1;
  pot+=elecIon;

  //cout << "elec-ion: " << elecIon << endl;
  Array1 <doublevar> r1(3), r2(3);
  doublevar IonIon=0;
  for(int i=0; i< nions; i++)
  {
    sample->getIonPos(i,r1);
    //cout << i << "   " << r1(2) << endl;
    for(int j=0; j<i; j++)
    {
      sample->getIonPos(j,r2);
      doublevar r=sqrt( (r1(0)-r2(0))*(r1(0)-r2(0))
                        + (r1(1)-r2(1))*(r1(1)-r2(1))
                        + (r1(2)-r2(2))*(r1(2)-r2(2)));
      IonIon+=sample->getIonCharge(i)*sample->getIonCharge(j)/r;
    }
  }
  pot+=IonIon;
  
  //cout << "Ion-ion: " << IonIon << endl;

  doublevar elecElec=0;
  Array1 <doublevar> R2(5);
  for(int i=0; i< nelectrons; i++)
  {
    for(int j=0; j<i; j++)
    {
      sample->getEEDist(j,i,R2);
      elecElec+= 1/R2(0);
    }
  }
  pot+=elecElec;
  //cout << "elec-elec: " << elecElec << endl;
  //cout << "pot " << pot << endl;

  doublevar fieldPot=0;
  Array1 <doublevar> pos(3);
  for(int e=0; e< nelectrons; e++) {
    sample->getElectronPos(e,pos);
    for(int d=0; d< 3; d++) 
      fieldPot-=electric_field(d)*pos(d);
  }
  for(int i=0; i< nions; i++) {
    sample->getIonPos(i,pos);
    for(int d=0; d< 3; d++) 
      fieldPot+=sample->getIonCharge(i)*electric_field(d)*pos(d);
  }

  pot+=fieldPot;

  return pot;
}




void Molecular_system::calcLocSeparated(Sample_point * sample, Array1<doublevar> & PotLoc)
{
  Array1 <doublevar> onebody; 
  Array2 <doublevar> twobody; 
  separatedLocal(sample, onebody, twobody); 
  int nelectrons = sample->electronSize();
  PotLoc = 0.0; 
  for (int e = 0; e < nelectrons; e++) {
    PotLoc(e) += onebody(e); 
    for (int ep =e+1; ep < nelectrons; ep ++) {
      PotLoc(e) += twobody(e, ep); 
      PotLoc(ep) += twobody(e, ep); 
    }
  }
  /*
  int nions=sample->ionSize();
  int nelectrons=sample->electronSize();
  
  //cout << "Calculating local energy\n";

  Array1 <doublevar> R(5);
  doublevar pot=0;

  doublevar elecIon=0;
  sample->updateEIDist();
  sample->updateEEDist();
  PotLoc = 0.0; 
  for(int e=0; e< nelectrons; e++)
  {

    // Electron-ion interaction
    for(int i=0; i < nions; i++)
    {
      sample->getEIDist(e,i, R);
      PotLoc(e) += -sample->getIonCharge(i)/R(0);
    }
    // Electron-electron interaction
    for(int i=0; i<e; i++) {
      sample->getEEDist(i,e, R); // noted that the first index should be smaller than the second index
      PotLoc(e) += 1/R(0);
      PotLoc(i) += 1/R(0);
    }

    // Local external field
    Array1 <doublevar> pos(3);
    sample->getElectronPos(e,pos);
    for(int d=0; d< 3; d++) 
      PotLoc(e) += -electric_field(d)*pos(d);
  }
  */
}

//----------------------------------------------------------------------
//Keeping this separate from calcLoc for efficiency reasons..
//It's most likely faster to add to a double than fill matrices..
//We also don't need to calculate the ion-ion interaction for this
void Molecular_system::separatedLocal(Sample_point * sample,Array1 <doublevar> & onebody,
    Array2 <doublevar> & twobody)
{
  int nions=sample->ionSize();
  int nelectrons=sample->electronSize();
  onebody.Resize(nelectrons);
  twobody.Resize(nelectrons,nelectrons);
  onebody=0.0; twobody=0.0;
  Array1 <doublevar> R(5);
  sample->updateEIDist();
  sample->updateEEDist();

  for(int e=0; e< nelectrons; e++) {
    for(int i=0; i < nions; i++) {
      sample->getEIDist(e,i, R);
      onebody(e)-=sample->getIonCharge(i)/R(0);
    }
  }
  Array1 <doublevar> R2(5);
  for(int i=0; i< nelectrons; i++) {
    for(int j=i+1; j<nelectrons; j++) {
      sample->getEEDist(i,j,R2);
      twobody(i,j)+= 1/R2(0);
    }
  }

  Array1 <doublevar> pos(3);
  for(int e=0; e< nelectrons; e++) {
    sample->getElectronPos(e,pos);
    for(int d=0; d< 3; d++) 
      onebody(e)-=electric_field(d)*pos(d);
  }
}

//----------------------------------------------------------------------

void Molecular_system::locDerivative(int ion, Sample_point * sample,
                                     Force_fitter & fitter,
                                     Array1 <doublevar> & der) {

  sample->updateEIDist();
  der.Resize(3);
  der=0;
  Array1 <doublevar> R(5);
  int ne=sample->electronSize();
  Array1 <doublevar> tmpder(3), tmpder_fit(3);
  for(int e=0; e< ne; e++) {
    sample->getEIDist(e,ion, R);
    for(int d=0; d< 3; d++) 
      //der(d)-=sample->getIonCharge(ion)*R(d+2)/(R(1)*R(0));
      tmpder(d)=-sample->getIonCharge(ion)*R(d+2)/(R(1)*R(0));
    fitter.fit_force(R,tmpder, tmpder_fit);
    for(int d=0; d< 3; d++)
      der(d)+=tmpder_fit(d);
  }

  Array1 <doublevar> vec(3);
  Array1 <doublevar> r_ion(3);
  sample->getIonPos(ion, r_ion);

  int nIon=sample->ionSize();

  for(int j=0; j< nIon; j++) {
    if(j!=ion) {
    sample->getIonPos(j,R);
    doublevar r2=0;
    for(int d=0; d <3; d++) 
      vec(d)=r_ion(d)-R(d);

    for(int d=0; d< 3; d++)
      r2+=vec(d)*vec(d);
    doublevar r=sqrt(r2);

    for(int d=0; d< 3; d++) 
      der(d)-=sample->getIonCharge(ion)*sample->getIonCharge(j)*vec(d)/(r2*r);
    }
  }

}

//------------------------------------------------------------------------


int Molecular_system::getBounds(Array2 <doublevar> & latticevec,Array1<doublevar> & origin) { 
  cout << "getBounds()" << endl;
  latticevec.Resize(3,3);
  latticevec=0.0;
  Array1 <doublevar> minmax(6);
  for(int d=0; d< 3; d++) { 
    minmax(2*d)=minmax(2*d+1)=0.0;
  }
  
  int nions=nIons();
  Array1 <doublevar> ionpos(3);
  for(int i=0; i< nions; i++) {
    getIonPos(i,ionpos);
    for(int d=0; d< 3; d++) {
      if(ionpos(d) < minmax(2*d)) minmax(2*d) = ionpos(d);
      if(ionpos(d) > minmax(2*d+1)) minmax(2*d+1)=ionpos(d);
    }
  }

  for(int d=0; d< 3; d++) {
    origin(d)=minmax(2*d)-8.0;
    latticevec(d,d)=minmax(2*d+1)-origin(d)+8.0;
  }
  return 1;
}
