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

#include "Molecular_sample.h"
#include "ulec.h"
#include "Wavefunction.h"




void Molecular_sample::generateStorage(Sample_storage * & store)
{
  int nions=parent->ions.size();
  store=new Sample_storage;
  store->iondist_temp.Resize(5,nions);
  store->pointdist_temp.Resize(5,nelectrons);
  store->pos_temp.Resize(3);
}

void Molecular_sample::saveUpdate(int e, Sample_storage * store) {
  //cout << "saveUpdate" << endl;
  getElectronPos(e, store->pos_temp);
  int nions=parent->ions.size();
  for(int d=0; d< 5; d++)
  {
    for(int i=0; i< nions; i++)
    {
      store->iondist_temp(d,i)=iondist(d,i,e);
    }
    for(int i=0; i<e; i++)
    {
      store->pointdist_temp(d,i)=pointdist(d,i,e);
    }
    for(int j=e; j < nelectrons; j++)
    {
      store->pointdist_temp(d,j)=pointdist(d,e,j);
    }
  }
}

void Molecular_sample::restoreUpdate(int e, Sample_storage * store)
{
  //cout << "restoreUpdate" << endl;
  for(int i=0; i<3; i++)
  {
    elecpos(e,i)=store->pos_temp(i);
  }
  int nions=parent->ions.size();
  for(int d=0; d< 5; d++)
  {
    for(int i=0; i< nions; i++)
    {
      iondist(d,i,e)=store->iondist_temp(d,i);
    }
    for(int i=0; i<e; i++)
    {
      pointdist(d,i,e)=store->pointdist_temp(d,i);
    }
    for(int j=e; j < nelectrons; j++)
    {
      pointdist(d,e,j)=store->pointdist_temp(d,j);
    }
  }
}




void Molecular_sample::updateEIDist()
{
  //cout << "updateEIDist" << endl;
  for(int e=0; e< nelectrons; e++)
  {
    if(ionDistStale(e))
    {
      ionDistStale(e)=0;
      int nions=parent->ions.size();
      for(int j=0; j<nions; j++)
      {
        iondist(1,j,e)=0;
      }

      for(int i=0; i<3; i++)
      {
        for(int j=0; j<nions; j++)
        {
          iondist(i+2,j,e)=elecpos(e,i)-parent->ions.r(i,j);
          iondist(1,j,e)+=iondist(i+2, j, e) * iondist(i+2, j, e);
        }
      }

      for(int j=0; j<nions; j++)
      {
        iondist(0,j, e)=sqrt(iondist(1,j, e));
      }
    }
  }
  //cout << "done" << endl;
}


/*!

 */
void Molecular_sample::updateEEDist()
{
  //cout << "elecelecDistupdate" << endl;
  for(int e=0; e< nelectrons; e++)
  {
    if(elecDistStale(e)==1)
    {
      elecDistStale(e)=0;

      for(int j=0; j<e; j++)
      {
        pointdist(1,j, e)=0;
      }
      for(int k=e+1; k<nelectrons; k++)
      {
        pointdist(1,e,k)=0;
      }

      for(int i=0; i<3; i++)
      {
        for(int j=0; j<e; j++)
        {
          pointdist(i+2,j, e)=elecpos(j,i)-elecpos(e,i);
          //cout << "elecpos  " << j << "  " << i << "  "
          //     << elecpos(j,i) << endl;
          pointdist(1,j, e)+=pointdist(i+2, j, e) * pointdist(i+2, j, e);
        }
        for(int k=e+1; k<nelectrons; k++)
        {
          pointdist(i+2,e, k)=elecpos(e,i)-elecpos(k,i);
          //cout << "elecpos " << k << " " << i << "  "
          //     << elecpos(k,i) << "   " << pointdist(i+2,e,k) << endl;
          pointdist(1,e, k)+=pointdist(i+2, e, k) * pointdist(i+2, e, k);
        }
      }

      for(int j=0; j<e; j++)
      {
        pointdist(0,j, e)=sqrt(pointdist(1,j, e));
        //cout << j << "   " << e << "   " <<  pointdist(0,j,e) << endl;
      }
      for(int k=e+1; k<nelectrons; k++)
      {
        pointdist(0,e, k)=sqrt(pointdist(1,e, k));
        //cout << e << "   " << k << "   " <<  pointdist(0,e,k) << endl;
      }
    }
  }
}


void Molecular_sample::init(System * sys) {
  assert(sys != NULL);
  recast(sys, parent); //assign sys to parent

  int nions=parent->ions.size();
  nelectrons=parent->nelectrons(0)+parent->nelectrons(1);

  elecpos.Resize(nelectrons, 3);
  elecpos=0.0;
  iondist.Resize(5,nions,nelectrons);
  iondist=0.0;
  pointdist.Resize(5, nelectrons, nelectrons);
  pointdist=0.0;
  elecDistStale.Resize(nelectrons);
  ionDistStale.Resize(nelectrons);
  elecDistStale=1;
  ionDistStale=1;
}



void Molecular_sample::randomGuess()
{
  //cout << "randomGuess " << endl;
  const doublevar range=3.0;  //range of the cube around atoms to fill
  Array1 <doublevar> trialPos(3);

  Array1 <doublevar> ioncenter(3); ioncenter=0;
  for(int ion=0; ion < parent->ions.size(); ion++) {
    for(int d=0;d < 3; d++) ioncenter(d)+=parent->ions.r(d,ion)/parent->ions.size();
  }
  
  
  int e=0;
  //Try to put them around the ions
  //cout << "nions " << ions.size() << endl;
  for(int spin=0; spin < 2; spin++)
  {
    for(int ion=0; ion< parent->ions.size(); ion++)
    {
      //cout << "charge " << parent->ions.charge(ion) << endl;
      int laste=e;

      //putting half electrons, rounded up for the first iteration
      //and rounded down for the second one.

      for(;
          e-laste < int(parent->ions.charge(ion)/2)+(1-spin)*int(parent->ions.charge(ion))%2;
          e++)
      {
        if(e< nelectrons) {
          //for(int d=0; d< 3; d++)
          //{
          for(int d=0; d< 3; d++)
          {
            trialPos(d)=parent->ions.r(d,ion)+2*range*(rng.ulec()-.5);
          }
          setElectronPos(e,trialPos);          
            //elecpos(e,d)=parent->ions.r(d,ion)+2*range*(rng.ulec()-.5);
            //cout << "electron position " << e << "  " << d
            //     << "  " << elecpos(e,d) << endl;
          //}
        }
      }
    }
  }


  //Throw the rest randomly
  for(;e<nelectrons; e++)
  {
    for(int d=0; d< 3; d++)
    {
      trialPos(d)=4*range*(rng.ulec()-.5)+ioncenter(d);
      //elecpos(e,d)=4*range*(rng.ulec()-.5);
    }
    setElectronPos(e,trialPos);
  }

  ionDistStale=1;
  elecDistStale=1;

  if(wfObserver)
    wfObserver->notify(all_electrons_move, 0);


}



void Molecular_sample::setElectronPos(const int e,
                                    const Array1 <doublevar> & position)
{
  //cout << "setElectronPos" << endl;
  assert( position.GetDim(0) == 3 );
  Array1 <doublevar> temp(position.GetDim(0));
  temp=position;
  if(parent->use_bounding_box) {
    parent->bounding_box.enforcePbc(temp);
  }

  for(int i=0; i<3; i++)
  {
    elecpos(e,i) = temp(i);
  }
  ionDistStale(e)=1;
  elecDistStale(e)=1;


  if(wfObserver)
    wfObserver->notify(electron_move, e);
}


void Molecular_sample::setElectronPosNoNotify(const int e, 
                                    const Array1 <doublevar> & position) {
  assert( position.GetDim(0) == 3 );
  Array1 <doublevar> temp(position.GetDim(0));
  temp=position;
  if(parent->use_bounding_box) {
    parent->bounding_box.enforcePbc(temp);
  }

  for(int i=0; i<3; i++)
  {
    elecpos(e,i) = temp(i);
  }
  ionDistStale(e)=1;
  elecDistStale(e)=1;
}

void Molecular_sample::getElectronPos(const int e, Array1 <doublevar> & R){
    
    assert( R.GetDim(0) >= 3 );
    
    for(int i=0; i< 3; i++) {
      R(i)=elecpos(e,i);
    }
    
  }


void Molecular_sample::moveIon(const int ion, const Array1 <doublevar> & r)
{
  assert(r.GetDim(0) >=3);
  Array1 <doublevar> temp(r.GetDim(0));
  temp=r;
  if(parent->use_bounding_box) {
    parent->bounding_box.enforcePbc(temp);
  }


  for(int i=0; i< 3; i++) {

    parent->ions.r(i,ion)=temp(i);
  }
  elecDistStale=1;
  ionDistStale=1;

  if(wfObserver)
    wfObserver->notify(ion_move, ion);
}




#include <iomanip>

void Molecular_sample::rawOutput(ostream & os)
{
  os << "Molecular_sample\n";
  os << "numIons  " << iondist.GetDim(1) << endl;
  for(int i=0; i< parent->ions.size(); i++)
  {
    for(int d=0; d< 3; d++)
    {
      os << setw(20) << setprecision(12) << parent->ions.r(d,i);
    }
    os << "\n";
  }
  os << "   numElectrons   " << nelectrons << endl;
  for(int e=0; e< nelectrons; e++) {
    for(int d=0; d< 3; d++) {
      os << setw(20) << setprecision(12) << elecpos(e,d);
    }
    os << "\n";
  }
}


void Molecular_sample::rawInput(istream & is)
{
  string text;
  int nions;
  is >> text;
  if(text != "Molecular_sample")
    error("expected Molecular_sample, got ", text);

  is >> text;
  if(text != "numIons")
    error("expected numIons, got ", text);

  is >> nions;

  Array2 <doublevar> old_ionpos(nions, 3);
  for(int i=0; i< nions; i++) {
    for(int d=0; d< 3; d++) {
      //Reading in the ionic coordinates is confusing and screws up
      //continuous sampling!
      //is >> parent->ions.r(d,i);
      //double dummy;
      is >> old_ionpos(i,d);
    }
  }

  is >> text;
  if(text != "numElectrons")
    error("expected numElectrons, got ", text);

  is >> nelectrons;
  for(int e=0; e< nelectrons; e++) {
    for(int d=0; d< 3; d++) {
      is >> elecpos(e,d);
    }
  }


  //shift the electronic positions if the ions have moved
  updateEIDist();
  for(int e=0; e< nelectrons; e++) {
    Array1 <doublevar> displace(3,0.0);
    doublevar norm=0.0;
    Array1 <doublevar> dist(5);
    for(int at=0; at < nions; at++) {
      getEIDist(e,at, dist);
      doublevar z=getIonCharge(at);
      doublevar func=z/(dist(1)*dist(1));
      norm+=func;      
      //cout << "func " << func << "  dist " << dist(1) << 
      //  "  z " << z << endl;
      for(int d=0; d< 3; d++) 
	displace(d) += func*(parent->ions.r(d,at)-old_ionpos(at,d));
    }
    //cout << "norm = " << norm << endl;
    if(norm>0) {
      for(int d=0; d< 3; d++) 
        elecpos(e,d)+=displace(d)/norm;
      //cout << "displacing electron " << e << " by "
      //     << displace(0)/norm << "   " << displace(1)/norm
      //     << "   " << displace(2)/norm << endl;
    }
  }
     


}

//-------------------------------------------------------------------------
