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


#include "ulec.h"
#include "Ring_sample.h"

#include "Wavefunction.h"




void Ring_sample::generateStorage(Sample_storage * & store) {
  //int nions=parent->ions.size();
  store=new Sample_storage;
  //store->iondist_temp.Resize(5,nions);
  store->pointdist_temp.Resize(5,nelectrons);
  store->pos_temp.Resize(3);
}

void Ring_sample::saveUpdate(int e, Sample_storage * store) {
  //cout << "saveUpdate" << endl;
  getElectronPos(e, store->pos_temp);
  //int nions=parent->ions.size();
  for(int d=0; d< 5; d++)
  {
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

void Ring_sample::restoreUpdate(int e, Sample_storage * store) {
  //cout << "restoreUpdate" << endl;
  for(int i=0; i<3; i++)
  {
    elecpos(e,i)=store->pos_temp(i);
  }
  //int nions=parent->ions.size();
  for(int d=0; d< 5; d++)
  {
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




void Ring_sample::updateEIDist()
{
  //cout << "updateEIDist" << endl;
  
  /*
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
  */
  //cout << "done" << endl;
}


/*!

 */
void Ring_sample::updateEEDist()
{
  //cout << "elecelecDistupdate" << endl;
  for(int e=0; e< nelectrons; e++) {
    if(elecDistStale(e)==1) {
      elecDistStale(e)=0;

      for(int j=0; j<e; j++) {
        pointdist(2,j, e)=elecpos(j,0)-elecpos(e,0);
        if(pointdist(2,j,e) > parent->length*.5) pointdist(2,j,e)-=parent->length;
        else if(pointdist(2,j,e) < -parent->length*.5) pointdist(2,j,e)+=parent->length;
        
        pointdist(0,j,e)=fabs(pointdist(2,j,e));
        
        if(pointdist(0,j,e) > parent->length*.5) cout << "distance too big " 
                    << pointdist(0,j,e) << endl;
        pointdist(1,j,e)=pointdist(0,j,e)*pointdist(0,j,e);        
      }
      for(int k=e+1; k<nelectrons; k++) {
        pointdist(2,e,k)=elecpos(e,0)-elecpos(k,0);
        if(pointdist(2,e,k) > parent->length*.5) pointdist(2,e,k)-=parent->length;
        else if(pointdist(2,e,k) < -parent->length*.5) pointdist(2,e,k)+=parent->length;
        
        pointdist(0,e,k)=fabs(pointdist(2,e,k));
        pointdist(1,e,k)=pointdist(0,e,k)*pointdist(0,e,k);

      }

    }
  }
}


void Ring_sample::init(System * sys) {
  assert(sys != NULL);
  recast(sys, parent); //assign sys to parent

  //int nions=parent->ions.size();
  nelectrons=parent->nelectrons(0)+parent->nelectrons(1);

  elecpos.Resize(nelectrons, 3);
  elecpos=0.0;
  //iondist.Resize(5,nions,nelectrons);
  //iondist=0.0;
  pointdist.Resize(5, nelectrons, nelectrons);
  pointdist=0.0;
  elecDistStale.Resize(nelectrons);
  //ionDistStale.Resize(nelectrons);
  elecDistStale=1;
  //ionDistStale=1;
}



void Ring_sample::randomGuess()
{
  //cout << "randomGuess " << endl;

  for(int e=0; e< nelectrons; e++) {
    elecpos(e,0)=rng.ulec()*parent->length;
    elecpos(e,1)=elecpos(e,2)=0;
  }

  //ionDistStale=1;
  elecDistStale=1;

  if(wfObserver) {
    wfObserver->notify(all_electrons_move, 0);
  }

}



void Ring_sample::setElectronPos(const int e,
                                    const Array1 <doublevar> & position)
{
  //cout << "setElectronPos" << endl;
  assert( position.GetDim(0) == 3 );

  //cout << "position " << position(0) << "  " 
  //     << position(1) << "   " << position(2) << endl;
  assert(fabs(position(1)) <1e-8);
  assert(fabs(position(2)) <1e-8);
  
  doublevar x=position(0);
  int count=0;
  while(x>parent->length) {
    x-=parent->length;
    if(count++ > 5000) {
      error("position way off ring ", x);
    }
  }
  while(x < 0) {
    x+=parent->length;
    if(count++ > 5000) {
      error("position way off ring ", x);
    }
  }
  elecpos(e,0)=x;

  //ionDistStale(e)=1;
  elecDistStale(e)=1;


  if(wfObserver)
  {
    wfObserver->notify(electron_move, e);
  }

}


void Ring_sample::moveIon(const int ion, const Array1 <doublevar> & r)
{
  error("Don't support moveIon in RING");
}

//----------------------------------------------------------------------


#include <iomanip>

void Ring_sample::rawOutput(ostream & os)
{
  os << "Ring_sample\n";
  //os << "numIons  " << iondist.GetDim(1) << endl;
  //for(int i=0; i< parent->ions.size(); i++)
  //{
  //  for(int d=0; d< 3; d++)
  //  {
  //    os << setw(20) << setprecision(12) << parent->ions.r(d,i);
  //  }
  //  os << "\n";
  //}
  os << "   numElectrons   " << nelectrons << endl;
  for(int e=0; e< nelectrons; e++)
  {
    for(int d=0; d< 3; d++)
    {
      os << setw(20) << setprecision(12) << elecpos(e,d);
    }
    os << "\n";
  }
}

//----------------------------------------------------------------------


void Ring_sample::rawInput(istream & is)
{
  string text;
  //int nions;
  is >> text;
  if(text != "Ring_sample")
    error("expected Ring_sample, got ", text);

  //is >> text;
  //if(text != "numIons")
  //  error("expected numIons, got ", text);

  //is >> nions;
  //for(int i=0; i< nions; i++)
  //{
  //  for(int d=0; d< 3; d++)
  //  {
      //Reading in the ionic coordinates is confusing and screws up
      //continuous sampling!
      //is >> parent->ions.r(d,i);
  //    double dummy;
  //    is >> dummy;
      //cout << ions.r(d,i) << "   ";
  //  }
    //cout << endl;
  //}

  is >> text;
  if(text != "numElectrons")
    error("expected numElectrons, got ", text);

  is >> nelectrons;
  for(int e=0; e< nelectrons; e++)
  {
    for(int d=0; d< 3; d++)
    {
      is >> elecpos(e,d);

    }
  }


}

//-------------------------------------------------------------------------
