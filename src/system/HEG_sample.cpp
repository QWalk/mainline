/*
 
Copyright (C) 2007 Lucas K. Wagner, Jindrich Kolorenc

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



#include "HEG_sample.h"
#include "ulec.h"
#include "Wavefunction.h"
#include "HEG_system.h"

//----------------------------------------------------------------------

void HEG_sample::generateStorage(Sample_storage * & store)
{
  store=new Sample_storage;
  store->pointdist_temp.Resize(nelectrons,5);
  store->pos_temp.Resize(3);
}

void HEG_sample::saveUpdate(int e, Sample_storage * store)
{
  getElectronPos(e, store->pos_temp);
  for(int d=0; d< 5; d++) {
    for(int i=0; i<e; i++)
    {
      store->pointdist_temp(i,d)=pointdist(i,e,d);
    }
    for(int j=e; j < nelectrons; j++)
    {
      store->pointdist_temp(j,d)=pointdist(e,j,d);
    }
  }
}

void HEG_sample::restoreUpdate(int e, Sample_storage * store){
  for(int i=0; i<3; i++) {
    elecpos(e,i)=store->pos_temp(i);
  }
  for(int d=0; d< 5; d++)  {

    for(int i=0; i<e; i++)
      pointdist(i,e,d)=store->pointdist_temp(i,d);

    for(int j=e; j < nelectrons; j++)
      pointdist(e,j,d)=store->pointdist_temp(j,d);
  }
}



//----------------------------------------------------------------------

/*!

 */
void HEG_sample::updateEEDist() {

#define ORTHOGONAL_CELL 1

#ifdef ORTHOGONAL_CELL

  // Fast evaluation of electron electron distances. Unfortunately,
  // it REQUIRES simple ORTHOGONAL simulation cell and thus cannot
  // be used in general crystalline environment (for homogeneous
  // gas such a loss of generality does not matter, does it?). It
  // is tempting to replace 1/latVec with recipLatVec and do proper
  // matrix multiplications, but as noted above, it is flawed as
  // a general concept.
  Array1 <doublevar> dist_lc(3);
  for(int e=0; e< nelectrons; e++) {
    if(elecDistStale(e)==1) {
      for (int d=0; d<3; d++)
	elecpos_lc(e,d)=elecpos(e,d)/parent->latVec(d,d);
    }
  }
  
  for(int e=0; e< nelectrons; e++) {
    if(elecDistStale(e)==1) {

      elecDistStale(e)=0;

      for(int j=0; j< e; j++) {
        for(int d=0; d< 3; d++) {
          dist_lc(d)=elecpos_lc(j,d)-elecpos_lc(e,d);
	  dist_lc(d)-=floor(dist_lc(d)+0.5);
	}
	pointdist(j,e,1)=0.0;
	for (int d=0; d<3; d++) {
	  pointdist(j,e,d+2)=parent->latVec(d,d)*dist_lc(d);
	  pointdist(j,e,1)+=pointdist(j,e,d+2)*pointdist(j,e,d+2);
	}
	pointdist(j,e,0)=sqrt(pointdist(j,e,1));
      }
      
      for(int k=e+1; k< nelectrons; k++) {
	for(int d=0; d< 3; d++) {
	  dist_lc(d)=elecpos_lc(e,d)-elecpos_lc(k,d);
	  dist_lc(d)-=floor(dist_lc(d)+0.5);
	}
	pointdist(e,k,1)=0.0;
	for (int d=0; d<3; d++) {
	  pointdist(e,k,d+2)=parent->latVec(d,d)*dist_lc(d);
	  pointdist(e,k,1)+=pointdist(e,k,d+2)*pointdist(e,k,d+2);
	}
	pointdist(e,k,0)=sqrt(pointdist(e,k,1));
      }

    } // if(elecDistStale(e)==1) {
  }  // loop over electrons

#else
  // Original Lucas version from Periodic_sample
  for(int e=0; e< nelectrons; e++) {
    if(elecDistStale(e)==1) {
      
      elecDistStale(e)=0;
      for(int j=0; j<e; j++) {
        pointdist(j, e,1)=0.0;
      }
      for(int k=e+1; k<nelectrons; k++) {
        pointdist(e,k,1)=0.0;
      }
      
      //Find the closest image
      Array1 <doublevar> tmpvec(3);
      doublevar tmpdis;
      doublevar height2=parent->smallestheight*parent->smallestheight*.25;
      
      for(int j=0; j< e; j++) {
        //first try within the primitive cell
        for(int d=0; d< 3; d++) 
          pointdist(j,e,d+2)=elecpos(j,d)-elecpos(e,d);
        for(int d=0; d< 3; d++)
          pointdist(j,e,1)+=pointdist(j,e,d+2)*pointdist(j,e,d+2);
        
        //then check outside if we haven't found the minimum image
        if(pointdist(j,e,1) > height2) {
          for(int a=0; a< 26; a++) {
            for(int d=0; d< 3; d++) 
              tmpvec(d)=elecpos(j,d)-elecpos(e,d)+tmplat(a,d);
	    
            tmpdis=0;
            for(int d=0; d<3; d++) tmpdis+=tmpvec(d)*tmpvec(d);
            
            if(tmpdis < pointdist(j,e,1)) {
              pointdist(j,e,1)=tmpdis;
              for(int d=0; d< 3; d++) {
                pointdist(j, e, d+2)=tmpvec(d);
              }
            }
            if(tmpdis < height2)
              break;
          }
        }
      }
      
      for(int k=e+1; k< nelectrons; k++) {
        //first try within the primitive cell
        for(int d=0; d< 3; d++) 
          pointdist(e,k,d+2)=elecpos(e,d)-elecpos(k,d);
        for(int d=0; d< 3; d++)
          pointdist(e,k,1)+=pointdist(e,k,d+2)*pointdist(e,k,d+2);
        
        //then check outside if we haven't found the minimum image
        if(pointdist(e,k,1) > height2) {
          for(int a=0; a< 26; a++) {
            for(int d=0; d< 3; d++) 
              tmpvec(d)=elecpos(e,d)-elecpos(k,d)+tmplat(a,d);
	    
            tmpdis=0;
            for(int d=0; d<3; d++) tmpdis+=tmpvec(d)*tmpvec(d);
            
            if(tmpdis < pointdist(e,k,1)) {
              pointdist(e,k,1)=tmpdis;
              for(int d=0; d< 3; d++) {
                pointdist(e, k, d+2)=tmpvec(d);
              }
            }
            if(tmpdis < height2)
              break;
          }
        }
      }     
      
      //---------------------------
      for(int j=0; j<e; j++) {
        pointdist(j, e,0)=sqrt(pointdist(j, e,1));
	//cout << j << " " << e << " " << pointdist(j,e,2) << endl;
      }
      for(int k=e+1; k<nelectrons; k++) {
        pointdist(e, k,0)=sqrt(pointdist(e, k,1));
	//cout << e << " " << k << " " << pointdist(e,k,2) << endl;
      }
      
    } // if(elecDistStale(e)==1) {
  }  // loop over electrons

#endif

}


//----------------------------------------------------------------------

void HEG_sample::init(System * sys) {
  assert(sys != NULL);
  recast(sys, parent); //assign sys to parent

  nelectrons=parent->nelectrons(0)+parent->nelectrons(1);

  elecpos.Resize(nelectrons, 3);
  pointdist.Resize(nelectrons, nelectrons,5);
  elecpos_lc.Resize(nelectrons, 3);

  elecDistStale.Resize(nelectrons);
  elecDistStale=1;

  tmplat.Resize(26,3);
  int counter=0;
  for(int aa=-1; aa <= 1; aa++) {
    for(int bb=-1; bb <= 1; bb++) {
      for(int cc=-1; cc <= 1; cc++) {    
	if(aa!=0 || bb !=0 || cc!=0) {
	  for(int d=0; d< 3; d++) {
	    tmplat(counter,d)=aa*parent->latVec(0,d)
	      +bb*parent->latVec(1,d)
	      +cc*parent->latVec(2,d);
	  }
	  counter++;
	}
      }
    }
  }
  
  // update_overall_sign is false for complex-valued wavefunctions,
  // i.e., for non-integer k-points
  update_overall_sign=true;
  for (int d=0; d<3; d++) {
    doublevar intpart;
    modf(parent->kpt(d),&intpart);
    if ( parent->kpt(d) != intpart ) update_overall_sign=false;
  }

}


//----------------------------------------------------------------------

void HEG_sample::randomGuess()
{
  Array1 <doublevar> trialPos(3);

  // homogeneous distribution of electrons in the simulation cell
  for(int e=0; e<nelectrons; e++)
  {
    doublevar rn0=rng.ulec();
    doublevar rn1=rng.ulec();
    doublevar rn2=rng.ulec();
    for(int d=0; d< 3; d++)
    {
      trialPos(d)=rn0*parent->latVec(0,d) +
	rn1*parent->latVec(1,d) + rn2*parent->latVec(2,d);
    }
    setElectronPos(e, trialPos);
  }

  elecDistStale=1;

  if(wfObserver)
  {
    wfObserver->notify(all_electrons_move, 0);
  }

}


//----------------------------------------------------------------------

void HEG_sample::setElectronPos(const int e,
                                   const Array1 <doublevar> & position)
{

  assert( position.GetDim(0) == 3 );
  Array1 <doublevar> temp(position.GetDim(0));
  temp=position;
  Array1 <int> nshift;
  parent->enforcePbc(temp, nshift);

  for(int i=0; i<3; i++)
  {
    //points.r(i,e) = temp(i);
    elecpos(e,i)=temp(i);
  }
  elecDistStale(e)=1;

  if(wfObserver)
  {
    wfObserver->notify(electron_move, e);
  }

}


void HEG_sample::translateElectron(const int e, const Array1 <doublevar> & trans) {
  Array1 <doublevar> temp(trans.GetDim(0));
  
  for(int d=0; d< 3; d++) 
    temp(d)=elecpos(e,d)+trans(d);
  Array1<int> nshift;
  parent->enforcePbc(temp, nshift);
  doublevar kdotr=0;
  for(int d=0; d< 3; d++) 
    kdotr+=parent->kpt(d)*nshift(d);
  if ( update_overall_sign ) overall_sign*=cos(pi*kdotr);
  // the sign of the phase here is critical for correct evaluation of
  // density matrices
  overall_phase-=pi*kdotr;
  // I suppose the value of overall_phase averages around zero, no
  // "modulo" operation should be needed
  //overall_phase=overall_phase - 2*pi * (int) (overall_phase/pi/2);
  
  //if(fabs(kdotr) > 1e-8) 
  //  cout << "shift " << nshift(0) << "   " 
  //      << nshift(1) << "   " << nshift(2) << endl;
  
  for(int i=0; i<3; i++)
    elecpos(e,i)=temp(i);
  
  elecDistStale(e)=1;

  if(wfObserver)
    wfObserver->notify(electron_move, e);
       
}

//----------------------------------------------------------------------

#include <iomanip>

void HEG_sample::rawOutput(ostream & os)
{
  os << "HEG_sample\n";
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


void HEG_sample::rawInput(istream & is)
{
  string text;
  is >> text;
  if(text != "HEG_sample")
    error("expected HEG_sample, got ", text);

  is >> text;
  if(text != "numElectrons")
    error("expected numElectrons, got ", text);
  //cout << "electrons " << endl;
  is >> nelectrons;
  for(int e=0; e< nelectrons; e++)
  {
    for(int d=0; d< 3; d++)
    {
      is >> elecpos(e,d);

    }
  }

  elecDistStale=1;
  
}

//-------------------------------------------------------------------------
