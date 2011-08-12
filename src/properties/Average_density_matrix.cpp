/*
 
 Copyright (C) 2011 Lucas K. Wagner
 
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

#include "Average_density_matrix.h"
#include "ulec.h"

void Average_tbdm_basis::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & avg) { 
  avg.type="tbdm_basis";
  avg.vals.Resize(4*nmo);
  avg.vals=0;
  wf->updateVal(wfdata,sample);

  Array2 <doublevar> movals1(nmo,1),movals2(nmo,1),movals1_old(nmo,1),movals2_old(nmo,1);
  Wf_return wfval_base(wf->nfunc(),2);
  Wf_return wfval_2b(wf->nfunc(),2),wfval_1b(wf->nfunc(),2); //
  wf->getVal(wfdata,0,wfval_base);
  int nelec_1b=sys->nelectrons(0);
  int npairs=sys->nelectrons(0)*sys->nelectrons(1);
  for(int i=0; i< npoints_eval ;i++) { 
    Array1 <doublevar> r1(3),r2(3),oldr1(3),oldr2(3);
    int k=int(rng.ulec()*sys->nelectrons(0));
    int l=int(rng.ulec()*sys->nelectrons(1))+sys->nelectrons(0);
    
    sample->getElectronPos(k,oldr1);
    sample->getElectronPos(l,oldr2);
    for(int d=0; d< 3; d++) {
      r1(d)=10*rng.ulec()-5;
      r2(d)=10*rng.ulec()-5;
    }
    //Calculate the orbital values for r1 and r2
    momat->updateVal(sample,k,0,movals1_old); 
    momat->updateVal(sample,l,0,movals2_old); 
    //
    sample->setElectronPos(k,r1);
    momat->updateVal(sample,k,0,movals1); 
    wf->updateVal(wfdata,sample);
    wf->getVal(wfdata,0,wfval_1b);
    sample->setElectronPos(l,r2);
    momat->updateVal(sample,l,0,movals2); 

    //calculate Psi(r1,r2, etc)
    wf->updateVal(wfdata,sample);
    wf->getVal(wfdata,0,wfval_2b);
    //Accumulate matrix elements 
    doublevar psiratio_2b=exp(wfval_2b.amp(0,0)-wfval_base.amp(0,0))
        *wfval_2b.sign(0)*wfval_base.sign(0);
    doublevar psiratio_1b=exp(wfval_1b.amp(0,0)-wfval_base.amp(0,0))
        *wfval_1b.sign(0)*wfval_base.sign(0);
    for(int orbnum=0; orbnum < nmo; orbnum++) { 
      avg.vals(nmo*2+orbnum)+=npairs*(movals1(orbnum,0)*movals1_old(orbnum,0)
          *movals2(orbnum,0)*movals2_old(orbnum,0) 
          *psiratio_2b)/npoints_eval ;
      avg.vals(nmo*3+orbnum)+=0.5*(movals1(orbnum,0)*movals1(orbnum,0)+movals2(orbnum,0)*movals2(orbnum,0))/npoints_eval;

      avg.vals(orbnum)+=nelec_1b*movals1(orbnum,0)*movals1_old(orbnum,0)*psiratio_1b/npoints_eval;
      avg.vals(nmo+orbnum)+=movals1(orbnum,0)*movals1(orbnum,0)/npoints_eval;

    }
    //Restore the electronic positions
    sample->setElectronPos(k,oldr1);
    sample->setElectronPos(l,oldr2);
  }

}

//----------------------------------------------------------------------

void Average_tbdm_basis::read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words) { 
  unsigned int pos=0;
  vector <string> mosec;
  if(!readsection(words,pos=0, mosec,"ORBITALS")) { 
    error("Need ORBITALS section in TBDM_BASIS");
  }
  allocate(mosec,sys,momat);


  Array1 <Array1 <int> > occupations(1);
  occupations[0].Resize(momat->getNmo());
  for(int i=0; i< momat->getNmo(); i++) { 
    occupations[0][i]=i;
  }
  npoints_eval=100;
  momat->buildLists(occupations);
  nmo=momat->getNmo();
  
}

//----------------------------------------------------------------------

void Average_tbdm_basis::write_init(string & indent, ostream & os) { 
  os << indent << "TBDM_BASIS" << endl;
  os << indent << "NMO " << nmo << endl;
  os << indent << "NPOINTS " << npoints_eval << endl;
  os << indent << "ORBITALS { \n";
  momat->writeinput(indent,os);
  os << indent << "}\n";
}
//----------------------------------------------------------------------
void Average_tbdm_basis::read(vector <string> & words) { 
  unsigned int pos=0;
  //if(!readsection(words,pos=0, mosec,"ORBITALS")) { 
  //  error("Need ORBITALS section in TBDM_BASIS");
  //}
  //allocate(mosec,sys,momat);
  //Array1 <Array1 <int> > occupations(1);
  //occupations[0].Resize(momat->getNmo());
  //for(int i=0; i< momat->getNmo(); i++) { 
  //  occupations[0][i]=i;
  //}

  readvalue(words, pos=0,nmo,"NMO");
  readvalue(words,pos=0,npoints_eval,"NPOINTS");
}
//----------------------------------------------------------------------
void Average_tbdm_basis::write_summary(Average_return &avg,Average_return &err, ostream & os) { 
  for(int orbnum=0; orbnum < nmo; orbnum++) { 
    os << "Average_obdm_basis " << orbnum << setw(17) 
      << avg.vals(orbnum) << " +/- " << setw(16) << err.vals(orbnum) 
      << setw(16) << avg.vals(nmo+orbnum) << " +/- "<< setw(16) << err.vals(nmo+orbnum) <<  setw(17) << avg.vals(orbnum)/avg.vals(nmo+orbnum) << " +/- " << err.vals(orbnum)/avg.vals(nmo+orbnum) <<  endl;
  }
  for(int orbnum=2*nmo; orbnum < 3*nmo; orbnum++) { 
    os << "Average_tbdm_basis " << orbnum-2*nmo << setw(16) 
      << avg.vals(orbnum) << " +/- " << setw(16) << err.vals(orbnum) 
      << setw(16) << avg.vals(nmo+orbnum) << " +/- "<< setw(16) << err.vals(nmo+orbnum) <<  setw(17) << avg.vals(orbnum)/avg.vals(nmo+orbnum)/avg.vals(nmo+orbnum) 
      << " +/- " << setw(16) << err.vals(orbnum)/avg.vals(nmo+orbnum)/avg.vals(nmo+orbnum) << endl;
  }
  
}

//----------------------------------------------------------------------
