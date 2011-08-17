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
  avg.vals.Resize(nmo+2*nmo*nmo);
  avg.vals=0;
  wf->updateVal(wfdata,sample);

  Array2 <doublevar> movals1(nmo,1),movals2(nmo,1),movals1_old(nmo,1),movals2_old(nmo,1);
  Wf_return wfval_base(wf->nfunc(),2);
  Wf_return wfval_2b(wf->nfunc(),2),wfval_1b(wf->nfunc(),2); //
  wf->getVal(wfdata,0,wfval_base);
  int nelec_1b=sys->nelectrons(0);
  int npairs=sys->nelectrons(0)*sys->nelectrons(1);
  doublevar test_avg_moval1=0, test_avg_moval2=0, test_avg_wfratio=0;
  for(int i=0; i< npoints_eval ;i++) { 
    Array1 <doublevar> r1(3),r2(3),oldr1(3),oldr2(3);
    int k=int(rng.ulec()*sys->nelectrons(0));
    int l=int(rng.ulec()*sys->nelectrons(1))+sys->nelectrons(0);
    //Calculate the orbital values for r1 and r2
    momat->updateVal(sample,k,0,movals1_old); 
    momat->updateVal(sample,l,0,movals2_old); 
    
    sample->getElectronPos(k,oldr1);
    sample->getElectronPos(l,oldr2);
    for(int d=0; d< 3; d++) {
      r1(d)=10*rng.ulec()-5;
      r2(d)=10*rng.ulec()-5;
    }
    sample->setElectronPos(k,r1);
    gen_sample(nstep_sample,.2,k,movals1, sample);
    wf->updateVal(wfdata,sample);

    wf->getVal(wfdata,0,wfval_1b);
    sample->setElectronPos(l,r2); 
    gen_sample(nstep_sample,.2,l,movals2,sample);
    wf->updateVal(wfdata,sample);
    wf->getVal(wfdata,0,wfval_2b);

    //Accumulate matrix elements 
    doublevar psiratio_2b=exp(wfval_2b.amp(0,0)-wfval_base.amp(0,0))
        *wfval_2b.sign(0)*wfval_base.sign(0);
    doublevar psiratio_1b=exp(wfval_1b.amp(0,0)-wfval_base.amp(0,0))
        *wfval_1b.sign(0)*wfval_base.sign(0);

    doublevar dist1=0,dist2=0;
    for(int orbnum=0; orbnum < nmo; orbnum++) {
      dist1+=movals1(orbnum,0)*movals1(orbnum,0);
      dist2+=movals2(orbnum,0)*movals2(orbnum,0);
    }

    //orbital normalization
    for(int orbnum=0; orbnum < nmo; orbnum++) { 
      avg.vals(orbnum)+=0.5*(movals1(orbnum,0)*movals1(orbnum,0)/dist1
            +movals2(orbnum,0)*movals2(orbnum,0)/dist2 )/npoints_eval;
    }
    for(int orbnum=0; orbnum < nmo; orbnum++) { 
      for(int orbnum2=0; orbnum2 < nmo; orbnum2++) { 
        avg.vals(nmo+orbnum*nmo+orbnum2)+=
            npairs*(movals1(orbnum,0)*movals1_old(orbnum,0)
            *movals2(orbnum2,0)*movals2_old(orbnum2,0) 
            *psiratio_2b)/npoints_eval/dist1/dist2 ;

        avg.vals(nmo+nmo*nmo+orbnum*nmo+orbnum2)+=
          nelec_1b*movals1(orbnum,0)*movals1_old(orbnum2,0)
          *psiratio_1b/npoints_eval/dist1;
        //if(orbnum==orbnum2 && orbnum==5) { 
        //  cout << movals1(orbnum,0) << " " << movals1_old(orbnum2,0) 
        //    << " " << psiratio_1b << endl;
//          test_avg_moval1+=movals1(orbnum,0)/npoints_eval;
//          test_avg_moval2+=movals1_old(orbnum,0)/npoints_eval;
//          test_avg_wfratio+=psiratio_1b/npoints_eval;

        //}
      }

    }
    //Restore the electronic positions
    sample->setElectronPos(k,oldr1);
    sample->setElectronPos(l,oldr2);
  }
//  cout << "avg 1bdm " << avg.vals(nmo+nmo*nmo+5*nmo+5) 
//    << " moval1 " << test_avg_moval1
//    << " moval2 " << test_avg_moval2 
 //   << " wfratio " << test_avg_wfratio << endl;

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
  momat->buildLists(occupations);
  nmo=momat->getNmo();

  if(!readvalue(words, pos=0,nstep_sample,"NSTEP_SAMPLE"))
    nstep_sample=10;
  if(!readvalue(words,pos=0,npoints_eval,"NPOINTS"))
    npoints_eval=100;
  
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

  /*
  for(int orbnum=0; orbnum < nmo; orbnum++) { 
    os << "Average_obdm_basis " << orbnum << setw(17) 
      << avg.vals(orbnum) << " +/- " << setw(16) << err.vals(orbnum) 
      << setw(16) << avg.vals(nmo+orbnum) << " +/- "<< setw(16) << err.vals(nmo+orbnum) <<  setw(17) << avg.vals(orbnum)/avg.vals(nmo+orbnum) << " +/- " << err.vals(orbnum)/avg.vals(nmo+orbnum) <<  endl;
  }
  for(int orbnum=2*nmo; orbnum < 3*nmo; orbnum++) { 
    os << "Average_tbdm_basis " << orbnum-2*nmo << setw(17) 
      << avg.vals(orbnum) << " +/- " << setw(16) << err.vals(orbnum) 
      << setw(16) << avg.vals(nmo+orbnum) << " +/- "<< setw(16) << err.vals(nmo+orbnum) <<  setw(17) << avg.vals(orbnum)/avg.vals(nmo+orbnum)/avg.vals(nmo+orbnum) 
      << " +/- " << setw(16) << err.vals(orbnum)/avg.vals(nmo+orbnum)/avg.vals(nmo+orbnum) << endl;
  }
  */
  for(int orbnum=0; orbnum < nmo; orbnum++) { 
    for(int orbnum2=0; orbnum2 < nmo; orbnum2++)  { 
      doublevar norm=sqrt(avg.vals(orbnum)*avg.vals(orbnum2));
      doublevar tobdm=avg.vals(nmo+nmo*nmo+orbnum*nmo+orbnum2)/norm;
      doublevar tobdm_err= err.vals(nmo+nmo*nmo+orbnum*nmo+orbnum2)/norm;
      doublevar ttbdm= avg.vals(nmo+orbnum*nmo+orbnum2)/(norm*norm);
      doublevar ttbdm_err=err.vals(nmo+orbnum*nmo+orbnum2)/(norm*norm);
      if(fabs(tobdm)> 3*tobdm_err or fabs(ttbdm) > 3*ttbdm_err) { 
        os << setw(10) << orbnum << setw(10) << orbnum2 
          << setw(17) << tobdm << setw(17) << tobdm_err
          << setw(17) << ttbdm << setw(17) << ttbdm_err << endl;
      }
    }
  }


  Array2 <doublevar> obdm(nmo,nmo);
  Array2 <doublevar> tbdm(nmo,nmo);
  
  for(int orbnum=0; orbnum < nmo; orbnum++) {
    for(int orbnum2=0; orbnum2<nmo; orbnum2++) {
      doublevar norm=sqrt(avg.vals(orbnum)*avg.vals(orbnum2));
      obdm(orbnum,orbnum2)=avg.vals(nmo+nmo*nmo+orbnum*nmo+orbnum2)/norm;
      tbdm(orbnum,orbnum2)= avg.vals(nmo+orbnum*nmo+orbnum2)/(norm*norm);
    }
  }

  doublevar trace=0;
  for(int i=0; i< nmo; i++) trace+=obdm(i,i);
  os << "trace of the obdm " << trace << endl;

  trace=0;
  for(int i=0; i< nmo; i++) trace+=tbdm(i,i);
  os << "trace of the tbdm " << trace << endl;

  //Symmetrize the matrix
  for(int i=0; i< nmo; i++) {
    for(int j=i+1; j< nmo; j++) { 
      obdm(i,j)=obdm(j,i)=0.5*(obdm(i,j)+obdm(j,i));
    }
  }
  Array1 <doublevar> evals(nmo);
  Array2 <doublevar> evecs(nmo,nmo);
  EigenSystemSolverRealSymmetricMatrix(obdm,evals,evecs);

  os << "evals ";
  for(int i=0; i< nmo; i++) os << evals(i) <<" " ;
  os << endl;


  os << "Evecs " << endl;
  for(int i=0; i< nmo; i++) { 
    for(int j=0; j< nmo; j++) { 
      os << setw(17) << evecs(i,j);
    }
    os << endl;
  }


}

//----------------------------------------------------------------------
//
//
//Note this needs to be changed for non-zero k-points!
void Average_tbdm_basis::gen_sample(int nstep, doublevar  tstep, 
    int e, Array2 <doublevar> & movals, Sample_point * sample) { 
  int nmo=momat->getNmo();
  int ndim=3;
  Array1 <doublevar> r(ndim),rold(ndim);
  Array2 <doublevar> movals_old(nmo,1);
  movals.Resize(nmo,1);

  sample->getElectronPos(e,rold);
  momat->updateVal(sample,e,0,movals_old);

  for(int step=0; step < nstep; step++) { 
    for(int d=0; d< ndim; d++) { 
      r(d)=rold(d)+sqrt(tstep)*rng.gasdev();
    }
    sample->setElectronPos(e,r);
    momat->updateVal(sample,e,0,movals); 

    doublevar sum_old=0,sum=0;
    for(int mo=0; mo < nmo; mo++) { 
      sum_old+=movals_old(mo,0)*movals_old(mo,0);
      sum+=movals(mo,0)*movals(mo,0);
    }
    if(rng.ulec() < sum/sum_old) { 
      movals_old=movals;
      rold=r;
    }
  }

  movals=movals_old;
  sample->setElectronPos(e,rold);

}

