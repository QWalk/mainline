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
  //orbital normalization, then 1bdm up, 1bdm down, 2bdm up up, 2bdm up down 2bdm down up 2bdm down down, with everything complex
  if(eval_tbdm) { 
    avg.vals.Resize(nmo+4*nmo*nmo+8*nmo*nmo);
  }
  else { 
    avg.vals.Resize(nmo+(4)*nmo*nmo);  
  }
  avg.vals=0;


  wf->updateVal(wfdata,sample);

  Array2 <dcomplex> movals1(nmo,1),movals2(nmo,1),movals1_old(nmo,1),movals2_old(nmo,1);
  Wf_return wfval_base(wf->nfunc(),2);
  Wf_return wfval_2b(wf->nfunc(),2),wfval_1b(wf->nfunc(),2); //
  wf->getVal(wfdata,0,wfval_base);
  int nup=sys->nelectrons(0);
  int ndown=sys->nelectrons(1);
  int nupup=0,nupdown=0,ndownup=0,ndowndown=0;
  for(int i=0; i< npoints_eval  ;i++) { 
    Array1 <doublevar> r1(3),r2(3),oldr1(3),oldr2(3);
    int k=0,l=0;
    while(k==l) { 
      k=int(rng.ulec()*(nup+ndown));
      l=int(rng.ulec()*(nup+ndown));
      if(nup==1 and ndown==1) { 
        k=0; l=1;
      }
      else if(nup==1 or ndown==1) { 
        error("Need to fix density_matrix");
      }
      else { 
        if(i%4==0) { 
          k=int(rng.ulec()*nup);
          l=int(rng.ulec()*nup);
        }
        else if(i%4==1) { 
          k=int(rng.ulec()*nup);
          l=nup+int(rng.ulec()*ndown);
        }
        else if(i%4==2) { 
          k=nup+int(rng.ulec()*ndown);
          l=int(rng.ulec()*nup);
        }
        else if(i%4==3) { 
          k=nup+int(rng.ulec()*ndown);
          l=nup+int(rng.ulec()*ndown);
        }
      }
    }
    //Calculate the orbital values for r1 and r2
    calc_mos(sample,k,movals1_old);
    calc_mos(sample,l,movals2_old);
    sample->getElectronPos(k,oldr1);
    sample->getElectronPos(l,oldr2);

    r1=saved_r(i);
    r2=saved_r(npoints_eval+i);
    sample->setElectronPos(k,r1);
    doublevar dist1=gen_sample(nstep_sample,1.0,k,movals1, sample);

    wf->updateVal(wfdata,sample);
    wf->getVal(wfdata,0,wfval_1b);


    dcomplex psiratio_1b=exp(dcomplex(wfval_1b.amp(0,0)-wfval_base.amp(0,0),
          wfval_1b.phase(0,0)-wfval_base.phase(0,0)));
   
    doublevar dist2; dcomplex psiratio_2b;
    if(eval_tbdm) { 

      sample->setElectronPos(l,r2); 
      dist2=gen_sample(nstep_sample,1.0,l,movals2,sample);
      wf->updateVal(wfdata,sample);
      wf->getVal(wfdata,0,wfval_2b);

      psiratio_2b=exp(dcomplex(wfval_2b.amp(0,0)-wfval_base.amp(0,0),
            wfval_2b.phase(0,0)-wfval_base.phase(0,0)));
    }

    //orbital normalization
    for(int orbnum=0; orbnum < nmo; orbnum++) { 
      //avg.vals(orbnum)+=0.5*(norm(movals1(orbnum,0))/dist1
      //      +norm(movals2(orbnum,0))/dist2 )/npoints_eval;
      avg.vals(orbnum)+=norm(movals1(orbnum,0))/(dist1*npoints_eval);

    }
    int which_tbdm=0;
    int which_obdm=0;
    int npairs=0;
    int nelec_1b=nup;
    if(k >= nup) { which_obdm=1; nelec_1b=ndown; } 
    if(k < nup and l < nup) { which_tbdm=2; npairs=nup*(nup-1);nupup++; } 
    else if(k < nup and l >= nup) { which_tbdm=3; npairs=nup*ndown; nupdown++;} 
    else if(k >= nup and l < nup) { which_tbdm=4; npairs=ndown*nup; ndownup++; } 
    else if(k >= nup and l >= nup) { which_tbdm=5; npairs=ndown*(ndown-1);ndowndown++; } 

    dcomplex tmp;
    int place=0;
    for(int orbnum=0; orbnum < nmo; orbnum++) { 
      for(int orbnum2=0; orbnum2 < nmo; orbnum2++) { 
        tmp=doublevar(nelec_1b)*movals1(orbnum,0)*conj(movals1_old(orbnum2,0))
            *psiratio_1b/dist1;

        avg.vals(nmo+2*which_obdm*nmo*nmo+place)+=tmp.real();
        avg.vals(nmo+2*which_obdm*nmo*nmo+place+1)+=tmp.imag();

        place+=2;
      }
    }
    if(eval_tbdm) { 
      place=0;
      for(int orbnum=0; orbnum < nmo; orbnum++) { 
        for(int orbnum2=0; orbnum2 < nmo; orbnum2++) { 
          tmp=doublevar(npairs)*conj(movals1(orbnum,0))*movals1_old(orbnum2,0)
            *conj(movals2(orbnum,0))*movals2_old(orbnum2,0)
            *psiratio_2b/dist1/dist2;
          avg.vals(nmo+2*which_tbdm*nmo*nmo+place)+=tmp.real();
          avg.vals(nmo+2*which_tbdm*nmo*nmo+place+1)+=tmp.imag();
          place+=2;
        }
      }
    }
    
    //Restore the electronic positions
    sample->getElectronPos(k,saved_r(i));
    sample->getElectronPos(l,saved_r(npoints_eval+i));
    sample->setElectronPos(k,oldr1);
    sample->setElectronPos(l,oldr2);
  }

  int place=0;
  for(int i=0;i < nmo; i++) { 
    for(int j=0; j<nmo; j++) { 
      for(int k=0; k< 2; k++) { 
        avg.vals(nmo+2*0*nmo*nmo+place)/=nupup+nupdown;
        avg.vals(nmo+2*1*nmo*nmo+place)/=ndownup+ndowndown;
        if(eval_tbdm) { 
          avg.vals(nmo+2*2*nmo*nmo+place)/=nupup;
          avg.vals(nmo+2*3*nmo*nmo+place)/=nupdown;
          avg.vals(nmo+2*4*nmo*nmo+place)/=ndownup;
          avg.vals(nmo+2*5*nmo*nmo+place)/=ndowndown;
        }
        place++;
      }
    }
  }
  //cout << nupup << " " << nupdown << " " << ndownup << " " << ndowndown << endl;


}

//----------------------------------------------------------------------

void Average_tbdm_basis::read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words) { 
  unsigned int pos=0;
  vector <string> orbs;

  vector <string> mosec;
  if(readsection(words,pos=0, mosec,"ORBITALS")) { 
    //error("Need ORBITALS section in TBDM_BASIS");
    complex_orbitals=false;
    allocate(mosec,sys,momat);
    nmo=momat->getNmo();
  }
  else if(readsection(words,pos=0,mosec,"CORBITALS")) { 
    complex_orbitals=true;
    allocate(mosec,sys,cmomat);
    nmo=cmomat->getNmo();
  }
  else { error("Need ORBITALS or CORBITALS in TBDM_BASIS"); } 

  Array1 <Array1 <int> > occupations(1);
  if(readsection(words,pos=0,orbs,"EVAL_ORBS")) { 
    nmo=orbs.size();
    occupations[0].Resize(nmo);
    for(int i=0; i< nmo; i++) 
      occupations[0][i]=atoi(orbs[i].c_str());
  }
  else { 
    occupations[0].Resize(nmo);
    for(int i=0; i< nmo; i++) { 
      occupations[0][i]=i;
    }
  }

  if(complex_orbitals) 
    cmomat->buildLists(occupations);
  else 
    momat->buildLists(occupations);



  if(!readvalue(words, pos=0,nstep_sample,"NSTEP_SAMPLE"))
    nstep_sample=10;
  if(!readvalue(words,pos=0,npoints_eval,"NPOINTS"))
    npoints_eval=100;

  eval_tbdm=true;
  if(haskeyword(words,pos=0,"ONLY_OBDM"))
    eval_tbdm=false;

  //Since we rotate between the different pairs of spin channels, make sure
  //that npoints_eval is divisible by 4
  if(npoints_eval%4!=0) {
    npoints_eval+=4-npoints_eval%4;
  }

  int ndim=3;
  saved_r.Resize(npoints_eval*2); //r1 and r2;
  for(int i=0; i< npoints_eval*2; i++) saved_r(i).Resize(ndim);
  
  int warmup_steps=1000;
  Sample_point * sample=NULL;
  sys->generateSample(sample);
  sample->randomGuess();
  Array2 <dcomplex> movals(nmo,1);
  for(int i=0; i< npoints_eval*2; i++) { 
    gen_sample(warmup_steps,1.0,0,movals,sample);
    sample->getElectronPos(0,saved_r(i));
  }
}

//----------------------------------------------------------------------

void Average_tbdm_basis::write_init(string & indent, ostream & os) { 
  os << indent << "TBDM_BASIS" << endl;
  os << indent << "NMO " << nmo << endl;
  os << indent << "NPOINTS " << npoints_eval << endl;
  if(!eval_tbdm) { 
    os << indent << "ONLY_OBDM" << endl;
  }
  if(complex_orbitals) { 
    os << indent << "CORBITALS { \n";
    cmomat->writeinput(indent,os); 
  }
  else { 
    os << indent << "ORBITALS { \n";
    momat->writeinput(indent,os);
  }
  os << indent << "}\n";
}
//----------------------------------------------------------------------
void Average_tbdm_basis::read(vector <string> & words) { 
  unsigned int pos=0;
  /*
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
  */

  readvalue(words, pos=0,nmo,"NMO");
  readvalue(words,pos=0,npoints_eval,"NPOINTS");
  if(haskeyword(words, pos=0,"ONLY_OBDM")) eval_tbdm=false;
  else eval_tbdm=true;
}
//----------------------------------------------------------------------
void Average_tbdm_basis::write_summary(Average_return &avg,Average_return &err, ostream & os) { 

  Array2 <dcomplex> obdm_up(nmo,nmo),obdm_down(nmo,nmo), 
         obdm_up_err(nmo,nmo),obdm_down_err(nmo,nmo),
         tbdm_uu(nmo,nmo),tbdm_uu_err(nmo,nmo),
         tbdm_ud(nmo,nmo),tbdm_ud_err(nmo,nmo),
         tbdm_du(nmo,nmo),tbdm_du_err(nmo,nmo),
         tbdm_dd(nmo,nmo),tbdm_dd_err(nmo,nmo);

  os << "tbdm: nmo " << nmo << endl;
  os << "Orbital normalization " << endl;
  for(int i=0; i< nmo; i++) { 
    os << avg.vals(i) << " +/- " << err.vals(i) << endl;
  }

  dcomplex i_c(0.,1.);
  int place=0;
  for(int i=0; i < nmo; i++) { 
    for(int j=0; j < nmo; j++)  { 
      doublevar norm=sqrt(avg.vals(i)*avg.vals(j));
      int index=nmo+place;
      obdm_up(i,j)=dcomplex(avg.vals(index),avg.vals(index+1))/norm;
      obdm_up_err(i,j)=dcomplex(err.vals(index),err.vals(index+1))/norm;

      index+=2*nmo*nmo;
      obdm_down(i,j)=dcomplex(avg.vals(index),avg.vals(index+1))/norm;
      obdm_down_err(i,j)=dcomplex(err.vals(index),err.vals(index+1))/norm;
      if(eval_tbdm) { 
        index+=2*nmo*nmo;
        tbdm_uu(i,j)=dcomplex(avg.vals(index),avg.vals(index+1))/(norm*norm);
        tbdm_uu_err(i,j)=dcomplex(err.vals(index),err.vals(index+1))/(norm*norm);
        index+=2*nmo*nmo;
        tbdm_ud(i,j)=dcomplex(avg.vals(index),avg.vals(index+1))/(norm*norm);
        tbdm_ud_err(i,j)=dcomplex(err.vals(index),err.vals(index+1))/(norm*norm);
        index+=2*nmo*nmo;
        tbdm_du(i,j)=dcomplex(avg.vals(index),avg.vals(index+1))/(norm*norm);
        tbdm_du_err(i,j)=dcomplex(err.vals(index),err.vals(index+1))/(norm*norm);
        index+=2*nmo*nmo;
        tbdm_dd(i,j)=dcomplex(avg.vals(index),avg.vals(index+1))/(norm*norm);
        tbdm_dd_err(i,j)=dcomplex(err.vals(index),err.vals(index+1))/(norm*norm);
      }
      place+=2;
    }
  }

  os << "One-body density matrix " << endl;
  int colwidth=40;
  os << setw(10) << " " << setw(10) << " "
    << setw(colwidth) << "up" << setw(colwidth) << "up err"
    << setw(colwidth) << "down" << setw(colwidth) << "down err"
    << endl;
  for(int i=0; i< nmo ; i++) { 
    for(int j=0; j<nmo; j++) { 
      os << setw(10) << i << setw(10) << j
        << setw(colwidth) << obdm_up(i,j) << setw(colwidth) << obdm_up_err(i,j)
        << setw(colwidth) << obdm_down(i,j) << setw(colwidth) << obdm_down_err(i,j)
        << endl;
    }

  }
  if(eval_tbdm) { 
    os << "two-body density matrix " << endl;
    os << setw(10) << " " << setw(10) << " "
      << setw(colwidth) << "upup" << setw(colwidth) << "upup err"
      << setw(colwidth) << "updown" << setw(colwidth) << "updown err"
      << setw(colwidth) << "downup" << setw(colwidth) << "downup err"
      << setw(colwidth) << "downdown" << setw(colwidth) << "downdown err"
      << endl;
    for(int i=0; i< nmo ; i++) { 
      for(int j=0; j<nmo; j++) { 
        os << setw(10) << i << setw(10) << j
          << setw(colwidth) << tbdm_uu(i,j) << setw(colwidth) << tbdm_uu_err(i,j)
          << setw(colwidth) << tbdm_ud(i,j) << setw(colwidth) << tbdm_ud_err(i,j)
          << setw(colwidth) << tbdm_du(i,j) << setw(colwidth) << tbdm_du_err(i,j)
          << setw(colwidth) << tbdm_dd(i,j) << setw(colwidth) << tbdm_dd_err(i,j)
          << endl;
      }

    }
  }


  dcomplex trace=0;
  for(int i=0;i< nmo; i++) trace+=obdm_up(i,i);
  os << "Trace of the obdm: up: " << trace;
  trace=0;
  for(int i=0; i< nmo; i++) trace+=obdm_down(i,i);
  os << " down: " << trace << endl;


}

//----------------------------------------------------------------------
//
//
//Note this needs to be changed for non-zero k-points!
doublevar Average_tbdm_basis::gen_sample(int nstep, doublevar  tstep, 
    int e, Array2 <dcomplex> & movals, Sample_point * sample) { 
  int ndim=3;
  Array1 <doublevar> r(ndim),rold(ndim);
  Array2 <dcomplex> movals_old(nmo,1);
  movals.Resize(nmo,1);

  sample->getElectronPos(e,rold);
  calc_mos(sample,e,movals_old);
  doublevar acc=0;

  for(int step=0; step < nstep; step++) { 
    for(int d=0; d< ndim; d++) { 
      r(d)=rold(d)+sqrt(tstep)*rng.gasdev();
    }
    sample->setElectronPos(e,r);
    calc_mos(sample,e,movals); 

    doublevar sum_old=0,sum=0;
    for(int mo=0; mo < nmo; mo++) { 
      sum_old+=norm(movals_old(mo,0));
      sum+=norm(movals(mo,0));
    }
    if(rng.ulec() < sum/sum_old) { 
      movals_old=movals;
      rold=r;
      acc++;
    }
  }
  //cout << "acceptance " << acc/nstep << endl;

  movals=movals_old;
  sample->setElectronPos(e,rold);

  doublevar sum=0;
  for(int mo=0; mo < nmo; mo++) {
    sum+=norm(movals(mo,0));
  }
  return sum;

}

//----------------------------------------------------------------------
void Average_tbdm_basis::calc_mos(Sample_point * sample, int e, Array2 <dcomplex> & movals) { 
  if(complex_orbitals) { 
    movals.Resize(nmo,1);
    cmomat->updateVal(sample,e,0,movals); 
  }
  else { 
    Array2 <doublevar> movals_d(nmo,1);
    momat->updateVal(sample,e,0,movals_d);
    movals.Resize(nmo,1);
    for(int i=0; i< nmo; i++) { 
      movals(i,0)=movals_d(i,0);
    }
  }
}

//----------------------------------------------------------------------
