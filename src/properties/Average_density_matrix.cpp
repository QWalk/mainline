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
  avg.vals.Resize(nmo+6*nmo*nmo); //orbital normalization, then 1bdm up, 1bdm down, 2bdm up up, 2bdm up down 2bdm down up 2bdm down down
  avg.vals=0;
  wf->updateVal(wfdata,sample);

  Array2 <doublevar> movals1(nmo,1),movals2(nmo,1),movals1_old(nmo,1),movals2_old(nmo,1);
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
    int which_tbdm=0;
    int which_obdm=0;
    int npairs=0;
    int nelec_1b=nup;
    if(k >= nup) { which_obdm=1; nelec_1b=ndown; } 
    if(k < nup and l < nup) { which_tbdm=2; npairs=nup*(nup-1);nupup++; } 
    else if(k < nup and l >= nup) { which_tbdm=3; npairs=nup*ndown; nupdown++;} 
    else if(k >= nup and l < nup) { which_tbdm=4; npairs=ndown*nup; ndownup++; } 
    else if(k >= nup and l >= nup) { which_tbdm=5; npairs=ndown*(ndown-1);ndowndown++; } 
    for(int orbnum=0; orbnum < nmo; orbnum++) { 
      for(int orbnum2=0; orbnum2 < nmo; orbnum2++) { 
        avg.vals(nmo+which_tbdm*nmo*nmo+orbnum*nmo+orbnum2)+=
            npairs*(movals1(orbnum,0)*movals1_old(orbnum,0)
            *movals2(orbnum2,0)*movals2_old(orbnum2,0) 
            *psiratio_2b)/dist1/dist2 ;

        avg.vals(nmo+which_obdm*nmo*nmo+orbnum*nmo+orbnum2)+=
          nelec_1b*movals1(orbnum,0)*movals1_old(orbnum2,0)
          *psiratio_1b/dist1;
      }

    }
    //Restore the electronic positions
    sample->setElectronPos(k,oldr1);
    sample->setElectronPos(l,oldr2);
  }

  for(int i=0;i < nmo; i++) { 
    for(int j=0; j<nmo; j++) { 
      int place=i*nmo+j;
      avg.vals(nmo+0*nmo*nmo+place)/=nupup+nupdown;
      avg.vals(nmo+1*nmo*nmo+place)/=ndownup+ndowndown;
      avg.vals(nmo+2*nmo*nmo+place)/=nupup;
      avg.vals(nmo+3*nmo*nmo+place)/=nupdown;
      avg.vals(nmo+4*nmo*nmo+place)/=ndownup;
      avg.vals(nmo+5*nmo*nmo+place)/=ndowndown;
    }
  }
  //cout << nupup << " " << nupdown << " " << ndownup << " " << ndowndown << endl;


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

  //Since we rotate between the different pairs of spin channels, make sure
  //that npoints_eval is divisible by 4
  if(npoints_eval%4!=0) {
    npoints_eval+=4-npoints_eval%4;
  }

  
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
  Array2 <doublevar> obdm_up(nmo,nmo),obdm_down(nmo,nmo), 
         obdm_up_err(nmo,nmo),obdm_down_err(nmo,nmo),
         tbdm_uu(nmo,nmo),tbdm_uu_err(nmo,nmo),
         tbdm_ud(nmo,nmo),tbdm_ud_err(nmo,nmo),
         tbdm_du(nmo,nmo),tbdm_du_err(nmo,nmo),
         tbdm_dd(nmo,nmo),tbdm_dd_err(nmo,nmo);
  for(int i=0; i < nmo; i++) { 
    for(int j=0; j < nmo; j++)  { 
      doublevar norm=sqrt(avg.vals(i)*avg.vals(j));
      int place=i*nmo+j;
      obdm_up(i,j)=avg.vals(nmo+place)/norm;
      obdm_up_err(i,j)=err.vals(nmo+place)/norm;
      obdm_down(i,j)=avg.vals(nmo+nmo*nmo+place)/norm;
      obdm_down_err(i,j)=err.vals(nmo+nmo*nmo+place)/norm;
      tbdm_uu(i,j)=avg.vals(nmo+2*nmo*nmo+place)/(norm*norm);
      tbdm_uu_err(i,j)=err.vals(nmo+2*nmo*nmo+place)/(norm*norm);
      tbdm_ud(i,j)=avg.vals(nmo+3*nmo*nmo+place)/(norm*norm);
      tbdm_ud_err(i,j)=err.vals(nmo+3*nmo*nmo+place)/(norm*norm);
      tbdm_du(i,j)=avg.vals(nmo+4*nmo*nmo+place)/(norm*norm);
      tbdm_du_err(i,j)=err.vals(nmo+4*nmo*nmo+place)/(norm*norm);
      tbdm_dd(i,j)=avg.vals(nmo+5*nmo*nmo+place)/(norm*norm);
      tbdm_dd_err(i,j)=err.vals(nmo+5*nmo*nmo+place)/(norm*norm);
    }
  }

  os << "One-body density matrix " << endl;
  os << setw(10) << " " << setw(10) << " "
    << setw(17) << "up" << setw(17) << "up err"
    << setw(17) << "down" << setw(17) << "down err"
    << endl;
  for(int i=0; i< nmo ; i++) { 
    for(int j=0; j<nmo; j++) { 
      os << setw(10) << i << setw(10) << j
        << setw(17) << obdm_up(i,j) << setw(17) << obdm_up_err(i,j)
        << setw(17) << obdm_down(i,j) << setw(17) << obdm_down_err(i,j)
        << endl;
    }

  }

  os << "two-body density matrix " << endl;
  os << setw(10) << " " << setw(10) << " "
    << setw(17) << "upup" << setw(17) << "upup err"
    << setw(17) << "updown" << setw(17) << "updown err"
    << setw(17) << "downup" << setw(17) << "downup err"
    << setw(17) << "downdown" << setw(17) << "downdown err"
    << endl;
  for(int i=0; i< nmo ; i++) { 
    for(int j=0; j<nmo; j++) { 
      os << setw(10) << i << setw(10) << j
        << setw(17) << tbdm_uu(i,j) << setw(17) << tbdm_uu_err(i,j)
        << setw(17) << tbdm_ud(i,j) << setw(17) << tbdm_ud_err(i,j)
        << setw(17) << tbdm_du(i,j) << setw(17) << tbdm_du_err(i,j)
        << setw(17) << tbdm_dd(i,j) << setw(17) << tbdm_dd_err(i,j)
        << endl;
    }

  }

  os << "Two-body density matrix minus product of one-body density matrices\n";
  os << setw(10) << " " << setw(10) << " "
    << setw(17) << "upup" << setw(17) << "upup err"
    << setw(17) << "updown" << setw(17) << "updown err"
    << setw(17) << "downup" << setw(17) << "downup err"
    << setw(17) << "downdown" << setw(17) << "downdown err"
    << endl;
  for(int i=0; i< nmo ; i++) { 
    for(int j=0; j<nmo; j++) { 
      os << setw(10) << i << setw(10) << j
        << setw(17) << tbdm_uu(i,j)-obdm_up(i,i)*obdm_up(j,j) << setw(17) << tbdm_uu_err(i,j)
        << setw(17) << tbdm_ud(i,j)-obdm_up(i,i)*obdm_down(j,j) << setw(17) << tbdm_ud_err(i,j)
        << setw(17) << tbdm_du(i,j)-obdm_down(i,i)*obdm_up(j,j) << setw(17) << tbdm_du_err(i,j)
        << setw(17) << tbdm_dd(i,j)-obdm_down(i,i)*obdm_down(j,j) << setw(17) << tbdm_dd_err(i,j)
        << endl;
    }
  }


  doublevar trace=0;
  for(int i=0;i< nmo; i++) trace+=obdm_up(i,i);
  os << "Trace of the obdm: up: " << trace;
  trace=0;
  for(int i=0; i< nmo; i++) trace+=obdm_down(i,i);
  os << " down: " << trace << endl;
  trace=0;
  for(int i=0; i< nmo; i++) trace+=tbdm_uu(i,i);
  os << "Trace of tbdm: upup: " << trace;
  trace=0;
  for(int i=0; i< nmo; i++) trace+=tbdm_ud(i,i);
  os << " updown: " << trace;
  trace=0;
  for(int i=0; i< nmo; i++) trace+=tbdm_du(i,i);
  os << " downup: " << trace;
  trace=0;
  for(int i=0; i< nmo; i++) trace+=tbdm_dd(i,i);
  os << " downdown: " << trace << endl;


  os << "OBDM_up" << endl;
  for(int i=0; i< nmo; i++) { 
    for(int j=0; j< nmo; j++)  { 
      os << setw(17) << obdm_up(i,j);
    }
    os << endl;
  }

  os << "TBDM_upup" << endl;
  for(int i=0; i< nmo; i++) { 
    for(int j=0; j< nmo; j++) { 
      os << setw(17) << tbdm_uu(i,j);
    }
    os << endl;
  }

  os << "TBDM_updown" << endl;
  for(int i=0; i< nmo; i++) { 
    for(int j=0; j< nmo; j++) { 
      os << setw(17) << tbdm_ud(i,j);
    }
    os << endl;
  }

  

  /*
  os << "Extra correlation" << endl;
  for(int i=0; i< nmo; i++) { 
    for(int j=0; j< nmo; j++) { 
      doublevar tmp= tbdm(i,j)-obdm(i,i)*obdm(j,j);
      doublevar tmp_err= tbdm_err(i,j)+obdm_err(i,i)+obdm_err(j,j); //wrong error!
      if(fabs(tmp) > 3*tmp_err) { 
        os << setw(10) << i << setw(10) << j
          << setw(17) << tmp
          << setw(17) << tmp_err
          << endl;
      }
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
*/
/*

  os << "temporary matrix output " << endl;
  vector <string> orbnames(13);
  orbnames[0]="Mn2s";
  orbnames[1]="Mnz";
  orbnames[2]="Mny+x";
  orbnames[3]="Mny-x";
  orbnames[4]="Os";
  orbnames[5]="Oy+x";
  orbnames[6]="Oy-x";
  orbnames[7]="Oz";
  orbnames[8]="Mnxy";
  orbnames[9]="Mnxx-yy";
  orbnames[10]="Mnzz-rr";
  orbnames[11]="Mnyz+xz";
  orbnames[12]="Mnyz-xz";

  os << setw(10) << " ";
  for(int i=0; i< 13; i++) os << setw(10) << orbnames[i];
  os << endl;
  for(int i=0; i< 13; i++) {
    os << setw(10) << orbnames[i];
    for(int j=0; j< 13; j++) { 
      if(fabs(obdm(i,j)) > 0.01) 
        os << setw(10) << setprecision(2) << obdm(i,j);
      else 
        os << setw(10) << 0.0;
    }
    os << endl;
  }


  Array2 <doublevar> obdm_tmp(13,13);
  for(int i=0; i< 13; i++) 
    for(int j=0; j< 13; j++) obdm_tmp(i,j)=obdm(i,j);
  EigenSystemSolverRealSymmetricMatrix(obdm,evals,evecs);
  nmo=13;
  os << "evals ";
  for(int i=0; i< nmo; i++) os << evals(i) <<" " ;
  os << endl;


  os << "Evecs " << endl;
  for(int i=0; i< nmo; i++) { 
    os << setw(10) << orbnames[i];
    for(int j=0; j< nmo; j++) { 
      os << setw(10) << evecs(i,j);
    }
    os << endl;
  }
*/

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

