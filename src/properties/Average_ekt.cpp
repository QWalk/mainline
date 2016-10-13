/*

   Copyright (C) 2015 Huihuo Zheng

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

#include "Average_ekt.h"
#include "ulec.h"
#include "Pseudopotential.h"
//#include "System.h"
void Average_ekt::randomize(Wavefunction_data * wfdata, Wavefunction * wf,
    System * sys, Sample_point * sample_tmp) { 
  int nup=sys->nelectrons(0);
  int ndown=sys->nelectrons(1);
  totnelectrons = nup + ndown; 
  Array2 <dcomplex> movals1(nmo,1),movals2(nmo,1);
  //Make a copy of the Sample so that it doesn't have to update the wave function.
  Sample_point * sample;
  sys->generateSample(sample);
  for(int i=0; i< npoints_eval; i++) { 
    int k=0,l=0;
    while(k==l) { 
      k=int(rng.ulec()*(nup+ndown));
      l=int(rng.ulec()*(nup+ndown));
      if(nup==1 and ndown==1) { 
        k=0; l=1;
      }
      else if(nup==1) { 
        if(i%2==0) { 
          k=0;
          l=nup+int(rng.ulec()*ndown);
        }
        else { 
          k=nup+int(rng.ulec()*ndown);
          l=nup+int(rng.ulec()*ndown);
        }
      }
      else if(ndown==1) { 
        if(i%2==0) { 
          k=int(rng.ulec()*nup);
          l=nup;
        }
        else { 
          k=int(rng.ulec()*nup);
          l=int(rng.ulec()*nup);
        }
      }
      else if(ndown==0) { 
        k=int(rng.ulec()*nup);
        l=int(rng.ulec()*nup);
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
    rk(i)=k;
    rk(npoints_eval+i)=l;
    Array1 <doublevar> r1=saved_r(i);
    Array1 <doublevar> r2=saved_r(npoints_eval+i);
    sample->setElectronPos(k,r1);
    sample->setElectronPos(l,r2);
    doublevar dist1=gen_sample(nstep_sample,1.0,k,movals1, sample);
    doublevar dist2=gen_sample(nstep_sample,1.0,l,movals2, sample);
    sample->getElectronPos(k,saved_r(i));
    sample->getElectronPos(l,saved_r(npoints_eval+i));

  }


  delete sample; 
}

//----------------------------------------------------------------------

void Average_ekt::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
    System * sys, Pseudopotential *psp, Sample_point * sample, Properties_point & pt, Average_return & avg) {
  evaluate(wfdata, wf, sys, psp, sample, avg);
  // for ekt, we do not need pt. 
}

void Average_ekt::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
    System * sys, Pseudopotential * psp, Sample_point * sample, Average_return & avg) { 
  /*
     For our first testing case, we implemented the three matrix independently, but in principle they can be grouped 
     together to reduced the computation cost. 
     */
  int psp_curr_deterministic=psp->getDeterministic();
  psp->setDeterministic(deterministic_psp);
  //int n = eval_conduction + eval_valence + eval_obdm; 
  int nelectrons = sys->nelectrons(0) + sys->nelectrons(1);
  int nup = sys->nelectrons(0); 
  int ndown = sys->nelectrons(1); 
  totnelectrons = nup + ndown; 
  //  avg.vals.Resize(nmo + 3*4*nmo*nmo + nelectrons + nelectrons); 
  avg.vals.Resize(nmo + 12*nmo*nmo); 
  avg.vals = 0.0; 
  avg.type="EKT";
  //  cout << "Test" << endl; 
  //  cout << "Evaluating EKT" << endl; 
  //if(eval_conduction)
  //evaluate_conduction(wfdata,wf,sys, psp, sample,avg); 
  evaluate_valence(wfdata,wf,sys, psp, sample,avg);
  //  evaluate_obdm(wfdata,wf,sys, sample, avg);
  psp->setDeterministic(psp_curr_deterministic);

}



void Average_ekt::evaluate_conduction(Wavefunction_data * wfdata, Wavefunction * wf,
    System * sys, Pseudopotential *psp, Sample_point * sample, Average_return & avg) {
  cout << "Not implemented " << endl; 
}


//----------------------------------------------------------------------

void Average_ekt::evaluate_valence(Wavefunction_data * wfdata, Wavefunction * wf,
    System * sys, Pseudopotential *psp, Sample_point * sample, Average_return & avg) { 
  /*
     This is to evaluate the v_ij in extended Koopmans' theorem -- for valence bands
     */
  wf->updateVal(wfdata,sample);
  Wf_return wfval_base(wf->nfunc(),2);
  wf->getVal(wfdata,0,wfval_base);
  int nup=sys->nelectrons(0);
  int ndown=sys->nelectrons(1);
  int nelectrons=nup+ndown;

  Array1 <doublevar> rn(3), rnp(3), rnpp(3); 
  Array1 <Array2 <dcomplex> > movals1_base(nelectrons);
  Array2 <dcomplex> movals_p(nmo, 1); 
  Array1 <doublevar> VLoc(nelectrons); 
  Array1 <doublevar> VLoc0(nelectrons);
  Array1 <doublevar> Vtest(nelectrons+1);
  Array1 <dcomplex> pseudo_t(nmo); 
  /***************************************
    Routines to get Kin and Vloc
   ****************************************/
  Array2<doublevar> Kin(nelectrons, wf->nfunc());
  //Array2<doublevar> Kin0(nelectrons, wf->nfunc());
  Array2<dcomplex> movals_lap(nmo, 5);
  for(int e=0; e< nelectrons; e++) { 
    movals1_base(e).Resize(nmo,1);
    //Kin(e).Resize(1); //The Kinetic energy 
    calc_mos(sample, e, movals1_base(e));
  }
//  sys->calcKineticSeparated(sample,saved_r(i),Kin);

  //***** calculate kinetic energy and potential energy
  sys->calcKineticSeparated(wfdata, sample, wf, Kin); 
  sys->calcLocSeparated(sample, VLoc);

  Array1 <Wf_return> wfs(nelectrons);

  Wavefunction_storage * store;
  wf->generateStorage(store);
  for(int i=0; i< npoints_eval; i++) {
    Array1 <doublevar> oldpos(3);
    sample->getElectronPos(0,oldpos);
    sample->setElectronPosNoNotify(0, saved_r(i));
    calc_mosLap(sample, 0, movals_lap);
    int nrandvar=psp->nTest();
    Array1 <doublevar> rand_num(nrandvar);
    for(int l=0; l< nrandvar; l++) 
      rand_num(l)=rng.ulec();
    calc_mos(sample,0,movals_p);

    calcPseudoMo(sys, sample, psp, rand_num, pseudo_t);
    sample->setElectronPosNoNotify(0,oldpos);
    Array1 <Wf_return> wf_eval;
    wf->evalTestPos(saved_r(i),sample,wfs);
    sys->calcLocWithTestPos(sample, saved_r(i), Vtest);

    doublevar vtot = 0.0;
    for (int e=0; e<nelectrons; e++) {
      vtot += Vtest(e);
    }
    doublevar dist1=0;
    for(int m=0; m < nmo; m++)
      dist1+=norm(movals_p(m,0));
    
    for(int orbnum=0; orbnum < nmo; orbnum++) {
      avg.vals(orbnum)+=norm(movals_p(orbnum,0))/(dist1*npoints_eval);//this is the normalization
    }
    
    Array2<doublevar> totalv(nelectrons, wf->nfunc());
    totalv = 0.0;
    psp->calcNonlocSeparated(wfdata, sys, sample, wf, totalv);
    
    int place = 0;
    for (int orbnum = 0; orbnum < nmo; orbnum++ ) {
      for (int orbnum2 = 0; orbnum2 < nmo; orbnum2++ ) {
        dcomplex tmp3=conj(movals_p(orbnum, 0))*
                    (pseudo_t(orbnum2) +
                     Vtest(nelectrons)*movals_p(orbnum2, 0)
                      -0.5*movals_lap(orbnum2, 4)
                     +movals_p(orbnum2, 0)*vtot)
                /(dist1*npoints_eval);
        //	tmp3=conj(movals_p(orbnum, 0))*(movals_p(orbnum2, 0))/(dist1*npoints_eval);
        //dcomplex tmp3=conj(movals_p(orbnum, 0))*(pseudo_t(orbnum2))/(dist1*npoints_eval);
        int which_obdm = 0;
        avg.vals(nmo+8*nmo*nmo + 2*which_obdm*nmo*nmo+place) += tmp3.real();
        avg.vals(nmo+8*nmo*nmo + 2*which_obdm*nmo*nmo+place+1) += tmp3.imag();
        
        which_obdm = 1;
        //tmp3=conj(movals_p(orbnum, 0))*(-0.5*movals_lap(orbnum2, 4))/(dist1*npoints_eval);
        avg.vals(nmo+8*nmo*nmo + 2*which_obdm*nmo*nmo+place) += tmp3.real();
        avg.vals(nmo+8*nmo*nmo + 2*which_obdm*nmo*nmo+place+1) += tmp3.imag();
        place += 2;
      }
    }
    ofstream dump;
    if(dump_data) {
      dump.open("EKT_DUMP",ios::app);
    }
    for(int e=0; e< nelectrons; e++) {
      dcomplex psiratio_1b=conj(exp(dcomplex(wfs(e).amp(0,0)
                                             -wfval_base.amp(0,0),
                                             wfs(e).phase(0,0)
                                             -wfval_base.phase(0,0))));
      int which_obdm=0;
      if(e >= nup) { which_obdm=1;  }
      dcomplex tmp, tmp2, tmp3;
      int place=0;
      dcomplex prefactor=psiratio_1b/(dist1*npoints_eval);
      //cout << "Local potential" << "   Kinetic" << "   Pseudo" << endl;
      //      cout << VLoc(e) << "  " << Kin(e, 0) << "  " << totalv(e, 0) << endl;
      
      for(int orbnum=0; orbnum < nmo; orbnum++) {
        for(int orbnum2=0; orbnum2 < nmo; orbnum2++) {
          //	  assert(wf->nfunc()==1);
          tmp = 0.5*(movals_p(orbnum,0)*conj(movals1_base(e)(orbnum2,0))*prefactor
                     + movals1_base(e)(orbnum, 0)*conj(movals_p(orbnum2, 0))*conj(prefactor));
          tmp2 = tmp*(VLoc(e) + Kin(e, 0) + totalv(e, 0));
          //rho_ij part the one body reduced density matrix
          doublevar tmp2_mag=abs(tmp2);
          if(tmp2_mag > ekt_cutoff) {
            tmp2*=ekt_cutoff/tmp2_mag;
          }
          avg.vals(nmo+2*which_obdm*nmo*nmo+place) += tmp.real();
          avg.vals(nmo+2*which_obdm*nmo*nmo+place+1) += tmp.imag();
          //v_ij^v part
          //tmp3 = 0.0;
          avg.vals(nmo+4*nmo*nmo + 2*which_obdm*nmo*nmo+place)+=tmp2.real();
          avg.vals(nmo+4*nmo*nmo + 2*which_obdm*nmo*nmo+place+1)+=tmp2.imag();
          if(dump_data) {
            if(orbnum==orbnum2 and orbnum==0) {
              sample->getElectronPos(e,oldpos);
              dump << which_obdm << "," << tmp2.real() << "," << tmp.real() << ","
                 << VLoc(e) << "," << Kin(e,0) << "," << totalv(e,0) <<
               "," << oldpos(0) << "," << oldpos(1) << "," << oldpos(2) << endl;
            }
          }
          if (eval_conduction) {
            //tmp3 = -1.0*psiratio_1b*conj(movals1_base(e)(orbnum, 0))*(vtot*movals_lap(orbnum2, 0)
            //						 + pseudo_t(orbnum2) - 0.5*movals_lap(orbnum2, 4) + Vtest(nelectrons)*movals_lap(orbnum2, 0))/(dist1*npoints_eval);
            tmp3 = -1.0*conj(movals1_base(e)(orbnum, 0))*movals_p(orbnum2, 0)*prefactor*(VLoc(e) + Kin(e, 0) + totalv(e, 0) + Vtest(e));
            doublevar tmp3_mag=abs(tmp3);
            if(tmp3_mag > ekt_cutoff) {
              tmp3*=ekt_cutoff/tmp3_mag;
            }
            avg.vals(nmo+8*nmo*nmo + 2*which_obdm*nmo*nmo+place) += tmp3.real();
            avg.vals(nmo+8*nmo*nmo + 2*which_obdm*nmo*nmo+place+1) += tmp3.imag();
          }
          place+=2;
        }
      }
    }
  }
  delete store;
}


void Average_ekt::evaluate_obdm(Wavefunction_data * wfdata, Wavefunction * wf,
    System * sys, Sample_point * sample, Average_return & avg) { 


  wf->updateVal(wfdata,sample);
  Wf_return wfval_base(wf->nfunc(),2);
  wf->getVal(wfdata,0,wfval_base);
  int nup=sys->nelectrons(0);
  int ndown=sys->nelectrons(1);
  int nelectrons=nup+ndown;

  Array1 <Array2 <dcomplex> > movals1_base(nelectrons);
  for(int e=0; e< nelectrons; e++) { 
    movals1_base(e).Resize(nmo,1);
    calc_mos(sample,e,movals1_base(e));
    /*! Huihuo
      Here, permutation symmetry has been applied, the real evaluation quantity is 
      sum_(e=1)^N phi_i* (r_e) phi_j(r') psi(r_1,... r',...r_N)/psi(r_1, ..., r_n, ..., r_N)
      Therefore, movals1_base[e, i] = phi_i(r_e), where e is the index of the electron. 
      !!!One need to be careful about the off-gamma point 1 RDM 
      */
  }
  //  avg.vals.Resize(nmo+4*nmo*nmo);
  //  avg.vals=0;

  Array2 <dcomplex> movals1(nmo,1);
  Array1 <Wf_return> wfs(nelectrons);

  Wavefunction_storage * store;
  wf->generateStorage(store);
  for(int i=0; i< npoints_eval; i++) { 
    Array1 <doublevar> oldpos(3);
    sample->getElectronPos(0,oldpos);
    sample->setElectronPosNoNotify(0,saved_r(i));
    calc_mos(sample,0,movals1);
    /*!
      This is simply to calculate phi_j(r'): Noted that r' has been given, therefore, movals1 is a one dimentional array, with matrix element to be 
      phi_j(r') -- j is the index 
      */
    sample->setElectronPosNoNotify(0,oldpos);

    Array1 <Wf_return> wf_eval;
    wf->evalTestPos(saved_r(i),sample,wfs);

    //Testing the evalTestPos
    //for(int e=0; e< nelectrons; e++) { 
    //  Wf_return test_wf(wf->nfunc(),2);
    //  sample->getElectronPos(e,oldpos);
    //  wf->saveUpdate(sample,e,store);
    //  sample->setElectronPos(e,saved_r(i));
    //  wf->updateVal(wfdata,sample);
    //  wf->getVal(wfdata,e,test_wf);
    //  sample->setElectronPos(e,oldpos);
    //  wf->restoreUpdate(sample,e,store);
    //  cout << "e " << e << " test " << test_wf.amp(0,0) << " evalTestPos " << wfs(e).amp(0,0) << endl;
    //}

    doublevar dist1=0;
    for(int m=0; m < nmo; m++) 
      dist1+=norm(movals1(m,0));

    for(int orbnum=0; orbnum < nmo; orbnum++) { 
      avg.vals(orbnum)+=norm(movals1(orbnum,0))/(dist1*npoints_eval);
    }

    for(int e=0; e< nelectrons; e++) { 
      dcomplex psiratio_1b=exp(dcomplex(wfs(e).amp(0,0)-wfval_base.amp(0,0),
            wfs(e).phase(0,0)-wfval_base.phase(0,0)));
      /*!
        psi(r_1,...,r',...,r_N)/psi(r_1,...,r_e,...,r_N), the ratio of the new with respect to the old if one change the position of the e-th electron from 
        r_e to r'
        */
      int which_obdm=0;
      if(e >= nup) { which_obdm=1;  } 
      dcomplex tmp;
      int place=0;
      dcomplex prefactor=psiratio_1b/(dist1*npoints_eval);
      for(int orbnum=0; orbnum < nmo; orbnum++) { 
        for(int orbnum2=0; orbnum2 < nmo; orbnum2++) { 
          //tmp=movals1(orbnum,0)*conj(movals1_base(e)(orbnum2,0))*prefactor;
          /*! Huihuo Zheng
            This is to based on my new definition: 
            \rho(r, r') = N\int dR psi*(r, R) psi(r', R) = \sum rho_{ij} phi_i*(r) phi_j(r')
            in this way, one will find that, the phase factor is cancelled
            */
          tmp=conj(movals1(orbnum,0))*movals1_base(e)(orbnum2,0)*prefactor;
          avg.vals(nmo+2*which_obdm*nmo*nmo+place)+=tmp.real();
          avg.vals(nmo+2*which_obdm*nmo*nmo+place+1)+=tmp.imag();

          place+=2;
        }
      }

    }
  }
  delete store;
}
//----------------------------------------------------------------------

void Average_ekt::read(System * sys, Wavefunction_data * wfdata, vector
    <string> & words) { 
  unsigned int pos=0;
  vector <string> orbs;
  vector <string> mosec;
  if(readsection(words,pos=0, mosec,"ORBITALS")) { 
    //error("Need ORBITALS section in ekt");
    complex_orbitals=false;
    allocate(mosec,sys,momat);
    nmo=momat->getNmo();
  }
  else if(readsection(words,pos=0,mosec,"CORBITALS")) { 
    complex_orbitals=true;
    allocate(mosec,sys,cmomat);
    nmo=cmomat->getNmo();
  }
  else { error("Need ORBITALS or CORBITALS in EKT"); } 
  occupations.Resize(1);

  if(readsection(words,pos=0,orbs,"STATES")) { 
    nmo=orbs.size();
    occupations[0].Resize(nmo);
    for(int i=0; i< nmo; i++) 
      occupations[0][i]=atoi(orbs[i].c_str())-1;
  }
  else if(readsection(words,pos=0,orbs,"EVAL_ORBS")) { 
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

  deterministic_psp=0;
  if(haskeyword(words,pos=0,"DETERMINISTIC_PSP"))
    deterministic_psp=1;

  dump_data=false;
  if(haskeyword(words,pos=0,"DUMP_DATA"))
    dump_data=true;
  if(!readvalue(words, pos=0,ekt_cutoff,"EKT_CUTOFF"))
    ekt_cutoff=1e3;
  
  
  if(complex_orbitals) 
    cmomat->buildLists(occupations);
  else 
    momat->buildLists(occupations);



  if(!readvalue(words, pos=0,nstep_sample,"NSTEP_SAMPLE"))
    nstep_sample=10;
  if(!readvalue(words,pos=0,npoints_eval,"NPOINTS"))
    npoints_eval=4;
  string mode;
  eval_obdm=false; 
  eval_conduction=false; 
  eval_valence=false; 
  if(readvalue(words,pos=0,mode,"MODE")) {
    if(caseless_eq(mode,"ALL") or caseless_eq(mode, "CONDUCTION")) {
      eval_obdm=true;
      eval_conduction=true;
      eval_valence=true; 
    }
    else if (caseless_eq(mode, "VALENCE")) {
      eval_obdm=true; 
      eval_conduction=true;
    }
    else if (caseless_eq(mode, "OBDM")) {
      eval_obdm=true; 
    }

  }

  eval_obdm=true; 
  eval_conduction=true; 
  eval_valence=true; 
  int ndim=3;
  saved_r.Resize(npoints_eval*2); //r1 and r2;
  for(int i=0; i< npoints_eval*2; i++) saved_r(i).Resize(ndim);
  rk.Resize(npoints_eval*2);

  int warmup_steps=1000;
  Sample_point * sample=NULL;
  sys->generateSample(sample);
  sample->randomGuess();
  Array2 <dcomplex> movals(nmo,1);
  for(int i=0; i< npoints_eval*2; i++) { 
    gen_sample(warmup_steps,1.0,0,movals,sample);
    sample->getElectronPos(0,saved_r(i));
  }
  delete sample;
}

//----------------------------------------------------------------------

void Average_ekt::write_init(string & indent, ostream & os) { 
  os << indent << "EKT" << endl;
  os << indent << "NMO " << nmo << endl;
  os << indent << "NPOINTS " << npoints_eval << endl;
  if(eval_valence) { 
    os << indent << "valence band" << endl;
  }
  if(eval_conduction) { 
    os << indent << "conduction band" << endl;
  }
  if (eval_obdm) {
    os << indent << "obdm only" << endl; 
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
  os << indent << "STATES { ";
  int nstates=occupations(0).GetDim(1);
  for(int i=0; i < nstates; i++) { 
    os << occupations(0)(i)+1 << " ";
    if( (i+1)%30==0) os << "\n" << indent << " ";
  }
  os << " } ";
}
//----------------------------------------------------------------------
void Average_ekt::read(vector <string> & words) { 
  unsigned int pos=0;

  readvalue(words, pos=0,nmo,"NMO");
  readvalue(words,pos=0,npoints_eval,"NPOINTS");
  if(haskeyword(words, pos=0,"OBDM")) 
    eval_obdm=true; 
  else 
    eval_obdm=false;
  if(haskeyword(words, pos=0,"CONDUCTION")) 
    eval_conduction=true; 
  else 
    eval_conduction=false;
  if(haskeyword(words,pos=0,"VALENCE")) 
    eval_valence=true;
  else 
    eval_valence=false;
  if(haskeyword(words,pos=0,"CORBITALS")) 
    complex_orbitals=true;
  else 
    complex_orbitals=false;

  eval_obdm=true; 
  eval_valence=true; 
  eval_conduction=true; 
  occupations.Resize(1);
  occupations(0).Resize(nmo);
  vector <string> orbs;
  if(!readsection(words,pos=0,orbs,"STATES")) { 
    cout << "WARNING: init section does not contain states. Numbering of the density matrix may be incorrect." << endl;
    for(int i=0; i < nmo; i++) 
      occupations[0][i]=i;
  }
  else { 
    nmo=orbs.size();
    occupations[0].Resize(nmo);
    for(int i=0; i< nmo; i++) 
      occupations[0][i]=atoi(orbs[i].c_str())-1;
  }

}
//----------------------------------------------------------------------
void Average_ekt::write_summary(Average_return &avg,Average_return &err, ostream & os) { 

  Array2 <dcomplex> obdm_up(nmo,nmo),obdm_down(nmo,nmo), 
         obdm_up_err(nmo,nmo),obdm_down_err(nmo,nmo), 
         cndc_up(nmo, nmo), cndc_down(nmo, nmo),
         vlnc_up(nmo, nmo), vlnc_down(nmo, nmo), 
         cndc_up_err(nmo, nmo), cndc_down_err(nmo, nmo), 
         vlnc_up_err(nmo, nmo), vlnc_down_err(nmo, nmo); 
  //Array4 <dcomplex>
  //       tbdm_uu(nmo,nmo),tbdm_uu_err(nmo,nmo),
  //      tbdm_ud(nmo,nmo),tbdm_ud_err(nmo,nmo),
  //      tbdm_du(nmo,nmo),tbdm_du_err(nmo,nmo),
  //      tbdm_dd(nmo,nmo),tbdm_dd_err(nmo,nmo);

  os << "EKT: nmo " << nmo << endl;
  os << "EKT: states { ";
  for(int i=0; i< nmo; i++) { 
    os << occupations[0][i]+1 << " " ;
  }
  os << " } \n";
  os << "Orbital normalization " << endl;
  for(int i=0; i< nmo; i++) { 
    os << avg.vals(i) << " +/- " << err.vals(i) << endl;
  }

  // for(int i=0; i < nmo; i++) { 
  //   avg.vals(i)=1.0/nmo; //assume that the orbitals are normalized.
  // }

  if (eval_obdm) {
    dcomplex i_c(0.,1.);
    int place=0;
    for(int i=0; i < nmo; i++) { 
      for(int j=0; j < nmo; j++)  { 
        doublevar norm=sqrt(avg.vals(i)*avg.vals(j));
        int index=nmo + place;
        obdm_up(i,j)=dcomplex(avg.vals(index),avg.vals(index+1))/norm;
        obdm_up_err(i,j)=dcomplex(err.vals(index),err.vals(index+1))/norm;

        index+=2*nmo*nmo;
        obdm_down(i,j)=dcomplex(avg.vals(index),avg.vals(index+1))/norm;
        obdm_down_err(i,j)=dcomplex(err.vals(index),err.vals(index+1))/norm;
        place+=2;
      }
    }

    os << "One-body density matrix " << endl;
    int colwidth=40;
    if(!complex_orbitals) colwidth=20;
    os << setw(10) << " " << setw(10) << " "
      << setw(colwidth) << "up" << setw(colwidth) << "up err"
      << setw(colwidth) << "down" << setw(colwidth) << "down err"
      << endl;
    for(int i=0; i< nmo ; i++) { 
      for(int j=0; j<nmo; j++) { 
        if(complex_orbitals) { 
          os << setw(10) << i << setw(10) << j
            << setw(colwidth) << obdm_up(i,j) << setw(colwidth) << obdm_up_err(i,j)
            << setw(colwidth) << obdm_down(i,j) << setw(colwidth) << obdm_down_err(i,j)
            << endl;
        }
        else { 
          os << setw(10) << i << setw(10) << j
            << setw(colwidth) << obdm_up(i,j).real() 
            << setw(colwidth) << obdm_up_err(i,j).real()
            << setw(colwidth) << obdm_down(i,j).real() 
            << setw(colwidth) << obdm_down_err(i,j).real()
            << endl;

        }
      }

    }
  }

  if (eval_valence) {
    dcomplex i_c(0.,1.);
    int place=0;
    for(int i=0; i < nmo; i++) { 
      for(int j=0; j < nmo; j++)  { 
        doublevar norm=sqrt(avg.vals(i)*avg.vals(j));
        int index=nmo+4*nmo*nmo + place;
        vlnc_up(i,j)=dcomplex(avg.vals(index),avg.vals(index+1))/norm;
        vlnc_up_err(i,j)=dcomplex(err.vals(index),err.vals(index+1))/norm;

        index+=2*nmo*nmo;
        vlnc_down(i,j)=dcomplex(avg.vals(index),avg.vals(index+1))/norm;
        vlnc_down_err(i,j)=dcomplex(err.vals(index),err.vals(index+1))/norm;
        place+=2;
      }
    }

    os << "Valence band matrix " << endl;
    int colwidth=40;
    if(!complex_orbitals) colwidth=20;
    os << setw(10) << " " << setw(10) << " "
      << setw(colwidth) << "up" << setw(colwidth) << "up err"
      << setw(colwidth) << "down" << setw(colwidth) << "down err"
      << endl;
    for(int i=0; i< nmo ; i++) { 
      for(int j=0; j<nmo; j++) { 
        if(complex_orbitals) { 
          os << setw(10) << i << setw(10) << j
            << setw(colwidth) << vlnc_up(i,j) << setw(colwidth) << vlnc_up_err(i,j)
            << setw(colwidth) << vlnc_down(i,j) << setw(colwidth) << vlnc_down_err(i,j)
            << endl;
        }
        else { 
          os << setw(10) << i << setw(10) << j
            << setw(colwidth) << vlnc_up(i,j).real() 
            << setw(colwidth) << vlnc_up_err(i,j).real()
            << setw(colwidth) << vlnc_down(i,j).real() 
            << setw(colwidth) << vlnc_down_err(i,j).real()
            << endl;

        }
      }

    }
  }

  if (eval_conduction) {
    dcomplex i_c(0.,1.);
    int place=0;
    for(int i=0; i < nmo; i++) { 
      for(int j=0; j < nmo; j++)  { 
        doublevar norm=sqrt(avg.vals(i)*avg.vals(j));
        int index=nmo + 8*nmo*nmo + place;
        cndc_up(i,j)=dcomplex(avg.vals(index),avg.vals(index+1))/norm;
        cndc_up_err(i,j)=dcomplex(err.vals(index),err.vals(index+1))/norm;

        index+=2*nmo*nmo;
        cndc_down(i,j)=dcomplex(avg.vals(index),avg.vals(index+1))/norm;
        cndc_down_err(i,j)=dcomplex(err.vals(index),err.vals(index+1))/norm;
        place+=2;
      }
    }

    os << "Conduction band matrix " << endl;
    int colwidth=40;
    if(!complex_orbitals) colwidth=20;
    os << setw(10) << " " << setw(10) << " "
      << setw(colwidth) << "up" << setw(colwidth) << "up err"
      << setw(colwidth) << "down" << setw(colwidth) << "down err"
      << endl;
    for(int i=0; i< nmo ; i++) { 
      for(int j=0; j<nmo; j++) { 
        if(complex_orbitals) { 
          os << setw(10) << i << setw(10) << j
            << setw(colwidth) << cndc_up(i,j) << setw(colwidth) << cndc_up_err(i,j)
            << setw(colwidth) << cndc_down(i,j) << setw(colwidth) << cndc_down_err(i,j)
            << endl;
        }
        else { 
          os << setw(10) << i << setw(10) << j
            << setw(colwidth) << cndc_up(i,j).real() 
            << setw(colwidth) << cndc_up_err(i,j).real()
            << setw(colwidth) << cndc_down(i,j).real() 
            << setw(colwidth) << cndc_down_err(i,j).real()
            << endl;

        }
      }

    }
  }
  /*
     totnelectrons = 2; 
     cout << "Ntot: " << totnelectrons << endl; 
     for(int i=0; i < totnelectrons; i++) { 
     cout << "Kin: " << i << " " << avg.vals(nmo + 12*nmo*nmo + i) << "+/-" << err.vals(nmo + 12*nmo*nmo + i) << endl; 
     cout << "Pot: " << i << " " << avg.vals(nmo + 12*nmo*nmo + totnelectrons + i) << "+/-" << err.vals(nmo + 12*nmo*nmo + totnelectrons + i) << endl; 
     }*/


  dcomplex trace=0;
  for(int i=0;i< nmo; i++) trace+=obdm_up(i,i);
  os << "Trace of the obdm: up: " << trace;
  trace=0;
  for(int i=0; i< nmo; i++) trace+=obdm_down(i,i);
  os << " down: " << trace << endl;
}


//------------------------------------------------------------------
void Average_ekt::joutput(Array2 <dcomplex> &obdm, ostream & os) {
  
  for(int i=0; i< nmo ; i++) {
    os << "[";
    for(int j=0; j<nmo; j++) {
      if(complex_orbitals) {
        if(j==nmo-1) os << "["<< obdm(i,j).real() <<","<< obdm(i,j).imag()<<"]";
        else  os << "["<< obdm(i,j).real() <<","<< obdm(i,j).imag()<<"],";
      }
      else {
        if(j==nmo-1) os << obdm(i,j).real();
        else os << obdm(i,j).real() << ",";
      }
    }
    if(i==nmo-1) os << "]" << endl;
    else os <<"]," << endl;
  }
}


//-----------------------------------------------------------------
void Average_ekt::jsonOutput(Average_return &avg,Average_return &err, ostream & os) {
  
  Array2 <dcomplex> obdm_up(nmo,nmo),obdm_down(nmo,nmo),
  obdm_up_err(nmo,nmo),obdm_down_err(nmo,nmo),
  cndc_up(nmo, nmo), cndc_down(nmo, nmo),
  vlnc_up(nmo, nmo), vlnc_down(nmo, nmo),
  cndc_up_err(nmo, nmo), cndc_down_err(nmo, nmo),
  vlnc_up_err(nmo, nmo), vlnc_down_err(nmo, nmo);
  
  os <<"\""<< avg.type << "\":{" << endl;
  os << "\"nmo\":" << nmo <<","<< endl;
  os << "\"states\":[";
  for(int i=0; i< nmo; i++) {
    if(i==nmo-1) os << occupations[0][i]+1 ;
    else os << occupations[0][i]+1 << "," ;
  }
  os << "]," << endl;
  os << "\"normalization\":{" << endl;
  os << "\"value\":[";
  for(int i=0; i< nmo; i++) {
    if(i< nmo-1) os << avg.vals(i) << ",";
    else os << avg.vals(i);
  }
  os << "],"<< endl;
  os << "\"error\":[";
  for(int i=0; i< nmo; i++) {
    if(i< nmo-1) os << err.vals(i) << ",";
    else os << err.vals(i);
  }
  os << "]"<< endl;
  os << "}," << endl;
  
  if(eval_obdm) {
    dcomplex i_c(0.,1.);
    int place=0;
    for(int i=0; i < nmo; i++) {
      for(int j=0; j < nmo; j++)  {
        doublevar norm=sqrt(avg.vals(i)*avg.vals(j));
        int index=nmo + place;
        obdm_up(i,j)=dcomplex(avg.vals(index),avg.vals(index+1))/norm;
        obdm_up_err(i,j)=dcomplex(err.vals(index),err.vals(index+1))/norm;
        
        index+=2*nmo*nmo;
        obdm_down(i,j)=dcomplex(avg.vals(index),avg.vals(index+1))/norm;
        obdm_down_err(i,j)=dcomplex(err.vals(index),err.vals(index+1))/norm;
        place+=2;
      }
    }
    os << "\"obdm\":{" << endl;
    int colwidth=40;
    if(!complex_orbitals) colwidth=20;
    
    os << "\"up\":["<< endl;
    joutput(obdm_up, os);
    os << "]," << endl;
    os << "\"up_err\":["<< endl;
    joutput(obdm_up_err, os);
    os << "]," << endl;
    os << "\"down\":["<< endl;
    joutput(obdm_down, os);
    os << "]," << endl;
    os << "\"down_err\":["<< endl;
    joutput(obdm_down_err, os);
    os << "]" << endl;
    os << "}," << endl;
  }
  
  if (eval_valence) {
    dcomplex i_c(0.,1.);
    int place=0;
    for(int i=0; i < nmo; i++) {
      for(int j=0; j < nmo; j++)  {
        doublevar norm=sqrt(avg.vals(i)*avg.vals(j));
        int index=nmo+4*nmo*nmo + place;
        vlnc_up(i,j)=dcomplex(avg.vals(index),avg.vals(index+1))/norm;
        vlnc_up_err(i,j)=dcomplex(err.vals(index),err.vals(index+1))/norm;
        
        index+=2*nmo*nmo;
        vlnc_down(i,j)=dcomplex(avg.vals(index),avg.vals(index+1))/norm;
        vlnc_down_err(i,j)=dcomplex(err.vals(index),err.vals(index+1))/norm;
        place+=2;
      }
    }
    
    os << "\"valence\":{" << endl;
    int colwidth=40;
    if(!complex_orbitals) colwidth=20;
    
    os << "\"up\":["<< endl;
    joutput(vlnc_up, os);
    os << "]," << endl;
    os << "\"up_err\":["<< endl;
    joutput(vlnc_up_err, os);
    os << "]," << endl;
    os << "\"down\":["<< endl;
    joutput(vlnc_down, os);
    os << "]," << endl;
    os << "\"down_err\":["<< endl;
    joutput(vlnc_down_err, os);
    os << "]" << endl;
    os << "}," << endl;
  }
  
  if (eval_conduction) {
    dcomplex i_c(0.,1.);
    int place=0;
    for(int i=0; i < nmo; i++) {
      for(int j=0; j < nmo; j++)  {
        doublevar norm=sqrt(avg.vals(i)*avg.vals(j));
        int index=nmo + 8*nmo*nmo + place;
        cndc_up(i,j)=dcomplex(avg.vals(index),avg.vals(index+1))/norm;
        cndc_up_err(i,j)=dcomplex(err.vals(index),err.vals(index+1))/norm;
        
        index+=2*nmo*nmo;
        cndc_down(i,j)=dcomplex(avg.vals(index),avg.vals(index+1))/norm;
        cndc_down_err(i,j)=dcomplex(err.vals(index),err.vals(index+1))/norm;
        place+=2;
      }
    }
    
    os << "\"conduction\":{" << endl;
    int colwidth=40;
    if(!complex_orbitals) colwidth=20;
    
    os << "\"up\":["<< endl;
    joutput(cndc_up, os);
    os << "]," << endl;
    os << "\"up_err\":["<< endl;
    joutput(cndc_up_err, os);
    os << "]," << endl;
    os << "\"down\":["<< endl;
    joutput(cndc_down, os);
    os << "]," << endl;
    os << "\"down_err\":["<< endl;
    joutput(cndc_down_err, os);
    os << "]" << endl;
    os << "}" << endl;
  }
  
  os << "}" << endl;
}

//----------------------------------------------------------------------
//
//
//Note this needs to be changed for non-zero k-points!
doublevar Average_ekt::gen_sample(int nstep, doublevar  tstep, 
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
/*
   Return the molecule orbital values
   */
void Average_ekt::calc_mos(Sample_point * sample, int e, Array2 <dcomplex> & movals) { 
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


/*
   Calculate the values, gradient and laplacians of molecule orbitals, 
   returned as: [val, grad, lap] 
   */
void Average_ekt::calc_mosLap(Sample_point * sample, int e, Array2 <dcomplex> & molaps) { 
  if(complex_orbitals) { 
    molaps.Resize(nmo,5);
    cmomat->updateLap(sample, e, 0, molaps);
  }
  else { 
    Array2 <doublevar> molaps_d(nmo,5);
    momat->updateLap(sample,e,0,molaps_d);
    molaps.Resize(nmo,5);
    for(int i=0; i< nmo; i++) { 
      for (int d=0; d<5; d++) 
        molaps(i,d)=molaps_d(i,d);
    }
  }
}
//----------------------------------------------------------------------

doublevar legendre(doublevar x, int n); 

void Average_ekt::calcPseudoMo(System * sys,
    Sample_point * sample,
    Pseudopotential *psp, 
    const Array1 <doublevar> & accept_var,
    Array1 <dcomplex> & totalv)//, 
{
  int natoms=sample->ionSize();
  assert(accept_var.GetDim(0) >= nTest());
  //assert(totalv.GetDim(0) >= nwf);
  assert(nelectrons == sample->electronSize());

  Array1 <doublevar> ionpos(3), oldpos(3), newpos(3);
  Array1 <doublevar> newdist(5), olddist(5);

  int accept_counter=0;
  //deriv.Resize(natoms, 3);
  //deriv=0;
  dcomplex nonlocal; 
  Array1 <doublevar> rDotR(psp->getMaxAIP());
  Array2 <dcomplex> movals_t(nmo, 1);
  Array2 <dcomplex> movals(nmo, 1);
  sample->getElectronPos(0, oldpos);
  calc_mos(sample, 0, movals);
  dcomplex local; 
  for (int jmo = 0; jmo < nmo; jmo ++) {
    nonlocal=0.0; 
    local = 0.0; 
    for(int at = 0; at< natoms; at++){
      if(psp->getNumL(at) != 0) {
        Array1 <doublevar> v_l(psp->getNumL(at));
        sample->getIonPos(at, ionpos);
        sample->updateEIDist();
        sample->getEIDist(0, at, olddist);
        int spin=1;
        //	if(e < sys->nelectrons(0)) spin=0;
        spin = 0;
        psp->getRadialOut(at, spin, sample, olddist, v_l);

        //----------------------------------------
        //Start integral

        int accept;
        if(psp->getDeterministic()) {
          accept= olddist(0) < psp->getCutoff(at);
        }
        else {
          doublevar strength=0;
          const doublevar calculate_threshold=10;

          for(int l=0; l<psp->getNumL(at)-1; l++) {
            strength+=calculate_threshold*(2*l+1)*fabs(v_l(l));
          }
          strength=min((doublevar) 1.0, strength);
          doublevar rand=accept_var(accept_counter++);
          //cout << at <<"  random number  " << rand
          //   << "  p_eval  " << strength  << endl;
          if ( strength > 0.0 ) {
            for(int l=0; l<psp->getNumL(at)-1; l++)
              v_l(l)/=strength;
          }
          accept=strength>rand;
        }

        //bool localonly = true;
        if(accept)  {
          for(int i=0; i< psp->getAIP(at); i++) {
            for(int d=0; d < 3; d++) 
              newpos(d)=psp->getIntegralPt(at,i,d)*olddist(0)-olddist(d+2);
            sample->translateElectron(0, newpos);
            sample->updateEIDist();
            sample->getEIDist(0,at,newdist);
            rDotR(i)=0;
            for(int d=0; d < 3; d++)
              rDotR(i)+=newdist(d+2)*olddist(d+2);

            rDotR(i)/=(newdist(0)*olddist(0));  //divide by the magnitudes
            calc_mos(sample, 0, movals_t);//update value
            //----
            //cout << "signs " << base_sign << "  " << new_sign << endl;;
            for(int l=0; l< psp->getNumL(at)-1; l++) {
              nonlocal+=(2*l+1)*v_l(l)*legendre(rDotR(i), l)*psp->getIntegralWeight(at, i)*movals_t(jmo, 0);
            }
            sample->setElectronPos(0, oldpos);
          } 
        }
        int localL=psp->getNumL(at)-1; //The l-value of the local part is
        local += v_l(localL); 
      } // if atom has psp
    }  //atom loop
    //    cout << "Local:" << local << endl; 
    //cout << "Nonlocal: " << nonlocal/movals(jmo, 0) << endl; 
    totalv(jmo) = local*movals(jmo, 0)+nonlocal;
    //totalv(jmo) =0.0; 
  }  //nmo loop
  //cout << "psp: local part " << accum_local
  // << "  nonlocal part " << accum_nonlocal << endl;
}
