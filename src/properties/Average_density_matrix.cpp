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

void Average_tbdm_basis::randomize(Wavefunction_data * wfdata, Wavefunction * wf,
                       System * sys, Sample_point * sample_tmp) { 
  int nup=sys->nelectrons(0);
  int ndown=sys->nelectrons(1);
  Array2 <dcomplex> movals1(nmo,1),movals2(nmo,1);
  //Make a copy of the Sample so that it doesn't have to update the wave function.
  Sample_point * sample=NULL;
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

void Average_tbdm_basis::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & avg) { 
  if(eval_old) evaluate_old(wfdata,wf,sys,sample,avg);
  else { 
    if(eval_tbdm) evaluate_tbdm(wfdata,wf,sys,sample,avg);
    else evaluate_obdm(wfdata,wf,sys,sample,avg);
  }
}


//----------------------------------------------------------------------

void Average_tbdm_basis::evaluate_obdm(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & avg) { 
  avg.type="tbdm_basis";
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
    /*! 
      Here, permutation symmetry has been applied, the real evaluation quantity is 
      sum_(e=1)^N phi_i* (r_e) phi_j(r') psi(r_1,... r',...r_N)/psi(r_1, ..., r_n, ..., r_N)
      Therefore, movals1_base[e, i] = phi_i(r_e), where e is the index of the electron. 
      !!!One need to be careful about the off-gamma point 1 RDM 
     */
  }
  avg.vals.Resize(nmo+4*nmo*nmo);
  avg.vals=0;

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
      // This is because we are doing the one body reduced density matrix
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
      /*
	shouldn't the normalization factor just /norm(movals1(orb, 0))
       */
      for(int orbnum=0; orbnum < nmo; orbnum++) { 
        for(int orbnum2=0; orbnum2 < nmo; orbnum2++) { 
          //tmp=movals1(orbnum,0)*conj(movals1_base(e)(orbnum2,0))*prefactor;
	  /*!
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

void Average_tbdm_basis::evaluate_old(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & avg) { 

  avg.type="tbdm_basis";
  
  //orbital normalization, then 1bdm up, 1bdm down, 2bdm up up, 2bdm up down 2bdm down up 2bdm down down, with everything complex
  if(eval_tbdm) { 
    avg.vals.Resize(nmo+4*nmo*nmo+8*nmo*nmo*nmo*nmo);
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
    int k=rk(i);
    int l=rk(npoints_eval+i);
    //Calculate the orbital values for r1 and r2
    calc_mos(sample,k,movals1_old);
    calc_mos(sample,l,movals2_old);
    sample->getElectronPos(k,oldr1);
    sample->getElectronPos(l,oldr2);

    r1=saved_r(i);
    r2=saved_r(npoints_eval+i);
    
    sample->setElectronPos(k,r1);
    //doublevar dist1=gen_sample(nstep_sample,1.0,k,movals1, sample);
    doublevar dist1=0,dist2=0;

    calc_mos(sample,k,movals1);
    for(int m=0; m < nmo; m++) 
      dist1+=norm(movals1(m,0));
    

    wf->updateVal(wfdata,sample);
    wf->getVal(wfdata,0,wfval_1b);


    dcomplex psiratio_1b=exp(dcomplex(wfval_1b.amp(0,0)-wfval_base.amp(0,0),
          wfval_1b.phase(0,0)-wfval_base.phase(0,0)));
   
     dcomplex psiratio_2b;
    if(eval_tbdm) { 

      sample->setElectronPos(l,r2); 
      //dist2=gen_sample(nstep_sample,1.0,l,movals2,sample);
      calc_mos(sample,l,movals2);
      for(int m=0; m < nmo; m++) 
        dist2+=norm(movals2(m,0));
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
    tbdm_t which_tbdm;
    int which_obdm=0;
    int npairs=0;
    int nelec_1b=nup;
    if(k >= nup) { which_obdm=1; nelec_1b=ndown; } 
    if(k < nup and l < nup) { which_tbdm=tbdm_uu; npairs=nup*(nup-1);nupup++; } 
    else if(k < nup and l >= nup) { which_tbdm=tbdm_ud; npairs=nup*ndown; nupdown++;} 
    else if(k >= nup and l < nup) { which_tbdm=tbdm_du; npairs=ndown*nup; ndownup++; } 
    else if(k >= nup and l >= nup) { which_tbdm=tbdm_dd; npairs=ndown*(ndown-1);ndowndown++; } 

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
      for(int oi=0; oi < nmo; oi++) { 
        for(int oj=0; oj < nmo; oj++) { 
          for(int ok=0; ok < nmo; ok++) { 
            for(int ol=0; ol < nmo; ol++) { 
              tmp=doublevar(npairs)*conj(movals1(oi,0))*movals1_old(ok,0)
                *conj(movals2(oj,0))*movals2_old(ol,0)
                *psiratio_2b/dist1/dist2;
              int ind=tbdm_index(which_tbdm,oi,oj,ok,ol);
              avg.vals(ind)+=tmp.real();
              avg.vals(ind+1)+=tmp.imag();
              //cout << "index " << ind << endl;
            }
          }
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
        place++;
      }
    }
  }
  

  if(eval_tbdm) { 
    place=0;
    for(int oi=0; oi < nmo; oi++) { 
      for(int oj=0; oj < nmo; oj++) { 
        for (int ok=0; ok < nmo; ok++) { 
          for(int ol=0; ol < nmo; ol++) { 
            int induu=tbdm_index(tbdm_uu,oi,oj,ok,ol);
            int indud=tbdm_index(tbdm_ud,oi,oj,ok,ol);
            int inddu=tbdm_index(tbdm_du,oi,oj,ok,ol);
            int inddd=tbdm_index(tbdm_dd,oi,oj,ok,ol);
            for(int k=0; k < 2; k++) { 
              avg.vals(induu+k)/=nupup;
              avg.vals(indud+k)/=nupdown;
              avg.vals(inddu+k)/=ndownup;
              avg.vals(inddd+k)/=ndowndown;
            }
          }
        }
      }
    }
  }


  //cout << nupup << " " << nupdown << " " << ndownup << " " << ndowndown << endl;


}

void Average_tbdm_basis::evaluate_tbdm(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & avg) { 
  avg.type="tbdm_basis";

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
  }
  
  //orbital normalization, then 1bdm up, 1bdm down, 2bdm up up, 2bdm up down 2bdm down up 2bdm down down, with everything complex
  if(tbdm_diagonal)
    avg.vals.Resize(nmo+16*nmo*nmo);
  else 
    avg.vals.Resize(nmo+4*nmo*nmo+8*nmo*nmo*nmo*nmo);


  avg.vals=0;

  Wavefunction_storage * store;
  wf->generateStorage(store);
 

  Wf_return wfval1(wf->nfunc(),2);
  Array1 <Wf_return> wfval2(nelectrons);
  Array2 <dcomplex> movals1(nmo,1),movals2(nmo,1);
  Array1 <doublevar> oldpos1(3),oldpos2(3);
  for(int i=0; i< npoints_eval; i++) { 
    for(int e1=0; e1 < nelectrons; e1++) { 


      //Evaluate all the quantities we will be needing here
      sample->getElectronPos(e1,oldpos1);
      wf->saveUpdate(sample,e1,store);
      sample->setElectronPos(e1,saved_r(i));
      calc_mos(sample,e1,movals1);
      wf->updateVal(wfdata,sample);
      wf->getVal(wfdata,e1,wfval1);
      wf->evalTestPos(saved_r(i+npoints_eval),sample,wfval2);

      sample->getElectronPos(0,oldpos2);
      sample->setElectronPosNoNotify(0,saved_r(i+npoints_eval));
      calc_mos(sample,0,movals2);
      sample->setElectronPosNoNotify(0,oldpos2);
     
      sample->setElectronPos(e1,oldpos1);
      wf->restoreUpdate(sample,e1,store);
      //Done evaluating

      //---Add to the OBDM
     
      doublevar dist1=0,dist2=0;
      for(int m=0; m < nmo; m++)  { 
        dist1+=norm(movals1(m,0));
        dist2+=norm(movals2(m,0));
      }
      for(int orbnum=0; orbnum < nmo; orbnum++) { 
        avg.vals(orbnum)+=norm(movals1(orbnum,0))/(dist1*npoints_eval*nelectrons);
      }

      dcomplex psiratio_1b=exp(dcomplex(wfval1.amp(0,0)-wfval_base.amp(0,0),
            wfval1.phase(0,0)-wfval_base.phase(0,0)));

      int which_obdm=0;
      if(e1 >= nup) { which_obdm=1;  } 
      dcomplex tmp;
      int place=0;
      for(int orbnum=0; orbnum < nmo; orbnum++) { 
        for(int orbnum2=0; orbnum2 < nmo; orbnum2++) { 
          tmp=conj(movals1(orbnum,0))*movals1_base(e1)(orbnum2,0)
            *psiratio_1b/dist1;
	  /*Corrected*/
          avg.vals(nmo+2*which_obdm*nmo*nmo+place)+=tmp.real()/doublevar(npoints_eval);
          avg.vals(nmo+2*which_obdm*nmo*nmo+place+1)+=tmp.imag()/doublevar(npoints_eval);

          place+=2;
        }
      }
      //---Add to the TBDM
      
      for(int e2=0; e2 < nelectrons; e2++)  { 
        
        if(e2!=e1)  { 
          dcomplex psiratio_2b=exp(dcomplex(wfval2(e2).amp(0,0)-wfval_base.amp(0,0),
                wfval2(e2).phase(0,0)-wfval_base.phase(0,0)));
          tbdm_t which_tbdm;
          if(e1 < nup and e2 < nup) { which_tbdm=tbdm_uu; } 
          else if(e1 < nup and e2 >= nup) { which_tbdm=tbdm_ud; } 
          else if(e1 >= nup and e2 < nup) { which_tbdm=tbdm_du;  } 
          else if(e1 >= nup and e2 >= nup) { which_tbdm=tbdm_dd; } 
          dcomplex orbind=psiratio_2b/(dist1*dist2*npoints_eval); 
          dcomplex tmp1,tmp2,tmp3; 
          if(tbdm_diagonal) { 
            int ind=nmo+4*nmo*nmo+2*which_tbdm*nmo*nmo;
            for(int oi=0; oi < nmo; oi++) { 
              tmp1=orbind*conj(movals1(oi,0));
              for(int oj=0; oj < nmo; oj++) { 
                int ok=oi; int ol=oj;
                tmp2=tmp1*conj(movals2(oj,0));
                tmp3=tmp2*movals1_base(e1)(ok,0);
                tmp=tmp3*movals1_base(e2)(ol,0);
                avg.vals.v[ind++]+=tmp.real();
                avg.vals.v[ind++]+=tmp.imag();
              }
            }
          } //TBDM_DIAGONAL
          else { 
            int ind=tbdm_index(which_tbdm,0,0,0,0);
            for(int oi=0; oi < nmo; oi++) { 
              tmp1=orbind*conj(movals1(oi,0));
              for(int oj=0; oj < nmo; oj++) { 
                tmp2=tmp1*conj(movals2(oj,0));

                for(int ok=0; ok < nmo; ok++) {
                  tmp3=tmp2*movals1_base(e1)(ok,0);
                  for(int ol=0; ol< nmo; ol++) { 
                    tmp=tmp3*movals1_base(e2)(ol,0);
                    avg.vals.v[ind++]+=tmp.real();
                    avg.vals.v[ind++]+=tmp.imag();
                  }
                }
              }
            }
          }//offdiagonal tbdm
        }
      }
    }
  }

  delete store;

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



  if(complex_orbitals) 
    cmomat->buildLists(occupations);
  else 
    momat->buildLists(occupations);



  if(!readvalue(words, pos=0,nstep_sample,"NSTEP_SAMPLE"))
    nstep_sample=10;
  if(!readvalue(words,pos=0,npoints_eval,"NPOINTS"))
    npoints_eval=4;


  eval_old=false;
  if(haskeyword(words,pos=0,"NUMERICAL_EVALUATION")) eval_old=true;

  string mode;
  eval_tbdm=true;
  tbdm_diagonal=false;
  if(readvalue(words,pos=0,mode,"MODE")) {
    if(caseless_eq(mode,"OBDM")) eval_tbdm=false;
    if(caseless_eq(mode,"TBDM_DIAGONAL")) tbdm_diagonal=true;
  }
  
  
  if(haskeyword(words,pos=0,"ONLY_OBDM")) 
    error("ONLY_OBDM is depreciated; use MODE instead.");
    //eval_tbdm=false;
  
  //if(haskeyword(words,pos=0,"TBDM_DIAGONAL")) 
  //  error("TBDM_DIAGONAL is depreciated; use MODE instead.");

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

void Average_tbdm_basis::write_init(string & indent, ostream & os) { 
  os << indent << "TBDM_BASIS" << endl;
  os << indent << "NMO " << nmo << endl;
  os << indent << "NPOINTS " << npoints_eval << endl;
  if(!eval_tbdm) { 
    os << indent << "ONLY_OBDM" << endl;
  }
  if(tbdm_diagonal) { 
    os << indent << "TBDM_DIAGONAL" << endl;
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
  int nstates=occupations(0).GetDim(0);
  for(int i=0; i < nstates; i++) { 
    os << occupations(0)(i)+1 << " ";
    if( (i+1)%30==0) os << "\n" << indent << " ";
  }
  os << " } ";
}
//----------------------------------------------------------------------
void Average_tbdm_basis::read(vector <string> & words) { 
  unsigned int pos=0;

  readvalue(words, pos=0,nmo,"NMO");
  readvalue(words,pos=0,npoints_eval,"NPOINTS");
  if(haskeyword(words, pos=0,"ONLY_OBDM")) eval_tbdm=false;
  else eval_tbdm=true;
  if(haskeyword(words,pos=0,"TBDM_DIAGONAL")) tbdm_diagonal=true;
  else tbdm_diagonal=false;
  if(haskeyword(words,pos=0,"CORBITALS")) complex_orbitals=true;
  else complex_orbitals=false;

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
void Average_tbdm_basis::write_summary(Average_return &avg,Average_return &err, ostream & os) { 

  Array2 <dcomplex> obdm_up(nmo,nmo),obdm_down(nmo,nmo), 
         obdm_up_err(nmo,nmo),obdm_down_err(nmo,nmo);
  //Array4 <dcomplex>
  //       tbdm_uu(nmo,nmo),tbdm_uu_err(nmo,nmo),
   //      tbdm_ud(nmo,nmo),tbdm_ud_err(nmo,nmo),
   //      tbdm_du(nmo,nmo),tbdm_du_err(nmo,nmo),
   //      tbdm_dd(nmo,nmo),tbdm_dd_err(nmo,nmo);

  os << "tbdm: nmo " << nmo << endl;
  os << "tbdm: states { ";
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
  if(eval_tbdm) { 
    os << "two-body density matrix " << endl;
    os << setw(10) << " i " << setw(10) << " j " << setw(10) << " k " 
      << " l " 
      << setw(colwidth) << "upup" << setw(colwidth) << "upup err"
      << setw(colwidth) << "updown" << setw(colwidth) << "updown err"
      << setw(colwidth) << "downup" << setw(colwidth) << "downup err"
      << setw(colwidth) << "downdown" << setw(colwidth) << "downdown err"
      << endl;
    dcomplex uu,uu_err,ud,ud_err,du,du_err,dd,dd_err;

    if(tbdm_diagonal) { 
      for(int i=0; i< nmo; i++) { 
        for(int j=0; j< nmo; j++) { 
          doublevar norm=sqrt(avg.vals(i)*avg.vals(i)*avg.vals(j)*avg.vals(j));
          int induu=nmo+4*nmo*nmo+2*i*nmo+2*j;
          int indud=induu+2*nmo*nmo;
          int inddu=indud+2*nmo*nmo;
          int inddd=inddu+2*nmo*nmo;
          uu=dcomplex(avg.vals(induu),avg.vals(induu+1))/norm;
          uu_err=dcomplex(err.vals(induu),err.vals(induu+1))/norm;
          ud=dcomplex(avg.vals(indud),avg.vals(indud+1))/norm;
          ud_err=dcomplex(err.vals(indud),err.vals(indud+1))/norm;
          du=dcomplex(avg.vals(inddu),avg.vals(inddu+1))/norm;
          du_err=dcomplex(err.vals(inddu),err.vals(inddu+1))/norm;
          dd=dcomplex(avg.vals(inddd),avg.vals(inddd+1))/norm;
          dd_err=dcomplex(err.vals(inddd),err.vals(inddd+1))/norm;


          if(complex_orbitals) { 
            os << setw(10) << i << setw(10) << j << setw(10) << i <<  setw(10) << j 
              << setw(colwidth) << uu << setw(colwidth) << uu_err
              << setw(colwidth) << ud << setw(colwidth) << ud_err
              << setw(colwidth) << du << setw(colwidth) << du_err
              << setw(colwidth) << dd << setw(colwidth) << dd_err
              << endl;
          }
          else { 
            os << setw(10) << i << setw(10) << j << setw(10) << i <<  setw(10) << j 
              << setw(colwidth) << uu.real() << setw(colwidth) << uu_err.real()
              << setw(colwidth) << ud.real() << setw(colwidth) << ud_err.real()
              << setw(colwidth) << du.real() << setw(colwidth) << du_err.real()
              << setw(colwidth) << dd.real() << setw(colwidth) << dd_err.real()
              << endl;
            
          }
        }
      }

    }
    else { 
      for(int i=0; i< nmo ; i++) { 
        for(int j=0; j<nmo; j++) { 
          for(int k=0; k< nmo; k++) {
            for(int l=0; l< nmo; l++) { 
              doublevar norm=sqrt(avg.vals(i)*avg.vals(j)*avg.vals(k)*avg.vals(l));
              int induu=tbdm_index(tbdm_uu,i,j,k,l);
              int indud=tbdm_index(tbdm_ud,i,j,k,l);
              int inddu=tbdm_index(tbdm_du,i,j,k,l);
              int inddd=tbdm_index(tbdm_dd,i,j,k,l);

              uu=dcomplex(avg.vals(induu),avg.vals(induu+1))/norm;
              uu_err=dcomplex(err.vals(induu),err.vals(induu+1))/norm;
              ud=dcomplex(avg.vals(indud),avg.vals(indud+1))/norm;
              ud_err=dcomplex(err.vals(indud),err.vals(indud+1))/norm;
              du=dcomplex(avg.vals(inddu),avg.vals(inddu+1))/norm;
              du_err=dcomplex(err.vals(inddu),err.vals(inddu+1))/norm;
              dd=dcomplex(avg.vals(inddd),avg.vals(inddd+1))/norm;
              dd_err=dcomplex(err.vals(inddd),err.vals(inddd+1))/norm;
              if(complex_orbitals) { 
                os << setw(10) << i << setw(10) << j << setw(10) << k <<  setw(10) << l 
                  << setw(colwidth) << uu << setw(colwidth) << uu_err
                  << setw(colwidth) << ud << setw(colwidth) << ud_err
                  << setw(colwidth) << du << setw(colwidth) << du_err
                  << setw(colwidth) << dd << setw(colwidth) << dd_err
                  << endl;
              }
              else { 
                os << setw(10) << i << setw(10) << j << setw(10) << k <<  setw(10) << l 
                  << setw(colwidth) << uu.real() << setw(colwidth) << uu_err.real()
                  << setw(colwidth) << ud.real() << setw(colwidth) << ud_err.real()
                  << setw(colwidth) << du.real() << setw(colwidth) << du_err.real()
                  << setw(colwidth) << dd.real() << setw(colwidth) << dd_err.real()
                  << endl;

              }



            }
          }
        }

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

// print out the diagonal elements of two body density matrix
void Average_tbdm_basis::jdgoutput(string s, int indx,  Average_return &avg1, Average_return &avg2, ostream & os){
  
  for(int i=0; i< nmo; i++) {
    for(int j=0; j< nmo; j++) {
      doublevar norm=sqrt(avg1.vals(i)*avg1.vals(i)*avg1.vals(j)*avg1.vals(j));
      vector<int> ind(4);
      ind[0]= nmo+4*nmo*nmo+2*i*nmo+2*j;
      ind[1]= ind[0]+2*nmo*nmo;
      ind[2]= ind[1]+2*nmo*nmo;
      ind[3]= ind[2]+2*nmo*nmo;
      dcomplex uu;
      uu=dcomplex(avg2.vals(ind[indx]),avg2.vals(ind[indx]+1))/norm;
      if(complex_orbitals) {
        if(i==nmo-1 && j==nmo-1) os << "[" << uu.real() <<","<< uu.imag()<<"]";
        else os <<"["<< uu.real() <<","<< uu.imag()<<"],";
      }
      else {
        if( i==nmo-1 && j==nmo-1) os << uu.real();
        else os << uu.real() << ",";
      }
    }
  }
}
//----------------------------------------------------------------------

// print out one body density matrix in json format
void Average_tbdm_basis::jobdmoutput(Array2 <dcomplex> &obdm, ostream & os) {
  
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
//----------------------------------------------------------------------

void Average_tbdm_basis::jsontbdmvalHelper(Average_return & avg,tbdm_t tbdm,ostream & os) { 
  dcomplex uu;
  for(int i=0; i< nmo ; i++) {
    os << "[";
    for(int j=0; j<nmo; j++) {
      os << "[";
      for(int k=0; k< nmo; k++) {
        os << "[";
        for(int l=0; l< nmo; l++) {
          doublevar norm=sqrt(avg.vals(i)*avg.vals(j)*avg.vals(k)*avg.vals(l));
          int ind=tbdm_index(tbdm,i,j,k,l);
          uu=dcomplex(avg.vals(ind),avg.vals(ind+1))/norm;
          if(complex_orbitals) {
            if(l==nmo-1) os << "[" << uu.real() <<","<< uu.imag()<<"]";
            else os <<"["<< uu.real() <<","<< uu.imag()<<"],";
          }
          else {
            os << uu.real();
            if(l<nmo-1) os << ",";
          }
        }
        if(k==nmo-1) os << "]" << endl;
        else os << "]," << endl;
      }
      if(j==nmo-1) os << "]" << endl;
      else os << "]," << endl;
    }
    if(i==nmo-1) os << "]" << endl;
    else os << "]," << endl;
  }
  
}

//----------------------------------------------------------------------


void Average_tbdm_basis::jsontbdmerrHelper(Average_return & avg,
                                           Average_return & err,
                                          tbdm_t tbdm,ostream & os) { 
  dcomplex uu;
  for(int i=0; i< nmo ; i++) {
    os << "[";
    for(int j=0; j<nmo; j++) {
      os << "[";
      for(int k=0; k< nmo; k++) {
        os << "[";
        for(int l=0; l< nmo; l++) {
          doublevar norm=sqrt(avg.vals(i)*avg.vals(j)*avg.vals(k)*avg.vals(l));
          int ind=tbdm_index(tbdm,i,j,k,l);
          uu=dcomplex(err.vals(ind),err.vals(ind+1))/norm;
          if(complex_orbitals) {
            if(l==nmo-1) os << "[" << uu.real() <<","<< uu.imag()<<"]";
            else os <<"["<< uu.real() <<","<< uu.imag()<<"],";
          }
          else {
            os << uu.real();
            if(l < nmo-1) os << ",";
          }
        }
        if(k==nmo-1) os << "]" << endl;
        else os << "]," << endl;
      }
      if(j==nmo-1) os << "]" << endl;
      else os << "]," << endl;
    }
    if(i==nmo-1) os << "]" << endl;
    else os << "]," << endl;
  }
  
}

//----------------------------------------------------------------------
void Average_tbdm_basis::jsonOutput(Average_return &avg,Average_return &err, ostream & os) {
  
  Array2 <dcomplex> obdm_up(nmo,nmo),obdm_down(nmo,nmo),
  obdm_up_err(nmo,nmo),obdm_down_err(nmo,nmo);
  //Array4 <dcomplex>
  //       tbdm_uu(nmo,nmo),tbdm_uu_err(nmo,nmo),
  //      tbdm_ud(nmo,nmo),tbdm_ud_err(nmo,nmo),
  //      tbdm_du(nmo,nmo),tbdm_du_err(nmo,nmo),
  //      tbdm_dd(nmo,nmo),tbdm_dd_err(nmo,nmo);
  
  os << "\"" << avg.type << "\":{" << endl;
  os << "\"nmo\":" << nmo << "," << endl;
  os << "\"states\":[";
  for(int i=0; i< nmo; i++) {
    if(i==nmo-1) os << occupations[0][i]+1;
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
      place+=2;
    }
  }
  
  os << "\"obdm\":{" << endl;
  int colwidth=40;
  if(!complex_orbitals) colwidth=20;
  
  os << "\"up\":["<< endl;
  jobdmoutput(obdm_up,os);
  os << "]," << endl;
  os << "\"up_err\":["<< endl;
  jobdmoutput(obdm_up_err,os);
  os << "]," << endl;
  os << "\"down\":["<< endl;
  jobdmoutput(obdm_down,os);
  os << "]," << endl;
  os << "\"down_err\":["<< endl;
  jobdmoutput(obdm_down_err,os);
  os << "]" << endl;
  os << "}," << endl;
  
  if(eval_tbdm) {
    os << "\"tbdm\":{" << endl;
    //dcomplex uu,uu_err,ud,ud_err,du,du_err,dd,dd_err;
    if(tbdm_diagonal) {
      os << "\"diagonal\":\"yes\"," << endl;
      
      os << "\"upup\":[" << endl;
      jdgoutput("upup",0,avg,avg,os);
      os << "]," << endl;
      os << "\"upup_err\":[" << endl;
      jdgoutput("upup_err",0,avg,err,os);
      os << "]," << endl;
      os << "\"updown\":[" << endl;
      jdgoutput("updown",1,avg,avg,os);
      os << "]," << endl;
      os << "\"updown_err\":[" << endl;
      jdgoutput("updown_err",1,avg,err,os);
      os << "]," << endl;
      os << "\"downup\":[" << endl;
      jdgoutput("downup",2,avg,avg,os);
      os << "]," << endl;
      os << "\"downup_err\":[" << endl;
      jdgoutput("downup_err",2,avg,err,os);
      os << "]," << endl;
      os << "\"downdown\":[" << endl;
      jdgoutput("downdown",3,avg,avg,os);
      os << "]," << endl;
      os << "\"downdown_err\":[" << endl;
      jdgoutput("downdown_err",3,avg,err,os);
      os << "]" << endl;
      
    }
    else {
      os << "\"diagonal\":\"no\"," << endl;
      
      os << "\"upup\":[" << endl;
      jsontbdmvalHelper(avg,tbdm_uu,os);
      os << "]," << endl;
      
      os << "\"upup_err\":[" << endl;
      jsontbdmerrHelper(avg,err,tbdm_uu,os);
      os << "]," << endl;
      
      os << "\"updown\":[" << endl;
      jsontbdmvalHelper(avg,tbdm_ud,os);
      os << "]," << endl;
      
      
      os << "\"updown_err\":[" << endl;
      jsontbdmerrHelper(avg,err,tbdm_ud,os);
      os << "]," << endl;
      
      os << "\"downup\":[" << endl;
      jsontbdmvalHelper(avg,tbdm_du,os);
      os << "]," << endl;
      
      os << "\"downup_err\":[" << endl;
      jsontbdmerrHelper(avg,err,tbdm_du,os);
      os << "]," << endl;
      
      
      os << "\"downdown\":[" << endl;
      jsontbdmvalHelper(avg,tbdm_dd,os);
      os << "]," << endl;
      
      
      os << "\"downdown_err\":[" << endl;
      jsontbdmerrHelper(avg,err,tbdm_dd,os);
      os << "]" << endl;
      
    }
    os << "}" << endl;
  }
  
  
  os << "}" << endl;
  
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
