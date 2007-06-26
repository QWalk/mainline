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
#include "Properties_gather.h"
#include "ulec.h"
#include "qmc_io.h"
#include "Wavefunction_data.h"
//----------------------------------------------------------------------

Properties_gather::~Properties_gather() {
  int naux=aux_sys.GetDim(0);
  for(int i=0; i< naux; i++) {
    if(aux_sys(i)) delete aux_sys(i);
    if(aux_wf(i)) delete aux_wf(i);
    if(aux_wfdata(i)) delete aux_wfdata(i);
    if(aux_sample(i)) delete aux_sample(i);
  }
}

//----------------------------------------------------------------------

void Properties_gather::read(vector < string> & words) {
  
  unsigned int pos=0;

  vector < vector < string> > aux_systxt;
  vector < vector < string > > aux_wftxt;
  vector < string> tmp, tmp2;

  pos=0;
  while(readsection(words, pos, tmp, "AUX_SYS")) {
    unsigned int pos2=0;
    if(readsection(tmp, pos2, tmp2, "SYSTEM"))
      aux_systxt.push_back(tmp2);
    else 
      aux_systxt.push_back(tmp);
  }

  pos=0;
  while(readsection(words, pos, tmp, "AUX_WF")) 
    aux_wftxt.push_back(tmp);
  

  if(aux_systxt.size() != aux_wftxt.size())
    error("the number of AUX_SYS and AUX_WF must be the same in properties");

  vector <string> warpsec;
  if(readsection(words, pos=0, warpsec,"WARPER")) {
    warper.read(warpsec);
  }

  if(haskeyword(words, pos=0, "NOWARP")) 
    warper.set_warp(0);
  
  zpol_manye=0;
  if(haskeyword(words, pos=0, "ZPOL_manye"))
    zpol_manye=1;


  int naux=aux_systxt.size();

  aux_wfdata.Resize(naux);
  aux_sys.Resize(naux);
  aux_wf.Resize(naux);
  aux_sample.Resize(naux);

  aux_wfdata=NULL; aux_sys=NULL;
  aux_wf=NULL; aux_sample=NULL;

  for(int i=0; i < naux; i++) {
    allocate(aux_systxt[i], aux_sys(i));
    allocate(aux_wftxt[i], aux_sys(i), aux_wfdata(i));
    aux_sys(i)->generateSample(aux_sample(i));
    aux_wfdata(i)->generateWavefunction(aux_wf(i));
    aux_sample(i)->attachObserver(aux_wf(i));
  }  
  

}


//----------------------------------------------------------------------
#include "Split_sample.h"

void Properties_gather::updateAuxFunctions(Sample_point * sample) {
  int naux=aux_sys.GetDim(0);
  Array1 <doublevar> pos(3); doublevar jacobian;

  Array1 <doublevar> ref_pos(3);
  for(int i=0; i< naux; i++) {
    for(int e=0; e < sample->electronSize() ; e++) {
      sample->getElectronPos(e,ref_pos);
      warper.space_warp(sample, aux_sample(i), e,ref_pos,pos, jacobian);
      aux_sample(i)->setElectronPos(e,pos);
    }
    aux_wf(i)->notify(all_electrons_move, 0);
    aux_wf(i)->updateLap(aux_wfdata(i), aux_sample(i));
  }
}


void Properties_gather::getGreensFunctions(Dynamics_generator * dyngen,
                                           Dynamics_info & dinfo,
                                           int e, Sample_point * sample,
                                           Guiding_function * guide,
                                           doublevar timestep,
                                           Array1 <Dynamics_info> & aux_dinfo,
                                           int eval_gf) {
  doublevar jacobian;
  //Array1 <doublevar> pos(3);
  int naux=aux_sys.GetDim(0);
  //gf.Resize(naux);
 
  //Array1 <doublevar> ref_pos(3);
  aux_dinfo.Resize(naux);
  for(int i=0; i< naux; i++) {
    
    Wavefunction_data * usewfdata=NULL;
    if(aux_wfdata(i)==NULL) 
      error("bad assumption in Properties_gather::getGreensFunctions");
    else usewfdata=aux_wfdata(i);
    int ndim=aux_sys(i)->ndim();
    
    Array1 <doublevar> test_pos(3);
    aux_sample(i)->getElectronPos(e,test_pos);
        
    Array1 <doublevar> new_pos(3);
    warper.space_warp(sample, aux_sample(i), e,dinfo.new_pos, new_pos, jacobian);
    
    Array1 <doublevar> ref_pos(3);
    warper.space_warp(sample, aux_sample(i), e,dinfo.orig_pos, ref_pos, jacobian);
    
    
    if(dinfo.accepted && eval_gf) {
      //cout << "moving " << e << " to " << new_pos(0) << endl;
    
      dyngen->greenFunction(aux_sample(i), aux_wf(i), usewfdata,
                            guide,e, new_pos, timestep, aux_dinfo(i), dinfo);
    
    }
    else { //cout << "not moving " << e<< " at " << test_pos(0) << endl; 
           aux_dinfo(i).green_forward=1;}
      
    
    aux_dinfo(i).accepted=dinfo.accepted;

    Array1 <doublevar> diffuse_start(3), diffuse_end(3);
    warper.space_warp(sample, aux_sample(i), e,dinfo.diffuse_start, diffuse_start, jacobian);
    warper.space_warp(sample, aux_sample(i), e,dinfo.diffuse_end, diffuse_end, jacobian);
    aux_dinfo(i).diffusion_rate=0;
    for(int d=0; d< ndim; d++)
      aux_dinfo(i).diffusion_rate+=(diffuse_end(d)-diffuse_start(d))
          *(diffuse_end(d)-diffuse_start(d));
                                      
  }
 
}

//----------------------------------------------------------------------
#include "Force_fitter.h"

void getZpol(System * sys, Sample_point * sample, Array1 <dcomplex> & zpol,
	     int zpol_manye) { 

  int nelectrons=sample->electronSize();
  Array2 <doublevar> gvec;
  Array1 <doublevar> pos(3);
  zpol.Resize(3);
  
  if(sys->getRecipLattice(gvec)) {
    zpol=dcomplex(0.0,0.0);
    if(zpol_manye) { 
      Array1 <doublevar> sum(3,0.0);
      //See Souza et al. PRB v 62 pp 1666
      //We're using the relation given in section
      //VII.B 
      for(int e=0; e< nelectrons; e++) {
        sample->getElectronPos(e,pos);
        for(int i=0; i< 3; i++) {
          for(int d=0; d< 3; d++) {
            sum(i)+=gvec(i,d)*pos(d);
          }
        }
      }
      
      //cout << "pos " << pos(2) << "  sum " << sum(2) << endl;
      for(int i=0; i< 3; i++) {
        zpol(i)=dcomplex(cos(2*pi*sum(i)),sin(2*pi*sum(i)));
      }
    }
    else { 
      sys->getPrimRecipLattice(gvec);
      for(int e=0; e< nelectrons; e++) { 
        sample->getElectronPos(e,pos);
        for(int i=0; i< 3; i++) {
          doublevar tmp=0;
          for(int d=0; d< 3; d++) {
            tmp+=gvec(i,d)*pos(d);
          }
          zpol(i)+=dcomplex(cos(2*pi*tmp),sin(2*pi*tmp))/doublevar(nelectrons);
        }
      }
    }

  }
  else {
    Array1 <doublevar> sum(3,0.0);
    for(int e=0; e< nelectrons; e++) {
      sample->getElectronPos(e,pos);
      for(int d=0; d< 3; d++) {
        sum(d)-=pos(d);
      }
    }
    int nions=sample->ionSize();
    for(int at=0; at < nions; at++) {
      sample->getIonPos(at,pos);
      doublevar charge=sample->getIonCharge(at);
      for(int d=0; d< 3; d++) {
        sum(d)+=charge*pos(d);
      }
    }
    for(int d=0; d< 3; d++)
      zpol(d)=dcomplex(sum(d),0.0);
  }
    
}

/*!

 */

void Properties_gather::extendedGather(Properties_point & myprop,
                                   Pseudopotential * psp, 
                                   System * sys, 
                                   Wavefunction_data * wfdata,
                                   Wavefunction * wf, 
                                   Sample_point * sample, 
                                   Guiding_function * guide, int n_converge,
				       int aux_updated,
				       Array2 <doublevar> & drift,
				       Array1 < Array2 <doublevar> > & aux_drift
,
				       Array1 <Array2 <doublevar> > & aux_positions) { 
  int nwf=wf->nfunc();
  int naux=aux_sys.GetDim(0);

  //updated all the aux wfs, etc.
  gatherData(myprop, psp, sys, wfdata, wf, sample, guide, n_converge,
	     aux_updated);
  Wf_return temp(nwf,5);
  
  int nelectrons=sample->electronSize();
  drift.Resize(nelectrons,3);
  for(int e=0; e< nelectrons; e++) { 
    wf->getLap(wfdata, e, temp);
    for(int d=0; d< 3; d++) { 
      drift(e,d)=temp.amp(0,d+1);
    }
  }
  aux_drift.Resize(naux);
  aux_positions.Resize(naux);
  Array1 <doublevar> apos(3);
  for(int a=0; a< naux; a++) { 
    aux_positions(a).Resize(nelectrons,3);
    aux_drift(a).Resize(nelectrons,3);
    for(int e=0;e < nelectrons; e++) { 
      aux_sample(a)->getElectronPos(e,apos);
      aux_wf(a)->getLap(aux_wfdata(a),e,temp);
      for(int d=0; d< 3; d++) { 
	aux_positions(a)(e,d)=apos(d);
	aux_drift(a)(e,d)=temp.amp(0,d+1);
      }
    }
  }
}


void Properties_gather::gatherData(Properties_point & myprop,
                                   Pseudopotential * psp, 
                                   System * sys, 
                                   Wavefunction_data * wfdata,
                                   Wavefunction * wf, 
                                   Sample_point * sample, 
                                   Guiding_function * guide, int n_converge,
                                   int aux_updated) {
  Force_fitter force_fitter;
  force_fitter.setup(1.0,5);
  
  int nwf=wf->nfunc();
  int naux=aux_sys.GetDim(0);
  
  myprop.setSize(nwf, naux, n_converge);
  

  Array2 <doublevar> deriv_temp; //We'll just toss the derivative..

  wf->updateLap(wfdata, sample);
  wf->getVal(wfdata, 0, myprop.wf_val);

  sys->calcKinetic(wfdata, sample, wf, myprop.kinetic);
  myprop.potential=sys->calcLoc(sample);

  for(int w=0; w< nwf; w++) {
    myprop.weight(w)=guide->getWeight(myprop.wf_val, 
                                      myprop.wf_val, w);
    if(square_weight) 
      myprop.weight(w)*=myprop.weight(w);
      
  }

  int nrandvar=psp->nTest();
  Array1 <doublevar> rand_num(nrandvar);
  for(int i=0; i< nrandvar; i++) 
    rand_num(i)=rng.ulec();

  psp->calcNonlocWithTest(wfdata, sample, wf,force_fitter,
                          rand_num,  myprop.nonlocal,
                          deriv_temp);

  
  getZpol(sys, sample, myprop.z_pol, zpol_manye);
  

  //Collect auxillary points

  Array1 <doublevar> a_kin(nwf);
  Array1 <doublevar> a_nonloc(nwf);
  doublevar a_loc;
  int nelectrons=sys->nelectrons(0)+sys->nelectrons(1);
  Array1 <doublevar> pos(3);
  
  
  
  
  
  //cout << "----------------------Auxillary points " << endl;

  for(int i=0; i< naux; i++) {
    Wavefunction_data * usewfdata;
    if(aux_wfdata(i)==NULL) usewfdata=wfdata;
    else usewfdata=aux_wfdata(i);

    doublevar tot_jacob=1;
    Array1 <doublevar> ref_drift(3), aux_drift(3);
    Array1 <doublevar> ref_pos(3), orig_pos(3);

    for(int e=0;e< nelectrons; e++) {
      doublevar jacobian=1;
      sample->getElectronPos(e,ref_pos);
      aux_sample(i)->getElectronPos(e,orig_pos);
      warper.space_warp(sample, aux_sample(i), 
                        e,ref_pos, pos, jacobian);
      tot_jacob*=jacobian;

      if(!aux_updated)
        aux_sample(i)->setElectronPos(e,pos);
      
      //for(int d=0; d< 3; d++) 
      //   cout << "e=" << e << "  diff " << pos(d) << "  " << orig_pos(d) << endl;
    }
    if(!aux_updated) {
      aux_wf(i)->notify(all_electrons_move, 0);
      aux_wf(i)->updateLap(usewfdata, aux_sample(i));
    }
    aux_wf(i)->getVal(usewfdata, 0, myprop.aux_wf_val(i));  
    //cout << "pg: new wfval " << myprop.aux_wf_val(i).amp(0,0) << endl;
    aux_sys(i)->calcKinetic(usewfdata, aux_sample(i), 
                            aux_wf(i), a_kin);
    psp->calcNonlocWithTest(usewfdata, aux_sample(i),
                            aux_wf(i), force_fitter,
                            rand_num, a_nonloc, deriv_temp);
    a_loc=aux_sys(i)->calcLoc(aux_sample(i));
        
    
    myprop.aux_jacobian(i)=tot_jacob;

    Array1<dcomplex> z_pol_tmp(3);
    getZpol(aux_sys(i), aux_sample(i),z_pol_tmp, zpol_manye);
    for(int d=0; d< 3; d++) 
      myprop.aux_z_pol(i,d)=z_pol_tmp(d);

    for(int w=0; w< n_converge; w++) {
      myprop.aux_energy(i,w)=a_kin(0)+a_nonloc(0)+a_loc;
      myprop.aux_weight(i,w)=guide->getWeight(myprop.aux_wf_val(i), 
					      myprop.wf_val, 0);
      myprop.aux_weight(i,w)*=myprop.aux_weight(i,w);
      myprop.aux_weight(i,w)*=tot_jacob;
    }
  }

  //for(int i=0; i< naux; i++) {
  //  cout << "aux en " << myprop.aux_energy(i,0) 
  //      << "  weight " << myprop.aux_weight(i,0) 
  //      << "  val " << myprop.aux_wf_val(i).amp(0,0)
  //      << "  jacobian " << myprop.aux_jacobian(i) << endl;
  //}
  
  myprop.count=1;


}
