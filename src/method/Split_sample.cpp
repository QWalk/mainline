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

#include "Split_sample.h"
#include "qmc_io.h"

int allocate(vector <string> & words, Dynamics_generator *& sam) {
  if(caseless_eq(words[0],"SPLIT"))
    sam=new Split_sampler;
  else if(caseless_eq(words[0],"RUNGE"))
    sam=new SRK_dmc;
  else if(caseless_eq(words[0],"UNR")) 
    sam=new UNR_sampler;
  else 
    error("unknown type of sampler: ", words[0]);

  sam->read(words);
  return 1;
}

//----------------------------------------------------------------------

void limDrift(Array1 <doublevar> & drift, doublevar tau, drift_type dtype)
{
  doublevar drift2=0;
  doublevar drifta=0;


  if(dtype==drift_cutoff) {
    const doublevar drmax=4;
    
    for(int d=0; d<3; d++)
        drift2+=drift(d)*drift(d);

    drifta=drmax/max( (doublevar) sqrt(drift2), drmax);
    //cout << "drifta " << drifta*tau << endl;
    for(int d=0; d<3; d++) {
        drift(d)*= tau*drifta;
      }
  }
  else if(dtype==drift_cyrus) {
    
    const doublevar acyrus=0.25;
    
    doublevar tau2=tau*2.;
    
    for(int d=0; d<3; d++)
      {
        drift2+=drift(d)*drift(d);
      }
    
    drift2*=acyrus;
    if(drift2 > 1e-16) 
      drifta=(sqrt(1.+tau2*drift2)-1)/(drift2);
    else drifta=0.0;
    
    //cout << "drifta " << drifta << endl;
    for(int d=0; d<3; d++)
      drift(d)*=drifta;
    
  }
  else {
    error("Unknown drift_type in limDrift");
  }
  
}

//-------------------------------------------------------------------

/*!
returns the exponent of the transition probability for
going from point1 to point2
 */
doublevar Split_sampler::transition_prob(int point1, int point2,
                                         doublevar timestep, 
                                         drift_type dtype) {
  doublevar prob=0;
  Array1 <doublevar> drift(3);
  //cout << "transition probability" << endl;

  drift=trace(point1).drift;

  limDrift(drift, timestep, dtype);

  for(int d=0; d< 3; d++) {
    
    prob-=(trace(point2).pos(d)-trace(point1).pos(d)-drift(d))
      *(trace(point2).pos(d)-trace(point1).pos(d)-drift(d));
  }
  prob/=(2.0*timestep);
  
  return prob;
}

doublevar transition_prob(Point &  point1, Point & point2,
                          doublevar timestep, drift_type dtype) {
  doublevar prob=0;
  Array1 <doublevar> drift(3);
  //cout << "transition probability" << endl;

  drift=point1.drift;

  
  limDrift(drift, timestep, dtype);

  for(int d=0; d< 3; d++) {
    prob-=(point2.pos(d)-point1.pos(d)-drift(d))
      *(point2.pos(d)-point1.pos(d)-drift(d));
  }
  prob/=(2.0*timestep);
  return prob;
}

//---------------------------------------------------------------------


doublevar runge_kutta_resamp(Point & p1, Point & p2,
                               doublevar timestep, drift_type dtype,
                          int ndim=3) {
  Array1 <doublevar> dr1(3), dr2(3);
  dr1=p1.drift; dr2=p2.drift;

  limDrift(dr1, timestep, dtype);
  limDrift(dr2, timestep, dtype);
  Array1 <doublevar> avg(3);
  for(int d=0; d< 3; d++) {
    avg(d)=(dr1(d)+dr2(d))/2.0;
  }
  doublevar green_forward=0;
  for(int d=0; d< ndim; d++) {
    green_forward+=(p2.pos(d)-p1.pos(d)-avg(d))*(p2.pos(d)-p1.pos(d)-avg(d));
    //green_forward+=(p2.pos(d)-p1.pos(d))*(p2.pos(d)-p1.pos(d))+avg(d)*avg(d);
  }
  return -green_forward/(2*timestep);
}


doublevar runge_kutta_symm(Point & p1, Point & p2,
                               doublevar timestep, drift_type dtype,
                          int ndim=3) {
  Array1 <doublevar> dr1(3), dr2(3);
  dr1=p1.drift; dr2=p2.drift;

  limDrift(dr1, timestep, dtype);
  limDrift(dr2, timestep, dtype);
  Array1 <doublevar> avg(3);
  for(int d=0; d< 3; d++) {
    avg(d)=(dr1(d)+dr2(d))/2.0;
  }
  //limDrift(avg,timestep, dtype);
  doublevar green_forward=0;
  for(int d=0; d< ndim; d++) {
    //green_forward+=(p2.pos(d)-p1.pos(d)-avg(d))*(p2.pos(d)-p1.pos(d)-avg(d));
    green_forward+=(p2.pos(d)-p1.pos(d))*(p2.pos(d)-p1.pos(d))+avg(d)*avg(d);
  }
  return -green_forward/(2*timestep);
}


doublevar linear_symm(Point & p1, Point & p2,
		      doublevar timestep, drift_type dtype,
		      int ndim=3) {
  Array1 <doublevar> dr1(3), dr2(3);
  dr1=p1.drift; dr2=p2.drift;

  limDrift(dr1, timestep,dtype);
  limDrift(dr2, timestep, dtype);

  doublevar green_forward=0;
  for(int d=0; d< ndim; d++) {
    
    green_forward+=(p2.pos(d)-p1.pos(d))*(p2.pos(d)-p1.pos(d))
      +(p2.pos(d)-p1.pos(d))*(dr2(d)-dr1(d))
      +.5*(dr2(d)*dr2(d)+dr1(d)*dr1(d));

    // Runge-kutta move
    //doublevar tmp=p2.pos(d)-p1.pos(d)-.5*(dr2(d)+dr1(d));
    //green_forward+=tmp*tmp;
    //green_forward+=(p2.pos(d)-p1.pos(d))*(p2.pos(d)-p1.pos(d))
    //      +.25*(dr2(d)+dr1(d))*(dr2(d)+dr1(d));
  }
  return -green_forward/(2*timestep);
}

//----------------------------------------------------------------------

void Split_sampler::read(vector <string> & words) {
  unsigned int pos=0;

  readvalue(words, pos=0, divide_, "DIVIDER");
  if(haskeyword(words, pos=0, "RESTRICT_NODES"))
    restrict_nodes=1;
  
  int depth;
  if(readvalue(words, pos=0, depth, "DEPTH")) 
    setRecursionDepth(depth);

  string drifttype;
  if(readvalue(words, pos=0, drifttype, "DRIFT_TYPE")) {
    if(drifttype=="CYRUS")
      setDriftType(drift_cyrus);
    else if(drifttype=="CUTOFF")
      setDriftType(drift_cutoff);
    else error("Didn't understand DRIFT_TYPE ", drifttype);
  }

  
}


doublevar Dynamics_generator::greenFunction(Sample_point * sample, Wavefunction * wf,
                                  Wavefunction_data * wfdata, 
                                  Guiding_function * guidingwf,
                                  int e,
                                  Array1 <doublevar> & newpos, 
                                  doublevar timestep,
                                  Dynamics_info & info, Dynamics_info & oldinfo) {
                                    
  //cout << "auxillary " << endl;
  drift_type dtype=drift_cyrus;
  assert(newpos.GetDim(0) >= 3);
  wf->updateLap(wfdata, sample);
  Point p1; p1.lap.Resize(wf->nfunc(), 5);
  p1.sign=sample->overallSign();
  int ndim=sample->ndim();

  sample->getElectronPos(e,p1.pos);
  wf->getLap(wfdata, e, p1.lap);
  guidingwf->getLap(p1.lap, p1.drift);  
  
  
  for(int d=0; d< ndim; d++) 
    p1.translation(d)=newpos(d)-p1.pos(d);
  
  sample->translateElectron(e,p1.translation);
  wf->updateLap(wfdata, sample);
  Point p2; p2.lap.Resize(wf->nfunc(), 5);
  p2.sign=sample->overallSign();
  p2.pos=newpos;
  wf->getLap(wfdata, e, p2.lap);
  guidingwf->getLap(p2.lap, p2.drift);  
  

  info.green_forward=exp(::transition_prob(p1, p2, timestep, dtype));
  info.green_backward=::transition_prob(p2,p1, timestep, dtype);
  info.symm_gf=exp(linear_symm(p1, p2, timestep, dtype));
  
  //cout << "gf:elec " << e << " new wfval " << p2.lap.amp(0,0) << endl;
  info.diffuse_start=p1.drift;
  limDrift(info.diffuse_start, timestep, dtype);
  for(int d=0; d< ndim; d++) {
    //cout << "secondary drift " << info.drift_pos(d) << endl;
    info.diffuse_start(d)+=p1.pos(d);
  }
  info.diffuse_end=newpos;
  info.orig_pos=p1.pos;
  info.new_pos=newpos;
  info.diffusion_rate=0;
  for(int d=0; d< ndim; d++) 
    info.diffusion_rate+=(info.diffuse_end(d)-info.diffuse_start(d))
                        *(info.diffuse_end(d)-info.diffuse_start(d));
  
  info.acceptance=min(1.0, (exp(info.green_backward)/info.green_forward)
                           *guidingwf->getTrialRatio(p2.lap, p1.lap)
                           *guidingwf->getTrialRatio(p2.lap, p1.lap));


  //trying a better gf
  //info.green_forward*=info.acceptance;///oldinfo.acceptance;
  //--------
  

  info.accepted=0;
  doublevar ratio=guidingwf->getTrialRatio(p1.lap,p2.lap)*p1.sign*p2.sign;
  if(ratio < 0) cout << "crossed node " << endl;
  
  return info.acceptance;
}


//----------------------------------------------------------------------
/*!
From x to y in the trace
 */
doublevar Split_sampler::get_acceptance(Guiding_function * guidingwf, 
                                        int x, int y) {
  //indent = indent + "  ";

  doublevar ratio=guidingwf->getTrialRatio(trace(y).lap, 
                                           trace(x).lap)*trace(x).sign
        *trace(y).sign;

  if(restrict_nodes && ratio < 0) return 0;

  ratio=ratio*ratio;
  //cout << indent << "*********" << endl;
  //cout << indent << "ratio x " << x << " y " << y << " : "  << ratio << endl;

  int dir=1;
  if(y-x < 0) {
    dir=-1;
  }


  doublevar prob_transition=0; //exponent of transition probability
  //We use the current move as well(thus y+1 or x-1)
  for(int i=x+dir; i != y+dir; i+= dir ) {
    int dist=abs(x-i);
    //cout << indent << "transition from " << y << " to " << y-dir*dist 
    //     << " and " << x << " to " << x+dir*dist << " using " << dist 
    //     << " (timestep " << timesteps(dist)
    //     <<  endl;

    doublevar num=transition_prob(y,y-dir*dist,timesteps(dist), dtype);
    doublevar den=transition_prob(x,x+dir*dist,timesteps(dist), dtype);
    prob_transition+=num-den;
    //cout << indent <<  "num " << num << " den " << den << " num-den " << num-den << endl;
  }
  
  doublevar reject_prob_den=1; //denominator(forward rejections)
  doublevar reject_prob_num=1; //numerator(backward rejection)
  for(int i=x+dir; i != y; i+=dir) {
    int dist=abs(x-i);
    //cout << indent << "->acceptance from " << y << " to " << y-dir*dist 
    //     << endl;
    reject_prob_num*=1-get_acceptance(guidingwf, y, y-dir*dist);
    //cout << indent << "->over " << x << " to " << x+dir*dist << endl;
    reject_prob_den*=1-get_acceptance(guidingwf, x,x+dir*dist);
  }

  doublevar acc=0;

  const doublevar tiny=1e-10;
  if(fabs(reject_prob_den) < tiny) {
    if(fabs(reject_prob_num) < tiny) acc=0;
    else {
      //If we have a zero on the denominator, we know that
      //alpha(y->x)=1/alpha(x->y) must be zero, so therefore the acceptance is 1
      acc=1;
      //cout << indent << "warn:denominator zero : " << reject_prob_den << "   num "
      //     << reject_prob_num << endl;
      //error("problem in split_sample; denominator zero");
    }
  }
  else {
    //cout << indent << "prob_transition " << prob_transition
    //    << "   " << reject_prob_num << "   " <<  reject_prob_den << endl;
    acc=ratio*exp(prob_transition)*reject_prob_num/reject_prob_den;
  }
  //cout << indent << "acceptance " << acc << endl;
  //indent.erase(indent.end()-1);
  //indent.erase(indent.end()-1);

  return min(1.0,acc);

}

#include "ulec.h"
#include "Wavefunction_data.h"


//----------------------------------------------------------------------
int Split_sampler::split_driver(int e,
                                Sample_point * sample,
                                Wavefunction * wf, 
                                Wavefunction_data * wfdata,
                                Guiding_function * guidingwf,
                                int depth,
                                Dynamics_info & info,
                                doublevar & efftimestep) {

  //cout << "primary " << endl;
  assert(trace.GetDim(0) >= depth);
  assert(recursion_depth_ <= timesteps.GetDim(0));

  if(depth > recursion_depth_) return 0;

  Array1 <doublevar> c_olddrift(3);
  Array1 <doublevar> c_newdrift(3);
  
  c_olddrift=trace(0).drift;  
  limDrift(c_olddrift, timesteps(depth), dtype);

  int ndim=sample->ndim();

  //cout << "drift " << c_olddrift(0) << "  " 
  //     << c_olddrift(1) << "  " << c_olddrift(2) << endl;
  //info.drift_pos.Resize(3);
  
  for(int d=0; d< ndim; d++) {
    trace(depth).gauss(d)=rng.gasdev();
    trace(depth).translation(d)=trace(depth).gauss(d)*sqrt(timesteps(depth))
        + c_olddrift(d);
    trace(depth).pos(d)=trace(0).pos(d)
        + trace(depth).translation(d);
    

  }

  doublevar diffusion_rate=0;
  for(int d=0; d< ndim; d++) 
    diffusion_rate+=trace(depth).gauss(d)*timesteps(depth)*trace(depth).gauss(d);;
  
  
  sample->translateElectron(e, trace(depth).translation);
  trace(depth).sign=sample->overallSign();
  
  if(wfdata->supports(laplacian_update) ) {
    wf->updateLap(wfdata, sample);
    wf->getLap(wfdata, e, trace(depth).lap);
  }
  else {
    wf->updateForceBias(wfdata, sample);
    wf->getForceBias(wfdata, e, trace(depth).lap);
  }
  
  guidingwf->getLap(trace(depth).lap, trace(depth).drift);

  //indent="";
  //cout << "#######################acceptance for " << depth << endl;
  doublevar acc=get_acceptance(guidingwf, 0,depth);
  //cout << "acceptance for " << depth << " : " <<  acc << endl;    

  info.green_forward=exp(transition_prob(0,depth,timesteps(depth), dtype));
  //cout << "green_forward " << info.green_forward << endl;
  info.green_backward=exp(transition_prob(depth,0,timesteps(depth),dtype));
  info.diffusion_rate=diffusion_rate;
  info.acceptance=acc;
  info.orig_pos=trace(0).pos;
  info.diffuse_start.Resize(3);
  for(int d=0; d< ndim; d++)
    info.diffuse_start(d)=trace(0).pos(d)+c_olddrift(d);
  info.diffuse_end=trace(depth).pos;
  info.new_pos=trace(depth).pos;
  info.gauss=trace(depth).gauss;
  
  info.symm_gf=exp(linear_symm(trace(0), trace(depth), timesteps(depth), dtype));
  info.resample_gf=info.symm_gf;
  //trying a better gf
  //info.green_forward*=info.acceptance;
  //--------
  
  
  if (acc+rng.ulec() > 1.0) {
    info.accepted=1;
    return depth;
  }
  else {
    info.accepted=0;
    
    Array1 <doublevar> rev(3,0.0);
    for(int d=0; d< 3; d++) rev(d)=-trace(depth).translation(d);
    sample->translateElectron(e,rev);

    depth++;
    
    return split_driver(e,sample, wf, wfdata, guidingwf, 
                        depth, info, efftimestep);
  }
  
}


//----------------------------------------------------------------------

int Split_sampler::sample(int e,
                          Sample_point * sample, 
                          Wavefunction * wf, 
                          Wavefunction_data * wfdata,
                          Guiding_function * guidingwf,
                          Dynamics_info & info,
                          doublevar & efftimestep) {

  if(! wfStore.isInitialized())
    wfStore.initialize(sample, wf);
  
  wf->updateLap(wfdata, sample);
  wfStore.saveUpdate(sample, wf, e);
  trace.Resize(recursion_depth_+1);

  for(int i=0; i < recursion_depth_+1; i++) {
    trace(i).lap.Resize(wf->nfunc(), 5);
  }

  timesteps.Resize(recursion_depth_+1);
  timesteps=efftimestep;

    
  for(int i=2; i< recursion_depth_+1; i++) {
    timesteps(i)=efftimestep/pow(divide_,i-1);
  }

  int depth=0;

  sample->getElectronPos(e,trace(depth).pos);
  wf->getLap(wfdata, e, trace(depth).lap);
  trace(depth).sign=sample->overallSign();

  guidingwf->getLap(trace(depth).lap, trace(depth).drift);
  depth++;

  
  int acc=split_driver(e, sample, wf, wfdata, guidingwf, depth,  
                      info, efftimestep);

  if(acc > 0) {
    acceptances(acc-1)++;
    for(int i=0; i< acc; i++) {
      tries(i)++;
    }
  }
  else {
    for(int i=0; i< recursion_depth_; i++) {
      tries(i)++;
    }
  }

  if(!acc) {
    wfStore.restoreUpdate(sample, wf, e);
  }
  //cout << "-----------split done" << endl;
  return acc;
}


//----------------------------------------------------------------------


void Split_sampler::showStats(ostream & os) {
  doublevar totacc=0;
  for(int i=0; i< recursion_depth_; i++) {
    totacc+=acceptances(i);
  }

  os << "Total acceptance " << totacc/tries(0) << endl;
  for(int i=0; i< recursion_depth_; i++) {
    os << "accept_level" << i << "   " << acceptances(i)/tries(i) 
       << " tries " << tries(i) <<  endl;
  }
}

//----------------------------------------------------------------------


void Split_sampler::resetStats() {
  acceptances=0;
  tries=0;
}
//----------------------------------------------------------------------
//######################################################################

void UNR_sampler::showStats(ostream & os) { 
  os << "acceptance " << acceptance/tries << endl;
}

void UNR_sampler::resetStats() { 
  acceptance=0;
  tries=0;
}


void UNR_sampler::getDriftEtc(Point & pt, Sample_point * sample,
			      doublevar tstep, int e,
			      doublevar & exponent, 
			      //!< exponent at which to do an exponential move
			      doublevar & probability,
			      //!< probability to do an exponential move
			      Array1 <doublevar> & rnuc
			      //!< distance array for closest ion
			      ) {
  Array1 <doublevar> z(5);
  Array1 <doublevar> dist(5);
  int closest_ion=-1;
  z(0)=1e8;
  int nions=sample->ionSize();
  sample->updateEIDist();
  for(int at=0; at< nions; at++) {
    sample->getEIDist(e,at,dist);
    if(dist(0) < z(0)) { 
      z=dist;
      closest_ion=at;
    }
  }

  assert(closest_ion!=-1);
  doublevar charge=sample->getIonCharge(closest_ion);
  rnuc.Resize(3);
  //Here we back out the nucleus position from the distance,
  //which is the best way to do it for periodic systems.
  //There may be a small bug here if the periodic cell is really
  //small and we can traverse it in one step..there may be extra
  //rejections.  It shouldn't generally be a concern, though.
  for(int d=0; d< 3; d++) { 
    rnuc(d)=pt.pos(d)-z(d+2);
  }

  //the exponent for the exponential moves
  exponent=sqrt(charge*charge+1/tstep);

  doublevar vhatdotzhat=0;
  doublevar driftmag=0;
  for(int d=0; d< 3; d++) driftmag+=pt.drift(d)*pt.drift(d);
  driftmag=sqrt(driftmag);
  for(int d=0;d < 3; d++) {
    vhatdotzhat+=z(d+2)*pt.drift(d)/(driftmag*z(0));
    //cout << "vhat " << vhatdotzhat << endl;
  }

  //magic numbers!!  taken from UNR as usual
  doublevar a=.5*(1+vhatdotzhat)
    +charge*charge*z(0)*z(0)/(10*(4+charge*charge*z(0)*z(0)));
  
  doublevar prefac=(-1+sqrt(1+2*a*driftmag*driftmag*tstep))
    /(a*driftmag*driftmag*tstep);

  driftmag=0;
  for(int d=0;d < 3; d++) {
    pt.drift(d)*=prefac;
    driftmag+=pt.drift(d)*pt.drift(d);
  }
  driftmag=sqrt(driftmag);
  
  doublevar vz=vhatdotzhat*driftmag;
  Array1 <doublevar> rho(3);
  rho=0;
  doublevar rhomag=0;
  for(int d=0; d< 3; d++) {
    rho(d)=pt.drift(d)-vz*z(d+2)/z(0);
    rhomag+=rho(d)*rho(d);
  }
  rhomag=sqrt(rhomag);
  doublevar zpp=max(z(0)+vz*tstep,0.0);
  doublevar rhoscale=2*tstep*zpp/(z(0)+zpp);
  for(int d=0;d < 3; d++) {
    pt.drift(d)=rnuc(d)+rhoscale*rho(d)+zpp*z(d+2)/z(0);
  }

  probability=0.5*erfc((z(0)+vz*tstep)/sqrt(2*tstep)); 

}

//----------------------------------------------------------------------

doublevar UNRgreenfunc(Point & pt1, Point & pt2,
		       doublevar tstep,
		       doublevar exponent, doublevar prob,
		       Array1 <doublevar> & rnuc) { 
  assert(rnuc.GetDim(0)==3);
  doublevar gaussmov=0;
  doublevar expmov=0;
  for(int d=0; d< 3; d++) {
    gaussmov+=(pt2.pos(d)-pt1.drift(d))*(pt2.pos(d)-pt1.drift(d));
    expmov+=(pt2.pos(d)-rnuc(d))*(pt2.pos(d)-rnuc(d));
    //cout << "rnuc " << rnuc(d) << "  pt2 " << pt2.pos(d) << " pt1 " << pt1.drift(d) << endl;
  }
  doublevar gaussnorm=sqrt(8*pi*pi*pi*tstep*tstep*tstep);
  gaussnorm=1/gaussnorm;
  

  return prob*exponent*exponent*exponent*exp(-2*exponent*sqrt(expmov))/pi
    + (1-prob)*exp(-gaussmov/(2*tstep))*gaussnorm;
}
//----------------------------------------------------------------------

int UNR_sampler::sample(int e, Sample_point * sample,
			Wavefunction * wf, Wavefunction_data * wfdata,
			Guiding_function * guidewf, 
			Dynamics_info & info,
			doublevar & tstep) {
    
  tries++;
  Point p1; p1.lap.Resize(wf->nfunc(), 5);
  
  wf->updateLap(wfdata, sample);
  sample->getElectronPos(e,p1.pos);
  wf->getLap(wfdata, e, p1.lap);
  p1.sign=sample->overallSign();
  guidewf->getLap(p1.lap, p1.drift);

  //Starting determination of move..
  //using (close to) the notation of UNR
  Array1 <doublevar> rnuc(3);
  doublevar prob,exponent;

  //cout << "p1.drift "; write_array(cout,p1.drift); cout << endl;
  getDriftEtc(p1,sample,tstep,e, exponent,prob,rnuc);
  
  //cout << "p1.pos ";
  //write_array(cout, p1.pos);
  //cout << endl;
  //cout << "p1.drift "; write_array(cout,p1.drift); cout << endl;
  

  Point p2; p2.lap.Resize(wf->nfunc(),5);
  //cout << "prob " << prob << endl;
  if(prob+rng.ulec() > 1.0) { 
    //cout << "doing exponential move "<< endl;
    doublevar r=ranr2exponential()/(2*exponent);
    doublevar theta=rancos();
    doublevar phi=2*pi*rng.ulec();
    p2.pos=rnuc;
    p2.pos(0)+=r*sin(theta)*cos(phi);
    p2.pos(1)+=r*sin(theta)*sin(phi);
    p2.pos(2)+=r*cos(theta);
    
  }
  else { 
    //cout << " do a gaussian move." << endl;
    for(int d=0; d< 3; d++) {
      p2.pos(d)=p1.drift(d)+sqrt(tstep)*rng.gasdev();
    }
  }

  //cout << "p2.pos "; write_array(cout,p2.pos); cout << endl;

  Array1 <doublevar> translate(3);
  for(int d=0; d< 3; d++) translate(d)=p2.pos(d)-p1.pos(d);
  sample->translateElectron(e,translate);
  wf->updateLap(wfdata, sample);
  wf->getLap(wfdata, e, p2.lap);
  p2.sign=sample->overallSign();
  guidewf->getLap(p2.lap, p2.drift);

  Array1 <doublevar> rnucp(3);
  doublevar probp,exponentp;
  getDriftEtc(p2,sample,tstep,e,exponentp,probp,rnucp);

  //cout << "p2.drift "; write_array(cout,p2.drift); cout << endl;


  doublevar ratio=guidewf->getTrialRatio(p2.lap,p1.lap)
    *p2.sign*p1.sign;

  doublevar acc=ratio*ratio*UNRgreenfunc(p2,p1,tstep,exponentp,probp,rnucp)
    /UNRgreenfunc(p1,p2,tstep,exponent,prob,rnuc);

  //cout << "forward " << UNRgreenfunc(p2,p1,tstep,exponentp,probp,rnucp) << endl;
  //cout << "acceptance " << acc << " ratio " << ratio 
  //  << "forward " << UNRgreenfunc(p2,p1,tstep,exponentp,probp,rnucp)
   // << " backward " << UNRgreenfunc(p1,p2,tstep,exponent,prob,rnuc) << endl;
  
  if(restrict_nodes && ratio < 0) {
    //cout << "stopping nodal crossing " << endl;
    acc=0;
  }

  //cout << "acceptance " << acc << endl;
  //cout <<"##############################\n";

  info.diffusion_rate=0;
  for(int d=0;d < 3; d++) {
    info.diffusion_rate+=(p2.pos(d)-p1.drift(d))*(p2.pos(d)-p1.drift(d))*tstep;
  }
  info.acceptance=acc;
  info.orig_pos=p1.pos;
  info.diffuse_start=p1.drift;
  info.diffuse_end=p2.pos;
  info.new_pos=p2.pos;
  
  if(acc+rng.ulec()>1.0) { 
    info.accepted=1;
    acceptance++;
    return 1;
  }
  else { 
    sample->setElectronPos(e,p1.pos);
    info.accepted=0;
    return 0;
  }
  
  
}

//----------------------------------------------------------------------
//######################################################################
void SRK_dmc::read(vector <string> & words) {
  unsigned int pos=0;
  if(haskeyword(words, pos=0, "RESAMPLE"))
    resample=1;
  else resample=0;
  
  readvalue(words, pos=0, resample_tol, "RESAMPLE_TOL");

  string drifttype;
  if(readvalue(words, pos=0, drifttype, "DRIFT_TYPE")) {
    if(drifttype=="CYRUS")
      setDriftType(drift_cyrus);
    else if(drifttype=="CUTOFF")
      setDriftType(drift_cutoff);
    else error("Didn't understand DRIFT_TYPE ", drifttype);
  }
  

}




doublevar SRK_dmc::greenFunction(Sample_point * sample, Wavefunction * wf,
                                  Wavefunction_data * wfdata, 
                                  Guiding_function * guidingwf,
                                  int e,
                                  Array1 <doublevar> & newpos, 
                                  doublevar timestep,
                                  Dynamics_info & info, Dynamics_info & oldinfo) {
                                    
  //cout << "auxillary " << endl;
  //drift_type dtype=drift_cyrus;
  assert(newpos.GetDim(0) >= 3);
  wf->updateLap(wfdata, sample);
  Point p1; p1.lap.Resize(wf->nfunc(), 5);
  p1.sign=sample->overallSign();
  int ndim=sample->ndim();

  sample->getElectronPos(e,p1.pos);
  wf->getLap(wfdata, e, p1.lap);
  guidingwf->getLap(p1.lap, p1.drift);  
  
  
  for(int d=0; d< ndim; d++) 
    p1.translation(d)=newpos(d)-p1.pos(d);
  
  sample->translateElectron(e,p1.translation);
  wf->updateLap(wfdata, sample);
  Point p2; p2.lap.Resize(wf->nfunc(), 5);
  p2.sign=sample->overallSign();
  p2.pos=newpos;
  wf->getLap(wfdata, e, p2.lap);
  guidingwf->getLap(p2.lap, p2.drift);  
  

  info.green_forward=exp(runge_kutta_resamp(p1, p2, timestep, dtype));

  info.symm_gf=exp(runge_kutta_symm(p1, p2, timestep, dtype));
  info.green_backward=exp(transition_prob(p2,p1, timestep, dtype));

  //cout << "gf:elec " << e << " new wfval " << p2.lap.amp(0,0) << endl;
  info.diffuse_start=p1.drift;
  limDrift(info.diffuse_start, timestep, dtype);
  for(int d=0; d< ndim; d++) {
    //cout << "secondary drift " << info.drift_pos(d) << endl;
    info.diffuse_start(d)+=p1.pos(d);
  }
  info.diffuse_end=newpos;
  info.orig_pos=p1.pos;
  info.new_pos=newpos;
  info.diffusion_rate=0;
  for(int d=0; d< ndim; d++) 
    info.diffusion_rate+=(info.diffuse_end(d)-info.diffuse_start(d))
                        *(info.diffuse_end(d)-info.diffuse_start(d));
  
  info.acceptance=min(1.0, (exp(info.green_backward)/info.green_forward)
                           *guidingwf->getTrialRatio(p2.lap, p1.lap)
                           *guidingwf->getTrialRatio(p2.lap, p1.lap));
  info.accepted=0;
  doublevar ratio=guidingwf->getTrialRatio(p1.lap,p2.lap)*p1.sign*p2.sign;
  if(ratio < 0) cout << "crossed node " << endl;
  
  return info.acceptance;
}

//----------------------------------------------------------

int SRK_dmc::rk_step(int e,
                          Sample_point * sample, 
                          Wavefunction * wf, 
                          Wavefunction_data * wfdata, 
                          Guiding_function * guidingwf,
                          Dynamics_info & info,
                          doublevar & timestep, Array1 <Point> & trace
                          ) {
 

  //drift_type dtype=drift_cyrus;
  assert(trace.GetDim(0)==3);          
  
  Array1 <doublevar> c_olddrift=trace(0).drift;
  limDrift(c_olddrift,timestep, dtype);
  int ndim=sample->ndim();

  for(int d=0; d< ndim; d++) {
    trace(1).translation(d)=trace(0).gauss(d)*sqrt(timestep)
        + c_olddrift(d);
    trace(1).pos(d)=trace(0).pos(d)
        + trace(1).translation(d);
  }

  info.diffusion_rate=0;
  for(int d=0; d< ndim; d++) 
    info.diffusion_rate+=trace(0).gauss(d)*timestep*trace(0).gauss(d);;
  
  sample->translateElectron(e, trace(1).translation);
  trace(1).sign=sample->overallSign();
  wf->updateLap(wfdata, sample);
  wf->getLap(wfdata, e, trace(1).lap);
  guidingwf->getLap(trace(1).lap, trace(1).drift);
  Array1 <doublevar> c_newdrift=trace(1).drift;
  limDrift(c_newdrift,timestep, dtype);

  for(int d=0; d< 3; d++) {
    trace(2).gauss(d)=trace(0).gauss(d);
    trace(2).translation(d)=.5*(c_newdrift(d)+c_olddrift(d))
        +trace(0).gauss(d)*sqrt(timestep)-trace(1).translation(d);
    trace(2).pos(d)=trace(1).pos(d)+trace(2).translation(d);
  }
  
  sample->translateElectron(e,trace(2).translation);
  trace(2).sign=sample->overallSign();
  wf->updateLap(wfdata, sample);
  wf->getLap(wfdata, e, trace(2).lap);
  guidingwf->getLap(trace(2).lap, trace(2).drift);
  
  info.symm_gf=exp(runge_kutta_symm(trace(0), trace(2), timestep, dtype));
  info.resample_gf=exp(runge_kutta_resamp(trace(0), trace(2), timestep, dtype));
  info.green_forward=exp(-info.diffusion_rate/(2*timestep));
  //cout << "symm " << info.symm_gf << "  forward " << info.green_forward
  //      << " ratio " << info.symm_gf/info.green_forward << endl;
  info.green_backward=transition_prob(trace(2),trace(0),timestep,dtype);
  info.acceptance=1.0;
  info.orig_pos=trace(0).pos;
  info.new_pos=trace(2).pos;
   
  return 1;                       
                           
 }

 
//----------------------------------------------------------------------

int SRK_dmc::sample(int e,
                          Sample_point * sample, 
                          Wavefunction * wf, 
                          Wavefunction_data * wfdata, 
                          Guiding_function * guidingwf,
                          Dynamics_info & info,
                          doublevar & timestep
                          ) {

  //drift_type dtype=drift_cyrus;
  tries++;
  int ntrace=40;
  Array1 <Point> trace(ntrace);
  for(int i=0; i < ntrace; i++) {
    trace(i).lap.Resize(wf->nfunc(), 5);
  }
  for(int d=0; d< 3; d++) {
    trace(0).gauss(d)=rng.gasdev();
  }

  
  wf->updateLap(wfdata, sample);
  sample->getElectronPos(e,trace(0).pos);
  wf->getLap(wfdata, e, trace(0).lap);
  trace(0).sign=sample->overallSign();
  guidingwf->getLap(trace(0).lap, trace(0).drift);
  
  //cout << "initial rk step " << endl;
  rk_step(e,sample, wf, wfdata, guidingwf, info, timestep, trace);
  
  /*  
  int pos=3;
  Array1 <doublevar> err(ntrace,0.0);
  err(2)=fabs(info.resample_gf/info.green_forward -1);
  while(resample && 
        fabs(info.resample_gf/info.green_forward-1) > resample_tol && 
        pos < ntrace) {
    //cout << "big error " << info.symm_gf/info.green_forward << endl;
    Array1 <doublevar> newdrift(3), orig_drift(3);
    orig_drift=trace(0).drift;
    newdrift=trace(pos-1).drift;
    limDrift(orig_drift, timestep, dtype);
    limDrift(newdrift, timestep, dtype);
    for(int d=0; d< 3; d++) {
      trace(pos).gauss(d)=trace(0).gauss(d);
      trace(pos).translation(d)=.5*(newdrift(d)+orig_drift(d))
          +trace(0).gauss(d)*sqrt(timestep);//-trace(pos-1).translation(d);
      trace(pos).pos(d)=trace(0).pos(d)+trace(pos).translation(d);
    }
    
    sample->setElectronPos(e,trace(pos).pos);
    trace(pos).sign=sample->overallSign();
    wf->updateLap(wfdata, sample);
    wf->getLap(wfdata, e, trace(pos).lap);
    guidingwf->getLap(trace(pos).lap, trace(pos).drift);
    
    info.symm_gf=exp(runge_kutta_symm(trace(0), trace(pos), timestep, dtype));
    info.resample_gf=exp(runge_kutta_resamp(trace(0), trace(pos), 
                                            timestep, dtype));
    info.new_pos=trace(pos).pos;

    err(pos)=fabs(info.resample_gf/info.green_forward-1);
    pos++;
  }
  
  if(pos==ntrace) { 
    nbottom++;
    int bestindex=2; doublevar bestratio=err(bestindex);
    for(int i=2; i < ntrace; i++) {
      if(err(i) < bestratio) { bestindex=i; bestratio=err(i); } 
    }
    if(bestindex!=pos-1) {
      //cout << " large remaining error "; write_array(cout, err); cout << endl;
      //cout << "best index " << bestindex << endl;
    
    
      info.symm_gf=exp(runge_kutta_symm(trace(0), trace(bestindex), 
                                        timestep, dtype));
      info.resample_gf=exp(runge_kutta_resamp(trace(0), trace(bestindex), 
                                        timestep, dtype));
      info.new_pos=trace(bestindex).pos;
      sample->setElectronPos(e,trace(bestindex).pos);
      pos=bestindex+1;
    }
  }
  */
  //if(pos > 3) cout << "final ratio " << info.symm_gf/info.green_forward << endl;
  
  doublevar ratio=guidingwf->getTrialRatio(trace(2).lap,
                                           trace(0).lap)
    *trace(0).sign*trace(2).sign;

  
  info.diffuse_start.Resize(3);
  for(int d=0; d< 3; d++) 
    info.diffuse_start(d)=trace(2).pos(d)-sqrt(timestep)*trace(0).gauss(d);
  
  info.diffuse_end=trace(2).pos;

  //retries+=pos-3;
  if(ratio < 0) {
    //cout << "rejecting node crossing " <<  endl;
    info.accepted=0;
    sample->setElectronPos(e, trace(0).pos);
    wf->updateLap(wfdata, sample);
    return 0;
  } 
  acceptances++;
  info.accepted=1;
  return 1;
}


