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

/* 

Homogeneous electron gas, created as a simplification of crystalline
Periodic_system.

*/

#include "HEG_system.h"
#include "Array.h"
#include "Wavefunction.h"
#include "Sample_point.h"
#include "HEG_sample.h"
#include "qmc_io.h"
#include "MatrixAlgebra.h"

void HEG_system::notify(change_type change, int n)
{
  switch(change)
  {
  default:
    cout << "WARNING: HEG system got a signal that it doesn't know: "
	 << change << endl;
  }
}


int HEG_system::generateSample(Sample_point *& samptr)
{
  assert(samptr==NULL);
  samptr=new HEG_sample;
  samptr->init(this);
  return 1;
}

int HEG_system::showinfo(ostream & os)
{
  const int ndim=3;
  os << "Homogeneous periodic system (HEG)" << endl;

  // Perhaps we do not need to be this general here, but since it is already
  // coded ...
  os << "Lattice vectors:" << endl;
  for(int i=0; i< ndim; i++) {
    os << i << " : ";
    for(int j=0; j < ndim; j++) {
      os << latVec(i,j) << "   " ;
    }
    os << endl;
  }
  os << "reciprocal vectors(multiplied by 2*pi) : " << endl;
  for(int i=0; i< ndim; i++) {
    os << i << " : ";
    for(int j=0; j< ndim; j++) {
      os << 2*pi*recipLatVec(i,j) << "    ";
    }
    os << endl;
  }

  os << endl;
  os << "Model for \"periodic\" e-e interaction used:" << endl;
  switch (eeModel) {
  case 1:
    os << "Ewald" << endl;
    break;
  case 2:
    os << "truncated Coulomb interaction, Fraser et al., PRB 53, 1814 (1994)"
       << endl;
    break;
  case 3:
    os << "Gaussian short-range interaction" << endl;
    os << "  amplitude      " << Gauss_a << endl;
    os << "  std. deviation " << Gauss_s << endl;
  }
  os << endl;

  if (!same_spin_int)
    os << "Interaction between like spins switched off." << endl << endl;
  if (!diff_spin_int)
    os << "Interaction between unlike spins switched off." << endl << endl;

  if ( eeModel == 1 ) {
    os << "total number of points in reciprocal ewald sum: "
       << ngpoints << endl;
    os << "Self e-e " << self_ee << endl;
    os << "xc correction " << xc_correction << endl;
  }

  doublevar rs=(3./(4*pi*totnelectrons/cellVolume));
  rs=pow(rs, 1.0/3.0);
  os << "rs "  << rs << endl;

  // mean-fields should be close for small rs (high densities)
  os << "approximate energies (for unpolarized Coulomb system of given rs)" << endl;
  doublevar ene=pow(9*pi/4,2.0/3)*3.0/(10*rs*rs) 
    - pow(9*pi/4,1.0/3)*3/(pi*4*rs);
  os << "  Hartree-Fock         "
     << totnelectrons*ene << " (" << ene << " per electron)" << endl;
  ene=ene + (0.0622*log(rs) - 0.096)/2.0;
  os << "  Gell-Mann & Bruecker "
     <<  totnelectrons*ene << " (" << ene << " per electron)" << endl;

  os << endl;

  return 1;
}


//----------------------------------------------------------------------
/*!
*/
int HEG_system::read(vector <string> & words,
                          unsigned int & pos)
{
  const int ndim=3;
  int startpos=pos;

  vector <string> latvectxt;

  vector <string> spintxt;
  if(!readsection(words, pos, spintxt, "NSPIN")) {
    error("Need NSPIN in HEG system");
  }
  nspin.Resize(2);
  nspin(0)=atoi(spintxt[0].c_str());
  nspin(1)=atoi(spintxt[1].c_str());
  totnelectrons=nspin(0)+nspin(1);

  vector <string> ktxt;
  if(readsection(words, pos=0, ktxt, "KPOINT")) {
    if(ktxt.size()!=3) error("KPOINT must be a section of size 3");
    kpt.Resize(3);
    for(int i=0; i< 3; i++) 
      kpt(i)=atof(ktxt[i].c_str());
  }
  else {
    kpt.Resize(3);
    kpt=0;
  }

  pos=startpos;
  if(!readsection(words, pos, latvectxt, "BOXSIZE"))
    error("BOXSIZE is required in HEG");
  if(latvectxt.size() != ndim)
    error("BOXSIZE must have exactly ",ndim, " values");

  vector< vector<string> > perturb_sec;
  vector <string> dumstring;
  pos=startpos;
  while(readsection(words, pos,dumstring, "PERTURB")) perturb_sec.push_back(dumstring);
  nperturb=perturb_sec.size();
  perturb_pos.Resize(nperturb, ndim);
  perturb_strength.Resize(nperturb);
  perturb_alpha.Resize(nperturb);
  perturb_spin.Resize(nperturb);
  for(int i=0; i< nperturb; i++) { 
    if(!readsection(perturb_sec[i], pos=0,dumstring,"POS")) error("Need POS in PERTURB");
    if(dumstring.size()!=ndim) error("wrong dimension in POS in PERTURB");
    for(int d=0; d< ndim; d++) perturb_pos(i,d)=atof(dumstring[d].c_str());
    if(!readvalue(perturb_sec[i], pos=0,perturb_strength(i),"STRENGTH")) error("Need STRENGTH in PERTURB");
    if(!readvalue(perturb_sec[i], pos=0,perturb_alpha(i),"SCALE")) error("Need SCALE in PERTURB");
    if(!readvalue(perturb_sec[i], pos=0,perturb_spin(i),"SPIN")) error("Need SPIN in PERTURB");
    
  }

  latVec.Resize(ndim, ndim);
  latVec=0.0;
  for(int i=0; i< ndim; i++)
    latVec(i,i)=atof(latvectxt[i].c_str());
  /*
  for(int i=0; i< ndim; i++) {
    for(int j=0; j< ndim; j++) {
      latVec(i,j)=atof(latvectxt[i*ndim+j].c_str());
    }
  }
  */

  origin.Resize(3);
  vector <string> origintxt;
  if(readsection(words, pos=0, origintxt, "ORIGIN")) {
    if(origintxt.size() < 3) error("ORIGIN section must have at least 3 elements.");
    for(int i=0; i< 3; i++) origin(i)=atof(origintxt[i].c_str());
  }
  else {
    origin=0;   //defaulting the origin to zero
  }

  vector <string> interactiontxt;
  // default is truncated Coulomb
  calcLocChoice=&HEG_system::calcLocTrunc;
  eeModel=2;
  same_spin_int=true;
  diff_spin_int=true;
  if(readsection(words, pos=0, interactiontxt, "INTERACTION")) {
    if( haskeyword(interactiontxt, pos=0, "EWALD") ) {
      calcLocChoice=&HEG_system::calcLocEwald;
      eeModel=1;
    }
    if( haskeyword(interactiontxt, pos=0, "TRUNCCOUL") ) {
      calcLocChoice=&HEG_system::calcLocTrunc;
      eeModel=2;
    }
    if( haskeyword(interactiontxt, pos=0, "GAUSS") ) {
      calcLocChoice=&HEG_system::calcLocGauss;
      eeModel=3;
      if(!readvalue(interactiontxt,pos=0,Gauss_a, "AMP"))
	error("Gauss interaction model requested, but no AMP given.");
      if(!readvalue(interactiontxt,pos=0,Gauss_s, "STDEV"))
	error("Gauss interaction model requested, but no STDEV given.");
      Gauss_s2=Gauss_s*Gauss_s;
    }
    if( haskeyword(interactiontxt, pos=0, "NO_SAME_SPIN") )
      same_spin_int=false;
    if( haskeyword(interactiontxt, pos=0, "NO_DIFF_SPIN") )
      diff_spin_int=false;
  }

  if ( eeModel==1 && ( !same_spin_int || !diff_spin_int ) )
    error("Spin-dependent Ewald interaction not supported.");

  

  //-------------cross products

  //cross product:  0->1x2, 1->2x0, 2->0x1
  Array2 <doublevar> crossProduct(ndim, ndim);
  crossProduct(0,0)=(latVec(1,1)*latVec(2,2)-latVec(1,2)*latVec(2,1));
  crossProduct(0,1)=(latVec(1,2)*latVec(2,0)-latVec(1,0)*latVec(2,2));
  crossProduct(0,2)=(latVec(1,0)*latVec(2,1)-latVec(1,1)*latVec(2,0));

  crossProduct(1,0)=(latVec(2,1)*latVec(0,2)-latVec(2,2)*latVec(0,1));
  crossProduct(1,1)=(latVec(2,2)*latVec(0,0)-latVec(2,0)*latVec(0,2));
  crossProduct(1,2)=(latVec(2,0)*latVec(0,1)-latVec(2,1)*latVec(0,0));

  crossProduct(2,0)=(latVec(0,1)*latVec(1,2)-latVec(0,2)*latVec(1,1));
  crossProduct(2,1)=(latVec(0,2)*latVec(1,0)-latVec(0,0)*latVec(1,2));
  crossProduct(2,2)=(latVec(0,0)*latVec(1,1)-latVec(0,1)*latVec(1,0));


  //------------ reciprocal cell  (used in setupEwald and showinfo)
  recipLatVec.Resize(ndim, ndim);
  doublevar det=Determinant(latVec, ndim);

  debug_write(cout, "cell volume ", det,"\n");
  cellVolume=det;

  for(int i=0; i< ndim; i++) {
    for(int j=0; j< ndim; j++) {
      recipLatVec(i,j)=crossProduct(i,j)/det;
    }
  }

  //------------ normal vectors (needed in enforcePbc)

  normVec.Resize(ndim, ndim);

  for(int i=0; i< ndim; i++) {
    for(int j=0; j < ndim; j++) {
      normVec(i,j)=crossProduct(i,j);
    }

    //Check to make sure the direction is facing out
    doublevar dotprod=0;
    for(int j=0; j < ndim; j++) {
      dotprod+=normVec(i,j)*latVec(i,j);
    }
    if(dotprod < 0) {
      for(int j=0; j< ndim; j++) {
        normVec(i,j)= -normVec(i,j);
      }
    }
  }

  corners.Resize(ndim, ndim);
  for(int i=0; i< ndim; i++) {
    for(int j=0; j< ndim; j++) {
      corners(i,j)=origin(j)+latVec(i,j);
    }
  }


  //------------ setup for interaction calculators (reciprocal part of Ewald sum etc.)
  if ( eeModel == 1 ) setupEwald(crossProduct);
  if ( eeModel == 2 ) setupTruncCoulomb();

  return 1;
}

//------------------------------------------------------------------------

int HEG_system::setupEwald(Array2 <doublevar> & crossProduct) {
  //----Set up reciprocal Ewald sum

  const int ndim=3;
  smallestheight=1e99;   //!< smallest distance that spans the cell

  for(int i=0; i< ndim; i++) {
    doublevar tempheight=0;
    doublevar length=0;
    for(int j=0; j< ndim; j++) {
      tempheight+=crossProduct(i,j)*latVec(i,j);
      length+=crossProduct(i,j)*crossProduct(i,j);
    }
    tempheight=fabs(tempheight)/sqrt(length);
    if(tempheight < smallestheight ) smallestheight=tempheight;
  }
  debug_write(cout, "elsu ", smallestheight,"\n");

  // We want only the simulation cell and its nearest neighbours
  // to contribute into the real-space part of Ewald sum; alpha
  // has to be chosen so that there are no contributions from distances
  // larger than smallestheight, i.e., we want
  //   erfc(alpha*smallestheight) \sim 0
  // note: erfc(5)=1.5e-12, erfc(6)=2.2e-17, erfc(6.5)=3.8e-20
  // ---
  // 5.0 chosen in PERIODIC system seems to be large enough; change
  // from 5.0 to 6.5 doubles the number of g-points 
  alpha=5.0/smallestheight;

  debug_write(cout, "alpha ", alpha, "\n");

  const int gmax=24;  //Could make this an option.
  const double qweight_cut=1e-18;

  ngpoints=(gmax+1)*(2*gmax+1)*(2*gmax+1)-gmax*(2*gmax+1)-gmax;
  Array2 <doublevar> gpointtemp(ngpoints,3);
  Array1 <doublevar> gweighttemp(ngpoints);
  int currgpt=0;
  int totgpt=0;
  int smallgpt=0;

  for(int ig=0; ig <= gmax; ig++) {
    int jgmin=-gmax;
    if(ig==0) jgmin=0;
    for(int jg=jgmin; jg <= gmax; jg++) {
      int kgmin=-gmax;
      if(ig==0 && jg==0) kgmin=0;
      for(int kg=kgmin; kg <= gmax; kg++) {
        totgpt++;
        for(int i=0; i< ndim; i++) {
          gpointtemp(currgpt, i)=2*pi*(ig*recipLatVec(0,i)
              +jg*recipLatVec(1,i)
              +kg*recipLatVec(2,i));
        }
        doublevar gsqrd=0;  // |g|^2
        for(int i=0; i< ndim; i++) {
          gsqrd+=gpointtemp(currgpt,i)*gpointtemp(currgpt,i);
        }

        if(gsqrd > 1e-8) {  // throw away (0,0,0)
          gweighttemp(currgpt)=4.0 * pi*exp(-gsqrd/(4*alpha*alpha))
            /(cellVolume*gsqrd);


          if(gweighttemp(currgpt) > qweight_cut) {
            currgpt++;
          }
        } else {
          smallgpt++;
        }
      }
    }
  }
  if ( ngpoints != totgpt ) {
    error("Wrong number of g-points while setting up Ewald summation.");
  }
  if ( smallgpt != 1 ) {
    error("More than 1 g-point detected as 0.");
  }
  if ( currgpt+smallgpt >= totgpt ) {
    error("Increase gmax in HEG_system.cpp, current value is too small.");
  }
  cout << "Reciprocal Ewald will use " << currgpt << " g-points out of "
       << totgpt << " examined." << endl;

  //Adjust to the correct number of kpoints..
  ngpoints=currgpt;
  gpoint.Resize(ngpoints, 3);
  gweight.Resize(ngpoints);

  for(int i=0; i< ngpoints; i++) {
    for(int d=0; d< ndim; d++) {
      gpoint(i,d)=gpointtemp(i,d);
    }
    gweight(i)=gweighttemp(i);
  }

  constEwald();

  return 1;

}

//------------------------------------------------------------------------

int HEG_system::setupTruncCoulomb() {

  // Setup for truncated Coulomb. It is correct only for cubic cell, I do not
  // see a reason to generalize it, since we are in a homogeneous system. Some
  // condition to test for cubic cell should be added, though.
  doublevar L=sqrt(latVec(0,0)*latVec(0,0)+latVec(0,1)*latVec(0,1)
		+latVec(0,2)*latVec(0,2));
  // backgr_trunc=1/V \int_{cell} 1/r
  // [ Mathematica helped with this one... ]
  doublevar backgr_trunc_single=-0.5*(pi+log(1351-780*sqrt(3.0)))/L;
  backgr_trunc=totnelectrons*totnelectrons/2.0*backgr_trunc_single;

  return 1;

}

//------------------------------------------------------------------------

int HEG_system::enforcePbc(Array1 <doublevar> & pos, Array1 <int> & nshifted) {
  assert(pos.GetDim(0) >=3);
  int shifted=0;
  nshifted.Resize(3);
  nshifted=0;

  int nshift=0;
  for(int i=0; i< 3; i++) {
    int shouldcontinue=1;
    while(shouldcontinue) {

      //Whether we're past the origin side
      doublevar tooshort=0;
      for(int j=0; j<3;j++) tooshort+=normVec(i,j)*(pos(j)-origin(j));

      //Whether we're past the lattice vector
      doublevar toofar=0;
      for(int j=0; j< 3; j++) toofar+=normVec(i,j)*(pos(j)-corners(i,j));

      //the 1e-12 seems to help avoid numerical problems, esp 
      //when integrating over a grid(which tends to hit the edges)
      if(tooshort < -1e-12) {
	//cout <<"tooshort " << tooshort << endl;
        for(int j=0; j< 3; j++) pos(j)+=latVec(i,j);
        shifted=1;
        nshifted(i)+=1;
        nshift++;
      }
      else if(toofar > 1e-12) {
	//cout << "toofar " << toofar << endl;
        for(int j=0; j< 3; j++) pos(j)-=latVec(i,j);
        shifted=1;
	// JK: here was +1, which works fine for real k-points that
	// correspond to standing waves. For general (complex) k-points the
	// wavefunction is "directional", however.
        nshifted(i)-=1;

        nshift++;
      }
      else {
        shouldcontinue=0;
      }
      if(nshift > 1000)
	
        error("Did over 1000 shifts and we're still out of the simulation cell."
	            "  There's probably something wrong.  Position : ", pos(i));
    }
  }

  return shifted;

}


//----------------------------------------------------------------------

doublevar HEG_system::calcLoc(Sample_point * sample)
{
  //cout << "Local energy" << endl;
  
  return (*this.*calcLocChoice)(sample)+calcLocPerturb(sample);
}
//---------
doublevar HEG_system::calcLocPerturb(Sample_point * sample) { 

  doublevar pot=0;
  int nd=ndim();
  for(int p=0; p < nperturb; p++) { 
    Array1 <doublevar> dist(nd), epos(nd);
    int s=perturb_spin(p);
    int estart, eend;
    doublevar r=0;
    if(s==0) { estart=0; eend=nspin(0); } 
    else { estart=nspin(0); eend=totnelectrons; } 
    for(int e=estart; e< eend; e++) { 
      sample->getElectronPos(e,epos);
      r=0;
      for(int d=0; d< nd; d++) { 
        doublevar del=epos(d)-perturb_pos(p,d);
        while(del > latVec(d,d)*0.5) del-=latVec(d,d);
        while(del < -latVec(d,d)*0.5) del+=latVec(d,d);
        r+=del*del;
      }
      pot+=perturb_strength(p)*exp(-perturb_alpha(p)*r);
    }
  }
  return pot;
}
//---------

doublevar HEG_system::calcLocEwald(Sample_point * sample)
{
  //cout << "Ewald local energy" << endl;
  doublevar ewalde=ewaldElectron(sample);
  return self_ee+ewalde+xc_correction;
}

//---------

doublevar HEG_system::calcLocTrunc(Sample_point * sample)
{
  //cout << "Truncated Coulomb local energy" << endl;
  
  sample->updateEEDist();
  
  doublevar ee_pot=0.0;
  Array1 <doublevar> eidist(5);

  if (diff_spin_int) {
    for(int e1=0; e1< nspin(0); e1++) {
      for(int e2 =nspin(0); e2 < totnelectrons; e2++) {
	sample->getEEDist(e1,e2, eidist);
	ee_pot+=1.0/eidist(0);
      }
    }
  }
  if (same_spin_int) {
    for(int e1=0; e1< nspin(0); e1++) {
      for(int e2 =e1+1; e2 < nspin(0); e2++) {
	sample->getEEDist(e1,e2, eidist);
	ee_pot+=1.0/eidist(0);
      }
    }
    for(int e1=nspin(0); e1< totnelectrons; e1++) {
      for(int e2 =e1+1; e2 < totnelectrons; e2++) {
	sample->getEEDist(e1,e2, eidist);
	ee_pot+=1.0/eidist(0);
      }
    }   
  }

  return ee_pot-backgr_trunc;
}

//---------

doublevar HEG_system::calcLocGauss(Sample_point * sample)
{
  //cout << "Gaussian toy local energy" << endl;
  
  sample->updateEEDist();

  doublevar ee_pot=0.0;
  Array1 <doublevar> eidist(5);

  if (diff_spin_int) {
    for(int e1=0; e1< nspin(0); e1++) {
      for(int e2 =nspin(0); e2 < totnelectrons; e2++) {
        sample->getEEDist(e1,e2, eidist);
        ee_pot+=exp(-eidist(1)/2.0/Gauss_s2);
      }
    }
  }
  if (same_spin_int) {
    for(int e1=0; e1< nspin(0); e1++) {
      for(int e2 =e1+1; e2 < nspin(0); e2++) {
        sample->getEEDist(e1,e2, eidist);
        ee_pot+=exp(-eidist(1)/2.0/Gauss_s2);
      }
    }
    for(int e1=nspin(0); e1< totnelectrons; e1++) {
      for(int e2 =e1+1; e2 < totnelectrons; e2++) {
        sample->getEEDist(e1,e2, eidist);
        ee_pot+=exp(-eidist(1)/2.0/Gauss_s2);
      }
    }   
  }


  return Gauss_a*ee_pot/Gauss_s/sqrt(2*pi);
}


//----------------------------------------------------------------------

void HEG_system::constEwald() {

  // that extra term not included in PERIODIC system;
  // contribution of the tripple loop seems to be essentially zero as
  // it should be due to our choice of alpha
  doublevar xi=0;
  int nlatvec=1;
  for (int i=-nlatvec; i <= nlatvec; i++) {
    for (int j=-nlatvec; j <= nlatvec; j++) {
      for (int k=-nlatvec; k <= nlatvec; k++) {
	if ( i==0 && j==0 && k==0 ) break;
	doublevar pos;
	doublevar pos2=0;
	for (int d=0; d < 3; d++) {
	  pos=i*latVec(0,d)+j*latVec(1,d)+k*latVec(2,d);
	  pos2+=pos*pos;
	}
	doublevar dxi=erfcm(alpha*sqrt(pos2))/sqrt(pos2);
	xi+=dxi;
      }
    }
  }
  // half of the self-interaction belongs to the simulation cell and
  // half to the given periodic image
  xi/=2;  
  xi-=alpha/sqrt(pi);
	
  cout << "correction in constEwald" << endl;
  cout << xi << " instead of " << -alpha/sqrt(pi) << endl;

  self_ee=-0.5*totnelectrons*totnelectrons*pi/(cellVolume*alpha*alpha)
    +totnelectrons*xi;

  //Correct for the exchange-correlation false interaction with
  //the uniform electron gas number.  cexc was fitted to several
  //systems..
  doublevar rs=(3./(4*pi*totnelectrons/cellVolume));
  rs=pow(rs, 1.0/3.0);
  doublevar cexc=0.36;
  xc_correction=cexc/rs;
  // let's make this zero for a moment
  xc_correction=0;

}


//----------------------------------------------------------------------

doublevar HEG_system::ewaldElectron(Sample_point * sample) {
  sample->updateEEDist();

  const int nlatvec=1;
  Array1 <doublevar> r1(3), r2(3);
  Array1 <doublevar> eidist(5);

  //-------------Electron-electron real-space part (pairs only,
  //   no self-interaction)         

  doublevar elecElec_real=0;

  for(int e1=0; e1< totnelectrons; e1++) {
    for(int e2 =e1+1; e2 < totnelectrons; e2++) {
      sample->getEEDist(e1,e2, eidist);
      for(int d=0; d< 3; d++) r1(d)=eidist(d+2);

      //----over  lattice vectors
      for(int kk=-nlatvec; kk <=nlatvec; kk++) {
        for(int jj=-nlatvec; jj <=nlatvec; jj++) {
          for(int ii=-nlatvec; ii <=nlatvec; ii++) {
            for(int d=0; d< 3; d++) {
              r2(d)=r1(d)+kk*latVec(0,d)+jj*latVec(1,d)+ii*latVec(2,d);
            }
            doublevar r=sqrt(r2(0)*r2(0)+r2(1)*r2(1)+r2(2)*r2(2));

            elecElec_real+=erfcm(alpha*r)/r;
          }
        }
      }
      //----done lattice vectors
    }
  }

  //---------electron reciprocal part (this DOES include the self-interaction)

  doublevar rdotg;
  Array2 <doublevar> elecpos(totnelectrons, 3);
  sample->getAllElectronPos(elecpos);
  doublevar elecElec_recip=0;
  for(int gpt=0; gpt < ngpoints; gpt++) {
    doublevar sum_sin=0, sum_cos=0;
    for(int e=0; e< totnelectrons; e++) {
      rdotg=gpoint(gpt, 0)*elecpos(e,0)
            +gpoint(gpt, 1)*elecpos(e,1)
            +gpoint(gpt, 2)*elecpos(e,2);

      sum_sin+=sin(rdotg);
      sum_cos+=cos(rdotg);

    }
    // NOTE: 1/2 from \sum_{e != e'} is cancelled with the fact that
    // we use only one g-point from g,-g pair
    elecElec_recip+=(sum_cos*sum_cos + sum_sin*sum_sin)*gweight(gpt);
  }

  //cout << "elecElec_real " << elecElec_real << endl;
  //cout << "elecElec_recip " << elecElec_recip << endl;

  return elecElec_real + elecElec_recip;
}


//------------------------------------------------------------------------
