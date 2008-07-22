//------------------------------------------------------------------------
//src/Backflow_wf_data.cpp

#include "Qmc_std.h"
#include "qmc_io.h"
#include "Backflow_wf_data.h"
#include "Wavefunction_data.h"
#include "Backflow_wf.h"
#include "Cslat_wf.h"
#include <algorithm>


void backflow_config(Sample_point * sample, 
		     int e,
		     Array3 <doublevar> & corr,
		     Array3 <doublevar> & onebody,
		     Array3 <doublevar> & threebody_diffspin,
		     Sample_point * temp_samp) {
  //cout << "bfconfig " << endl;
  int nelectrons=sample->electronSize();
  Array1 <doublevar> newpos(3);
  sample->updateEEDist();
  sample->updateEIDist();
  sample->getElectronPos(e,newpos);

  Array1 <doublevar> dist(5);
  for(int j=0; j< e; j++) { 
    sample->getEEDist(j,e,dist);
    for(int d=0; d< 3; d++) { 
      newpos(d)+=corr(j,e,0)*dist(d+2);
    }
  }
  for(int k=e+1; k< nelectrons; k++) { 
    sample->getEEDist(e,k,dist);
    for(int d=0; d< 3; d++) { 
      newpos(d)-=corr(e,k,0)*dist(d+2);
    }
  }

  //cout << "here" << e << endl;
  
  int natoms=sample->ionSize();
  for(int at=0; at< natoms; at++) { 
    sample->getEIDist(e,at,dist);
    for(int d=0; d< 3; d++) { 
      //start: added for ei back-flow
      newpos(d)+=(onebody(e,at,0)+threebody_diffspin(e,at,0))*dist(d+2);
      //end: added for ei back-flow
    }
  }

 
  //Array1 <doublevar> oldpos(3);
  //sample->getElectronPos(e,oldpos);
  //for(int d=0; d< 3; d++) { 
  //  newpos(d)+=onebody(e,0)*oldpos(d);
  //}

  temp_samp->setElectronPos(0,newpos);

  
}

//start: added for ei back-flow
void Backflow_wrapper::updateLapjastgroup(Sample_point * sample, int e, Array3 <doublevar> & threebody_diffspin){
  //cout <<"updateLapjastgroup:start"<<endl;
  int nelectrons=sample->electronSize();
  int natoms=sample->ionSize();
  int eibasis_max=jgroupdata.maxEIBasis();
  int eebasis_max=jgroupdata.nEEBasis();
  threebody_diffspin.Resize(nelectrons, natoms, 5);
  threebody_diffspin=0.0;

  if(has_electron_ion_bf){
    Array4 <doublevar> eibasis_save(nelectrons, natoms, eibasis_max ,5);
    Array3 <doublevar> eebasis(nelectrons, eebasis_max, 5);
    Array3 <doublevar> eibasis(natoms, eibasis_max ,5);
    for(int i=0; i< nelectrons; i++) {
      jgroupdata.updateEIBasis(i,sample,eibasis);
      for(int at=0; at< natoms; at++) {
	for(int j=0; j< eibasis_max; j++) {
	  for(int d=0; d< 5; d++) {
	    eibasis_save(i,at,j,d)=eibasis(at,j,d);
	  }
	}
      }
    }
    
    
    jgroupdata.updateEEBasis(e,sample,eebasis);
    jgroupdata.three_body_diffspin.updateLap_E_I(e,eibasis_save,eebasis,threebody_diffspin);
  }
  //cout <<"done"<<endl;
}

void Backflow_wrapper::updateValjastgroup(Sample_point * sample, int e, Array3 <doublevar> & threebody_diffspin){
  int nelectrons=sample->electronSize();
  int natoms=sample->ionSize();
  int eibasis_max=jgroupdata.maxEIBasis();
  int eebasis_max=jgroupdata.nEEBasis();

  threebody_diffspin.Resize(nelectrons, natoms, 5);
  threebody_diffspin=0.0;
  if(has_electron_ion_bf){
    Array4 <doublevar> eibasis_save(nelectrons, natoms, eibasis_max ,5);
    Array3 <doublevar> eebasis(nelectrons, eebasis_max, 5);
    Array3 <doublevar> eibasis(natoms, eibasis_max ,5);
    for(int i=0; i< nelectrons; i++) {
      jgroupdata.updateEIBasis(i,sample,eibasis);
      for(int at=0; at< natoms; at++) {
	for(int j=0; j< eibasis_max; j++) {
	  for(int d=0; d< 5; d++) {
	    eibasis_save(i,at,j,d)=eibasis(at,j,d);
	  }
	}
      }
    }
    
    
    jgroupdata.updateEEBasis(e,sample,eebasis);
    jgroupdata.three_body_diffspin.updateVal_E_I(e,eibasis_save,eebasis,threebody_diffspin);
  }
}
//end: added for ei back-flow


void Backflow_wrapper::updateVal(Sample_point * sample, 
				 Jastrow2_wf & jast,
				 int e,
				 int listnum, 
				 Array2 <doublevar> & newvals) { 

  int nelectrons=sample->electronSize();
  Array3<doublevar> jast_corr;
  jast.updateVal(&jdata,sample);
  jast.get_twobody(jast_corr);
  Array3 <doublevar> onebody;
  jast.get_onebody(onebody);
  //start: added for ei back-flow
  Array3 <doublevar> threebody_diffspin;
  updateValjastgroup(sample,e,threebody_diffspin);
  //end: added for ei back-flow
  backflow_config(sample,e,jast_corr,onebody,threebody_diffspin,temp_samp);
  temp_samp->updateEIDist();
  molecorb->updateVal(temp_samp,0,listnum,newvals);
}

void Backflow_wrapper::getNeighbors(Sample_point * sample,
				    Jastrow2_wf & jast, 
				    int e, Array1 <int> & list, 
				    int & nlist) {

  //needs to be changed in the future
  //does not check for neighbors of ei back-flow term: threebody_diffspin.
  //which in principle is using different ee, ei basis.
  nlist=0;
  int nelectrons=sample->electronSize();
  list.Resize(nelectrons);
  Array3 <doublevar> jast_corr;
  jast.updateLap(&jdata,sample);
  jast.get_twobody(jast_corr);
  doublevar threshold=1e-10;
  for(int j=0; j< 3; j++) { 
    if(fabs(jast_corr(j,e,0)) > threshold) { 
      list(nlist++)=j;
    }
  }
  list(nlist++)=e;
  for(int k=e+1; k< nelectrons; k++) { 
    if(fabs(jast_corr(e,k,0)) > threshold) { 
      list(nlist++)=k;
    }
  }
}


void Backflow_wrapper::updateLap(Sample_point * sample, 
				 Jastrow2_wf & jast,
				 int e, 
				 int listnum, 
				 Array2 <doublevar> & newvals, 
				 //!<(mo,[val,grad,hess])
				 Array3 <doublevar>& coor_deriv, 
				 //!< (i,alpha,beta)
				 Array2 <doublevar> & coor_laplacian 
				 //!< (i,alpha)
				 ) { 

  //cout << "Backflow_wrapper: updateLap " << endl;
  int nelectrons=sample->electronSize();
  int natoms=sample->ionSize();
  
  Array3 <doublevar> jast_corr;
  jast.updateLap(&jdata,sample);
  jast.get_twobody(jast_corr);
  Array3 <doublevar> onebody;
  jast.get_onebody(onebody);
  //start: added for ei back-flow
  Array3 <doublevar> threebody_diffspin;
  updateLapjastgroup(sample,e,threebody_diffspin);
  //end: added for ei back-flow
  backflow_config(sample,e,jast_corr,onebody,threebody_diffspin,temp_samp);

  //Array1 <doublevar> tmp_pos(3);
  //cout << "Backflow_wrapper::updateLap"<<endl;
  //temp_samp->getElectronPos(0,tmp_pos);
  //cout <<"xyz: "<<tmp_pos(0)<<",  "<<tmp_pos(1)<<",  "<<tmp_pos(2)<<endl;

  temp_samp->updateEIDist();
  molecorb->updateHessian(temp_samp,0,listnum,newvals);
 
  coor_deriv.Resize(nelectrons,3,3);
  coor_laplacian.Resize(nelectrons,3);
  coor_deriv=0;
  coor_laplacian=0;
  //self-derivative..
  
  for(int d=0; d< 3;d ++) 
    coor_deriv(e,d,d)=1;

  Array1 <doublevar> dist(5);
  
  for(int at=0; at < natoms; at++) { 
    sample->getEIDist(e,at,dist);
    for(int a=0;a < 3; a++) { 
      for(int b=0; b< 3; b++) { 	
	coor_deriv(e,a,b)+=onebody(e,at,a+1)*dist(b+2);
      }
      coor_deriv(e,a,a)+=onebody(e,at,0);
      coor_laplacian(e,a)+=2.0*onebody(e,at,a+1)
	                  +onebody(e,at,4)*dist(a+2);
    }
  }



  for(int k=0; k < e; k++) { 
    sample->getEEDist(k,e,dist);
    for(int a=0; a< 3; a++) {
      coor_deriv(e,a,a)-=jast_corr(k,e,0);
      for(int b=0; b< 3; b++) { 
	coor_deriv(e,a,b)+=jast_corr(e,k,a+1)*dist(b+2);
      }
      coor_laplacian(e,a)+=jast_corr(e,k,4)*dist(a+2)
	-2*jast_corr(e,k,a+1);
    }
  }

  for(int k=e+1; k < nelectrons; k++) { 
    sample->getEEDist(e,k,dist);
    for(int a=0; a< 3; a++) {
      coor_deriv(e,a,a)-=jast_corr(e,k,0);
      for(int b=0; b< 3; b++) { 
	coor_deriv(e,a,b)-=jast_corr(e,k,a+1)*dist(b+2);
      }
      coor_laplacian(e,a)-= jast_corr(e,k,4)*dist(a+2)
	+2*jast_corr(e,k,a+1);
    }
  }


  //cout << "middle " <<endl;
  //for(int d=0; d< 3; d++) { 
  //  cout << "lap " << coor_laplacian(e,d) << endl;
  //}
  
  /*
  cout << "jast corr " << endl;
  for(int i=0; i< nelectrons; i++) { 
    for(int j=0; j< nelectrons; j++) { 
      cout << jast_corr(i,j,4) << "   ";
    }
    cout << endl;
  }
  */


  for(int i=0; i< e; i++) { 
    sample->getEEDist(i,e,dist);
    for(int a=0; a< 3; a++) { 
      for(int b=0; b< 3; b++) { 
	coor_deriv(i,a,b)+=jast_corr(i,e,a+1)*dist(b+2);
      }
      coor_deriv(i,a,a)+=jast_corr(i,e,0);
      coor_laplacian(i,a)= jast_corr(i,e,4)*dist(a+2)
	+2*jast_corr(i,e,a+1);
    }
  }
  for(int j=e+1; j<nelectrons; j++) { 
    sample->getEEDist(e,j,dist);
    for(int a=0; a< 3; a++) { 
      for(int b=0; b< 3; b++) { 
	coor_deriv(j,a,b)-=jast_corr(j,e,a+1)*dist(b+2);
      }
      coor_deriv(j,a,a)+=jast_corr(e,j,0);
      coor_laplacian(j,a)= -jast_corr(j,e,4)*dist(a+2)+2*jast_corr(j,e,a+1);
    }
  }

  //start: added for ei back-flow
  //extra three body stuff 
  //diagonal contribution
  for(int at=0; at < natoms; at++) { 
    sample->getEIDist(e,at,dist);
    for(int a=0;a < 3; a++) { 
      for(int b=0; b< 3; b++) { 	
	coor_deriv(e,a,b)+=threebody_diffspin(e,at,a+1)*dist(b+2);
      }
      coor_deriv(e,a,a)+=threebody_diffspin(e,at,0); 
      coor_laplacian(e,a)+=2.0*threebody_diffspin(e,at,a+1)+threebody_diffspin(e,at,4)*dist(a+2);
    }
  }

  //off-diagonal contribution
  for(int at=0; at < natoms; at++) {
    sample->getEIDist(e,at,dist);
    for(int i=0; i< e; i++) {
      for(int a=0;a < 3; a++) { 
	for(int b=0; b< 3; b++) { 	
	  coor_deriv(i,a,b)+=threebody_diffspin(i,at,a+1)*dist(b+2);
	}
	coor_laplacian(i,a)+=threebody_diffspin(i,at,4)*dist(a+2);
      }
    }
    for(int j=e+1; j<nelectrons; j++) { 
      for(int a=0;a < 3; a++) { 
	for(int b=0; b< 3; b++) { 	
	  coor_deriv(j,a,b)+=threebody_diffspin(j,at,a+1)*dist(b+2);
	}
	coor_laplacian(j,a)+=threebody_diffspin(j,at,4)*dist(a+2);
      }
    }
  }
  //end: added for ei back-flow
  

  //for(int d=0; d< 3; d++) { 
  //  cout << "lap " << coor_laplacian(e,d) << endl;
  //}
  //cout << "done " << endl;

}

//----------------------------------------------------------------------

void Backflow_wrapper::readOrbitals(System * sys, 
			    vector <string> & words) {
  unsigned int pos=0;
  vector <string> mowords;
  if(!readsection(words, pos=0, mowords, "ORBITALS"))
    error("Need ORBITALS section");
  allocate(mowords,sys, molecorb);
}

//----------------------------------------------------------------------

void Backflow_wrapper::init(System * sys, 
			    Array1 <Array1 <int> > & totoccupation,
			    vector <string> & words) {

  //cout <<"Backflow_wrapper::init"<<endl;
  unsigned int pos=0;
  //vector <string> mowords;
  //if(!readsection(words, pos=0, mowords, "ORBITALS"))
  //  error("Need ORBITALS section");
  //allocate(mowords,sys, molecorb);

  molecorb->buildLists(totoccupation);
  sys->generateSample(temp_samp);

  vector <string> jwords;
  if(!readsection(words,pos=0,jwords, "EE_BF"))
    error("Need EE_BF section");
  pos=0;
  jdata.read(jwords,pos,sys);
  //start: added for ei back-flow
  vector <string> jgroupwords;
  if(readsection(words,pos=0,jgroupwords, "EI_BF")){
    has_electron_ion_bf=1;
    jgroupdata.set_up(jgroupwords, sys);
    if(jgroupdata.hasThreeBodySpin()==0)
      error("EI_BF works only with Three body spin jastrow group!");
  }
  else
    has_electron_ion_bf=0;
  //end: added for ei back-flow
  //cout <<"end of: Backflow_wrapper::init"<<endl;
}

//----------------------------------------------------------------------

void Backflow_wrapper::showinfo(ostream & os) { 
  os << "Molecular Orbital object : ";
  molecorb->showinfo(os);
  os << "\n";

  os<< "Electron-electron backflow : "<<endl;
  jdata.showinfo(os);
  //start: added for ei back-flow
  if(has_electron_ion_bf){
    os<< "Electron-ion backflow : "<<endl;
    string indent="  ";
    jgroupdata.showinfo(indent,os);
  }
  //end: added for ei back-flow
}


//----------------------------------------------------------------------

void Backflow_wrapper::writeinput(string & indent, ostream & os) { 
  os << indent << "ORBITALS {\n";
  string indent2=indent+"  ";
  molecorb->writeinput(indent2, os);
  os << indent << "}\n";

  os << indent << "EE_BF { \n";
  jdata.writeinput(indent2,os);
  os << indent << "}" << endl;
  //start: added for ei back-flow
  if(has_electron_ion_bf){
    os << indent << "EI_BF { \n";
    jgroupdata.writeinput(indent2,os);
    os << indent << "}" << endl;
  }
  //end: added for ei back-flow
}

//######################################################################
//----------------------------------------------------------------------
//######################################################################

void Determinant_keeper::read(vector <string> & words, System * sys) {

  unsigned int pos=0;
  optimize_det=haskeyword(words, pos=0, "OPTIMIZE_DET");

  vector <string> strdetwt;
  readsection(words, pos=0, strdetwt, "DETWT");
  ndet=strdetwt.size();
  detwt.Resize(ndet);
  for(int det=0; det < ndet; det++)
    detwt(det)=atof(strdetwt[det].c_str());

  vector <string>  statevec;
  if(!readsection(words, pos=0, statevec, "STATES"))
    error("Need STATES in Backflow wave function");


  nelectrons.Resize(2);
  nelectrons(0)=sys->nelectrons(0);
  nelectrons(1)=sys->nelectrons(1);

  //pos=startpos;
  unsigned int canonstates=ndet*(nelectrons(0)+nelectrons(1));
  if( canonstates != statevec.size())  {
    error("in STATES section, expecting to find ", canonstates,
	  " states(as calculated from NSPIN), but found ",
	  statevec.size(), " instead.");
  }
  

  int tote=nelectrons(0)+nelectrons(1);

  occupation_orig.Resize(ndet, 2);
  for(int det=0; det < ndet; det++) {
    for(int s=0; s<2; s++) {
      occupation_orig(det,s).Resize(nelectrons(s));
    }
  }

  int counter=0;
  for(int det=0; det<ndet; det++) {
    for(int s=0; s<2; s++) {
      for(int e=0; e<nelectrons(s); e++) {
	occupation_orig(det,s)(e)=atoi(statevec[counter].c_str())-1;

	counter++;
      }
    }
  }

  //Fill up calculation helpers

  spin.Resize(tote);
  opspin.Resize(tote);
  rede.Resize(tote);

  int eup=nelectrons(0);
  for(int e=0; e<eup; e++) {
    spin(e)=0;
    opspin(e)=1;
    rede(e)=e;
  }
  for(int e=eup; e<tote; e++) {
    spin(e)=1;
    opspin(e)=0;
    rede(e)=e-eup;
  }
	
}

//----------------------------------------------------------------------

void Determinant_keeper::getOccupation(Array1 <Array1 <int> > & totoccupation,
				       Array2 <Array1 <int> > & occupation) { 
  occupation.Resize(ndet, 2);

  for(int det=0; det < ndet; det++) {
    for(int s=0; s<2; s++) {
        //cout << det << " "  << s << endl;
        occupation(det,s).Resize(nelectrons(s));
    }
  }
  




  //Find what MO's are necessary for each spin

  totoccupation.Resize(2);
  for(int s=0; s<2; s++)  {
    vector <int> totocctemp;
    for(int det=0; det<ndet; det++)  {
      for(int mo=0; mo < nelectrons(s); mo++) {
	int place=-1;
	int ntot=totocctemp.size();
	for(int i=0; i< ntot; i++) {
	  if(occupation_orig(det,s)(mo)==totocctemp[i]) {
	    place=i;
	    break;
	  }
	}
	if(place==-1) { //if we didn't find the MO
	  occupation(det,s)(mo)=totocctemp.size();
	  totocctemp.push_back(occupation_orig(det,s)(mo));
	}
	else {
	  occupation(det,s)(mo)=place;
	}
      }
    }
    
    //cout << "done assignment " << endl;

    totoccupation(s).Resize(totocctemp.size());
    for(int i=0; i<totoccupation(s).GetDim(0); i++)
    {
      totoccupation(s)(i) = totocctemp[i];
      //cout << "total occupation for " << s<< " : "
      // << totoccupation(s)(i) << endl;
    }
  }
  
}


//----------------------------------------------------------------------

void Determinant_keeper::showinfo(ostream & os ) { 
  if(ndet > 1)
    os << ndet << " Determinants\n";
  else
    os << "Determinant" << endl;

  int nfunc=1;
  for(int f=0; f< nfunc; f++)
  {
    if(nfunc > 1)
      os << "For function " << f << endl;
    for(int det=0; det<ndet; det++)
      {
	if(ndet > 1) {
	  os << "Determinant " << det << ":\n";
	  os << "Weight: " << detwt(det) << endl;
	}

	os << "State: \n";
	for(int s=0; s<2; s++)
	  {
	    if(s==0)
	      os << "spin up:\n";
	    if(s==1)
	      os << "spin down: \n";

	    os << "  ";
	    for(int e=0; e<nelectrons(s); e++)
	      {
		os << occupation_orig(det,s)(e)+1 << " ";
		if((e+1)%10 == 0)
		  os << endl << "  ";
	      }
	    os << endl;
	  }
      }
  }
}

//----------------------------------------------------------------------

void Determinant_keeper::writeinput(string & indent, ostream & os ) { 
  os << indent << "DETWT { ";
  for(int det=0; det < ndet; det++) {
    os << detwt(det) << "  ";
  }
  os << "}" << endl;

  os << indent << "STATES { " << endl << indent <<"  ";
  for(int det=0; det < ndet; det++) {
    for(int s=0; s<2; s++) {
      for(int e=0; e< nelectrons(s); e++) {
	os << occupation_orig(det,s)(e)+1 << " ";
	if((e+1)%20 ==0)
	  os << endl << indent << "  ";
      }
      os << endl << indent << "  ";
    }
  }
  os << "}" << endl;

}

//######################################################################
//----------------------------------------------------------------------
//######################################################################

/*!
*/
void Backflow_wf_data::read(vector <string> & words, unsigned int & pos,
                        System * sys)
{
  unsigned int startpos=pos;
  //optimize_backflow=haskeyword(words, pos=startpos,"OPTIMIZE_BACKFLOW");
  vector <string> detwords;
  if(!readsection(words, pos=0, detwords,"DETERMINANT"))
    error("Couldn't find DETERMINANT section");
  vector <string> bfwords;
  if(!readsection(words,pos=0, bfwords,"BFLOW"))
    error("Couldn't find BACKFLOW section");

  bfwrapper.readOrbitals(sys,bfwords);
  dkeeper.read(detwords,sys);
  Array1 <Array1 <int> > totoccupation;
  dkeeper.getOccupation(totoccupation, occupation);
  bfwrapper.init(sys, totoccupation, bfwords);

}

//----------------------------------------------------------------------

int Backflow_wf_data::supports(wf_support_type support) {
  switch(support) {
  case laplacian_update:
    return 1;
  case density:
    return 0;
  case parameter_derivatives:
    return 0;
  default:
    return 0;
  }
}

//----------------------------------------------------------------------



void Backflow_wf_data::generateWavefunction(Wavefunction *& wf)
{
  assert(wf==NULL);

  wf=new Backflow_wf;
  Backflow_wf * slatwf;
  recast(wf, slatwf);
  slatwf->init(this);
  attachObserver(slatwf);
}

int Backflow_wf_data::showinfo(ostream & os)
{
  os << "Backflow-Slater" << endl;
  dkeeper.showinfo(os);
  bfwrapper.showinfo(os);
  return 1;
}

//----------------------------------------------------------------------

int Backflow_wf_data::writeinput(string & indent, ostream & os)
{


  os << indent << "BACKFLOW" << endl;
  os << indent << "BFLOW { " << endl;
  string indent2=indent+"  ";
  bfwrapper.writeinput(indent2,os);
  os << indent << "}" << endl;
  os << indent << "DETERMINANT { " << endl;
  dkeeper.writeinput(indent2,os);
  os << indent << "}" << endl;

  return 1;
}

//------------------------------------------------------------------------
void Backflow_wf_data::getVarParms(Array1 <doublevar> & parms)
{
  //cout <<"start getVarParms"<<endl;
  
  /*
  if(optimize_backflow) {
    error("need to do parms for backflow");
  }
  else if(optimize_det) {
    parms.Resize(detwt.GetDim(0)-1);
    for(int i=1; i< detwt.GetDim(0); i++) {
      parms(i-1)=detwt(i);
    }
  }
  else {
    parms.Resize(0);
  }
  */
  //error("need to do parameter optimization");
  //cout <<"done getVarParms"<<endl;
  
  parms.Resize(nparms());
  Array1 <doublevar> parms_tmp1;
  Array1 <doublevar> parms_tmp2;
  bfwrapper.jdata.getVarParms(parms_tmp1);
  //start: added for ei back-flow
  if(bfwrapper.has_electron_ion_bf){
    bfwrapper.jgroupdata.getVarParms(parms_tmp2);
  }
  //end: added for ei back-flow
  int sizetotal=bfwrapper.nparms();
  int sizejdata=bfwrapper.jdata.nparms();
  parms.Resize(sizetotal);
  //start: added for ei back-flow
  for(int i=0;i<sizetotal;i++){
    if(i<sizejdata)
      parms(i)=parms_tmp1(i);
    else
      parms(i)=parms_tmp2(i-sizejdata);
  }
  //end: added for ei back-flow
  //cout <<"done getVarParms"<<endl;
}

void Backflow_wf_data::setVarParms(Array1 <doublevar> & parms)
{
  //cout <<"start setVarParms"<<endl;
  assert(parms.GetDim(0)==nparms());
  /*
  if(optimize_backflow) {

  }
  else if(optimize_det) {
    for(int i=1; i< detwt.GetDim(0); i++) {
      detwt(i)=parms(i-1);
    }
  }
  else {
    parms.Resize(0);
  }
  */
  //error("need to do parameter optimization");
  Array1 <doublevar> parms_tmp1(bfwrapper.jdata.nparms());
  //start: added for ei back-flow
  Array1 <doublevar> parms_tmp2(bfwrapper.jgroupdata.nparms());
  //end: added for ei back-flow

  int sizetotal=bfwrapper.nparms();
  int sizejdata=bfwrapper.jdata.nparms();
  for(int i=0;i<sizetotal;i++){
    if(i<sizejdata)
      parms_tmp1(i)=parms(i);
    else
      parms_tmp2(i-sizejdata)=parms(i);
  }
  bfwrapper.jdata.setVarParms(parms_tmp1);
  //start: added for ei back-flow
  if(bfwrapper.has_electron_ion_bf){
    bfwrapper.jgroupdata.setVarParms(parms_tmp2);
  }
  //end: added for ei back-flow
  int max=wfObserver.size();
  //cout << "slatmax " << max << endl;
  for(int i=0; i< max; i++) {
    wfObserver[i]->notify(all_wf_parms_change, 0);
  }
  //cout <<"done setVarParms"<<endl;
}
