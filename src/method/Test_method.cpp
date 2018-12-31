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

#include "Test_method.h"
#include "qmc_io.h"
#include "Program_options.h"
#include "System.h"
#include "Basis_function.h"
void Test_method::read(vector <string> words,
                       unsigned int & pos,
                       Program_options & options)
{


  //--Set up variables
  
  if(readvalue(words, pos=0, readconfig, "READCONFIG")) 
    canonical_filename(readconfig);
  
  sysprop=NULL;
  allocate(options.systemtext[0], sysprop);

  nelectrons=sysprop->nelectrons(0)+sysprop->nelectrons(1);


  wfdata=NULL;
  if(options.twftext.size() > 0) {
    cout << "wfdata allocate" << endl;
    allocate(options.twftext[0], sysprop, wfdata);
  }

  sysprop->generatePseudo(options.pseudotext, psp);

  vector <string> btxt;
  basis=NULL;
  if(readsection(words, pos=0, btxt, "BASIS")) {
    allocate(btxt,basis);
  }

  vector <string> backtxt;
  test_backflow=0;
  if(readsection(words,pos=0, backtxt,"BACKFLOW_TEST")) { 
    test_backflow=1;
    Array1 <Array1 <int> >  occupation(1);
    occupation(0).Resize(sysprop->nelectrons(0));
    for(int i=0; i< sysprop->nelectrons(0); i++) { 
      occupation(0)(i)=i;
    }
    backflow.readOrbitals(sysprop,backtxt);
    backflow.init(sysprop,occupation,backtxt);
  }
  
  if(haskeyword(words, pos=0, "PLOT_EE_CUSP"))
    plot_cusp=1;
  else
    plot_cusp=0;
  
  vector <string> derstxt;
  if(readsection(words, pos=0, derstxt, "TEST_PARMS_DERS")){
    parms_ders=1;
    if(haskeyword(derstxt, pos=0, "HESSIAN"))
      testhessian=1;
    else
      testhessian=0;
  }
  else{
    parms_ders=0;
    testhessian=0;
  }

  cout << "done setup " << endl;
}

int Test_method::showinfo(ostream & os)
{
  os << "#############Testing#################\n";
  wfdata->showinfo(os);
  //sysprop->showinfo(os);

  return 1;
}


#include "Split_sample.h"
#include "Guiding_function.h"


void check_numbers(doublevar num1,doublevar num2, ostream & os,
		   doublevar tol=1e-2) {
  string ok="OK";
  
  //make sure you are not divinding by zero
  //doublevar tol2=tol*tol;
  doublevar norm=num2;
  if(fabs(norm)<tol){
    if(fabs(num1)<tol)
      norm=1.0;
    else
      norm=num1;
  }
   
  if(fabs((num1-num2)/norm) > tol) ok="FAILED";
  os << num1 << "    " << num2 <<  "   " << ok << endl;
}

/*!

*/
void Test_method::run(Program_options & options, ostream & output)
{

  cout << "in run " << endl;

  cout.precision(15);
  
  
  Wavefunction * mywf=NULL;
  Sample_point * sample=NULL;
  
  wfdata->generateWavefunction(mywf);
  sysprop->generateSample(sample);

  /*
  if(readconfig!="") {
    ifstream checkfile(readconfig.c_str());
    if(!checkfile) 
      error("Couldn't open ", readconfig);
    string dummy;
    while(checkfile >> dummy) {
      if(read_config(dummy, checkfile, sample)) break;
    }
  }
  */

  Array1 <Config_save_point> config_pos;
  config_pos.Resize(0);
  if(readconfig!="") { 
    read_configurations(readconfig, config_pos);
    config_pos(0).restorePos(sample);
  }
  else sample->randomGuess();
    

  
  sample->attachObserver(mywf);
  Array1 <Wf_return> first_calc(nelectrons);
  cout <<"**********updateLap " << endl;
  mywf->updateLap(sample);
  cout << "******done " << endl;
  string indent="";
  for(int i=0; i< nelectrons; i++) {
    first_calc(i).Resize(mywf->nfunc(),5);
    mywf->getLap(i, first_calc(i));
    //first_calc(i).write(indent,cout);
  }
  
  cout << "#############checking derivatives" << endl;
  const doublevar del=1e-6;
  Array1 <doublevar> epos(3);
  Array1 <doublevar> new_epos(3);
  Wf_return test_wf(mywf->nfunc(), 5);
  for(int e=0; e < nelectrons; e++) {
    cout << "#######################\n";
    cout << "####electron " << e << " #####" << endl;
    cout << "#######################\n";
    sample->getElectronPos(e,epos);
    //cout <<" position "<<epos(0)<<" "<<epos(1)<<" "<<epos(2)<<endl;
    doublevar lap=0;
    doublevar phaselap=0;
    for(int d=0;d < 3; d++) {
      new_epos=epos;
      new_epos(d)+=del;
      sample->setElectronPos(e, new_epos);
      //mywf->notify(all_electrons_move,0);
      mywf->updateLap(sample);
      mywf->getLap(e,test_wf);
      doublevar ratio=exp(test_wf.amp(0,0)-first_calc(0).amp(0,0));
      doublevar derivative=(ratio-1)/del;
      doublevar phasederivative=(exp(test_wf.phase(0,0)-first_calc(0).phase(0,0))-1)/del;
      lap+=(test_wf.amp(0,d+1)*ratio-first_calc(e).amp(0,d+1))/del;
      phaselap+=(test_wf.phase(0,d+1)-first_calc(e).phase(0,d+1))/del;
      
      cout << "amplitude ";
      check_numbers(derivative,first_calc(e).amp(0,d+1),cout);
      cout << "phase ";
      check_numbers(phasederivative,first_calc(e).phase(0,d+1),cout);
    }
    cout << "amplitude laplacian ";
    check_numbers(lap,first_calc(e).amp(0,4),cout);
    cout << "phase laplacian ";
    check_numbers(phaselap,first_calc(e).phase(0,4),cout);
    sample->setElectronPos(e,epos);
    mywf->updateVal(sample);
  }
  
  

  if(basis!=NULL) { 
    for(double r=0; r < 10.0; r+=.05) { 
      Array1 <double> rscan(5); rscan=0.0;
      rscan(0)=r; rscan(1)=r*r; rscan(2)=r;
      Array2 <doublevar> vals(basis->nfunc(),5);
      basis->calcLap(rscan,vals);
      cout << "basisplot " << r << "  " << vals(0,0)  << " " << vals(0,4) << endl;
    }
    
    Array2 <doublevar> finite_der(basis->nfunc(),3,0.0);
    Array3 <doublevar> finite_hessian(basis->nfunc(),3,3,0.0);
    Array2 <doublevar> hessian(basis->nfunc(),10);
    Array2 <doublevar> tmp_hess(basis->nfunc(),10);
    sample->updateEIDist();
    Array1 <doublevar> dist(5);
    sample->getEIDist(0,0,dist);
    basis->calcHessian(dist,hessian);
    Array1 <doublevar> epos(3);
    Array1 <doublevar> npos(3);
    sample->getElectronPos(0,epos);
    for(int d=0; d< 3; d++) { 
      npos=epos;
      npos(d)+=del;
      sample->setElectronPos(0,npos);
      sample->updateEIDist();      
      sample->getEIDist(0,0,dist);
      basis->calcHessian(dist,tmp_hess);
      for(int f=0; f< basis->nfunc(); f++) { 
	  
        finite_der(f,d)=(tmp_hess(f,0)-hessian(f,0))/del;
        for(int d1=0; d1 < 3; d1++) { 
          finite_hessian(f,d,d1)=(tmp_hess(f,d1+1)-hessian(f,d1+1))/del;
        }
      }
    }

    for(int f=0; f< basis->nfunc(); f++) { 
      cout << "val " << hessian(f,0) << endl;    
      cout << "analytic    finite diff\n";
      check_numbers(hessian(f,1),finite_der(f,0),cout);
      check_numbers(hessian(f,2),finite_der(f,1),cout);
      check_numbers(hessian(f,3),finite_der(f,2),cout);
      check_numbers(hessian(f,4),finite_hessian(f,0,0), cout);
      check_numbers(hessian(f,5),finite_hessian(f,1,1), cout);
      check_numbers(hessian(f,6),finite_hessian(f,2,2), cout);
      check_numbers(hessian(f,7),finite_hessian(f,0,1), cout);
      check_numbers(hessian(f,8),finite_hessian(f,0,2), cout);
      check_numbers(hessian(f,9),finite_hessian(f,1,2), cout);
    }
      
  }

  

  if(plot_cusp){
    plotCusp(mywf, sample);
  }
  
  if(parms_ders){
    testParmDeriv(mywf, sample);
  }

  

  delete mywf; mywf=NULL;
  delete sample;
  sample=NULL;
  if(basis) delete basis;
}



void Test_method::testParmDeriv(Wavefunction * mywf, Sample_point * sample){
  int nparms=wfdata->nparms();
  cout <<"#######################################################\n";
  cout <<" Checking ParmDeriv for "<<nparms<<" parameters"<<endl;
  cout <<"#######################################################\n";
  doublevar del=1e-8;
  
  if(!wfdata->supports(parameter_derivatives))
    error("Wavefunction needs to supports analytic parameter derivatives");
  

  Array1 <Wf_return> base_wfval(nelectrons),test_wfval(nelectrons);
  Parm_deriv_return wfders;
  if(testhessian)
    wfders.need_hessian=1;
  else
    wfders.need_hessian=0;
  Parm_deriv_return wfders_tmp;
  wfders_tmp.need_hessian=0;
  mywf->updateLap(sample);
  for(int e=0; e< nelectrons; e++) { 
    base_wfval(e).Resize(1,5);
    test_wfval(e).Resize(1,5);
    mywf->getLap(e, base_wfval(e));
  }
  cout << "parmderiv" << endl;
  mywf->getParmDeriv(wfdata, sample,  wfders);
  cout << "done " << endl;
  Array1 <doublevar> Psi(nparms); 
  Array1 <doublevar> temp_parms(nparms),orig_parms(nparms);
  wfdata->getVarParms(temp_parms);
  for (int i=0;i<nparms;i++)
    orig_parms(i)=temp_parms(i);

  cout << "finite diff   analytic\n";
  for (int i=0;i<nparms;i++){
    temp_parms(i)+=del;
    //set new parameters
    wfdata->setVarParms(temp_parms);
    mywf->updateLap(sample);
    mywf->notify(all_electrons_move,0);
    for(int e=0; e< nelectrons; e++) 
      mywf->getLap(e, test_wfval(e));
    mywf->getParmDeriv(wfdata, sample,  wfders_tmp);
    doublevar psi=base_wfval(0).sign(0)*test_wfval(0).sign(0)
      *exp(test_wfval(0).amp(0,0)-base_wfval(0).amp(0,0));
    Psi(i)=psi;
    doublevar derivative=(psi-1)/del;
    cout <<" derivative check for parameter "<< i<< endl;
    check_numbers(derivative,wfders.gradient(i),cout);

    cout << "Gradient and laplacian check: size " << wfders.gradderiv.GetDim(0)
      << " " << wfders.gradderiv.GetDim(1) << " " << wfders.gradderiv.GetDim(2) << endl;
    for(int e=0; e< nelectrons; e++) { 
      for(int d=1; d< 5; d++) { 
        doublevar der=(test_wfval(e).amp(0,d)-base_wfval(e).amp(0,d))/del;
        check_numbers(der,wfders.gradderiv(i,e,d-1),cout);
      }
    }
    
    if(testhessian){
      cout <<" hessian check1"<<endl;
      for(int j=i;j<nparms;j++){
        doublevar hessian_wf=(psi*wfders_tmp.gradient(j)-wfders.gradient(j))/del;
        check_numbers(hessian_wf,wfders.hessian(i,j),cout);
      }
    }
    temp_parms(i)-=del;
  }//end of loop over i
  

}




//----------------------------------------------------------------------

void Test_method::plotCusp(Wavefunction * mywf, Sample_point * sample){
  
  Array1 < Array1 <doublevar> >orig_walker(nelectrons);
  for(int e=0; e < nelectrons; e++){
    orig_walker(e).Resize(3);
    sample->getElectronPos(e,orig_walker(e));
  }
  ofstream cusp("cusp_r_12.dat");
  int nsteps=500;
  cusp << "#### Cusp plot: r_12  , phi(r_12), 1/2(lap(r_1)+lap(r_2)) -1/r_12" <<endl;
  const doublevar ndel=1e-7;
  Wf_return test_wf(mywf->nfunc(), 5);
  Wf_return test_wf1(mywf->nfunc(), 5);
  Array1 <doublevar> new_pos(3), old_pos(3);
  sample->getElectronPos(1,old_pos);
  new_pos=old_pos;
  new_pos(0)+=nsteps*ndel;
  old_pos(0)-=nsteps*ndel;
  sample->setElectronPos(0,new_pos);
  sample->setElectronPos(1,old_pos);
  for(int f=0;f<nsteps-1;f++){
    sample->getElectronPos(0,new_pos);
    new_pos(0)-=ndel;
    sample->getElectronPos(1,old_pos);
    old_pos(0)+=ndel;
    sample->setElectronPos(0,new_pos);
    sample->setElectronPos(1,old_pos);
    mywf->updateLap(sample);
    mywf->getLap(0,test_wf);
    mywf->getLap(1,test_wf1);
    doublevar distance=new_pos(0)-old_pos(0);
    if(fabs(distance)>0)
      cusp <<  distance << "  "<<test_wf.sign(0)*exp(test_wf.amp(0,0))
	//<<"  "<<test_wf.amp(0,4)<<"  "<<test_wf1.amp(0,4)<<"  "<<1/fabs(distance) 
	   <<"  "<<0.5*(test_wf.amp(0,4)+test_wf1.amp(0,4))-1.0/fabs(distance)<<endl;
  }
  cusp <<endl;
  cusp.close();
  
  for(int e=0; e < nelectrons; e++){
    sample->setElectronPos(e,orig_walker(e));
  }
  sample->getElectronPos(0,old_pos);
  ofstream cusp2("cusp_r_1N.dat");
  cusp2 << "#### Cusp plot: r_1N  , phi(r_1N), 1/2(lap(r_1)+lap(r_N)) -1/r_1N" <<endl;
  new_pos=old_pos;
  new_pos(0)+=nsteps*ndel;
  old_pos(0)-=nsteps*ndel;
  sample->setElectronPos(0,new_pos);
  sample->setElectronPos(nelectrons-1,old_pos);
  for(int f=0;f<nsteps-1;f++){
    sample->getElectronPos(0,new_pos);
    new_pos(0)-=ndel;
    sample->getElectronPos(nelectrons-1,old_pos);
    old_pos(0)+=ndel;
    sample->setElectronPos(0,new_pos);
    sample->setElectronPos(nelectrons-1,old_pos);
    mywf->updateLap(sample);
    mywf->getLap(0,test_wf);
    mywf->getLap(nelectrons-1,test_wf1);
    doublevar distance=new_pos(0)-old_pos(0);
    if(fabs(distance)>0)
      cusp2 <<  distance << "  "<<test_wf.sign(0)*exp(test_wf.amp(0,0))
	//<<"  "<<test_wf.amp(0,4)<<"  "<<test_wf1.amp(0,4)<<"  "<<1/fabs(distance) 
	    <<"  "<<0.5*(test_wf.amp(0,4)+test_wf1.amp(0,4))-1.0/fabs(distance)<<endl;
  }
  cusp2 <<endl;
  cusp2.close();
  for(int e=0; e < nelectrons; e++){
    sample->setElectronPos(e,orig_walker(e));
  }
  
}

