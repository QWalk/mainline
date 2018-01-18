#include "Average_ldots.h"
#include "ulec.h"

void Average_ldots::read(System * sys, Wavefunction_data * wfdata,
                      vector <string> & words){  
}

void Average_ldots::read(vector <string> & words){
}

void Average_ldots::write_init(string & indent, ostream & os){
  os << indent << "L dot S\n";
}

void Average_ldots::evaluate(Wavefunction_data * wfdata, Wavefunction * wf, System * sys,
                          Sample_point * sample, Average_return & avg) {
  avg.type="ldots";
  int ndim=3;
 // int nwf=wf->nfunc();
  avg.vals.Resize(1);
  avg.vals=0.0;
  int nelectrons=sample->electronSize();
  int nions=sample->ionSize();
  Array1 <doublevar> e_pos(ndim);    
  Array1 <doublevar> ion_pos(ndim);
  //get the gradient of wavefunction
  Parm_deriv_return derivatives;
  derivatives.need_hessian = 0;
  wf->updateVal(wfdata, sample);
  if(!wfdata->supports(parameter_derivatives))
    error("Wavefunction needs to supports analytic parameter derivatives");
 // if(!wf->getParmDeriv(wfdata, sample,derivatives)) {
 //   error("WF needs to support parmderivs for now.");
//  }
  wf->getParmDeriv(wfdata, sample, derivatives);  
  for(int e=0; e< nelectrons; e++) {
    int spin = 1;
    if(e >= sys->nelectrons(0)) spin = -1;
    sample->getElectronPos(e,e_pos);
    cout << "get here!!" << endl;
    for(int at=0; at < nions; at++) {
      sample->getIonPos(at,ion_pos);
      avg.vals(0) += spin*(e_pos(0)-ion_pos(0))*derivatives.val_gradient(e,1) - spin*(e_pos(1)-ion_pos(1))*derivatives.val_gradient(e,0); 
    }  
  }    
}

void Average_ldots::write_summary(Average_return & avg, Average_return & err,
                               ostream & os) {
  int ndim=avg.vals.GetDim(0);
  assert(ndim <= err.vals.GetDim(0));
  os << "Expectation value of L dot S \n";

  for(int i=0;i<avg.vals.GetDim(0);i++)
    os << avg.vals(i) << " +/- " << err.vals(i) << endl;

}

void Average_ldots::jsonOutput(Average_return & avg, Average_return & err,
                           ostream & os){

  os << "\"" << avg.type << "\:{" << endl;

  for(int i=0;i<avg.vals.GetDim(0);i++){
    os << "\"value\":[" << "  " <<  avg.vals(i) << " ]," << endl;
    os << "\"error\":[" << "  " << err.vals(i) << " ]," << endl;
  }
}
