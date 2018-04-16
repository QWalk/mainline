#include "Average_ldots.h"
#include "ulec.h"

void Average_ldots::read(System * sys, Wavefunction_data * wfdata,
                      vector <string> & words){  
  unsigned int pos=0;
//  Array1 <doublevar> direction(3);
  direction.Resize(3);
  if(!readvalue(words, pos=0, direction(0), "ZAXIS_X")) direction(0) = 0;
  if(!readvalue(words, pos=0, direction(1), "ZAXIS_Y")) direction(1) = 0;
  if(!readvalue(words, pos=0, direction(2), "ZAXIS_Z")) direction(2) = 1;
  
  doublevar dir_amp;
  dir_amp=direction(0)*direction(0)+direction(1)*direction(1)+direction(2)*direction(2);
  dir_amp=sqrt(dir_amp);
  for(int i=0;i<3;i++) direction(i) = direction(i)/dir_amp;

}

void Average_ldots::read(vector <string> & words){
}

void Average_ldots::write_init(string & indent, ostream & os){
  os << indent << "LdotS\n";
}

void Average_ldots::evaluate(Wavefunction_data * wfdata, Wavefunction * wf, System * sys,
                          Sample_point * sample, Average_return & avg) {
  avg.type="ldots";
  int ndim=3;
//  int nwf=wf->nfunc();
  avg.vals.Resize(2);
  int nelectrons=sample->electronSize();
  avg.vals=0.0;
  int nions=sample->ionSize();
  Array1 <doublevar> e_pos(ndim);    
  Array1 <doublevar> ion_pos(ndim);
  
  wf->updateLap(wfdata, sample);
  Wf_return temp_lap;
  temp_lap.Resize(1,5);

  Array1 <doublevar> ls_real(3);
  Array1 <doublevar> ls_imag(3);
  ls_real = 0.0;
  ls_imag = 0.0;
 // for(int i=0; i< nwf; i++){
  for(int e=0; e< nelectrons; e++) {
    int spin = 1;
    if(e >= sys->nelectrons(0)) spin = -1;
    sample->getElectronPos(e,e_pos);
    wf->getLap(wfdata, e, temp_lap);
    for(int at=0; at < nions; at++) {
      sample->getIonPos(at,ion_pos);
      ls_real(0) += (e_pos(2)-ion_pos(2))*temp_lap.amp(0,2)*spin - (e_pos(1)-ion_pos(1))*temp_lap.amp(0,3)*spin;
      ls_imag(0) += (e_pos(2)-ion_pos(2))*temp_lap.phase(0,2)*spin - (e_pos(1)-ion_pos(1))*temp_lap.phase(0,3)*spin; //Lx\cdot S
      ls_real(1) += (e_pos(2)-ion_pos(2))*temp_lap.amp(0,1)*spin - (e_pos(0)-ion_pos(0))*temp_lap.amp(0,3)*spin;
      ls_imag(1) += (e_pos(2)-ion_pos(2))*temp_lap.phase(0,1)*spin - (e_pos(0)-ion_pos(0))*temp_lap.phase(0,3)*spin; //Ly \cdot S
      ls_real(2) += (e_pos(0)-ion_pos(0))*temp_lap.amp(0,2)*spin - (e_pos(1)-ion_pos(1))*temp_lap.amp(0,1)*spin;
      ls_imag(2) += (e_pos(0)-ion_pos(0))*temp_lap.phase(0,2)*spin - (e_pos(1)-ion_pos(1))*temp_lap.phase(0,1)*spin; //Lz \cdot S
  
    }  
  }
  avg.vals(0) = ls_real(0)*direction(0)+ls_real(1)*direction(1)+ls_real(2)*direction(2);
  avg.vals(1) = ls_imag(0)*direction(0)+ls_imag(1)*direction(1)+ls_imag(2)*direction(2);

  //}

}

void Average_ldots::write_summary(Average_return & avg, Average_return & err,
                               ostream & os) {
  int ndim=avg.vals.GetDim(0);
  assert(ndim <= err.vals.GetDim(0));
  os << "Expectation value of L dot S\n"; 
//  os << "Expectation value of L dot S along " << "(" << direction(0) << " " << direction(1) << " " << direction(2) << ")\n";

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
