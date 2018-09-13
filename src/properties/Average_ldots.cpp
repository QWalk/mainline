#include "Average_ldots.h"
#include "ulec.h"

void Average_ldots::read(System * sys, Wavefunction_data * wfdata,
                      vector <string> & words){  
  unsigned int pos=0;
//  Array1 <doublevar> direction(3);
  init_grid.Resize(2);
  del_grid.Resize(2);  // \delta_\theta, \delta_\phi
  num_grid.Resize(2);  // num_\theta, num_\phi
  
  int nions=sys->nIons();

  if(!readvalue(words, pos=0, init_grid(0), "INIT_THETA")) init_grid(0) = 0;
  if(!readvalue(words, pos=0, init_grid(1), "INIT_PHI")) init_grid(1) = -3.14159;

  if(!readvalue(words, pos=0, del_grid(0), "DELTA_THETA")) del_grid(0) = 0;
  if(!readvalue(words, pos=0, del_grid(1), "DELTA_PHI")) del_grid(1) = 0;

  if(!readvalue(words, pos=0, num_grid(0), "NUM_THETA")) num_grid(0) = 1;
  if(!readvalue(words, pos=0, num_grid(1), "NUM_PHI")) num_grid(1) = 1;
  
  if(!readvalue(words, pos=0, at_i, "ATOM_I")) at_i = 1;
  if(!readvalue(words, pos=0, at_f, "ATOM_F")) at_f = nions;
  
}

void Average_ldots::read(vector <string> & words){
}

void Average_ldots::write_init(string & indent, ostream & os){
  os << indent << "LdotS\n";
  os << indent << "from atom " << at_i << " to atom " << at_f << "\n";
}

void Average_ldots::evaluate(Wavefunction_data * wfdata, Wavefunction * wf, System * sys,
                          Sample_point * sample, Average_return & avg) {
  avg.type="ldots";
  int ndim=3;
  int num_dir=num_grid(0)*num_grid(1);
//  int nwf=wf->nfunc();
  avg.vals.Resize(6*num_dir);
  int nelectrons=sample->electronSize();
  avg.vals=0.0;
  //int nions=sample->ionSize();
  Array1 <doublevar> e_pos(ndim);    
  Array1 <doublevar> ion_pos(ndim);
  
  wf->updateLap(wfdata, sample);
  Wf_return temp_lap;
  temp_lap.Resize(1,5);

  Array1 <doublevar> ls_real(3);
  Array1 <doublevar> ls_imag(3);

  Array1 <doublevar> direction(3);

  for (int i=0; i<num_grid(0); i++){
    doublevar spin_theta = init_grid(0)+i*del_grid(0);
    for (int j=0; j< num_grid(1); j++){
      doublevar spin_phi = init_grid(1)+j*del_grid(1);
      ls_real = 0.0;
      ls_imag = 0.0;
      for(int e=0; e< nelectrons; e++) {
        int spin = 1;
        if(e >= sys->nelectrons(0)) spin = -1;
        sample->getElectronPos(e,e_pos);
        wf->getLap(wfdata, e, temp_lap);
        //for(int at=0; at < nions; at++) {
        for(int at=at_i-1; at < at_f; at++) {
          sample->getIonPos(at,ion_pos);
          ls_real(0) += (e_pos(1)-ion_pos(1))*temp_lap.amp(0,3)*spin - (e_pos(2)-ion_pos(2))*temp_lap.amp(0,2)*spin;
          ls_imag(0) += (e_pos(1)-ion_pos(1))*temp_lap.phase(0,3)*spin - (e_pos(2)-ion_pos(2))*temp_lap.phase(0,2)*spin; //Lx\cdot S, yz-zy
          ls_real(1) += (e_pos(2)-ion_pos(2))*temp_lap.amp(0,1)*spin - (e_pos(0)-ion_pos(0))*temp_lap.amp(0,3)*spin;
          ls_imag(1) += (e_pos(2)-ion_pos(2))*temp_lap.phase(0,1)*spin - (e_pos(0)-ion_pos(0))*temp_lap.phase(0,3)*spin; //Ly \cdot S, zx-xz??
          ls_real(2) += (e_pos(0)-ion_pos(0))*temp_lap.amp(0,2)*spin - (e_pos(1)-ion_pos(1))*temp_lap.amp(0,1)*spin;
          ls_imag(2) += (e_pos(0)-ion_pos(0))*temp_lap.phase(0,2)*spin - (e_pos(1)-ion_pos(1))*temp_lap.phase(0,1)*spin; //Lz \cdot S, xy-yx
  
        }  
      }
      direction(0) = sin(spin_theta)*cos(spin_phi);
      direction(1) = sin(spin_theta)*sin(spin_phi);
      direction(2) = cos(spin_theta);

      avg.vals(6*i*num_grid(1)+6*j) = ls_real(0)*direction(0);//+ls_real(1)*direction(1)+ls_real(2)*direction(2);
      avg.vals(6*i*num_grid(1)+6*j+1) = ls_imag(0)*direction(0);//+ls_imag(1)*direction(1)+ls_imag(2)*direction(2);
      avg.vals(6*i*num_grid(1)+6*j+2) = ls_real(1)*direction(1);
      avg.vals(6*i*num_grid(1)+6*j+3) = ls_imag(1)*direction(1);
 
      avg.vals(6*i*num_grid(1)+6*j+4) = ls_real(2)*direction(2);
      avg.vals(6*i*num_grid(1)+6*j+5) = ls_imag(2)*direction(2);
    }
  }

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
