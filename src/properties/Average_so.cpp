#include "Average_so.h"
#include "Pseudopotential_so.h"
#include "ulec.h"

void Average_so::read(System * sys, Wavefunction_data * wfdata, 
                      vector <string> & words){
  unsigned int pos=0;
  vector <string> psp_file;
  psp_file.resize(0);
  vector <string> words_pspfile;
  words_pspfile.reserve(50);
  vector <string> pseudosec;
  vector < vector <string> > psp_txt;
  while(readsection(words,pos,psp_file,"READPSP_SO")!=0){
    //ifstream psp_soFile(psp_file[0]);
    ifstream psp_soFile(psp_file[0].c_str());
    parsefile(psp_soFile,words_pspfile);
  }
  pos=0;
  while( readsection(words_pspfile, pos, pseudosec, "PSEUDO") != 0) {
    psp_txt.push_back(pseudosec);
  } 
  psp_so = new Pseudopotential_so;
  psp_so->read(psp_txt,sys);

  init_grid.Resize(2);
  del_grid.Resize(2);  // \delta_\theta, \delta_\phi
  num_grid.Resize(2);  // num_\theta, num_\phi
  
  if(!readvalue(words, pos=0, init_grid(0), "INIT_THETA")) init_grid(0) = 0;
  if(!readvalue(words, pos=0, init_grid(1), "INIT_PHI")) init_grid(1) = -3.14159;

  if(!readvalue(words, pos=0, del_grid(0), "DELTA_THETA")) del_grid(0) = 0;  
  if(!readvalue(words, pos=0, del_grid(1), "DELTA_PHI")) del_grid(1) = 0;  
  
  if(!readvalue(words, pos=0, num_grid(0), "NUM_THETA")) num_grid(0) = 1;  
  if(!readvalue(words, pos=0, num_grid(1), "NUM_PHI")) num_grid(1) = 1;  

}

void Average_so::read(vector <string> & words){
}

void Average_so::write_init(string & indent, ostream & os){
  os << indent << "psp_so\n";
  os << del_grid(0) << " " << del_grid(1) << " " << num_grid(0) << " " << num_grid(1) << "\n";
 // psp_so->showinfo(os);
}

/*
void Average_so::evaluate(Wavefunction_data * wfdata, Wavefunction * wf, System * sys,
                          Sample_point * sample, Average_return & avg) {
  avg.type="psp_so";
  int nwf=wf->nfunc();
  avg.vals.Resize(nwf);
  avg.vals=0.0;
  int nelectrons=sample->electronSize();
 
  Array2 <doublevar> totalv;
  totalv.Resize(nelectrons,nwf);

  psp_so->calcNonlocSeparated(wfdata,sys,sample,wf,totalv);
  for(int w=0;w<nwf;w++){
    for(int e=0;e<nelectrons;e++){
      avg.vals(w)+=totalv(e,w);
    }
  } 
 
}
*/
void Average_so::evaluate(Wavefunction_data * wfdata, Wavefunction * wf, System * sys,
                          Sample_point * sample, Average_return & avg) {
  avg.type="psp_so";
  
  
  int num_dir=num_grid(0)*num_grid(1);
  //cout << "total number of spin directions" << num_dir << endl;  

  int nwf=wf->nfunc();
  avg.vals.Resize(nwf*num_dir); // how many soi energies do we compute
  avg.vals=0.0;
  int nelectrons=sample->electronSize();
 
  Array2 <doublevar> totalv(nwf,3); //nwf=1
  //totalv.Resize(nwf,3);
  
  int nrandvar=psp_so->nTest();
  Array1 <doublevar> rand_num(nrandvar);
  for(int i=0; i< nrandvar; i++)
    rand_num(i)=rng.ulec();
 
  psp_so->calcNonlocWithTest(wfdata,sys,sample,wf,rand_num,totalv);

  for (int i=0; i<num_grid(0); i++){
    doublevar spin_theta = init_grid(0)+i*del_grid(0);
    for (int j=0; j< num_grid(1); j++){
      doublevar spin_phi = init_grid(1)+j*del_grid(1);
      for (int iwf=0; iwf<nwf; iwf++){
        avg.vals(iwf*nwf+i*num_grid(1)+j) = totalv(iwf,0)*sin(spin_theta)*cos(spin_phi)
                                        +totalv(iwf,1)*sin(spin_theta)*sin(spin_phi)
	      		                +totalv(iwf,2)*cos(spin_theta);
      }
    }
  }

}

void Average_so::randomize(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample) {
  psp_so->randomize();
}

void Average_so::write_summary(Average_return & avg, Average_return & err, 
                               ostream & os) {
  psp_so->showinfo(os);
  os << "First order energy correction due to spin-orbit interaction\n";
 
  for(int i=0;i<avg.vals.GetDim(0);i++)
    os << avg.vals(i) << " +/- " << err.vals(i) << endl;

 /*
  for (int i=0; i<num_grid(0); i++){
    doublevar spin_theta = 0.0+i*del_grid(0);
    for (int j=0; j< num_grid(1); j++){
      doublevar spin_phi = 0.0+j*del_grid(1);
      for (int iwf=0; iwf<1; iwf++){
        doublevar spin_x = sin(spin_theta)*cos(spin_phi);
        doublevar spin_y = sin(spin_theta)*sin(spin_phi);
        doublevar spin_z = cos(spin_theta);
        os << "spin direction" << spin_theta << "  " << spin_phi <<  "   ";
        os << spin_x << "  " << spin_y << "  " << spin_z << endl; 
        os << avg.vals(iwf*1+i*num_grid(0)+j) << " +/- " << err.vals(i) << endl;
      }
    }
  }

 */

  
}

void Average_so::jsonOutput(Average_return & avg, Average_return & err,
                           ostream & os){
  
  os << "\"" << avg.type << "\:{" << endl;

  for(int i=0;i<avg.vals.GetDim(0);i++){
    os << "\"value\":[" << "  " <<  avg.vals(i) << " ]," << endl;
    os << "\"error\":[" << "  " << err.vals(i) << " ]," << endl; 
  }
}



