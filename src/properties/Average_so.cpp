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
}

void Average_so::read(vector <string> & words){
}

void Average_so::write_init(string & indent, ostream & os){
  os << indent << "psp_so\n";
  psp_so->showinfo(os);
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
  int nwf=wf->nfunc();
  avg.vals.Resize(nwf);
  avg.vals=0.0;
  int nelectrons=sample->electronSize();
 
  Array1 <doublevar> totalv;
  totalv.Resize(nwf);

  int nrandvar=psp_so->nTest();
  Array1 <doublevar> rand_num(nrandvar);
  for(int i=0; i< nrandvar; i++)
    rand_num(i)=rng.ulec();

  psp_so->calcNonlocWithTest(wfdata,sys,sample,wf,rand_num,totalv);
  avg.vals=totalv;

}

void Average_so::randomize(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample) {
  psp_so->randomize();
}

void Average_so::write_summary(Average_return & avg, Average_return & err, 
                               ostream & os) {
  int ndim=avg.vals.GetDim(0);
  assert(ndim <= err.vals.GetDim(0));
  psp_so->showinfo(os);
  os << "First order energy correction due to spin-orbit interaction\n";

  for(int i=0;i<avg.vals.GetDim(0);i++)
    os << avg.vals(i) << " +/- " << err.vals(i) << endl;
  
}

void Average_so::jsonOutput(Average_return & avg, Average_return & err,
                           ostream & os){
  
  os << "\"" << avg.type << "\:{" << endl;

  for(int i=0;i<avg.vals.GetDim(0);i++){
    os << "\"value\":[" << "  " <<  avg.vals(i) << " ]," << endl;
    os << "\"error\":[" << "  " << err.vals(i) << " ]," << endl; 
  }
}



