#include "Wavefunction_data.h"
#include <vector>

class Concatenate_wf;

class Concatenate_wf_data: public Wavefunction_data {
public:
  Concatenate_wf_data(vector<Wavefunction_data*> wf_data_in):wf_datas(wf_data_in){}

  virtual ~Concatenate_wf_data(){
    for(int i=0;i<wf_datas.size();i++){
      if(wf_datas[i]!=NULL)
        delete wf_datas[i];
    }
  }

  int supports(wf_support_type support){
    for(int i=0;i<wf_datas.size();i++){
      if(!wf_datas[i]->supports(support))
        return 0;
    }
    return 1;
  }
 
  //Unimplemented and unused methods
  void read(vector <string> & words, unsigned int & pos, System * sys){error("read not impl in Concatenate");}
  void getVarParms(Array1 <doublevar>&){error("getVarParms not impl in Concatenate");}
  void setVarParms(Array1 <doublevar>&){error("setVarParms not impl in Concatenate");}
  int nparms(){error("nparms not impl in concatenate");return 0;}
  int valSize(){error("valSize not impl in concatenate");return 0;}
  int writeinput(string &, ostream &){error("writeinput not imp in cocatenate");return 0;}
  void generateWavefunction(Wavefunction * &){error("generateWaveFunction not impl in cocatenate");}
  Concatenate_wf_data * clone() const{error("clonse not impl in concatenate");return NULL;}

private:
  vector <Wavefunction_data *> wf_datas;
  friend class Concatenate_wf;
};

