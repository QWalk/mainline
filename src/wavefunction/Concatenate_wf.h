#include "Wavefunction.h"
#include <vector>

class Wavefunction_data;
class Concatenate_wf_data;

class Concatenate_wf_storage: public Wavefunction_storage{
  Concatenate_wf_storage(int nwfs){
     wf_stores.resize(nwfs);
  }
  ~Concatenate_wf_storage(){
    for(int i=0;i<wf_stores.size();i++){
      if(wf_stores[i]){
        delete wf_stores[i];
      }
    }
  }

private:
  friend class Concatenate_wf;
  vector<Wavefunction_storage *> wf_stores;
};

class Concatenate_wf: public Wavefunction {
public:
  Concatenate_wf(vector<Wavefunction*> wf_in):wfs(wf_in){}

  virtual ~Concatenate_wf(){
    for(int i=0;i<wfs.size();i++){
      if(wfs[i]!=NULL)
        delete wfs[i];
    }
  }
  
  void updateLap(Wavefunction_data * wfdata, Sample_point * point){
    Concatenate_wf_data * dataptr;
    recast(wfdata,dataptr);
    for(int i=0;i<wfs.size();i++){
      wfs[i]->updateLap(dataptr->wf_datas[i],point);  
    }
  }
  void getLap(Wavefunction_data * wfdata, int e, Wf_return & vals){
    Concatenate_wf_data * dataptr;
    recast(wfdata,dataptr);
    vals.Resize(wfs.size(),5);
    Wf_return tmp;
    tmp.Resize(1,5);
    for(int i=0;i<wfs.size();i++){
      wfs[i]->getLap(dataptr->wf_datas[i],e,tmp);
      for(int d=0;d<5;d++){
        vals.amp(i,d)=tmp.amp(0,d);
        vals.phase(i,d)=tmp.phase(0,d);
      }
    }
  }

  //Unimplemented and unused methods
  void notify(change_type change,int num){
    for(int i=0;i<wfs.size();i++)
      wfs[i]->notify(change,num);
  }

  void updateVal(Wavefunction_data *, Sample_point *){error("updateVal not impl in Concatenate");}
  void getVal(Wavefunction_data *, int, Wf_return &){error("getVal not impl in Concatenate");}
  void getSymmetricVal(Wavefunction_data *, int, Wf_return &){error("getSymmetricVal not impl in Concatenate");}
  
  void generateStorage(Wavefunction_storage * & wfstore){
    wfstore=new Concatenate_wf_storage(wfs.size());
    Concatenate_wf_storage * store;
    recast(wfstore,store);
    for(int i=0;i<wfs.size();i++)
      wfs[i]->generateStorage(store->wf_stores[i]);
  }

  void saveUpdate(Sample_point * sample, int e, Wavefunction_storage * wfstore){
    Concatenate_wf_storage * store;
    recast(wfstore,store);
    for(int i=0;i<wfs.size();i++)
      wfs[i]->saveUpdate(sample,e,store->wf_stores[i]);
  }

  void restoreUpdate(Sample_point * sample, int e, Wavefunction_storage * wfstore){
    Concatenate_wf_storage * store;
    recast(wfstore,store);
    for(int i=0;i<wfs.size();i++)
      wfs[i]->restoreUpdate(sample,e,store->wf_stores[i]);
  }

private:
  vector<Wavefunction*> wfs;

};


