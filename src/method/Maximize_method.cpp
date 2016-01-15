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


#include "Maximize_method.h"
#include "qmc_io.h"
#include "ulec.h"
#include <iomanip>
#include "Program_options.h"
#include "System.h"
#include "Split_sample.h"
#include "Properties.h"
#include "Generate_sample.h"

//----------------------------------------------------------------------

void Maximize_method::read(vector <string> words,
            unsigned int & pos,
            Program_options & options) { 
  allocate(options.systemtext[0],  sys);
  sys->generatePseudo(options.pseudotext, pseudo);

  debug_write(cout, "wfdata allocate\n");
  wfdata=NULL;
  if(options.twftext.size() < 1) error("Need TRIALFUNC section for OPTIMIZE");
  allocate(options.twftext[0], sys, wfdata);

  if(!readvalue(words,pos=0,nconfig,"NCONFIG"))
    nconfig=100;
  
}

//----------------------------------------------------------------------

void Maximize_method::run(Program_options & options, ostream & output) { 
  Wavefunction * wf=NULL;
  Sample_point * sample=NULL;
  sys->generateSample(sample);
  wfdata->generateWavefunction(wf);
  sample->attachObserver(wf);
  
  Wf_return wfret(1,2);

  Array1 <Config_save_point> config_pos(nconfig);
  Array1 <doublevar> epos(3);
  
  Primary guidewf;
  generate_sample(sample,wf,wfdata,&guidewf,nconfig,config_pos);
  int nelectrons=sample->electronSize();
  Properties_gather mygather;

  string tablename=options.runid+".table";
  ofstream tableout(tablename.c_str());

  for(int e=0; e< nelectrons; e++) {
    sample->getElectronPos(e,epos);
    for(int d=0; d< 3; d++) {
      tableout << "e" << e << "_"<< d << " ";
    }
  }
  tableout << "elocal" << " ";
  tableout << "logpsi" << endl;
  
  for(int i=0; i < nconfig; i++) { 
    maximize(sample,wf,config_pos(i));
    for(int e=0; e< nelectrons; e++) {
      sample->getElectronPos(e,epos);
      for(int d=0; d< 3; d++) {
        tableout << epos(d) << " ";
      }
    }
    Properties_point pt;
    mygather.gatherData(pt, pseudo, sys, wfdata, wf, 
                            sample, &guidewf);
    
    tableout << pt.energy(0) << " ";
    wf->getVal(wfdata,0,wfret);
    tableout << wfret.amp(0,0) << endl;
     
    config_pos(i).write(cout);
  }

  delete wf;
  delete sample;
}
//----------------------------------------------------------------------

#include "macopt.h"
class Maximize_fn:public Macopt {
public:
  Maximize_fn(int _n,int _verbose, double _tol,
              int _itmax,int _rich):Macopt(_n,_verbose,_tol,_itmax,_rich) { 
              } 
  ~Maximize_fn(){ };


  //-----------------------------------------
  void update_positions(double * _p) { 
    int count=1; 
    Array1 <doublevar> epos(3);
    int nelec=sample->electronSize();
    for(int e=0; e< nelec; e++) { 
      for(int d=0; d< 3; d++) { 
        epos(d)=_p[count++];
      }
      sample->setElectronPos(e,epos);
    }
    
  }
  //-----------------------------------------

  double func(double * _p) { 
    update_positions(_p);
    wf->updateVal(wfdata,sample);
    Wf_return wfret(1,2);
    wf->getVal(wfdata,0,wfret);
    return -wfret.amp(0,0);
  }
  //-----------------------------------------
  
  double dfunc(double * _p, double * _g) { 
    update_positions(_p);
    wf->updateLap(wfdata,sample);
    int count=0;
    int nelec=sample->electronSize();
    Wf_return lap(1,5);
    wf->getVal(wfdata,0,lap);
    _g[count++]=-lap.amp(0,0);
    for(int e=0; e< nelec; e++) { 
      wf->getLap(wfdata,e,lap);
      for(int d=0; d< 3; d++) { 
        _g[count++]=-lap.amp(0,d+1);
      }
    }
    return _g[0];
  }
  //-----------------------------------------
  
  void iteration_print(double f, double gg, double tol,  int itn) {
    cout << "It " << endl;
  }
  //-----------------------------------------

  Wavefunction * wf;
  Wavefunction_data * wfdata;
  Sample_point * sample;
  System * system;
};

//----------------------------------------------------------------------


void Maximize_method::maximize(Sample_point * sample,Wavefunction * wf,Config_save_point & pt) { 
  int nelectrons=sample->electronSize();
  pt.restorePos(sample);
  Maximize_fn maximizer(nelectrons*3,1,1e-12,1000,1);
  maximizer.wf=wf;
  maximizer.wfdata=wfdata;
  maximizer.sample=sample;
  maximizer.system=sys;

  Array1 <doublevar> allpos(nelectrons*3+1);
  Array1 <doublevar> epos(3);
  int count=1;
  for(int e=0; e< nelectrons; e++) {
    sample->getElectronPos(e,epos);
    for(int d=0; d< 3; d++) { 
      allpos(count++)=epos(d);
    }
  }

  maximizer.maccheckgrad(allpos.v,nelectrons*3,0.001,nelectrons*3);
  maximizer.macoptII(allpos.v,nelectrons*3);  
  


  pt.savePos(sample);
}

//----------------------------------------------------------------------
  
int Maximize_method::showinfo(ostream & os) { 
  return 1;
}
//----------------------------------------------------------------------

