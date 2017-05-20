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

#include "System.h"
#include "Molecular_system.h"
#include "Periodic_system.h"
#include "ulec.h"
#include "SHO_system.h"
#include "HEG_system.h"
#include "qmc_io.h"

int allocate(vector <string> & syswords,
             System * & sysptr)
{
  assert(sysptr==NULL);
  if(caseless_eq(syswords[0],"MOLECULE")
     || caseless_eq(syswords[0],"MOLECULAR")) 
    sysptr=new Molecular_system;
  
  else if(caseless_eq(syswords[0],"PERIODIC"))
    sysptr=new Periodic_system;
  else if(caseless_eq(syswords[0],"SHO"))
    sysptr=new SHO_system;
  else if(caseless_eq(syswords[0],"HEG"))
    sysptr=new HEG_system;
  else
    error("Couldn't understand ", syswords[0], " in system section.");
  
  unsigned int pos=0;
  sysptr->read(syswords, pos);
  return 1;
}

//---------------------------------------------------------------------

void System::calcKinetic(Wavefunction_data * wfdata,
                         Sample_point * sample,
                         Wavefunction * wf,
                         Array1 <doublevar> & lap)
{
  Array2<doublevar> Kin;
  calcKineticSeparated(wfdata, sample, wf,Kin);
  int nelectrons = sample->electronSize();
  int nwf = wf->nfunc(); 
  for (int w=0; w<nwf; w++) {
    lap(w) = 0.0; 
    for(int e=0; e<nelectrons; e++) {
      lap(w) += Kin(e, w); 
    }
  }
  //cout << "Calculating kinetic energy \n";
  /*
  assert(lap.GetDim(0)>= wf->nfunc());
  int nelectrons=sample->electronSize();
  int nwf=wf->nfunc();
  lap=0;
  Wf_return temp(nwf,5);
  for(int w=0; w< nwf; w++)
  {
    for(int e=0; e< nelectrons; e++)
    {

      wf->getLap(wfdata, e, temp);
      lap(w)+=temp.amp(w,4);
      if ( temp.is_complex==1 ) {
        lap(w)-=(  temp.phase(w,1)*temp.phase(w,1)
            +temp.phase(w,2)*temp.phase(w,2)
            +temp.phase(w,3)*temp.phase(w,3) );
      }
      //cout << "total laplacian: " << lap(0) <<  " amp  " 
      //  << temp.amp(w,4) << endl;
        
    }
    lap(w)*=-0.5;
  }
  */
  //cout << "laplacian " << lap(0) << endl;
  //cout << "Calculating kinetic energy done \n";
}


//----------------------------------------------------------------------

/*
   Calculate the kinetic energy for all electrons
   into an array
   */
void System::calcKineticSeparated(Wavefunction_data * wfdata,
    Sample_point * sample,
    Wavefunction * wf, Array2<doublevar> & Kin)
{
  //cout << "Calculating kinetic energy \n";

  //  assert(lap.GetDim(0)>= wf->nfunc());
  int nelectrons=sample->electronSize();
  int nwf=wf->nfunc();
  wf->updateLap(wfdata,sample);
  Wf_return temp(nwf,5);
  Kin.Resize(nelectrons, nwf);

  for(int w=0; w< nwf; w++) {
    for(int e=0; e< nelectrons; e++) {
      wf->getLap(wfdata, e, temp);
      //  lap(e, w)+=temp.amp(w,4);
      Kin(e, w) = temp.amp(w, 4); 
      if ( temp.is_complex==1 ) {
        // lap(e, w)-=(  temp.phase(w,1)*temp.phase(w,1)
        //   +temp.phase(w,2)*temp.phase(w,2)
        //   +temp.phase(w,3)*temp.phase(w,3) );
        Kin(e, w)-=(  temp.phase(w,1)*temp.phase(w,1)
            +temp.phase(w,2)*temp.phase(w,2)
            +temp.phase(w,3)*temp.phase(w,3) );
      }
      //cout << "total laplacian: " << lap(0) <<  " amp  " 
      //  << temp.amp(w,4) << endl;
      //lap(e, w)*=-0.5;
      Kin(e, w)*=-0.5; 
    }
  }

  //cout << "laplacian " << lap(0) << endl;
  //cout << "Calculating kinetic energy done \n";
}


//----------------------------------------------------------------------

void System::generatePseudo(vector <vector <string> > & words,
                            Pseudopotential * & pseudo) {
  pseudo=new Pseudopotential;
  pseudo->read(words, this);
}



//----------------------------------------------------------------------

int write_xyz(System * sys, ostream & os) {
  
  vector <string> atomlabels;
  sys->getAtomicLabels(atomlabels);
  os<<atomlabels.size()<<endl;
  os << endl; //need a empty line here!!!
  Array1 <doublevar> pos(3);
  for(unsigned int i=0; i<atomlabels.size(); i++) {
    sys->getIonPos(i, pos);
    os<<atomlabels[i] <<" "<< pos(0)
        <<" "<<pos(1)
        <<" "<< pos(2)<<endl;
  }
  return 1;
}

//----------------------------------------------------------------------
