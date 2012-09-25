/*
 
Copyright (C) 2007 Zachary Helms
 with further modifications by Lucas K. Wagner

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

#include "Program_options.h"
#include "Determinant_select_method.h"
#include "qmc_io.h"
#include "System.h"
#include "MatrixAlgebra.h"
#include "ulec.h"
/*!
Read the "words" from the method section in the input file
via doinput() parsing, and store section information in private
variables orbs, resolution, and minmax.
Set up MO_matrix and Sample_point objects for wavefunction from input
*/
void Determinant_select_method::read(vector <string> words,
                       unsigned int & pos,
                       Program_options & options)
{
  cout << "rrreeead " << endl;
  
  sys=NULL;
  allocate(options.systemtext[0],  sys);

  vector <vector < string> > orbgroup_txt;
  pos=0;
  vector < string> dummy;
  while( readsection(words,pos,dummy,"ORB_GROUP")) { 
    orbgroup_txt.push_back(dummy);
  }
  cout << "aaaahree " << endl;
  orbital_groups.Resize(orbgroup_txt.size());
  for(int i=0; i< orbital_groups.GetDim(0); i++) { 
    int norb_this=orbgroup_txt[i].size();
    orbital_groups[i].Resize(norb_this);
    for(int j=0; j< norb_this; j++) { 
      orbital_groups[i][j]=atoi(orbgroup_txt[i][j].c_str())-1;
    }
  }

  cout << "here " << endl;


  vector <string> orbtext;
  if(!readsection(words, pos=0, orbtext, "ORBITALS")){
      error("Need ORBITALS");
  }
  allocate(orbtext, sys, mymomat);

  vector <string> occupied_up, virtu_up, occupied_down, virtu_down;
  if(!readsection(words, pos=0, occupied_up, "OCCUPIED_UP"))
    error("Need OCCUPIED section");
  if(!readsection(words, pos=0, virtu_up, "VIRTUAL_UP"))
    error("Need VIRTUAL section");
  if(!readsection(words, pos=0, occupied_down, "OCCUPIED_DOWN"))
    error("Need OCCUPIED section");
  if(!readsection(words, pos=0, virtu_down, "VIRTUAL_DOWN"))
    error("Need VIRTUAL section");
  nocc.Resize(2); nvirt.Resize(2);
  nocc(0)=occupied_up.size();
  nocc(1)=occupied_down.size();
  nvirt(0)=virtu_up.size();
  nvirt(1)=virtu_down.size();
  occ.Resize(2,max(nocc(0),nocc(1)));
  virt.Resize(2,max(nvirt(0),nvirt(1)));
  for(int i=0; i< nocc(0); i++) occ(0,i)=atoi(occupied_up[i].c_str())-1;
  for(int i=0; i< nocc(1); i++) occ(1,i)=atoi(occupied_down[i].c_str())-1;

  for(int i=0; i< nvirt(0); i++) virt(0,i)=atoi(virtu_up[i].c_str())-1;
  for(int i=0; i< nvirt(1); i++) virt(1,i)=atoi(virtu_down[i].c_str())-1;
  cout << "done read " << endl; 
  mywalker=NULL;
  sys->generateSample(mywalker);

  mymomat->buildLists(orbital_groups);
  cout << " done build " << endl;
}

//----------------------------------------------------------------------
int Determinant_select_method::showinfo(ostream & os) { 
  os << "Determinant_select method " << endl;
  return 1;
}

#include "Generate_sample.h"

//----------------------------------------------------------------------
void Determinant_select_method::run(Program_options & options, ostream & output) {
  int nsamples=1000;
  for(int g=0; g< orbital_groups.GetDim(0); g++) { 
    int norb=orbital_groups[g].GetDim(0);
    Array1 <Array1 <doublevar> > r_sample(nsamples);
    generate_mo_sample(mywalker,sys,mymomat,g,nsamples,r_sample);
    Array2 <doublevar> vals(nsamples,norb);
    Array2 <doublevar> tmp_val(norb,2);
    Array1 <doublevar> orb_norm(norb,0.0);
    for(int s=0;s < nsamples;s++) {
      mywalker->setElectronPos(0,r_sample(s));
      mymomat->updateVal(mywalker,0,g,tmp_val);
      for(int i=0; i< norb; i++) {
        vals(s,i)=tmp_val(i,0);
        cout << s << i << " vals " << vals(s,i) << endl;
        orb_norm(i)+=vals(s,i)*vals(s,i);
      }
    }
    for(int i=0; i< norb; i++) 
      orb_norm(i)/=nsamples;
    
    Array4 <doublevar> voverlap(norb,norb,norb,norb);
    voverlap=0.0;
    for(int e1=0; e1 < nsamples; e1++) { 
      cout << "sample " << e1 << endl;
      for(int e2=e1+1; e2< nsamples; e2++) { 
        mywalker->setElectronPos(0,r_sample(e1));
        mywalker->setElectronPos(1,r_sample(e2));
        mywalker->updateEEDist();
        Array1 <doublevar> ree(5);
        mywalker->getEEDist(0,1,ree);
        doublevar rij=1.0/ree(0);
        //cout << "rij " << rij << endl;

        for(int i=0; i< norb; i++) { 
          for(int j=0; j< norb; j++) { 
            for(int k=0; k< norb; k++) { 
              for(int l=0; l< norb; l++) { 
                voverlap(i,j,k,l)+=vals(e1,i)*vals(e2,j)*rij*vals(e1,k)*vals(e2,l);
              }
            }
          }
        }
      }
    }
    Array2 <int> lookup_occ=occ, lookup_virt=virt;
    for(int s=0; s< 2; s++) { 
      for(int i=0; i< occ.GetDim(1); i++) { 
        for(int j=0; j< norb; j++) { 
          if(occ(s,i)==orbital_groups(g)(j))
            lookup_occ(s,i)=j;
        }
      }
      for(int i=0; i< virt.GetDim(1); i++) { 
        for(int j=0; j< norb; j++) { 
          if(virt(s,i)==orbital_groups(g)(j))
            lookup_virt(s,i)=j;
        }
      }
    }
    doublevar renorm=1.0/(nsamples*(nsamples-1)/2.0);
    for(int s1=0; s1 < 2; s1++) { 
    for(int s2=s1; s2 < 2; s2++) { 
      output << "----------Between channels " << s1 << " and " << s2 << endl;
      for(int i=0; i< nocc(s1); i++) { 
      for(int j=0; j< nocc(s2); j++) { 
        for(int k=0; k< nvirt(s1); k++) { 
          for(int l=0; l< nvirt(s2); l++) { 
            int oi=lookup_occ(s1,i);
            int oj=lookup_occ(s2,j);
            int ok=lookup_virt(s1,k);
            int ol=lookup_virt(s2,l);
            //cout << "i " << i << " j " << j << " k " << k << " l " << l <<endl;
            //cout << "oi " << oi << " oj " << oj << " ok " << ok << " ol " << ol << endl;
            output << occ(s1,i)+1 << " " << occ(s2,j)+1 << " " << virt(s1,k)+1 << " "
              << virt(s2,l)+1 << " " 
              << voverlap(oi,oj,ok,ol)*renorm/sqrt(orb_norm(oi)*orb_norm(oj)*orb_norm(ok)*orb_norm(ol)) << endl;
          }
        }
      }
    }
    }
    }
    
            
  }
}
