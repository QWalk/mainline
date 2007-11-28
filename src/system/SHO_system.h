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

#ifndef SHO_SYSTEM_H_INCLUDED
#define SHO_SYSTEM_H_INCLUDED
#include "System.h"
#include "Pseudopotential.h"
#include "Particle_set.h"

/*!
The Hamiltonian is:
\f$H= \omega*x^2 \f$
*/
class SHO_system:public System {
public:
  int showinfo(ostream & os);
  int generateSample(Sample_point * &);
  virtual void makeCopy(System *& ptr) {
    ptr=new SHO_system(*this);
  }

  SHO_system() { ee_interaction=0;
    origin.Resize(3); origin=0;
  };
  void notify(change_type, int);
  int read(vector <string> & words, unsigned int & pos);
  doublevar calcLoc(Sample_point *);

  virtual int nelectrons(int s) {
    assert(s==0 || s==1);
    return nspin(s);
  }


  virtual void getAtomicLabels(vector <string> & labels) {
    labels.clear();
  }

  virtual int ndim() { return omega.GetDim(0);}
  virtual int nIons() { return 0; }
  virtual void getIonPos(int i, Array1 <doublevar> & pos) {
    error("getIonPos not implemented");
  }
  virtual doublevar getIonCharge(const int number) 
  { error("getIonCharge not implemented"); return 0; } 
  virtual void setIonPos(int i, Array1 <doublevar> & pos) 
    { error("setIonPos not implemented"); }
  
  virtual void getEquivalentCenters(Array2 <int> & equiv_centers,
                                    Array1 <int> & ncenters_atom,
                                    Array2 <int> & displacements) {
    error("getEquivalentCenters() not implemented");
  }


private:
  
  friend class SHO_sample;
  Array1 <int> nspin;
  int ee_interaction;
  Array1 <doublevar> omega;
  Array1 <doublevar> origin;

};

#endif //SHO_SYSTEM_H_INCLUDED
//----------------------------------------------------------------------
