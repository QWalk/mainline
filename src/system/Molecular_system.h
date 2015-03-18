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

#ifndef MOLECULAR_SYSTEM_H_INCLUDED
#define MOLECULAR_SYSTEM_H_INCLUDED
#include "System.h"
#include "Pseudopotential.h"
#include "Particle_set.h"
#include "Pbc_enforcer.h"

class Molecular_system:public System
{
public:


  int showinfo(ostream & os);
  int generateSample(Sample_point * &);
  virtual void makeCopy(System *& ptr) {
    ptr=new Molecular_system(*this);
  }
  Molecular_system() {};

  Molecular_system(Molecular_system & sys) {
    ions=sys.ions;
    nspin.Resize(sys.nspin.GetDim(0));
    nspin=sys.nspin;
    atomLabels=sys.atomLabels;
    bounding_box=sys.bounding_box;
    use_bounding_box=sys.use_bounding_box;
  }

  void notify(change_type, int);
  int read(vector <string> & words, unsigned int & pos);



  doublevar calcLoc(Sample_point *);
  void calcLocWithTestPos(Sample_point * sample, Array1 <doublevar> & tpos, Array1 <doublevar> & Vtest); 
  void calcLocSeparated(Sample_point *, Array1<doublevar> & ); 
  //  doublevar calcPotNoNotify(const int e, Sample_point * sample); //added by Huihuo Zheng for extended Koopmans' theorem
  void locDerivative(int ion, Sample_point * sample, 
                     Force_fitter &,Array1 <doublevar> & der);
  virtual void separatedLocal(Sample_point *,Array1 <doublevar> & onebody,
                     Array2 <doublevar> & twobody); 
  
  virtual int nelectrons(int s) {
    assert(s==0 || s==1);
    return nspin(s);
  }
  virtual void getAtomicLabels(vector <string> & labels) {
    labels=atomLabels;
  }
  virtual int nIons() { return atomLabels.size(); }
  virtual void getIonPos(int i, Array1 <doublevar> & pos) {
    assert(pos.GetDim(0)==3);
    for(int d=0; d< 3; d++) {
      pos(d)=ions.r(d,i);
    }
  }
  virtual void setIonPos(int i, Array1 <doublevar> & pos);
  doublevar getIonCharge(const int ion)
  {
    return ions.charge(ion);
  }
  virtual void getEquivalentCenters(Array2 <int> & equiv_centers,
                                    Array1 <int> & ncenters_atom,
                                    Array2 <int> & displacements) {
    int nat=ions.size();
    equiv_centers.Resize(nat, 1);
    ncenters_atom.Resize(nat);
    displacements.Resize(nat,3);
    displacements=0;
    for(int at=0; at< nat; at++) {
      equiv_centers(at,0)=at;
      ncenters_atom(at)=1;
    }
  }

  virtual int getBounds(Array2<doublevar> & latvec,Array1 <doublevar> & origin);


private:
  friend class Molecular_sample;
  Particle_set ions;
  Array1 <int> nspin;

  vector <string> atomLabels;
  Pbc_enforcer bounding_box;
  int use_bounding_box;

  Array1 <doublevar> electric_field; 
  doublevar inirange;
};

#endif //MOLECULAR_SYSTEM_H_INCLUDED
//----------------------------------------------------------------------
