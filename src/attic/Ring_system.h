/*
 
Copyright (C) 2007 Lucas K. Wagner
 with modifications by Pavel Vagner

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


#ifndef RING_SYSTEM_H_INCLUDED
#define RING_SYSTEM_H_INCLUDED
#include "System.h"
#include "Pseudopotential.h"
#include "Particle_set.h"
#include "Pbc_enforcer.h"
/*!
The Hamiltonian is:
\f$H= \Pi^2/2m^* + V_{ee}(|r_i-r_j|)\f$

where \f$ \Pi=P^2+2A*P+A^2 \f$, and 
\f$ A=\frac{\Phi}{\Phi_0} \frac{\hbar}{R_0} \f$

\f$ V_{ee}(|r_i-r_j|) = \frac{1}{|r_i-r_j|+r_0}\f$

Input variables are:

\f$m^*\f$: effective mass of the electron

\f$r_0\f$: cutoff of e-e interaction

\f$R_0\f$: radius of ring

\f$\frac{\Phi}{\Phi_0} \f$: quantum of magnetic flux through the ring.

*/
class Ring_system:public System {
public:
  int showinfo(ostream & os);
  int generateSample(Sample_point * &);
  virtual void makeCopy(System *& ptr) {
    ptr=new Ring_system(*this);
  }

  Ring_system() {
    exp_interaction=0;
    fixed_phase=0;
  };
  void notify(change_type, int);
  int read(vector <string> & words, unsigned int & pos);
  doublevar calcLoc(Sample_point *);

  virtual void calcKinetic(Wavefunction_data *, 
                           Sample_point *,
                           Wavefunction *,
                           Array1 <doublevar> &);
  void calcLocSeparated(Sample_point *, Array1<doublevar> &) {
    cout << "Not implemented!" << endl; 
  }
  
  virtual int nelectrons(int s) {
    assert(s==0 || s==1);
    return nspin(s);
  }


  virtual void getAtomicLabels(vector <string> & labels) {
    labels.clear();
  }

  virtual int ndim() { return 1;}
  virtual int nIons() { error("nIons() not implemented"); return 1; }
  virtual void getIonPos(int i, Array1 <doublevar> & pos) {
    error("getIonPos not implemented");
  }
  virtual void setIonPos(int i, Array1 <doublevar> & pos) 
    { error("setIonPos not implemented"); }
  virtual doublevar getIonCharge(const int number) 
  { error("getIonCharge not implemented"); return 0; } 

    
    
  virtual void getEquivalentCenters(Array2 <int> & equiv_centers,
                                    Array1 <int> & ncenters_atom,
                                    Array2 <int> & displacements) {
    error("getEquivalentCenters() not implemented");
  }


private:
  
  friend class Ring_sample;
  Array1 <int> nspin;
  doublevar length;
  doublevar radius;
  doublevar effmass;
  doublevar flux;
  /*
    e-e interactions:
    1. dist=2*radius * sin(pi*fabs(x1-x2)/length)
       V_ee= 1/(dist+ee_0)

    2. dist=fabs( x1-x2 )
       V_ee=ee_0*exp(dist/ee_range)
  */
  doublevar ee_0;
  doublevar ee_range;

  int ee_interaction;
  int exp_interaction;
  int fixed_phase;

  int nbarrier;
  Array1 <doublevar> barrier_strength;
  Array1 <doublevar> barrier_lower;
  Array1 <doublevar> barrier_upper;

};

#endif //RING_SYSTEM_H_INCLUDED
//----------------------------------------------------------------------
