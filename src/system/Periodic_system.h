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

#ifndef PERIODIC_SYSTEM_H_INCLUDED
#define PERIODIC_SYSTEM_H_INCLUDED

#include "Qmc_std.h"
#include "System.h"
#include "Pseudopotential.h"
#include "Particle_set.h"
#include "Pbc_enforcer.h"
class Periodic_sample;

/*!
\brief 
Represents a periodic system. Keyword: PERIODIC

For the ewald interaction, we use the formula:

\f[
V=.5 \sum_{l}\sum_{i \ne j} q_i q_j \frac{erfc(\alpha |{\bf r_{ij}+l}|)}{r_{ij}}
  + \frac{4\pi}{V} 
\sum_{{\bf k} \ne 0}\frac{e^{-|{\bf k}|^2/4\alpha^2}}{|{\bf k}|^2} {1 \over 2}
((\sum_i q_i cos({\bf r_i \cdot k}))^2+(\sum_i q_i sin({\bf r_i \cdot k}))^2 )
-.5 \sum_{i \ne j} q_i q_j \frac{\pi}{V \alpha^2} 
+ .5 \sum_i q_i^2 (\sum_{l \ne 0} \frac{erfc(\alpha |{\bf l}|)}{|{\bf l}|}
-\frac{2\alpha}{\sqrt{2}}-\frac{\pi}{V\alpha^2})

\f]

\todo
Stop using Particle_set for the ion positions.  Use Pbc_enforcer instead of 
having the member function.
 */
class Periodic_system:public System
{
public:
  int showinfo(ostream & os);
  int generateSample(Sample_point * &);
  virtual void makeCopy(System *& ptr) {
    ptr=new Periodic_system(*this);
  }

  Periodic_system() {}

  void notify(change_type, int);
  int read(vector <string> & words, unsigned int & pos);

  doublevar calcLoc(Sample_point *);


  /*!
    Enforce the periodic boundary conditions on the position.
    Returns one if the vector was changed.
   */
  int enforcePbc(Array1 <doublevar> & pos, Array1 <int> & nshifted);


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
  
  virtual void getCenterLabels(vector <string> & cenlabels) {
    int ncen=equiv_atom.GetDim(0);
    cenlabels.resize(ncen);
    for(int cen=0; cen < ncen; cen++) {
      cenlabels[cen]=atomLabels[equiv_atom(cen)];
    }
  }


  virtual void getEquivalentCenters(Array2 <int> & equiv_centers_,
                                    Array1 <int> & ncenters_atom_, 
                                    Array2 <int> & displacements) {
    equiv_centers_.Resize(equiv_centers.GetDim(0), equiv_centers.GetDim(1));
    ncenters_atom_.Resize(ncenters_atom.GetDim(0));
    equiv_centers_=equiv_centers;
    ncenters_atom_=ncenters_atom;
    displacements=center_displacement;
  }

  int getBounds(Array2 <doublevar> & latvec) {
    latvec=latVec;
    return 1;
  }

  int getRecipLattice(Array2 <doublevar> & gvec) {
    gvec=recipLatVec;
    return 1;
  }
  
  int getPrimRecipLattice(Array2 <doublevar> & gvec) { 
    gvec=prim_recip_vec;
    return 1;
  }
    
  int getPrimLattice(Array2 <doublevar> & gvec) { 
    gvec=primlat;
    return 1;
  }

  void kpoint(Array1 <doublevar> & kp) {
    kp=kpt;
  }

  virtual void getorigin (Array1 <doublevar> & o) {
    o=origin;
  }

private:
  friend class Periodic_sample;

  Array1 <int> nspin;
  vector <string> atomLabels;


  Particle_set ions;

  Array2 <doublevar> centerpos;
  Array1 <int> equiv_atom; //!< For a given center, what atom it corresponds to
  Array2 <int> equiv_centers; //!< for a given atom, what centers correspond to it
  Array1 <int> ncenters_atom; //!< number of centers on a given atom
  Array2 <int> center_displacement; //!< in what direction the centers have been displaced
  

  Array1 <doublevar> origin;  //!< the origin of the simulation cell

  Array2 <doublevar> latVec;
  //!< lattice vectors, first index is a,b,c, second is x,y,z
  Array1 <doublevar> kpt; //!< the k-point we're simulating at.    
  Array2 <doublevar> recipLatVec; //!<reciprocal lattice vectors
  Array2 <doublevar> prim_recip_vec; //!< lattice vectors of the reciprocal primitive cell
  Array2 <doublevar> primlat; //!< lattice vectors of the primitive cell
  Array2 <doublevar> normVec;  //!< normal vectors to the sides, pointing out
  Array2 <doublevar> corners; //!< the position of the corner by moving one lattice vector

  Array2 <doublevar> gpoint; //!< A list of non-zero g points in the ewald sum
  Array1 <doublevar> gweight;
  //!< A list of the weights(\f$4\pi exp(|g|^2/4 \alpha^2) \over V_{cell}|g|^2\f$)

  Array1 <doublevar> ion_cos, ion_sin;
  //!< The sums for each g point over the ions.
  doublevar ion_ewald; //!< ionic ewald energy

  int ngpoints; //!< number of k points in ewald sum
  int totnelectrons; //!< number of electrons
  doublevar smallestheight; //!< smallest distance that spans the cell
  doublevar alpha; //!< the ewald parameter
  doublevar cellVolume; //!< Simulation cell volume

  doublevar self_ii; //!< self ion-ion energy
  doublevar self_ei; //!< self electron-ion energy
  doublevar self_ee; //!< self electron-electron energy
  doublevar xc_correction; //!<exchange-correlation correction


  /*!
    Set the constant ewald terms
   */
  void constEwald();

  /*!
    Calculate the ion-ion interaction from the ions stored here.
   */
  doublevar ewaldIon();

  /*!
    electron-ion interaction and electron-electron interaction
   */
  doublevar ewaldElectron(Sample_point * sample);
  
  Array1 <doublevar> ion_polarization;
};



/*!
  \brief
  An approximation to erfc(complimentary error function).  

  Supposedly has an absolute error of less than 1.5e-07.
 */
inline doublevar erfcm(doublevar x) {
  const doublevar p=.3275911, a=.254829592, b=-.284496736, c=1.421413741;
  const doublevar d=-1.453152027, e=1.061405429;

  doublevar y=1./(1.+p*x);
  return y*(a+y*(b+y*(c+y*(d+y*e))))*exp(-x*x);
}

#endif //PERIODIC_SYSTEM_H_INCLUDED
//----------------------------------------------------------------------
