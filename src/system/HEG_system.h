/*
 
Copyright (C) 2007 Lucas K. Wagner, Jindrich Kolorenc

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

#ifndef HEG_SYSTEM_H_INCLUDED
#define HEG_SYSTEM_H_INCLUDED

#include "Qmc_std.h"
#include "System.h"
#include "Particle_set.h"
#include "Pbc_enforcer.h"
class HEG_sample;

/*!
\brief 
Represents a homogeneus periodic system (homogeneous electron gas), created
as a simplification of crystalline Periodic_system.

Keyword: HEG

There are two choices for e-e interaction:
  1) Ewald
  2) truncated Coulomb a.k.a. MPC [Fraser et al., PRB 53, 1814 (1994)]
At present, the default is 2.

*/

class HEG_system:public System
{
public:
  int showinfo(ostream & os);
  int generateSample(Sample_point * &);
  virtual void makeCopy(System *& ptr) {
    ptr=new HEG_system(*this);
  }

  HEG_system() {}

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
    labels.clear();
  }

  virtual int nIons() { return 0; }
  virtual void getIonPos(int i, Array1 <doublevar> & pos) {
    error("getIonPos not implemented");
  }
  virtual void setIonPos(int i, Array1 <doublevar> & pos) {
    error("setIonPos not implemented");
  }
  doublevar getIonCharge(const int ion) {
    error("getIonCharge not implemented");
    return 0.0;
  }
  
  virtual void getEquivalentCenters(Array2 <int> & equiv_centers_,
                                    Array1 <int> & ncenters_atom_, 
                                    Array2 <int> & displacements) {
    error("getEquivalentCenters() not implemented");
  }

  int getBounds(Array2 <doublevar> & latvec) {
    latvec=latVec;
    return 1;
  }

  int getRecipLattice(Array2 <doublevar> & gvec) {
    gvec=recipLatVec;
    return 1;
  }
    
  void kpoint(Array1 <doublevar> & kp) {
    kp=kpt;
  }


private:

  friend class HEG_sample;

  Array1 <int> nspin;

  Array2 <doublevar> centerpos;

  Array1 <doublevar> origin;  //!< the origin of the simulation cell

  Array2 <doublevar> latVec;
  //!< lattice vectors, first index is a,b,c, second is x,y,z
  Array1 <doublevar> kpt; //!< the k-point we're simulating at.    
  Array2 <doublevar> recipLatVec; //!<reciprocal lattice vectors
  Array2 <doublevar> normVec;  //!< normal vectors to the sides, pointing out
  Array2 <doublevar> corners; //!< the position of the corner by moving one lattice vector
  doublevar smallestheight;   //!< smallest distance that spans the cell

  int totnelectrons; //!< number of electrons
  doublevar cellVolume; //!< Simulation cell volume

  Array2 <doublevar> gpoint; //!< A list of non-zero g points in the ewald sum
  Array1 <doublevar> gweight;
  //!< A list of the weights(\f$4\pi exp(|g|^2/4 \alpha^2) \over V_{cell}|g|^2\f$)

  int ngpoints; //!< number of k points in ewald sum
  doublevar alpha; //!< the ewald parameter
  doublevar self_ee; //!< self electron-electron energy
  doublevar xc_correction; //!< exchange-correlation correction

  doublevar backgr_trunc; //!< constant to adjust truncated Coulomb for positive background charge


  /*!
    Initializations for Ewald sum evaluation
  */
  int setupEwald(Array2 <doublevar> & crossProduct);
    
  /*!
    Set the constant ewald terms
   */
  void constEwald();

  /*!
    electron-electron interaction (Ewald formula)
   */
  doublevar ewaldElectron(Sample_point * sample);

  
  /*!
    Initializations for truncated Coulomb evaluation
  */
  int setupTruncCoulomb(); 


  /*!
    local energy evaluators and pointer to the selected one
  */
  doublevar calcLocEwald(Sample_point *);
  doublevar calcLocTrunc(Sample_point *);
  doublevar (HEG_system::*calcLocChoice)(Sample_point *);
  /*!
    description for chosen model (for status info)
    0  default (truncated Coulomb)
    1  Ewald
    2  truncated Coulomb
  */
  int eeModel;  
  
};



/*!
  \brief
  An approximation to erfc(complimentary error function).  

  Supposedly has an absolute error of less than 1.5e-07.
 */
#ifndef PERIODIC_SYSTEM_H_INCLUDED
inline doublevar erfcm(doublevar x) {
  const doublevar p=.3275911, a=.254829592, b=-.284496736, c=1.421413741;
  const doublevar d=-1.453152027, e=1.061405429;

  doublevar y=1./(1.+p*x);
  return y*(a+y*(b+y*(c+y*(d+y*e))))*exp(-x*x);
}
#endif //PERIODIC_SYSTEM_H_INCLUDED

#endif //HEG_SYSTEM_H_INCLUDED
//----------------------------------------------------------------------
