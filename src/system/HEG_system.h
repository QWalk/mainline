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

There are two choices for Coulomb e-e interaction:
\li Ewald \n
       <tt>interaction { Ewald }</tt>
\li truncated Coulomb a.k.a. MPC [Fraser et al., PRB 53, 1814 (1994)] \n
       <tt>interaction { truncCoul }</tt>
  
and some others for model calculations:
\li Gaussian e-e well (needs overall amplitude and standard deviation,
     i.e., interaction range, to be given in the \c system input deck) \n
       <tt>interaction { Gauss amp -5. stdev .5 }</tt>
     
At present, the default is the truncated Coulomb.

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
    labels.push_back("ORIGIN");
  }

  virtual int nIons() { return 1; }
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

  Array1 <int> nspin;             //!< number of spin-up and -down electrons
  int totnelectrons;              //!< total number of electrons
  
  Array2 <doublevar> centerpos;
  //!< position of centers to which basis is bound, could perhaps be simplified in this system

  Array1 <doublevar> origin;      //!< the origin of the simulation cell

  Array2 <doublevar> latVec;
  //!< lattice vectors, first index is a,b,c, second is x,y,z
  Array1 <doublevar> kpt;         //!< the k-point we're simulating at.    
  Array2 <doublevar> recipLatVec; //!<reciprocal lattice vectors
  Array2 <doublevar> normVec;     //!< normal vectors to the sides, pointing out
  Array2 <doublevar> corners;     //!< the position of the corner by moving one lattice vector
  doublevar smallestheight;       //!< smallest distance that spans the cell
  doublevar cellVolume;           //!< Simulation cell volume

  Array2 <doublevar> gpoint;      //!< A list of non-zero g points in the Ewald sum
  Array1 <doublevar> gweight;
  //!< A list of the weights(\f$4\pi exp(|g|^2/4 \alpha^2) \over V_{cell}|g|^2\f$)
  int ngpoints;                   //!< number of k points in ewald sum
  doublevar alpha;                //!< the Ewald parameter
  doublevar self_ee;              //!< self electron-electron energy
  doublevar xc_correction;        //!< exchange-correlation correction

  doublevar backgr_trunc;
  //!< constant to adjust truncated Coulomb for positive background charge
  
  doublevar Gauss_a;              //!< amplitude of Gaussian interaction
  doublevar Gauss_s;              //!< std. deviation (i.e., range) of Gaussian interaction
  doublevar Gauss_s2;             //!< \c Gauss_s squared

  int eeModel;
  /*!< \brief
   *  description for chosen e-e interaction model (used in showinfo())
   *  \li 0  default (truncated Coulomb)
   *  \li 1  Ewald
   *  \li 2  truncated Coulomb
   *  \li 3  short-range gaussian well
  */
  
  bool same_spin_int;             //!< switch for same spin interaction
  bool diff_spin_int;             //!< switch for different spin interaction

  doublevar (HEG_system::*calcLocChoice)(Sample_point *);
  //!< pointer to the selected local energy evaluator
  
  /*! \brief
    Initializations for Ewald sum evaluation
  */
  int setupEwald(Array2 <doublevar> & crossProduct);
    
  /*! \brief
    Set the constant ewald terms
   */
  void constEwald();

  /*! \brief
    electron-electron interaction (Ewald formula)
   */
  doublevar ewaldElectron(Sample_point * sample);

  
  /*! \brief
    Initializations for truncated Coulomb evaluation
  */
  int setupTruncCoulomb(); 


  /*! \brief
    evaluates local energy for Ewald model
  */
  doublevar calcLocEwald(Sample_point *);

  /*! \brief
    evaluates local energy for truncated Coulomb model
  */
  doublevar calcLocTrunc(Sample_point *);

  /*! \brief
    evaluates local energy for Gaussian short-range interaction
  */
  doublevar calcLocGauss(Sample_point *);


  //These variables are for putting a perturbing potential on the system
  doublevar calcLocPerturb(Sample_point *);
  int nperturb;
  Array2 <doublevar> perturb_pos;
  Array1 <doublevar> perturb_strength;
  Array1 <doublevar> perturb_alpha;
  Array1 <doublevar> perturb_spin;

  
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
