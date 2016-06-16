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


#ifndef SYSTEM_H_INCLUDED
#define SYSTEM_H_INCLUDED

#include "Qmc_std.h"
#include "Wavefunction.h"
#include "Force_fitter.h"
#include "Pseudopotential.h"
class Sample_point;
class Wavefunction_data;
class Wavefunction;
//class Pseudopotential;
/*!
  \brief
  Holds the information about the system; knows how to calculate energies
  and generate Sample_points.


 */
class System
{
public:
  virtual ~System()
  {}
  virtual int showinfo(ostream & os)=0;
  virtual int read(vector <string> & words, unsigned int & pos)=0;
  virtual int generateSample(Sample_point * &)=0;
  //virtual void notify(change_type, int)=0;
  virtual void generatePseudo(vector < vector < string > > &,
                              Pseudopotential * & pseudo);
  virtual doublevar calcLoc(Sample_point *)=0;
  /*!
    Added by Huihuo Zheng for extended Koopmans' theorem; 
   */
  virtual void calcLocSeparated(Sample_point *, Array1 <doublevar> &)=0;
  virtual void calcLocWithTestPos(Sample_point * sample, Array1 <doublevar> & tpos, Array1 <doublevar> & Vtest) {};
  /*!
    Return the local potential separated into onebody and twobody 
    terms.  Onebody is per quantum particle, and twobody is a nxn upper
    triangular matrix (i.e., filled in (i,j) for i < j)
    */
  virtual void separatedLocal(Sample_point *,Array1 <doublevar> & onebody,
                     Array2 <doublevar> & twobody) { 
    error("separatedLocal not available for this system");
  }
  virtual void locDerivative(int ion, Sample_point *, Force_fitter &,
                             Array1 <doublevar> & der) {
    error("this system doesn't support locDerivative");
  }


  virtual void calcKinetic(Wavefunction_data *,
                           Sample_point *,
                           Wavefunction *,
                           Array1 <doublevar> &);
  /*!
    \brief 
    This is to return the kinetic energy of each electron
     ------------ lap(e, iw)
    e (electron index), iw (wavefunction index)
  */
  virtual void calcKineticSeparated(Wavefunction_data *,
				    Sample_point *,
				    Wavefunction *, Array2<doublevar> & Kin);

  /*!
    \brief
     give the reciprocal lattice.  return 1 on success and 0 on error.
  */
  virtual int getRecipLattice(Array2 <doublevar> & gvec) {
    //gvec.Resize(3,3);
    //gvec=0.0;
    //for(int d=0; d< 3; d++) gvec(d,d)=1.0;
    return 0;
  }
  /*!
    \brief 
   return the primitive lattice (usually given in input)
    */
  virtual int getPrimLattice(Array2 <doublevar> & gvec) { 
    return 0;
  }

  /*!
    \brief 
   return the primitive reciprocal lattice (usually given in input)
    */
  virtual int getPrimRecipLattice(Array2 <doublevar> & gvec) { 
    return 0;
  }
  
  /*!
    \brief 
    Get a representative box for the simulation cell
   */
  virtual int getBounds(Array2 <doublevar> & latvec,Array1 <doublevar> & origin) { return 0;};

  
  virtual void kpoint(Array1 <doublevar> & kp) {
    kp.Resize(3); kp=0;
  }
  virtual void getorigin (Array1 <doublevar> & o) {
    //o.Resize(3);o=0;
    error("Origin not implemented for this method");
  }
  /*!
    \brief
    Number of electrons for spin s
   */
  virtual int nelectrons(int s)=0;

  virtual int nIons()=0;
  virtual void getIonPos(int i, Array1 <doublevar> & pos)=0;
  /*!
    \brief
    Get an ion's charge(in a.u.)
  */
  virtual doublevar getIonCharge(const int number)=0;

  virtual int ndim() { return 3;}

  /*!
    \brief
    Gives a listing of the labels of the atomic coordinates,
    in order.
   */
  virtual void getAtomicLabels(vector <string> & labels)=0;

  virtual void getCenterLabels(vector <string> & labels) {
    getAtomicLabels(labels);
  }
  virtual void getEquivalentCenters(Array2 <int> & equiv_centers,
                                    Array1 <int> & ncenters_atom, 
                                    Array2 <int> & displacements)=0;

  virtual void makeCopy(System *& ptr)=0;

};

int allocate(vector <string> &, System * &);


int write_xyz(System * sys, ostream & os);

#endif //SYSTEM_H_INCLUDED

//----------------------------------------------------------------------
