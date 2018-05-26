/*
 
Original Copyright (C) 2007 Lucas K. Wagner

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

#ifndef SLAT_WF_DATA_H_INCLUDED
#define SLAT_WF_DATA_H_INCLUDED

#include "Qmc_std.h"
#include "Wavefunction_data.h"
#include "MO_matrix.h"
#include "clark_updates.h"
#include "Orbital_rotation.h"
template <class T> class Slat_wf;
/*!
\brief
A Slater determinant or several determinants.  Also can use the same
molecular orbitals to create several separate combinations of determinants,
and return them as a vector.

<h3> Options </h3> <br>
<b> Required </b> <br>
<ul>

<li>
<b> DETWT </b> List of the weights of the determinants. 
</li>

<li>
<b> STATES </b> List of the occupations of the molecular orbitals, first
spin up, then spin down.  For example, a RHF determinant might be <br>
1 2 3 <br>
1 2 3 <br>
and a UHF might be <br>
1 2 3 <br>
4 5 6 <br>
If there are more than one determinants, simply continue listing up and 
down occupations.<br>  Also, list several STATES functions, and they
will be added to the vector in order. 
</li>

<li>
<b> ORBITALS </b> This section is an input deck for an MO_matrix
</li>
</ul>

<b> Optional </b>
<ul>
<li>
<b> NSPIN </b>  Override the global parameter.
</li>
</ul>

 */
class Slat_wf_data : public Wavefunction_data
{
public:

  Slat_wf_data():molecorb(NULL) {}

  ~Slat_wf_data()
  {
    if(molecorb != NULL ) delete molecorb;
  }

  int optimize_mo; //!< whether to optimize MO basis
  virtual int valSize() {
    //if(optimize_mo) return 0;
    if(optimize_det) return 2*detwt.GetDim(0);
    else return nfunc*2;
  }
  virtual void lockInParms(Array1<doublevar> & parms){
    if(optimize_mo && optimize_det){
      Array1<doublevar> orbital_parms;
      orbital_parms.Resize(orbrot->nparms());
      for(int i=0;i<orbrot->nparms();i++){
        orbital_parms[i]=parms[ncsf-1+i];
      }
      orbrot->lockIn(orbital_parms);
      for(int i=0;i<orbrot->nparms();i++){
        parms[ncsf-1+i]=orbital_parms[i];
      }
    }else if(optimize_mo){
      orbrot->lockIn(parms);
    }
    return; 
  }
  virtual void getVarParms(Array1 <doublevar> & parms);
  virtual void setVarParms(Array1 <doublevar> & parms);
  virtual void linearParms(Array1 <bool> & is_linear);

  virtual void renormalize();

  virtual void read(vector <string> & words,
                    unsigned int & pos,
                    System * sys
                   );
  virtual int supports(wf_support_type );
  void generateWavefunction(Wavefunction *&);


  int showinfo(ostream & os);

  int writeinput(string &, ostream &);

  int nparms()
  {
    if(optimize_mo && optimize_det) return ncsf - 1 + orbrot->nparms();
    else if(optimize_mo) return orbrot->nparms();
    else if(optimize_det){ return ncsf-1; } 
    else return 0;
  }

  Slat_wf_data * clone() const{
    return new Slat_wf_data(*this);
  }

private:
  void init_mo();
  friend class Slat_wf<doublevar>;
  friend class Slat_wf<dcomplex>;
  friend class Cslat_wf;

  Array3 < Array1 <int> > occupation;
  //!< The place in totoccupation for the molecular orbitals for (fn, det, spin) Used to look up the right MO from MO_matrix.
  Array3 < Array1 <int> > occupation_orig;
  //!< The molecular orbitals numbers for (function, determinant, spin)
  Array1 <log_value<doublevar> > detwt;
  Array1 < Array1 <int> > totoccupation; //!< all the molecular orbitals for a given spin

  int max_occupation_changes;

  Array1 <int> nelectrons; //!< number of electrons as a function of spin

  int nmo;        //!<Number of molecular orbitals
  int ndet;       //!<Number of determinants
  int ndim;       //!<Number of (spacial) dimensions each electron has
  int nfunc;      //!<Number of separate functions this represents

  Array1 <int> orbitals_for_optimize_mo; //!< which orbitals should be optimized 
  string mo_place; //!< where to place the new mo's
  int optimize_det; //!< whether to optimize determinant coefficients
  Array1 <Array1 <doublevar> > CSF;
  int ncsf; 
  int sort; //whether to sort det. weights by size
  
  
  General_MO_matrix * genmolecorb; 
  MO_matrix * molecorb;
  int use_complexmo;
  bool use_clark_updates; //!<Use Bryan Clark's updates.
  Excitation_list excitations;
  Complex_MO_matrix * cmolecorb;
  Orbital_rotation * orbrot;

};

#endif //SLAT_WF_DATA_H_INCLUDED
//------------------------------------------------------------------------
