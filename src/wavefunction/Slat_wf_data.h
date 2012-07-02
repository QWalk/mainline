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
//#include "Slat_wf.h"
//#include "Cslat_wf.h"
template <class T> class Slat_wf;
//class Cslat_wf;
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


  virtual int valSize() {
    if(optimize_mo) return 0;
    if(optimize_det) return 2*detwt.GetDim(0);
    else return nfunc*2;
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
    if(optimize_mo) 
      return orbitals_for_optimize_mo.GetSize()*(molecorb->nMoCoeff()/molecorb->getNmo());
    else if(optimize_det){ return ncsf-1; } 
    /*
      if(use_csf)
        if(all_weights)
          return ncsf;
        else
          return ncsf-1;
      else
        if(all_weights)
          return ndet;
        else
          return ndet-1;
    }
    */
    else return 0;
  }

private:
  void init_mo();
  friend class Slat_wf<doublevar>;
  friend class Slat_wf<dcomplex>;
  friend class Cslat_wf;
  friend class FSlat_wf;

  void calc_determinant_evaluation_order(const int & f, const int & s, 
					 const Array3<Array1<int> > & changes,
					 const Array3 <int> & nchanges,
					 const Array3< Array1 <int> > & occ,
					 Array1 <int> & order,
					 int dstart, int dend, int dlevel );
  //!< Compute best order for determinant evaluation for determinants numbered inclusively between dstart,dend and for excitation level dlevel


  int cost_iterative_determinant_order(const int & func, const int & spin, 
				       const Array3<Array1<int> > & changes,
				       const Array3 <int> & nchanges,
				       const Array3< Array1 <int> > & occ,
				       const Array1 <int> & order);
  //!< Helper function: sum cost of updating determinants using iterative algorithm in O(N) units


  Array3 < Array1 <int> > occupation;
  //!< The place in totoccupation for the molecular orbitals for (fn, det, spin) Used to look up the right MO from MO_matrix.
  Array3 < Array1 <int> > occupation_orig;
  //!< The molecular orbitals numbers for (function, determinant, spin)
  Array1 <doublevar> detwt;
  Array1 < Array1 <int> > totoccupation; //!< all the molecular orbitals for a given spin

  Array3 <int> occupation_nchanges;           //!< number of changes to occupations for (fn,det,spin)
  Array3 <int> evaluation_order;              //!< Order to evaluate determinants (fn,idx,spin)
  Array3 <int> occupation_first_diff_change_from_last_det;           //!< change index of first change different from previous determinant (fn,det,spin)
  Array3 < Array1 <int> > occupation_changes; //!< list of changes to occupations for (fn,det,spin)
  int max_occupation_changes;

  Array1 <int> spin;  //!< spin as a function of electron
  Array1 <int> opspin;//!< the opposing spin for an electron
  Array1 <int> rede;  //<! The number of the electron within its spin
  Array1 <int> nelectrons; //!< number of electrons as a function of spin

  int nmo;        //!<Number of molecular orbitals
  int ndet;       //!<Number of determinants
  int ndim;       //!<Number of (spacial) dimensions each electron has
  int nfunc;      //!<Number of separate functions this represents

  int optimize_mo; //!< whether to optimize MO basis
  Array1 <int> orbitals_for_optimize_mo; //!< which orbitals should be optimized 
  string mo_place; //!< where to place the new mo's
  int optimize_det; //!< whether to optimize determinant coefficients
  Array1 <Array1 <doublevar> > CSF;
//  int use_csf; //whether to use csfs instead of pure dets
  int ncsf; 
  int sort; //whether to sort det. weights by size
//  int all_weights; 
  
  
  General_MO_matrix * genmolecorb; 
  MO_matrix * molecorb;
  int use_complexmo;
  int use_iterative_updates; //!<Whether to use "fast"/low-memory Nukala-Kent iterative updates for multideterminants
  Complex_MO_matrix * cmolecorb;

};

#endif //SLAT_WF_DATA_H_INCLUDED
//------------------------------------------------------------------------
