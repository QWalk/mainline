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

#ifndef JASTROW2_WF_H_INCLUDED
#define JASTROW2_WF_H_INCLUDED
#include "Qmc_std.h"
#include "Array45.h"

#include "Jastrow2_one.h"
#include "Jastrow2_two.h"
#include "Jastrow2_three.h"
#include "Jastrow2_three_diffspin.h"

#include "Wavefunction.h"
#include "Wavefunction_data.h"
#include "Basis_function.h"

//######################################################################


/*!

*/
class Jastrow_group {
public:
  Jastrow_group() {
    has_one_body=0;
    has_two_body=0;
    two_body=NULL;

  }

  ~Jastrow_group() {
    if(two_body) delete two_body;
    two_body=NULL;

    for(int i=0; i< eibasis.GetDim(0); i++)
      if(eibasis(i)) delete eibasis(i);
    eibasis=NULL;

    for(int i=0; i< eebasis.GetDim(0); i++)
      if(eebasis(i)) delete eebasis(i);
    eebasis=NULL;
  }

  void set_up(vector <string> & words, System * sys);
  void updateEIBasis(int e, Sample_point * sample, Array3 <doublevar> & );
  void updateEEBasis(int e, Sample_point * sample, Array3 <doublevar> & );

  int maxEIBasis() { return maxbasis_on_center; }
  int nEEBasis() { return n_eebasis; }

  int hasOneBody() { return has_one_body; }
  int hasTwoBody() { return has_two_body; }
  int hasThreeBody() { return has_three_body; } 
  int hasThreeBodySpin() { return has_three_body_diffspin; } 
  int optimizeBasis() { return optimize_basis; }
  Jastrow_onebody_piece one_body;
  Jastrow_twobody_piece * two_body;
  Jastrow_threebody_piece three_body;
  Jastrow_threebody_piece_diffspin three_body_diffspin;

  int writeinput(string &, ostream &);
  int showinfo(string &, ostream &);
  int valSize();
  int nparms();

  void getVarParms(Array1 <doublevar> & parms);
  void setVarParms(Array1 <doublevar> & parms);

  void getEEbasisPlot(Array1 <doublevar> &, Array1 <doublevar> &);
  void getEIbasisPlotInfo(vector <string> &, Array1 <int> &);
  void getEIbasisPlot(int, Array1 <doublevar> &, Array1 <doublevar> &);

private:
  int check_consistency();
  vector <string> atomnames;
  
  int has_one_body;
  int has_two_body;
  int has_three_body;
  int has_three_body_diffspin;
  int optimize_basis;
  int have_diffspin;
  //Stuff to handle the electron-ion basis

  Array1 <Basis_function * > eibasis;
  Array1 <int> nfunc_eib; //!< number of functions/basis object
  int toteibasis; //!< total number of basis functions
  int maxeibasis; //!< maximum number of basis functions/object
  //Bug--need something that has max number/atom for sizing the basis output
  Array2 <int> atom2basis; //!< which basis object is on each atom
  Array1 <int> nbasis_at; //!< number of basis objects/atom
  int maxbasis_on_center; //!< maximum number of basis functions on a center


  //Electron-electron basis

  Array1 <Basis_function * > eebasis;
  int n_eebasis; //!< Total number of electron-electron
  Array1 <int> nfunc_eeb; //!< number of functions/object
  int maxeebasis;
  int nelectrons;
  int n_spin_up;

};

//######################################################################

class Jastrow2_wf;

/*!
  \brief data object for a Jastrow that supports Laplacian updates
 */
class Jastrow2_wf_data:public Wavefunction_data {
public:
  virtual void read(vector <string> & words,
                    unsigned int & pos,
                    System * sys
                   );
  virtual int supports(wf_support_type );
  virtual void getVarParms(Array1 <doublevar> &);
  virtual void setVarParms(Array1 <doublevar> &);
  virtual int nparms();
  virtual int valSize();
  virtual int writeinput(string &, ostream &);
  virtual int showinfo(ostream &);
  virtual void generateWavefunction(Wavefunction * &);
private:
  friend class Jastrow2_wf;
  Array1 < Jastrow_group > group;
  int nelectrons;
  int natoms;
};


//######################################################################

class Jastrow2_storage : public Wavefunction_storage {
    Jastrow2_storage(int nelectrons) {
      one_body_part.Resize(5);
      two_body_part_e.Resize(nelectrons,5);
      two_body_part_others.Resize(nelectrons,5);

      one_body_part_2.Resize(5);
      two_body_part_e_2.Resize(nelectrons,5);
      two_body_part_others_2.Resize(nelectrons,5);
    }
  private:
    friend class Jastrow2_wf;
    vector <Wavefunction *> wfObserver;
    Array1 <doublevar> one_body_part;
    Array2 <doublevar> two_body_part_e; //!< parts with this electron as first index
    Array2 <doublevar> two_body_part_others; //!< parts with this electron as second index
    Array1 <Array3<doublevar> > eibasis;
	
    Array1 <doublevar> one_body_part_2;
    Array2 <doublevar> two_body_part_e_2; //!< parts with this electron as first index
    Array2 <doublevar> two_body_part_others_2; //!< parts with this electron as second index
    Array1 <Array3<doublevar> > eibasis_2;
};

//######################################################################

class Jastrow2_wf : public Wavefunction {
public:

  void init(Wavefunction_data *);

  virtual void notify(change_type , int );
  virtual void updateVal(Wavefunction_data *, Sample_point *);
  virtual void updateLap(Wavefunction_data *, Sample_point *);
  virtual void updateForceBias(Wavefunction_data *, Sample_point *);

  virtual void getVal(Wavefunction_data *, int, Wf_return &);
  virtual void getLap(Wavefunction_data *, int, Wf_return &);

  virtual void getForceBias(Wavefunction_data *, int, Wf_return &);

  virtual void getDensity(Wavefunction_data *,int,  Array2 <doublevar> &);


  virtual void generateStorage(Wavefunction_storage *&);
  virtual void saveUpdate(Sample_point *, int e, Wavefunction_storage *);
  virtual void restoreUpdate(Sample_point *, int e, Wavefunction_storage *);

  virtual void saveUpdate(Sample_point *, int e1, int e2, Wavefunction_storage *);
  virtual void restoreUpdate(Sample_point *, int e1, int e2, Wavefunction_storage *);

  virtual void storeParmIndVal(Wavefunction_data *, Sample_point *,
                               int, Array1 <doublevar> & );
  virtual void getParmDepVal(Wavefunction_data *,
                             Sample_point *,
                             int,
                             Array1 <doublevar> &,
                             Wf_return &);
 
  virtual int getParmDeriv(Wavefunction_data *, Sample_point *,
                           Parm_deriv_return & );

  virtual void getSymmetricVal(Wavefunction_data *, 
			       int, 
			       Wf_return &);

  virtual void plot1DInternals(Array1 <doublevar> &,
			       vector <Array1 <doublevar> > &,
			       vector <string> &,
			       string );


  // JK: to implement analytical derivatives in BCS_wf (which uses two-body
  // Jastrow piece as a pair orbital), "electron-resolved" ParamDeriv is needed
  int get_twobody_ParmDeriv(Sample_point *, Array3<doublevar> & );
    
    
  void get_twobody(Array3 <doublevar>& twobody) { 
    twobody=two_body_save;
  }

  void get_onebody(Array3 <doublevar> & onebody) {
    onebody=one_body_ion;
  }

  void keep_ion_dep() { 
    keep_ion_dependent=1;
    one_body_ion.Resize(nelectrons, parent->natoms,5);
    one_body_ion=0.0;
    updateEverythingVal=1;
    updateEverythingLap=1;
    
  }

  void get_onebody_save(Array2 <doublevar> & onebody) { 
    onebody=one_body_save;
  }

  

private:

  //Control flags for lazy evaluation
  Array1 <int> electronIsStaleVal;
  Array1 <int> electronIsStaleLap;
  int updateEverythingVal;
  int updateEverythingLap;
  int sampleAttached;
  int dataAttached;
  int staticSample;
  int parmChanged;

  
  //stored variables
  Array2 <doublevar> one_body_save; 
  //!< save one-body parts(electron, valgradlap)

  Array3 <doublevar> two_body_save;
  //!< save two-body parts(see Jastrow_twobody_piece for explanation)
  doublevar u_twobody;
  //!< save for the total value of the two-body terms

  Array1 <  Array4 <doublevar> > eibasis_save; 
  //!< first array is group, 4d array is (electron, ion, basis#, valgradlap)
  


  //These are for backflow, which needs per-ion information
  int keep_ion_dependent;
  Array3 <doublevar> one_body_ion;
  //!< (electron,ion, valgradlap)


  Jastrow2_wf_data * parent;
  int nelectrons;
  
  int maxeibasis;
  int maxeebasis;

};
#endif //JASTROW2_WF_H_INCLUDED

//--------------------------------------------------------------------------
