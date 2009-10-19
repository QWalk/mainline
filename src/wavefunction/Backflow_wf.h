//--------------------------------------------------------------------------
// include/Backflow_wf.h
//
//
#ifndef BACKFLOW_WF_H_INCLUDED

#define BACKFLOW_WF_H_INCLUDED

#include "Qmc_std.h"
#include "Array.h"
#include "Wavefunction.h"
#include "Jastrow2_wf.h"
class Wavefunction_data;
class Backflow_wf_data;
class System;


class Backflow_wf_storage : public Wavefunction_storage
{
public:
  virtual ~Backflow_wf_storage()
  {}
private:
  friend class Backflow_wf;

   Array2 <doublevar> gradlap;
  Array2 <doublevar> detVal;

};


/*!
A slater wavefunction; \f$\Psi=\sum_i det_i(\Phi_1\Phi_2...)\f$
where the \f$\Phi\f$'s are one-particle molecular orbitals.
Also supports multiple states as a vector of wavefunction values.  Just
specify multiple STATE keywords.
*/
class Backflow_wf : public  Wavefunction
{

public:

  Backflow_wf()
  {}


  virtual int nfunc() {
    return 1;
  }


  virtual void notify(change_type , int );

  virtual void updateVal(Wavefunction_data *, Sample_point *);
  virtual void updateLap(Wavefunction_data *, Sample_point *);

  virtual void getVal(Wavefunction_data *, int, Wf_return &);
  virtual void getLap(Wavefunction_data *, int, Wf_return &);

  virtual void getDensity(Wavefunction_data *,int, Array2 <doublevar> &);

  virtual void saveUpdate(Sample_point *, int e, Wavefunction_storage *);
  virtual void restoreUpdate(Sample_point *, int e, Wavefunction_storage *);


  virtual void storeParmIndVal(Wavefunction_data *, Sample_point *,
                               int, Array1 <doublevar> & );
  virtual void getParmDepVal(Wavefunction_data *,
                             Sample_point *,
                             int,
                             Array1 <doublevar> &,
                             Wf_return &);

  virtual int getParmDeriv(Wavefunction_data *, 
			   Sample_point *,
			   Parm_deriv_return & );

  virtual void getSymmetricVal(Wavefunction_data *, 
			       int, 
			       Wf_return &);

  void generateStorage(Wavefunction_storage * & wfstore);


  void init(Wavefunction_data *);

  //--
private:


  void calcVal(Sample_point *);
  void calcLap(Sample_point *);
  void updateVal(int e, Sample_point *);
  Array1 <int> electronIsStaleVal;
  Array1 <int> electronIsStaleLap;
  int updateEverythingVal;
  int updateEverythingLap;

  int sampleAttached;
  int dataAttached;
  Backflow_wf_data * parent;

  //Saved variables for electron updates
  Array3 <doublevar>  moVal;//(electron,mo,[val grad hess])
  Array2 <doublevar> updatedMoVal;//(mo,[val grad hess])
  Array2 < Array2 <doublevar> > inverse;
  //!<inverse of the value part of the mo_values array transposed
  //(det,spin)(elec,elec)
  Array2 <doublevar> detVal;

  //Backflow-specific
  Array4 <doublevar> coor_grad; 
  //!< gradient of the coordinate transform (i,grad,alpha,beta)
  //where:
  // i is for the i'th quasiparticle coordinate
  // grad is the electron we're taking the gradient to
  // alpha is the direction of the gradient
  // beta is the direction of the quasiparticle

  Array3 <doublevar> coor_lap;
  //!< laplacian of the coordinate transform (i,grad,alpha)
  // where the indices are the same as for the gradient,
  // except alpha is the direction of the quasiparticle



  Array2 <doublevar> gradlap;
  //!< (elec,[grad lap])


  int staticSample;

  int nmo;        //!<Number of molecular orbitals
  int ndet;       //!<Number of determinants
  int ndim;       //!<Number of (spacial) dimensions each electron has
  Array1 <int> nelectrons; //!< 2 elements, for spin up and down
  Array1 <int> spin;       //!< lookup table for the spin of a given electron


  Jastrow2_wf jast;
};

#endif //SLAT_WF_H_INCLUDED
//--------------------------------------------------------------------------
