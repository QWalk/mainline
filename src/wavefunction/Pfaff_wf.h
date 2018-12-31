//--------------------------------------------------------------------------
// include/Pfaff_wf.h
//
//
#ifndef PFAFF_WF_H_INCLUDED

#define PFAFF_WF_H_INCLUDED

#include "Qmc_std.h"
#include "Array.h"
#include "Wavefunction.h"
class Wavefunction_data;
class Pfaff_wf_data;
class System;




void UpdatePfaffianRowLap( Array1 < Array1 <doublevar> > & mopfaff_row, 
                           int & e,  
                           const Array3 <doublevar> &  moVal, 
                           const Array1 < Array1 <int> > & occupation_pos,
                           const Array1 <int> & npairs, 
                           const Array2 <int> & order_in_pfaffian,
                           const Array1 < Array1 <doublevar> > & tripletorbuu,
			   const Array1 < Array1 <doublevar> > & tripletorbdd,
                           const Array1 < Array1 <doublevar> > & singletorb,
                           const Array1 < Array1 <doublevar> > & unpairedorb,
                           const Array1 < Array1 <doublevar> > & alpha_angle,
                           doublevar coef_eps
                           );

void UpdatePfaffianRowHess( Array1 < Array2 <doublevar> > & mopfaff_row, 
			    int & e,  
			    const Array3 <doublevar> &  moVal, 
			    const Array1 < Array1 <int> > & occupation_pos,
			    const Array1 <int> & npairs, 
			    const Array2 <int> & order_in_pfaffian,
			    const Array1 < Array1 <doublevar> > & tripletorbuu,
			    const Array1 < Array1 <doublevar> > & tripletorbdd,
			    const Array1 < Array1 <doublevar> > & singletorb,
			    const Array1 < Array1 <doublevar> > & unpairedorb,
			    const Array1 < Array1 <doublevar> > & alpha_angle,
                            doublevar coef_eps
			    );

void FillPfaffianMatrix( Array2 <doublevar> & mopfaff_tot,  
                         const Array3 <doublevar> &  moVal,
                         const Array1 < Array1 <int> > & occupation_pos,
                         const Array1 <int> & npairs, 
                         const Array2 <int> & order_in_pfaffian,
                         const Array1 < Array1 <doublevar> > & tripletorbuu,
			 const Array1 < Array1 <doublevar> > & tripletorbdd,
                         const Array1 < Array1 <doublevar> > & singletorb,
                         const Array1 < Array1 <doublevar> > & unpairedorb,
                         const Array1 < Array1 <doublevar> > & alpha_angle,
                         doublevar coef_eps
                         );


/*!
A Pfaffian wavefunction; \f$\Psi=Pf(\Phi(1,2)\Phi(3,4)...)\f$
where the \f$\Phi(i,j)\f$'s are two-particle pairing orbitals.
*/
class Pfaff_wf : public  Wavefunction
{

public:

  Pfaff_wf()
  {}

  virtual int nfunc()
  {
    return 1;
  }

  virtual void notify(change_type , int );

  virtual void updateVal(Sample_point * const);
  virtual void updateLap(Sample_point * const);

  virtual void getVal(Wf_return &);
  virtual void getLap(int, Wf_return &);



  virtual int getParmDeriv(Wavefunction_data *, 
			   Sample_point *,
			   Parm_deriv_return & );

  virtual void getSymmetricVal(Wavefunction_data *, 
			       int, 
			       Wf_return &);


  void init(Wavefunction_data *);
  doublevar coef_eps;

  //--
private:

 // void save_for_static();

  void calcVal(Sample_point * const );
  void updateVal(Sample_point * const, int);
  void calcLap(Sample_point * const);
  void updateLap(Sample_point * const, int);

  Array1 <doublevar> electronIsStaleVal;
  Array1 <doublevar> electronIsStaleLap;
  int updateEverythingVal;
  int updateEverythingLap;
  int updateMoLap;

  int sampleAttached;
  int dataAttached;
  Pfaff_wf_data * parent;

  //Saved variables for electron updates
  Array3 <doublevar> moVal;

  Array2 <doublevar> updatedMoVal;

  Array1 < Array2 <doublevar> > inverse;
  //!<inverse of the value part of the mo_values array transposed

  Array1 <doublevar> pfaffVal;

  //Variables for a static(electrons not moving) calculation
  Array3 <doublevar> saved_laplacian;
  //!<Saved laplacian for a static calculation (electron, function, [val grad lap])

  Array1 < Array2 <doublevar> > mopfaff_tot;



  int nmo;        //!<Number of molecular orbitals
  int npf;        //!<Number of pfaffians
  int nsfunc;     //!<Number of pfaffian functions such singlets, triplets and etc separately
  int ndim;       //!<Number of (spacial) dimensions each electron has
  int npairs;     //!<Number of all pairs icluding no pairs (also total dim of pfaffian matrix)
  bool firstime;
  
  
  Array1 <int> nelectrons; //!< 2 elements, for spin up and down
  
};

#endif //PFAFF_WF_H_INCLUDED
//------------------------------------------------------------

