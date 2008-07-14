//------------------------------------------------------------------------
//include/Wavefunction_data.h

#ifndef PFAFF_WF_DATA_H_INCLUDED
#define PFAFF_WF_DATA_H_INCLUDED

#include "Qmc_std.h"
#include "Wavefunction_data.h"
#include "MO_matrix.h"
#include "Pfaff_wf.h"
#include "MatrixAlgebra.h"


//const int noccupied=4; //4 occupied orbitals for Mn 
void Renormalization(Array1 <doublevar> & coef, doublevar & norm, doublevar coef_eps);
void Getnormalization(Array1 <doublevar> & coef, doublevar & norm);

class Pfaff_wf_data : public Wavefunction_data
{
public:

  Pfaff_wf_data():molecorb(NULL) {}

  ~Pfaff_wf_data()
  {
    if(molecorb != NULL && is_my_mo) delete molecorb;
  }

  
  virtual int valSize()
  {
    if(optimize_pf) 
      return totoccupation(0).GetSize();
    else if(optimize_pfwt)
      return pfwt.GetSize();
    else return 2;
  }
  
  /*
  doublevar Triplet_to_Singlet_ratio( Array1 <doublevar> & triplet_coef,
                                      Array1 <doublevar> & singlet_coef,
                                      int nmo
                                      );
  */

  //void Renormalization( Array1 <doublevar> & coef,doublevar & norm );
  void allocate_pairing_orbital(vector <string> & pf_place, int orbint);
  void Pfaffian_optimize_read(vector <string> & pf_place, int pf, unsigned int pos);

  virtual void getVarParms(Array1 <doublevar> & parms);
  virtual void setVarParms(Array1 <doublevar> & parms);


  virtual void read(vector <string> & words,
                    unsigned int & pos,
                    System * sys
                   );
  virtual int supports(wf_support_type );

  void generateWavefunction(Wavefunction *&);

  int nparms();

  int showinfo(ostream & os);

  int writeinput(string &, ostream &);

  int useMo(MO_matrix * mo);
  
  doublevar coef_eps;
    
private:
  
  void init_mo();
  void getPFCoeff(Array1 <doublevar> & parms);
  void setPFCoeff(Array1 <doublevar> & parms);
  friend class Pfaff_wf;
  //vector <Wavefunction *> wfObserver;
  //!< The children of this _data

  //Array1 < Array1 <int> > occupation;
  //!< The place in totoccupation for the molecular orbitals for (fn, det, spin) Used to look up the right MO from MO_matrix.
  //Array1 < Array1 <int> > occupation_orig;
  //!< The molecular orbitals numbers for (function, determinant, spin)
  Array1 < Array1 <int> >  totoccupation; //!< all the molecular orbitals 
  Array1 < Array1 <int> >  occupation; //!< which orbs are occupied for each pfaffian
  Array1 < Array1 <int> >  occupation_pos; //!< position in totoccupation for each orb ccupied for each pfaffian
  Array1 < Array2 <int> >  order_in_pfaffian; //!< which pairing function should be used for which pfaffian 

  // Array1 < Array1 <int> >  nopairs_occupation; //!< molecular orbitals used for unpaired part of pfaffian matrix

  Array1 <int> spin;  //!< spin as a function of electron
  Array1 <int> opspin;//!< the opposing spin for an electron
  Array1 <int> rede;  //<! The number of the electron within its spin
  Array1 <int> nelectrons; //!< number of electrons as a function of spin
  Array1 <int> npairs; //!<number of triplet up,up pairs, number uf triplet down down 
  //and finally number of unpaired electrons spin up
  
  

  Array1 < Array1 <doublevar> > tripletorbuu; //!< coeficients for triplet pairing orbital
  Array1 < Array1 <doublevar> > tripletorbdd; //!< coeficients for triplet pairing orbital
  Array1 < Array1 <doublevar> > singletorb; //!< coeficients for singlet pairing orbital 
  Array1 < Array1 <doublevar> > unpairedorb; //!< coeficients for unpaired orbitals
  Array1 < Array1 <doublevar> > normalization; //!< coeficients for mixing Singlet and Unpaired w Triplet
  
  Array1 <doublevar> pfwt; //!< pfaffian weights
 
  int nmo;        //!<Number of molecular orbitals
  int npf;        //!<Number of pfaffians  
  int nsfunc;     //!<Number of pfaffian functions such singlets, triplets and etc separately
  int tote;       //!<Total number of electrons
  int ndim;       //!<Number of (spacial) dimensions each electron has
  Array1 < int > ntote_pairs; //!<Size of array for totoccupation
  Array1 < int > noccupied; //first n orbitals to be excluded from opt when VIRTUAL is set 
  

  int optimize_pf; //!< whether to optimize Pairing orbitals
  int optimize_pfwt; //!< whether to optimize pfaffian weights
  int check_pfwt_sign; //!< whether to check sign of weihght of multipfaffian
  Array1 < Array1 <string>  > optimize_string; //!<acctually counts all the coeficients to be optimized
  Array1 < Array1 < Array1 <int> > > optimize_total; //! optimize_total(pf)(uu,dd,singlet,unpaired,alpha)(poss)=1/0;
  //Array1 <int> optimize_bonds; //! optimize_bonds(....,1,...,2,....,2,...,1,...,3,3,3,3....) each ineteger connects one bond


  MO_matrix * molecorb;
  int is_my_mo; //!< Whether we allocated the MO_matrix
  


};

#endif //PFAFF_WF_DATA_H_INCLUDED
//------------------------------------------------------------------------
