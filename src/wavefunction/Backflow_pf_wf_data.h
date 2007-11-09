//------------------------------------------------------------------------
//include/Backflow_pf_wf_data.h

#ifndef BACKFLOW_PF_WF_DATA_H_INCLUDED
#define BACKFLOW_PF_WF_DATA_H_INCLUDED

#include "Qmc_std.h"
#include "Wavefunction_data.h"
#include "MO_matrix.h"
#include "Backflow_pf_wf.h"
#include "Backflow_wf_data.h"
#include "Pfaff_wf_data.h"
#include "Jastrow2_wf.h"


class Pfaffian_keeper { 

 public:
  void read(vector <string> & words, System * sys );
  void getOccupation(Array1 <Array1 <int> > & totoccupation, 
		     int nmo);

  void writeinput(string & indent, ostream & os);
  void showinfo(ostream & os );


  int optimize_pf; //!< whether to optimize Pairing orbitals
  int optimize_pfwt; //!< whether to optimize pfaffian weights
  int npf;        //!<Number of pfaffians 
    
  void allocate_pairing_orbital(vector <string> & pf_place, int orbint);
  void Pfaffian_optimize_read(vector <string> & pf_place, int pf, unsigned int pos);

  void getVarParms(Array1 <doublevar> & parms);
  void setVarParms(Array1 <doublevar> & parms);
  int nparms();
  Array1 <int> nelectrons; //!< number of electrons as a function of spin
  Array1 < Array1 <int> >  totoccupation; //!< all the molecular orbitals 
  Array1 < Array1 <int> >  occupation; //!< which orbs are occupied for each pfaffian
  Array1 < Array1 <int> >  occupation_pos; //!< position in totoccupation for each orb occupied for each pfaffian
  Array1 < Array2 <int> >  order_in_pfaffian; //!< which pairing function should be used for which pfaffian 

  Array1 <int> npairs; //!<number of triplet up,up pairs, number uf triplet down down 
  //and finally number of unpaired electrons spin up
  Array1 < Array1 <doublevar> > tripletorbuu; //!< coeficients for triplet pairing orbital
  Array1 < Array1 <doublevar> > tripletorbdd; //!< coeficients for triplet pairing orbital
  Array1 < Array1 <doublevar> > singletorb; //!< coeficients for singlet pairing orbital 
  Array1 < Array1 <doublevar> > unpairedorb; //!< coeficients for unpaired orbitals
  Array1 < Array1 <doublevar> > normalization; //!< normalization of pairing functions
  Array1 <doublevar> pfwt; //!< pfaffian weight
  doublevar coef_eps;

private:
  void getPFCoeff(Array1 <doublevar> & parms);
  void setPFCoeff(Array1 <doublevar> & parms);
  friend class Pfaff_wf;
  
  
 
  int nmo;        //!<Number of molecular orbitals 
  int nsfunc;     //!<Number of pfaffian functions such singlets, triplets and etc separately
  int tote;       //!<Total number of electrons
  int ndim;       //!<Number of (spacial) dimensions each electron has
  Array1 < int > ntote_pairs; //!<Size of array for totoccupation
  Array1 < int > noccupied; //first n orbitals to be excluded from opt when VIRTUAL is set 
  

  
  int check_pfwt_sign; //!< whether to check sign of weihght of multipfaffian
  Array1 < Array1 <string>  > optimize_string; //!<acctually counts all the coeficients to be optimized
  Array1 < Array1 < Array1 <int> > > optimize_total; //! optimize_total(pf)(uu,dd,singlet,unpaired,alpha)(poss)=1/0;
  //Array1 <int> optimize_bonds; //! optimize_bonds(....,1,...,2,....,2,...,1,...,3,3,3,3....) each ineteger connects one bond

};


//######################################################################

/*!
 */
class Backflow_pf_wf_data : public Wavefunction_data
{
public:

  Backflow_pf_wf_data() {}

  ~Backflow_pf_wf_data()
  {
  }


  virtual int valSize() {
    // if(optimize_pf) 
    //  return totoccupation(0).GetSize();
    //else if(optimize_pfwt)
    //  return pfwt.GetSize();
    //else return 2;
    return 0;
  }

  virtual void getVarParms(Array1 <doublevar> & parms);
  virtual void setVarParms(Array1 <doublevar> & parms);

  virtual void read(vector <string> & words,
                    unsigned int & pos,
                    System * sys
                   );
  virtual int supports(wf_support_type );
  void generateWavefunction(Wavefunction *&);


  int showinfo(ostream & os);

  int writeinput(string &, ostream &);

  int nparms(){
    return bfwrapper.nparms()+pfkeeper.nparms();
  }
 

private:
  friend class Backflow_pf_wf;
  
  Backflow_wrapper bfwrapper;
  Pfaffian_keeper pfkeeper;

  string mo_place; //!< where to place the new mo's
  int optimize_pf; //!< whether to optimize determinant coefficients
  

};

#endif //BACKFLOW_PF_WF_DATA_H_INCLUDED
//------------------------------------------------------------------------
