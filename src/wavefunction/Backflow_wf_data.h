//------------------------------------------------------------------------
//include/Backflow_wf_data.h

#ifndef BACKFLOW_WF_DATA_H_INCLUDED
#define BACKFLOW_WF_DATA_H_INCLUDED

#include "Qmc_std.h"
#include "Wavefunction_data.h"
#include "MO_matrix.h"
#include "Backflow_wf.h"
#include "Jastrow2_wf.h"


//######################################################################

class Backflow_wrapper {
 public:
  void updateVal(Sample_point * sample, Jastrow2_wf & jast,int e, 
		 int listnum, Array2 <doublevar> & newvals);
  //first index is mo#, then derivatives
  void updateLap(Sample_point * sample, Jastrow2_wf & jast,int e, 
		 int listnum, Array2 <doublevar> & newvals, //!<(mo,[val,grad,hess])
		 Array3 <doublevar>& coor_deriv, //!< (i,alpha,beta)
		 Array2 <doublevar> & coor_laplacian //!< (i,alpha)
		 );

  //start: added for ei back-flow
  //gives array of threebody_diffspin(i, a, d) i-electron index, a-atom index, d-dimension 
  //which is val/grad/lap w.r.t i-th electron coordinate of e-th electron & a-th atom three-body diff spin jastrow.
  void updateLapjastgroup(Sample_point * sample, int e, Array3 <doublevar> & threebody_diffspin);
  void updateValjastgroup(Sample_point * sample, int e, Array3 <doublevar> & threebody_diffspin);
  //end: added for ei back-flow

  
  void readOrbitals(System * sys, vector <string> & words);

  void init(System * sys, Array1 <Array1 <int> > & occupation,
	    vector <string> & words);
 
  int nparms() { 
    //start: added for ei back-flow
    if(has_electron_ion_bf) 
      return jdata.nparms()+jgroupdata.nparms();
    //end: added for ei back-flow
    else
      return jdata.nparms();
  }

  int nmo(){ return molecorb->getNmo(); } 
  void showinfo(ostream & os);
  void writeinput(string & indent, ostream & os);
  ~Backflow_wrapper() { 
    if(temp_samp != NULL) delete temp_samp;
    if(molecorb != NULL) delete molecorb;
  }
  Backflow_wrapper() { 
    temp_samp=NULL;
    molecorb=NULL;
  }

  void getNeighbors(Sample_point * sample,
		    Jastrow2_wf & jast,
		    int e, Array1 <int> & list,
		    int & nlist);


  Jastrow2_wf_data jdata;
  //start: added for ei back-flow
  Jastrow_group jgroupdata;
  int has_electron_ion_bf;
  //end: added for ei back-flow
 private:
  Sample_point * temp_samp;
  MO_matrix * molecorb;

  
};
void backflow_config(Sample_point * sample, 
		     int e,
		     Array3 <doublevar> & corr,
		     Array3 <doublevar> & onebody,
		     Array3 <doublevar> & threebody_diffspin,
		     Sample_point * temp_samp);

//######################################################################

class Determinant_keeper { 

 public:

  void read(vector <string> & words, System * sys );
  void getOccupation(Array1 <Array1 <int> > & totoccupation, 
		     //!< all the molecular orbitals for a given spin
		     Array2 <Array1 <int> > & occupation
		     //!< The place in totoccupation for the molecular orbitals for (det, spin) Used to look up the right MO from MO_matrix after we cull the non-calculated orbitals
		     );

  int nparms() { 
    if(optimize_det) return detwt.GetDim(0);
    else return 0;
  }

  void writeinput(string & indent, ostream & os);
  void showinfo(ostream & os );


  Array2 < Array1 <int> > occupation_orig;
  //!< The molecular orbitals numbers for (determinant, spin) as entered by the user
  Array1 <doublevar> detwt;

  Array1 <int> spin;  //!< spin as a function of electron
  Array1 <int> opspin;//!< the opposing spin for an electron
  Array1 <int> rede;  //<! The number of the electron within its spin
  Array1 <int> nelectrons; //!< number of electrons as a function of spin

  //int nmo;        //!<Number of molecular orbitals
  int ndet;       //!<Number of determinants
  //int ndim;       //!<Number of (spacial) dimensions each electron has

  int optimize_det;
};


//######################################################################

/*!
 */
class Backflow_wf_data : public Wavefunction_data
{
public:

  Backflow_wf_data() {}

  ~Backflow_wf_data()
  {
  }


  virtual int valSize() {
    //error("redo valsize()");
    //if(optimize_backflow) return 0;
    //if(optimize_det) return 2*dkeeper.detwt.GetDim(0);
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
    return bfwrapper.nparms()+dkeeper.nparms();
  }

private:
  friend class Backflow_wf;
  
  Backflow_wrapper bfwrapper;
  Determinant_keeper dkeeper;

  Array2 < Array1 <int> > occupation;
  string mo_place; //!< where to place the new mo's
  int optimize_det; //!< whether to optimize determinant coefficients

};

#endif //BACKFLOW_WF_DATA_H_INCLUDED
//------------------------------------------------------------------------
