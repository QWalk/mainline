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
//------------------------------------------------------------------------
//include/Jastrow_wf_data.h

#ifndef JASTROW_WF_DATA_H_INCLUDED
#define JASTROW_WF_DATA_H_INCLUDED

#include "Qmc_std.h"
#include "Wavefunction_data.h"
#include "Basis_function.h"

class Jastrow_wf;

class Jastrow_wf_data : public Wavefunction_data
{
public:

  Jastrow_wf_data():elecIonBasis(NULL), elecElecBasis(NULL)
  {}

  virtual void read(vector <string> & words,
                    unsigned int & pos,
                    System * sys
                   );
                   
  virtual int supports(wf_support_type );
  void generateWavefunction(Wavefunction *&);

  virtual int valSize();
  void getVarParms(Array1 <doublevar> & );
  void setVarParms(Array1 <doublevar> & );
  int nparms()
  {
    return cuspBasis(0)->nparms()+cuspBasis(1)->nparms()+ eeCorrelation.GetDim(0)
           +eiCorrelation.GetDim(0) +eeiCorrelation.GetDim(0)+nbasisparms;
  }
  int showinfo(ostream & os);
  int writeinput(string &, ostream &);

  void renormalize();
  void resetNormalization();


  void calcEECuspVal(int i, int j, Array1 <doublevar> & eedist);

  //Helper functions for calcLap()
  void calcEECuspLap(int i, int j,Array1 <doublevar> & eedist);
  void sumElecIon(Array1 <doublevar> &,Array2 <doublevar> &);
  void sumElecElec(Array2 <doublevar> & , Array2 <doublevar> &);


  /*!
    Get the contribution from one electron to the value.
   */
  void updateVal(Array2 < Array1 <doublevar> > & a_kval, //!< The a_k's that we need to update
                 Sample_point * sample, //!< sample point we're using
                 int e ,                //!< electron to update
                 Array1 <doublevar> & ionPartialSum,
                 //!< contribution from this electron to the value
                 Array2 <doublevar> & elecPartialSum

                );



  /*!
    Fill the parameter independent values.
   */
  void fillParmInd(Sample_point * sample,
                   int e,
                   Array1 <doublevar> & vals
                  );
  /*!
    Use the values from fillParmInd to update the total value.
   */
  void updateParmsVal(Sample_point * sample,
                      int e, doublevar & value,
                      Array1 <doublevar> & oldvalue);
  /*!
    makes the saved variables for a static calculation.
   */
  void makeStaticSave(Sample_point * sample,
                      Array1 <doublevar> & values,
                      Array2 <doublevar> & valPartial,
                      Array3 <doublevar> & derivatives);
  /*!
    uses the variables from makeStaticSave to update value, valPartialSum, derivatives
   */
  void updateParms(Sample_point * sample,
                   Array1 <doublevar> & oldvalue,
                   Array2 <doublevar> & oldValPartial,
                   Array3 <doublevar> & oldderiv,
                   doublevar & value,
                   Array1 <doublevar> & valPartialSum,
                   Array2 <doublevar> & derivatives);

private:
  friend class Jastrow_wf;

  int optimize_basis;
  //!< Whether or not to optimize the basis
  int nbasisparms;
  //!< Number of parameters in the basis

  //! the \f$ \gamma \f$ in \f$ \frac{-c}{\gamma} e^{-\gamma r_{ij} \f$
  //Array1 <doublevar> gamma;

  //! the \f$ c \f$ in \f$ \frac{-c}{\gamma} e^{-\gamma r_{ij} \f$.  Can be diff for diff spins
  //Array1 <doublevar> cusp;

  //c_klm's
  Array1 <doublevar> eeCorrelation;
  Array1 <doublevar> eiCorrelation;
  Array1 <doublevar> eeiCorrelation;

  //indexing to the wavefunction parameters.
  Array1 <int> eiCorr_start, eiCorr_end;
  Array1 <int> eeiCorr_start, eeiCorr_end;
  Array1 <int> c_k, c_l, c_m;


  //Calculation info

  Array1 <Basis_function * > cuspBasis;

  Basis_function * elecIonBasis;
  Basis_function * elecElecBasis;
  doublevar maxeecutoff; //!< maximum electron-electron cutoff
  doublevar maxeicutoff; //!< maximum electron-ion cutoff

  int nelectrons;
  int nions;
  int eenum;
  int einum;
  int nUniqueAtoms;
  int basisoffset; //!<Number of basis functions we have defined as a constant
  Array1 <int> firste, laste;
  Array1 <int> spin;
  doublevar faccp;
  doublevar facco;
  doublevar normalization;

  vector <string> atomname;


  //Temporary storage for the basis function values.
  Array2 < Array2 <doublevar> > a_k;
  //Array2 <Array1 <doublevar> > a_kval;


  //! (electron#, electron#) for outer, (basis, valgradlap) for inner
  Array2 <Array2 <doublevar> > b_m;
  Array3 <doublevar> eecusp;

  Array1 <Array1 <doublevar> > bupdate;

  Array1 <doublevar> elecPartial_temp;

};

#endif //JASTROW_WF_DATA_H_INCLUDED
//------------------------------------------------------------------------
