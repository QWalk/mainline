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

#ifndef JASTROW2_ONE_H_INCLUDED
#define JASTROW2_ONE_H_INCLUDED
#include "Qmc_std.h"
#include "Array45.h"
#include "Wavefunction.h"
void get_onebody_parms(vector<string> & words, vector<string> & atomnames,
                        Array2 <doublevar> & unique_parameters,
                        Array1 <int> & nparms,
                        vector <string> & parm_labels,
                        Array2 <int> & linear_parms,
                        Array1 <int> & parm_centers);

/*!

\brief 
One-body jastrow term

Basis format:
Array of electron-ion basis : (ion, basis, valgradlap)  

 */
class Jastrow_onebody_piece {
public:
  void set_up(vector <string> & words,
              vector <string> & atomnames //!< A list in order of center names
              );

  int writeinput(string &, ostream &);

  int showinfo(string & indent, ostream & os);

  /*!

  */
  void updateLap(int e,
                 const Array3 <doublevar> & eionbasis,
                 Array1 <doublevar> & lap);


  void updateLap_ion(int e, const Array3 <doublevar> & eionbasis,
			Array3 <doublevar> & lap);
  /*!

  */
  void updateVal(int e,
                 const Array3 <doublevar> & eionbasis,
                 doublevar & updated_val);

  void getParmDeriv(int e, const Array3 <doublevar> & eionbasis,
                    Parm_deriv_return & );
  /*!
    Number of parameters for an atom
  */
  int nparms(int at) {
    return _nparms(parm_centers(at));
  }
  int nparms();
  void getParms(Array1 <doublevar> & parms);
  void setParms(Array1 <doublevar> & parms);

                  
  /*!
    updated_val = sum(parm(i) * save(i))
  */
  void parameterSaveVal(int e,
                        const Array3 <doublevar> & eionbasis,
                        Array1 <doublevar> & save, int begin_fill);
  //void parameterUpdateVal(Array1 <doublevar> & save, doublevar & val, int begin_fill);


  /*!
  take all the electron basis functions and sum over them,
  returning only the parameter-dependent parts
  total_lap(j) = sum(parm(i)*lap(i,j))
  */
  void parameterSaveLap(const Array3 <doublevar> & eionbasis,
                        Array2 <doublevar> & save_lap, int begin_fill);
  //void parameterCalcLap(const Array2 <doublevar> & save_lap,
  //                      Array2 <doublevar> & lap, int begin_fill);


  /*! \brief
    these come in handy when we need parameters visible, e.g. in
    Jastrow2_wf::plot1DInternals
  */
  void unfreeze() { freeze=0; }
  void refreeze() { freeze=1; }

  
  /*! \brief
    need to export this, so that Jastrow2_wf::plot1DInternals can map
    from "uniquie parameters" to atomic basis
  */
  int atom_kind(int at) { return parm_centers(at); }

  Jastrow_onebody_piece():freeze(0) { }
private:
  Array2 <doublevar> unique_parameters;
  Array1 <int> _nparms;        //!<Number of parameters in each row above
  vector <string> parm_labels;
  Array2 <int> linear_parms; //!< map the parameters onto a serial array
  Array1 <int> parm_centers; //!<Look up the parameter from the atom
  int freeze;

};


#endif //JASTROW2_ONE_H_INCLUDED

//--------------------------------------------------------------------------
