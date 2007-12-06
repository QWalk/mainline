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

#ifndef JASTROW2_THREE_DIFFSPIN_H_INCLUDED
#define JASTROW2_THREE_DIFFSPIN_H_INCLUDED
#include "Qmc_std.h"
#include "Array45.h"
#include "Jastrow2_one.h"

/*!

\brief 
Three-body jastrow term


now, in this case, uij=uji, gradi uij  !=  gradj uij, and it's also not 
antisymmetric, so we have to have two vectors: one filled with gradi uij and
one with gradj uij.  Taking the derivative matrix as given by 

<pre>
  grad1    grad2     grad3    grad4

    0      g2u12     g3u13    g4u14
  g1u12      0       g3u23    g4u24
  g1u13    g2u23       0      g4u23
  g1u14    g2u24     g3u34      0
</pre>

we return the i'th column and i'th row.

 */
class Jastrow_threebody_piece_diffspin {
public:
  void set_up(vector <string> & words,
              System * sys //!< A list in order of center names
              );

  int writeinput(string &, ostream &);

  int showinfo(string & indent, ostream & os);

  /*!

  */

 
  void updateLap(int e,
                 const Array4 <doublevar> & eionbasis,
                 const Array3 <doublevar> & eebasis,
                 Array3 <doublevar> & lap);

  void updateVal(int e,
                 const Array4 <doublevar> & eionbasis,
                 const Array3 <doublevar> & eebasis,
                 Array1 <doublevar> & updated_val);

  void getParmDeriv(const Array3 <doublevar> & eionbasis, //i,at, basis
                    const Array3 <doublevar> & eebasis, // i,j,basis, with i<j
                    Parm_deriv_return & deriv);


  //start: added for ei back-flow
  void updateLap_E_I(int e,
		     const Array4 <doublevar> & eibasis,
		     const Array3 <doublevar> & eebasis,
		     Array3 <doublevar> & lap); 

  void updateVal_E_I(int e,
		     const Array4 <doublevar> & eibasis,
		     const Array3 <doublevar> & eebasis,
		     Array3 <doublevar> & lap);
  //end: added for ei back-flow
  
  /*!
    Number of parameters for an atom
  */
  int nparms(int at) {
    return _nparms(parm_centers(at));
  }
  int nparms();
  void getParms(Array1 <doublevar> & parms);
  void setParms(Array1 <doublevar> & parms);


  int eibasis_needed(int at) { 
    /*
    if(nparms(at) <=3) return 2;
    else if(nparms(at) <=7) return 3;
    else if(nparms(at) <=12) return 4;
    else return 4;
    */
    return eibasis_max(at);
  }
  int eebasis_needed() { 
    return eebasis_max;
		       
  }
  Jastrow_threebody_piece_diffspin():freeze(0),eebasis_max(0) {}
  
private:

  int make_default_list();

  //Array2 <doublevar> unique_parameters;
  Array1 < Array2 <doublevar> > unique_parameters_spin;
  Array1 <int> _nparms;        //!<Number of parameters in each row above
  vector <string> parm_labels;
  Array1 < Array2 <int> > linear_parms; //!< map the parameters onto a serial array
  Array1 <int> parm_centers; //!<Look up the parameter from the atom
  int freeze;

  Array1 <int> eibasis_max; //!< (atoms): number of basis functions used in the sum for each atom
  int eebasis_max; //number of ee basis functions used.
  int nspin_up;
  int nelectrons;
  
  
  Array2 <int> klm; //!< which basis functions to use for each parameter in the list(parm#, klm);
};


#endif //JASTROW2_THREE_DIFFSPIN_H_INCLUDED

//--------------------------------------------------------------------------
