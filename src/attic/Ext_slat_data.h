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
//include/Ext_slat_data.h

#ifndef EXT_SLAT_DATA_H_INCLUDED
#define EXT_SLAT_DATA_H_INCLUDED

#include "Qmc_std.h"
#include "Wavefunction_data.h"
#include "MO_matrix.h"
#include "Ext_slat.h"

/*!


<li>
<b> ORBITALS </b> This section is an input deck for an MO_matrix
</li>
</ul>

PAIRS

 */
class Ext_slat_data : public Wavefunction_data
{
public:

  Ext_slat_data():molecorb(NULL) {}

  ~Ext_slat_data()
  {
    if(molecorb != NULL ) delete molecorb;
  }


  virtual int valSize()
  {
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


  int nparms()
  {

    //We have nelectrons coefficients and nelectrons amplitudes
    return amplitude.GetDim(1);
  }

private:

  void propagate_irreducible();

  friend class Ext_slat;
  MO_matrix * molecorb;

  Array3 <int> pairs; //!< which pairs to combine
  Array3 <doublevar> pair_coeff; //!< their coefficients(NB: c_0=sqrt(1-c_1^2)
  Array2 <doublevar> amplitude; //!< amplitude of the MO
  
  int nmo;
  Array1 <int> spin;

};

#endif //SLAT_WF_DATA_H_INCLUDED
//------------------------------------------------------------------------
