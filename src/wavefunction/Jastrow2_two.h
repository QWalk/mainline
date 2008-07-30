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

#ifndef JASTROW2_TWO_H_INCLUDED
#define JASTROW2_TWO_H_INCLUDED
#include "Qmc_std.h"
#include "Array45.h"

#include "Wavefunction.h"
class System;

/*!

<b> Computational description </b>

For the value, imagine we have U= sum_i sum_j>i u_ij
So, for 4 electrons, for example, we are summing over the matrix
<pre>
0   u12   u13    u14
0    0    u23    u24
0    0     0     u34
0    0     0      0
</pre>
Note that the u's must be symmetric with respect to  the vector r_ij.
This implies that the gradient is antisymmetric for only electron-electron
interactions.  We can write the
matrix as:
<pre>
  grad1    grad2     grad3    grad4

    0      g2u12     g3u13    g4u14
  g1u12      0       g3u23    g4u24
  g1u13    g2u23       0      g4u23
  g1u14    g2u24     g3u34      0
</pre>
(where gi means operation with the gradient with respect to coordinate i)
The gradient is given by a sum down the column, and since we have
antisymmetry, gi uij = - gj uij

The laplacian is similar, except li uij = lj uij, so it is symmetric.

So, for our updates, we take an array that is (electron, function, valgradlap)
the electron index is as follows(for electron number j)
<pre>
1j
2j
...
jj(ignored)
j(j+1)
j(j+2)
...
j(nelectrons)
</pre>

For value, this does not matter, since it is symmetric, but for the gradient,
it is important to get the ordering(=signs) right.

We return in the same ordering as the basis functions, so for the gradient,
we return

<pre>
g1u1j
g2u2j
...
0
gjuj(j+1)
...
gjuj(nelectrons)
</pre>

*/

class Jastrow_twobody_piece {
public:
  virtual void set_up(vector <string> & words, System * sys);
  virtual int writeinput(string &, ostream &);
  virtual int showinfo(string &, ostream &);

  virtual int nparms() {if(freeze) return 0;  return parameters.GetDim(0);}
  virtual void getParms(Array1 <doublevar> & parms);
  virtual void setParms(Array1 <doublevar> & parms);

  /*! Number of basis functions we need */
  virtual int nbasis_needed() { return parameters.GetDim(0); }


  virtual void updateLap(int e,
                 const Array3 <doublevar> & eebasis,
                 Array2 <doublevar> & lap);


  virtual void updateVal(int e,
                 const Array3 <doublevar> & eebasis,
                 Array1 <doublevar> & val);

  virtual void getParmDeriv(const Array3 <doublevar> & eebasis, //expects in form i,j,basis, with i<j
                            Parm_deriv_return &);

  virtual void parameterSaveVal(int e,
                        const Array3 <doublevar> & eebasis,
                        Array1 <doublevar> & save, int begin_fill);

  virtual void parameterSaveLap(const Array4 <doublevar> & eebasis,
                        Array2 <doublevar> & save_lap, int begin_fill);

  /*! \brief
    these come in handy when we need parameters visible, e.g. in
    Jastrow2_wf::plot1DInternals
  */
  virtual void unfreeze() { freeze=0; }
  virtual void refreeze() { freeze=1; }

  virtual ~Jastrow_twobody_piece() {}
private:
  Array1 <doublevar>  parameters;
  int nelectrons;
  int freeze;
};

/*!
\brief two-body for different spins
 */
class Jastrow_twobody_piece_diffspin:public Jastrow_twobody_piece {
public:
  virtual void set_up(vector <string> & words, System * sys);
  virtual int writeinput(string &, ostream &);
  virtual int showinfo(string&, ostream &);
  /*! Number of variational parameters */
  virtual int nparms() { if(freeze) return 0; else return 2*spin_parms.GetDim(1);}

  /*! Number of basis functions we need */
  virtual int nbasis_needed() { return spin_parms.GetDim(1); }
  virtual void getParms(Array1 <doublevar> & parms);
  virtual void setParms(Array1 <doublevar> & parms);


  virtual void updateLap(int e,
                 const Array3 <doublevar> & eebasis,
                 Array2 <doublevar> & lap);


  virtual void updateVal(int e,
                 const Array3 <doublevar> & eebasis,
                 Array1 <doublevar> & val);

  virtual void getParmDeriv(const Array3 <doublevar> & eebasis, //expects in form i,j,basis, with i<j
                            Parm_deriv_return &);

  virtual void parameterSaveVal(int e,
                        const Array3 <doublevar> & eebasis,
                        Array1 <doublevar> & save, int begin_fill);

  virtual void parameterSaveLap(const Array4 <doublevar> & eebasis,
                        Array2 <doublevar> & save_lap, int begin_fill);

  /*! \brief
    these come in handy when we need parameters visible, e.g. in
    Jastrow2_wf::plot1DInternals
  */
  virtual void unfreeze() { freeze=0; }
  virtual void refreeze() { freeze=1; }

private:
  Array2 <doublevar> spin_parms;
  int nspin_up;
  int freeze;
};



#endif //JASTROW2_TWO_H_INCLUDED

//--------------------------------------------------------------------------
