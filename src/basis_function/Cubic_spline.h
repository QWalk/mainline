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

#ifndef CUBIC_SPLINE_H_INCLUDED
#define CUBIC_SPLINE_H_INCLUDED

#include "Basis_function.h"
#include "Spline_fitter.h"

/*!
This class represents a 1-D cubic spline interpolation.  It is meant
to be used for the radial part of an atomic orbital, and will return
symmetry-adjusted values for the calculation methods.

*/


class Cubic_spline: public Basis_function
{
private:

  //Overall symmetry type of a spline
  // added an extra term sym_7F_crystal   //CRYSTAL F orbital
  // added sym_9G_pyscf
  enum symmetry_type { sym_S, sym_P, sym_6D, sym_10F, sym_5D, sym_7F, sym_7F_crystal, sym_15G, sym_9G, sym_P_siesta, sym_D_siesta, sym_F_siesta, sym_9G_pyscf};
  
  string symmetry_lookup(symmetry_type );
  symmetry_type symmetry_lookup(string &);
  void assign_indiv_symmetries();
  int symmetry_lvalue(symmetry_type s);

  //symmetry type of the expansion(all functions.)
  enum indiv_symm_type {
    isym_S, isym_Px, isym_Py, isym_Pz,
    isym_Dxx, isym_Dyy, isym_Dzz, isym_Dxy, isym_Dxz, isym_Dyz,
    isym_Dz2r2, isym_Dx2y2,
    isym_Fxxx, isym_Fyyy, isym_Fzzz, isym_Fxxy, isym_Fxxz, //GAMESS spherical harmonics 
    isym_Fyyx, isym_Fyyz, isym_Fzzx, isym_Fzzy, isym_Fxyz,
    //these are linear combinations of spherical harmonics.  m2 is the imaginary part of m=2, 
    //p2 is the real part, etc.
    isym_Fm3, isym_Fm1, isym_F0, isym_Fp1, isym_Fp2, isym_Fp3, isym_Fp3mod,
    isym_Gxxxx, isym_Gyyyy, isym_Gzzzz, isym_Gxxxy, isym_Gxxxz,
    isym_Gyyyx, isym_Gyyyz, isym_Gzzzx, isym_Gzzzy, isym_Gxxyy,
    isym_Gxxzz, isym_Gyyzz, isym_Gxxyz, isym_Gyyxz, isym_Gzzxy,
    isym_G0,isym_G1,isym_G2,isym_G3,isym_G4,isym_G5,isym_G6,isym_G7,isym_G8         
};

  Array1 <symmetry_type> symmetry;
//!< The overall symmetry of a spline

  Array1 <indiv_symm_type> indiv_symmetry;
  //!< The orbital expansion of each function(S, Px, Py, etc..)

  Array1 <int> nfuncspline;
  //!< Number of functions represented by each spline

  doublevar threshold;
  //!< the upper limit of the spline

  int nsplines;
  //!< Number of splines this object represents

  double requested_cutoff;
  //!< If the user asked for a cutoff, this will have it.

  vector < vector < doublevar > > exponent;
  //!< Save the exponents for a gaussian expansion to print again

  vector < vector < doublevar > > coefficient;
  //!< Save the coefficients to reprint

  Array1 <doublevar> rcut;
  //!< actual cutoff for each function

  int nfunctions;
  //!<Number of functions that this class represents(taking into account 1 S, 3 P, etc)

  string atomname; //!< Name of the atom that this belongs to.
  string norm_type; //!< normalization type(CRYSTAL, GAMESS, etc.)
  bool renormalize; //!< Whether or not to renormalize the basis functions

  bool enforce_cusp; //!< whether or not to force the spline to have a cusp (only works for SPLINE inputs for now)
  bool match_sto; //!< whether to replace the basis function with a Slater-type orbital on the interval [o:rc]
  bool zero_derivative; //!< whether or not force zero derivative at the nucleus (strictly no cusp), if we want to treat the cusp in Jastrow
  double cusp;
  double cusp_matching; //!< the radius at which the cusp should take over from the rest of the basis function
  doublevar customspacing;

  Array1 <Spline_fitter> splines;

  void findCutoffs();

  /*!
    Read the spline fit points.  
  */
  int readspline(vector <string> & words);

public:

  Cubic_spline()
  {}
  ~Cubic_spline()
  {}

  //-----------------------------------------------------

  int showinfo(string & indent, ostream & os);
  int writeinput(string &, ostream &);

  virtual int read(
    vector <string> & words,
    //!< The words from the basis section that will create this basis function
    unsigned int & pos
    //!< The current position in the words(important if one basis section makes several functions); will be incremented as the Basis_function reads the words.
  );


  int readbasis(vector <string> & words, unsigned int & pos,
                Array1 <double> & parms) ;


  int nfunc();
  doublevar cutoff(int);
  virtual string label()
  {
    return atomname;
  }
  void raw_input(ifstream & input);

  void calcVal(const Array1 <doublevar> & r,
               Array1 <doublevar> & symvals,
               const int startfill=0);

  /*!
     All returned values are of the form \f$f(r),
     \frac{1}{r}\frac{df(r)}{dr}, \frac{d^2f(r)}{dr^2} \f$.
  */
  void calcLap(
    const Array1 <doublevar> & r,
    //!< in form r, r^2, x, y, z
    Array2 <doublevar> & symvals,
    //!< The values of the spline propogated through symmetry.  For example, a p state would be a 3x5 matrix, and an s state a 1x5.
    const int startfill=0
  );


  virtual void calcHessian(const Array1 <doublevar> & r,
			   Array2 <doublevar> & symvals,
			   //!< (func, [val,grad,d2f/dx2,d2f/dy2,d2f/dz2
			   //,d2f/dxdy,d2f/dxdz,d2f/dydz]) (nfunc,10)
			   const int startfill=0
			   );
};

#endif // CUBIC_SPLINE_H_INCLUDED
//--------------------------------------------------------------------------
