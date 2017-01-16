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
#if 0
#ifndef CBASIS_FUNCTION_H_INCLUDED
#define CBASIS_FUNCTION_H_INCLUDED

#include "Array.h"
#include "Qmc_std.h"

/*!
This class represents a basis set around a single center.  
*/
class CBasis_function
{
public:


  virtual ~CBasis_function()
  {}

  /*!
    \brief
    Return only the value of the function in symvals
  */
  virtual void calcVal(const Array1 <doublevar> & r,
                       //!< in form r, r^2, x, y, z
                       Array1 <dcomplex> & symvals,
                       //!< The values and derivatives of the basis function
                       const int startfill=0
                       //!< The place in the array to start the fill
                       )=0;

  /*!
    \brief
    Returns the value and derivatives of the basis functions in the
    form [function, [val, df/dx, df/dy, df/dz, \f$\nabla^2f \f$] ]
   */
  virtual void calcLap(
    const Array1 <doublevar> & r,
    //!< in form r, r^2, x, y, z
    Array2 <dcomplex> & symvals,
    //!< The values and derivatives of the basis function
    const int startfill=0
    //!< The place in the array to start the fill
  ) =0;


  virtual void calcHessian(const Array1 <doublevar> & r,
			   Array2 <dcomplex> & symvals,
			   //!< (func, [val,grad,d2f/dx2,d2f/dy2,d2f/dz2
			   //,d2f/dxdy,d2f/dxdz,d2f/dydz]) (nfunc,10)
			   const int startfill=0
			   ) { 
    error("This function doesn't support Hessians");
  }


  virtual int read(
    vector <string> & words,
    //!< The words from the basis section that will create this basis function
    unsigned int & pos
    //!< The current position in the words(should be obseleted, since one section -> one set); will be incremented as the Basis_function reads the words.
  )=0;


  /*!
    \brief
    Return the number of functions that this represents
  */
  virtual int nfunc()=0;

  /*!
    \brief
    Return the label of the center that this function is around
  */
  virtual string label()=0;

  /*!
    \brief
    Return the cutoff of a given function, at which point it is
    safe to assume the function and its derivatives are zero.
  */
  virtual doublevar cutoff(int )=0;

  /*!
    \brief
    Show some summary information for the output file.  Return
    one on success.
  */
  virtual int showinfo(string & indent, ostream & os)
  {
    os << indent 
       <<  "Showing basis function info..this basis function doesn't"
      " have the showinfo function\n";
    return 1;
  }

  /*!
    \brief
    Write a valid input section for rereading later.  Note that
    this also serves as some documentation for the section.
  */
  virtual int writeinput(string &, ostream &)=0;

  /*!
    \brief
    put the set of variational parameters into the argument
  */
  virtual void getVarParms(Array1 <doublevar> & parms) {
    parms.Resize(0);
  }

  /*!
    \brief
     After changing the variational parameters, set them within
     the basis function.

     NB: the argument of this function must be exactly the same
     length as getVarParms()
  */
  virtual void setVarParms(Array1 <doublevar> & parms) {
    assert(parms.GetDim(0)==0);
  }

  /*!
    \brief
    Number of parameters this function has.
  */
  virtual int nparms() {
    return 0;
  }
};


#endif //CBASIS_FUNCTION_H_INCLUDED
//--------------------------------------------------------------------------
#endif
