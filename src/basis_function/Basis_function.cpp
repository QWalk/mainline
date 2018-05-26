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


#include "Basis_function.h"
#include "qmc_io.h"
#include "Cubic_spline.h"
#include "Pade_function.h"
#include "Poly_pade_function.h"
#include "Gaussian_function.h"
#include "Exponent_cusp.h"
#include "Cutoff_cusp.h"
#include "Planewave_function.h"
#include "Rgaussian_function.h"
#include "Gen_pade_function.h"
#include "Group_function.h"
#include "Cosine_function.h"
#include "Spherical_bessel_function.h"	

int allocate(vector <string> & basistext, Basis_function * & bptr)
{
  assert( bptr == NULL );
  unsigned int pos=0;
  string type;
  readnext(basistext, pos, type);
  if(caseless_eq(type ,"AOSPLINE"))
    bptr= new Cubic_spline;
  else if(caseless_eq(type, "PADE"))
    bptr= new Pade_function;
  else if(caseless_eq(type, "POLYPADE"))
    bptr= new Poly_pade_function;
  else if(caseless_eq(type, "GAUSSIAN"))
    bptr= new Gaussian_function;
  else if(caseless_eq(type,"EXPONENTIAL_CUSP"))
    bptr= new Exponent_cusp;
  else if(caseless_eq(type,"CUTOFF_CUSP"))
    bptr=new Cutoff_cusp;
  else if(caseless_eq(type,"PLANEWAVE"))
    bptr=new Planewave_function;
  else if(caseless_eq(type,"RGAUSSIAN"))
    bptr=new Rgaussian_function;
  else if(caseless_eq(type,"GENPADE"))
    bptr=new Gen_pade_function;
  else if(caseless_eq(type,"BASIS_GROUPS"))
    bptr=new Group_function;
  else if(caseless_eq(type,"COSINE"))
    bptr=new Cosine_function;
  else if(caseless_eq(type,"SBESSEL"))
    bptr=new Spherical_bessel_function;
  else
    error("Didn't understand the basis section type ", type,
          " in basis section ", basistext[0]);

  bptr->read(basistext, pos);
  return 1;

}

int deallocate(Basis_function * & bptr)
{
  if(bptr == NULL)
    return 0;

  delete bptr;
  bptr=NULL;
  return 1;
}
//--------------------------------------------------------------------------
