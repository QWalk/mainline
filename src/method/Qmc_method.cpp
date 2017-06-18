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
#include "Qmc_method.h"
#include "Vmc_method.h"
#include "Optimize_method.h"
#include "Shdmc_method.h"
#include "Dmc_method.h"
#include "Rndmc_method.h"
#include "Test_method.h"
#include "Plot_method.h"
#include "Localize_method.h"
#include "Nodes_method.h"
#include "Reptation_method.h"
#include "Postprocess_method.h"
#include "Plot1d_method.h"
#include "Wannier_method.h"
#include "Linear_optimization.h"
#include "Rotate_orbs.h"
#include "Lowdin_method.h"
#include "Maximize_method.h"

int allocate(vector <string> & words,
             Program_options & options,
             Qmc_method * & methptr)
{
  assert(methptr == NULL);

  if(caseless_eq(words[0],"VMC"))
    methptr=new Vmc_method;
  
  else if(caseless_eq(words[0],"OPTIMIZE"))
    methptr=new Optimize_method;
  
  else if(caseless_eq(words[0],"SHDMC"))
    methptr=new Shdmc_method;
  
  else if(caseless_eq(words[0],"DMC"))
    methptr=new Dmc_method;

  else if(caseless_eq(words[0],"RNDMC"))
    methptr=new Rndmc_method;
    
  else if(caseless_eq(words[0],"TEST"))
    methptr=new Test_method;
  
  else if(caseless_eq(words[0],"PLOT"))
    methptr=new Plot_method;

   else if(caseless_eq(words[0],"LOCALIZE"))
    methptr=new Localize_method;
 
  else if(caseless_eq(words[0],"NODES"))
    methptr=new Nodes_method;
  
  else if(caseless_eq(words[0],"REPTATION"))
    methptr=new Reptation_method;
  
  else if(caseless_eq(words[0],"POSTPROCESS"))
    methptr=new Postprocess_method;

  else if(caseless_eq(words[0],"PLOT1D"))
    methptr=new Plot1D_method;

  else if(caseless_eq(words[0],"WANNIER"))
    methptr=new Wannier_method;
  else if(caseless_eq(words[0],"LINEAR"))
    methptr=new Linear_optimization_method;
  else if(caseless_eq(words[0],"ROTATE_ORBS"))
    methptr=new Rotate_orbs_method;
  else if(caseless_eq(words[0],"LOWDIN"))
    methptr=new Lowdin_method;
  else if(caseless_eq(words[0],"MAXIMIZE"))
    methptr=new Maximize_method;
  
  
  else
    error("Error parsing the method section; unknown keyword ",
          words[0]);


  unsigned int pos=1;

  methptr->read(words,pos,options);
  return 1;
}

int allocate(vector <string> & words,
             Program_options & options,
             Qmc_avg_method * & methptr)
{
  assert(methptr == NULL);
  if(words[0] == "VMC") {
     methptr=new Vmc_method;
  }
  else if(words[0] == "DMC") {
     methptr=new Dmc_method;
  }
  else if(words[0] == "RNDMC") {
     methptr=new Rndmc_method;
  }
  else error(" ",words[0], " is not a valid averaging method");

  unsigned int pos=1;

  methptr->read(words,pos,options);
  return 1;
}

		


int deallocate(Qmc_method * & methptr)
{
  if(methptr == NULL)
    return 0;

  delete methptr;

  methptr=NULL;
  return 1;
}
