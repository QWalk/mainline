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
#include "Optimize_method2.h"
#include "Newton_opt_method.h"
#include "Dmc_method.h"
#include "Test_method.h"
#include "Plot_method.h"
#include "Localize_method.h"
#include "Nodes_method.h"
#include "Reptation_method.h"
#include "Postprocess_method.h"

int allocate(vector <string> & words,
             Program_options & options,
             Qmc_method * & methptr)
{
  assert(methptr == NULL);

  if(caseless_eq(words[0],"VMC"))
    methptr=new Vmc_method;
  
  else if(caseless_eq(words[0],"OPTIMIZE"))
    methptr=new Optimize_method;

  else if(caseless_eq(words[0],"OPTIMIZE2"))
    methptr=new Optimize_method2;
  
  else if(caseless_eq(words[0],"NEWTON_OPT"))
    methptr=new Newton_opt_method;
  
  else if(caseless_eq(words[0],"DMC"))
    methptr=new Dmc_method;
  
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
