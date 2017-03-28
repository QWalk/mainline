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

#include "Wavefunction_data.h"
#include "Slat_wf_data.h"
#include "Program_options.h"
#include "Slat_Jastrow_data.h"
#include "Jastrow2_wf.h"
#include "Backflow_wf_data.h"
#include "Backflow_pf_wf_data.h"
#include "Pfaff_wf_data.h"
#include "BCS_wf_data.h"
#include "qmc_io.h"

int allocate(vector <string> & wftext, System * sys, Wavefunction_data * & wfptr)
{
  assert(wfptr == NULL);
  if(wftext.size() < 1) error("Empty wavefunction section");
  if(caseless_eq(wftext[0],"SLATER"))
    wfptr=new Slat_wf_data;
  else if(caseless_eq(wftext[0],"JASTROW2") 
	  || caseless_eq(wftext[0], "JASTROW")) 
    wfptr=new Jastrow2_wf_data;
  else if(caseless_eq(wftext[0],"SLATER-JASTROW"))
    wfptr=new Slat_Jastrow_data;
  else if(caseless_eq(wftext[0],"BACKFLOW"))
    wfptr=new Backflow_wf_data;
  else if(caseless_eq(wftext[0],"BACKFLOW-PFAFFIAN"))
    wfptr=new Backflow_pf_wf_data;
  else if(caseless_eq(wftext[0],"PFAFFIAN")) 
    wfptr=new Pfaff_wf_data;
  else if(caseless_eq(wftext[0],"BCS"))
    wfptr=new BCS_wf_data;
  else
    error("Error parsing the trial wavefunction section; unknown keyword ",
          wftext[0]);


  unsigned int pos=1;
  wfptr->read(wftext,pos,sys);
  return 1;
}

int deallocate(Wavefunction_data * & wfptr)
{
  if(wfptr == NULL)
    return 0;

  delete wfptr;
  wfptr=NULL;
  return 1;
}

Wavefunction_data * duplicate(Wavefunction_data * wfdata){
  return wfdata->clone();
}

//------------------------------------------------------------------------
