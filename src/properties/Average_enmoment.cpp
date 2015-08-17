/*
 
 Copyright (C) 2012 Lucas K. Wagner
 
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


#include "Average_enmoment.h"
#include <cmath>
#include "Properties_point.h"
//######################################################################


void Average_enmoment::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                       System * sys, Sample_point * sample, Properties_point & pt,
                       Average_return & avg) { 
  Average_return tmpavg;
  avg_gen->evaluate(wfdata,wf,sys,sample,pt,tmpavg);
  avg.type="enmoment_"+tmpavg.type;
  int navg=tmpavg.vals.GetDim(0);
  
  avg.vals.Resize(nmoment*(navg+1));
  Array1 <doublevar> enmom(nmoment);
  for(int i=0; i< nmoment; i++)  { 
    enmom(i)=pow(pt.energy(0),i);
    avg.vals(i)=enmom(i);
  }
  int count=nmoment;
  for(int i=0; i< nmoment; i++) {
    for(int j=0; j< navg; j++) {
      avg.vals(count)=enmom(i)*tmpavg.vals(j);
      count++;
    }
  }
}
//-----------------------------------------------------------------------------
void Average_enmoment::read(System * sys, Wavefunction_data * wfdata, vector
                   <string> & words) { 
  unsigned int pos=0;
  vector <string> avgsec;
  if(!readsection(words, pos=0,avgsec,"AVERAGE"))
    error("Need AVERAGE section in ENMOMENT");
  allocate(avgsec,sys,wfdata,avg_gen);

  if(readvalue(words,pos=0,nmoment, "NMOMENT"))
    nmoment++;
  else 
    nmoment=4;
}
//-----------------------------------------------------------------------------
void Average_enmoment::write_init(string & indent, ostream & os) { 
  os << indent << "ENMOMENT" << endl;
  os << indent << "AVERAGE { ";
  avg_gen->write_init(indent,os);
  os << indent << "}" << endl;
  os << indent << "NMOMENT " << nmoment << endl;
}
//-----------------------------------------------------------------------------
void Average_enmoment::read(vector <string> & words) { 
  vector<string> avgsec;
  unsigned int pos=0;
  if(!readsection(words, pos=0,avgsec,"AVERAGE"))
    error("Need AVERAGE section in ENMOMENT");
  allocate(avgsec,avg_gen);
  if(!readvalue(words, pos=0,nmoment,"NMOMENT"))
    error("need NMOMENT");
  
}
//-----------------------------------------------------------------------------
void Average_enmoment::write_summary(Average_return &avg ,Average_return & err,
    ostream & os) { 
  os << "energy moments" << endl;
  os << "Avg energy " << avg.vals[1] << " +/- " << err.vals[1] << endl;
  int navg=(avg.vals.GetDim(0)-nmoment)/nmoment;
  Average_return tmpavg, tmperr;
  tmpavg.vals.Resize(navg);
  tmperr.vals.Resize(navg);
  int count=nmoment;
  for(int i=0; i < nmoment; i++) {
    for(int j=0; j< navg; j++) {
      tmpavg.vals(j)=avg.vals(count);
      tmperr.vals(j)=err.vals(count);
      tmpavg.vals(j)/=avg.vals[i];
      tmperr.vals(j)/=avg.vals[i];
      
      count++;
    }
    os << "Average for moments of the local energy:  Moment " << i << endl;
    avg_gen->write_summary(tmpavg,tmperr,os);
  } 

   
}
//-----------------------------------------------------------------------------

