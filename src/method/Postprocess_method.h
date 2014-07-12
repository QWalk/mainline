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


#ifndef POSTPROCESS_METHOD_H_INCLUDED
#define POSTPROCESS_METHOD_H_INCLUDED

#include <vector>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include "Qmc_method.h"
#include "qmc_io.h"
#include "Array.h"
#include "MO_matrix.h"
#include "Sample_point.h"
#include "Properties.h"
class System;
class Wavefunction;
class Program_options;

/*!
In DMC, use the option SAVE_TRACE filename to save the walker positions and their weights to a binary file (called 'filename') once every block.
These walkers should be decorrelated, so they can now be treated independently.

Now, you can use 
method { POSTPROCESS
  READCONFIG filename
  AVERAGE { ... } 
  DENSITY { ... }
  }
  include sys
  trialfunc { .. }

to take any averages that you'd like. 

There are some tradeoffs with this method. First, since we only save the walkers occasionally, the error bars will b
e larger than if you had put AVERAGE into the DMC section directly. However, if the AVERAGE takes some time to evalu
ate, such as the 2-RDM, then you will also not waste time evaluating the object on highly correlated walkers. You sh
ould use this scheme when one of the following is true:
a) You're not sure what averaging variables you want when the DMC calculation is run.
b) The averaging variables you want dominate the cost of the DMC calculation (2-RDM!)

Otherwise, evaluation on the fly (AVERAGE in the DMC section) is better.
 
 */
class Postprocess_method : public Qmc_method
{
public:
  void read(vector <string> words,
            unsigned int & pos,
            Program_options & options);
  void run(Program_options & options, ostream & output);
  int showinfo(ostream & os);
  Postprocess_method():sys(NULL),pseudo(NULL),wfdata(NULL) {}
  ~Postprocess_method()
  {
    if(sys != NULL) delete sys;
    if(pseudo!=NULL) delete pseudo;
    if(wfdata!=NULL) delete wfdata;
    for(int i=0; i< densplt.GetDim(0); i++) delete densplt(i);
    for(int i=0; i< average_var.GetDim(0); i++) delete average_var(i);
  }

private:
  int gen_point(Wavefunction * wf, Sample_point * sample,
    Config_save_point & configpos, doublevar weight, Properties_point & pt);
  int worker(Wavefunction * wf, Sample_point * sample);
  int master(Wavefunction * wf, Sample_point * sample,FILE * f,ostream & );
  
  
  
  System * sys;
  Pseudopotential * pseudo;
  Wavefunction_data * wfdata;
  Array1 < Local_density_accumulator *> densplt;
  Array1 < Average_generator * > average_var;
  string configfile;
};

//----------------------------------------------------------------------


inline void weighted_update_average(doublevar weight, doublevar xnew,doublevar &  xavg, doublevar &  xvar, doublevar & totweight)
  {
  doublevar nweight=totweight+weight;
  doublevar navg=(weight*xnew+totweight*xavg)/nweight;
  doublevar nvar=(weight*xnew*xnew+totweight*(xavg*xavg+xvar))/nweight-navg*navg;
  xavg=navg;
  xvar=nvar;
}
//----------------------------------------------------------------------

class Postprocess_average { 
public:
  Postprocess_average(int navg) { 
    avgavg.Resize(navg);
    varavg.Resize(navg);
    avgen=varen=avgwt=varwt=0.0;
    npoints=0;
  }


  void update_average(Properties_point & pt) { 
   // cout << "here" << endl;
    weighted_update_average(pt.weight(0),pt.energy(0),avgen,varen,avgwt);
    for(int i=0; i< avgavg.GetDim(0); i++) { 
      //cout << "resizing avg " <<i << " avgavg " << avgavg.GetDim(0) 
       // << " varavg " << varavg.GetDim(0) << " pt.avgrets "<< pt.avgrets.GetDim(0) << " " 
       // << pt.avgrets.GetDim(1) << endl;
      if(avgavg(i).vals.GetDim(0)==0) { 
        avgavg(i)=pt.avgrets(0,i);
        varavg(i)=pt.avgrets(0,i);
        avgavg(i).vals=0;
        varavg(i).vals=0;
      }
      //cout << "updaing average " << endl;
      for(int j=0; j< pt.avgrets(0,i).vals.GetDim(0); j++) { 
        weighted_update_average(pt.weight(0),pt.avgrets(0,i).vals(j),avgavg(i).vals(j),varavg(i).vals(j),avgwt);
      }
    }
    avgwt+=pt.weight(0);
    npoints++;
  }

  void print(Array1<Average_generator *> average_var,ostream & output) { 
    output << "Averages" << endl;
    output << "total_energy " << avgen << " +/- " << sqrt(varen/npoints) << " (sigma " << sqrt(varen) << ") "<< endl;
    output << "weight " << avgwt/npoints << endl;
    for(int i=0; i< avgavg.GetDim(0); i++) { 
      for(int j=0; j< varavg(i).vals.GetDim(0); j++) 
        varavg(i).vals(j)=sqrt(varavg(i).vals(j)/npoints);
      average_var(i)->write_summary(avgavg(i),varavg(i),output);
  }
    
  }
private:
  doublevar avgen,varen,avgwt,varwt;
  long int npoints;
  Array1 <Average_return> avgavg, varavg;
  
};

#endif //POSTPROCESS_METHOD_H_INCLUDED
//------------------------------------------------------------------------
