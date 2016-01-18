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

#ifndef AVERAGE_ENMOMENT_H_INCLUDED
#define AVERAGE_ENMOMENT_H_INCLUDED

#include "Average_generator.h"
#include "MO_matrix.h"

//Evaluate the two-body density matrix on a basis.
class Average_enmoment:public Average_generator { 
 public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & ) { 
    error("need to use Properties_point version of evaluate");
  }
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Properties_point & pt,
                        Average_return &avg);
  
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void randomize(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample) { 
    avg_gen->randomize(wfdata,wf,sys,sample);
  } 
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
  virtual ~Average_enmoment() {
    delete avg_gen;
  }
  Average_enmoment() { 
    avg_gen=NULL;
  }
 private:
  Average_generator * avg_gen;
  int nmoment;

};




#endif //AVERAGE_ENMOMENT_H_INCLUDED
