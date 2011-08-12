/*
 
 Copyright (C) 2011 Lucas K. Wagner
 
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

#ifndef AVERAGE_DENSITY_MATRIX_H_INCLUDED
#define AVERAGE_DENSITY_MATRIX_H_INCLUDED

#include "Average_generator.h"
#include "MO_matrix.h"

//Evaluate the two-body density matrix on a basis.
class Average_tbdm_basis:public Average_generator { 
 public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & );
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
  virtual ~Average_tbdm_basis() { 
    if(momat) delete momat;
  }
  Average_tbdm_basis() { 
    momat=NULL;
  }
 private:
  MO_matrix * momat;
  int nmo;
  int npoints_eval; //Number of integration points to use
};

#endif //AVERAGE_DENSITY_MATRIX_H_INCLUDED
