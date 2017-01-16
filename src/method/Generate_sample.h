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


#ifndef GENERATE_SAMPLE_H_INCLUDED
#define GENERATE_SAMPLE_H_INCLUDED


#include "Split_sample.h"

//The objective of this function is to generate a converged set of samples as 
//efficiently as possible.
//It will use the provided system objects to generate nconfig samples, stored in 
//the configs variable.

void generate_sample(Sample_point * sample,
    Wavefunction * wf,
    Wavefunction_data * wfdata,
    Guiding_function * guidewf,
    int nconfig,
    Array1 <Config_save_point> & configs);
#include "MO_matrix.h"

void generate_mo_sample(Sample_point * sample, System * sys,
    MO_matrix * mo, int list, int nconfig, Array1 <Array1 <doublevar> > & r);

#endif //GENERATE_SAMPLE_H_INCLUDED
