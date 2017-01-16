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

#include "Array.h"
#include "Qmc_std.h"



/*!
\brief
Averaging with weights
 */
void average(int start, //!<start average at point #
             int end,   //!< stop average at point #
             Array2 <doublevar> &, //!< values to average, points second index
             Array2 <doublevar> &, //!< weights to average, points second index
             Array1 <doublevar> &, //!< return average
             Array1 <doublevar> &, //!< return variance
             int paravg=0          //!< whether to sync across process or not
            );

/*!
\brief
  Averaging with weights for a single vector
 */
void average(int start,
             int end,
             Array1 <doublevar> &,
             Array1 <doublevar> &,
             doublevar & avg,
             doublevar & var,
             int paravg=0
            );

/*!
\brief
Averaging without weights
 */
void average(int start,
             int end,
             Array2 <doublevar> & , //!<values
             Array1 <doublevar> & , //!< average
             Array1 <doublevar> &,  //!< variance
             int paravg=0           //!< sync across processes
            );

/*!
\brief
Get autocorrelation time of a vector

 */
void autocorrelation(int start,
		     int end,
		     const Array2 <doublevar> &, //!< values
		     const Array1 <doublevar> &, //!< average in
		     const Array1 <doublevar> &, //!< variance in
		     Array2 <doublevar> &, //!< autocorrelation in (wf,n-1 steps)
		     int autocorr_depth,
		     int paravg=0 //!< sync across processes
);

//----------------------------------------------------------------------
