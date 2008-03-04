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

#ifndef SPLINE_FITTER_H_INCLUDED
#define SPLINE_FITTER_H_INCLUDED

#include "Qmc_std.h"

class Spline_fitter {
 public:

  doublevar findCutoff();
  void splinefit(
    Array1 <doublevar>& x,
    //!< x interpolation points
    Array1 <doublevar>& y,
    //!< f(x) interpolation points
    double yp1,
    //!< df/dx at x_0
    double ypn
    //!< df/dx at x_n, with n=number of points
     );

  //Modify the spline so that the function goes smoothly to 
  //zero at 'cut'
  void enforceCutoff(doublevar cut);
  
  void raw_input(ifstream & input);
  
  int getInterval(doublevar r) {
    int interval=int(r*invspacing);
    //if(r >= threshold) cout << "r " << r << " threshold " << threshold << endl;
    assert(r < threshold);
    //cout << "interval " << interval << " spacing " << spacing << " r " << r << endl;
    //assert(interval*spacing <=r);
    //assert((interval+1)*spacing >= r);    
    return interval;
  }

  doublevar getVal(doublevar r, int i) {
    doublevar height=r-i*spacing;
    return coeff(i,0)+height*(coeff(i,1)
                              +height*(coeff(i,2)
                                       +height*coeff(i,3)));
  }


  void getDers(doublevar r, int i, 
               doublevar & func, doublevar & fdir,
               doublevar & f2dir) {
    //int i=int(r/spacing);
    doublevar height=r-i*spacing;    
    func=coeff(i,0)+height*(coeff(i,1)
                            +height*(coeff(i,2)
                                     +height*coeff(i,3)));
    
    fdir=(coeff(i,1)+height*(2*coeff(i,2)
                             +height*(3*coeff(i,3))))/r;
    
    f2dir=2*coeff(i,2)+6*height*coeff(i,3);    
  }

  int readspline(vector <string> & words, bool enforce_cusp, doublevar cusp);
  int writeinput(string & indent, ostream & os); 


  /*!
    \brief
    Whether this spline matches another..
   */
  int match(Spline_fitter & s) {

    //cout << "spacing " << spacing << " s " << s.spacing << endl;
    //cout << "threshold " << threshold << " s " << s.threshold << endl;
    if(fabs(spacing-s.spacing) > 1e-8)
      return 0;
    if(fabs(threshold-s.threshold) > 1e-8)
      return 0;


    return 1;
  }
  
  //pad with zeros out to thresh
  //if we already have support past thresh, does nothing.
  void pad(doublevar thresh);

 private:
    doublevar invspacing; //1/spacing, so we can do multiplications instead of divisions.
  doublevar spacing;
  doublevar threshold;
  Array2 <doublevar> coeff;
};
  

//----------------------------------------------------------------------

class Spline_fitter_nonuniform {
 public:

  doublevar findCutoff();
  void splinefit(
    Array1 <doublevar>& x,
    //!< x interpolation points
    Array1 <doublevar>& y,
    //!< f(x) interpolation points
    double yp1,
    //!< df/dx at x_0
    double ypn
    //!< df/dx at x_n, with n=number of points
     );

  
  int getInterval(doublevar r);


  doublevar getVal(doublevar r, int i) {
    doublevar height=r-pos(i);
    return coeff(i,0)+height*(coeff(i,1)
                              +height*(coeff(i,2)
                                       +height*coeff(i,3)));
  }


  void getDers(doublevar r, int i, 
               doublevar & func, doublevar & fdir,
               doublevar & f2dir) {
    //int i=int(r/spacing);
    doublevar height=r-pos(i);    
    func=coeff(i,0)+height*(coeff(i,1)
                            +height*(coeff(i,2)
                                     +height*coeff(i,3)));
    
    fdir=(coeff(i,1)+height*(2*coeff(i,2)
                             +height*(3*coeff(i,3))))/r;
    
    f2dir=2*coeff(i,2)+6*height*coeff(i,3);    
  }

  int readspline(vector <string> & words);


  /*!
    \brief
    Whether this spline matches another..
   */
  int match(Spline_fitter_nonuniform & s) {
    return 1;
  }

 private:
  doublevar threshold;
  Array1 <doublevar> pos;
  Array2 <doublevar> coeff;
};

#endif //SPLINE_FITTER_H_INCLUDED
