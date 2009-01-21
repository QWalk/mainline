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


#ifndef PERIODIC_SAMPLE_H_INCLUDED
#define PERIODIC_SAMPLE_H_INCLUDED


#include "Sample_point.h"
#include "Periodic_system.h"
class Wavefunction;

class Sample_storage;

/*!
 
*/
class Periodic_sample : public Sample_point
{
public:

  Periodic_sample() { overall_sign=1.0; overall_phase=0.0; } 
  ~Periodic_sample()
  {}

  void init(System * sys);

  void randomGuess();

  int electronSize()
  {
    return nelectrons;
  }
  int ionSize()
  {
    return parent->ions.size();
  }

  int centerSize()
  {
    return ionSize();
  }
  doublevar getIonCharge(const int ion)
  {
    return parent->ions.charge(ion);
  }

  virtual int getBounds(Array2 <doublevar> & latvec, Array1 <doublevar> & origin) {
    latvec=parent->latVec;
    origin=parent->origin;
    return 1;
  }
  /*!
    Moves electron e
    position should be an array of length at least 3

  */
  void setElectronPos(const int e,const Array1 <doublevar> & position);
  void translateElectron(const int e, const Array1 <doublevar> & trans);
  
  void getElectronPos(const int e, Array1 <doublevar> & R)
  {
    assert( R.GetDim(0) >= 3 );
    for(int i=0; i< 3; i++)
    {
      //R(i) = points.r(i,e);
      R(i) = elecpos(e,i);
    }
  }
  void getAllElectronPos(Array2 <doublevar> & pos) {
    assert(pos.GetDim(0) >= nelectrons);
    assert(pos.GetDim(1) >= 3);

    pos=elecpos;
  }

  void moveIon(const int ion, const Array1 <doublevar> & r);
  
  void getIonPos(const int ion, Array1 <doublevar> & r)
  {
    assert(r.GetDim(0) >= 3);

    for(int i=0; i<3; i++)
    {
      r(i)=parent->ions.r(i,ion);
    }
  }



  void updateEIDist();
  void updateEEDist();


  void updateECDist();
  //{
  //  updateEIDist();
  //}

  void getECDist(const int e, const int cent,
                 Array1 <doublevar> & distance)
  {
    assert(!cenDistStale(e));
    assert(distance.GetDim(0) >=5);
    for(int i=0; i< 5; i++)
      distance(i)=cendist(e, cent, i);
  }

  void getEIDist(const int e,const int ion, Array1 <doublevar> & distance)
  {
    assert( distance.GetDim(0) >= 5 );
    assert( ! ionDistStale(e));
    for(int i=0; i<5; i++)
    {
      distance(i)=iondist(ion,e,i);
    }
  }
  /*!
  Returns the vector pointing from e1 to e2.
  */
  void getEEDist(const int e1,const int e2, Array1 <doublevar> & distance)
  {
    assert(distance.GetDim(0) >= 5);
    assert( ! elecDistStale(e1));
    assert( e1 < e2 );
    for(int i=0; i< 5; i++)
    {
      distance(i)=pointdist(e1,e2,i);
    }
  }

  void rawOutput(ostream &);
  void rawInput(istream &);

  void generateStorage(Sample_storage *& );
  void saveUpdate(int, Sample_storage * );
  void restoreUpdate(int, Sample_storage *);

  /*!
  Note that this is not bulletproof; it only works when translateElectron() is used,
  so it doesn't work when you start doing complicated things with setElectronPos(),
  etc.  It's enough for the sampling method, and that's about it.
   */
  doublevar overallSign() { return overall_sign; }
  doublevar overallPhase() { return overall_phase; }
private:

  int nelectrons;

  Array2 <doublevar> elecpos; //electron positions

  Array3 <doublevar> cendist; //!< center distances
  Array1 <int> cenDistStale;

  Array3 <doublevar> iondist;
  Array1 <int> elecDistStale;
  Array1 <int> ionDistStale;
  Array3 <doublevar> pointdist;  //this is an upper triangular matrix;
  //the lower is currently wasted
  
  doublevar overall_sign;
  doublevar overall_phase;
  // false for complex-valued wavefunctions, i.e., for non-integer k-points
  bool update_overall_sign;

  Periodic_system * parent;     //The System that created this object
};



#endif //PERIODIC_SAMPLE_H_INCLUDED
//-------------------------------------------------------------------------
