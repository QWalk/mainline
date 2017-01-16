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


#ifndef MOLECULAR_SAMPLE_H_INCLUDED
#define MOLECULAR_SAMPLE_H_INCLUDED


#include "Sample_point.h"
class Wavefunction;
#include "Molecular_system.h"
class Sample_storage;


class Molecular_sample : public Sample_point
{
public:

  ~Molecular_sample()
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


  /*!
    Moves electron e
    position should be an array of length at least 3

  */
  void setElectronPos(const int e,const Array1 <doublevar> & position);
  virtual void setElectronPosNoNotify(const int e, 
                                    const Array1 <doublevar> & pos);

  void getElectronPos(const int e, Array1 <doublevar> & R);

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
  void updateECDist()
  {
    updateEIDist();
  }




  void getECDist(const int e, const int cent,
                 Array1 <doublevar> & distance)
  {
    getEIDist(e,cent, distance);
  }

  void getEIDist(const int e,const int ion, Array1 <doublevar> & distance)
  {
    assert( distance.GetDim(0) >= 5 );
    assert( ! ionDistStale(e));
    for(int i=0; i<5; i++)
    {
      distance(i)=iondist(i,ion,e);
    }
  }

  virtual int getEIDist_temp(const int e, int ion,
			      Array1 <doublevar> & distance) {
    assert(distance.GetDim(0)>=5);
    distance(1)=0;
    for(int d=0; d< 3; d++) {
      distance(d+2)=elecpos(e,d)-parent->ions.r(d,ion);
      distance(1)+=distance(d+2)*distance(d+2);
    }
    distance(0)=sqrt(distance(1));
    return 1;
  }

  void getEEDist(const int e1,const int e2, Array1 <doublevar> & distance)
  {
    //cout << "getEEDist" << endl;
    assert(distance.GetDim(0) >= 5);
    assert( ! elecDistStale(e1));
    assert( e1 < e2 );
    for(int i=0; i< 5; i++)
    {
      distance(i)=pointdist(i,e1,e2);
    }
  }

  void rawOutput(ostream &);
  void rawInput(istream &);

  void generateStorage(Sample_storage *& );
  void saveUpdate(int, Sample_storage * );
  void restoreUpdate(int, Sample_storage *);

private:

  int nelectrons;

  Array2 <doublevar> elecpos; //electron positions

  Array3 <doublevar> iondist;
  Array1 <int> elecDistStale;
  Array1 <int> ionDistStale;
  Array3 <doublevar> pointdist;  //this is an upper triangular matrix;
                                 //the lower is currently wasted

  Molecular_system * parent;

};



#endif //MOLECULAR_SAMPLE_H_INCLUDED
//-------------------------------------------------------------------------
