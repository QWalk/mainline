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

//--------------------------------------------------------------------------
// Cat_wf.h
//

#ifndef CAT_WF_H_INCLUDED
#define CAT_WF_H_INCLUDED

#include "Wavefunction.h"
#include "Cat_wf_data.h"

class Cat_wf_storage : public Wavefunction_storage
{
  Cat_wf_storage()
  {}
  ~Cat_wf_storage()
  {
    for(int i=0; i< store_array.GetDim(0); i++)
    {
      if(store_array(i))
        delete store_array(i);
    }
  }
private:
  friend class Cat_wf;
  Array1 <Wavefunction_storage *> store_array;

};

/*!
 
*/
class Cat_wf : public Wavefunction
{
public:

  Cat_wf()
  {}
  ~Cat_wf()
  {
    for(int i=0; i< wfarray.GetDim(0); i++)
      deallocate(wfarray(i));
  }



  void generateStorage(Wavefunction_storage * & wfstore);


  int nfunc()
  {
    return nFunc;
  }

  virtual void notify(change_type , int );


  virtual void updateVal(Wavefunction_data *, Sample_point *);
  virtual void updateLap(Wavefunction_data *, Sample_point *);
  virtual void updateForceBias(Wavefunction_data *, Sample_point *);

  virtual void getVal(Wavefunction_data *, int, Wf_return &,
                      int startpos=0);
  virtual void getLap(Wavefunction_data *, int, Wf_return &,
                      int startpos=0) ;
  virtual void getForceBias(Wavefunction_data *, int, Wf_return &,
                            int startpos=0);

  virtual void getDensity(Wavefunction_data *,int,  Array2 <doublevar> &);

  virtual void saveUpdate(Sample_point *, int e, Wavefunction_storage *);

  virtual void restoreUpdate(Sample_point *, int e, Wavefunction_storage *);


  virtual void storeParmIndVal(Wavefunction_data *, Sample_point *,
                               int, Array1 <doublevar> & );
  virtual void getParmDepVal(Wavefunction_data *,
                             Sample_point *,
                             int,
                             Array1 <doublevar> &,
                             Wf_return  &);


private:
  friend class Cat_wf_data;
  Array1 <Wavefunction *> wfarray;
  Array1 <int> ordering;
  int nFunc;
  int nUniqueFunc;
};

#endif //CAT_WF_H_INCLUDED
//--------------------------------------------------------------------------
