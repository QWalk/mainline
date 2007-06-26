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

//------------------------------------------------------------------------
//Cat_wf_data.h

#ifndef CAT_WF_DATA_H_INCLUDED
#define CAT_WF_DATA_H_INCLUDED

#include "Qmc_std.h"
#include "Wavefunction_data.h"

class Cat_wf;

class Cat_wf_data : public Wavefunction_data
{
public:

  virtual void read(vector <string> & words,
                    unsigned int & pos,
                    System * sys
                   );
  virtual int supports(wf_support_type );
  void getVarParms(Array1 <doublevar> & );
  void setVarParms(Array1 <doublevar> & );

  void attachObserver(Wavefunction * wf);
  void generateWavefunction(Wavefunction *&);
  int nparms();


  virtual int valSize();


  int showinfo(ostream & os);
  int writeinput(string &, ostream &);


  void renormalize()
  {
    for(int i=0; i< nUniqueFunc; i++)
    {
      data_array(i)->renormalize();
    }
  }

  void resetNormalization()
  {
    for(int i=0; i< nUniqueFunc; i++)
    {
      data_array(i)->resetNormalization();
    }
  }

  Cat_wf_data()
  {}
  ~Cat_wf_data()
  {
    for(int i=0; i< data_array.GetDim(0); i++)
      deallocate(data_array(i));
  }


private:
  friend class Cat_wf;
  int nFunc, nUniqueFunc;
  Array1 <int> ordering;
  Array1 <Wavefunction_data *> data_array;
};

#endif //CAT_WF_DATA_H_INCLUDED
//------------------------------------------------------------------------
