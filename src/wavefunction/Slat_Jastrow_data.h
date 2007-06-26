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


#ifndef SLAT_JASTROW_DATA_H_INCLUDED
#define SLAT_JASTROW_DATA_H_INCLUDED

#include "Qmc_std.h"

#include "Wavefunction_data.h"

class Slat_Jastrow;

class Slat_Jastrow_data : public Wavefunction_data
{
public:

  void read(vector <string> & words,
            unsigned int & pos,
            System * sys
           );
  virtual int supports(wf_support_type );
  void generateWavefunction(Wavefunction *&);
  int valSize();
  void getVarParms(Array1 <doublevar> & );
  void setVarParms(Array1 <doublevar> & );

  void attachObserver(Wavefunction * wf);

  int nparms()
  {
    return slater->nparms()+jastrow->nparms();
  }



  int showinfo(ostream & os);
  int writeinput(string &, ostream &);

  void renormalize()
  {
    slater->renormalize();
    jastrow->renormalize();
  }
  void resetNormalization()
  {
    slater->resetNormalization();
    jastrow->resetNormalization();
  }

  Slat_Jastrow_data():slater(NULL), jastrow(NULL)
  {}
  ~Slat_Jastrow_data()
  {
    deallocate(slater);
    deallocate(jastrow);
  }

  virtual void clearObserver() {
    wfObserver.clear();
    slater->clearObserver();
    jastrow->clearObserver();
  }
  

private:
  friend class Slat_Jastrow;

  Wavefunction_data * slater;
  Wavefunction_data * jastrow;

};

#endif //SLAT_JASTROW_DATA_H_INCLUDED
//------------------------------------------------------------------------
