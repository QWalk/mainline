/*
 
Copyright (C) 2007 Lucas K. Wagner, 2008 Jindrich Kolorenc

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

#ifndef BCS_WF_DATA_H_INCLUDED
#define BCS_WF_DATA_H_INCLUDED

#include "Qmc_std.h"
#include "Wavefunction_data.h"
#include "BCS_wf.h"
#include "Jastrow2_wf.h"
#include "System.h"

//######################################################################

class BCS_jastrow_cofactor { 
 public:
  virtual void valGradLap(const Array3 <doublevar> & twobody,
			  const Array2 <doublevar> & onebody,
			  int i, int j,
			  Array2 <doublevar> & valgradlap)=0;
  virtual ~BCS_jastrow_cofactor() { }
};

class BCS_jastrow_u:public BCS_jastrow_cofactor { 
 public:
  virtual void valGradLap(const Array3 <doublevar> & twobody,
			  const Array2 <doublevar> & onebody,
			  int i, int j,
			  Array2 <doublevar> & valgradlap) { 
    valgradlap.Resize(2,5);
    //valgradlap(0,0)=twobody(i,j,0)+onebody(i,0)*onebody(j,0);
    valgradlap(0,0)=twobody(i,j,0);
    for(int d=1; d< 5; d++) { 
      //valgradlap(0,d)=twobody(i,j,d)+onebody(i,d)*onebody(j,0);
      //valgradlap(1,d)=twobody(j,i,d)+onebody(i,0)*onebody(j,d);
      valgradlap(0,d)=twobody(i,j,d);
      valgradlap(1,d)=twobody(j,i,d);
    }
  }
};

//######################################################################


/*!
 */
class BCS_wf_data : public Wavefunction_data
{
public:

  BCS_wf_data() {
    //  molecorb=NULL;
  }

  ~BCS_wf_data()
  {
    //if(molecorb) delete molecorb;
  }


  virtual int valSize() {
    return 0;
  }

  virtual void getVarParms(Array1 <doublevar> & parms);
  virtual void setVarParms(Array1 <doublevar> & parms);

  virtual void read(vector <string> & words,
                    unsigned int & pos,
                    System * sys
                   );
  virtual int supports(wf_support_type );
  void generateWavefunction(Wavefunction *&);


  int showinfo(ostream & os);

  int writeinput(string &, ostream &);

  int nparms(){
    int tot=jastdata.nparms();
    return tot;
  }

private:
  friend class BCS_wf;
  
  Jastrow2_wf_data jastdata;

  //MO_matrix * molecorb;

  Array1 <int> nelectrons; //!< 2 elements, for spin up and down
  Array1 <int> spin;       //!< lookup table for the spin of a given electron
  Array1 <int> rede;

  doublevar magnification_factor;

};

#endif //BCS_WF_DATA_H_INCLUDED
//------------------------------------------------------------------------
