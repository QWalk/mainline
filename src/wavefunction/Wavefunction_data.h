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

#ifndef WAVEFUNCTION_DATA_H_INCLUDED
#define WAVEFUNCTION_DATA_H_INCLUDED

#include "Qmc_std.h"
#include "Array.h"
#include "Wavefunction.h"



enum wf_support_type { laplacian_update, density, parameter_derivatives };

/*!
\brief
Base class for wavefunctions..
 */
class Wavefunction_data
{
public:

  /*!
    Parse an input section in the first argument, starting from the
    second argument.
   */
  virtual void read(vector <string> & words,
                    unsigned int & pos,
                    System * sys
                   ) =0;


  virtual int supports(wf_support_type )=0;
  /*!
    Attach a wavefunction to be observed.  This object will
    now notify the wavefunction of any parameter changes.
   */
  virtual void attachObserver(Wavefunction * wf)
  {
    assert(wf != NULL);
    wfObserver.push_back(wf);
    wf->notify(data_attach,0);
  }

  /*!
    Get variational parameters.
   */
  virtual void getVarParms(Array1 <doublevar> &)=0;

  /*!
    Set variational parameters
   */
  virtual void setVarParms(Array1 <doublevar> &)=0;

  /*! 
   * Return true if the wave function in linear in the parameter,
   * false if not.
   * */
  virtual void linearParms(Array1 <bool> & is_linear) { 
    is_linear.Resize(nparms());
    is_linear=false;
  }

  /*!
    Number of variational parameters
   */
  virtual int nparms()=0;

  /*!
    The size of the array that is stored when doing
    Wavefunction::storeDepVal()
   */
  virtual int valSize()=0;

  /*!
    Print out a pretty output that is informative about the 
    wave function.
   */
  virtual int showinfo(ostream & os)
  {
    os << "This is the default showinfo.  The inherited "
    "class hasn't implemented it yet.\n";
    return 1;
  }


  /*!
    Use the wave function objects that this _data is linked to
    to find a normalization. 
   */
  virtual void renormalize()
  {}
  ;

  /*!
    Resets the normalization back to the original value.  
    Useful if a parameter has changed, and we want to restart
    the normalization procedure.
   */
  virtual void resetNormalization()
  {}
  ;

  /*!
    Write an input section with the current parameters for this 
    wave function.  Note: if the parameters don't change, 
    this can just be the original input echoed back.
    The first argument is the amount of indentation to 
    append to each line, and the second is the output stream.
   */
  virtual int writeinput(string &, ostream &)=0;

  virtual void generateWavefunction(Wavefunction * &)=0;

  virtual ~Wavefunction_data()
  {}

  virtual void clearObserver() {
    wfObserver.clear();
  }


protected:
  vector <Wavefunction *> wfObserver;
};




int allocate(vector <string> & wftext, System * sys,
             Wavefunction_data * & wfptr);
int deallocate(Wavefunction_data * & wfptr);

#endif //WAVEFUNCTION_DATA_H_INCLUDED
//------------------------------------------------------------------------
