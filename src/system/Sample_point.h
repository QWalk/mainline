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


#ifndef SAMPLE_POINT_H_INCLUDED
#define SAMPLE_POINT_H_INCLUDED

#include "Qmc_std.h"
class Wavefunction;
class Sample_storage;
class System;


/*!
\brief
Manages the coordinates and distances of a single point in 3N
dimensional space.

Whenever a coordinate is moved, it keeps track of whether the
distances have been updated.  When you request them to be updated,
it checks, and only does it if necessary.

All distances are returned in the form \f$ [r, r^2, x, y,z] \f$

*/
class Sample_point
{
public:

  Sample_point():wfObserver(NULL)
  {}
  virtual ~Sample_point()
  {}
  ;

  //Set up
  /*!
    \brief
    Initialize with a parent
  */
  virtual void init(System * parent) {
    error("init isn't implemented for this type of sample_point.");
  }

  /*!
    \brief
    Randomly place the electrons
  */
  virtual void randomGuess()=0;

  //General information

  /*!
    \brief
    Number of electrons
  */
  virtual int electronSize()=0;

  /*!
    Number of ions
  */
  virtual int ionSize()=0;
  virtual int centerSize()=0;

  virtual int ndim();
  //Position and property functions

  virtual int getBounds(Array2 <doublevar> & latvec, Array1 <doublevar> & origin) {
    return 0;
  }
  
  /*!
    \brief
    Set an electron's position to pos
  */
  virtual void setElectronPos(const int e, const Array1 <doublevar> & pos)=0;

  /*!
    \brief
     Move an electron by trans.  This should be used when there is a concept of 
     an electron path, since it will correctly change the sign with regards to 
     kpoints, and possibly other instances where some translational symmetry is
     broken.
   */
  virtual void translateElectron(const int e, const Array1 <doublevar> & trans) {
    Array1 <doublevar> pos(3);
    getElectronPos(e,pos);
    for(int d=0; d< 3; d++) 
      pos(d)+=trans(d);
    setElectronPos(e,pos);
  }
   
  /*!
  return any prefactor due to k-point sampling or some such thing
   */
  virtual doublevar overallSign() {
    return 1.0;
  }        
  
  
  /*!
    \brief
    Dangerous!  move an electron without notifying observing wave functions
   */
  virtual void setElectronPosNoNotify(const int e, const Array1 <doublevar> & pos)
    { error("setElectronPosNoNotify not implemented");}

  /*!
    \brief
    Get an electron position
  */
  virtual void getElectronPos(const int e, Array1 <doublevar> & pos)=0;

  /*!
    \brief
    Copy all the electron positions to the argument
  */
  virtual void getAllElectronPos(Array2 <doublevar> & pos) {
    error("getAllElectronPos not implemented for this sample point");
  }

  /*!
    \brief
    Change an ion's position
  */
  virtual void moveIon(const int ion, const Array1 <doublevar> & pos)=0;

  /*!
    \brief
    Get an ion's position
  */
  virtual void getIonPos(const int ion, Array1 <doublevar> & pos)=0;

  /*!
    \brief
    Get an ion's charge(in a.u.)
  */
  virtual doublevar getIonCharge(const int number)=0;


  //Distance functions

  //Electron-electron

  /*!
    \brief
    Update the electron-electron distances
  */
  virtual void updateEEDist()=0;

  /*!
    \brief
    Returns the vector pointing from e1 to e2(e2 must be greater
    than e1)
    
    Each Sample_point implementation should have an assert(e2 > e1) within 
    the function body.
  */
  virtual void getEEDist(const int e1, const int e2,
                         Array1 <doublevar> & distance)=0;

  //Electron-ion
  /*!
    \brief
    Update the electron-ion distances
  */
  virtual void updateEIDist()=0;

  /*!
    \brief
    Get the vector pointing from the ion to the electron
  */
  virtual void getEIDist(const int e, const int ion,
                         Array1 <doublevar> & distance)=0;

  /*! 
    \brief
    Get the vector from the ion to the electron.

    This is a speed optimization; it does not update the
    state of the Sample_point; it just returns the distance
    for a single ion.  This should be used only when 
    you need just one distance from an ion, and you're
    going to move the electron again soon.
    It's also optional for implementations, so it will
    return 0 if it's not supported.(used in pseudopotential
    when doing optimization of the wave function)
  */
  virtual int getEIDist_temp(const int e, int ion,
			      Array1 <doublevar> & distance) {
    return 0;
  }

  //Electron-center
  virtual void updateECDist()=0;
  virtual void getECDist(const int e, const int cent,
                         Array1 <doublevar> & distance)=0;


  //I/O

  /*!
    \brief
    Save the coordinates to a file
  */
  virtual void rawOutput(ostream &)=0;

  /*!
    Invert rawOutput(), and read coordinates from a file
  */
  virtual void rawInput(istream &)=0;


  virtual void attachObserver(Wavefunction * wfptr);

  //Update functions

  virtual void generateStorage(Sample_storage *&)=0;
  virtual void saveUpdate(int e, Sample_storage *)=0;
  virtual void restoreUpdate(int e, Sample_storage *)=0;

protected:
  Wavefunction * wfObserver;
};


//--------------------------------------------------


class Sample_storage
{
private:
  friend class Sample_point;
  friend class Molecular_sample;
  friend class Periodic_sample;
  friend class Ring_sample;
  friend class SHO_sample;
  friend class HEG_sample;
  Array2 <doublevar> iondist_temp;
  Array2 <doublevar> pointdist_temp;
  Array1 <doublevar> pos_temp;
};

//int restore_configs(Array1 <Sample_point *> & sample_array,
 //                   string & startconfig);

//int initialize_configs(vector <string> words, Array1 <Sample_point * > & sample_array, Program_options & options);


/*!
  returns 1 if it read into the sample and 0 if not
 */
int read_config(string & last_read, istream & is, 
                Sample_point * sample);

void write_config(ostream & os, 
                  Sample_point * sample);



class Config_save_point {
 public:
  void savePos(Sample_point * sample);
  void restorePos(Sample_point * sample);
  void mpiSend(int node);
  void mpiReceive(int node);
  
 private:
  Array1 <Array1 <doublevar> > electronpos;
};
#endif
//-------------------------------------------------------------------------
