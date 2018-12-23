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

//
#ifndef AVERAGE_GENERATOR_H_INCLUDED
#define AVERAGE_GENERATOR_H_INCLUDED
#include "Qmc_std.h"
#include "Wavefunction.h"
#include "Wavefunction_data.h"
#include "Basis_function.h"
#include "System.h"
#include "Basis_function.h"
struct Properties_point;
struct Average_return {
  string type;
  Array1 <doublevar> vals;
};

/*!
 The average generators are meant to be averages which are not *too* big(ie, < 1000 doubles), and should
 go into the Properties_point structure.  They will be averaged with error bars, and can do post-processing
 either in the output of gosling or the main code(they use the same functions, so it's automatic).
 
 The interface gives the code writer a lot of power, including the power to mess up the VMC/RMC/DMC calculation.
 There are thus a few rules of engagement.  
 1) Leave the state of sample and wf as you found them.  The final electronic coordinates should never change.
 2) Try to keep the amount of data saved to the bare minimum you need, although if you know that in the future 
 you'll need some internal state to make nice output, then go ahead and store it in the logfile in the init section.
 
 If you find that you need more information than is provided via the evaluate() function(for example, maybe 
                                                                                         there are special weights for RMC or something), you can do the following:
 1) Subclass Average_generator with a different evaluate() function.  Try to keep it general so more than 
 one type can be included with the same parameters.
 2) Make an allocate() function for your subclass and add the subclass to the allocate() for Average_generator
 3) Keep an array of Subclass *'s in the averaging function and have them evaluate into the Properties_points, 
 appending to the normal list of Average_returns.
 That's about all.  Assuming that you write the write_init, read, and write_summaries, correctly, everything
 else should be taken care of.
 
 Futher hints:
 -Averages with weights can be done.  Pack the vals array with val1*weight1, weight1, val2*weight2, weight2, etc.
 The values and weights will be averaged properly, and in write_summary(), you can report <val1*weight1>/<weight1>
 as the average value.
 
 */
class Average_generator {
public:
  virtual ~Average_generator() {};
  //these are used in the main program.  evaluate() packs the
  //averaging data into the 'vals' array.  read() does as normal
  //write_init will write a section in the log file with relevant data
  //for gosling
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & )=0;

  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Properties_point & pt,Average_return &avg) {
    evaluate(wfdata,wf,sys,sample,avg);//in most of the case, just let it be the same; 
  }
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
			System * sys, Pseudopotential * psp, Sample_point * sample, Average_return &avg ) {
    evaluate(wfdata,wf,sys,sample,avg);//In most of the case, just let it be the same; 
  }
  
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Pseudopotential *psp, Sample_point * sample,
                        Properties_point & pt,Average_return &avg) {
    evaluate(wfdata,wf,sys, sample, pt, avg);
  }

  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words)=0;
  virtual void write_init(string & indent, ostream & os)=0;
  //This is used for Average_generators that use some kind of random sampling
  //in addition what's done already in the larger Monte Carlo stuff.  It's separated
  //from evaluate() so that correlated sampling can work efficiently
  virtual void randomize(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample) { } 
  
  //these are used in gosling: read will read in the write_init stuff
  //from above and set any relevant
  //internal variables.  gosling will average together all the
  // Average_returns and give a final averaged one
  //to write_summary, which will give a nice interpretation of the data.
  virtual void read(vector <string> & words)=0;
  virtual void write_summary(Average_return &,Average_return &, ostream & os)=0;
  virtual void jsonOutput(Average_return &,Average_return &, ostream & os) {
    error("jsonOutput not implemented for this Average_generator");
  }
    
};

int allocate(vector<string> & words, System * sys, Wavefunction_data * wfdata, Average_generator *& avg);
int allocate(vector<string> & words, Average_generator * & avg);

//##########################################################################

class Average_dipole:public Average_generator { 
public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & );
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
 private:
};


//----------------------------------------------------------------------------

class Average_structure_factor:public Average_generator { 
public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & );
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
  virtual void jsonOutput(Average_return &,Average_return &, ostream & os);
private:
  int npoints;
  Array2 <doublevar> kpts;
  
};

//----------------------------------------------------------------------------
class Average_magnetic_structure_factor:public Average_generator { 
public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & );
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
  virtual void jsonOutput(Average_return &,Average_return &, ostream & os);
private:
  int npoints;
  Array2 <doublevar> kpts;
  
};

//----------------------------------------------------------------------------

class Average_fourier_density:public Average_generator { 
public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & );
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
private:
  int npoints;
  Array2 <doublevar> kpts;
  
};
//----------------------------------------------------------------------------


class Average_twobody_correlation:public Average_generator { 
public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & );
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
  virtual void jsonOutput(Average_return &,Average_return &, ostream & os);
private:
  double resolution;
  int npoints;
  int direction; //we can look only in the x, y, or z direction if we like 0=r, 2=x,3=y,4=z.
};

//----------------------------------------------------------------------------

class Average_manybody_polarization:public Average_generator { 
public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & );
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
  virtual void jsonOutput(Average_return &,Average_return &, ostream & os);

private:
    Array2 <doublevar> gvec;
};  


//----------------------------------------------------------------------------

/*!
Spherically averaged one-body density matrix.
*/
class Average_obdm:public Average_generator { 
public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & );
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
private:
  int npoints;            //!< number of points in the plotting grid
  doublevar dR;           //!< grid spacing 
  int np_aver;            //!< number of points in quadrature for spherical average
  Array2 <doublevar> ptc; //!< cartesian coordinates of quadrature points
  Array1 <doublevar> wt;  //!< weights for quadrature points
};  

//----------------------------------------------------------------------------

/*!
Projected two-body density matrix, spherically averaged.
*/
class Average_tbdm:public Average_generator { 
public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & );
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
private:
  int nelectrons;
  int npoints;            //!< number of points in the plotting grid
  doublevar dR;           //!< grid spacing 
  int np_aver;            //!< number of points in quadrature for spherical average
  Array2 <doublevar> ptc; //!< cartesian coordinates of quadrature points
  Array1 <doublevar> wt;  //!< weights for quadrature points
};  

//----------------------------------------------------------------------------

/*!
Local magnetic (spin) moments and charges obtained by integrating electron
density in muffin-tin spheres.
*/
class Average_local_moments:public Average_generator { 
public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & );
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
private:
  int nelectrons;            //!< number of electrons
  int nup;                   //!< number of spin-up electrons
  int natoms;                //!< number of atoms
  Array1 <doublevar> rMT;    //!< radii of muffin-tin spheres
  vector <string> atomnames; //!< atom labels
};  

//----------------------------------------------------------------------------

//##########################################################################

class Average_density_moments:public Average_generator { 
public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & );
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
};


//----------------------------------------------------------------------------

//##########################################################################
class Average_linear_derivative:public Average_generator { 
public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & );
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
private:  
  doublevar tau;
  int unr;
  
};
//----------------------------------------------------------------------------

//##########################################################################
class Average_linear_delta_derivative:public Average_generator { 
public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & );
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
private:
  doublevar tau;
  int unr;
  
  
};
//----------------------------------------------------------------------------
/*!
Spherically averaged density on the basis
*/
class Average_spherical_density:public Average_generator { 
public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & );
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
private:
  int npoints;            //!< number of points in the plotting grid
  doublevar dR;           //!< grid spacing 
  doublevar cutoff;        //!< grid cutoff 
  int nup,ndown;         //!< number of electrons
  Basis_function * basis;
  int nfunc;
  vector <string> basistext;
  //int np_aver;            //!< number of points in quadrature for spherical average
  //Array2 <doublevar> ptc; //!< cartesian coordinates of quadrature points
  //Array1 <doublevar> wt;  //!< weights for quadrature points
};  

//----------------------------------------------------------------------------



//----------------------------------------------------------------------------

/*!
Spherically averaged density on the grid
*/
class Average_spherical_density_grid:public Average_generator { 
public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & );
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
private:
  int npoints;            //!< number of points in the plotting grid
  doublevar dR;           //!< grid spacing 
  int nup,ndown;         //!< number of elec
};  

//----------------------------------------------------------------------------

/*!
Density along a line.
  */
class Average_line_density:public Average_generator { 
public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & );
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
 private:
   Array1 <doublevar> vec;     //!< Line goes in this direction
   Array1 <doublevar> origin;  //!< Line starts here
   doublevar resolution;  //!< Bin size
   int npoints; //Number of points on the line
};

//----------------------------------------------------------------------------

/*!
 \brief 
 Accumulate the derivative of the wave function with respect to the parameters. 
 */
class Average_wf_parmderivs:public Average_generator { 
public:
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & avg);
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Properties_point & pt, Average_return & avg);
  virtual void evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
			System * sys, Pseudopotential * psp, Sample_point * sample, Properties_point & pt, Average_return &avg );
  
  virtual void read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words);
  virtual void write_init(string & indent, ostream & os);
  virtual void read(vector <string> & words);
  virtual void write_summary(Average_return &,Average_return &, ostream & os);
  virtual void jsonOutput(Average_return &,Average_return &, ostream & );
  
private:
  bool evaluate_pseudopotential;
  doublevar nodal_cutoff;
  
};





#endif //AVERAGE_GENERATOR_H_INCLUDED

