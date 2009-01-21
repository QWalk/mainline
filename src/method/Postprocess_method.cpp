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

#include "Program_options.h"
#include "Postprocess_method.h"
#include "qmc_io.h"
#include "System.h"
#include "Sample_point.h"
#include "ulec.h"

/*!
Read the "words" from the method section in the input file
via doinput() parsing, and store section information in private
variables orbs, resolution, and minmax.
Set up MO_matrix and Sample_point objects for wavefunction from input
*/
void Postprocess_method::read(vector <string> words,
                       unsigned int & pos,
                       Program_options & options)
{
  sys=NULL;
  allocate(options.systemtext[0],  sys);
  sample=NULL;
  sys->generateSample(sample);

  if(!readvalue(words, pos=0, configfile, "READCONFIG"))
    error("Need READCONFIG in Postprocess");

  canonical_filename(configfile);
}

class Postprocess_accumulator {
public:
  Postprocess_accumulator(System * sys) { 
    npoints=0;
    nelectrons=sys->nelectrons(0)+sys->nelectrons(1);
    
    vector <double> cvg_tmp;
    for(double i=1; i <= nelectrons; i++) 
      cvg_tmp.push_back(i);
    

    cvg.Resize(cvg_tmp.size());
    for(int i=0; i< cvg.GetDim(0); i++) cvg(i)=cvg_tmp[i];
   

    sys->getAtomicLabels(atom_names);
    
    zpol_cvg_avg.Resize(cvg.GetDim(0),3);
    zpol_cvg_var.Resize(cvg.GetDim(0),3);
    zpol_cvg_avg=0;
    zpol_cvg_var=0;
    
    zpol_single.Resize(3);
    zpol_single_var.Resize(3);
    zpol_single=dcomplex(0.0,0.0);zpol_single_var=dcomplex(0.0,0.0);
  }


  //--------------------------------------------------
  void accumulate(System * sys, Sample_point * sample) {


    int nelectrons=sample->electronSize();
    Array2 <doublevar> gvec;
    //Array2 <doublevar> latvec;
    Array1 <doublevar> pos(3);
    sys->getRecipLattice(gvec);
    //sys->getBounds(latvec);

    dcomplex zpol;
    Array1 <doublevar> sum(3,0.0);
    for(int e=0; e< nelectrons; e++) {
      sample->getElectronPos(e,pos);
      for(int i=0; i< 3; i++) {
        for(int d=0; d< 3; d++) 
          sum(i)+=gvec(i,d)*pos(d);
      }
    }
    

    for(int j=0; j< cvg.GetDim(0); j++) { 
      for(int i=0; i< 3; i++) {
        double temp_sum=sum(i);
        zpol=dcomplex(cos(2*pi*temp_sum/doublevar(cvg(j))),
                      sin(2*pi*temp_sum/doublevar(cvg(j))));

        dcomplex old_av=zpol_cvg_avg(j,i);
        dcomplex new_av=old_av+(zpol-old_av)/doublevar(npoints+1);
        dcomplex old_var=zpol_cvg_var(j,i);
        dcomplex new_var=old_var;
        if(npoints>1)
          new_var=dcomplex(doublevar(1-1.0/npoints)*old_var.real()
                           +doublevar(npoints+1)*(new_av.real()-old_av.real())
                           *(new_av.real()-old_av.real()),
                           doublevar(1-1.0/npoints)*old_var.imag()
                           +doublevar(npoints+1)*(new_av.imag()-old_av.imag())
                           *(new_av.imag()-old_av.imag()));
        zpol_cvg_avg(j,i)=new_av;
        zpol_cvg_var(j,i)=new_var;

      }
    }
    

    for(int i=0; i< 3; i++) {
      zpol=dcomplex(0.0,0.0);
      for(int e=0; e< nelectrons; e++) {
        doublevar ov=0;
        for(int d=0; d< 3; d++) { 
          ov+=gvec(i,d)*pos(d);
        }
        zpol+=dcomplex(cos(2*pi*ov),sin(2*pi*ov));
        dcomplex old_av=zpol_single(i);
        dcomplex new_av=old_av+(zpol-old_av)/doublevar(npoints+1);
        dcomplex old_var=zpol_single_var(i);
        dcomplex new_var=old_var;
        if(npoints>1)
          new_var=dcomplex(doublevar(1-1.0/npoints)*old_var.real()
                           +doublevar(npoints+1)*(new_av.real()-old_av.real())
                           *(new_av.real()-old_av.real()),
                           doublevar(1-1.0/npoints)*old_var.imag()
                           +doublevar(npoints+1)*(new_av.imag()-old_av.imag())
                           *(new_av.imag()-old_av.imag()));
        zpol_single(i)=new_av;
        zpol_single_var(i)=new_var;        
      }
    }
    npoints++;
  }

  //----------------------------------------
  void printout(ostream & os) { 

    os << "Single-point zpol " << endl;
    for(int d=0; d< 3; d++) {
      os  << "sing " <<  d<< "   " <<zpol_single(d).real()<< "  " 
      << zpol_single(d).imag() 
      << "  " << sqrt(zpol_single_var(d).real()/npoints)
      << "  " << sqrt(zpol_single_var(d).imag()/npoints) <<  "   " 
      <<  endl;      
    } 
    
    os << "Convergence: " << endl;
    os << "d  n    z +/- zerr phase(without error bars) \n";
    
    for(int j=0; j< cvg.GetDim(0); j++) {
      for(int d=0; d< 3; d++) {
        
        os <<cvg(j) << "   " <<  d<< "   " <<zpol_cvg_avg(j,d).real()<< "  " 
        << zpol_cvg_avg(j,d).imag() 
        << "  " << sqrt(zpol_cvg_var(j,d).real()/npoints)
        << "  " << sqrt(zpol_cvg_var(j,d).imag()/npoints) <<  "   " 
        <<  endl;
      }
    }

  }
  //----------------------------------------------------

private:
  long int npoints;
  int nelectrons;
  Array2 <dcomplex>  zpol_cvg_avg; //!< n-body zpol
  Array2 <dcomplex> zpol_cvg_var;
  Array1 <doublevar> cvg;

  Array1 <doublevar> atomcharges;
  vector <string> atom_names;
  
  
  Array1 <dcomplex> zpol_single; //!< single-body zpol
  Array1 <dcomplex> zpol_single_var;
};

//----------------------------------------------------------------------

void Postprocess_method::run(Program_options & options, ostream & output) {

  ifstream checkfile(configfile.c_str());
  if(!checkfile) error("Couldn't open ", configfile);
  string dummy;
  int nread=0;
  Postprocess_accumulator accumulator(sys);
  while(checkfile >> dummy) {
    if(read_config(dummy, checkfile,sample)) {
      accumulator.accumulate(sys,sample);
      nread++;
      if(nread%1000==0) output << nread << " configs read" << endl;
    }
  }
  checkfile.close();
  accumulator.printout(output);
  

}


/*!

*/
int Postprocess_method::showinfo(ostream & os)
{
  os<<"#############Postprocess_method#################\n";
  sys->showinfo(os);

  return 1;
}
