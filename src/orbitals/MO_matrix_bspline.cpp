/*
 
Copyright (C) 2007 Jindrich Kolorenc, Michal Bajdich

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



#include "Qmc_std.h"
#include "MO_matrix_bspline.h"
#include "Sample_point.h"
#include "qmc_io.h"

//--------------------------------------------------------------------------

void MO_matrix_bspline::init() {

  single_write(cout, "Bspline MO\n");
  Einspline.Resize(nmo);
  for (int m=0;m< nmo;m++){
    string orbfile_0=valfiles[m];
    
    ifstream ORB_0(orbfile_0.c_str());
    if(!ORB_0){
      error("couldn't find orb file ", orbfile_0);
    }
    Einspline[m]=new EinsplineOrb;

    if(periodic)
      Einspline[m]->periodic=1;
    else
      Einspline[m]->periodic=0;

    Einspline[m]->read(ORB_0);
    ORB_0.close();
  }

}


//----------------------------------------------------------------------


void MO_matrix_bspline::buildLists(Array1 < Array1 <int> > & occupations) {
  int numlists=occupations.GetDim(0);
  moLists.Resize(numlists);
  for(int lis=0; lis < numlists; lis++) {
    int nmo_list=occupations(lis).GetDim(0);
    moLists(lis).Resize(nmo_list);
    for(int mo=0; mo < nmo_list; mo++) {
      moLists(lis)(mo)=occupations(lis)(mo);
    }
  }
}

//----------------------------------------------------------------------

int MO_matrix_bspline::showinfo(ostream & os)
{
  os << "Bspline Molecular Orbital\n";
  string indent="  ";
  os << "Using orbital value plotfiles : \n";
  for(int i=0;i<valfiles.size();i++){
    os << indent << valfiles[i] <<endl; 
    os << indent <<"origin:   "<<Einspline[i]->origin(0)<<"  "<< Einspline[i]->origin(1)<<"  "<<Einspline[i]->origin(2)<<endl;
    os << indent <<"box_size: "<<Einspline[i]->box_size(0)<<"  "<<Einspline[i]->box_size(1)<<"  "<<Einspline[i]->box_size(2)<<endl;
    os << indent <<"spacing:  "<<Einspline[i]->spacing(0)<<"  "<<Einspline[i]->spacing(1)<<"  "<<Einspline[i]->spacing(2)<<endl;
  }
  os << "Number of molecular orbitals: " << nmo << endl;
  
  return 1;
}

int MO_matrix_bspline::writeinput(string & indent, ostream & os)
{
  os << indent << "BSPLINE_MO" << endl;
  os << indent << "NMO " << nmo << endl;
  os << indent << "MAGNIFY " << magnification_factor << endl;
  string indent2=indent+"  ";
  os << indent2 << "VALPLT {" <<endl;

  for(int i=0;i<valfiles.size();i++)
    os << indent2 << valfiles[i] <<endl; 

  os << indent2 << "}"<<endl;
  os << indent << "}" << endl;
  return 1;
}

void MO_matrix_bspline::read(vector <string> & words, unsigned int & startpos, System * sys){
   unsigned int pos=startpos;

   if(!readvalue(words, pos, nmo, "NMO"))
     {
       error("Need NMO in molecular orbital section");
     }
   
   if(nmo > 40000) 
     error("You have entered more than 40,000 for NMO.  This seems a bit big; we most likely"
	   " can't handle it.");
   
   
   pos=0;
   if(!readvalue(words, pos, magnification_factor, "MAGNIFY")) {
     magnification_factor=1;
   }
   pos=0;

   pos=0;
   if(readsection(words, pos, valfiles, "VALPLT")) {
     if(valfiles.size()!=nmo){
       cout << "#plotfiles " << valfiles.size()<<" inside VALPLT" << endl;
       error("Must have the same number as NMO");
     }
   }
   else{
     error("Need VALPLT section in BSPLINE_MO ORBITALS");
   }
   if(haskeyword(words, pos=0, "PERIODIC") )
     periodic=1;
   else periodic=0; 
   
   init();
}

//------------------------------------------------------------------------

void MO_matrix_bspline::updateVal(Sample_point * sample, int e,
                                   int listnum,
                                   //!< which list to use
                                   Array2 <doublevar> & newvals
                                   //!< The return: in form (MO, val)
) {
  newvals=0;
  int nmo_list=moLists(listnum).GetDim(0);
  Array1 <doublevar> xyz(3);
  Array1 <doublevar> vals(1);
  sample->getElectronPos(e,xyz);

  for(int m=0; m < nmo_list; m++) {
     int mo=moLists(listnum)(m);
     Einspline[mo]->evaluate_spline (xyz, vals);
     for(int d=0;d<vals.GetSize();d++)
       newvals(mo,d)= magnification_factor*vals(d);
  }
}

//------------------------------------------------------------------------


void MO_matrix_bspline::updateLap(Sample_point * sample, int e,
				  int listnum,
				  //!< which list to use
				  Array2 <doublevar> & newvals
				  //!< The return: in form (MO, [val, grad, lap])
) {
  newvals=0;
  Array1 <doublevar> vals(5);
  Array1 <doublevar> xyz(3);
  int nmo_list=moLists(listnum).GetDim(0);
  sample->getElectronPos(e,xyz);

  for(int m=0; m < nmo_list; m++) {
     int mo=moLists(listnum)(m);
     Einspline[mo]->evaluate_spline_vgl (xyz, vals);
     for(int d=0;d<vals.GetSize();d++)
       newvals(mo,d)= magnification_factor*vals(d);
  }
}

//--------------------------------------------------------------------------


void MO_matrix_bspline::updateHessian(
  Sample_point * sample,
  int e,
  int listnum,
  //const Array1 <int> & occupation,
  //!<A list of the MO's to evaluate
  Array2 <doublevar> & newvals
  //!< The return: in form (MO, [val, grad, dxx,dyy,...])
)
{
  newvals=0;
  Array1 <doublevar> vals(10);
  Array1 <doublevar> xyz(3);
  int nmo_list=moLists(listnum).GetDim(0);
  sample->getElectronPos(e,xyz);

  for(int m=0; m < nmo_list; m++) {
     int mo=moLists(listnum)(m);
     Einspline[mo]->evaluate_spline_vgh (xyz, vals);
     for(int d=0;d<vals.GetSize();d++)
       newvals(mo,d)= magnification_factor*vals(d);
  }
}
//--------------------------------------------------------------------------
