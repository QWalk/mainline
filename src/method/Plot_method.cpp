/*
 
Copyright (C) 2007 Zachary Helms
 with further modifications by Lucas K. Wagner

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

#include "Plot_method.h"
#include "qmc_io.h"
#include "System.h"
#include "Program_options.h"
/*!
Read the "words" from the method section in the input file
via doinput() parsing, and store section information in private
variables orbs, resolution, and minmax.
Set up MO_matrix and Sample_point objects for wavefunction from input
*/
void Plot_method::read(vector <string> words,
                       unsigned int & pos,
                       Program_options & options)
{
  pos=0;	//always start from first word
  vector <string> Torbs;
  vector <string> Tminmax;
  vector <string> orbtext;


  if(!readsection(words, pos=0, orbtext, "ORBITALS"))
    error("Need ORBITALS in PLOT section");

  if(! readvalue(words,pos=0,resolution,"RESOLUTION"))
    resolution=.2;

  if(! readsection(words,pos=0,Tminmax,"MINMAX"))
    error("Need MINMAX in METHOD section");
  if(Tminmax.size() != 6)
    error("MINMAX needs 6 values");
  minmax.Resize(6);
  for(unsigned int i=0; i<Tminmax.size(); i++) {
    minmax(i)=atof(Tminmax[i].c_str());
  }


  sysprop=NULL;
  allocate(options.systemtext[0],  sysprop);
  allocate(orbtext, sysprop, mymomat);
  int maxorb=0;
  if(readsection(words,pos=0,Torbs,"PLOTORBITALS")) {
    orbs.Resize(Torbs.size());
    for(unsigned int i=0; i<Torbs.size(); i++) {
      orbs(i)=atoi(Torbs[i].c_str());
      if(orbs(i)> maxorb) maxorb=orbs(i);
    }
  }
  else {
    int totmo=mymomat->getNmo();
    orbs.Resize(totmo);
    for(int i=0; i< totmo; i++) {
      orbs(i)=i+1;
    }
  }
  if(maxorb > mymomat->getNmo()) {
    error("Too high orbital requested in PLOTORBITALS.  Try increasing " 
	  "NORB in the ORBITALS section.");
  }

  Array1 <Array1 <int> > orblist(1);
  orblist(0).Resize(orbs.GetDim(0));
  for(int i=0; i< orbs.GetDim(0); i++) {
    //cout << "i " << i << endl;
    //cout << "orbs " << orbs(i) << endl;
    orblist(0)(i)=orbs(i)-1;
  }

  mymomat->buildLists(orblist);

  mywalker=NULL;
  sysprop->generateSample(mywalker);
  mymovals.Resize(mymomat->getNmo(),5);
}

/*!
 
*/
void Plot_method::run(Program_options & options, ostream & output) {
  ofstream os; //for writing to *.plt files
  unsigned int electron=0; //# of the electron that will move through grid
  string pltfile; //name of plotfile being written
  Array1 <doublevar> xyz(3),resolution_array(3); //position of electron "in" MO
  Array1 <int> D_array1(3); //dummy array1
  D_array1=0; //sets all 3 components to 0. use as counter for gridpoints

  D_array1(0)=roundoff((minmax(1)-minmax(0))/resolution);
  D_array1(1)=roundoff((minmax(3)-minmax(2))/resolution);
  D_array1(2)=roundoff((minmax(5)-minmax(4))/resolution);

  resolution_array(0)=(minmax(1)-minmax(0))/(D_array1(0)-1);  
  resolution_array(1)=(minmax(3)-minmax(2))/(D_array1(1)-1);
  resolution_array(2)=(minmax(5)-minmax(4))/(D_array1(2)-1);

  int npts=D_array1(0)*D_array1(1)*D_array1(2);
  Array2 <doublevar> grid(orbs.GetSize(),npts);
  Array1 <doublevar> density(npts);
  //generate .xyz file for gOpenMol to view coordinates
  pltfile=options.runid + ".xyz";
  os.open(pltfile.c_str());
  write_xyz(sysprop,os);
  os.close();


  //calculate value of each molecular orbital at each grid point and store in an Array1
  // grid values with x=fastest running variable, and z=slowest
  cout<<"calculating "<<D_array1(0)*D_array1(1)*D_array1(2) <<" grid points"<<endl;
  cout<<"for "<< orbs.GetDim(0) <<" molecular orbitals"<<endl;
  int count=0;
  xyz(0)=minmax(0);
  xyz(1)=minmax(2);
  xyz(2)=minmax(4); //init elec probe to xmin ymin zmin
   for(int xx=0;xx<D_array1(0);xx++){
     xyz(0)=minmax(0)+xx*resolution_array(0); //move forward on z axis one resolution unit
     cout << "x " << xyz(0) << endl;
     for(int yy=0; yy<D_array1(1);yy++){
       xyz(1)=minmax(2)+yy*resolution_array(1);  
       for(int zz=0; zz<D_array1(2);zz++){
         xyz(2)=minmax(4)+zz*resolution_array(2);
         mywalker->setElectronPos(electron,xyz); //move elec#1 to point specified by xyz
         mymomat->updateVal(mywalker,electron,0,mymovals); //recalculate MO value for elec#1
         density(count)=0;
         for(int i=0; i<orbs.GetSize(); i++) {
             grid(i,count)=mymovals(i,0);
             density(count)+=mymovals(i,0)*mymovals(i,0);
         }
         count++; //index for cycling through grid points
       }
     }
   }

  //Loop through and generate plot files for each orbital requested
  if(orbs.GetSize()<=0)
    error("number of orbitals requested is not a positive number");
  cout<<"saving data for "<<orbs.GetSize()<<" molecular orbitals"<<endl;
  
  
  for(int i=0; i<orbs.GetSize(); i++) {
    //output to file with orbital number in it
    string basename;
    char strbuff[40];
    sprintf(strbuff, "%d", orbs(i));
    basename = options.runid;
    basename += ".orb";
    basename += strbuff;
    
    int dopltfile=0;
    int cubefile=1;
    //Note that the pltfile will be rotated, since it requires z,y,x, 
    //but cube requires x,y,z.
    if(dopltfile) {
      pltfile = basename+".plt"; 
      os.open(pltfile.c_str());
      cout<<"writing to "<<pltfile<<endl;

      // http://www.csc.fi/gopenmol/developers/plt_format.phtml
      os<<"3 "; //rank=3 always
      os<<"2\n"; //dummy variable => "Orbital/density surface"
      //number of grid points for x, y, & z direction
      os <<D_array1(2)<<" "<<D_array1(1)<<" "<<D_array1(0)<<endl;
      os <<minmax(4)<<" "<<minmax(5)<<" "<<minmax(2)<<" "<<minmax(3)
         <<" "<<minmax(0)<<" "<<minmax(1)<<endl;

      for(int j=0; j<(D_array1(0)*D_array1(1)*D_array1(2)); j++) {
        os<<grid(i,j)<<endl;
      }
      os.close();
    }
    if(cubefile) {
      string cubename=basename+".cube";
      os.open(cubename.c_str());
      os << "GOS plot output\n";
      os << "Molecular orbital " << orbs(i) << endl;
      int natoms=sysprop->nIons();
      os << "  " << natoms << "   " << minmax(0) << "   "
          << minmax(2) << "   " << minmax(4) << endl;
      os << D_array1(0) << "   " << resolution_array(0) << "  0.0000   0.0000" << endl;
      os << D_array1(1) << "   0.0000   " << resolution_array(1) << "  0.0000" << endl;
      os << D_array1(2) << "   0.0000    0.0000    " << resolution_array(2) << endl;
      Array1 <doublevar> pos(3);
      for(int at=0; at< natoms; at++) {
        mywalker->getIonPos(at,pos);
        os << "   " << mywalker->getIonCharge(at) << "   0.0000    " << pos(0) 
            <<"    " << pos(1) << "   " << pos(2) << endl;
      }
      
      for(int j=0; j< npts; j++) {
        os << grid(i,j) <<  "    ";
        if(j%6 ==5) os << endl;
      }
      os << endl;
      os.close();
      
    }
  }

  
  string cubename=options.runid+".dens.cube";
  os.open(cubename.c_str());
  os << "GOS plot output\n";
  os << "Electron density" << endl;
  int natoms=sysprop->nIons();
  os << "  " << natoms << "   " << minmax(0) << "   "
      << minmax(2) << "   " << minmax(4) << endl;
  os << D_array1(0) << "   " << resolution_array(0) << "  0.0000   0.0000" << endl;
  os << D_array1(1) << "   0.0000   " << resolution_array(1) << "  0.0000" << endl;
  os << D_array1(2) << "   0.0000    0.0000    " << resolution_array(2) << endl;
  Array1 <doublevar> pos(3);
  for(int at=0; at< natoms; at++) {
    mywalker->getIonPos(at,pos);
    os << "   " << mywalker->getIonCharge(at) << "   0.0000    " << pos(0) 
        <<"    " << pos(1) << "   " << pos(2) << endl;
  }
      
  for(int j=0; j< npts; j++) {
    os << density(j) <<  "    ";
    if(j%6 ==5) os << endl;
  }
  os << endl;
  os.close(); 

}


/*!
Print information about private variables {orbs,resolution,minmax}
*/
int Plot_method::showinfo(ostream & os)
{
  os<<"#############Plot_method#################\n";
  os<<"orbs="<<orbs(0);
  for(int i=1;i<orbs.GetSize();i++)
    os<<", "<<orbs(i);
  os<<endl;
  os<<"resolution="<<resolution<<endl;
  os<<"xmin="<<minmax(0)<<" xmax="<<minmax(1)<<endl;
  os<<"ymin="<<minmax(2)<<" ymax="<<minmax(3)<<endl;
  os<<"zmin="<<minmax(4)<<" zmax="<<minmax(5)<<endl;
  os<<"done"<<endl;
  return 1;
}
