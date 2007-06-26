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
//src/Mesh_function.cpp
#include "Mesh_function.h"
#include "qmc_io.h"
#include "Cubic_spline3D.h"
#include "jeep_utils.h"
/*!
\todo
Rewrite the code to do the edges.  As it is, it's quite inefficient
and confusing..

*/
int Mesh_function::read(
  vector <string> & words,
  unsigned int & pos)
{

  centername=words[0];

  rcut=1e99;

  vector <string> valfiles, dxfiles, dyfiles, dzfiles, lapfiles;
  readsection(words, pos=0, valfiles, "VALPLT");
  readsection(words, pos=0, dxfiles,  "DXPLT");
  readsection(words, pos=0, dyfiles,  "DYPLT");
  readsection(words, pos=0, dzfiles,  "DZPLT");
  readsection(words, pos=0, lapfiles, "LAPPLT");

  nmax=valfiles.size();
  if((valfiles.size()!=dxfiles.size())
      || (valfiles.size() != dyfiles.size())
      || (valfiles.size() !=dzfiles.size())
      || (valfiles.size() !=lapfiles.size()))

  {
    cout << "num of vals " << valfiles.size() << endl;
    cout << "num of dx " << dxfiles.size() << endl;
    cout << "num of dy " << dyfiles.size() << endl;
    cout << "num of dz " << dzfiles.size() << endl;
    cout << "num of lap " << lapfiles.size() << endl;
    error("Must have the same number of files for the values and derivatives");
  }

  //int count,nfunctions,natoms;
  Array1 <doublevar> origin(3),box_size(3);
  Array1 <int> indie_points(3),grid_points(3);
  string orbfile_0,orbfile_1,orbfile_2,orbfile_3,orbfile_lap;
  //doublevar dummy1,dummy2,dummy3,dummy4,dummy5;
  single_write(cout, "Mesh Basis\n");

  //iterations for all orbitals
  grid.Resize(nmax,5);
  minmax.Resize(nmax);
  resolution.Resize(nmax);

  for (int m=0;m< nmax;m++){
    orbfile_0=valfiles[m];
    orbfile_1=dxfiles[m];
    orbfile_2=dyfiles[m];
    orbfile_3=dzfiles[m];
    orbfile_lap=lapfiles[m];


    ifstream ORB_0(orbfile_0.c_str());
    if(!ORB_0)
	{
	  error("couldn't find orb file ", orbfile_0);
	}
    ifstream ORB_1(orbfile_1.c_str());
    if(!ORB_1)
	{
	  error("couldn't find orb file ", orbfile_1);
	}
    ifstream ORB_2(orbfile_2.c_str());
    if(!ORB_2)
	{
	  error("couldn't find orb file ", orbfile_2);
	}
    ifstream ORB_3(orbfile_3.c_str());
    if(!ORB_3)
	{
	  error("couldn't find orb file ", orbfile_3);
	}
    ifstream ORB_lap(orbfile_lap.c_str());
    if(!ORB_lap)
	{
	  error("couldn't find orb file ", orbfile_lap);
	}


    Array1 <int> padding(3);

    padding=7;
    Array1 <int> startfill(3);
    startfill=3;

    get_grid_from_plot(ORB_0, grid(m,0), box_size, origin, padding, startfill);
    get_grid_from_plot(ORB_1, grid(m,1), box_size, origin, padding, startfill);
    get_grid_from_plot(ORB_2, grid(m,2), box_size, origin, padding, startfill);
    get_grid_from_plot(ORB_3, grid(m,3), box_size, origin, padding, startfill);
    get_grid_from_plot(ORB_lap, grid(m,4), box_size, origin, padding, startfill);


    ORB_0.close();
    ORB_1.close();
    ORB_2.close();
    ORB_3.close();
    ORB_lap.close();


    Array1 <int> indie_points(3);

    for(int i=0; i< 3; i++) indie_points(i)=grid(m,0).GetDim(i)-padding(i);

    minmax(m).Resize(6);
    resolution(m).Resize(3);

    cout << "resolution of grid ";
    for (int i=0;i<3;i++) {
      resolution(m)(i)=box_size(i)/(indie_points(i)-1);
      cout << resolution(m)(i) << "    ";
      minmax(m)(2*i)=origin(i)-3.0*resolution(m)(i);
      minmax(m)(2*i+1)=origin(i)+box_size(i)+4.0*resolution(m)(i);
    }
    cout << endl;

    cout << "origin ";
    for(int i=0; i< 3; i++) cout << origin(i) << "   " ;
    cout << endl;

    cout << "box_size ";
    for(int i=0; i< 3; i++) cout << box_size(i) << "   ";
    cout << endl;

    cout << "npoints ";
    for(int i=0; i< 3; i++) cout << indie_points(i) << "    ";
    cout << endl;

    for(int i=0; i< 3; i++) grid_points(i)=indie_points(i)+padding(i);

    for(int i=0;i<grid_points(0);i++)
      for(int j=0;j<grid_points(1);j++)
        for(int k=0;k<grid_points(2);k++){
         if (i<3){
           if (j<3){
             if(k<3) for (int l=0;l<5;l++){
               grid(m,l)(i,j,k)=grid(m,l)(i+indie_points(0),
                                          j+indie_points(1),
                                          k+indie_points(2));
		         }
             else  for (int l=0;l<5;l++){
		           grid(m,l)(i,j,k)=grid(m,l)(i+indie_points(0),j+indie_points(1),k);
            }
          }
	        else {
            if (k<3) for (int l=0;l<5;l++){
              grid(m,l)(i,j,k)=grid(m,l)(i+indie_points(0),j,k+indie_points(2));
            }
            else for (int l=0;l<5;l++){
             grid(m,l)(i,j,k)=grid(m,l)(i+indie_points(0),j,k);
           }
	       }
	    }
	    else {
	      if (j<3){
          if(k<3) for (int l=0;l<5;l++){
            grid(m,l)(i,j,k)=grid(m,l)(i,j+indie_points(1),k+indie_points(2));
          }
          else for (int l=0;l<5;l++){
           grid(m,l)(i,j,k)=grid(m,l)(i,j+indie_points(1),k);
          }
	      }
	      else if (k<3) for (int l=0;l<5;l++){
          grid(m,l)(i,j,k)=grid(m,l)(i,j,k+indie_points(2));
	      }
	    }
	    if (i>indie_points(0)+2){
	      if (j>indie_points(1)+2){
          if(k>indie_points(2)+2) for (int l=0;l<5;l++){
            grid(m,l)(i,j,k)=grid(m,l)(i-indie_points(0),
                                       j-indie_points(1),
                                       k-indie_points(2));
		      }
		      else for (int l=0;l<5;l++){
            grid(m,l)(i,j,k)=grid(m,l)(i-indie_points(0),j-indie_points(1),k);
          }
	      }
	      else {
          if (k>indie_points(2)+2)for (int l=0;l<5;l++){
            grid(m,l)(i,j,k)=grid(m,l)(i-indie_points(0),j,k-indie_points(2));
          }
          else for (int l=0;l<5;l++){
            grid(m,l)(i,j,k)=grid(m,l)(i-indie_points(0),j,k);
          }
	      }
	    }
	    else {
	      if (j>indie_points(1)+2){
          if(k>indie_points(2)+2) for (int l=0;l<5;l++){
            grid(m,l)(i,j,k)=grid(m,l)(i,j-indie_points(1),k-indie_points(2));
          }
          else for (int l=0;l<5;l++){
            grid(m,l)(i,j,k)=grid(m,l)(i,j-indie_points(1),k);
          }
	      }
	      else if (k>indie_points(2)+2) for (int l=0;l<5;l++){
          grid(m,l)(i,j,k)=grid(m,l)(i,j,k-indie_points(2));
	      }
      }
	  }

    single_write(cout,"Reading of",m+1,"-th orbital done ...");
    cout <<endl;
  }


  return 0;
}

//------------------------------------------------------------------------

void Mesh_function::getVarParms(Array1 <doublevar> & parms) {
  //cout << "getVarParms " << endl;
  parms.Resize(0);
}
//------------------------------------------------------------------------
void Mesh_function::setVarParms(Array1 <doublevar> & parms) {
  //cout << "setVArParms " << endl;
  assert(parms.GetDim(0)==0);

}

//------------------------------------------------------------------------
int Mesh_function::nfunc()
{
  return nmax;
}

//------------------------------------------------------------------------
int Mesh_function::showinfo(string & indent, ostream & os)
{
  os << indent << "Mesh function\n";
  os << indent << "Number of functions " << nmax << endl;
  return 1;
}

//------------------------------------------------------------------------
int Mesh_function::writeinput(string & indent, ostream & os)
{
  os << indent << centername << endl;
  os << indent << "MESH\n";
  os << indent << "CUTOFF " << rcut << endl;
  return 1;
}

//------------------------------------------------------------------------
void Mesh_function::raw_input(ifstream & input)
{error("Raw input not supported by Mesh_function");}


//------------------------------------------------------------------------

void Mesh_function::calcVal(const Array1 <doublevar> & r,
                                Array1 <doublevar> & symvals,
                                const int startfill)
{
  static Array1 <doublevar> xyz(3);
  static Array1 <doublevar> rmin(3);
  for(int d=0; d< 3; d++) {
    xyz(d)=r(d+2);
  }

  int index=startfill;
  for(int m=0; m < nmax; m++) {
    rmin(0)=minmax(m)(0);
    rmin(1)=minmax(m)(2);
    rmin(2)=minmax(m)(4);
    splin3_my(rmin, resolution(m), grid(m,0), xyz, symvals(index));
    index++;
  }
}

//------------------------------------------------------------------------

void Mesh_function::calcLap(
  const Array1 <doublevar> & r,
  Array2 <doublevar> & symvals,
  const int startfill
)
{
  static Array1 <doublevar> xyz(3);
  static Array1 <doublevar> rmin(3);
  for(int d=0; d< 3; d++) {
    xyz(d)=r(d+2);
  }

  int index=startfill;
  for(int m=0; m < nmax; m++) {
    rmin(0)=minmax(m)(0);
    rmin(1)=minmax(m)(2);
    rmin(2)=minmax(m)(4);
    //cout <<m<<endl;
    splin3_my(rmin, resolution(m), grid(m,0), xyz, symvals(index,0));
    splin3_my(rmin, resolution(m), grid(m,1), xyz, symvals(index,1));
    splin3_my(rmin, resolution(m), grid(m,2), xyz, symvals(index,2));
    splin3_my(rmin, resolution(m), grid(m,3), xyz, symvals(index,3));
    splin3_my(rmin, resolution(m), grid(m,4), xyz, symvals(index,4));
    //cout << "function value" <<newvals(m,0)<<endl;
    index++;
  }



}

//------------------------------------------------------------------------
