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

#include "MO_1d.h"
#include "qmc_io.h"

void MO_1d::buildLists(Array1 < Array1 <int> > & occupations) {
  moLists=occupations;
}

int MO_1d::showinfo(ostream & os) {
  os << "MO_1D" << endl;
  os << nmo << " MO's from " << readfile << endl;
  return 1;
}

int MO_1d::writeinput(string & indent, ostream & os) {
  os << endl;
  os << indent << "MO_1D" << endl;
  os << indent << "NMO " << nmo << endl;
  os << indent << "MESH_FILE " << readfile << endl;
  os << indent << "NPOINTS " << npoints << endl;
  return 1;

}

//----------------------------------------------------------------------

void MO_1d::read(vector <string> & words, unsigned int & startpos, 
                 System * sys) {

  if(!readvalue(words, startpos=0, nmo, "NMO")) 
    error("Need NMO in MO_1D");

  splines.Resize(nmo, 2);//one each for real and complex
  
  //Read the wave function file directly
  if(!readvalue(words, startpos=0, readfile, "MESH_FILE")) 
    error("Need MESH_FILE in MO_1D");

  int n;
  if(!readvalue(words, startpos=0, n, "NPOINTS"))
    error("Need NPOINTS IN MO_1D");

  doublevar magnify;
  if(!readvalue(words, startpos=0, magnify, "MAGNIFY"))
    magnify=1;

  npoints=n;

  Array1 <doublevar> x(n), y(n), yc(n);
  ifstream is(readfile.c_str());

  if(!is) error("Couldn't open ", readfile);

  int  mo=0;
  string dumline;
  while(mo < nmo && is) {
    
    //is.ignore(180, '\n');
    //is.ignore(180, '\n');
    getline(is, dumline);
    getline(is, dumline);
    doublevar dummy;

    doublevar max=0;
    for(int i=0; i< n; i++) {
      is >> dummy >> x(i) >> y(i) >> yc(i);
      if(max < fabs(y(i)))
        max=fabs(y(i));
    }

    max*=.1;
    for(int i=0; i< n; i++) {
      y(i)*=magnify/max;
      yc(i)*=magnify/max;
    }

    is.ignore(180,'\n');
    doublevar yp1=(y(1)-y(0))/(x(1)-x(0));
    doublevar ypn=(y(n-1)-y(n-2))/(x(n-1)-x(n-2));
    splines(mo,0).splinefit(x,y,yp1, ypn);
    
    yp1=(yc(1)-yc(0))/(x(1)-x(0));
    ypn=(yc(n-1)-yc(n-2))/(x(n-1)-x(n-2));
    splines(mo,1).splinefit(x,yc,yp1, ypn);
    mo++;
  }
  
  for(int i=0; i< nmo; i++) {
    for(int j=0; j< 2; j++) {
      if(!splines(0,0).match(splines(i,j))) {
        error("Spline ", i, " doesn't match 0 ");
      }
    }
  }

  //if(mpi_info.node==0 && 0) {
  //  for(int i=0; i< nmo; i++) {
  //    string spline="spline";
  //    append_number(spline, i);
  //    spline+=".dat";
  //    ofstream splnout(spline.c_str());
  //    for(doublevar r=x(0); r< x(n-1); r+=.005) {
//	int interval=splines(0,0).getInterval(r);
//	splnout << r << "   " << splines(i,0).getVal(r,interval) << "   "
//		<< splines(i,1).getVal(r,interval)
//		<< endl;
//      }
//      splnout.close();
//    }
//  }
  
}


//----------------------------------------------------------------------
void MO_1d::updateVal(Sample_point * sample,int e, int listnum,
                      Array2 <dcomplex> & newvals) {

  Array1 <doublevar> r(3);
  sample->getElectronPos(e,r);
  doublevar x=r(0);
  
  int interval=splines(0,0).getInterval(x);
  
  int nmo_list=moLists(listnum).GetDim(0);
  for(int m=0; m< nmo_list; m++) {
    int mo=moLists(listnum)(m);
    newvals(m,0)=dcomplex(splines(mo,0).getVal(x,interval), 
                           splines(mo,1).getVal(x,interval) );
  }

}

//----------------------------------------------------------------------
void MO_1d::updateLap(Sample_point * sample,int e,int listnum,
                      Array2 <dcomplex> & newvals) {

  Array1 <doublevar> r(3);
  sample->getElectronPos(e,r);
  doublevar x=r(0);
  
  int interval=splines(0,0).getInterval(x);
  
  int nmo_list=moLists(listnum).GetDim(0);
  doublevar f,fdir, f2dir, fi, fidir, fi2dir;
  for(int m=0; m< nmo_list; m++) {
    int mo=moLists(listnum)(m);
    splines(mo,0).getDers(x,interval, f,fdir, f2dir);
    splines(mo,1).getDers(x,interval, fi, fidir, fi2dir);
    newvals(m,0)=dcomplex(f,fi);
    newvals(m,1)=x*dcomplex(fdir, fidir);
    newvals(m,2)=newvals(m,3)=0;
    newvals(m,4)=dcomplex(f2dir, fi2dir);
  }
}

//----------------------------------------------------------------------
