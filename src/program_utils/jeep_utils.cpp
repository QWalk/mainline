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
//----------------------------------------------------------------------

#include "jeep_utils.h"
#include "qmc_io.h"

//----------------------------------------------------------------------

void get_grid_from_plot(istream & plotfile,
                        Array3 <doublevar> & grid,
                        Array1 <doublevar> & box_size,
                        Array1 <doublevar> & origin ) {
  Array1 <int> padding(3), startfill(3);
  padding=0;startfill=0;
  get_grid_from_plot(plotfile, grid, box_size, origin, padding, startfill);
}


void get_grid_from_plot(istream & plotfile,
                        Array3 <doublevar> & grid,
                        Array1 <doublevar> & box_size,
                        Array1 <doublevar> & origin,
                        Array1 <int> & padding,
                        Array1 <int> & startfill) {

  Array1 <int> npoints(3);
  box_size.Resize(3);
  origin.Resize(3);


  int nfunctions;
  get_global_header(plotfile, nfunctions);
  Array1 <int> start(3);
  Array1 <int> end(3);

  for(int function=0; function < nfunctions; function++) {

    get_function_header(plotfile, npoints, box_size, origin);

    if(function==0) {
      grid.Resize(npoints(0)+padding(0), npoints(1)+padding(1), npoints(2)+padding(2));
      grid=0;
    }

    for(int i=0; i< 3; i++) {
      start(i)=startfill(i);
      end(i)=npoints(i)+startfill(i);
      assert(end(i) <= npoints(i)+padding(i));
    }

    for(int x=start(0); x < end(0); x++) {
    for(int y=start(1); y < end(1); y++) {
    for(int z=start(2); z < end(2); z++) {
      doublevar mo_val;

      if(!(plotfile >> mo_val)) {
        error("error reading plot file");
      }


      if(function==0) {
      }
      else if(function==1) {
        mo_val=-mo_val;
      }
      else {
        error("too many functions..");
      }
      grid(x,y,z)+= mo_val;


    }
    }
    }

  }
}

//----------------------------------------------------------------------

void get_function_header(istream & plotfile,
                Array1 <int> & npoints,
                Array1 <doublevar> & box_size,
                Array1 <doublevar> & origin) {


  plotfile >> origin(0) >> origin(1) >> origin(2);
  plotfile.ignore(180, '\n');
  plotfile >> box_size(0) >> box_size(1) >> box_size(2);
  plotfile.ignore(180, '\n');
  plotfile.ignore(180, '\n');
  plotfile.ignore(180, '\n');
  plotfile >> npoints(0) >> npoints(1) >> npoints(2);
  plotfile.ignore(180, '\n');
  /*
  cout << "origin ";
  for(int i=0; i< 3; i++) cout << origin(i) << "   " ;
  cout << endl;

  cout << "box_size ";
  for(int i=0; i< 3; i++) cout << box_size(i) << "   ";
  cout << endl;

  cout << "npoints ";
  for(int i=0; i< 3; i++) cout << npoints(i) << "    ";
  cout << endl;
  */


}

//----------------------------------------------------------------------

void get_global_header(istream & plotfile, int & nfunctions) {

  for(int i=0; i< 8; i++) plotfile.ignore(180, '\n');
  int natoms;
  plotfile >> natoms; plotfile.ignore(180, '\n');
  //cout << "natoms " << natoms << endl;
  for(int i=0; i< natoms; i++) plotfile.ignore(180, '\n');

  plotfile >> nfunctions;
  //cout << "found " << nfunctions << " functions " << endl;
  plotfile.ignore(180, '\n'); //nfunctions

}

//----------------------------------------------------------------------

doublevar get_cutoff_radius(Array3 <doublevar> & grid, Array1 <doublevar> & box_size,
                       Array1 <doublevar> & origin, Array1 <doublevar> & center) {

  assert(box_size.GetDim(0)==3);
  assert(origin.GetDim(0)==3);
  assert(center.GetDim(0)==3);

  Array1 <int> npoints(3);
  for(int i=0; i< 3; i++) npoints(i)=grid.GetDim(i);

  Array2 <doublevar> latVec(3,3);
  for(int i=0; i< 3; i++)
    for(int j=0; j< 3; j++) latVec(i,j)=0;

  for(int i=0; i< 3; i++) latVec(i,i)=box_size(i);

  int nnorm=6;
  Array1 <doublevar> norm_radius(nnorm);
  norm_radius(0)=3.0;
  norm_radius(1)=4.0;
  norm_radius(2)=5.0;
  norm_radius(3)=7.0;
  norm_radius(4)=11.0;
  norm_radius(5)=6.0;

  Array1 <doublevar> partial_norm(nnorm);
  partial_norm=0;

  doublevar tot_norm=0;

  Array1 <doublevar> pos(3);


  doublevar delta0=box_size(0)/(npoints(0)-1);
  doublevar delta1=box_size(1)/(npoints(1)-1);
  doublevar delta2=box_size(2)/(npoints(2)-1);
  doublevar mo_val;
  doublevar norm;
  doublevar threshold=1e-15;
  for(int x=0; x< npoints(0); x++) {
    cout << "."; cout.flush();
    pos(0)=origin(0)+delta0*x;
    for(int y=0; y< npoints(1); y++) {
      pos(1)=origin(1)+delta1*y;
      for(int z=0; z< npoints(2); z++) {
        pos(2)=origin(2)+delta2*z;
        mo_val=grid(x,y,z);
        norm=mo_val*mo_val;
        if(norm> threshold) {
          doublevar dist=get_periodic_distance(pos, center, latVec);
          //cout << "dist " << dist << endl;
          tot_norm+=norm;
          for(int i=0; i< nnorm; i++) {
            //cout << "norm_radius " << norm_radius(i) << endl;
            if(dist < norm_radius(i)) {
              //cout << "included " << endl;
              partial_norm(i)+=norm;
            }
          }
        }
      }
    }
  }

  cout << endl;
  doublevar min_radius=box_size(0);

  cout << "total norm " << tot_norm << endl;
  for(int i=0; i < nnorm; i++) {
    doublevar percent=partial_norm(i)/tot_norm;
    cout << "radius " << norm_radius(i) << " fraction of norm " << percent << endl;
    if(percent > .996 && min_radius > norm_radius(i) ) min_radius=norm_radius(i);
  }
  cout << "choosing cutoff of  " << min_radius << endl;
  return min_radius;
}

//----------------------------------------------------------------------


void get_wannier_centers(istream & is, Array2 <doublevar> & centers) {

  vector < vector < doublevar> > centers_temp;

  string dummy;
  vector <doublevar> temp_vec;
  temp_vec.resize(3);
  while(is >> dummy) {
    if( dummy == "MLWF") {
      centers_temp.clear();
      is.ignore(180, '\n');
      is.ignore(180, '\n');
      is >> dummy;
      while(dummy == "&&") {
        is >> dummy;
        is >> temp_vec[0] >> temp_vec[1] >> temp_vec[2];
        is.ignore(180, '\n');
        is >> dummy;
        centers_temp.push_back(temp_vec);
      }

    }
  }

  int ncenters=centers_temp.size();
  centers.Resize(ncenters, 3);
  for(int i=0; i< ncenters; i++) {
    cout << "wannier center " << i << "   ";
    for(int d=0; d< 3; d++) {
      centers(i,d)=centers_temp[i][d];
      cout << centers(i,d) << "  ";
    }
    cout << endl;
  }

}


//----------------------------------------------------------------------

doublevar get_periodic_distance(const Array1 <doublevar> & r1,
				const Array1 <doublevar> & r2,
				const Array2 <doublevar> & latVec) {
  assert(r1.GetDim(0)==3);
  assert(r2.GetDim(0)==3);
  assert(latVec.GetDim(0)==3);
  assert(latVec.GetDim(1)==3);

  Array1 <doublevar> delta_r(3);
  Array1 <doublevar> delta_r2(3);

  for(int d=0; d< 3; d++) {
    delta_r(d)=r2(d)-r1(d);
  }
  int nlatvec=1;

  doublevar r=1e99;
  for(int kk=-nlatvec; kk <=nlatvec; kk++) {
    for(int jj=-nlatvec; jj <=nlatvec; jj++) {
      for(int ii=-nlatvec; ii <=nlatvec; ii++) {
        for(int d=0; d< 3; d++) {
          delta_r2(d)=delta_r(d)+kk*latVec(0,d)+jj*latVec(1,d)+ii*latVec(2,d);
        }
        r=min(sqrt(delta_r2(0)*delta_r2(0)+
                   delta_r2(1)*delta_r2(1)+
                   delta_r2(2)*delta_r2(2)), r);
      }
    }
  }

  return r;
}

//----------------------------------------------------------------------

//----------------------------------------------------------------------
#include "jeep_utils.h"

/*
  (format comment copied from T. Ogitsu and B. Militzer)
Jeep wf format

spin up data

int = 0
int nst
double[nst]             = occ[n], n=0, .. , nst-1   --number of electrons in each orbital
int ngw
complex<double>[nst*ngw]  = wf(t)
complex<double>[nst*ngw]  = wf(t-dt)

next lines only if nspin == 2:

spin down data

int = 0
int nst
double[nst]             = occ[n], n=0, .. , nst-1
int ngw
complex<double>[nst*ngw]  = wf(t)
complex<double>[nst*ngw]  = wf(t-dt)

*/
int read_int(FILE * file) {
  int i;
  fread(&i, sizeof(int), 1, file);
  return i;
}

double read_double(FILE * file) {
  double f;
  fread(&f, sizeof(double), 1, file);
  return f;
}

void skip(FILE * file, const int bytes) {
  int pos=ftell(file);
  fseek(file, pos+bytes, SEEK_SET);
}

void get_jeep_header(FILE * wfin, int & nst, int & ngw, Array1 <doublevar> & occupation) {
   int test=read_int(wfin);
   if(test != 0) error("First int in wf file should be zero");
   nst=read_int(wfin);
   occupation.Resize(nst);
   for(int i=0; i< nst; i++) occupation(i)=read_double(wfin);
   ngw=read_int(wfin);
}

//----------------------------------------------------------------------

void summarize_jeep_wf(string & filename, int & nst, int & ngw, 
                       Array1 <int> & mo_places) {

  FILE * wffile= fopen(filename.c_str(), "r");
  if(wffile==0) { 
    cout << "Couldn't open " << filename << endl;
    exit(1);
  }
  
  int startpos=ftell(wffile);
  fseek(wffile, 0L, SEEK_END);
  int size=ftell(wffile);
  fseek(wffile, startpos, SEEK_SET);

  int test=read_int(wffile); //should be zero
  nst=read_int(wffile);

  for(int i=0; i< nst; i++ )
    read_double(wffile); //clear out the occupation; we don't need it
  
  ngw=read_int(wffile);

  int size_of_mo=2*ngw*sizeof(double);
  int nskip=4*nst*ngw*sizeof(double);
  //cout << "nskip " << nskip << " size " << size << endl;
  skip(wffile, nskip);
  int pos=ftell(wffile);
  //cout << "position " << pos << endl;

  if(pos==size) {
    fclose(wffile);
    mo_places.Resize(nst);
    int start=3*sizeof(int)+nst*sizeof(double);
    for(int i=0; i< nst; i++) {
      mo_places(i)=start+i*size_of_mo;
    }
    return;
  }

  if(pos > size) {
    cout << "WARNING: count seems off in wf file " << endl;
  }



  //if we're this far, then there's a second spin
  test=read_int(wffile);
  assert(test==0);
  int nst2=read_int(wffile);

  mo_places.Resize(nst+nst2);
  int start=3*sizeof(int)+nst*sizeof(double);
  int start2=3*sizeof(int)+nst2*sizeof(double);
  for(int i=0; i< nst; i++) {
    mo_places(i)=start+i*size_of_mo;
    //cout << "mo_places " << mo_places(i) << endl;
  }  
  for(int i=nst; i< nst+nst2; i++) {
    mo_places(i)=start //first header
      + nst*size_of_mo //the second wf of spin 0
      + start2         //second header
      + i*size_of_mo;  //first wf of spin 0+current mo
    //cout << "mo_places " << mo_places(i) << endl;
  }
 
  nst+=nst2;

  fclose(wffile);
  return;
  
}


//----------------------------------------------------------------------

void get_next_pw_mo(FILE * wfin, int ngw,  Array2 <doublevar> & coeff) {
  coeff.Resize(ngw,2);
  for(int j=0; j < ngw; j++) {
    coeff(j,0)=read_double(wfin);
    coeff(j,1)=read_double(wfin);
  }
}


