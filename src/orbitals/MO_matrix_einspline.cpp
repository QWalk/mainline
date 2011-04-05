/*
 
Copyright (C) 2011 Lucas K. Wagner

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

#include "MO_matrix_einspline.h"
#include "qmc_io.h"
#include "MatrixAlgebra.h"

//----------------------------------------------------------------------
void MO_matrix_einspline::buildLists(Array1 <Array1 <int> > & occupations) { 
#ifdef USE_EINSPLINE
  
  int nsplines=occupations.GetDim(0);
  spline.Resize(occupations.GetDim(0));
  int ngridpts=npoints(0)*npoints(1)*npoints(2);
  Array1 <Ugrid> grids(ndim);
  for(int d=0; d<ndim; d++) { 
    grids(d).start=0;
    grids(d).end=1.0;
    grids(d).num=npoints(d);
  }
  BCtype_d bc;
  bc.lCode=PERIODIC; bc.rCode=PERIODIC;
  nmo_lists.Resize(nsplines);
  for(int s=0; s< nsplines; s++) { 
    nmo_lists(s)=occupations(s).GetDim(0);
    spline(s)=create_multi_UBspline_3d_d(grids(0),grids(1),grids(2),bc,bc,bc,nmo_lists(s));
    for(int i=0; i < nmo_lists(s); i++) { 
      set_multi_UBspline_3d_d(spline(s),i,modata.v+occupations(s)(i)*ngridpts);
    }
  }

#endif
}
//----------------------------------------------------------------------

void MO_matrix_einspline::read(vector <string> & words, unsigned int & startpos, System * sys) { 
#ifdef USE_EINSPLINE
  unsigned int pos=startpos;
  ndim=3;
  string orbfile;
  if(!readvalue(words,pos=startpos,orbfile,"ORBFILE")) 
    error("Need keyword ORBFILE..");
  if(!readvalue(words,pos=startpos,nmo,"NMO"))
    error("Need keyword NMO");
  double magnify=1.0;
  readvalue(words,pos=startpos,magnify,"MAGNIFY");
  //Should probably just make node0 read this and send to others over MPI..
  ifstream is(orbfile.c_str());
  if(!is) error("Couldn't open ",orbfile);
  string dummy;
  is.ignore(180,'\n');
  is >> dummy; 
  int nmo_file;
  is >> nmo_file;
  is.ignore(180,'\n');
  is.ignore(180,'\n');
  latvec.Resize(ndim,ndim);
  for(int i=0; i< ndim; i++) { 
    for(int j=0; j< ndim; j++) {
      is >> latvec(i,j);
    }
  }
  latvecinv.Resize(ndim,ndim);
  InvertMatrix(latvec,latvecinv,ndim);
  is >> dummy;
  resolution.Resize(ndim);
  for(int i=0; i< ndim; i++) is >> resolution(i);
  is >> dummy;
  npoints.Resize(ndim);
  for(int i=0; i< ndim; i++) is >> npoints(i);
  is.ignore(180,'\n'); is.ignore(180,'\n');
  modata.Resize(nmo,npoints(0),npoints(1),npoints(2));
  int totpts=nmo*npoints(0)*npoints(1)*npoints(2);
  for(doublevar *p=modata.v; p!=modata.v+totpts; p++){
    is >> *p;
  }
  is.close();

#endif //USE_EINSPLINE
}
//----------------------------------------------------------------------

int MO_matrix_einspline::showinfo(ostream & os) { 
  return 1;
}
//----------------------------------------------------------------------

int MO_matrix_einspline::writeinput(string &, ostream &) { 
  return 1;
}
//----------------------------------------------------------------------


void MO_matrix_einspline::updateVal(Sample_point * sample,int e,int listnum,Array2 <doublevar> & newvals) { 
#ifdef USE_EINSPLINE
  Array1 <doublevar> vals(nmo_lists(listnum));
  Array1 <doublevar> pos(ndim),u(ndim);
  sample->getElectronPos(e,pos);
  u=0;
  for(int d=0; d< ndim; d++) { 
    for(int d1=0; d1 < ndim; d1++) { 
      u(d)+=pos(d1)*latvecinv(d,d1);
    }
    u(d)-=floor(u(d));
  }
  //cout << "pos " << pos(0) << " " << pos(1) << " " << pos(2) 
  //  << " u " << u(0) << " " << u(1) << " " << u(2) << endl;
  eval_multi_UBspline_3d_d(spline(listnum),u(0),u(1),u(2),vals.v);

  for(int i=0; i< nmo_lists(listnum); i++) { 
    newvals(i,0)=vals(i);
    //cout << "vals " << vals(i) << endl;
  }

#endif
}
//----------------------------------------------------------------------

  
void MO_matrix_einspline::updateLap(Sample_point * sample,int e,int listnum,Array2 <doublevar> & newvals) { 
#ifdef USE_EINSPLINE
  Array1 <doublevar> vals(nmo_lists(listnum));
  Array2 <doublevar> grad(nmo_lists(listnum),ndim);
  Array3 <doublevar> hess(nmo_lists(listnum),ndim,ndim);
  Array1 <doublevar> pos(ndim),u(ndim);
  sample->getElectronPos(e,pos);
  u=0;
  for(int d=0; d< ndim; d++) { 
    for(int d1=0; d1 < ndim; d1++) { 
      u(d)+=pos(d1)*latvecinv(d,d1);
    }
    u(d)-=floor(u(d));
  }
  eval_multi_UBspline_3d_d_vgh(spline(listnum),u(0),u(1),u(2),vals.v,grad.v,hess.v);
  for(int i=0; i< nmo_lists(listnum); i++) { 
    newvals(i,0)=vals(i);
    //cout << "lap vals " << vals(i) << endl;
  }
  for(int i=0; i< nmo_lists(listnum); i++) { 
    for(int d=0; d< ndim; d++) { 
      newvals(i,d+1)=0.0;
      for(int d1=0; d1 < ndim; d1++) { 
        newvals(i,d+1)+=latvecinv(d1,d)*grad(i,d1);
      }
    }
  }

  for(int i=0; i< nmo_lists(listnum); i++) { 
    doublevar lap=0;
    for(int d1=0; d1< ndim; d1++) { 
      for(int d2=0; d2 < ndim; d2++) { 
        for(int d3=0; d3 < ndim; d3++) { 
          lap+=hess(i,d2,d3)*latvecinv(d2,d1)*latvecinv(d3,d1);
        }
      }
    }
    newvals(i,ndim+1)=lap;
  }


#endif
}
//----------------------------------------------------------------------



