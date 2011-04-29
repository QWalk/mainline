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

#ifdef COMPLEX_WF
#define EINSPLINE_NAME MO_matrix_Ceinspline
#else
#define EINSPLINE_NAME MO_matrix_einspline
#endif

//----------------------------------------------------------------------
void EINSPLINE_NAME::buildLists(Array1 <Array1 <int> > & occupations) { 
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
#ifdef COMPLEX_WF 
  BCtype_z bc;
#else
  BCtype_d bc;
#endif
  bc.lCode=PERIODIC; bc.rCode=PERIODIC;
  nmo_lists.Resize(nsplines);
  for(int s=0; s< nsplines; s++) { 
    nmo_lists(s)=occupations(s).GetDim(0);
#ifdef COMPLEX_WF
    spline(s)=create_multi_UBspline_3d_z(grids(0),grids(1),grids(2),
        bc,bc,bc,nmo_lists(s));
#else
    spline(s)=create_multi_UBspline_3d_d(grids(0),grids(1),grids(2),
        bc,bc,bc,nmo_lists(s));
#endif
    for(int i=0; i < nmo_lists(s); i++) { 
#ifdef COMPLEX_WF
      set_multi_UBspline_3d_z(spline(s),i,modata.v+occupations(s)(i)*ngridpts);
#else
      set_multi_UBspline_3d_d(spline(s),i,modata.v+occupations(s)(i)*ngridpts);
#endif
    }
  }

#endif
}
//----------------------------------------------------------------------

void EINSPLINE_NAME::read(vector <string> & words, unsigned int & startpos, System * sys) { 
#ifdef USE_EINSPLINE
  unsigned int pos=startpos;
  ndim=3;
//  string orbfile;
  sys->kpoint(kpoint);
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
#ifdef COMPLEX_WF
  for(complex_double *p=modata.v;p!=modata.v+totpts; p++) {
    is >> *p;
  }
#else
  for(doublevar *p=modata.v; p!=modata.v+totpts; p++){
    is >> *p;
  }
#endif
  is.close();

#endif //USE_EINSPLINE
}
//----------------------------------------------------------------------

int EINSPLINE_NAME::showinfo(ostream & os) { 
  os << "Einspline" << endl;
  os << "NMO " << nmo << endl;
  os << "ORBFILE " << orbfile << endl;
  return 1;

}
//----------------------------------------------------------------------

int EINSPLINE_NAME::writeinput(string & indent, ostream &os ) { 
  os << indent << "EINSPLINE_MO" << endl;
  os<< indent << "NMO " << nmo << endl;
  os<< indent << "ORBFILE " << orbfile << endl;
  return 1;
}
//----------------------------------------------------------------------

#ifdef COMPLEX_WF
void EINSPLINE_NAME::updateVal(Sample_point * sample,int e,int listnum,Array2 <dcomplex> & newvals) { 
#else
void EINSPLINE_NAME::updateVal(Sample_point * sample,int e,int listnum,Array2 <doublevar> & newvals) { 
#endif
#ifdef USE_EINSPLINE
  Array1 <doublevar> pos(ndim),u(ndim);
  sample->getElectronPos(e,pos);
  u=0;
  for(int d=0; d< ndim; d++) { 
    for(int d1=0; d1 < ndim; d1++) { 
      u(d)+=pos(d1)*latvecinv(d,d1);
    }
    u(d)-=floor(u(d));
  }
 
#ifdef COMPLEX_WF
  Array1 <dcomplex> vals(nmo_lists(listnum));
  eval_multi_UBspline_3d_z(spline(listnum),u(0),u(1),u(2),vals.v);
  doublevar kr=0;
  for(int d=0; d< ndim; d++) { 
    kr+=kpoint(d)*u(d);
  }
  dcomplex eikr=dcomplex(cos(pi*kr),sin(pi*kr));
  //cout << u(0) << " " << u(1) << " " << u(2) << " " << vals(0) << " " <<  eikr << endl;
  //cout << u(2) << " " << sqrt(vals(15).real()*vals(15).real()+vals(15).imag()*vals(15).imag())<<  " " << vals(15).real() << " " << vals(15).imag() << " ";
  
  for(int i=0; i< nmo_lists(listnum); i++) {
    vals(i)*=eikr;
  }
  //cout << vals(15).real() << " " << vals(15).imag() << endl;
 // cout << pos(2) << " " << vals(15).real() << " " << vals(15).imag() << endl;
#else 
  Array1 <doublevar> vals(nmo_lists(listnum));
  eval_multi_UBspline_3d_d(spline(listnum),u(0),u(1),u(2),vals.v);
#endif

  for(int i=0; i< nmo_lists(listnum); i++) { 
    newvals(i,0)=vals(i);
  }

#endif
}
//----------------------------------------------------------------------

#ifdef COMPLEX_WF
void EINSPLINE_NAME::updateLap(Sample_point * sample,int e,int listnum,Array2 <dcomplex> & newvals) { 
#else
void EINSPLINE_NAME::updateLap(Sample_point * sample,int e,int listnum,Array2 <doublevar> & newvals) { 
#endif
#ifdef USE_EINSPLINE
  Array1 <doublevar> pos(ndim),u(ndim);
  sample->getElectronPos(e,pos);
  u=0;
  for(int d=0; d< ndim; d++) { 
    for(int d1=0; d1 < ndim; d1++) { 
      u(d)+=pos(d1)*latvecinv(d,d1);
    }
    u(d)-=floor(u(d));
  }
#ifdef COMPLEX_WF
  Array1 <dcomplex> vals(nmo_lists(listnum));
  Array2 <dcomplex> grad(nmo_lists(listnum),ndim);
  Array3 <dcomplex> hess(nmo_lists(listnum),ndim,ndim);
  eval_multi_UBspline_3d_z_vgh(spline(listnum),u(0),u(1),u(2),vals.v,grad.v,hess.v);

  doublevar kr=0;
  for(int d=0; d< ndim; d++) { 
    kr+=kpoint(d)*u(d);
  }
  dcomplex eikr=dcomplex(cos(pi*kr),sin(pi*kr));
  dcomplex tmp_val;
  Array1 <dcomplex> tmp_grad(ndim);
  dcomplex I(0,1.0);
  for(int i=0; i< nmo_lists(listnum); i++) { 
    tmp_val=vals(i);
    for(int d=0; d< ndim; d++) tmp_grad(d)=grad(i,d);
    vals(i)=eikr*tmp_val;
    for(int d=0; d< ndim; d++) 
      grad(i,d)=eikr*(I*pi*kpoint(d)*tmp_val+tmp_grad(d));
    for(int d1=0; d1 < ndim; d1++) 
      for(int d2=0; d2< ndim; d2++) 
        hess(i,d1,d2)=eikr*(hess(i,d1,d2)
            +I*pi*kpoint(d1)*tmp_grad(d2)
            +I*pi*kpoint(d2)*tmp_grad(d1)
            -pi*pi*kpoint(d1)*kpoint(d2)*tmp_val);
  }

#else
  Array1 <doublevar> vals(nmo_lists(listnum));
  Array2 <doublevar> grad(nmo_lists(listnum),ndim);
  Array3 <doublevar> hess(nmo_lists(listnum),ndim,ndim);
  eval_multi_UBspline_3d_d_vgh(spline(listnum),u(0),u(1),u(2),vals.v,grad.v,hess.v);
#endif
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
#ifdef COMPLEX_WF
    dcomplex lap(0.,0.);
#else
    doublevar lap=0;
#endif
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


#undef EINSPLINE_NAME
