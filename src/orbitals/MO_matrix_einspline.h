/*
 
Copyright (C) 2011 Lucas K. Wagner (based on work by Michal Bajdich)

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

#ifndef MO_MATRIX_EINSPLINE_H_INCLUDED
#define MO_MATRIX_EINSPLINE_H_INCLUDED

#include "MO_matrix.h"
#include "Array45.h"
// USE_RESTRICT with -std=gnu++98 flag for g++; needed for bspline.h 
#ifdef USE_RESTRICT
#define restrict __restrict__
#else
#define restrict 
#endif

// USE_EINSPLINE to add einspline library; tested on einspline-0.8.2
#ifdef USE_EINSPLINE 
#include <bspline.h>
#include <multi_bspline.h>
#endif

/*! This is a simple adaptor class to translate between the template language
 * and the C-type language of einspline
 * */
#ifdef USE_EINSPLINE
template <class T> class Spline_evaluator { 
  public:
    void create(Ugrid & gridx, Ugrid & gridy, Ugrid & gridz, int nspline) { } 
    void set(int i, T * data) { }
    void val(doublevar x, doublevar y, doublevar z, T * vals) { } 
    void hess(doublevar x, doublevar y, doublevar z, T* vals, T * grad, T * hess) { } 
};

template <> class Spline_evaluator<doublevar> { 
private:
  multi_UBspline_3d_d * spline;
public:
  Spline_evaluator() { spline=NULL; } 
  void create(Ugrid &  gridx, Ugrid & gridy, Ugrid & gridz,int nspline) { 
    BCtype_d bc;
    bc.lCode=PERIODIC; bc.rCode=PERIODIC;
    spline=create_multi_UBspline_3d_d(gridx,gridy,gridz,bc, bc, bc,nspline);
  }
  void set(int i, doublevar * data) { 
    set_multi_UBspline_3d_d(spline,i,data);
  }
  void val(doublevar x, doublevar y, doublevar z,doublevar * vals) { 
    eval_multi_UBspline_3d_d(spline,x,y,z,vals);
  }
  void hess(doublevar x, doublevar y, doublevar z, doublevar * vals, 
      doublevar * grad, doublevar * hess) { 
    eval_multi_UBspline_3d_d_vgh(spline,x,y,z,vals,grad,hess);
  }
  ~Spline_evaluator() { 
    if (spline) destroy_Bspline(spline);
  }
};


template <> class Spline_evaluator<dcomplex> { 
private:
  multi_UBspline_3d_z * spline;
public:
  Spline_evaluator() { spline=NULL; } 
  void create(Ugrid &  gridx, Ugrid & gridy, Ugrid & gridz,int nspline) { 
    BCtype_z bc;
    bc.lCode=PERIODIC; bc.rCode=PERIODIC;
    spline=create_multi_UBspline_3d_z(gridx,gridy,gridz,bc, bc, bc,nspline);
  }
  void set(int i, dcomplex * data) { 
    set_multi_UBspline_3d_z(spline,i,data);
  }
  void val(doublevar x, doublevar y, doublevar z,dcomplex * vals) { 
    eval_multi_UBspline_3d_z(spline,x,y,z,vals);
  }
  void hess(doublevar x, doublevar y, doublevar z, dcomplex * vals, 
      dcomplex * grad, dcomplex * hess) { 
    eval_multi_UBspline_3d_z_vgh(spline,x,y,z,vals,grad,hess);
  }
  ~Spline_evaluator() { 
    if (spline) destroy_Bspline(spline);
  }
};

#endif //USE_EINSPLINE


/*!
Represents a periodic set of orbitals using Ken Esler's EINSPLINE library. 
 */
template <class T> class MO_matrix_einspline:public Templated_MO_matrix<T> { 
protected:
  void init() { }
  using Templated_MO_matrix<T>::nmo;
  using Templated_MO_matrix<T>::orbfile;
private:
#ifdef USE_EINSPLINE
  Array1 <Spline_evaluator<T> > spline;
#endif
  Array2 <doublevar> latvec; //lattice vectors for the cell on which the function is defined
  Array2 <doublevar> latvecinv;
  Array1 <int> npoints;
  Array1 <doublevar> resolution;
  Array2 <doublevar> kpoint; //k-point for each orbital
  //Array4 <T > modata;
  Array1 <int> nmo_lists;
  int ndim;
  Array1 <Array1 <int> > occ;
public:
  virtual void buildLists(Array1 <Array1 <int> > & occupations);
  virtual void read(vector <string> & words, unsigned int & startpos, 
                    System * sys);
  virtual int showinfo(ostream & os);
  virtual int writeinput(string &, ostream &);
  virtual void writeorb(ostream &, Array2 <doublevar> & rotation, Array1 <int> & tmp) { } 
  virtual void updateVal(Sample_point *,int e,int listnum,Array2<T>&);
  virtual void updateLap(Sample_point *,int e, int listnum, Array2<T>&);
  virtual void updateHessian(Sample_point * sample,
			     int e, int listnum,Array2<T>&);

  MO_matrix_einspline() { } 
};


#include "qmc_io.h"
#include "MatrixAlgebra.h"


//----------------------------------------------------------------------
template <class T> void MO_matrix_einspline<T>::buildLists(Array1 <Array1 <int> > & occupations) { 
#ifdef USE_EINSPLINE
  ifstream is(orbfile.c_str());
  string dummy;
  is >> dummy;
  while(dummy != "orbitals") is>>dummy;
  is.ignore(180,'\n');
  Array1 <T> orb_data(npoints(0)*npoints(1)*npoints(2));

  occ=occupations; 
  int nsplines=occupations.GetDim(0);
  spline.Resize(occupations.GetDim(0));
  int ngridpts=npoints(0)*npoints(1)*npoints(2);
  Array1 <Ugrid> grids(ndim);
  for(int d=0; d<ndim; d++) { 
    grids(d).start=0;
    grids(d).end=1.0;
    grids(d).num=npoints(d);
  }
  nmo_lists.Resize(nsplines);
  for(int s=0; s< nsplines; s++) { 
    nmo_lists(s)=occupations(s).GetDim(0);
    spline(s).create(grids(0),grids(1), grids(2),nmo_lists(s));
  }

  for(int mo=0; mo < nmo; mo++) { 
    if(mpi_info.node==0) {
      is.read((char*)(orb_data.v),sizeof(T)*ngridpts);
    }

    for(int s=0; s< nsplines; s++) { 
      for(int i=0; i < nmo_lists(s); i++) { 
        if(occupations(s)(i)==mo) { 
#ifdef USE_MPI
          int nvals=sizeof(T)/sizeof(double);
          if(mpi_info.node==0) { 
            for(int proc=1; proc < mpi_info.nprocs; proc++) { 
              MPI_Send(orb_data.v,ngridpts*nvals,MPI_DOUBLE,proc,0,MPI_Comm_grp);
            }
          }
          else { 
            MPI_Status status;
            MPI_Recv(orb_data.v,ngridpts*nvals,MPI_DOUBLE,0,0,MPI_Comm_grp,&status);
          }          
#endif 
          spline(s).set(i,orb_data.v);
        }
      }
    }
  }

#endif
}
//----------------------------------------------------------------------

template <class T> void MO_matrix_einspline<T>::read(vector <string> & words, unsigned int & startpos, System * sys) { 
#ifdef USE_EINSPLINE
  unsigned int pos=startpos;
  ndim=3;
  if(!readvalue(words,pos=startpos,orbfile,"ORBFILE")) 
    error("Need keyword ORBFILE..");
  double magnify=1.0;
  readvalue(words,pos=startpos,magnify,"MAGNIFY");
  //Should probably just make node0 read this and send to others over MPI..
  ifstream is(orbfile.c_str());
  if(!is) error("Couldn't open ",orbfile);
  string dummy;
  is.ignore(180,'\n'); //header
  is >> dummy;  //nmo
  int nmo_file;
  is >> nmo_file;
  nmo=nmo_file;
  //cout << "nmo " << nmo << endl;
  is.ignore(180,'\n'); //clear nmo line
  is.ignore(180,'\n'); //K-point line
  kpoint.Resize(nmo_file,ndim);
  Array2 <doublevar> tmp_kpt(nmo_file,ndim);
  for(int i=0; i< nmo_file; i++) { 
    for(int d=0; d< ndim; d++)  { 
      is >> tmp_kpt(i,d);
    }

  }
  is.ignore(180,'\n');

  is.ignore(180,'\n');
  latvec.Resize(ndim,ndim);
  for(int i=0; i< ndim; i++) { 
    for(int j=0; j< ndim; j++) {
      is >> latvec(i,j);
      //cout << latvec(i,j) << " ";
    }
  }
  latvecinv.Resize(ndim,ndim);
//  cout << "inverting " << endl;
  InvertMatrix(latvec,latvecinv,ndim);
//  cout << "done " << endl;

  kpoint=0.0;
  for(int i=0; i< nmo_file; i++) { 
    for(int d1=0; d1 < ndim; d1++) { 
      for(int d2=0; d2 < ndim; d2++) { 
        kpoint(i,d2)+=tmp_kpt(i,d1)*latvecinv(d2,d1);
      }
    }
  }
  is >> dummy;
  resolution.Resize(ndim);
  for(int i=0; i< ndim; i++) is >> resolution(i);
  is >> dummy;
  npoints.Resize(ndim);
  for(int i=0; i< ndim; i++) is >> npoints(i);
  is.ignore(180,'\n'); is.ignore(180,'\n');

  is.close();
#else
  error("Not compiled with EINSPLINE support!");
#endif //USE_EINSPLINE

}
//----------------------------------------------------------------------

template <class T> int MO_matrix_einspline<T>::showinfo(ostream & os) { 
  os << "Einspline" << endl;
  os << "NMO " << nmo << endl;
  os << "ORBFILE " << orbfile << endl;
  return 1;

}
//----------------------------------------------------------------------

template <class T>int MO_matrix_einspline<T>::writeinput(string & indent, ostream &os ) { 
  os << indent << "EINSPLINE_MO" << endl;
  os<< indent << "NMO " << nmo << endl;
  os<< indent << "ORBFILE " << orbfile << endl;
  return 1;
}
//----------------------------------------------------------------------

template <class T> void MO_matrix_einspline<T>::updateVal(Sample_point * sample,
    int e, int listnum, Array2 <T> & newvals) { 
#ifdef USE_EINSPLINE
  Array1 <doublevar> pos(ndim),u(ndim);
  sample->getElectronPos(e,pos);
  u=0; 
  for(int d=0; d< ndim; d++) { 
    for(int d1=0; d1 < ndim; d1++) { 
      u(d)+=pos(d1)*latvecinv(d1,d);
    }
    u(d)-=floor(u(d));
  }
 
  Array1 <T> vals(nmo_lists(listnum));
  spline(listnum).val(u(0),u(1),u(2),vals.v);
  
  for(int i=0; i< nmo_lists(listnum); i++)  { 
    doublevar kr=0;
    for(int d=0; d< ndim; d++) { 
      kr+=kpoint(occ(listnum)(i),d)*pos(d);
    }
    T kfac=eval_kpoint_fac<T>(kr);
    vals(i)*=kfac;
  }
  
  for(int i=0; i< nmo_lists(listnum); i++) { 
    newvals(i,0)=vals(i);
  }
#endif
}
//----------------------------------------------------------------------

template <class T> void MO_matrix_einspline<T>::updateLap(Sample_point * sample,
    int e,int listnum,Array2 <T> & newvals) {
#ifdef USE_EINSPLINE
  Array2 <T> hessvals(nmo_lists(listnum),1+ndim+ndim*(ndim+1)/2);
  updateHessian(sample,e,listnum,hessvals);
  for(int i=0; i < nmo_lists(listnum); i++) {
    for (int d=0; d <= ndim; d++) {
      newvals(i,d)=hessvals(i,d);
    }
    newvals(i,ndim+1)=hessvals(i,ndim+1);
    for (int d=ndim+2; d <= 2*ndim; d++) {
      newvals(i,ndim+1)+=hessvals(i,d);
    }
  }
#endif
}


//----------------------------------------------------------------------

template <class T> void MO_matrix_einspline<T>::updateHessian(Sample_point * sample,
    int e,int listnum,Array2 <T> & newvals) { 
#ifdef USE_EINSPLINE
  Array1 <doublevar> pos(ndim),u(ndim);
  sample->getElectronPos(e,pos);
  u=0;
  //cout << "pos " << pos(0) << " " << pos(1) << " " << pos(2) << endl;
  for(int d=0; d< ndim; d++) { 
    for(int d1=0; d1 < ndim; d1++) { 
      u(d)+=pos(d1)*latvecinv(d1,d);
    }
    u(d)-=floor(u(d));
  }
  //cout << "u " << u(0) << " " << u(1) << " " << u(2) << endl;

  Array1 <T> vals(nmo_lists(listnum));
  Array2 <T> grad(nmo_lists(listnum),ndim);
  Array3 <T> hess(nmo_lists(listnum),ndim,ndim);
  spline(listnum).hess(u(0),u(1),u(2),vals.v,grad.v,hess.v);


  Array1 <T> tmp_grad(ndim);
  Array2 <T> tmp_hess(ndim,ndim);
  Array1 <doublevar> tmp_kpt(ndim);
  for(int i=0; i< nmo_lists(listnum); i++) { 
    doublevar kr=0;
    for(int d=0; d< ndim; d++) { 
      kr+=kpoint(occ(listnum)(i),d)*pos(d);
      tmp_kpt(d)=kpoint(occ(listnum)(i),d);
    }
   // cout << i << " kpt " << tmp_kpt(0) << " " << tmp_kpt(1) << " " << tmp_kpt(2) << endl;
    
    tmp_grad=T(0.0);
    tmp_hess=T(0.0);
    for(int d1=0; d1 < ndim; d1++) { 
      for(int d2=0; d2< ndim; d2++) { 
        tmp_grad(d1)+=latvecinv(d1,d2)*grad(i,d2);
      }
    }
    for(int d1=0; d1 < ndim; d1++) {
      for(int d2=0; d2 < ndim; d2++) {
        for(int d3=0; d3 < ndim; d3++) { 
          for(int d4=0; d4 < ndim; d4++) { 
            //tmp_hess(d1,d2)+=hess(i,d1,d2)*latvecinv(d3,d1)*latvecinv(d4,d2);
            //tmp_hess(d1,d2)+=hess(i,d3,d4)*latvecinv(d3,d1)*latvecinv(d4,d2);
            tmp_hess(d1,d2)+=hess(i,d3,d4)*latvecinv(d1,d3)*latvecinv(d2,d4);
 
          }
        }
      }
    }
    //-----------testing
    
    /*
    T lap=0.0;
    for(int d1=0; d1< ndim; d1++) { 
      for(int d2=0; d2 < ndim; d2++) { 
        for(int d3=0; d3 < ndim; d3++) { 
          lap+=hess(i,d2,d3)*latvecinv(d2,d1)*latvecinv(d3,d1);
        }
      }
    }

    T lap2=0;
    for(int d=0; d < ndim; d++) lap2+=tmp_hess(d,d);
    cout << i << " lap " << lap << " lap2 " << lap2 << endl;
    */
    //-----------

    eval_kpoint_deriv(tmp_kpt,kr,vals(i),tmp_grad,tmp_hess);
    for(int d1=0; d1 < ndim; d1++) {
      grad(i,d1)=tmp_grad(d1);
      for(int d2=0; d2 < ndim; d2++) 
        hess(i,d1,d2)=tmp_hess(d1,d2);
    }
  }
  /*
  for(int i=0; i< nmo_lists(listnum); i++) { 
    newvals(i,0)=vals(i);
  }
  for(int i=0; i< nmo_lists(listnum); i++) { 
    for(int d=0; d< ndim; d++) { 
      newvals(i,d+1)=0.0;
      for(int d1=0; d1 < ndim; d1++) { 
        newvals(i,d+1)+=latvecinv(d1,d)*grad(i,d1);
      }
    }
  }
*/
 
  //for(int i=0; i< nmo_lists(listnum); i++) { 
    //newvals(i,ndim+1)=lap;
  //}
  
  for(int i=0; i< nmo_lists(listnum); i++) { 
    newvals(i,0)=vals(i);
    int j=2*ndim+1;
    for(int d1=0; d1 < ndim; d1++) {
      newvals(i,d1+1)=grad(i,d1);
      newvals(i,d1+1+ndim)=hess(i,d1,d1);
      for(int d2=d1+1; d2 < ndim; d2++) {
        newvals(i,j)=hess(i,d1,d2);
        j++;
      }
    }
  }
#endif
}
//----------------------------------------------------------------------





#endif //MO_MATRIX_EINSPLINE_H_INCLUDED

