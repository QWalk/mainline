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

#include "integrals.h"
#include "qmc_io.h"
#include "MatrixAlgebra.h"
#include "Pbc_enforcer.h"
#include "jeep_utils.h"
#include "gaussian_set.h"
#include "overlaps.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

#ifdef USE_LAPACK
extern "C" {
#include "clapack.h"
#include "cblas.h"
}
#endif



#include "D3vector.h"
#include "JeepBasis.h"

void get_g_vec(double ecut, Array2 <doublevar> & latvec,
               Array2 <doublevar> & gvec) {

  if(latvec(0,1) != 0 || latvec(0,2) != 0 || latvec(1,0) != 0
     || latvec(1,2) != 0 || latvec(2,0) != 0 || latvec(2,1) !=0) {
     error("must have an orthogonal cell");
  }

  D3vector cell(latvec(0,0), latvec(1,1), latvec(2,2));

  RealBasis jeepbasis;
  bool status=jeepbasis.resize(cell, cell, ecut);
  if(! status) error("Basis allocation failed");

  int ngw=jeepbasis.size();
  gvec.Resize(ngw,3);
  single_write(cout, ngw, " plane waves from ecut \n");
  for(int i=0; i< ngw; i++) {
    gvec(i,0)=jeepbasis.gx[i];
    gvec(i,1)=jeepbasis.gx[ngw+i];
    gvec(i,2)=jeepbasis.gx[ngw+ngw+i];
  }

}




//----------------------------------------------------------------------

void divide_job(int node, int nprocs, int npoints, int & start, int & end) {
  start=node*npoints/nprocs;
  end=(node+1)*npoints/nprocs;
}

//------------------------------------------------------------------------


//send a matrix to all nodes assuming that this one is
//filled in places given by divide_job
void communicate_matrix(Array2 <doublevar> & matrix) {
  int tot=matrix.GetDim(0);
  int indx2=matrix.GetDim(1); 
#ifdef USE_MPI
  for(int i=0;i < tot;  i++) {
    int bcast_node=0;
    for(int n=0; n< mpi_info.nprocs; n++) {
      int node_start;
      int node_end;
      divide_job(n, mpi_info.nprocs, tot, node_start, node_end);
      if(i >= node_start && i < node_end) {
        bcast_node=n;
        break;
      }
    }

    //if(bcast_node==mpi_info.node)
    //  cout << "broadcasting " << mo << " "  << bcast_node << endl;
    MPI_Bcast(matrix.v+i*indx2, indx2, MPI_DOUBLE, bcast_node,
              MPI_COMM_WORLD);
  }
#endif
}

//send a matrix to all nodes assuming that this one is
//filled in places given by divide_job
void communicate_matrix(Array2 <doublevar> & partial, 
                        Array2 <doublevar> & full, int tot) {

  int start, end;
  divide_job(mpi_info.node, mpi_info.nprocs, tot, start, end);

  int n2=partial.GetDim(1);
  full.Resize(tot, partial.GetDim(1));

  for(int i=start; i < end; i++) {
    for(int j=0; j < n2; j++) 
      full(i,j)=partial(i-start, j);
  }

  communicate_matrix(full);

}

//----------------------------------------------------------------------


#include <ctime>

//----------------------------------------------------------------------
void write_orb_file(string & orbname,
                    Array1 <Center> & centers,
                    Gaussian_lookups & lookup,
                    Array2 <doublevar> & mo_coeff) {
  int nst=mo_coeff.GetDim(0);
  int ncenters=centers.GetDim(0);

  ofstream orbout(orbname.c_str());
  int totcounter=0;
  for(int mo=0; mo < nst; mo++) {
    for(int cen=0; cen < ncenters; cen++) {
      int nfunc=centers(cen).basis.GetDim(0);
      for(int f=0; f < nfunc; f++) {
        orbout << "   " << mo+1 << "   " << f+1
               << "  "  << cen+1 << "   " << totcounter+1 << endl;
        totcounter++;
      }
    }
  }


  orbout << "COEFFICIENTS" << endl;
  for(int mo=0; mo < nst; mo++) {
    int counter=0;
    for(int cen=0; cen < ncenters; cen++) {
      int ecen=lookup.equivalent_center(cen);
      for(int i=lookup.center_start(ecen); i < lookup.center_end(ecen); i++)  {
       
        orbout << mo_coeff(mo,i) << "   ";
        if(counter%5==4) orbout << endl;
        counter++;
      }
    }
    orbout << endl;
  }
  orbout.close();

}

//----------------------------------------------------------------------

void serial_solve(Array2 <doublevar> & lcao_overlap, 
                  Array2 <doublevar> & pw_lap,
                  Array1 <int> & mos, int nst, //_total_ number of states
                  Array2 <doublevar> & mo_coeff) {

  int totbasis=lcao_overlap.GetDim(1);
  //--------------------------------------serial LU decomp
  clock_t invert_start=clock();

  Array2 <doublevar> lcao_full;
  communicate_matrix(lcao_overlap, lcao_full, totbasis);
  //Enforce symmetric matrix
  for(int b1=0; b1 < totbasis; b1++) {
    for(int b2=0; b2 < b1; b2++) {
      lcao_full(b1,b2)=lcao_full(b2,b1);
    }
  }

  Array1 <int> index(totbasis);
  double d;

#ifdef USE_LAPACK
  //----LAPACK stuff
  clapack_dgetrf(CblasRowMajor, totbasis, totbasis, lcao_full.v,
                 totbasis, index.v);
  //-----------------
#else 
  ludcmp(lcao_full, totbasis, index, d);
#endif

  clock_t invert_end=clock();
  single_write(cout,"inverting matrix took "
               , ((double) invert_end-invert_start)/( (double) CLOCKS_PER_SEC)
               ," seconds\n");

  mo_coeff.Resize(nst, totbasis);


  Array1 <doublevar> tmp(totbasis);
  for(int mo=mos(0); mo < mos(1); mo++) {
    int moi=mo-mos(0);
    for(int b1=0; b1 < totbasis; b1++) 
      tmp(b1)=pw_lap(moi, b1);

#ifdef USE_LAPACK 
    int nrhs=1;
    clapack_dgetrs(CblasRowMajor, CblasNoTrans, totbasis, nrhs,
                   lcao_full.v, totbasis, index.v, tmp.v,totbasis);
#else 
    lubksb(lcao_full, totbasis, index,tmp);
#endif                   


    doublevar overlap=0;
    for(int b1=0; b1 < totbasis; b1++)
      overlap+=tmp(b1)*pw_lap(moi, b1);
    cout << mo << " :overlap " << overlap << endl;

    for(int b1=0; b1 < totbasis; b1++)
      mo_coeff(mo, b1)=tmp(b1);
  }    

  communicate_matrix(mo_coeff);
}



//----------------------------------------------------------------------
void fit_mo(string & wfin, Array2 <doublevar> & gvec,
            Array2 <doublevar> & latvec, Array1 <doublevar> & origin,
            Array1 <Center> & centers, Array1 <Contracted_gaussian> & basis,
            Array2 <doublevar> & mo_coeff) {

  Gaussian_lookups lookup;
  lookup.set_lookup_tables(latvec, origin, centers);

  int totbasis=lookup.totbasis2cen.GetDim(0);

  //-----------------------lcao overlaps with each other

  Array2 <doublevar> lcao_overlap(totbasis, totbasis);
  clock_t overlap_start=clock();
  calculate_overlap(latvec, origin, centers, basis,
                    lookup,lcao_overlap);
  clock_t overlap_end=clock();
  cout << mpi_info.node << " : gaussian overlaps took " 
       << ((double) overlap_end-overlap_start)/( double(CLOCKS_PER_SEC))
       << " seconds " << endl;


  //-----------------------lcao overlaps with the plane waves
  Array2 <doublevar> pw_lap;
  int nst, ngw; Array1 <int> mo_places;
  summarize_jeep_wf(wfin, nst, ngw, mo_places);
  Array1 <int> mos(2);
  mos(0)=0; mos(1)=nst;
#ifdef USE_MPI
  divide_job(mpi_info.node, mpi_info.nprocs, nst, mos(0), mos(1));
#endif

  pw_overlap(wfin, gvec, centers, basis, lookup, pw_lap, mos);


  //----------------------now solve the linear equation

  serial_solve(lcao_overlap, pw_lap, mos, nst, mo_coeff);
}



//----------------------------------------------------------------------

int main(int argc, char ** argv) {
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &(mpi_info.nprocs));
  MPI_Comm_rank(MPI_COMM_WORLD, &(mpi_info.node));
  cout << "processor " << mpi_info.node << " alive " << endl;
#endif

  if(argc < 2)
    error("usage: pw2lcao <filename>");

  ifstream in(argv[1]);

  string centerin;;
  string basisin;
  string wfin;

  double ecut=-1;
  Array2 <doublevar> latvec(3,3);
  latvec=0.0;
  Array1 <doublevar> origin(3);
  origin=0;
  string orbout="pw2lcao.orb";
  string test;
  string compare_orb;
  Array2 <doublevar> ref_latvec;


  while(in >> test) {
    if(test=="latvec") {
      for(int i=0; i<3; i++)
        for(int j=0; j< 3; j++) {
          if(!(in >> latvec(i,j)))
            error("Error reading in lattice vectors");
        }
    }
    else if(test=="ref_latvec") {
      ref_latvec.Resize(3,3);
      for(int i=0; i<3; i++)
        for(int j=0; j< 3; j++) {
          if(!(in >> ref_latvec(i,j)))
            error("Error reading in reference lattice vectors");
        }
    }
    else if(test=="origin") {
      for(int i=0; i< 3; i++) in >> origin(i);
    }
    else if(test=="ecut")
      in >> ecut;
    else if(test=="center_file")
      in >> centerin;
    else if(test=="basis_file")
      in >> basisin;
    else if(test=="wf_file")
      in >> wfin;
    else if(test=="orb_file")
      in >> orbout;
    else if(test=="compare_orb")
      in >> compare_orb;
    else
      error("Didn't understand ", test);
  }

  if(ref_latvec.GetDim(0)==0) {
    ref_latvec.Resize(3,3);
    for(int i=0; i< 3; i++) 
      for(int j=0; j< 3; j++)
        ref_latvec(i,j)=latvec(i,j);
  }

  Array1 <Center> centers;
  Array1 <Contracted_gaussian> basis;

  //cout << mpi_info.node << " Creating basis " << endl;
  create_local_basis(centerin, basisin, centers, basis);

  //int nst, ngw;
  //cout << mpi_info.node << " getting g-vectors " << endl;
  Array2 <doublevar> gvec;
  get_g_vec(ecut, ref_latvec, gvec);


  //-----------------------
  Array2 <doublevar> mo_coeff;
  fit_mo(wfin, gvec, latvec, origin, centers, basis, mo_coeff);


  Gaussian_lookups lookup;
  lookup.set_lookup_tables(latvec, origin, centers);
  if(mpi_info.node==0) {
    write_orb_file(orbout, centers,lookup,mo_coeff);
  }


  if(compare_orb != "" && mpi_info.node ==0) {
    Array2 <doublevar> compare_coeff;
    int nmo=mo_coeff.GetDim(0);
    int nfunc=mo_coeff.GetDim(1);
    read_mo(compare_orb, nmo,nfunc,compare_coeff);
    Array1 <int> compare_list(nmo);

    for(int m=0; m< nmo; m++) compare_list(m)=m;

    //Assume diagonal overlap for the moment.
    Array2 <doublevar> lcao_overlap(nfunc, nfunc);
    lcao_overlap=0.0;
    for(int i=0; i< nfunc; i++)
      lcao_overlap(i,i)=1.0;

    compare_mo(compare_coeff, mo_coeff, lcao_overlap, compare_list);
   }
#ifdef USE_MPI
  MPI_Finalize();
#endif

}


//----------------------------------------------------------------------

