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


#include "Qmc_std.h"
#include "Center_set.h"
#include "Basis_function.h"


//------------------------------------------------------------------------

void divide_job(int node, int nprocs, int npoints, int & start, int & end) {
  start=node*npoints/nprocs;
  end=(node+1)*npoints/nprocs;
}


//------------------------------------------------------------------------
void fit_molecorb(Center_set & centers,
                  Array1 <Basis_function * > & basis,
                  Sample_point * sample,
                  Array1 <doublevar> & origin,
                  Array1 <doublevar> & box_size,
                  Array1 <int> npoints,
                  Array2 <doublevar> & basis_matrix_inverse,
                  vector <string> & file_list,
                  Array2 <doublevar> & coefficients) {

  int totbasis=0;
  for(int i=0; i< centers.size(); i++)
  {
    int basiscent=0;
    for(int j=0; j< centers.nbasis(i); j++)
    {
      basiscent+=basis(centers.basis(i,j))->nfunc();
      //cout << "basiscent " << basiscent << endl;
    }
    totbasis+=basiscent;
  }

  int nmo_fit=file_list.size();


  Array1 <doublevar> projection(totbasis);
  coefficients.Resize(nmo_fit, totbasis);
  coefficients=0;

  Array3 <doublevar> mo_grid;



  Array1 <doublevar> pos(3);
  Basis_function * tempbasis=NULL;
  Array1 <doublevar> basisvals(totbasis);
  Array1 <doublevar> R(5);
  doublevar delta0=box_size(0)/(npoints(0)-1);
  doublevar delta1=box_size(1)/(npoints(1)-1);
  doublevar delta2=box_size(2)/(npoints(2)-1);

  int mo_start=0;
  int mo_end=nmo_fit;

#ifdef USE_MPI
  divide_job(mpi_info.node, mpi_info.nprocs, nmo_fit, mo_start, mo_end);
  cout << mpi_info.node << " mos: " << mo_start << " - " << mo_end << endl;
#endif

  for(int mo=mo_start; mo < mo_end; mo++) {
    cout << mpi_info.node << ":fitting mo " << mo << endl;

    ifstream orbin(file_list[mo].c_str());
    //cout << "reading grid " << endl;
    get_grid_from_plot(orbin, mo_grid, box_size, origin);

    projection=0;

    int pointnum=0;

    doublevar threshold=1e-6;

    for(int x=0; x< npoints(0); x++) {
      cout << "."; cout.flush();
      pos(0)=origin(0)+delta0*x;

      for(int y=0; y< npoints(1); y++) {
      pos(1)=origin(1)+delta1*y;

      for(int z=0; z< npoints(2); z++) {
      pos(2)=origin(2)+delta2*z;



      if(fabs(mo_grid.v[pointnum]) > threshold) {

        sample->setElectronPos(0, pos);
        centers.updateDistance(0, sample);
        int currfunc=0;
        for(int cent=0; cent< centers.size(); cent++)
        {
          centers.getDistance(0,cent , R);
          for(int j=0; j< centers.nbasis(cent); j++)
          {
            tempbasis=basis(centers.basis(cent,j));
            tempbasis->calcVal(R, basisvals, currfunc);
            currfunc+=tempbasis->nfunc();

          }
        }

        for(int bas=0; bas< totbasis; bas++) {
          projection(bas)+=mo_grid.v[pointnum]*basisvals(bas);
        }
      }
      pointnum++;
     }
     }
   }

   cout << endl;


   for(int bas=0; bas < totbasis; bas++) {
    cout << "projection " << bas << "   " << projection(bas)*delta0*delta1*delta2 << endl;
   }

    //now multiply the projection by the inverse to get the coefficients

    for(int bas=0; bas < totbasis; bas++) {
      for(int bas2=0; bas2 < totbasis; bas2++) {
        coefficients(mo, bas)+=basis_matrix_inverse(bas, bas2)*projection(bas2);
      }
    }


    orbin.close();


  }

#ifdef USE_MPI
  for(int mo=0; mo < nmo_fit; mo++) {
  
    int bcast_node=0;
    for(int n=0; n< mpi_info.nprocs; n++) {
      int node_start;//=n*nmo_fit/mpi_info.nprocs;
      int node_end;//=(n+1)*nmo_fit/mpi_info.nprocs;
      divide_job(n, mpi_info.nprocs, nmo_fit, node_start, node_end);
      if(mo >= node_start && mo < node_end) {
        bcast_node=n;
        break;
      }
    }

    if(bcast_node==mpi_info.node)
      cout << "broadcasting " << mo << " "  << bcast_node << endl;
    MPI_Bcast(coefficients.v+mo*totbasis, totbasis, MPI_DOUBLE, bcast_node,
              MPI_COMM_WORLD);
  }


#endif


}

//------------------------------------------------------------------------


void calculate_basis_overlap(Center_set & centers,
                             Array1 <Basis_function* > & basis,
                             Sample_point * sample,
                             Array1 <doublevar> & origin,
                             Array1 <doublevar> & box_size,
                             Array1 <int> & npoints,
                             Array2 <doublevar> & basis_matrix) {


  int totbasis=0;
  for(int i=0; i< centers.size(); i++)
  {
    int basiscent=0;
    for(int j=0; j< centers.nbasis(i); j++)
    {
      basiscent+=basis(centers.basis(i,j))->nfunc();
      //cout << "basiscent " << basiscent << endl;
    }
    totbasis+=basiscent;
  }

  Array1 <doublevar> pos(3);
  Basis_function * tempbasis=NULL;
  Array1 <doublevar> basisvals(totbasis);
  Array1 <doublevar> R(5);
  int pointnum=0;
  doublevar delta0=box_size(0)/(npoints(0)-1);
  doublevar delta1=box_size(1)/(npoints(1)-1);
  doublevar delta2=box_size(2)/(npoints(2)-1);

  cout << "spacing (x, y, z) :  " << delta0 << "   " << delta1 << "   " << delta2 << endl;

  doublevar basis_threshold=1e-8;
  basis_matrix.Resize(totbasis, totbasis);
  basis_matrix=0;
  cout << "calculating basis functions";
  cout.flush();


  int xstart=0;
  int xend=npoints(0);
#ifdef USE_MPI
  divide_job(mpi_info.node, mpi_info.nprocs, npoints(0), xstart, xend);
  //cout << "xstart for " << mpi_info.node << " : " << xstart << endl;
  //cout << "xend for " << mpi_info.node << " : " << xend << endl;
#endif
  for(int x=xstart; x< xend; x++) {
    cout << "."; cout.flush();
    pos(0)=origin(0)+delta0*x;
  //  cout << "x " << pos(0) << endl;
  for(int y=0; y< npoints(1); y++) {
    pos(1)=origin(1)+delta1*y;
   // cout << "y " << pos(1) << endl;
  for(int z=0; z< npoints(2); z++) {
    pos(2)=origin(2)+delta2*z;

    sample->setElectronPos(0, pos);
    centers.updateDistance(0, sample);
    int currfunc=0;
    for(int cent=0; cent< centers.size(); cent++)
    {
      centers.getDistance(0,cent , R);
      for(int j=0; j< centers.nbasis(cent); j++)
      {
        tempbasis=basis(centers.basis(cent,j));
        tempbasis->calcVal(R, basisvals, currfunc);
        currfunc+=tempbasis->nfunc();

      }
    }

    //Get overlap matrix
    for(int b1=0; b1 < totbasis; b1++) {
      if(fabs(basisvals(b1)) > basis_threshold)
      for(int b2=b1; b2 < totbasis; b2++) {
        basis_matrix(b1, b2)+=basisvals(b1)*basisvals(b2);
      }
    }
    pointnum++;

  }
  }
  }
  cout << endl;


  //this is a symmetric matrix..
  for(int b1=0; b1 < totbasis; b1++) {
    for(int b2=b1+1; b2 < totbasis; b2++) {
      basis_matrix(b2, b1)=basis_matrix(b1, b2);
    }
  }


#ifdef USE_MPI
  Array2 <doublevar> tempmatrix(totbasis, totbasis);
  tempmatrix=basis_matrix;
  MPI_Allreduce(tempmatrix.v, basis_matrix.v, totbasis*totbasis,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

}

//----------------------------------------------------------------------
