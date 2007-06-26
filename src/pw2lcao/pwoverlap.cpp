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
#include "Qmc_std.h"
#include "jeep_utils.h"
#include "gaussian_set.h"
#include "integrals.h"
#include "qmc_io.h"

void pw_overlap(string & wfin, Array2 <doublevar> & gvec,
                Array1 <Center> & centers,
                Array1 <Contracted_gaussian> & basis,
                Gaussian_lookups & lookup,
                Array2 <doublevar> & pw_lap, 
                Array1 <int> & mos) {

  assert(mos.GetDim(0)==2);

  int totbasis=lookup.totbasis2bas.GetDim(0);
  int ncenters=centers.GetDim(0);

  int nst, ngw;
  Array1 <int> mo_places;
  summarize_jeep_wf(wfin, nst, ngw, mo_places);

  FILE * wffile=fopen(wfin.c_str(), "r");

  single_write(cout, ngw, " g-vectors from wf file \n");
  if(ngw != gvec.GetDim(0) )
    error("Need the same number of g-vectors in wf file and from ecut");

  int mo_start=mos(0);
  int mo_end=mos(1);

  Array2 <doublevar> pw_coeff(ngw,2);

  pw_lap.Resize(mo_end-mo_start, totbasis);
  pw_lap=0;

  for(int mo=mo_start; mo < mo_end; mo++) {
    int moi=mo-mo_start;

    clock_t mo_start_time=clock();
    fseek(wffile, mo_places(mo), SEEK_SET);
    get_next_pw_mo(wffile, ngw, pw_coeff);

    Array1 <doublevar> sin_arr(ngw), cos_arr(ngw);
    for(int cen=0; cen < ncenters; cen++) {
     for(int i=1; i < ngw; i++) {
      doublevar scalarproduct=gvec(i,0)*centers(cen).pos(0)
                             +gvec(i,1)*centers(cen).pos(1)
                             +gvec(i,2)*centers(cen).pos(2);
#ifdef __USE_GNU
      sincos(scalarproduct, &(sin_arr(i)), &(cos_arr(i)));
#else
      cos_arr(i)=cos(scalarproduct);
      sin_arr(i)=sin(scalarproduct);
#endif
     }

      for(int b=lookup.center_start(cen); b < lookup.center_end(cen); b++) {
        int bas=lookup.totbasis2bas(b);
        int nalpha=basis(bas).alpha.GetDim(0);
        for(int a=0; a < nalpha; a++) {
          pw_lap(moi,b)+= basis(bas).coeff(a)
                          * P_x_c(basis(bas).lvals, centers(cen).pos,
                                  basis(bas).alpha(a), pw_coeff,
                                  gvec, sin_arr, cos_arr);
        }
      }
    }


    cout << mpi_info.node << ":calculating pw overlap for mo " << mo << "....";
    clock_t mo_end_time=clock();
    cout << "took "
         << ((double) mo_end_time-mo_start_time)/( (double) CLOCKS_PER_SEC)
         << " seconds" << endl;

  }

   fclose(wffile);
 }

//----------------------------------------------------------------------
