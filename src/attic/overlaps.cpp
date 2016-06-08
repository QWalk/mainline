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
#include "overlaps.h"
#include "integrals.h"

//----------------------------------------------------------------------

//Divide the job for the overlaps..
//NOTE: this is not the optimal division, since we're doing a triangular
//matrix.  Should be easy to change.
void divide_job_ov(int node, int nprocs, int npoints, int & start, int & end) {
  start=node*npoints/nprocs;
  end=(node+1)*npoints/nprocs;
}

//----------------------------------------------------------------------
#include <vector>
using namespace std;


void read_mo(string & orbin, int nmo, int nfunc, Array2 <doublevar> & moCoeff) {
  ifstream in(orbin.c_str());
  if(!in) error("Couldn't open ", orbin);
  string dummy;
  in >> dummy;
  while(dummy != "COEFFICIENTS")
    in >> dummy;

  moCoeff.Resize(nmo, nfunc);
  for(int m=0; m< nmo; m++) {
    for(int f=0; f< nfunc; f++) {
      in >> moCoeff(m,f);
     }
  }
}


/*!
Compare the sets of molecular orbitals, and
return true if they're the same, and false if not
Currently it doesn't use the lcao_overlap matrix, but
should be modified to do that; it's approximate otherwise.
*/
bool compare_mo(Array2 <doublevar> & oldMOCoeff,
                Array2 <doublevar> & newMOCoeff,
                Array2 <doublevar> & lcao_overlap,
                Array1 <int> compare_list) {

 assert(oldMOCoeff.GetDim(1)==newMOCoeff.GetDim(1));

 int nfunctions=oldMOCoeff.GetDim(1);

 int ncompare=compare_list.GetDim(0);
 if(newMOCoeff.GetDim(0) < ncompare)
   error("the punch file doesn't have enough molecular orbitals to compare.");
 if(oldMOCoeff.GetDim(0) < ncompare)
   error("the old punch file doesn't have enough MO's to compare.");

  vector <int> unresolved_mos;

  //First check to see if the mo's are in the same place
  //(most should be)
  for(int i=0; i< ncompare; i++) {
    int mo=compare_list(i);
    double dot=0, mag_old=0, mag_new=0;
    for(int f=0; f< nfunctions; f++) {
      dot+=newMOCoeff(mo,f)*oldMOCoeff(mo,f);
      mag_old+=oldMOCoeff(mo,f)*oldMOCoeff(mo,f);
      mag_new+=newMOCoeff(mo,f)*newMOCoeff(mo,f);
    }
    dot /= sqrt(mag_old*mag_new);
    dot =fabs(dot);
    cout << "mo " << mo << "  dot " << dot << endl;
    if(fabs(dot-1) > .01) {
      unresolved_mos.push_back(mo);
    }
  }

  int nunresolved=unresolved_mos.size();
  for(int i=0; i< nunresolved; i++) {
    cout << "not matched: " << unresolved_mos[i] << endl;
  }


  bool are_same=true;
  //See if any just swapped..
  for(int i=0; i< nunresolved; i++) {
    int mo1=unresolved_mos[i];
    bool resolved_swapping=false;

    for(int j=0; j< nunresolved; j++) {
      int mo2=unresolved_mos[j];

      double dot=0, mag_old=0, mag_new=0;
      for(int f=0; f< nfunctions; f++) {
        dot+=newMOCoeff(mo1,f)*oldMOCoeff(mo2,f);
        mag_old+=oldMOCoeff(mo2,f)*oldMOCoeff(mo2,f);
        mag_new+=newMOCoeff(mo1,f)*newMOCoeff(mo1,f);
      }
      dot /= sqrt(mag_old*mag_new);
      dot=fabs(dot);
      if(fabs(dot-1) < .01) {
        cout << "switched orbital: mo " << mo2 << " went to " << mo1
             << " dot product " << dot
             << endl;
        resolved_swapping=true;
      }
    }

    if(!resolved_swapping) {
      cout << "Unresolvable change in mo " << mo1 << endl;
      are_same=false;
    }

  }

  return are_same;
}


//----------------------------------------------------------------------

void calculate_overlap(Array2 <doublevar> & latvec,
                       Array1 <doublevar> & origin,
                       Array1 <Center> & centers,
                       Array1 <Contracted_gaussian> & basis,
                       Gaussian_lookups & lookup,
                       Array2 <doublevar> & lcao_overlap) {

  assert(lookup.totbasis2cen.GetDim(0)==lookup.totbasis2cen.GetDim(0));

  int totbasis=lookup.totbasis2cen.GetDim(0);
  //int ncenters=centers.GetDim(0);
  int nbasis=basis.GetDim(0);

  //Find the ranges of the basis functions
  //We calculate the range from viewing it as two
  //primitive S orbitals with parameter equal to
  //the smallest alpha orbital in the contraction
  const doublevar basis_tolerance=1e-10;
  Array1 <doublevar> minalpha(nbasis);
  for(int bas=0; bas < nbasis; bas++) {
    int nalpha=basis(bas).alpha.GetDim(0);
    minalpha(bas)=basis(bas).alpha(0);
    for(int a=1; a < nalpha; a++) {
      if(basis(bas).alpha(a) < minalpha(bas))
        minalpha(bas)=basis(bas).alpha(a);
    }
  }

  Array2 <doublevar> basis_range(nbasis, nbasis);
  for(int bas1=0; bas1< nbasis; bas1++) {
    for(int bas2=bas1; bas2 < nbasis; bas2++) {
      doublevar gamma=minalpha(bas1)+minalpha(bas2);
      basis_range(bas1, bas2)=-gamma
                           *log( (gamma/pi)*sqrt(gamma/pi)*basis_tolerance)
                           /(minalpha(bas1)*minalpha(bas2));
      //cout << sqrt(basis_range(bas1, bas2)) << "   ";
    }
    //cout << endl;
  }
  for(int bas1=0; bas1 < nbasis; bas1++) {
    for(int bas2=0; bas2 < bas1; bas2++) {
      basis_range(bas1, bas2)= basis_range(bas2, bas1);
    }
  }

  int basis_start=0;
  int basis_end=totbasis;
#ifdef USE_MPI
  divide_job_ov(mpi_info.node, mpi_info.nprocs, 
                totbasis, basis_start, basis_end);
#endif

  //lcao_overlap.Resize(totbasis, totbasis);
  lcao_overlap.Resize(basis_end-basis_start, totbasis);
  lcao_overlap=0;



  const int doublefac_size=10;
  //offset by one, so we can map the range -1 to inf to 0 to inf
  Array1 <int> doublefac_cache(doublefac_size);
  for(int i=0; i< doublefac_size; i++) {
    doublefac_cache(i)=doublefactorial(i-1);
  }

  const int binomial_size=8;
  Array2 <int> binomial_cache(binomial_size, binomial_size);
  for(int i=0; i< binomial_size; i++) {
    for(int j=0; j<=i; j++) {
      binomial_cache(i,j)=binomial2(i,j);
    }
  }

  //  int basis_start=0;
  //int basis_end=totbasis;
  //#ifdef USE_MPI
  //divide_job_ov(mpi_info.node, mpi_info.nprocs, totbasis, basis_start, basis_end);
  //#endif

  Array1 <doublevar> move_pos(3);
  for(int b1=basis_start; b1 < basis_end; b1++) {
    //if(b1%5==0) cout << "."; cout.flush();
    int offb=b1-basis_start;

    int bas1=lookup.totbasis2bas(b1);
    int cen1=lookup.totbasis2cen(b1);
    int nalpha1=basis(bas1).alpha.GetDim(0);

    for(int b2=b1; b2 < totbasis; b2++) {
      int bas2=lookup.totbasis2bas(b2);
      int cen2=lookup.totbasis2cen(b2);
      int nalpha2=basis(bas2).alpha.GetDim(0);


      doublevar overlap_tot=0;
      const int depth=4;
      for(int ii=-depth; ii <= depth; ii++) {
      for(int jj=-depth; jj <= depth; jj++) {
      for(int kk=-depth; kk <= depth; kk++) {

      for(int d=0; d< 3; d++)
        move_pos(d)=centers(cen2).pos(d)-centers(cen1).pos(d)
                    +ii*latvec(0,d)
                    +jj*latvec(1,d)
                    +kk*latvec(2,d);

      doublevar AB2=0.0;
      for(int d=0; d< 3; d++) AB2+=move_pos(d)*move_pos(d);

      if(AB2 < basis_range(bas1, bas2))
      for(int a1=0; a1 < nalpha1; a1++) {
        for(int a2=0; a2 < nalpha2; a2++) {
          //doublevar overlap=S12_f(basis(bas1).lvals, basis(bas2).lvals,
          //                    basis(bas1).alpha(a1), basis(bas2).alpha(a2),
          //                    centers(cen1).pos, move_pos);
          doublevar alpha1=basis(bas1).alpha(a1);
          doublevar alpha2=basis(bas2).alpha(a2);
          doublevar gamma=alpha1+alpha2;

          doublevar invgamma=1.0/gamma;
          doublevar PA, PB, dummyI=1.0;
          doublevar sqrtgamma=sqrt(pi/gamma);
          doublevar diff, tmpint, fsum;
          for(int d=0; d< 3; d++) {
            int l1=basis(bas1).lvals(d);
            int l2=basis(bas2).lvals(d);
            diff=move_pos(d);
            //AB2+=diff*diff;
            PA=diff*alpha2*invgamma;
            PB=-diff*alpha1*invgamma;
            tmpint=0.0;
            for(int j=0; j<= 0.5*(l1+l2); j++) {
              fsum=0.0;
              int qstart=max(-2*j, 2*(j-l2));
              int qend=min(2*j, 2*(l1-j));
              for(int q=qstart; q<= qend; q+=2) {
                int ti=int(j+0.5*q);
                int tj=int(j-0.5*q);
                //fsum+=binomial2(l1,ti)*binomial2(l2,tj)*pow(PA, l1-ti)
                //  *pow(PB, l2-tj);
                fsum+=binomial_cache(l1,ti)*binomial_cache(l2,tj)
                      *pow(PA, l1-ti)*pow(PB, l2-tj);
              }
              tmpint+=fsum
                      *doublefac_cache(2*j)
                      *pow(2.0*gamma, -j);
            }
            dummyI*=tmpint*sqrtgamma;
          }
          doublevar overlap=exp(-alpha1*alpha2*AB2*invgamma)*dummyI;
          overlap_tot+=basis(bas1).coeff(a1)*basis(bas2).coeff(a2)*overlap;
        }
      }
      }
      }
      }
      //cout << "overlap_tot " << overlap_tot << endl;
      lcao_overlap(offb, b2)=overlap_tot;

    }
  }
  //cout << endl;

  /*
#ifdef USE_MPI
  for(int b1=0; b1 < totbasis; b1++) {

    int bcast_node=0;
    for(int n=0; n< mpi_info.nprocs; n++) {
      int node_start;//=n*nmo_fit/mpi_info.nprocs;
      int node_end;//=(n+1)*nmo_fit/mpi_info.nprocs;
      divide_job_ov(n, mpi_info.nprocs, totbasis, node_start, node_end);
      if(b1 >= node_start && b1 < node_end) {
        bcast_node=n;
        break;
      }
    }

    //if(bcast_node==mpi_info.node)
    //  cout << "broadcasting " << mo << " "  << bcast_node << endl;
    MPI_Bcast(lcao_overlap.v+b1*totbasis, totbasis, MPI_DOUBLE, bcast_node,
              MPI_COMM_WORLD);
  }
#endif


  //Enforce symmetric matrix
  for(int b1=0; b1 < totbasis; b1++) {
    for(int b2=0; b2 < b1; b2++) {
      lcao_overlap(b1,b2)=lcao_overlap(b2,b1);
    }
  }
  */


}
//----------------------------------------------------------------------

