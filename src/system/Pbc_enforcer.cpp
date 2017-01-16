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

#include "Pbc_enforcer.h"

void Pbc_enforcer::init(Array2 <doublevar> & latVec) {
  assert(latVec.GetDim(0)==ndim);
  assert(latVec.GetDim(1)==ndim);

  _latVec.Resize(ndim, ndim);

  for(int i=0; i< ndim; i++) {
    for(int j=0; j< ndim; j++) {
      _latVec(i,j)=latVec(i,j);
    }
  }

  /*
  cout << "Lattice vectors:" << endl;
  for(int i=0; i< ndim; i++) {
    cout << i << " : ";
    for(int j=0; j < ndim; j++) {
      cout << latVec(i,j) << "   " ;
    }
    cout << endl;
  }
  */

  _origin.Resize(3);
  _origin=0;   //defaulting the origin to zero


  //-------------cross products

  //cross product:  0->1x2, 1->2x0, 2->0x1
  Array2 <doublevar> crossProduct(ndim, ndim);
  crossProduct(0,0)=(latVec(1,1)*latVec(2,2)-latVec(1,2)*latVec(2,1));
  crossProduct(0,1)=(latVec(1,2)*latVec(2,0)-latVec(1,0)*latVec(2,2));
  crossProduct(0,2)=(latVec(1,0)*latVec(2,1)-latVec(1,1)*latVec(2,0));

  crossProduct(1,0)=(latVec(2,1)*latVec(0,2)-latVec(2,2)*latVec(0,1));
  crossProduct(1,1)=(latVec(2,2)*latVec(0,0)-latVec(2,0)*latVec(0,2));
  crossProduct(1,2)=(latVec(2,0)*latVec(0,1)-latVec(2,1)*latVec(0,0));

  crossProduct(2,0)=(latVec(0,1)*latVec(1,2)-latVec(0,2)*latVec(1,1));
  crossProduct(2,1)=(latVec(0,2)*latVec(1,0)-latVec(0,0)*latVec(1,2));
  crossProduct(2,2)=(latVec(0,0)*latVec(1,1)-latVec(0,1)*latVec(1,0));

  //-------------------normal vectors

  _normVec.Resize(ndim, ndim);

  for(int i=0; i< ndim; i++) {
    for(int j=0; j < ndim; j++) {
      _normVec(i,j)=crossProduct(i,j);
    }

    //Check to make sure the direction is facing out
    doublevar dotprod=0;
    for(int j=0; j < ndim; j++) {
      dotprod+=_normVec(i,j)*latVec(i,j);
    }
    if(dotprod < 0) {
      for(int j=0; j< ndim; j++) {
        _normVec(i,j)= -_normVec(i,j);
      }
    }
  }

  _corners.Resize(ndim, ndim);
  for(int i=0; i< ndim; i++) {
    for(int j=0; j< ndim; j++) {
      _corners(i,j)=_origin(j)+latVec(i,j);
    }
  }


  init_called=true;
}

//----------------------------------------------------------------------
void Pbc_enforcer::setOrigin(Array1 <doublevar> & origin) {

  assert(origin.GetDim(0)==ndim);
  assert(init_called);
  _origin.Resize(ndim);
  for(int i=0; i< ndim; i++) _origin(i)=origin(i);

  for(int i=0; i< ndim; i++) {
    for(int j=0; j< ndim; j++) {
      _corners(i,j)=_origin(j)+_latVec(i,j);
    }
  }
}

//----------------------------------------------------------------------

int Pbc_enforcer::isInside(Array1 <doublevar>   & pos) {
  assert(pos.GetDim(0) >= 3);
  assert(init_called);
  for(int i=0; i< 3; i++) {
      //Whether we're past the origin side
      doublevar tooshort=0;
      for(int j=0; j<3;j++) tooshort+=_normVec(i,j)*(pos(j)-_origin(j));

      //Whether we're past the lattice vector
      doublevar toofar=0;
      for(int j=0; j< 3; j++) toofar+=_normVec(i,j)*(pos(j)-_corners(i,j));

      if(tooshort <0 || toofar > 0) return 0;
  }
  return 1;
}

//----------------------------------------------------------------------
int Pbc_enforcer::enforcePbc(Array1 <doublevar> & pos) {
  assert(pos.GetDim(0) >=3);
  assert(init_called);
  int shifted=0;
  for(int i=0; i< 3; i++) {
    int shouldcontinue=1;
    while(shouldcontinue) {

      //Whether we're past the origin side
      doublevar tooshort=0;
      for(int j=0; j<3;j++) tooshort+=_normVec(i,j)*(pos(j)-_origin(j));

      //Whether we're past the lattice vector
      doublevar toofar=0;
      for(int j=0; j< 3; j++) toofar+=_normVec(i,j)*(pos(j)-_corners(i,j));

      if(tooshort < 0) {
        for(int j=0; j< 3; j++) pos(j)+=_latVec(i,j);
      }
      else if(toofar >0) {
        for(int j=0; j< 3; j++) pos(j)-=_latVec(i,j);
      }
      else {
        shouldcontinue=0;
      }
    }
  }

  return shifted;

}

//----------------------------------------------------------------------
#include "vecmath.h"
double distance_edge(const Array1 <doublevar> & pos, 
                     const Array1 <doublevar> & a,
                     const Array1 <doublevar> & b,
                     const Array1 <doublevar> & origin, 
                     int jj, int kk) {

  double aproj=projection(pos, a);
  double bproj=projection(pos, b);
  
  Array1 <doublevar> corner(3,0.0);
  Array1 <doublevar> dis(3,0.0);
  for(int i=0; i< 3; i++) {
    corner[i]+=(jj+1)*a[i]/2.0;
    corner[i]+=(kk+1)*b[i]/2.0;
  }
  
  double origaproj=projection(origin, a);
  double origbproj=projection(origin, b);
  double lna=length_vec(a);
  double lnb=length_vec(b);
  
  for(int i=0; i< 3; i++) {
    corner[i]+=origaproj*a[i]/lna;
    corner[i]+=origbproj*b[i]/lnb;
  }
  
  for(int i=0; i< 3; i++) {
    dis[i]+=a[i]*aproj/lna;
    dis[i]+=b[i]*bproj/lnb;
    dis[i]-=corner[i];
  }
  return length_vec(dis);
}
//--------------------------------------------------------------------


#include "qmc_io.h"

double find_centers(const Array1 <doublevar> & origin_,
                    const Array2 <doublevar> & latvec_,
                    const Array2 <doublevar> & atompos,
                    Array2 <doublevar> & centerpos,
                    Array1 <int> & equiv_atom,
                    Array2 <int> & center_displacements,
                    doublevar cutoff_divider) {

  //const double cutoff_divider=1.000001;

  //Move the lattice vector into a vector of
  //vectors..  Makes it easier to translate this routine
  Array1 <Array1 <doublevar> > latvec(3);
  for(int i=0; i< 3; i++) {
    latvec(i).Resize(3);
    for(int j=0; j< 3; j++) {
      latvec(i)(j)=latvec_(i,j);
      //cout << latvec(i)(j) << "  " << latvec_(i,j) << endl;
    }
  }
  
  //Copy the input origin so we can manipulate it.
  Array1 <doublevar> origin(3);
  for(int i=0; i< 3; i++)
    origin(i)=origin_(i);

  // original (real) origin
  Array1 <doublevar> rorigin(3);
  rorigin=origin;

  Array1 <double> cross01(3), cross12(3), cross02(3);

  cross(latvec[0], latvec[1], cross01);
  cross(latvec[1], latvec[2], cross12);
  cross(latvec[0], latvec[2], cross02);

  Array1 <doublevar> height(3);
  height(0)=fabs(projection(latvec[0], cross12));
  height(1)=fabs(projection(latvec[1], cross02));
  height(2)=fabs(projection(latvec[2], cross01));
  double cutoff=min(min(height(0), height(1)), height(2))/cutoff_divider;

  double basis_cutoff=cutoff;

  //cout << "heights " << height(0) << "  " << height(1)
  //     << "  " << height(2) << endl;
  single_write(cout, "cutoff ", cutoff, "\n");

  // origin and "lattice vectors" cutoff_piped define a large cell whose
  // sides are cutoff apart from the sides of the simulation cell
  Array1 <Array1 <double > > cutoff_piped(3);
  for(int i=0; i< 3; i++) {
    cutoff_piped(i).Resize(3);
    for(int j=0; j< 3; j++) {
      origin[j]-=cutoff*latvec[i][j]/height(i);
      cutoff_piped(i)(j)=latvec[i][j]*(1+2*cutoff/height(i));
    }
  }

  Array2 < doublevar > cutoff_piped2d(3,3);
  for(int i=0; i< 3; i++) 
    for(int j=0; j< 3; j++)
      cutoff_piped2d(i,j)=cutoff_piped(i)(j);

  // center of mass of the simulation cell
  Array1 < doublevar > cm(3);
  cm=rorigin;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      cm(i)+=0.5*latvec[j][i];
    }
  }
  // radius of the sphere that contains ghost atoms is half of the longest
  // body diagonal plus the cutoff
  Array1 <doublevar> diag1(3), diag2(3);
  diag1=0.0;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      diag1[i]+=latvec[j][i];
    }
  }
  doublevar radius=sqrt(dot(diag1,diag1));
  for (int j=0; j<3; j++) {
    diag2=diag1;
    for (int i=0; i<3; i++) {
      diag2[i]-=2*latvec[j][i];
    }
    doublevar rad2=sqrt(dot(diag2,diag2));
    if ( rad2>radius ) radius=rad2;
  }
  radius=0.5*radius+cutoff;

  // number of adjacent cells to search
  // JK: was 1 but that is insufficient for long and thin simulation cells
  // (observed for MnO's antiferromagnetic primitive cell)
  const int nsearch=4;
  int natoms=atompos.GetDim(0);

  Array1 <doublevar> pos(3), temp(3);
  int total=0;
  //int ncorner=0;
  //int nedge=0;
  //int nside=0;

  Pbc_enforcer pbc;
  pbc.init(cutoff_piped2d);
  pbc.setOrigin(origin);

  vector < vector < doublevar> > tmpcenterpos;
  vector <int> tmpequiv_atom;
  vector < vector < int> > cell_displacement;
  
  for(int ii=-nsearch; ii < nsearch+1; ii++) {
    for(int jj=-nsearch; jj < nsearch+1; jj++) {
      for(int kk=-nsearch; kk < nsearch+1; kk++) {
        for(int at=0; at < natoms; at++) {
          for(int i=0; i< 3; i++) {
            pos[i]=atompos(at, i)
                   +(latvec[0][i])*ii
                   +(latvec[1][i])*jj
                   +(latvec[2][i])*kk;
          }
       
	  // JK: this branching of distance evaluation according to corners,
	  // edges and sides works only for nsearch=1 and strictly speaking
	  // only for orthogonal simulation cells. I will comment it out.
	  // The possible extra centers that will be selected by testing only
	  // for the cutoff_piped and the sphere around center of mass should
	  // not matter that much (speed-wise) in production runs utililizing
	  // MO_Cutoff orbitals.
	  /*
	  int use=0;
	  if(pbc.isInside(pos) ) {
            //if all three are non-zero, we're in a corner
            if(ii*jj*kk !=0) {
              Array1 <double> corner(3);
              corner=0;
              //find the closest corner
              for(int i=0; i< 3; i++) {
                corner[i]+=(ii+1)*latvec[0][i]/2;
                corner[i]+=(jj+1)*latvec[1][i]/2;
                corner[i]+=(kk+1)*latvec[2][i]/2;
                corner[i]+=rorigin[i];
              }
              if( distance_vec(pos, corner) < cutoff ) {
                use=1;
                ncorner++;
              }
            }
            //if we've moved two indices at the same time,
            //we're near the edge of the cell
            else if(ii*jj != 0) { //a and b move
              double dis=distance_edge(pos, latvec[0], latvec[1],
                                       rorigin, ii,jj);
              if(dis < cutoff) {
                use=1;
                nedge++;
              }
              else {
                use=0;
              }
            }
            else if(jj*kk != 0) { //b and c mov
              double dis=distance_edge(pos, latvec[1], latvec[2],
                                       rorigin, jj,kk);
              if(dis < cutoff) {
                use=1;
                nedge++;
              }
              else {
                use=0;
              }
            }
            else if(ii*kk != 0) { //a and c move
              double dis=distance_edge(pos, latvec[0], latvec[2],
                                       rorigin, ii,kk);
              if(dis < cutoff) {
                use=1;
                nedge++;
              }
              else {
                use=0;
              }
            }
            else {
              nside++;
              use=1;
            }
          }
	  */
	  if ( pbc.isInside(pos) 
	       && ( distance_vec(pos,cm) < radius ) ) {
            vector <doublevar> tpos;
            for(int i=0; i< 3; i++) tpos.push_back(pos(i));
            tmpcenterpos.push_back(tpos);
            tmpequiv_atom.push_back(at);
            vector <int> displ(3);
            displ[0]=ii; displ[1]=jj; displ[2]=kk;
            cell_displacement.push_back(displ);         
            total++;
          }
        }
      }
    }
  }

  assert(total==tmpcenterpos.size());

  centerpos.Resize(total,3);
  equiv_atom.Resize(total);
  center_displacements.Resize(total, 3);
  for(int i=0; i< total; i++) {
    for(int j=0; j< 3; j++) {
      centerpos(i,j)=tmpcenterpos[i][j];
      center_displacements(i,j)=cell_displacement[i][j];
      //cout << "pos " << centerpos(i,j) << " displacement "
      //   << cell_displacement[i][j] << endl;
    }
    equiv_atom(i)=tmpequiv_atom[i];
  }
  //debug_write(cout, "nedge ", nedge, " nside ", nside);
  //debug_write(cout, "  ncorner ", ncorner, "\n");
  single_write(cout, " total centers ", total, "\n");

  return basis_cutoff;
}


//----------------------------------------------------------------------
