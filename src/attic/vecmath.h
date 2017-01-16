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
//utils/vecmath.h

//This is a relatively slow library for vectors using the STL 
//vector <double>'s. 


#ifndef VECMATH_H_INCLUDED
#define VECMATH_H_INCLUDED

#include <vector>
#include <cmath>
/*!
 */
inline std::vector <double> cross(std::vector <double> & a, std::vector <double> & b) {
  assert(a.size()==3);
  assert(b.size()==3);
  std::vector <double> c;
  c.push_back(a[1]*b[2]-a[2]*b[1]);
  c.push_back(a[2]*b[0]-a[0]*b[2]);
  c.push_back(a[0]*b[1]-a[1]*b[0]);
  return c;
}

inline double dot(std::vector <double> & a, std::vector <double> & b) {
  assert(a.size()==3);
  assert(b.size()==3);
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

inline double length_vec(std::vector <double> & a) {
  assert(a.size()==3);
  return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

inline std::vector <double> operator+(std::vector <double> & a, std::vector <double> & b) {
  assert(a.size()==3);
  assert(b.size()==3);
  std::vector <double> c;
  for(int i=0; i< 3; i++) {
    c.push_back(a[i]+b[i]);
  }
  return c;
}

inline double distance_vec(std::vector <double> & a, std::vector <double> & b) {
  assert(a.size()==3);
  assert(b.size()==3);
  return sqrt(
	      (a[0]-b[0])*(a[0]-b[0])
	      +(a[1]-b[1])*(a[1]-b[1])
	      +(a[2]-b[2])*(a[2]-b[2]));
}

inline std::vector <double> operator-(std::vector <double> & a, std::vector <double> & b) {
  assert(a.size()==3);
  assert(b.size()==3);
  std::vector <double> c;
  for(int i=0; i< 3; i++) {
    c.push_back(a[i]-b[i]);
  }
  return c;
}

inline std::vector <double> operator*(std::vector <double> & v, int & i) {
  assert(v.size()==3);
  std::vector <double> c;
  for(int i=0; i< 3; i++) c.push_back(i*v[i]);
  return c;
}


//projection of a onto b (returns a dot b/abs(b))
inline double projection(std::vector <double> & a, std::vector <double> & b) {
  assert(a.size()==3);
  assert(b.size()==3);
  double c;
  c=dot(a,b);
  double bsize=dot(b,b);
  return c/sqrt(bsize);
}



class Shifter {

public:
  Shifter() {
    for(int i=0; i< 3; i++) origin.push_back(0);
  }
  std::vector <double> origin;
  int enforcepbc(std::vector <double> & x, std::vector <std::vector <double> > & latvec ) {

    //Find normal vectors..
    std::vector <std::vector <double> > norm;
    std::vector <double> temp;
    for(int i=0; i<3; i++) temp.push_back(0);
    for(int i=0; i<3; i++) norm.push_back(temp);

    norm[0]=cross(latvec[1], latvec[2]);
    norm[1]=cross(latvec[0], latvec[2]);
    norm[2]=cross(latvec[0], latvec[1]);

    /*
    for(int i=0; i< 3; i++) {
      cout << "normal before " << i << "   ";
      for(int j=0; j< 3; j++) {
	cout << norm[i][j] << "  ";
      }
      cout << endl;
    }
    */

    //Make sure the normal vectors are pointing out of the box
    for(int i=0; i< 3; i++) {
      //cout << "dot " << i << dot(norm[i], latvec[i]) << endl;
      if(dot(norm[i], latvec[i]) < 0) {
        for(int j=0; j< 3; j++) {
          norm[i][j]=-norm[i][j];
        }
      }
    }

    /*
    for(int i=0; i< 3; i++) {
      cout << "normal " << i << "   ";
      for(int j=0; j< 3; j++) {
	cout << norm[i][j] << "  ";
      }
      cout << endl;
    }
    */

    //find corners
    std::vector <std::vector <double> > corners;
    for(int i=0; i <3; i++) corners.push_back(temp);

    for(int i=0; i< 3; i++) {
      corners[i]=origin+latvec[i];
    }

    int shifted=0;
    //Initialization done, now do the pbc..
    for(int i=0; i< 3; i++) { //loop over lattice vectors
      int shouldcontinue=1;
      while(shouldcontinue) {
        temp=x-origin;
        double sm=dot(norm[i], temp);
        temp=x-corners[i];
        double sp=dot(norm[i], temp);
	      //	cout << "sm " << sm << "  sp " << sp << endl;

        if(sm < 0) {
          x=x+latvec[i];
          shifted=1;
        }
        else if(sp >0 ) {
          x=x-latvec[i];
          shifted=1;
        }
        else {
          shouldcontinue=0;
        }
      }
    }
    return shifted;
  }

};

//----------------------------------------------------------------------

//check whether the vector is enclosed in the parallelpiped defined by
//the three vectors in box(inclusive edges)
inline bool is_enclosed(std::vector <double> & x,
		 std::vector <std::vector <double> > & box) {

     /*
  //This is quite possibly wrong when the box edges are not orthogonal. 
  //Using the shifter instead
  assert(box.size()==3);
  assert(box[0].size()==3);
  assert(box[1].size()==3);
  assert(box[2].size()==3);
  assert(x.size()==3);
  std::vector <double > lengths;
  for(int i=0; i< 3; i++) {
    lengths.push_back(sqrt(dot(box[i], box[i])));
  }

  bool in_box=true;
  for(int i=0; i< 3; i++) {
    double proj=projection(x,box[i]);
    //cout << "projection " << i << "   " << proj
    // << "  length " << lengths[i]
    // << endl;
    if(proj > lengths[i] || proj < 0 ) {
      in_box=false;
    }
  }

  //return in_box;
  */
  Shifter shft;
  return !shft.enforcepbc(x,box);
}



#endif //VECMATH_H_INCLUDED

//----------------------------------------------------------------------
