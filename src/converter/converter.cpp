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
#include "converter.h"
#include <iomanip>
#include <cstdlib>
#include <algorithm>

using namespace std;

void split
(string & text, string & separators, vector<string> & words)
{
  int n = text.length();
  int start, stop;

  start = text.find_first_not_of(separators);
  while ((start >= 0) && (start < n)) {
    stop = text.find_first_of(separators, start);
    if ((stop < 0) || (stop > n)) stop = n;
    words.push_back(text.substr(start, stop - start));
    start = text.find_first_not_of(separators, stop+1);
  }
}

void append_number(string & str, int num)
{
  char strbuff[40];
  sprintf(strbuff, "%d", num);
  str+=strbuff;
}

//--------------------------------------------------------



void print_orbitals(ostream & os,
                    std::vector <Center> & centers,
                    std::vector <int> & nbasis,
                    std::vector < std::vector < double > > & moCoeff) {

  cout << "writing orb file.  This could take a few seconds." << endl;
  os.precision(12);
  const double threshold=1e-7;
  int ncenters=centers.size();
  int totcoeff=0;
  int totmo=moCoeff.size();

  vector <int> totfunc_start;
  int natoms=nbasis.size();
  int totbasis=0;

  for(int at=0; at< natoms; at++) {
    totfunc_start.push_back(totbasis);
    totbasis+=nbasis[at];
  }

  for(int mo=0; mo < totmo; mo++) {
    for(int cent=0; cent < ncenters; cent++) {
      int at=centers[cent].equiv_atom;
      for(int ao=0; ao < nbasis[at]; ao++) {
      //if(fabs(moCoeff[mo][totfunc_start[at]+ao]) > threshold) {
        os << setw(10) << mo+1
           << setw(10) << ao+1
           << setw(10) << cent+1
           << setw(10) << totcoeff+1 << endl;
        totcoeff++;
  //  }
      }
    }
  }

  totcoeff=0;
  os << "COEFFICIENTS\n";
  for(int mo=0; mo < totmo; mo++) {
    for(int cent=0; cent < ncenters; cent++) {
      int at=centers[cent].equiv_atom;
      for(int ao=0; ao < nbasis[at]; ao++) {
     //if(fabs(moCoeff[mo][totfunc_start[at]+ao]) > threshold) {
        os << moCoeff[mo][totfunc_start[at]+ao] << "  ";
        totcoeff++;
        if(totcoeff%5==0) os << endl;
  //  }
      }
    }
  }

}


void print_orbitals(ostream & os,
                    std::vector <Center> & centers,
                    std::vector <int> & nbasis,
                    std::vector < std::vector < dcomplex > > & moCoeff) {

  cout << "writing orb file.  This could take a few seconds." << endl;
  os.precision(12);
  const double threshold=1e-7;
  int ncenters=centers.size();
  int totcoeff=0;
  int totmo=moCoeff.size();

  vector <int> totfunc_start;
  int natoms=nbasis.size();
  int totbasis=0;

  for(int at=0; at< natoms; at++) {
    totfunc_start.push_back(totbasis);
    totbasis+=nbasis[at];
  }

  for(int mo=0; mo < totmo; mo++) {
    for(int cent=0; cent < ncenters; cent++) {
      int at=centers[cent].equiv_atom;
      for(int ao=0; ao < nbasis[at]; ao++) {
      //if(fabs(moCoeff[mo][totfunc_start[at]+ao]) > threshold) {
        os << setw(10) << mo+1
           << setw(10) << ao+1
           << setw(10) << cent+1
           << setw(10) << totcoeff+1 << endl;
        totcoeff++;
  //  }
      }
    }
  }

  totcoeff=0;
  os << "COEFFICIENTS\n";
  for(int mo=0; mo < totmo; mo++) {
    for(int cent=0; cent < ncenters; cent++) {
      int at=centers[cent].equiv_atom;
      for(int ao=0; ao < nbasis[at]; ao++) {
     //if(fabs(moCoeff[mo][totfunc_start[at]+ao]) > threshold) {
        os << moCoeff[mo][totfunc_start[at]+ao] << "  ";
        totcoeff++;
        if(totcoeff%5==0) os << endl;
  //  }
      }
    }
  }

}


//###########################################################################
void find_unique_atoms(const vector<Atom> & atoms, vector<string> & unique_atoms) { 
  unique_atoms.clear();
  for(vector<Atom>::const_iterator at=atoms.begin(); at != atoms.end();
      at++) { 
    if(find(unique_atoms.begin(), unique_atoms.end(),at->name)==unique_atoms.end()){
      unique_atoms.push_back(at->name);
    }
  }
}

//######################################################################


void print_centers(ostream & os, vector <Center> & centers) {
  int ncenters=centers.size();
  os.precision(12);
  os << ncenters << endl;
  for(int cent=0; cent < ncenters; cent++) {
    os << centers[cent].name
       << "  ";
    for(int i=0; i< 3; i++) os << centers[cent].pos[i] << "   ";
    os << endl;
  }

}


//######################################################################
void print_vmc_section( ostream & os, string & outputname, double eref) {
  os << "METHOD {\n"
     << "  VMC\n"
     << "  NBLOCK   10\n"
     << "  NSTEP    10\n"
     << "  NDECORR   4\n"
     << "  TIMESTEP 1.0\n"
     << "  NCONFIG  100\n"
     << "  STORECONFIG  " << outputname << ".config\n\n"
     << "  #uncomment the following to read a configuration\n"
     << "  #READCONFIG  " << outputname << ".config\n"
    //    << "  EREF " << eref << endl
     << "}\n\n";

}


void print_opt_section( ostream & os, string & outputname, double eref) {
  os << "METHOD {\n"
     << "  OPTIMIZE2\n\n"
     << "  #Number of optimization steps to do.  Should be roughly \n"
        "  #30 times the number of variational parameters\n"
     << "  ITERATIONS 30\n"
     << "  READCONFIG " << outputname << ".config\n"
     << "  NCONFIG  1000\n"
    //     << "  EREF " << eref << endl
     << "  MINFUNCTION MIXED \n"
    //<< "  PSEUDOTEMP " << outputname << ".pseudo\n"
     << "}\n\n";

}

void print_dmc_section( ostream & os, string & outputname, double eref) {
  os << "METHOD {\n"
     << "  DMC\n"
     << "  NBLOCK   10\n"
     << "  NSTEP    100\n"
     << "  NDECORR   4\n"
     << "  TIMESTEP 0.01\n"
     << "  NCONFIG  100\n"
     << "  STORECONFIG  " << outputname << ".config\n\n"
     << "  READCONFIG  " << outputname << ".config\n"
    //  << "  EREF " << eref << endl
     << "}\n\n";

}




//######################################################################
#include "vecmath.h"

double find_basis_cutoff(std::vector <std::vector <double> > & latvec) {
  double cutoff_divider=2.000001;
  assert(latvec.size()==3);
  assert(latvec[0].size()==3);
  assert(latvec[1].size()==3);
  assert(latvec[2].size()==3);


  //vector <double> origin;
  vector <double> cross01, cross12, cross02;
  for(int i=0; i< 3; i++) {
    //origin.push_back(0);
    cross01.push_back(0);
    cross12.push_back(0);
    cross02.push_back(0);
  }


  //TODO-Check the cutoff determination..
  cross01=cross(latvec[0], latvec[1]);
  cross12=cross(latvec[1], latvec[2]);
  cross02=cross(latvec[0], latvec[2]);

  double height0=fabs(projection(latvec[0], cross12));
  double height1=fabs(projection(latvec[1], cross02));
  double height2=fabs(projection(latvec[2], cross01));
  return min(min(height0, height1), height2)/cutoff_divider;
}

//--------------------------------------------------------------------

double find_centers(vector <double> & origin,
                    vector <vector <double> > & latvec,
                    vector <Atom> & atoms,
                    vector <Center> & centers ) {
  const double cutoff_divider=1.000001;

  assert(origin.size() == 3);
  assert(latvec.size()==3);
  assert(latvec[0].size()==3);
  assert(latvec[1].size()==3);
  assert(latvec[2].size()==3);


  //vector <double> origin;
  vector <double> cross01, cross12, cross02;
  for(int i=0; i< 3; i++) {
    //origin.push_back(0);
    cross01.push_back(0);
    cross12.push_back(0);
    cross02.push_back(0);
  }


  //TODO-Check the cutoff determination..
  cross01=cross(latvec[0], latvec[1]);
  cross12=cross(latvec[1], latvec[2]);
  cross02=cross(latvec[0], latvec[2]);

  double height0=fabs(projection(latvec[0], cross12));
  double height1=fabs(projection(latvec[1], cross02));
  double height2=fabs(projection(latvec[2], cross01));
  double cutoff=min(min(height0, height1), height2)/cutoff_divider;

  double basis_cutoff=cutoff;

  //cout << "heights " << height0 << "  " << height1 << "  " << height2 << endl;
  //cout << "cutoff " << cutoff << endl;

  vector <vector <double > > cutoff_piped;
  for(int i=0; i< 3; i++) {
    vector <double> tmp;
    double norm=sqrt(dot(latvec[i], latvec[i]));
    for(int j=0; j< 3; j++) {
      origin[j]-=cutoff*latvec[i][j]/norm;
      tmp.push_back(latvec[i][j]*(1+2*cutoff/norm));
    }
    cutoff_piped.push_back(tmp);
  }

  //the origin and cutoff_piped now define the parallelpiped that
  //we want to find centers in.

  //cout << "origin";
  //for(int i=0; i< 3; i++) cout << origin[i] << "  ";
  //cout << endl;

  //cout << "box " << endl;
  //for(int i=0; i< 3; i++) {
   // for(int j=0; j< 3; j++) {
   //   cout << cutoff_piped[i][j] << "  ";
   // }
  //  cout << endl;
  //}


  const int nsearch=1;//number of adjacent cells to search
  int natoms=atoms.size();
  vector <double > pos;
  for(int i=0; i< 3; i++) pos.push_back(0);
  vector <double> effpos; //position shifted by origin
  for(int i=0; i< 3; i++) effpos.push_back(0);
  vector <double> temp; //temporary vector
  for(int i=0; i< 3; i++) temp.push_back(0);
  Center tmpcenter;

  int total=0;

  int ncorner=0;
  int nedge=0;
  int nside=0;


  for(int ii=-nsearch; ii < nsearch+1; ii++) {
    for(int jj=-nsearch; jj < nsearch+1; jj++) {
      for(int kk=-nsearch; kk < nsearch+1; kk++) {
        for(int at=0; at < natoms; at++) {
          for(int i=0; i< 3; i++) {
            pos[i]=atoms[at].pos[i]
                   +(latvec[0][i])*ii
                   +(latvec[1][i])*jj
                   +(latvec[2][i])*kk;
          }
          effpos=pos-origin;
          int use=0;
          if(is_enclosed(effpos, cutoff_piped) ) {

            //if all three are non-zero, we're in a corner
            if(ii*jj*kk !=0) {
              vector <double> corner;
              for(int i=0; i< 3; i++) corner.push_back(0);

              //find the closest corner
              for(int i=0; i< 3; i++) {
                corner[i]+=(ii+1)*latvec[0][i]/2;
                corner[i]+=(jj+1)*latvec[1][i]/2;
                corner[i]+=(kk+1)*latvec[2][i]/2;
              }
              if( distance_vec(pos, corner) < cutoff ) {
                use=1;
                ncorner++;
              }
            }

            //if we've moved two indices at the same time,
            //we're near the edge of the cell
            else if(ii*jj != 0) { //a and b move
              double aproj=projection(pos, latvec[0]);
              double len=length_vec(latvec[0]);
              if(aproj > len) {
                aproj-=len;
              }
              double bproj=projection(pos, latvec[1]);
              len=length_vec(latvec[1]);
              if(bproj > len) {
                bproj -= len;
              }
              if(sqrt(aproj*aproj+bproj*bproj) < cutoff) {
                use=1;
                nedge++;
              }
              else {
                use=0;
              }
            }
            else if(jj*kk != 0) { //b and c move
              double aproj=projection(pos, latvec[1]);
              double len=length_vec(latvec[1]);
              if(aproj > len) {
                aproj-=len;
              }
              double bproj=projection(pos, latvec[2]);
              len=length_vec(latvec[2]);
              if(bproj > len) {
                bproj -= len;
              }
              if(sqrt(aproj*aproj+bproj*bproj) < cutoff) {
                use=1;
                nedge++;
              }
              else {
                use=0;
              }
            }
            else if(ii*kk != 0) { //a and c move
              double aproj=projection(pos, latvec[0]);
              double len=length_vec(latvec[0]);
              if(aproj > len) {
                aproj-=len;
              }
              double bproj=projection(pos, latvec[2]);
              len=length_vec(latvec[2]);
              if(bproj > len) {
                bproj -= len;
              }
              if(sqrt(aproj*aproj+bproj*bproj) < cutoff) {
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

          if(use ) {
            tmpcenter.pos=pos;
            tmpcenter.name=atoms[at].name;
            tmpcenter.equiv_atom=at;
            tmpcenter.basis=atoms[at].basis;
            centers.push_back(tmpcenter);
            total++;
            //cout << atoms[at].name
            //     << "  enclosed " << pos[0] << "  " << pos[1]
            //     << "   " << pos[2] << "  " << at+1 << endl;
          }
        }
      }
    }
  }

  cout << "nedge " << nedge << " nside " << nside << "  ncorner " << ncorner << endl;
  cout << "total of " << total << " centers " << endl;
  return basis_cutoff;
}


//#####################################################################


bool compare_mo(vector <vector < double> > & oldMOCoeff,
                vector <vector < double> > & moCoeff,
                vector <int> & compare_list  ) {


 int nfunctions=oldMOCoeff[0].size();

 int ncompare=compare_list.size();
 if(moCoeff.size() < ncompare) {
   cout << "the punch file doesn't have enough molecular orbitals to compare." << endl;
   exit(1);
 }
 if(oldMOCoeff.size() < ncompare) {
   cout << "the old punch file doesn't have enough MO's to compare." << endl;
 }

  if(moCoeff[0].size() != nfunctions) {
    cout << "Number of functions don't match between new and old mo\n";
    exit(1);
  }
  vector <int> unresolved_mos;

  //First check to see if the mo's are in the same place
  //(most should be)
  for(int i=0; i< ncompare; i++) {
    int mo=compare_list[i];
    double dot=0, mag_old=0, mag_new=0;
    for(int f=0; f< nfunctions; f++) {
      dot+=moCoeff[mo][f]*oldMOCoeff[mo][f];
      mag_old+=oldMOCoeff[mo][f]*oldMOCoeff[mo][f];
      mag_new+=moCoeff[mo][f]*moCoeff[mo][f];
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
        dot+=moCoeff[mo1][f]*oldMOCoeff[mo2][f];
        mag_old+=oldMOCoeff[mo2][f]*oldMOCoeff[mo2][f];
        mag_new+=moCoeff[mo1][f]*moCoeff[mo1][f];
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
