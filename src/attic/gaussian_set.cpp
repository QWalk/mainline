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
#include "gaussian_set.h"
#include "qmc_io.h"

/*!
Uses flat file of form:
ncenters
Label x y z
Label x y z
...
 */
void read_centerpos(string & filename, Array2 <doublevar> & position, vector <string> & labels)
{
  int ncenters;
  ifstream centin(filename.c_str());
  if(!centin)
    error("Couldn't open ", filename);
  centin >> ncenters;
  labels.clear();
  position.Resize(ncenters,3);
  string labeltemp;
  for(int i=0; i< ncenters; i++)
  {
    centin >> labeltemp;
    labels.push_back(labeltemp);
    for(int d=0; d< 3; d++)
    {
      centin >> position(i,d);
    }
  }
  centin.close();
}

//----------------------------------------------------------------------


void read_basis(string & basisin, Array1 <Contracted_gaussian> & basis) {
  ifstream basis_f(basisin.c_str());
  if(!basis_f)
    error("Couldn't open ", basisin);
  vector <string> words;
  parsefile(basis_f, words);
  unsigned int pos=0;
  vector <vector  <string> > sections;
  vector <string> section_temp;
  while(readsection(words, pos, section_temp, "BASIS") ) {
    sections.push_back(section_temp);
  }

  int nsections=sections.size();

  string label;
  vector <string> gamesstxt;
  vector < Contracted_gaussian > basis_tmp;


  Contracted_gaussian tempgauss;
  int totgauss=0;

  for(int sec=0; sec < nsections; sec++) {
    label=sections[sec][0];
    //cout << "label " << label << endl;
    if(!readsection(sections[sec], pos=0, gamesstxt, "GAMESS") ) {
      error("Need a GAMESS section");
    }
    unsigned int pos2=0;
    unsigned int txtsize=gamesstxt.size();
    string symm;
    int nexpand=0;

    while(pos2 < txtsize) {
      symm=gamesstxt[pos2];
      //cout << "symm " << symm << endl;
      
      nexpand=atoi(gamesstxt[++pos2].c_str());
      pos2+= 3*nexpand;
      pos2++;
      if(symm=="S") totgauss++;
      else if(symm=="P") totgauss+=3;
      else if(symm=="D" || symm=="6D") totgauss+=6;
      else error("Don't understand symmetry ", symm);
      //cout << "totgauss " << totgauss << endl;
    }
  }

  basis.Resize(totgauss);


  int gaussnum=0;
  for(int sec=0; sec < nsections; sec++) {
    label=sections[sec][0];
    //cout << "label " << label << endl;
    if(!readsection(sections[sec], pos=0, gamesstxt, "GAMESS") ) {
      error("Need a GAMESS section");
    }
    unsigned int pos2=0;
    unsigned int txtsize=gamesstxt.size();
    string symm;
    int nexpand=0;



    while(pos2 < txtsize) {
      symm=gamesstxt[pos2];
      nexpand=atoi(gamesstxt[++pos2].c_str());

      tempgauss.alpha.Resize(nexpand);
      tempgauss.coeff.Resize(nexpand);
      tempgauss.center_name=label;

      for(int i=0; i< nexpand; i++) {
        if(atoi(gamesstxt[++pos2].c_str()) != i+1) error("Some misformatting in basis section");
        tempgauss.alpha(i)=atof(gamesstxt[++pos2].c_str());
        tempgauss.coeff(i)=atof(gamesstxt[++pos2].c_str());
      }


      if(symm=="S") {
        tempgauss.lvals(0)=tempgauss.lvals(1)=tempgauss.lvals(2)=0;
        basis(gaussnum++)=tempgauss;
      }
      else if(symm=="P")  {
        for(int i=0; i< 3; i++) {
          tempgauss.lvals=0;
          tempgauss.lvals(i)=1;
          basis(gaussnum++) = tempgauss;
        }
      }
      else if(symm=="D" || symm=="6D") {

        //xx, yy, zz
        for(int i=0; i< 3; i++) {
          tempgauss.lvals=0;
          tempgauss.lvals(i)=2;
          basis(gaussnum++) = tempgauss;
        }
        //xy
        tempgauss.lvals(0)=1; tempgauss.lvals(1)=1; tempgauss.lvals(2)=0;
        basis(gaussnum++) = tempgauss;

        //xz
        tempgauss.lvals(0)=1; tempgauss.lvals(1)=0; tempgauss.lvals(2)=1;
        basis(gaussnum++) = tempgauss;

        //yz
        tempgauss.lvals(0)=0; tempgauss.lvals(1)=1; tempgauss.lvals(2)=1;
        basis(gaussnum++) = tempgauss;
      }
      else {
        error("Unsupported symmetry ", symm);
      }
      ++pos2;
    }


  }




/*
  for(int i=0; i< totgauss; i++) {
    cout << "basis " << i << "  center " << basis(i).center_name
	 << endl;
    cout << "  lvals "; for(int j=0; j< 3; j++) cout << basis(i).lvals(j) << "  ";
    cout << endl;
    for(int j=0; j< basis(i).alpha.GetDim(0); j++) {
      cout << "   " << basis(i).alpha(j) << "   " << basis(i).coeff(j) << endl;
    }

  }
*/

}

//----------------------------------------------------------------------

void create_local_basis(string & centerin, string & basisin,
			Array1 <Center> & centers,
                        Array1 <Contracted_gaussian >  & basis ) {

  Array2 <doublevar> centerpos;
  vector <string> labels;


  read_centerpos(centerin, centerpos, labels);

  read_basis(basisin, basis);


  int ncenters=centerpos.GetDim(0);
  int nbasis=basis.GetDim(0);

  centers.Resize(ncenters);
  for(int cen=0; cen < ncenters; cen++) {
    for(int d=0; d< 3; d++) centers(cen).pos(d)=centerpos(cen,d);
    centers(cen).label=labels[cen];
  }



  vector <vector <int> > cenbasis, basiscen;
  cenbasis.resize(ncenters);
  basiscen.resize(nbasis);
  for(int cen=0; cen < ncenters; cen++) {
    int foundbas=0;
    for(int bas=0; bas < nbasis; bas++) {

      if(centers(cen).label==basis(bas).center_name) {
        foundbas=1;
        cenbasis[cen].push_back(bas);
        basiscen[bas].push_back(cen);
      }
    }
    if(!foundbas) 
      cout << "WARNING!  Couldn't find a basis set for center "
           << centers(cen).label << endl;
  }

  

  for(int bas=0; bas < nbasis; bas++) {
    int ncen=basiscen[bas].size();
    basis(bas).center.Resize(ncen);
    for(int c=0; c< ncen; c++) {
      basis(bas).center(c)=basiscen[bas][c];
      //cout << "basis " << bas << " ->  center " << basis(bas).center(c) << endl;
    }
  }

  for(int cen=0; cen < ncenters; cen++) {
    int nbas=cenbasis[cen].size();
    centers(cen).basis.Resize(nbas);
    for(int b=0; b< nbas; b++) {
      centers(cen).basis(b)=cenbasis[cen][b];
      //cout << "center " << cen << " -> basis " << centers(cen).basis(b) << endl;
    }
  }



  for(int bas=0; bas < nbasis; bas++) {
    int nalpha=basis(bas).alpha.GetDim(0);
    //cout << "basis " << bas << endl;
    for(int a =0; a < nalpha; a++) {
        int totL=basis(bas).lvals(0)+basis(bas).lvals(1)+basis(bas).lvals(2);
        doublevar exponent=basis(bas).alpha(a);
        doublevar fac=sqrt(2.*exponent/pi);
        doublevar feg=4.*exponent;
        doublevar feg2=feg*feg;
        doublevar norm=0;
        switch(totL)
        {
        case 0:
          norm=sqrt(2.*feg*fac);
          break;
        case 1:
          norm=sqrt(2.*feg2*fac/3.);
          break;
        case 2:
          norm=sqrt(2.*feg*feg2*fac/15.);
          break;
        case 3:
          norm=sqrt(2.*feg2*feg2*fac/105.);
          break;
        default:
          norm=0;
          error("Unknown symmetry in Cubic_spline::readbasis! Shouldn't be here!");
        }
        basis(bas).coeff(a)*=norm;
    }
  }


}
//------------------------------------------------------------------------

#include "Pbc_enforcer.h"
void Gaussian_lookups::set_lookup_tables(Array2 <doublevar> & latvec,
                                         Array1 <doublevar> & origin,
                                         Array1 <Center> & centers) {

  int ncenters=centers.GetDim(0);
  // int nbasis=basis.GetDim(0);

  //Find the centers that are within the cell
  Pbc_enforcer pbc;
  pbc.init(latvec);
  pbc.setOrigin(origin);

  Array1 <int> center_in_cell(ncenters);
  int ncenters_inequivalent=0;
  Array1 <int> inequiv_centers(ncenters);
  for(int cen=0; cen< ncenters; cen++) {
    center_in_cell(cen)=pbc.isInside(centers(cen).pos);
    if(center_in_cell(cen)==1) {
      inequiv_centers(ncenters_inequivalent++) = cen;
    }
  }


  //cout << "There are " << ncenters << " total centers, of which "
  //     << ncenters_inequivalent << " are in the simulation cell." << endl;

  //For the centers that aren't in the cell, find their equivalent
  //one within.
  equivalent_center.Resize(ncenters);
  Array1 <doublevar> temp_pos(3);
  for(int cen=0; cen < ncenters; cen++) {
    if(!center_in_cell(cen)) {
      temp_pos=centers(cen).pos;

      pbc.enforcePbc(temp_pos);

      //cout << "reduced position ";
      //for(int d=0; d < 3; d++) cout << temp_pos(d) << "  ";
      //cout << endl;
      int found_ecen=0;
      for(int j=0; j< ncenters_inequivalent; j++) {
        int cen2=inequiv_centers(j);
        doublevar test_sum=0;
        for(int d=0; d< 3; d++)
          test_sum+=fabs(temp_pos(d)-centers(cen2).pos(d));
        //cout << cen2 << " test sum " << test_sum << endl;
        const doublevar threshold=1e-3;
        if(test_sum < threshold) {
          equivalent_center(cen)=cen2;
          found_ecen=1;
          //cout << "center " << cen << "  equiv to " << cen2 << endl;
          break;
        }
      }
      if(!found_ecen) error("Couldn't find equivalent center for center ", cen);
    }
    else {
      equivalent_center(cen)=cen;
      //cout << "center " << cen << " in the box " << endl;
    }
  }

  Array1 <int> cenbasis(ncenters);
  int totbasis=0;
  for(int c=0; c < ncenters_inequivalent; c++) {
    int cen=inequiv_centers(c);
    cenbasis(cen)=centers(cen).basis.GetDim(0);
    totbasis+=cenbasis(cen);
  }
  //cout << totbasis << " total basis functions " << endl;



  //Lookup tables to get the center and basis for
  //the absolute basis function order
  //Array1 <int> totbasis2cen(totbasis);
  //Array1 <int> totbasis2bas(totbasis);
  //Array1 <int> center_start(ncenters);
  //Array1 <int> center_end(ncenters);
  totbasis2cen.Resize(totbasis);
  totbasis2bas.Resize(totbasis);
  center_start.Resize(ncenters);
  center_end.Resize(ncenters);
  center_start=0; center_end=0;
  int counter=0;
  for(int c=0; c < ncenters_inequivalent; c++) {
    int cen=inequiv_centers(c);
    center_start(cen)=counter;
    for(int bas=0; bas < centers(cen).basis.GetDim(0); bas++) {
      totbasis2bas(counter)=centers(cen).basis(bas);
      totbasis2cen(counter)=cen;
      counter++;
    }
    center_end(cen)=counter;
    //cout << "center " << cen << " start " << center_start(cen)
    //     << "   end  " << center_end(cen)
    //     << endl;
  }

}

