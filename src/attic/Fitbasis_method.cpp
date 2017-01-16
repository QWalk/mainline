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
//------------------------------------------------------------------------
//src/Fitbasis_method.cpp
#include "Fitbasis_method.h"
#include "qmc_io.h"
#include "Program_options.h"
#include "MatrixAlgebra.h"

void Fitbasis_method::read(vector <string> words,
                           unsigned int & pos,
                           Program_options & options)
{

  vector <vector < string> > basis_sections;

  vector <string> basistmp;
  pos=0;
  while(readsection(words, pos,basistmp, "BASIS")) {
    basis_sections.push_back(basistmp);
  }

  int nbasis=basis_sections.size();

  if(nbasis<= 0) error("Need at least one BASIS section in FITBASIS.");

  basis.Resize(nbasis);
  basis=NULL;
  for(int i=0; i< nbasis; i++) {
    allocate(basis_sections[i], basis(i));
  }


  allocate(options.systemtext[0], sys);


  vector <string> centertext;
  if(readsection(words, pos=0, centertext, "CENTERS"))
  {}
  else
  {
    single_write(cout, "Defaulting to using the atoms as centers\n");
    string temp="USEATOMS";
    centertext.push_back(temp);
  }
  unsigned int newpos=0;
  centers.read(centertext, newpos, sys);
  centers.assignBasis(basis);



  totbasis=0;
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

  cout << "Fitting " << totbasis << " total basis functions." << endl;


  if(!readsection(words, pos=0, file_list, "MO_MESHES") ) {
    error("need list of files in MO_MESHES");
  }


  readvalue(words, pos=0, localize_output, "JEEP_CENTERS");

  if(!readvalue(words, pos=0, orb_out_file, "ORBFILE") ) {
    orb_out_file=options.runid+".orb";
  }

}

//------------------------------------------------------------------------

int Fitbasis_method::showinfo(ostream & os)
{
  os << "###################Fit basis#########################\n";
  return 1;
}

//------------------------------------------------------------------------



//------------------------------------------------------------------------

#include "jeep_utils.h"
#include "Basis_fitter.h"
void Fitbasis_method::run(Program_options & options, ostream & output)
{


  Array1 <doublevar> origin(3);
  Array1 <doublevar> box_size(3);
  Array1 <int> npoints(3);

  if(file_list.size() < 1) error("no files to open");
  ifstream orbin(file_list[0].c_str());
  if(!orbin) error("couldn't open ",file_list[0]);

  int nfunctions;
  get_global_header(orbin, nfunctions);
  get_function_header(orbin, npoints, box_size, origin);
  orbin.close();



  Sample_point * sample=NULL;
  sys->generateSample(sample);
  Array2 <doublevar> basis_matrix;
  calculate_basis_overlap(centers,basis, sample,
                          origin,box_size, npoints,
                          basis_matrix);




  //Note that this is a symmetric matrix, so we could use a better
  //inverter.
  cout << "inverting matrix " << endl;
  Array2 <doublevar> basis_matrix_inverse(totbasis, totbasis);

  doublevar delta0=box_size(0)/(npoints(0)-1);
  doublevar delta1=box_size(1)/(npoints(1)-1);
  doublevar delta2=box_size(2)/(npoints(2)-1);

  

  for(int b1=0; b1 < totbasis; b1++) {
    for(int b2=0; b2 < totbasis; b2++) {
      basis_matrix(b1,b2) *= delta0*delta1*delta2;
    }
  }



  InvertMatrix(basis_matrix, basis_matrix_inverse, totbasis);


  cout << "----overlap matrix-------" << endl;
  for(int b1=0; b1 < totbasis; b1++) {
    for(int b2=0; b2 < totbasis; b2++) {

      cout << setw(14) << basis_matrix(b1,b2);//*delta0*delta1*delta2;
    }
    cout << endl;
  }

  cout << "----inverse matrix-------" << endl;
  for(int b1=0; b1 < totbasis; b1++) {
    for(int b2=0; b2 < totbasis; b2++) {
      cout << setw(14) << basis_matrix_inverse(b1,b2);
    }
    cout << endl;
  }


  //loop over MO's to fit
  //Array1 <doublevar> projection(totbasis);
  Array2 <doublevar> coefficients;

  fit_molecorb(centers, basis, sample, origin, box_size, npoints,
               basis_matrix_inverse, file_list, coefficients);



  if(mpi_info.node==0) {
    cout << "writing orb file " << endl;
    int nmo_fit=file_list.size();

    //print header to the orb file
    ofstream orbout(orb_out_file.c_str());
    orbout.precision(14);
    int totmo=0;
    for(int mo=0; mo < nmo_fit; mo++) {
      for(int cent=0; cent< centers.size(); cent++)
      {
        int currfunc=0;
        for(int j=0; j< centers.nbasis(cent); j++)
        {
          //cout << "basi " << centers.basis(cent,j) << endl;
          int nfunc=basis(centers.basis(cent,j))->nfunc();
          for(int i=0; i< nfunc; i++) {
            currfunc++;
            totmo++;
            orbout << "     " << mo+1 << "   " << currfunc << "    " << cent+1
                  << "     " << totmo << endl;
          }

        }
      }
    }

    orbout << "COEFFICIENTS " << endl;


   for(int mo=0; mo < nmo_fit; mo++) {
      for(int bas=0; bas< totbasis; bas++) {
      orbout << coefficients(mo, bas) << "   ";
      if(bas%5==4) orbout << endl;
      }
      orbout << endl;
    }


    orbout.close();
  }
  delete sample;
  sample=NULL;


}



//------------------------------------------------------------------------
