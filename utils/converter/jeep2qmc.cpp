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
#include "wf_writer.h"
#include "basis_writer.h"
#include "Pseudo_writer.h"
#include <fstream>
#include <cstdio>
#include "vecmath.h"
#include "elements.h"
#include "JeepBasis.h"

using namespace std;

void read_jeep_sys(std::istream & is, vector < Atom > & atoms,
                   vector <Spline_pseudo_writer> & jeep_pseudo,
                   vector <double> & origin,
                   vector < vector < double> > & latvec,
                   vector < vector < double> > & ref_latvec);

void read_jeep_header(string & filename, vector <double> & occupation,
                      int & nspin, int & nst,int & nst2, int & ngw);


void read_jeep_wf(string & filename, RealBasis & jeepbasis,
     vector < vector <double> > & moCoeff);

void write_jeep_basis(ostream & os, RealBasis & jeepbasis) {
    int size=jeepbasis.size();
    cout << "size " << size << endl;
    os << "BASIS {\n";
    os << "  PLANEWAVE_ORIGIN\n";
    os << "  PLANEWAVE\n";
    os << "  GVECTOR { \n";
    for(int i=0; i< size; i++) {
      os << "      " << jeepbasis.gx[i] << "   "
      << jeepbasis.gx[size+i] << "   "
      << jeepbasis.gx[size+size+i] << endl;
    }
    os << "  }\n}\n";

}

void usage(const char * name) {
  cout << "usage: " << name <<   " <options> -ecut <cutoff energy> <output> " << endl;
  cout << "Where options can be: \n";
  cout << "-o          Base name for your run case\n";
  exit(1);
}


int main(int argc, char ** argv) {
  double ecut=-1;
  string lcao_basis, infilename, outputname;
  bool use_lcao=false;

  if(argc >= 2) {
    infilename=argv[argc-1];
  }
  else { usage(argv[0]); }
  use_lcao=false;
  for(int i=1; i< argc-1; i++) {
    if(!strcmp(argv[i], "-o") && argc > i+1) {
      outputname=argv[++i];
    }
    else if(!strcmp(argv[i], "-ecut") && argc > i+1) {
      ecut=atof(argv[++i]);
    }
    else if(!strcmp(argv[i], "-lcao") && argc > i+1) {
      use_lcao=true;
      lcao_basis=argv[++i];
    }
    else {
      usage(argv[0]);
    }
  }

  if(ecut==-1) {
    cout << "Must give -ecut <double> " << endl;
    exit(1);
  }



  if(outputname == "") {
    outputname=infilename;
  }


  ///////////////////////////////////////////////////////
  vector <Spline_pseudo_writer> pseudo;
  vector < vector <double> > latvec, ref_latvec;
  vector <double> origin;
  vector <Atom> atoms;
  Slat_wf_writer slwriter;
  vector <double> occupation; //electronic occupation

  RealBasis jeepbasis;
  if(use_lcao) slwriter.mo_matrix_type="CUTOFF_MO";
  else slwriter.mo_matrix_type="STANDARD_MO";
  slwriter.calctype="RHF";  //How do we know what the calculation type is?
  double eref=0;
  string sysfilename=infilename+".sys";
  string wffilename=infilename+".wf";

  ifstream sysfile(sysfilename.c_str());
  if(!sysfile) {
    cout << "couldn't open " << sysfilename << endl;
    exit(1);
  }
  read_jeep_sys(sysfile, atoms, pseudo, origin, latvec, ref_latvec);
  sysfile.close();
  int nspin, nst,nst2, ngw;
  read_jeep_header(wffilename, occupation, nspin, nst,nst2, ngw);

  //We should have the atoms and pseudopotentials..
  int nelectrons=0;
  for(unsigned int at=0; at < atoms.size(); at++) {
    nelectrons+=int(atoms[at].charge);
  }

  cout << "number of mo's " << occupation.size() << endl;
  slwriter.nup=0; slwriter.ndown=0;

  if(nspin==1) {
    for(int i=0; i< nst; i++) {
      if(occupation[i] >= .99) 
        slwriter.nup++;
      if(occupation[i] >= 1.05) 
        slwriter.ndown++;
    }
  }
  else if(nspin==2) {
    slwriter.calctype="UHF";
    slwriter.spin_dwn_start=nst;
    for(int i=0; i< nst; i++) {
      if(occupation[i] >=.99) 
        slwriter.nup++;
    }
    for(int i=nst; i< nst+nst2; i++) {
      if(occupation[i] >=.99) 
        slwriter.ndown++;
    }
  }
      


  //-------------------------------
  //print out the qmc input file

  string orboutname=outputname+".orb";
  slwriter.orbname=orboutname;
  string centeroutname;



  vector <Center> centers;
  vector < vector < double> > moCoeff;
  if(use_lcao) {
    centeroutname=lcao_basis+".centers";
    slwriter.write_centers=false;
    slwriter.use_global_centers=true;
    slwriter.basisname=lcao_basis;
    vector <double> origin_temp;
    for(int i=0; i< 3; i++) origin_temp.push_back(origin[i]);


    Shifter shft;
    for(int i=0; i< 3; i++)
      shft.origin[i]=origin[i];

    int natoms=atoms.size();
    vector<int> nshift;
    for(int a=0; a< natoms; a++) {
      shft.enforcepbc(atoms[a].pos, latvec, nshift);
    }

    centers.resize(atoms.size());
    for(int at=0; at < natoms; at++) {
      for(int i=0; i< 3; i++) centers[at].pos[i]=atoms[at].pos[i];
      centers[at].equiv_atom=at;
      centers[at].name=atoms[at].name;
    }
    string lcao_str=outputname+".pw2lcao";
    ofstream lcao(lcao_str.c_str());
    lcao << "latvec\n";
    for(int i=0; i< 3; i++) {
      for(int j=0; j< 3; j++) lcao << latvec[i][j] << "  ";
      lcao << endl;
    }
    lcao << "ref_latvec\n";
    for(int i=0; i< 3; i++) {
      for(int j=0; j< 3; j++) lcao << ref_latvec[i][j] << "  ";
      lcao << endl;
    }

    lcao << "origin ";
    for(int i=0; i< 3; i++) lcao <<  origin[i] << "   ";
    lcao << endl;
    lcao << "ecut " << ecut << endl;
    lcao << "center_file " << centeroutname << endl;
    lcao << "basis_file  " << lcao_basis << endl;
    lcao << "wf_file     " << wffilename << endl;
    lcao << "orb_file    " << orboutname << endl;
    lcao.close();

  }
  else {
    bool status;
    D3vector cell(latvec[0][0], latvec[1][1], latvec[2][2]);

    slwriter.write_centers=true;

    status = jeepbasis.resize(cell,cell,ecut);
    if ( !status )
    {
      cout << "basis allocation failed" << endl;
      exit(1);
    }

    centeroutname=outputname+".centers";
    slwriter.centername=centeroutname;
    string basisoutname=outputname+".basis";
    slwriter.basisname=basisoutname;



    //one center at the origin..
    Center tempcenter;
    tempcenter.name="PLANEWAVE_ORIGIN";
    tempcenter.equiv_atom=0;
    tempcenter.basis=0;
    centers.push_back(tempcenter);

    read_jeep_wf(wffilename, jeepbasis, moCoeff);

    vector <int> nbasis_list;
    nbasis_list.resize(1);
    nbasis_list[0]=2*jeepbasis.size(); //real and imaginary
    ofstream orbout(orboutname.c_str());
    print_orbitals(orbout, centers, nbasis_list, moCoeff);
    orbout.close();

    ofstream basisout(basisoutname.c_str());
    write_jeep_basis(basisout, jeepbasis);
    basisout.close();
  }


  ofstream centerout(centeroutname.c_str());
  centerout << centers.size() << endl;
  for(vector < Center>:: iterator i=centers.begin();
      i != centers.end(); i++) {
      i->print_center(centerout);
  }
  centerout.close();

  string slateroutname=outputname+".slater";
  ofstream slaterout(slateroutname.c_str());
  slwriter.print_wavefunction(slaterout);
  slaterout.close();

  //--------------------------Jastrow 2 output


  double basis_cutoff=find_basis_cutoff(latvec);
  string jast2outname=outputname+".jast2";
  Jastrow2_wf_writer jast2writer;
  jast2writer.set_atoms(atoms);
  
  
  ofstream jast2out(jast2outname.c_str());
  print_std_jastrow2(jast2writer, jast2out, basis_cutoff);
  jast2out.close();
  

  //-----------------------------------System output

  string sysoutname=outputname+".sys";
  ofstream sysout(sysoutname.c_str());
  sysout << "SYSTEM { PERIODIC \n";
  sysout << "  NSPIN { " << slwriter.nup << "  "
         << slwriter.ndown << " } \n";

  //If we're using the LCAO fits, we usually need at
  //least one unit cell layer of ghost centers surrounding the
  //simulation cell.  If we're using plane-waves, the
  //centers are unnecessary, and they're not used, so
  //we put a big divider in.
  //It would be nice to have an automatic determination
  //of this, like in crystal2qmc, but we need the gaussian
  //coefficients to do it.
  if(use_lcao) {
    sysout << "  CUTOFF_DIVIDER  1.00001" << endl;
  }
  else {
    sysout << "  CUTOFF_DIVIDER  5.00001" << endl;
  }

  sysout << "LATTICEVEC { \n";
  for(int i=0; i< 3; i++) {
    for(int j=0; j< 3; j++) {
      sysout << latvec[i][j] << "   ";
    }
    sysout << endl;
  }
  sysout << " } " << endl;

  sysout << "ORIGIN { ";
  for(int i=0; i< 3; i++) sysout << origin[i] << "  ";
  sysout << "}\n\n";

  for(vector <Atom>::iterator at=atoms.begin(); at != atoms.end(); at++) {
    at->print_atom(sysout);
  }
  sysout << "}\n\n\n";

  int npsp=pseudo.size();
  for(int psp=0; psp < npsp; psp++) {
    pseudo[psp].print_pseudo(sysout);
  }
  sysout.close();

  //---------------------------HF and OPT outputs

  string hfoutname=outputname+".hf";
  ofstream hfout(hfoutname.c_str());
  print_vmc_section(hfout, outputname, eref);
  hfout << "\n\n";
  hfout << "INCLUDE " << sysoutname << "  \n";
  hfout << "TRIALFUNC { INCLUDE " << slateroutname << "}\n\n";
  hfout.close();

  string optoutname=outputname+".opt";
  ofstream optout(optoutname.c_str());
  print_opt_section(optout, outputname, eref);
  optout << "\n\n";
  optout << "INCLUDE " << sysoutname << " \n";
  optout << "TRIALFUNC { \n  SLATER-JASTROW \n"
         << "  WF1 { INCLUDE " << slateroutname << " } \n"
         << "  WF2 { INCLUDE " << jast2outname   << " } \n"
         << "}\n\n";
  optout.close();


}


//-------------------------------------------------------------------


void read_jeep_sys(std::istream & is, vector < Atom > & atoms,
                   vector <Spline_pseudo_writer> & jeep_pseudo,
                   vector <double> & origin,
                   vector < vector < double> > & latvec,
                   vector < vector < double> > & ref_latvec) {
  vector <string> uniquenames;
  string dummy;
  vector <double> cell;
  vector <double> ref_cell;
  while(is >> dummy) {
    if(dummy == "set_cell") {
      double temp;
      cell.clear();
      for(int i=0; i< 3; i++) {
        if(!(is >> temp)) {
	  cout << "Unexpected end of sys file!  There should be three numbers after"
	       << " set_cell." << endl;
          exit(1);
	}
        cell.push_back(temp);
      }
    }
    else if(dummy=="set_ref_cell") {
      double temp;
      ref_cell.clear();
      for(int i=0; i< 3; i++) {
        if(!(is >> temp)) {
	  cout << "Unexpected end of sys file!  There should be three numbers after"
	       << " set_cell." << endl;
          exit(1);
	}
        ref_cell.push_back(temp);
      }
    }
    else if(dummy=="atom") {
      is >> dummy; //unique name; we don't need it..
      Atom tempatom;
      is >> tempatom.name;
      is >> tempatom.pos[0] >> tempatom.pos[1] >> tempatom.pos[2];
      atoms.push_back(tempatom);
      int nunique=uniquenames.size();
      bool unique=true;
      for(int i=0; i< nunique; i++) {
        if(tempatom.name==uniquenames[i]) {
          unique=false;
          break;
        }
      }
      if(unique) {
        uniquenames.push_back(tempatom.name);
      }

    }
  }


  if(cell.size() < 3) {
    cout << "Couldn't find set_cell in sys file!  Exiting" << endl;
    exit(1);
  }

  if(ref_cell.size() < 3) {
    ref_cell.clear();
    for(int i=0; i< 3; i++) ref_cell.push_back(cell[i]);
  }

  for(int i=0; i< 3; i++) origin.push_back(-cell[i]/2.0);

  latvec.resize(3);
  for(int i=0; i< 3; i++)
    for(int j=0; j< 3; j++)
       latvec[i].push_back(0.0);

  for(int i=0; i< 3; i++) latvec[i][i]=cell[i];


  ref_latvec.resize(3);
  for(int i=0; i< 3; i++)
    for(int j=0; j< 3; j++)
       ref_latvec[i].push_back(0.0);

  for(int i=0; i< 3; i++) ref_latvec[i][i]=ref_cell[i];

  //for(unsigned int at=0; at < atoms.size(); at++) {
  //  cout << atoms[at].name << "  " << atoms[at].pos[0] << "  " << atoms[at].pos[1]
  //  << "   " << atoms[at].pos[2] << endl;
  //}


  //There should be pseudopotential files for each of the atom names..
  string space=" ";
  const string comment="#";
  for(unsigned int file=0; file < uniquenames.size(); file++) {
    string pspfilename=uniquenames[file];
    ifstream pspfile(pspfilename.c_str());
    if(!pspfile) {
      cout << "Couldn't open " << pspfilename << endl;
      exit(1);
    }
    string line;
    Spline_pseudo_writer temppsp;
    temppsp.label=uniquenames[file];
    while(getline(pspfile, line)) {

      //Ignore comments until we find the first non-comment line
      unsigned int pos=line.find(comment);
      if( pos < line.size() )
      {
        line.erase(pos,line.size()-pos);
      }
      vector <string> words;
      split(line, space, words);
      if(words.size() >=2 ) {
        //cout << line << endl;
        int npoints=atoi(words[0].c_str());  //Number of points in this expansion
        double effcharge=atof(words[1].c_str()); //effective charge
        int numL=atoi(words[4].c_str());

        //Go ahead and assign the effective charge to the atoms.
        for(unsigned int at=0; at< atoms.size(); at++) {
          if(atoms[at].name==temppsp.label)
            atoms[at].charge=effcharge;
        }

        vector <double> positions; //for this l-value
        vector <double> values;
        positions.resize(npoints);
        values.resize(npoints);

        for(int i=0; i< numL; i++) {
          for(int point=0; point < npoints; point++) {
            pspfile >> positions[point] >> values[point];
          }
          temppsp.psp_pos.push_back(positions);
          temppsp.psp_val.push_back(values);

          if(i < numL-1) { //on all but the last one
            int npoints_temp=0;
            pspfile >> npoints_temp;
            string dummy;
            //pspfile >> dummy;
            //cout << "skipping " << dummy << endl;
            //if(dummy[0] == 'p') { //phi..  We just throw it away
              for(int p=0; p < npoints_temp+1; p++) { //number of lines + the newline left
                pspfile.ignore(180, '\n');
              }
              pspfile >> npoints_temp;
              //pspfile >> dummy;
            //}
            //if(dummy[0] != 'v') {
            //  cout << "expected v(something).  Instead, I got " << dummy << endl;
            //}
            pspfile.ignore(180, '\n');
            if(npoints != npoints_temp) {
              cout << "need the same number of points for all l-values" << endl;
              exit(1);
            }
            npoints=npoints_temp;

          }

        }

        //To get the psp's in the normal form, we:
        //Subtract the last l-value from the first ones
        //Add Z/r to the last one--this is done in the core code via a switch

        int lastL=numL-1;
        for(int i=0; i < lastL; i++) {
          for(int p=0; p < npoints; p++) {
            temppsp.psp_val[i][p]-= temppsp.psp_val[lastL][p];
          }
        }

        break;
      }
    }
    jeep_pseudo.push_back(temppsp);
  }



}

//----------------------------------------------------------------------


int read_int(FILE * file) {
  int i;
  fread(&i, sizeof(int), 1, file);
  return i;
}

double read_double(FILE * file) {
  double f;
  fread(&f, sizeof(double), 1, file);
  return f;
}

void skip(FILE * file, const int bytes) {
  int pos=ftell(file);
  fseek(file, pos+bytes, SEEK_SET);
}

/*
  (format comment copied from T. Ogitsu and B. Militzer)
Jeep wf format

spin up data

int = 0
int nst
double[nst]             = occ[n], n=0, .. , nst-1   --number of electrons in each orbital
int ngw
complex<double>[nst*ngw]  = wf(t)
complex<double>[nst*ngw]  = wf(t-dt)

next lines only if nspin == 2:

spin down data

int = 0
int nst
double[nst]             = occ[n], n=0, .. , nst-1
int ngw
complex<double>[nst*ngw]  = wf(t)
complex<double>[nst*ngw]  = wf(t-dt)

*/



void read_jeep_header(string & filename, vector <double> & occupation,
                      int & nspin, int & nst,int & nst2, int & ngw) {
  FILE * wffile= fopen(filename.c_str(), "r");
  if(wffile==0) { 
    cout << "Couldn't open " << filename << endl;
    exit(1);
  }
  
  int startpos=ftell(wffile);
  fseek(wffile, 0L, SEEK_END);
  int size=ftell(wffile);
  fseek(wffile, startpos, SEEK_SET);

  int test=read_int(wffile); //should be zero
  nst=read_int(wffile);
  occupation.resize(nst);
  for(int i=0; i< nst; i++ )
    occupation[i]=read_double(wffile);
  
  ngw=read_int(wffile);
  int nskip=4*nst*ngw*sizeof(double);
  cout << "nskip " << nskip << " size " << size << endl;
  skip(wffile, nskip);
  int pos=ftell(wffile);
  cout << "position " << pos << endl;
  if(pos==size) {
    fclose(wffile);
    nspin=1;
    return;
  }
  if(pos > size) {
    cout << "WARNING: count seems off in wf file " << endl;
  }

  //if we're this far, then there's a second spin
  test=read_int(wffile);
  assert(test==0);
  nst2=read_int(wffile);
  for(int i=0; i< nst2; i++) 
    occupation.push_back(read_double(wffile));

  nspin=2;
  fclose(wffile);
  return;

}


void read_jeep_wf(string & filename, RealBasis & jeepbasis,
     vector < vector <double> > & moCoeff) {

  vector <double> occupation;
  int nspin, ngw, nst, nst2;
  read_jeep_header(filename, occupation, nspin, nst,nst2, ngw);
  
  FILE * wffile=fopen(filename.c_str(), "r");
  int headerskip=2*sizeof(int) //nst, ngw
    +sizeof(int) //extra zero
    +nst*sizeof(double); //occupation numbers

  
  moCoeff.resize(nst*nspin);
  int currmo=0;
  for(int s=0; s< nspin; s++) {
    skip(wffile, headerskip);
    int nstt;
    if(s==0) 
      nstt=nst;
    else if(s==1)
      nstt=nst2;

    for(int i=0; i< nstt; i++) {
      for(int j=0; j< ngw; j++) {
        moCoeff[currmo].push_back(read_double(wffile)); //cos coefficient
        moCoeff[currmo].push_back(-read_double(wffile)); //sin coefficient
      }
      //Multiply the non-zero coefficients by 2, since we have only half the 
      //space with a real function
      for(int j=2; j< ngw*2; j++) {
        moCoeff[currmo][j]*=2;
      }
      currmo++;
    }

    int nskip=2*nst*ngw*sizeof(double);
    skip(wffile, nskip); //skip over the t-dt part
  }
  
  fclose(wffile);


  int nmo=moCoeff.size();
  //int ngw=jeepbasis.size();
  double sumkinetic=0;
  for(int mo=0; mo < nmo; mo++) {
    double kinetic=0;
    int counter=0;

    for(int i=0; i< ngw; i++) {
      double gx=jeepbasis.gx[i];
      double gy= jeepbasis.gx[ngw+i];
      double gz=jeepbasis.gx[ngw+ngw+i];
      double real=moCoeff[mo][counter++];
      double imag=moCoeff[mo][counter++];
      kinetic+=(real*real+imag*imag)*(gx*gx+gy*gy+gz*gz);

    }
    //cout << "kinetic energy for mo " << mo << "   :   " << .5*kinetic << endl;
    sumkinetic+=.5*kinetic;
  }

  cout << "Total kinetic energy : " << sumkinetic << "  (check in the jeep output) \n";
}
