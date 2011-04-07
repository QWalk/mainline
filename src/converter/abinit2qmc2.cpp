/*
 
Copyright (C) 2009 Lucas Wagner

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
#include "Pseudo_writer.h"
#include <fstream>
#include <cmath>
#include "elements.h"
using namespace std;


class Abinit_converter { 
  public:
    void readfile(string filename);
    void write_orbitals(string  filename);
    void write_sys(string filename);
    void write_slater(string filename);
  private:
    void prune_coefficients();  //Cancel out the G and -G degeneracies
    void add_psp(string filename);
    void evaluate_orbital(vector<double> & r, vector <complex <double> > & orbitals);
    void reassign_z(); 
    Slat_wf_writer slater;
    Jastrow_wf_writer jast;
    vector <Atom> atoms;
    vector <vector<double> > latvec;
    vector <Spline_pseudo_writer> pspspline; //one for each atom type
    vector < vector <double> > gvec;
    vector < vector <complex < double> > > coeff;
    vector<double> kpoint;
};

int main(int argc, char ** argv) { 
  Abinit_converter abconverter;
  abconverter.readfile(argv[1]);
  return 0;
}

//--------------------------------------------------------------

//Only read in the information from the converter output file.  All other processing goes in other
//functions.
void Abinit_converter::readfile(string  filename) { 
  ifstream is(filename.c_str());
  string line;
  vector <string> words;
  int nelectrons=0;
  string sep=" ";
  while(getline(is,line)) { 
    words.clear();
    if(line=="Number of electrons per primitive cell") {
      getline(is,line);
      nelectrons=atoi(line.c_str());
      cout <<  "number of electrons " <<  nelectrons << endl;
      slater.nup=nelectrons/2;
      slater.ndown=nelectrons/2;
    }
    else if(line=="GEOMETRY") { 
      getline(is,line); getline(is,line); getline(is,line);
      int nat=atoi(line.c_str());
      cout << "Number of atoms " << nat << endl;
      getline(is,line);
      for(int i=0; i< nat; i++) { 
        words.clear();
        getline(is,line);
        cout << "line " << line; 
        split(line,sep,words);
        Atom at;
        at.charge=atoi(words[0].c_str());
        for(int d=0;d < 3; d++) at.pos[d]=atof(words[d+1].c_str());
        at.name=element_lookup_caps[int(at.charge)];
        atoms.push_back(at);
      }
    }
    else if(line=="Primitive lattice vectors (au)") { 
      latvec.resize(3);
      for(int d1=0; d1< 3; d1++) { 
        latvec[d1].resize(3);
        words.clear();
        getline(is,line);
        split(line,sep,words);
        for(int d2=0; d2 < 3; d2++) latvec[d1][d2]=atof(words[d2].c_str());
      }
    }
    else if(line=="G VECTORS") { 
      getline(is,line); getline(is,line); getline(is,line);
      int ngvec=atoi(line.c_str());
      getline(is,line); cout << line << endl; assert(line=="Gx Gy Gz (au)");
      gvec.resize(ngvec);
      for(vector< vector<double> >::iterator i=gvec.begin(); i!=gvec.end(); i++) { 
        getline(is,line);
        words.clear(); split(line,sep,words);
        for(vector<string>::iterator j=words.begin(); j!=words.end(); j++) 
          i->push_back(atof(j->c_str()));
      }
    }
    else if(line=="WAVE FUNCTION") { 
      getline(is,line); getline(is,line); getline(is,line);
      int nkpts=atoi(line.c_str());
      assert(nkpts==1);
      getline(is,line); getline(is,line);
      words.clear(); split(line,sep,words);
      int kptnum=atoi(words[0].c_str());
      int nbandsup=atoi(words[1].c_str());
      int nbandsdown=atoi(words[2].c_str());
      for(int d=0; d< 3; d++) kpoint.push_back(atof(words[3+d].c_str()));
      assert(nbandsdown==0);
      coeff.resize(nbandsup);
      for(vector<vector < complex < double> > >::iterator band=coeff.begin(); band!=coeff.end(); band++) {
        getline(is,line); getline(is,line); getline(is,line);
        cout << line << endl;
        assert(line=="Eigenvector coefficients");
        band->resize(gvec.size());
        for(vector< complex <double> >::iterator c=band->begin(); c!=band->end(); c++) { 
          is >> *c;
        }
        is.ignore(180,'\n');
      }
    }
  }

  prune_coefficients();
  add_psp("pp.Si.burkatzki.fhi");
  reassign_z();
  write_orbitals("tmp.cube");
  write_sys("qwalk.sys");
  write_slater("qwalk.slater");
}


//--------------------------------------------------------------

void Abinit_converter::prune_coefficients() { 
   //In principle, we can save a bit of time by taking out the G,-G pairs.
   //For the moment, I won't worry about it.
   //We may be able to save a lot of time by making the wave functions real
   //instead of detecting it post facto, as well.
}

//--------------------------------------------------------------

void Abinit_converter::evaluate_orbital(vector<double> & r, 
    vector <complex< double> > & orbitals) { 
  assert(orbitals.size()==coeff.size());
  int ngvec=gvec.size();
  vector < complex<double> > basis(ngvec);
  vector <complex <double> >::iterator b=basis.begin();
  for(vector<vector<double> >::iterator g=gvec.begin(); g!=gvec.end(); g++,b++) { 
    double dot=(*g)[0]*r[0]+(*g)[1]*r[1]+(*g)[2]*r[2];
    *b=complex<double>(cos(dot),sin(dot));
  }
  
  vector< complex < double> >::iterator o=orbitals.begin();
  vector<vector<complex <double> > >::iterator co=coeff.begin();
  for(;o!=orbitals.end() && co!=coeff.end(); co++,*o++) { 
    *o=0.0;
    vector <complex<double> >::iterator co_j=co->begin();
    b=basis.begin();
    for(;b!=basis.end() && co_j!=co->end(); b++,co_j++) { 
      *o+=(*co_j)*(*b);
    }
  }
}

//--------------------------------------------------------------

void Abinit_converter::write_orbitals(string  filename) { 
  vector <int> npts(3);
  
  npts[0]=npts[1]=npts[2]=80;
  vector <complex < double> > orbitals(coeff.size());
  vector <double> res(3);
  for(int d=0;d < 3; d++) {
    res[d]=1.0/npts[d];
  }

  vector<double> r(3),prop(3);
  int count=0;
  long int ntotalpts=npts[0]*npts[1]*npts[2]*orbitals.size();
  cout << "size of all orbitals: " << ntotalpts*8.0/1024.0/1024.0 << " megabytes ("<< ntotalpts << ") total points" << endl;
  complex<double> * allorbitals=new complex<double>[ntotalpts];
  vector <complex<double> > orbnorm(coeff.size());
  for(vector<complex<double> >::iterator i=orbnorm.begin(); i!=orbnorm.end(); i++) 
    *i=0.0;

  for(int ii=0; ii < npts[0]; ii++) { 
    cout << "ii " << ii << endl;
    for(int jj=0; jj< npts[1]; jj++) { 
      for(int kk=0; kk< npts[2]; kk++) { 
        for(int i=0; i< 3; i++) r[i]=0.0;
        prop[0]=ii*res[0]; prop[1]=jj*res[1]; prop[2]=kk*res[2];
        for(int i=0; i< 3; i++) { 
          for(int j=0; j< 3; j++) { 
            r[i]+=prop[j]*latvec[i][j];
            //os << r[j] << " ";
          }
        }
        evaluate_orbital(r,orbitals);
        complex<double> * ptr=allorbitals+kk+jj*npts[2]+ii*npts[2]*npts[1];
        vector< complex <double> >::iterator i=orbitals.begin(),orbn=orbnorm.begin();
        for(; i!= orbitals.end(); i++,orbn++) { 
          *ptr=*i;
          *orbn+=complex<double>(i->real()*i->real(),i->imag()*i->imag());
          ptr+=npts[1]*npts[2]*npts[0];
        }
        //cout << r[0] <<  " " << r[1] << " " << r[2] << " "
        //  << orbitals[15] << endl;
        //cout << orbitals[0] << " " ;
        //if(count%6==5) cout << endl;
        count++;
      }
    }
  }

  for(vector<complex<double> >::iterator orbn=orbnorm.begin(); orbn!=orbnorm.end();
      orbn++) cout <<"norm " << *orbn << endl;

  ofstream os(filename.c_str());
  os << "Orbitals file" << endl;
  os << "norbitals " << orbitals.size() << endl;
  os << "Lattice vectors " << endl;
  for(int i=0; i< 3; i++) { 
    for(int j=0; j< 3; j++) { 
      os << latvec[i][j] << " ";
    }
    os << endl;
  }
  os << "resolution ";
  for(int i=0; i< 3; i++) os << res[i] << " ";
  os << endl;
  os << "npoints ";
  for(int i=0; i< 3; i++) os << npts[i] << " ";
  os << endl;
  os << "orbitals follow (orbital,x,y,z indices) " << endl;
  complex<double> * ptr=allorbitals;
  for(int orb=0; orb < orbnorm.size(); orb++) { 
    int nxyz=npts[0]*npts[1]*npts[2];
    if(orbnorm[orb].real() > orbnorm[orb].imag())  { 
      cout << "orb " <<orb << " real " << endl;
      for(int i=0; i< nxyz; i++) {
        os << ptr->real() << "\n";
        ptr++;
      }
    }
    else for(int i=0; i< nxyz; i++) { 
      os << ptr->imag() << "\n";
      ptr++;
    }
  }
  delete [] allorbitals;
  
}
//--------------------------------------------------------------
void Abinit_converter::write_sys(string filename) { 
  ofstream os(filename.c_str());
  os << "System { periodic \n";
  os << "  nspin { " << slater.nup << " " << slater.ndown << " } \n";
  for(vector<Atom>::iterator at=atoms.begin(); at!=atoms.end(); at++) { 
    at->print_atom(os);
  }
  os << "LatticeVec { \n";
  for(int i=0; i< 3; i++) { 
    for(int j=0; j< 3; j++) { 
      os << latvec[i][j] << " ";
    }
    os << endl;
  }
  os <<  "}\n";
  os << "   kpoint { 0 0 0 } \n";
  
  os << "} ";


  for(vector<Spline_pseudo_writer>::iterator i=pspspline.begin(); i!=pspspline.end();i++) { 
    i->print_pseudo(os);
  }
}


//--------------------------------------------------------------
void Abinit_converter::write_slater(string filename) { 
  ofstream os(filename.c_str());
  slater.orbname="tmp.cube";
  slater.calctype="RHF";
  slater.print_wavefunction(os);

}

//--------------------------------------------------------------

void Abinit_converter::add_psp(string filename) { 
  ifstream is(filename.c_str());
  string line;
  getline(is,line);
  getline(is,line);
  vector<string> words; string sep=" ";
  split(line,sep,words);
  int zatom=atoi(words[0].c_str());
  int zion=atoi(words[1].c_str());
  getline(is,line); words.clear(); split(line,sep,words);
  int npoints=atoi(words[4].c_str());
  int lmax=atoi(words[2].c_str());
  int lloc=atoi(words[3].c_str());
  cout << "lmax " << lmax << " lloc " << lloc << " zatom " << zatom << " zion " << zion << endl;
  assert(lloc == lmax); //for now, we'll assume the highest is  the local one
                        //I think that QWalk also assumes this at the time of this writing, so 
                        //we'd have to change it there too.
  Spline_pseudo_writer psp;
  for(int l=0; l< lmax; l++) { 
    //find the start of the values
    double rcut=0;
    while(getline(is,line)) { 
      words.clear(); split(line,sep,words);
      if(atoi(words[0].c_str())==npoints) { 
        rcut=atof(words[1].c_str());
        break;
      }
    }
    //
    vector<double> pos(npoints),val(npoints);
    vector<double>::iterator p=pos.begin(),v=val.begin();
    for(;p!=pos.end() && v!=val.end(); p++,v++) { 
      words.clear(); getline(is,line); split(line,sep,words);
      *p=atof(words[1].c_str());
      *v=atof(words[3].c_str());
    }
    vector <double> upos,uval;
    make_uniform(pos,val,upos,uval,0.01);
    psp.psp_pos.push_back(upos);
    psp.psp_val.push_back(uval);
  }
  psp.label=element_lookup_caps[int(zatom)];
  psp.effcharge=zion;
  int newnpoints=psp.psp_pos[0].size();
  for(int i=0; i< lloc-1; i++) { 
    assert(psp.psp_pos[i].size()==psp.psp_pos[lloc-1].size());
    for(int j=0; j< newnpoints; j++) {
      assert(fabs(psp.psp_pos[i][j]-psp.psp_pos[lloc-1][j]) < 1e-5);
      psp.psp_val[i][j]-=psp.psp_val[lloc-1][j];
    }
  }

  psp.print_pseudo(cout);
  pspspline.push_back(psp);

}
//--------------------------------------------------------------
void Abinit_converter::reassign_z() { 
  for(vector<Spline_pseudo_writer>::iterator p=pspspline.begin(); 
      p!=pspspline.end(); p++) { 
    for(vector<Atom>::iterator a=atoms.begin(); 
        a!=atoms.end(); a++) {
      if(a->name==p->label) { 
        a->charge=p->effcharge;
        a->print_atom(cout);
      }
    }
  }
}
//--------------------------------------------------------------
