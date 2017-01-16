/*
 
Copyright (C) 2008 Jindrich Kolorenc

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

#include "Blochwave_function.h"
#include "qmc_io.h"

/*!
   Note:  We may need to change this to read directly from a file, for memory
   constraints..
*/
int Blochwave_function::read(
  vector <string> & words,
  unsigned int & pos)
{
  // cout << "making Blochwave basis " << endl;
  unsigned int startpos=pos;

  centername=words[0];

  pos=startpos;
  vector <string> gvectxt;
  if(!readsection(words, pos, gvectxt, "GVECTOR")) {
    error("Need GVECTOR in BLOCHWAVE/PLANEWAVE input");
  }
  if(gvectxt.size() % 3 != 0) {
    error("Bad count on GVECTOR."
          "  There must be sets of three numbers.");
  }
  nmax=gvectxt.size()/3;
  // cout << nmax << " g-vectors " << endl;
  g_vector.Resize(nmax, 3);  
  int counter=0;
  for(int i=0; i< nmax; i++) {
    for(int d=0; d< 3; d++) {
      g_vector(i,d)=atof(gvectxt[counter].c_str());
      counter++;
    }
  }
  kmax=0;
  int b_max=0;
  vector < vector <string> > ckgtxt;      //text to store Ck(G) for all bands, all k points
  vector <string> ckgtmp;
  vector <string> cgtmp;
  vector <string> kvectxt;                //text to store k vectors        
  while(readsection(words, pos, ckgtmp, "KVECTOR")){
    for(int d=0; d< 3; d++)
      kvectxt.push_back(ckgtmp[d]);
    bmax.push_back(0);
    unsigned int pos_band=0;
    while(readsection(ckgtmp, pos_band, cgtmp, "BAND")){
      if( (cgtmp.size() !=  gvectxt.size()/3*2) ){ 
      // cout << kmax << "th k vector, " << bmax[kmax] << "th band:" << endl;
      error("Bad count of planewave coefficients Ck(G)"
            "  There must be same number of complex coefficients and G vectors.");
      }
      ckgtxt.push_back(cgtmp);
      bmax[kmax]++;
    }
    if( bmax[kmax] > b_max ) 
      b_max= bmax[kmax];       // keeping track of the largest number of bands for all k points
    if( bmax[kmax] == 0 ) {
      cout << "No BAND specifier for the " << kmax << "th k vector." << endl;
      error("Need BAND from input");
    }
   kmax++;
  }
  if( kmax == 0 )
    error("Need KVECTOR from input");
 
  k_vector.Resize(kmax, 3);
  ckg.Resize(kmax,b_max,nmax);  
  bn_total=0;
  int kvec_i=0;   // string counter for kvectxt
  for(int kn=0; kn< kmax; kn++) {
    for(int d=0; d< 3; d++) {
          k_vector(kn,d)=atof(kvectxt[kvec_i].c_str());
          kvec_i++;
    }
     // cout << kn+1 << "th numerical k vector stored" << endl;
    for(int bn=0; bn< bmax[kn]; bn++) {
      int ckg_i=0;    // string counter for ckgtxt, reset for every band
      for(int i=0; i< nmax; i++) {
        doublevar ckgr,ckgi; // real and imaginary parts of the coefficient
	//NB. some orbitals in Abinit have extra factor I
        ckgr=atof(ckgtxt[bn_total+bn][ckg_i].c_str()); ckg_i++;
        ckgi=atof(ckgtxt[bn_total+bn][ckg_i].c_str()); ckg_i++;
        ckg(kn,bn,i)= dcomplex(ckgr,ckgi);
      }
    }
    bn_total+=bmax[kn];
  }
  // cout << "reading done " << endl;
  return 0;  
}


void Blochwave_function::getVarParms(Array1 <doublevar> & parms) {
  //cout << "getVarParms " << endl;
  parms.Resize(0);
}

void Blochwave_function::setVarParms(Array1 <doublevar> & parms) {
  //cout << "setVArParms " << endl;
  assert(parms.GetDim(0)==0);
}

int Blochwave_function::nfunc()
{
  // cout << "Blochwave::nfunc() " << nmax << endl;
  return bn_total;
}

int Blochwave_function::showinfo(string & indent, ostream & os)
{
  os << indent << "Blochwave function\n";
  os << indent << kmax << " k points" << endl;
  os << indent << bn_total << " Blochwaves"  << endl;
  os << indent << nmax << " plane waves for each Bloch wave" << endl;
  return 1;
}

int Blochwave_function::writeinput(string & indent, ostream & os)
{
  os << indent << centername << endl;
  os << indent << "BLOCHWAVE\n";
  os << indent << "GVECTOR { \n";
  for(int i=0; i< nmax; i++) {
    os << indent << "  " << g_vector(i,0) << "  " << g_vector(i,1)
        << "  " << g_vector(i,2) << endl;
  }
  os << indent << "}" << endl;
  for(int kn=0; kn< kmax; kn++) {
    os << indent << "KVECTOR {\n";
    os << indent << "  " << k_vector(kn,0) << "  " << k_vector(kn,1)
          << "  " << k_vector(kn,2) << endl;
    for(int bn=0; bn< bmax[kmax]; bn++) {
      os << indent << "BAND {\n";
      for(int i=0; i< nmax; i++)
        os << indent << ckg(kn,bn,i).real() << "  " << ckg(kn,bn,i).imag() << endl;
      os << indent << "}" << endl;
    }
    os << indent << "}" << endl;
  }
  return 1;
}


void Blochwave_function::calcVal(const Array1 <doublevar> & r,
                                 Array1 <dcomplex> & symvals,
                                 const int startfill)
{
  //cout << "calcVal " << endl;
  assert(r.GetDim(0) >= 5);
  assert(symvals.GetDim(0) >= bn_total+startfill);
  int index=startfill;
  symvals=dcomplex(0.0,0.0);
  for(int kn=0; kn< kmax; kn++)
    for(int bn=0; bn< bmax[kn]; bn++) {
      for(int i=0; i< nmax; i++) {
        doublevar gdotr=0.0;
        for(int d=0; d< 3; d++)
          gdotr+=( g_vector(i,d)+k_vector(kn,d) )*r(d+2);
        symvals(index)+=ckg(kn,bn,i)*exp(I*gdotr);
      }
      index++;
   }
  //cout << "done" << endl;
}

void Blochwave_function::calcLap(
  const Array1 <doublevar> & r,
  Array2 <dcomplex> & symvals,
  const int startfill
)
{
  // cout << "calcLap " << endl;
  assert(r.GetDim(0) >= 5);
  assert(symvals.GetDim(0) >= bn_total+startfill);
  assert(symvals.GetDim(1) >= 5);

  int index=startfill;
  dcomplex t_exp;
  symvals=dcomplex(0.,0.);  
  for(int kn=0; kn< kmax; kn++)
    for(int bn=0; bn< bmax[kn]; bn++) {
      for(int fn=0; fn< nmax; fn++) {
        doublevar gdotr=0.0;
        doublevar gsquared=0.0;
        for(int d=0; d< 3; d++) {
          gdotr+=( g_vector(fn,d)+k_vector(kn,d) )*r(2+d);
          //Should probably store this one...
          gsquared+=(g_vector(fn,d)+k_vector(kn,d))*(g_vector(fn,d)+k_vector(kn,d));
        }
        t_exp=exp(I*gdotr);
	//t_exp=dcomplex(1.0,0.0);
        symvals(index,0)+=ckg(kn,bn,fn)*t_exp;
        for(int i=1; i< 4; i++) {
          symvals(index,i)+=ckg(kn,bn,fn)*I*(g_vector(fn, i-1)+k_vector(kn, i-1))*t_exp;
        }
        symvals(index,4)+=-ckg(kn,bn,fn)*gsquared*t_exp;
      }
      // cout << index << symvals(index,0) << " Ck(0) " << ckg(kn,bn,0) << endl;
      index++;
    }
  //cout << "done" << endl;
}



//------------------------------------------------------------------------
