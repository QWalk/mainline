/*
 
 Copyright (C) 2011 Lucas K. Wagner
 
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

#include "Average_region_fluctuation.h"
#include "Sample_point.h"
#include "ulec.h"

void Average_region_fluctuation::randomize(Wavefunction_data * wfdata, Wavefunction * wf,
                       System * sys, Sample_point * sample_tmp) { 
}

//----------------------------------------------------------------------

void Average_region_fluctuation::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & avg) { 
  avg.type="region_fluctuation";
  avg.vals.Resize(nregion*nregion*maxn*maxn*4);//the four is for the spins
  avg.vals=0.0;
  Array2 <int> n_in_regions;
  count_regions(sys,sample,n_in_regions);
  //for(int s=0; s< 2; s++) { 
  //  for(int r=0; r< nregion; r++) { 
  //    cout << "jj " << s << " " << r <<  " " << n_in_regions(s,r) << endl;
  //  }
  //}

  int maxn2=maxn*maxn;
  int nregion2=nregion*nregion;
  int nspin=2;
  for(int s1=0; s1 < nspin; s1++) { 
    for(int s2=0; s2< nspin; s2++) { 
      for(int r1=0; r1 < nregion; r1++) { 
        for(int r2=0; r2 < nregion; r2++) { 
          int n1=n_in_regions(s1,r1);
          int n2=n_in_regions(s2,r2);
          int index=n2+n1*maxn+r2*maxn2+r1*maxn2*nregion+s2*maxn2*nregion2+s1*maxn2*nregion2*nspin;
          //cout << "s1 " << s1 << " s2 " << s2 
          //  << " r1 " << r1 << " r2 " << r2 << " n1 " << n1 << " n2 " << n2 << " index "<< index << endl;
          avg.vals(index)=1.0;

        }
      }
    }
  }
}




//----------------------------------------------------------------------

void Average_region_fluctuation::read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words) { 
  unsigned int pos=0;

  int natoms=sys->nIons();
  if(!readvalue(words,pos=0,maxn,"MAXN"))
      maxn=20;
  nregion=natoms;

  
}

//----------------------------------------------------------------------

void Average_region_fluctuation::write_init(string & indent, ostream & os) { 
  os << indent << "REGION_FLUCTUATION" << endl;
  os << indent << "nregion " << nregion << endl;
  os << indent << "maxn " << maxn <<endl;
}
//----------------------------------------------------------------------
void Average_region_fluctuation::read(vector <string> & words) { 
  unsigned int pos=0;
  readvalue(words,pos=0,nregion,"NREGION");
  readvalue(words,pos=0,maxn, "MAXN");
  
}
//----------------------------------------------------------------------
void Average_region_fluctuation::write_summary(Average_return &avg,Average_return &err, ostream & os) { 
  int nspin=2;
  int maxn2=maxn*maxn;
  int nregion2=nregion*nregion;
  os << "Region fluctuation" << endl;
  os << "nspin " << nspin << " maxn " << maxn << " nregion " << nregion << endl;
  for(int s1=0; s1 < nspin; s1++) { 
    for(int s2=0; s2< nspin; s2++) { 
      for(int r1=0; r1 < nregion; r1++) { 
        for(int r2=0; r2 < nregion; r2++) { 
          os << "Spin " << s1 << " " << s2 << " Region " << r1 << " " << r2 << endl;
          for(int n1=0; n1 < maxn; n1++) { 
            for(int n2=0; n2 < maxn; n2++) { 
              int index=n2+n1*maxn+r2*maxn2+r1*maxn2*nregion+s2*maxn2*nregion2+s1*maxn2*nregion2*nspin;
              os <<  avg.vals(index) << " ";
            }
            for(int n2=0; n2 < maxn; n2++) { 
              int index=n2+n1*maxn+r2*maxn2+r1*maxn2*nregion+s2*maxn2*nregion2+s1*maxn2*nregion2*nspin;
              os << err.vals(index) << " ";
            }
            os << endl;
          }
          
        }
      }
    }
  }

}
//--------------------------------------------------------------------------

void Average_region_fluctuation::jsonOutput(Average_return &avg,Average_return &err, ostream & os) {
  int nspin=2;
  int maxn2=maxn*maxn;
  int nregion2=nregion*nregion;
  os << "\""<< avg.type <<"\":{" << endl;
  os << "\"nspin\":" << nspin << "," << endl;
  os << "\"maxn\":" << maxn <<"," << endl;
  os << "\"nregion\":" << nregion << "," << endl;
  os << "\"fluctuation data\":[" << endl;
  for(int s1=0; s1 < nspin; s1++) {
    for(int s2=0; s2< nspin; s2++) {
      for(int r1=0; r1 < nregion; r1++) {
        for(int r2=0; r2 < nregion; r2++) {
          os << "{" << endl;
          os << "\"spin\":[" << s1 << "," << s2 << "]," << endl;
          os << "\"region\":[" << r1 << "," << r2 << "]," << endl;
          os << "\"value\":[" << endl;
          for(int n1=0; n1 < maxn; n1++) {
            os << "[";
            for(int n2=0; n2 < maxn; n2++) {
              int index=n2+n1*maxn+r2*maxn2+r1*maxn2*nregion+s2*maxn2*nregion2+s1*maxn2*nregion2*nspin;
              os <<  avg.vals(index);
              if(n2 < maxn-1) os << ",";
            }
            if(n1==maxn-1) os << "]" << endl;
            else os << "]," << endl;
          }
          os << "]," << endl;
          
          os << "\"error\":[" << endl;
          for(int n1=0; n1 < maxn; n1++) {
            os << "[";
            for(int n2=0; n2 < maxn; n2++) {
              int index=n2+n1*maxn+r2*maxn2+r1*maxn2*nregion+s2*maxn2*nregion2+s1*maxn2*nregion2*nspin;
              os <<  err.vals(index);
              if(n2 < maxn-1) os << ",";
            }
            if(n1==maxn-1) os << "]" << endl;
            else os << "]," << endl;
          }
          os << "]" << endl;
          if(s1==nspin-1 && s2==nspin-1 && r1==nregion-1 && r2==nregion-1) os<<"}" << endl;
          else os << "}," << endl;
        }
      }
    }
  }
  os << "]" << endl;
  os << "}" << endl;
  
}


//--------------------------------------------------------------------------
//
void Average_region_fluctuation::count_regions(System * sys, 
    Sample_point * sample, Array2 <int> & n_in_regions) {
  n_in_regions.Resize(2,nregion);
  n_in_regions=0;

  int natoms=sys->nIons();
  assert(nregion==natoms+1);
  Array1 <int> nelec(2);
  for(int s=0;s < 2; s++) nelec(s)=sys->nelectrons(s);
  sample->updateEIDist();

  Array1 <doublevar> r(5);
  for(int s=0; s< 2; s++) { 
    for(int i=0; i< nelec(s); i++) { 
      int e=i+s*nelec(0);
      int atmin=0;
      doublevar rmin=1e99;
      for(int at=0; at < natoms; at++) { 
        
        sample->getEIDist(e,at,r);
        if(r(0)< rmin) {
          rmin=r(0);
          atmin=at;
        }
      }
      n_in_regions(s,atmin)++;
      
    }
  }

}


//#####################################################################
//----------------------------------------------------------------------

void Average_region_density_matrix::randomize(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample) {
  randomize(sys);
}

void Average_region_density_matrix::randomize(System * sys) { 
   int nsample=saved_r.GetDim(0);
   Array1 <doublevar> ran(3);
   for(int i=0; i < nsample; i++) { 
     saved_r(i)=0.0;

     for(int d=0; d< 3; d++) ran(d)=rng.ulec();
     for(int d1=0;d1 < 3; d1++) { 
       for(int d2=0; d2 < 3; d2++) { 
         saved_r(i)(d1)+=origin(d1)+latvec(d2,d1)*ran(d2);
       }
     }
   }
}
//----------------------------------------------------------------------

void Average_region_density_matrix::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                        System * sys, Sample_point * sample, Average_return & avg) { 
  avg.type="region_density_matrix";
  avg.vals.Resize(nmax*nmax*nregion*nregion*2*2);//the two is for the spins, two more for complex
  avg.vals=0.0;
  wf->updateVal(sample);
  Wf_return wfval_base(wf->nfunc(),2);
  wf->getVal(wfval_base);
  int nup=sys->nelectrons(0);
  int ndown=sys->nelectrons(1);
  int nelectrons=nup+ndown;

  Sample_point * sample_tmp=NULL;
  sys->generateSample(sample_tmp);
  
  Array1 <int> base_region(nelectrons);
  for(int e=0; e< nelectrons; e++) { 
    base_region(e)=which_region(sys,sample,e);
  }
  Array1<int> tot_n(nregion,0);
  for(int e=0; e< nelectrons; e++) {
    tot_n(base_region(e))++;
  }
  int nsample=saved_r.GetDim(0);
  for(int i=0; i < nsample; i++) { 
    Array1 <Wf_return> wf_eval;
    wf->evalTestPos(saved_r(i),sample,wf_eval);


    sample_tmp->setElectronPosNoNotify(0,saved_r(i));
    int region=which_region(sys,sample_tmp,0);


    for(int e=0; e< nelectrons; e++) { 
      dcomplex psiratio=exp(dcomplex(wf_eval(e).amp(0,0)-wfval_base.amp(0,0),
            wf_eval(e).phase(0,0)-wfval_base.phase(0,0)));
      int s=0;
      //if(e >=nup) s=1;
      
      int indx=2*(s*nregion*nregion*nmax*nmax
          +region*nregion*nmax*nmax+base_region(e)*nmax*nmax
          +tot_n(region)*nmax+tot_n(base_region(e)));
      avg.vals(indx)+=psiratio.real()/nsample;
      avg.vals(indx+1)+=psiratio.imag()/nsample;

    }

  }
  delete sample_tmp; 
}
//----------------------------------------------------------------------
int Average_region_density_matrix::which_region(System * sys, 
    Sample_point * sample,int e) { 

  int natoms=sys->nIons();
  Array1 <doublevar> r(5);
  int atmin=0;
  doublevar rmin=1e99;
  sample->updateEIDist();
  
  for(int at=0; at < natoms; at++) { 

    sample->getEIDist(e,at,r);
    if(r(0)< rmin) {
      rmin=r(0);
      atmin=at;
    }
  }

  return atmin;
}

//----------------------------------------------------------------------

void Average_region_density_matrix::read(System * sys, Wavefunction_data * wfdata, vector
                    <string> & words) { 
  unsigned int pos=0;
  int nsample=100;
  saved_r.Resize(nsample);
  for(int i=0; i< nsample; i++) saved_r(i).Resize(3);

  nmax=4;
  int natoms=sys->nIons();
  nregion=natoms;
  latvec.Resize(3,3);
  origin.Resize(3);
  if(sys->getBounds(latvec,origin)) { 
  }
  else { 
    Array1 <doublevar> mi(3),ma(3),pos(3);
    mi=1e99;
    ma=-1e99;
    for(int i=0; i< natoms; i++) { 
      sys->getIonPos(i,pos);
      for(int d=0;d < 3; d++) { 
        if(pos(d)< mi(d)) mi(d)=pos(d);
        if(pos(d)> ma(d)) ma(d)=pos(d);
      }
    }

    doublevar extra_space=2.0;
    latvec=0.0;
    for(int d=0; d< 3; d++) {
      origin(d)=mi(d)-extra_space;
      latvec(d,d)=ma(d)+extra_space;
    }
  }
  cout << "randomize " << endl;
  randomize(sys);
  cout << "done " << endl;
}
//----------------------------------------------------------------------


void Average_region_density_matrix::write_init(string & indent, ostream & os) { 
  os << indent << "REGION_DENSITY_MATRIX" << endl;
  os << indent << "nregion " << nregion << endl;
  os << indent << "nmax " << nmax << endl;
}
//----------------------------------------------------------------------
void Average_region_density_matrix::read(vector <string> & words) { 
  unsigned int pos=0;
  readvalue(words,pos=0,nregion,"NREGION");
  readvalue(words,pos=0,nmax,"NREGION");
  
}



void Average_region_density_matrix::write_summary(Average_return &avg,Average_return &err, ostream & os) { 


  for(int s=0; s< 2; s++) { 
    for(int rp=0; rp < nregion; rp++) { 
      for(int r=0; r< nregion; r++) { 
        for(int np=0; np < nmax; np++) { 
          for(int n=0; n< nmax; n++) { 
            int indx=2*(s*nregion*nregion*nmax*nmax
                +rp*nregion*nmax*nmax+r*nmax*nmax
                +np*nmax+n);
            os << "s " << s << " regions " << r << " -> " << rp << " numbers " << n << " " << np << " val " << avg.vals(indx) << " +/-  " << err.vals(indx) << endl;
          }
        }
      }
    }
  }
  
}
