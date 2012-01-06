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
  os.precision(3);
  cout << "avg.vals " << avg.vals(0) << endl;
  for(int s1=0; s1 < nspin; s1++) { 
    for(int s2=0; s2< nspin; s2++) { 
      for(int r1=0; r1 < nregion; r1++) { 
        for(int r2=0; r2 < nregion; r2++) { 
          os << "Spin " << s1 << " " << s2 << " Region " << r1 << " " << r2 << endl;
          for(int n1=0; n1 < maxn; n1++) { 
            for(int n2=0; n2 < maxn; n2++) { 
              int index=n2+n1*maxn+r2*maxn2+r1*maxn2*nregion+s2*maxn2*nregion2+s1*maxn2*nregion2*nspin;
              os << setw(12) << avg.vals(index);
            }
            for(int n2=0; n2 < maxn; n2++) { 
              int index=n2+n1*maxn+r2*maxn2+r1*maxn2*nregion+s2*maxn2*nregion2+s1*maxn2*nregion2*nspin;
              os << setw(12) << err.vals(index);
            }
            os << endl;
          }
          
        }
      }
    }
  }

}

//----------------------------------------------------------------------
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

  /*
  Array1 <int> nelec_found(2);
  nelec_found=0;
  for(int s=0;s < 2; s++) {
    for(int at=0; at< natoms; at++) { 
      nelec_found(s)+=n_in_regions(s,at);
    }
    n_in_regions(s,natoms)=nelec(s)-nelec_found(s);
  }
  */


}

//----------------------------------------------------------------------

