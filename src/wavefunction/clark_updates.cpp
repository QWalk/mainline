/*
 
Copyright (C) 2012 Lucas K. Wagner

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

#include "clark_updates.h"
void build_excitation_list(Array3 <Array1 <int> > & occupation,int f,//(function,det,spin) (orb #)
    Array1 <Excitation> & ex) { 
  int ndet=occupation.GetDim(1);
  int ns=occupation.GetDim(2);
  Array1 <int> nspin(ns);
  for(int s=0; s< ns; s++) nspin(s)=occupation(0,0,s).GetDim(0);
  ex.Resize(ndet);
  //Assume that the base determinant is the first one!
  int base=0;
  for(int d=0; d< ndet; d++) { 
    //Find the differences with the base determinant
    vector <vector <int> > tot_missing(ns);
    vector <vector <int> > tot_additional(ns);
    //This search could likely be improved, but be careful with ordering,
    for(int s=0; s< ns; s++) { 
      for(int e1=0; e1 < nspin(s); e1++) { 
        int o1=occupation(0,base,s)(e1);
        bool found=false;
        for(int e2=0; e2 < nspin(s); e2++) { 
          if(occupation(0,d,s)(e2)==o1) {
            found=true;
            break;
          }
        }
        if(!found)tot_missing[s].push_back(o1);
        int o2=occupation(0,d,s)(e1);
        found=false;
        for(int e2=0; e2 < nspin(s); e2++) {
          if(occupation(0,base,s)(e2)==o2) { 
            found=true;
            break;
          }
        }
        if(!found) { 
          tot_additional[s].push_back(o2);
        }

      }
    }
    //The update formulation works in terms of replacing an occupied orbital with a
    //virtual one.  By default, QWalk converters order determinants by orbital number.
    //These representations are equivalent up to a sign, which we determine here in 
    //a fairly general way

    ex[d].sign.Resize(ns);
    ex[d].sign=1;
    for(int s=0; s< ns; s++) { 
      Array1<int> newocc=occupation(0,base,s);
      int nex=tot_additional[s].size();
      for(int e=0; e< nspin(s); e++) { 
        for(int ex=0; ex< nex; ex++)  {
          if(newocc[e]==tot_missing[s][ex])
            newocc[e]=tot_additional[s][ex];
        }
      }
      int count=0;
      for(int e1=0; e1 < nspin(s); e1++) {
        for(int e2=e1+1; e2 < nspin(s); e2++) { 
          if(newocc[e2]<newocc(e1)) count++;
          if(occupation(0,d,s)(e2) < occupation(0,d,s)(e1)) count--;
        }
      }
      cout << "d " << d << " count " << count << endl;
      if(count%2==1) ex[d].sign(s)*=-1;
    }


    //tot_missing and tot_additional should be filled now.
    ex[d].g.Resize(ns);
    ex[d].e.Resize(ns);

    for(int s=0;s < ns; s++) { 
      int nex_s=tot_missing[s].size();
      if(nex_s != tot_additional[s].size()) { 
        error("In build_excitation_list: different numbers of holes and particles");
      }
      ex[d].g[s].Resize(nex_s);
      ex[d].e[s].Resize(nex_s);
      for(int i=0; i< nex_s; i++) {
        ex[d].g[s][i]=tot_missing[s][i];
        ex[d].e[s][i]=tot_additional[s][i];
        cout << "d " << d << " s " << s << " i " << i << " : " << ex[d].g[s][i] << " -> " << ex[d].e[s][i] << " sign " << ex[d].sign(s) <<  endl;
      }
    }

  }

}


