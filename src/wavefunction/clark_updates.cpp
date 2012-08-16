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
void Excitation_list::build_excitation_list(Array3 <Array1 <int> > & occupation,int f) { //(function,det,spin) (orb #)
  int ndet=occupation.GetDim(1);
  int ns=occupation.GetDim(2);
  Array1 <int> nspin(ns);
  for(int s=0; s< ns; s++) nspin(s)=occupation(0,0,s).GetDim(0);
  excitations.Resize(ndet);
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

    excitations[d].sign.Resize(ns);
    excitations[d].sign=1;
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
      if(count%2==1) excitations[d].sign(s)*=-1;
    }


    //tot_missing and tot_additional should be filled now.
    excitations[d].g.Resize(ns);
    excitations[d].e.Resize(ns);

    for(int s=0;s < ns; s++) { 
      int nex_s=tot_missing[s].size();
      if(nex_s != tot_additional[s].size()) { 
        error("In build_excitationscitation_list: different numbers of holes and particles");
      }
      excitations[d].g[s].Resize(nex_s);
      excitations[d].e[s].Resize(nex_s);
      for(int i=0; i< nex_s; i++) {
        excitations[d].g[s][i]=tot_missing[s][i];
        excitations[d].e[s][i]=tot_additional[s][i];
        cout << "d " << d << " s " << s << " i " << i << " : " << excitations[d].g[s][i] << " -> " << excitations[d].e[s][i] << " sign " << excitations[d].sign(s) <<  endl;
      }
    }

  }


  //Now we build some of the lookup variables
  int nex=excitations.GetDim(0);
  alle.resize(ns);
  allg.resize(ns);
  for(int e=0; e< nex; e++) { 
    for(int s=0; s< ns; s++) {  
      for(int i=0; i< excitations(e).g(s).GetDim(0); i++) { 
        int g=excitations(e).g(s)(i);
        bool found=false;
        for(vector<int>::iterator gi=allg[s].begin(); gi!=allg[s].end(); gi++) {
          if(g==*gi) {
            found=true;
            break;
          }
        }
        if(!found) allg[s].push_back(g);
        int ex=excitations(e).e(s)(i);
        found=false;
        for(vector<int>::iterator ei=alle[s].begin(); ei!=alle[s].end(); ei++) {
          if(ex==*ei) {
            found=true;
            break;
          }
        }
        if(!found) alle[s].push_back(ex);
      }
    }
  }

  

  //Array1 <Excitation> remap(nex);
  remap.Resize(nex);
  for(int e=0; e< nex; e++) { 
    remap(e).g.Resize(ns);
    remap(e).e.Resize(ns);
    for(int s=0; s< ns; s++) {
      int n=excitations(e).g(s).GetDim(0);

      remap(e).g(s).Resize(n);
      remap(e).e(s).Resize(n);

      int ng=allg[s].size();
      int ne=alle[s].size();
      for(int i=0; i< n; i++) { 
        for(int j=0; j< ng; j++) { 
          if(allg[s][j]==excitations(e).g(s)(i)) {
            remap(e).g(s)(i)=j;
            break;
          }
        }
        for(int j=0; j< ne; j++) { 
          if(alle[s][j]==excitations(e).e(s)(i)) { 
            remap(e).e(s)(i)=j;
            break;
          }
        }
      }
    }
  }

  
}


