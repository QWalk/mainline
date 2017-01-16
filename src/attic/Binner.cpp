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
//src/Binner.cpp
//Author: Lucas Wagner


#include "Binner.h"
#include "Sample_point.h"
#include "qmc_io.h"
#include "System.h"

void allocate(vector <string> & words, Binner * & bin, System * sys) {
  if(words[0]=="EEBIN") 
    bin=new EE_bin;
  else error("unknown bin type: ", words[0]);
  bin->read(words, sys);
}

//----------------------------------------------------------------------

void EE_bin::read(vector <string> & words, System * sys) {
  unsigned int pos=0;
  if(!readvalue(words, pos=0, range, "RANGE"))
    error("Need RANGE variable in EEBIN");
  if(!readvalue(words, pos=0, width, "WIDTH"))
    error("Need WIDTH in EEBIN");
  if(!readvalue(words, pos=0, output_file, "OUTPUT"))
    error("Need OUTPUT in EEBIN");

  int nbins=int(ceil(range/width));
  eebin.Resize(nbins);
  eebin=0;
  eebin_like.Resize(nbins);
  eebin_unlike.Resize(nbins);
  eebin_like=0;
  eebin_unlike=0;
  npoints=0;
  high_up_spin=sys->nelectrons(0);
}

//----------------------------------------------------------------------

void EE_bin::accumulate(Sample_point * sample) {
  Array1 <doublevar> dis(5);
  sample->updateEEDist();
  nelec=sample->electronSize();
  npoints++;
    
  for(int i=0; i< nelec; i++) {
    for(int j=i+1; j< nelec; j++) {
      sample->getEEDist(i,j,dis);
      if(dis(0) < range) {
        int bin=int(dis(0)/width);
        eebin(bin)++;
        if( (i < high_up_spin) == (j < high_up_spin) ) 
          eebin_like(bin)++;
        else
          eebin_unlike(bin)++;
      }
    }
  }
}

//----------------------------------------------------------------------

void EE_bin::write() {
  ofstream os(output_file.c_str());
  int nbins=eebin.GetDim(0);

  
  int ntot=0; //!< Number of pair distances
  int nltot=0; //!< number of like pair distances
  int nutot=0; //!< number of unlike pair distances
  
  //Ok, this is a stupid, but foolproof way to figure out the 
  //normalization.  There's probably a general way to do it..
  for(int i=0; i< nelec; i++) {
    for(int j=i+1; j< nelec; j++) {
      ntot++;
      if( (i < high_up_spin) == (j < high_up_spin) ) 
        nltot++;
      else
        nutot++;
    }
  }
      
  
  
  doublevar norm=ntot*3.0/(4.0*pi*range*range*range);
  for(int i=0; i< nbins; i++) {
    os << i*width << "   " << eebin(i)/(npoints*norm) << endl;
  }
  os.close();
  
  
  string likeoutname=output_file+"l";
  ofstream osl(likeoutname.c_str());
  norm=nltot*3.0/(4.0*pi*range*range*range);
  for(int i=0; i< nbins; i++) {
    osl << i*width << "   " << eebin_like(i)/(npoints*norm) << endl;
  }
  osl.close();
  
  string unlikeoutname=output_file+"u";
  ofstream osu(unlikeoutname.c_str());
  norm=nutot
       *3.0/(4.0*pi*range*range*range);
  for(int i=0; i< nbins; i++) {
    osu << i*width << "   " << eebin_unlike(i)/(npoints*norm) << endl;
  }
  osu.close();
  
}

//----------------------------------------------------------------------

