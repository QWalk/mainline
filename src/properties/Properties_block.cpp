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
#include "Properties_block.h"
#include <iomanip>
#include "qmc_io.h"

namespace Properties_types {
const char*  names[]={"total_energy", 
            "kinetic",
            "potential",
            "nonlocal",
            "weight"};
};

void print_array_sec(ostream & os, string & indent, 
                     const char * label, Array1 <doublevar>  arr) {
  os << indent << label << " { ";
  for(int i=0; i< arr.GetDim(0); i++) {
    os << arr(i) << "   ";
  }
  os << "} " << endl;
}


void read_array(vector<string> & words, unsigned int pos, 
                Array1 <doublevar> & array, const char * label) {
  vector <string> section;
  readsection(words, pos, section, label);
  array.Resize(section.size());
  int ns=section.size();
  for(int i=0; i< ns; i++) array(i)=atof(section[i].c_str());
}

int read_into_array(vector <string> & words, unsigned int pos,
                Array2 <doublevar> & array, const char * label, int place) {
  vector <string> section;
  if(!readsection(words, pos, section, label)) return 0;
  if(array.GetDim(1) <  int(section.size()) )
    error("In read_into_array, the section array is too large");
  if(array.GetDim(0) < int(place)) 
    error("in read_into array, place is too large");

  int ns=section.size();
  for(int i=0; i< ns; i++) array(place,i)=atof(section[i].c_str());
  return 1;
}

void Properties_block::storeToLog(string & indent, ostream & os,
                                    string & label) {

  using namespace Properties_types;

  int nwf=avg.GetDim(1);
  int naux=aux_energy.GetDim(0);
  //int autocorr_depth=autocorr.GetDim(1);

  os << indent << "label " << label << endl;
  os << indent << "totweight " << totweight << endl;
  os << indent << "nautocorr " << autocorr.GetDim(1) << endl;
  print_array_sec(os, indent, "total_energy", avg(total_energy));
  print_array_sec(os, indent, "total_energyvar", var(total_energy));
  print_array_sec(os, indent, "weight", avg(weight));
  print_array_sec(os, indent, "weightvar", var(weight));
  print_array_sec(os, indent, "kinetic_energy",avg(kinetic) );
  print_array_sec(os, indent, "kinetic_energyvar", var(kinetic));
  print_array_sec(os, indent, "potential_energy", avg(potential));
  print_array_sec(os, indent, "potential_energyvar", var(potential));
  print_array_sec(os, indent, "nonlocal_energy", avg(nonlocal));
  print_array_sec(os, indent, "nonlocal_energyvar", var(nonlocal));
  
  int navg_vals=avgrets.GetDim(1);
  for(int w=0; w< nwf; w++) { 
    for(int i=0; i< navg_vals; i++) { 
      os << indent <<  "average_generator { " << avgrets(w,i).type << " ";
      for(int j=0; j< avgrets(w,i).vals.GetDim(0); j++) { 
        os << avgrets(w,i).vals(j) << " ";
      }
      os << " } " << endl;
    }
  }
  
  cout << "aux " << endl; 
  string indent2=indent+"  ";
  for(int i=0; i< naux; i++) {
    os << indent << "aux_function  { " << endl;
    os << indent2 << "aux_size " << aux_size(i) << endl;
    int n_cvg=aux_diff.GetDim(1);
    os << indent2 << "nconverge " << n_cvg << endl;
    print_array_sec(os, indent2, "aux_energy", aux_energy(i));
    print_array_sec(os, indent2, "aux_energyvar", aux_energyvar(i));
    Array1 <doublevar> aux_diff_tmp(n_cvg), diffvar(n_cvg);
    
    
    
    for(int w=0; w< n_cvg; w++) {
      aux_diff_tmp(w)=//(aux_energy(i,w)-avg(total_energy,w))/aux_size(i);
          aux_diff(i,w)/aux_size(i);
      diffvar(w)=aux_diffvar(i,w)/aux_size(i);
    }
    print_array_sec(os, indent2, "aux_diff", aux_diff_tmp);
    print_array_sec(os, indent2, "aux_diffvar", diffvar);
    print_array_sec(os, indent2, "aux_weight", aux_weight(i));
    print_array_sec(os, indent2, "aux_weightvar", aux_weightvar(i));
    print_array_sec(os, indent2, "aux_autocorr", aux_autocorr(i));

    os << indent << "}" << endl;
  }
  cout << "autocorr " << endl;
  print_array_sec(os, indent, "autocorr_energy", autocorr(0));
  if(fabs(aux_weight_correlation) > 1e-10)
    os << indent << "aux_weight_correlation " << aux_weight_correlation 
       << endl;
  cout << "done " << endl;
}

//----------------------------------------------------------------------

void Properties_block::restoreFromLog(vector <string> & words) {
  
  vector <string> tmp;
  vector <vector <string > > auxtext;
  int nwf, naux;
  unsigned int pos=0;
  using namespace Properties_types;

  if(!readsection(words, pos=0,tmp, "total_energy"))
    error("error reading total_energy from log file");
  nwf=tmp.size();


  

  pos=0;
  while(readsection(words, pos, tmp, "aux_function")) 
    auxtext.push_back(tmp);
  naux=auxtext.size();
  int n_cvg=1;
  if(naux > 0) {
    readvalue(auxtext[0], pos=0, n_cvg, "NCONVERGE");
  }
  
  setSize(nwf, naux, n_cvg);

  if(!readvalue(words, pos=0, totweight, "totweight"))
    totweight=1;  

  read_into_array(words, pos=0, avg, "total_energy", total_energy);
  read_into_array(words, pos=0, var, "total_energyvar", total_energy);
  read_into_array(words, pos=0, avg, "weight", weight);
  read_into_array(words, pos=0, var, "weightvar", weight);
  read_into_array(words, pos=0, avg, "kinetic_energy", kinetic);
  read_into_array(words, pos=0, var, "kinetic_energyvar", kinetic);
  read_into_array(words, pos=0, avg, "potential_energy", potential);
  read_into_array(words, pos=0, var, "potential_energyvar", potential);
  read_into_array(words, pos=0, avg, "nonlocal_energy", nonlocal);
  read_into_array(words, pos=0, var, "nonlocal_energyvar", nonlocal);

  
  pos=0;
  vector <vector < string> > avgsecs;
  vector <string> avgwords;
  while(readsection(words, pos, avgwords, "average_generator")) { 
    avgsecs.push_back(avgwords);
  }
  if(avgsecs.size()%nwf !=0) { error("need average sections to equal a multiple of the number of wave functions"); }
  int navg=avgsecs.size()/nwf;
  avgrets.Resize(nwf, navg);
  int count=0;
  for(int w=0; w< nwf; w++) { 
    for(int i=0; i< avgrets.GetDim(1); i++) { 
      int n=avgsecs[i].size()-1;
      assert(n>=1);
      avgrets(w,i).vals.Resize(n);
      avgrets(w,i).type=avgsecs[count][0];
      for(int j=0; j< n; j++) { 
        avgrets(w,i).vals(j)=atof(avgsecs[count][j+1].c_str());
      }
      count++;
    }
  }
  
  
  int nautocorr=0;
  if(!readvalue(words, pos=0, nautocorr, "nautocorr"))
    error("No nautocorr in log file");
  if(nautocorr<0) nautocorr=0;

  autocorr.Resize(1,nautocorr);
  read_into_array(words, pos=0, autocorr, "autocorr_energy",0);

  aux_autocorr.Resize(naux, nautocorr);

  for(int i=0; i< naux; i++) {
    readvalue(auxtext[i], pos=0, aux_size(i), "aux_size");
    read_into_array(auxtext[i], pos=0, aux_energy, "aux_energy", i);
    read_into_array(auxtext[i], pos=0, aux_energyvar, "aux_energyvar",i);
    if(!read_into_array(auxtext[i], pos=0, aux_diff, "aux_diff", i)) {
      for(int w=0; w< n_cvg; w++) 
        aux_diff(i,w)=(aux_energy(i,w)-avg(total_energy,w))/aux_size(i);
    }
    Array1 <doublevar> diffvar(n_cvg);
    read_array(auxtext[i], pos=0, diffvar, "aux_diffvar");
    for(int w=0; w< n_cvg; w++) {
      aux_diff(i,w)*=aux_size(i);
      aux_diffvar(i,w)=diffvar(w)*aux_size(i);
    }

    read_into_array(auxtext[i], pos=0, aux_weight, "aux_weight", i);
    read_into_array(auxtext[i], pos=0, aux_weightvar, "aux_weightvar", i);
    read_into_array(auxtext[i], pos=0, aux_autocorr, "aux_autocorr", i);

  }  

  
}


//----------------------------------------------------------------------
void Properties_block::printBlockSummary(ostream & os) {
  int nwf=avg.GetDim(1);
  int naux=aux_energy.GetDim(0);

  using namespace Properties_types;
  int nprops=4;
  if(nwf>1) nprops++;

  assert(Properties_types::NUM_QUANTITIES >=5);
  os.setf(ios::showpoint);
  
  int field=os.precision()+10;
  int field1=4;


  os << setw(field1) << "  ";
  for(int p=0; p < nprops; p++) 
    os << setw(field) << Properties_types::names[p];
  os << endl;
  for(int w=0; w< nwf; w++) {
    os << setw(field1) << "##" << w;
    for(int p=0; p< nprops; p++)
      os << setw(field) << avg(p,w);
    os << "   (value)" << endl;
    os << setw(field1) << "&&" << w;
    for(int p=0; p< nprops; p++)
      os << setw(field) << sqrt(var(p,w));
    os << "   (sigma)" << endl;
  }


  if(naux > 0) { 
    os << "\nAuxillary functions : " << endl;
  
    os << setw(field1) << "  " 
       << setw(field) << "energy" 
       << setw(field) << "difference"
       << setw(field) << "weight" << endl;
  }
  int n_cvg=aux_diff.GetDim(1);
  for(int i=0; i< naux; i++) {
    
    for(int w=0; w< n_cvg; w++) {
      os << setw(field1) << "%%" << i << "-" << w
         << setw(field) << aux_energy(i,w) 
          << setw(field) << aux_diff(i,w)/aux_size(i) //(aux_energy(i,w)-avg(total_energy,w))/aux_size(i)
         << setw(field) << aux_weight(i,w) 
         << "   (value)" << endl;     
      
      os << setw(field1) << "$$" << i << "-" << w
         << setw(field) << sqrt(aux_energyvar(i,w))
         << setw(field) << sqrt(aux_diffvar(i,w))/aux_size(i)
         << setw(field) << sqrt(aux_weightvar(i,w))
         << "   (sigma)" << endl;     
    }
  }

  os.unsetf(ios::showpoint);
}

//----------------------------------------------------------------------

void Properties_block::reduceBlocks(Array1 <Properties_block> & blocks, 
                  int start, int end) {
  int nblocks=end-start;
  assert(nblocks >0);
  using namespace Properties_types;

  int nwf=blocks(0).avg.GetDim(1);
  int naux=blocks(0).aux_energy.GetDim(0);
  int n_cvg=blocks(0).aux_diff.GetDim(1);
  setSize(nwf, naux, n_cvg);

  avg=0.0; var=0.0;
  aux_energy=0; aux_weight=0; aux_energyvar=0; aux_weightvar=0;
  aux_diffvar=0; aux_size=0; aux_diff=0;
  int nautocorr=blocks(0).autocorr.GetDim(1);
  autocorr.Resize(1,nautocorr);
  autocorr=0.0;
  totweight=0.0;
  avgrets=blocks(start).avgrets;
  for(int w=0;w < nwf; w++) {
    for(int i=0; i< avgrets.GetDim(1); i++) { 
      avgrets(w,i).vals=0;
    }
  }
  
  for(int a=0; a< naux; a++) {
    aux_size(a)=blocks(start).aux_size(a);
  }

  for(int b=start; b< end; b++) {
    totweight+=blocks(b).totweight;
    
    for(int a=0; a< nautocorr; a++) {
      autocorr(0,a)+=blocks(b).autocorr(0,a)/nblocks;
    }
    
    for(int p=0; p < NUM_QUANTITIES; p++) {
      for(int w=0; w< nwf; w++) {
        avg(p,w)+=blocks(b).avg(p,w)/nblocks;
        var(p,w)+=blocks(b).var(p,w)*blocks(b).var(p,w);
      }
    }
    for(int w=0; w< nwf; w++) { 
      for(int i=0; i< avgrets.GetDim(1); i++) { 
        for(int j=0; j< avgrets(w,i).vals.GetDim(0); j++) { 
          avgrets(w,i).vals(j)+=blocks(b).avgrets(w,i).vals(j)/doublevar(nblocks);
        }
      }
    }

    for(int a=0; a< naux; a++) {
      for(int w=0; w< n_cvg; w++) {
        aux_energy(a,w)+=blocks(b).aux_energy(a,w)/nblocks;
        aux_weight(a,w)+=blocks(b).aux_weight(a,w)/nblocks;
        aux_energyvar(a,w)+=blocks(b).aux_energyvar(a,w)
          *blocks(b).aux_energyvar(a,w);
        aux_weightvar(a,w)+=blocks(b).aux_weightvar(a,w)
          *blocks(b).aux_weightvar(a,w);
        aux_diff(a,w)+=blocks(b).aux_diff(a,w)/nblocks;
        aux_diffvar(a,w)+=blocks(b).aux_diffvar(a,w)
          *blocks(b).aux_diffvar(a,w);
      }
    }
  }


  for(int p=0; p< NUM_QUANTITIES; p++) {
    for(int w=0; w< nwf; w++) {
      var(p,w)=sqrt(var(p,w))/nblocks;
    }
  }

  for(int a=0; a< naux; a++) {
    for(int w=0; w< n_cvg; w++) {
      aux_energyvar(a,w)=sqrt(aux_energyvar(a,w))/nblocks;
      aux_weightvar(a,w)=sqrt(aux_weightvar(a,w))/nblocks;
      aux_diffvar(a,w)=sqrt(aux_diffvar(a,w))/nblocks;
    }
  }


    
}
