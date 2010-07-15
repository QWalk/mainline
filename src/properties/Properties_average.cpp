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
#include "Properties_average.h"


void label_print(ostream & os,string label, Array1 <doublevar> & val,
                 Array1 <doublevar> & variance) {
  int pres=os.precision();
  for(int w=0; w< val.GetDim(0); w++) {
    os << label << w << setw(pres+5) << val(w) 
       << " dis " << setw(pres+5) << sqrt(variance(w)) 
       << endl;
  }
}


//######################################################################

void Properties_final_average::setSize(int nwf, int naux, int n_cvg) {

  using namespace Properties_types;

  avg.Resize(NUM_QUANTITIES, nwf);
  err.Resize(NUM_QUANTITIES, nwf);
  avgvar.Resize(NUM_QUANTITIES, nwf);
  avg=0.0; err=0.0;
  avgvar=0.0;

  diff_energy.Resize(nwf);
  diff_energyerr.Resize(nwf);
  aux_energy.Resize(naux, n_cvg);
  aux_energyerr.Resize(naux, n_cvg);
  aux_energyvar.Resize(naux, n_cvg);

  aux_diff.Resize(naux, n_cvg);
  aux_differr.Resize(naux, n_cvg);

  avgvar=0.0;
  aux_energy=0; aux_energyerr=0; aux_diff=0; aux_differr=0; aux_energyvar=0;  
 
}

//----------------------------------------------------------------------

void Properties_final_average::blockReduce(Array1 <Properties_block> & block_avg,
                                           int start_block, int end_block, 
                                           int find_equil) {
  
  using namespace Properties_types;
  assert(end_block > start_block);

  int nwf=block_avg(start_block).avg.GetDim(1);
  int naux=block_avg(start_block).aux_energy.GetDim(0);
  int n_cvg=block_avg(start_block).aux_energy.GetDim(1);
  setSize(nwf, naux, n_cvg);


  int autocorr_depth=block_avg(start_block).autocorr.GetDim(1);
  autocorr.Resize(nwf, autocorr_depth); autocorr=0;
  autocorrerr.Resize(nwf, autocorr_depth); autocorrerr=0;

  aux_autocorr.Resize(naux, autocorr_depth); aux_autocorr=0;
  aux_autocorrerr.Resize(naux, autocorr_depth); aux_autocorrerr=0;

  aux_size.Resize(naux);
  aux_size=block_avg(start_block).aux_size;

  int nblock=end_block-start_block;  
  avgavg=block_avg(start_block).avgrets;

  for(int w=0; w< nwf; w++)
    for(int i=0; i< avgavg.GetDim(1); i++) 
      avgavg(w,i).vals=0;
  avgerr=avgavg;
  
  //Stealing Jeff's algorithm for finding the throw-away blocks:
  //average over the last 4/5 of the blocks, and find the block
  //with energy within the standard deviation of that average..
  threw_out=0;
  //if we don't have more than five blocks, there aren't enough to determine 
  //equilibrium

  if(nblock > 5 && find_equil) { 
    doublevar tmpavg=0;
    int fkstart=start_block+nblock/5;
    for(int block=fkstart; block < end_block; block++) {
      tmpavg+=block_avg(block).avg(total_energy,0)/(end_block-fkstart);
    }
    doublevar tmpvar=0;
    for(int block=fkstart; block < end_block; block++) {
      tmpvar+=(block_avg(block).avg(total_energy,0)-tmpavg)
        *(block_avg(block).avg(total_energy,0)-tmpavg)/(end_block-fkstart);
    }
    tmpvar=sqrt(tmpvar);
    //cout << "tmpavg " << tmpavg << "  tmpvar " << tmpvar << endl;
    for(int block=start_block; block < end_block; block++) {
      if(fabs(block_avg(block).avg(total_energy,0)-tmpavg) < tmpvar) {
        threw_out=block-start_block;
        start_block=block;
        nblock=end_block-start_block;
        break;
      }
    }
  }

  totweight=0;
  for(int block=start_block; block < end_block; block++) {
    totweight+=block_avg(block).totweight;
  }
  

  for(int block=start_block; block < end_block; block++) {
    doublevar thisweight=block_avg(block).totweight/totweight;
    for(int w=0; w< nwf; w++) {
      for(int i=0; i< NUM_QUANTITIES; i++) {
        avg(i,w)+=block_avg(block).avg(i,w)*thisweight;
        avgvar(i,w)+=block_avg(block).var(i,w)*thisweight;
      }
      if(w==0) {
        for(int i=0; i< autocorr_depth; i++) {
          autocorr(w,i)+=block_avg(block).autocorr(w,i)*thisweight;
        }
      }
      for(int i=0; i< avgavg.GetDim(1); i++) { 
        for(int j=0; j< avgavg(w,i).vals.GetDim(0); j++) { 
          avgavg(w,i).vals(j)+=block_avg(block).avgrets(w,i).vals(j)*thisweight;
        }
      }
      if(w > 0) 
        diff_energy(w)+=(block_avg(block).avg(total_energy,w)
                         -block_avg(block).avg(total_energy,0))*thisweight;
      
    }
    
    
    for(int i=0; i< naux; i++) {
      for(int w=0; w< n_cvg; w++) {
        aux_diff(i,w)+=block_avg(block).aux_diff(i,w)*thisweight;
        aux_energy(i,w)+=block_avg(block).aux_energy(i,w)*thisweight;
        aux_energyvar(i,w)+=block_avg(block).aux_energyvar(i,w)*thisweight;
      }
    }
  }
  //int n2=nblock*nblock;
  for(int block=start_block; block < end_block; block++) {
    //doublevar wt=block_avg(block).totweight/totweight;
    doublevar thisweight=block_avg(block).totweight/(nblock*totweight);
    
    for(int w=0; w< nwf; w++) {
      for(int i=0; i< NUM_QUANTITIES; i++) {
        err(i,w)+=(block_avg(block).avg(i,w)-avg(i,w))
        *(block_avg(block).avg(i,w)-avg(i,w))*thisweight;
      }
      
      if(w==0) { 
        for(int i=0; i< autocorr_depth; i++) {
          autocorrerr(w,i)+=(block_avg(block).autocorr(w,i)-autocorr(w,i))
          *(block_avg(block).autocorr(w,i)-autocorr(w,i))*thisweight;
        }
      }
      for(int i=0; i< avgavg.GetDim(1); i++) { 
        for(int j=0; j< avgavg(w,i).vals.GetDim(0); j++) { 
          avgerr(w,i).vals(j)+=(block_avg(block).avgrets(w,i).vals(j)-avgavg(w,i).vals(j))
          *(block_avg(block).avgrets(w,i).vals(j)-avgavg(w,i).vals(j))*thisweight;
        }
      }
      
      if(w > 0) {
        doublevar ediff=block_avg(block).avg(total_energy,w)
        -block_avg(block).avg(total_energy,0);
        diff_energyerr(w)+=(ediff-diff_energy(w))*(ediff-diff_energy(w))*thisweight;
      }
    }
    
    
    for(int i=0; i< naux; i++) {
      for(int w=0; w< n_cvg; w++) {
        doublevar diff=block_avg(block).aux_diff(i,w);
        aux_differr(i,w)+=(diff-aux_diff(i,w))*(diff-aux_diff(i,w))*thisweight;
        aux_energyerr(i,w)=(block_avg(block).aux_energy(i,w)-aux_energy(i,w))
          *(block_avg(block).aux_energy(i,w)-aux_energy(i,w))*thisweight;
      }
    }
  }
  
  for(int w=0; w< nwf; w++) {
    for(int i=0; i< avgavg.GetDim(1); i++) { 
      for(int j=0; j< avgavg(w,i).vals.GetDim(0); j++) { 
        avgerr(w,i).vals(j)=sqrt(avgerr(w,i).vals(j));
      }
    }
  }
  
}

//----------------------------------------------------------------------

void Properties_final_average::averageReduce(Array1 <Properties_final_average> & final_avg,
					     int start, int end) {
  
  using namespace Properties_types;
  assert(end > start);

  int nwf=final_avg(start).avg.GetDim(1);
  int naux=final_avg(start).aux_energy.GetDim(0);
  int n_cvg=final_avg(start).aux_diff.GetDim(1);
  //cout << " n_cvg " << n_cvg << endl;
  //cout << " first avg " << endl;

  setSize(nwf, naux, n_cvg);

  avgvar=0.0;
  aux_energy=0; aux_energyerr=0; aux_diff=0; aux_differr=0;
  int autocorr_depth=final_avg(start).autocorr.GetDim(1);
  autocorr.Resize(nwf, autocorr_depth); autocorr=0;
  autocorrerr.Resize(nwf, autocorr_depth); autocorrerr=0;

  aux_autocorr.Resize(naux, autocorr_depth); aux_autocorr=0;
  aux_autocorrerr.Resize(naux, autocorr_depth); aux_autocorrerr=0;

  aux_size.Resize(naux);
  aux_size=final_avg(start).aux_size;

  int nblock=end-start;  
  
  //cout << "summing weights" << endl;
  
  totweight=0;
  for(int block=start; block < end; block++) {
    totweight+=final_avg(block).totweight;
    for(int w=0; w< nwf; w++) {
      for(int i=0; i< NUM_QUANTITIES; i++) {
        avg(i,w)+=final_avg(block).avg(i,w)/nblock;
        avgvar(i,w)+=final_avg(block).avgvar(i,w)/nblock;
      }
      
      for(int i=0; i< autocorr_depth; i++) {
        autocorr(w,i)+=final_avg(block).autocorr(w,i)/nblock;
      }
      if(w > 0) 
        diff_energy(w)+=(final_avg(block).avg(total_energy,w)
                         -final_avg(block).avg(total_energy,0))/nblock;
      
      
    }
    
    
    for(int i=0; i< naux; i++) {
      for(int w=0; w< n_cvg; w++) 
        aux_diff(i,w)+=final_avg(block).aux_diff(i,w)/nblock;
    }
  }
  
  
  int n2=nblock*nblock;
  for(int block=start; block < end; block++) {
    doublevar thisweight=1.0/n2;
    
    for(int w=0; w< nwf; w++) {
      for(int i=0; i< NUM_QUANTITIES; i++) {
        err(i,w)+=final_avg(block).err(i,w)*thisweight;
      }
      
      for(int i=0; i< autocorr_depth; i++) {
        autocorrerr(w,i)+=final_avg(block).autocorrerr(w,i)*thisweight;
      }
      if(w > 0) {
        diff_energyerr(w)+=final_avg(block).diff_energyerr(w)*thisweight;
      }
      
      
    }
    
    for(int i=0; i< naux; i++) {
      for(int w=0; w< n_cvg; w++) {
        aux_differr(i,w)+=final_avg(block).aux_differr(i,w)*thisweight;
      }
    }
  }
  
}

//--------------------------------------------------

#include "qmc_io.h"

//assuming that the imaginary and real parts are uncorrelated;
//probably not a great approximation, but can fix this later.
//We also return the bias from a nonlinear function.  This should
//scale like 1/M, where M is the number of sample points, but
//is definitely necessary to have under control

void complex_angle(dcomplex & val, dcomplex & err, 
		   double & angle, double & bias, double & angle_err) { 
  doublevar ratio=val.imag()/val.real();
  doublevar sigma_ratio=err.imag()*err.imag()/(val.real()*val.real())
    +val.imag()*val.imag()*err.real()*err.real()/(pow(val.real(),4));
  
  doublevar atander=1/(1+ratio*ratio); //first derivative of atan func
  doublevar atander2=-2*ratio*atander*atander; //second derivative
  
  bias=.5*atander2*sigma_ratio;
  angle_err=fabs(atander)*sqrt(sigma_ratio);
  
  angle=atan(fabs(ratio));
  if(val.real() <0 && val.imag() >=0)
    angle=pi-angle;
  else if(val.real() >=0 && val.imag() <0)
    angle=2*pi-angle;
  else if(val.real() <0 && val.imag() <0)
    angle+=pi;

}

//----------------------------------------------------------------------

void Properties_final_average::showSummary(ostream & os, Array1 <Average_generator*> avg_gen) { 
  int nwf=avg.GetDim(1);
  if(avg_gen.GetDim(0)%nwf != 0 ) error("Something wrong in the showSummary translation");
  int navg=avg_gen.GetDim(0)/nwf;
  Array2 <Average_generator*> avg2(nwf,navg);
  int count=0;
  for(int w=0; w< nwf; w++) {
    for(int i=0; i< navg; i++) avg2(w,i)=avg_gen(count++);
  }
  showSummary(os, avg2);
}
//----------------------------------------------------------------------

void Properties_final_average::showSummary(ostream & os, Array2 <Average_generator*>  avg_gen) {
  int nwf=avg.GetDim(1);
  int naux=aux_energy.GetDim(0);
  int autocorr_depth=autocorr.GetDim(1);
  os << "Threw out the first " << threw_out << " blocks in the averaging." 
     << endl;

  if(show_autocorr) {
    os << "Autocorrelation of energy " << endl;
    for(int i=0; i< autocorr_depth; i++) {
      os << "  en_ac_d" << i << "   " << autocorr(0,i) 
         << " +/- " << sqrt(autocorrerr(0,i)) << endl;
    }
    
    os << endl << endl;
  }
   
  using namespace  Properties_types;

  doublevar indep_points=sqrt(avgvar(0,0))/sqrt(err(0,0));
  int field=os.precision()+5;
  for(int w=0; w< nwf; w++) {
    for(int i=0; i < NUM_QUANTITIES; i++) {
      string tmp=names[i];
      append_number(tmp, w);
      os <<  setw(15)<< left  <<tmp << "    " 
         << right << setw(field) << avg(i,w)
         << " +/- " << setw(field) << sqrt(err(i,w))  
	 << " (sigma " << setw(field) << sqrt(avgvar(i,w)) << " ) " << endl;
    }

  }

  if(nwf > 1) {
    os << "Energy differences " << endl;
    for(int w =1; w < nwf; w++) 
      os << "  diff_" << w << "-0    " << diff_energy(w)
         << "  +/-   " << sqrt(diff_energyerr(w)) << endl;
  }
  
  
  if(naux > 0) 
    os << endl << endl << "Auxillary differences " << endl;
  assert(aux_size.GetDim(0)==naux);
  int n_cvg=aux_diff.GetDim(1);
  for(int i=0; i< naux; i++) {
    for(int w=0; w< n_cvg; w++) {
      os << "final_auxdiff" << i << "-" << w << "   "
      << aux_diff(i,w)/aux_size(i) << "  +/-   " 
      << sqrt(aux_differr(i,w))/aux_size(i) 
      << "   (sigma " << sqrt(aux_differr(i,w))*indep_points << " ) "
      << endl;
    }
  }
  os << endl;
  
 
  assert(avg_gen.GetDim(1)==avgavg.GetDim(1) && avgavg.GetDim(1)==avgerr.GetDim(1));
  for(int w=0; w< nwf; w++) { 
    for(int i=0; i< avg_gen.GetDim(1); i++) { 
      avg_gen(w,i)->write_summary(avgavg(w,i),avgerr(w,i), os);
    }
  }
  
  doublevar totpoints=indep_points*indep_points;
  os << "approximate number of independent points: " 
     << totpoints << endl;
  if(totpoints> totweight) {
    os << "Warning!  The estimated number of independent points is \n"
      " _greater_ than the points in the log file.  Error bars are \n"
      "probably unreliable!\n";
  }


}



void Properties_final_average::twoPointForces(Array2 <doublevar> & forces 
                                         ) {
  int nwf=avg.GetDim(1);
  int naux=aux_energy.GetDim(0);

  if(naux%2!=0) 
    error("there can't be two-point forces, since"
          " there's an odd number of differences");

  forces.Resize(naux/2, 2);
  for(int i=0; i < naux; i+=2) {
    assert(aux_size(i)==aux_size(i+1));
    forces(i/2,0)=(aux_diff(i,0)-aux_diff(i+1,0))/(aux_size(i)*2);
    forces(i/2,1)=sqrt( (aux_differr(i,0)+aux_differr(i+1,0))/4)/aux_size(i);
  }
}
