#include "Average_generator.h"
#include "qmc_io.h"
#include "gesqua.h"
#include "ulec.h"
#include "Average_density_matrix.h"
#include "Average_ekt.h"
#include "Average_region_fluctuation.h"
#include "Average_enmoment.h"
#include "Average_derivative_dm.h"
#include "Average_quadrupole.h"
#include "Average_so.h"
#include "Properties_point.h"
#include "Average_so.h"
#include "Average_ldots.h"
#include "jsontools.h"

//-----------------------------------------------------------------------------
int decide_averager(string & label, Average_generator *& avg) { 
  if(caseless_eq(label, "DIPOLE") ) 
    avg=new Average_dipole;
  else if(caseless_eq(label,"SK")) 
    avg=new Average_structure_factor;
  else if(caseless_eq(label,"MAGSK")) 
    avg=new Average_magnetic_structure_factor;
  
  else if(caseless_eq(label, "GR"))
    avg=new Average_twobody_correlation;
  else if(caseless_eq(label, "MANYBODY_POL"))
    avg=new Average_manybody_polarization;
  else if(caseless_eq(label, "OBDM"))
    avg=new Average_obdm;
  else if(caseless_eq(label, "EKT"))
    avg=new Average_ekt;
  else if(caseless_eq(label, "TBDM"))
    avg=new Average_tbdm;
  else if(caseless_eq(label, "LM"))
    avg=new Average_local_moments;
  else if(caseless_eq(label,"DENSITY_MOMENTS")) 
    avg=new Average_density_moments;
  else if(caseless_eq(label, "LINEAR_DER"))
    avg=new Average_linear_derivative;
  else if(caseless_eq(label, "LINEAR_DELTA_DER"))
    avg=new Average_linear_delta_derivative;
   else if(caseless_eq(label, "SPHERICAL_DENSITY"))
    avg=new Average_spherical_density;
  else if(caseless_eq(label, "SPHERICAL_DENSITY_GRID"))
    avg=new Average_spherical_density_grid;
  else if(caseless_eq(label, "LINE_DENSITY"))
    avg=new Average_line_density;
  else if(caseless_eq(label, "TBDM_BASIS"))
    avg=new Average_tbdm_basis;
  else if(caseless_eq(label,"REGION_FLUCTUATION"))
    avg=new Average_region_fluctuation;
  else if(caseless_eq(label,"REGION_DENSITY_MATRIX"))
    avg=new Average_region_density_matrix;
  else if(caseless_eq(label,"WF_PARMDERIV")) 
    avg=new Average_wf_parmderivs;
  else if(caseless_eq(label,"ENMOMENT"))
    avg=new Average_enmoment;
  else if(caseless_eq(label,"FOURIER_DENSITY"))
    avg=new Average_fourier_density;
  else if(caseless_eq(label,"AVERAGE_DERIVATIVE_DM"))
    avg=new Average_derivative_dm;
  else if(caseless_eq(label,"QUADRUPOLE"))
    avg=new Average_quadrupole;
  else if(caseless_eq(label,"JOINT_PBC"))
    avg=new Average_joint_pbc;
  else if(caseless_eq(label,"PSP_SO"))
    avg=new Average_so;  
  else if(caseless_eq(label,"LDOTS"))
    avg=new Average_ldots;
  else 
    error("Didn't understand ", label, " in Average_generator.");
  
  return 1;
}


//-----------------------------------------------------------------------------


int allocate(vector<string> & words, System * sys, Wavefunction_data * wfdata, Average_generator * & avg) { 
  decide_averager(words[0], avg);
  avg->read(sys, wfdata, words);
  return 1;
}
//-----------------------------------------------------------------------------


int allocate(vector<string> & words, Average_generator * & avg) { 
  decide_averager(words[0], avg);
  avg->read(words);
  return 1;
}
//############################################################################
//Average_dipole


void Average_dipole::read(System * sys, Wavefunction_data * wfdata, vector <string> & words) { 
}

void Average_dipole::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                              System * sys, Sample_point * sample, Average_return & avg ) { 
  avg.type="dipole";
  int ndim=3;
  avg.vals.Resize(ndim);
  avg.vals=0.0;
  int nelectrons=sample->electronSize();
  Array1 <doublevar> pos(ndim);
  for(int e=0; e< nelectrons; e++) {
    sample->getElectronPos(e,pos);
    for(int d=0; d< ndim; d++) {
      avg.vals(d)-=pos(d);
    }
  }
  int nions=sample->ionSize();
  for(int at=0; at < nions; at++) {
    sample->getIonPos(at,pos);
    doublevar charge=sample->getIonCharge(at);
    for(int d=0; d< ndim; d++) {
      avg.vals(d)+=charge*pos(d);
    }
  }
}

void Average_dipole::write_init(string & indent, ostream & os) { 
  os << indent << "dipole \n";
}

void Average_dipole::read(vector <string> & words) { 
  
}

void Average_dipole::write_summary(Average_return & avg, Average_return & err, ostream & os) {
  int ndim=avg.vals.GetDim(0);
  assert(ndim <= err.vals.GetDim(0));
  //Could put this in Debye if we want to be nice.
  os << "Dipole moment (a.u.) \n";
  for(int d=0; d< ndim; d++) { 
    if(d==0) os << "x ";
    else if(d==1) os << "y ";
    else if(d==2) os << "z ";
    os << avg.vals(d) << " +/- " << err.vals(d) << endl;
  }
    
}
//############################################################################


void Average_structure_factor::read(System * sys, Wavefunction_data * wfdata, vector <string> & words) {
  unsigned int pos=0;
  int np_side;
  if(!readvalue(words, pos=0, np_side, "NGRID")) 
    np_side=5;
  
  Array2 <doublevar> gvec;
  vector <string> gvec_sec;
  if(readsection(words, pos=0, gvec_sec, "GVEC")) {
    int count=0;
    gvec.Resize(3,3);
    for(int i=0; i< 3; i++) { 
      for(int j=0; j< 3; j++) { 
        gvec(i,j)=atof(gvec_sec[count++].c_str());
      }
    }
  }
  else { 
    if(!sys->getRecipLattice(gvec)) 
      error("You don't have a periodic cell and you haven't specified GVEC for S(k)");
  }
  
  int direction;
  if(readvalue(words, pos=0, direction, "DIRECTION")) { 
    npoints=np_side;
    kpts.Resize(npoints, 3);
    kpts=0;
    for(int ix=0; ix < np_side; ix++) { 
      for(int i=0; i< 3; i++) 
        kpts(ix,i)=2*pi*gvec(direction,i)*ix;
    }
  }
  else { 
    npoints=np_side*np_side*np_side; 

    kpts.Resize(npoints, 3);
    kpts=0;
    int c=0;
    for(int ix=0; ix < np_side; ix++) {
      for(int iy=0; iy < np_side; iy++) {
        for(int iz=0; iz < np_side; iz++) {
          for(int i=0; i< 3; i++)
            kpts(c,i)+=2*pi*(gvec(0,i)*ix+gvec(1,i)*iy+gvec(2,i)*iz);
          c++;
        }
      }
    }
  }
}  

//-----------------------------------------------------------------------------

void Average_structure_factor::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                                        System * sys, Sample_point * sample, Average_return & avg) {
  avg.type="sk";
  int nelectrons=sample->electronSize();
  
  Array1 <doublevar> pos(3);
  avg.vals.Resize(npoints);
  avg.vals=0;
  for(int p=0; p < npoints; p++) { 
    doublevar sum_cos=0, sum_sin=0;
    
    for(int e=0; e< nelectrons; e++) {
      sample->getElectronPos(e,pos);
      doublevar dot=0;
      for(int d=0; d< 3; d++) dot+=pos(d)*kpts(p,d);
      sum_cos+=cos(dot);
      sum_sin+=sin(dot);
    }
    avg.vals(p)+=(sum_cos*sum_cos+sum_sin*sum_sin)/nelectrons;
    
  }
  
}

//-----------------------------------------------------------------------------

void Average_structure_factor::write_init(string & indent, ostream & os) { 
  os << indent << "sk\n";
  for(int i=0; i< npoints; i++) { 
    os << indent << "kpoint { " << kpts(i,0) << " " << kpts(i,1) << " " 
    << kpts(i,2) << " } " <<endl;
  }
}
//-----------------------------------------------------------------------------

void Average_structure_factor::read(vector <string> & words) { 
  vector <vector <string> > kpttext;
  vector <string> tmp;
  unsigned int pos=0;
  while(readsection(words, pos, tmp, "KPOINT")) kpttext.push_back(tmp);
  npoints=kpttext.size();
  kpts.Resize(npoints,3);
  for(int i=0; i< npoints; i++) { 
    for(int d=0; d< 3; d++) { 
      kpts(i,d)=atof(kpttext[i][d].c_str());
    }
  }
}
//-----------------------------------------------------------------------------

void Average_structure_factor::write_summary(Average_return & avg, Average_return & err, ostream & os) { 
  os << "Structure factor \n";
  os << "    k  s(k) s(k)err  kx  ky  kz" << endl;
  assert(avg.vals.GetDim(0) >=npoints);
  assert(err.vals.GetDim(0) >=npoints);
  
  for(int i=0; i< npoints; i++) {
    doublevar sk=avg.vals(i);
    doublevar skerr=err.vals(i);
    doublevar r=0;
    for(int d=0; d< 3; d++) r+=kpts(i,d)*kpts(i,d);
    r=sqrt(r);
    os << "sk_out " <<  r << "   " << sk << "   " << skerr 
      << "  " << kpts(i,0) << "   " << kpts(i,1) << "   " << kpts(i,2) << endl;
  }
  
}

//--------------------------------------------------------------------------------
void Average_structure_factor::jsonOutput(Average_return & avg, Average_return & err, ostream & os) {
    
  assert(avg.vals.GetDim(0) >=npoints);
  assert(err.vals.GetDim(0) >=npoints);
    
  os << "\"" << avg.type << "\":{" << endl;
  os << "\"k\":[" << endl;
  for(int i=0; i< npoints; i++) {
    doublevar r=0;
    for(int d=0; d< 3; d++) r+=kpts(i,d)*kpts(i,d);
    r=sqrt(r);
    os << r;
    if(i<npoints-1) os << ",";
  }
  os << endl;
  os << "]," << endl;
    
  os << "\"value\":[" << endl;
  for(int i=0; i< npoints; i++) {
    doublevar sk=avg.vals(i);
    os << sk;
    if( i<npoints-1) os << ",";
  }
  os << endl;
  os << "]," << endl;
    
  os << "\"error\":[" << endl;
  for(int i=0; i< npoints; i++) {
    doublevar skerr=err.vals(i);
    os << skerr;
    if(i< npoints-1) os << ",";
  }
  os << endl;
  os << "]," << endl;
    
  os << "\"kx\":[" << endl;
  for(int i=0; i< npoints; i++) {
    os << kpts(i,0);
    if(i< npoints-1) os << ",";
  }
  os << endl;
  os << "]," << endl;
    
  os << "\"ky\":[" << endl;
  for(int i=0; i< npoints; i++) {
    os << kpts(i,1);
    if(i< npoints-1) os << ",";
  }
  os << endl;
  os << "]," << endl;
    
  os << "\"kz\":[" << endl;
  for(int i=0; i< npoints; i++) {
    os << kpts(i,2);
    if(i< npoints-1) os << ",";
  }
  os << endl;
  os << "]" << endl;
    
  os << "}" << endl;
}

//############################################################################

void Average_magnetic_structure_factor::read(System * sys, Wavefunction_data * wfdata, vector <string> & words) {
  unsigned int pos=0;
  int np_side;
  if(!readvalue(words, pos=0, np_side, "NGRID")) 
    np_side=5;
  
  Array2 <doublevar> gvec;
  vector <string> gvec_sec;
  if(readsection(words, pos=0, gvec_sec, "GVEC")) {
    int count=0;
    gvec.Resize(3,3);
    for(int i=0; i< 3; i++) { 
      for(int j=0; j< 3; j++) { 
        gvec(i,j)=atof(gvec_sec[count++].c_str());
      }
    }
  }
  else { 
    if(!sys->getRecipLattice(gvec)) 
      error("You don't have a periodic cell and you haven't specified GVEC for S(k)");
  }
  
  int direction;
  if(readvalue(words, pos=0, direction, "DIRECTION")) { 
    npoints=np_side;
    kpts.Resize(npoints, 3);
    kpts=0;
    for(int ix=0; ix < np_side; ix++) { 
      for(int i=0; i< 3; i++) 
        kpts(ix,i)=2*pi*gvec(direction,i)*ix;
    }
  }
  else { 
    npoints=np_side*np_side*np_side; 

    kpts.Resize(npoints, 3);
    kpts=0;
    int c=0;
    for(int ix=0; ix < np_side; ix++) {
      for(int iy=0; iy < np_side; iy++) {
        for(int iz=0; iz < np_side; iz++) {
          for(int i=0; i< 3; i++)
            kpts(c,i)+=2*pi*(gvec(0,i)*ix+gvec(1,i)*iy+gvec(2,i)*iz);
          c++;
        }
      }
    }
  }
}  

//-----------------------------------------------------------------------------

void Average_magnetic_structure_factor::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                                        System * sys, Sample_point * sample, Average_return & avg) {
  avg.type="magsk";
  int nelectrons=sample->electronSize();
  
  Array1 <doublevar> pos1(3),pos2(3);
  avg.vals.Resize(npoints*2);
  avg.vals=0;
  int nup=sys->nelectrons(0);
  for(int p=0; p < npoints; p++) { 
    doublevar sum_cos=0, sum_sin=0;
    
    for(int e1=0; e1< nelectrons; e1++) {
      sample->getElectronPos(e1,pos1);
      int s1=e1 < nup?1:-1;
      for(int e2=0; e2 < nelectrons; e2++) {
        sample->getElectronPos(e2,pos2);
        doublevar dot=0;
        for(int d=0; d< 3; d++) dot+=(pos1(d)-pos2(d))*kpts(p,d);
        int s2=e2 < nup?1:-1;
        sum_cos+=cos(dot)*s1*s2;
        sum_sin+=sin(dot)*s1*s2;
      }
    }
    avg.vals(2*p)+=sum_cos;
    avg.vals(2*p+1)+=sum_sin;
    
  }
  
}

//-----------------------------------------------------------------------------

void Average_magnetic_structure_factor::write_init(string & indent, ostream & os) { 
  os << indent << "magsk\n";
  for(int i=0; i< npoints; i++) { 
    os << indent << "kpoint { " << kpts(i,0) << " " << kpts(i,1) << " " 
    << kpts(i,2) << " } " <<endl;
  }
}
//-----------------------------------------------------------------------------

void Average_magnetic_structure_factor::read(vector <string> & words) { 
  vector <vector <string> > kpttext;
  vector <string> tmp;
  unsigned int pos=0;
  while(readsection(words, pos, tmp, "KPOINT")) kpttext.push_back(tmp);
  npoints=kpttext.size();
  kpts.Resize(npoints,3);
  for(int i=0; i< npoints; i++) { 
    for(int d=0; d< 3; d++) { 
      kpts(i,d)=atof(kpttext[i][d].c_str());
    }
  }
}
//-----------------------------------------------------------------------------

void Average_magnetic_structure_factor::write_summary(Average_return & avg, Average_return & err, ostream & os) { 
  os << "Structure factor \n";
  os << "    k  s(k) s(k)err  kx  ky  kz" << endl;
  assert(avg.vals.GetDim(0) >=2*npoints);
  assert(err.vals.GetDim(0) >=2*npoints);
  
  for(int i=0; i< npoints; i++) {
    doublevar skr=avg.vals(2*i), ski=err.vals(2*i+1);
    doublevar skrerr=err.vals(2*i),skierr=err.vals(2*i+1);
    doublevar r=0;
    for(int d=0; d< 3; d++) r+=kpts(i,d)*kpts(i,d);
    r=sqrt(r);
    os << "sk_out " <<  r << "   " << skr << "   " << skrerr 
       << "   " << ski << "   " << skierr 
      << "  " << kpts(i,0) << "   " << kpts(i,1) << "   " << kpts(i,2) << endl;
  }
  
}

//--------------------------------------------------------------------------------
void Average_magnetic_structure_factor::jsonOutput(Average_return & avg, Average_return & err, ostream & os) {
    
  assert(avg.vals.GetDim(0) >=npoints);
  assert(err.vals.GetDim(0) >=npoints);
    
  os << "\"" << avg.type << "\":{" << endl;
  os << "\"k\":";
  jsonarray(os,kpts);
  os << ",\n";

  os << "\"value\":";jsonarray(os,avg.vals);
  os << ",\n";
  os << "\"error\":"; jsonarray(os,err.vals); 
    
  os << "}" << endl;
}

//############################################################################


void Average_fourier_density::read(System * sys, Wavefunction_data * wfdata, vector <string> & words) {
  unsigned int pos=0;
  int np_side;
  if(!readvalue(words, pos=0, np_side, "NGRID")) 
    np_side=5;
  
  Array2 <doublevar> gvec;
  vector <string> gvec_sec;
  if(readsection(words, pos=0, gvec_sec, "GVEC")) {
    int count=0;
    gvec.Resize(3,3);
    for(int i=0; i< 3; i++) { 
      for(int j=0; j< 3; j++) { 
        gvec(i,j)=atof(gvec_sec[count++].c_str());
      }
    }
  }
  else { 
    if(!sys->getRecipLattice(gvec)) 
      error("You don't have a periodic cell and you haven't specified GVEC for Fourier_density");
  }
  
  npoints=np_side*np_side*np_side; 
  kpts.Resize(npoints, 3);
  kpts=0;
  int c=0;
  for(int ix=0; ix < np_side; ix++) {
    for(int iy=0; iy < np_side; iy++) {
      for(int iz=0; iz < np_side; iz++) {
        for(int i=0; i< 3; i++)
          kpts(c,i)+=2*pi*(gvec(0,i)*ix+gvec(1,i)*iy+gvec(2,i)*iz);
        c++;
      }
    }
  }
}  

//-----------------------------------------------------------------------------

void Average_fourier_density::write_init(string & indent, ostream & os) { 
  os << indent << "fourier_density\n";
  for(int i=0; i< npoints; i++) { 
    os << indent << "kpoint { " << kpts(i,0) << " " << kpts(i,1) << " " 
    << kpts(i,2) << " } " <<endl;
  }
}
//-----------------------------------------------------------------------------

void Average_fourier_density::read(vector <string> & words) { 
  vector <vector <string> > kpttext;
  vector <string> tmp;
  unsigned int pos=0;
  while(readsection(words, pos, tmp, "KPOINT")) kpttext.push_back(tmp);
  npoints=kpttext.size();
  kpts.Resize(npoints,3);
  for(int i=0; i< npoints; i++) { 
    for(int d=0; d< 3; d++) { 
      kpts(i,d)=atof(kpttext[i][d].c_str());
    }
  }
}
//-----------------------------------------------------------------------------


void Average_fourier_density::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                                        System * sys, Sample_point * sample, Average_return & avg) {
  avg.type="fourier_density";
  int nelectrons=sample->electronSize();
  int start[2],end[2],ne[2];
  ne[0]=sys->nelectrons(0);
  ne[1]=sys->nelectrons(1);
  start[0]=0; start[1]=ne[0];
  end[0]=ne[0]; end[1]=ne[0]+ne[1];
  int nspin=2;
  int stride=10;
  Array1 <doublevar> pos(3),pos2(3);
  avg.vals.Resize(stride*npoints); //(upcos,upsin,downcos,downsin,upupcos,upupsin)
  avg.vals=0;

  for(int p=0; p < npoints; p++) {
    //----------Here we do the total density 
    for(int s=0; s<nspin; s++) {  
      doublevar sum_cos=0, sum_sin=0;
      for(int e=start[s]; e< end[s]; e++) {
        sample->getElectronPos(e,pos);
        doublevar dot=0;
        for(int d=0; d< 3; d++) dot+=pos(d)*kpts(p,d);
        sum_cos+=cos(dot);
        sum_sin+=sin(dot);
      }
      avg.vals(stride*p+s*2)+=sum_cos/ne[s];
      avg.vals(stride*p+s*2+1)+=sum_sin/ne[s];
    }
    //------And here we do the pair density
    
    for(int s1=0; s1 < nspin; s1++) { 
      for(int e1=start[s1]; e1 < end[s1]; e1++) { 
        sample->getElectronPos(e1,pos);
        for(int s2=0; s2 < nspin; s2++) { 
          doublevar sum_cos=0, sum_sin=0;
          for(int e2=start[s2]; e2 < end[s2]; e2++) { 
            sample->getElectronPos(e2,pos2);
            doublevar dot=0;
            for(int d=0; d< 3; d++) dot+=kpts(p,d)*(pos2(d)-pos(d));
            sum_cos+=cos(dot);
            sum_sin+=sin(dot);
          }
          int offset=stride*p+nspin*2+2*(s1+s2);
          avg.vals(offset)+=sum_cos/ne[s2];
          avg.vals(offset+1)+=sum_sin/ne[s2];
        }
      }
    }
  }

  
}

void Average_fourier_density::write_summary(Average_return & avg, Average_return & err, ostream & os) { 
  //To get an output in table form, do 
  //gosling log  | grep fourier_out | awk '{$1=""; print $0}'
  int stride=avg.vals.GetDim(0)/npoints;
  os << "fourier_out kx ky kz upcos upcoserr upsin upsinerr downcos downcoserr downsin downsinerr uucos uucoserr uusin uusinerr udcos udcoserr udsin udsinerr  ddcos ddcoserr ddsin ddsinerr \n";
  for(int p=0; p < npoints; p++) { 
     os << "fourier_out ";
     for(int d=0; d< 3; d++) os << " " << kpts(p,d);
     for(int i=0; i< stride; i++) 
       os << " " << avg.vals(stride*p+i) << " " << err.vals(stride*p+i);
     os << endl;
  }

}

//############################################################################

void Average_twobody_correlation::read(System * sys, Wavefunction_data * wfdata, 
                                       vector <string> & words) { 
  doublevar range;
  unsigned int pos=0;
  if(!readvalue(words,pos=0,resolution, "RESOLUTION"))
    resolution=0.1;
  if(!readvalue(words,pos=0,range,"RANGE"))
    range=10.0;
  if(readvalue(words, pos=0, direction, "DIRECTION")) { 
    direction+=2;
  }
  else
    direction=0;
  
  npoints=int(range/resolution)+1;  
  
}

//-----------------------------------------------------------------------------

void Average_twobody_correlation::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                                        System * sys, Sample_point * sample, Average_return & avg) {
  avg.type="gr";
  int nelectrons=sample->electronSize();
  avg.vals.Resize(2*npoints); //for like and unlike
  avg.vals=0;
  int nup=sys->nelectrons(0);
  sample->updateEEDist();
  Array1 <doublevar> dist(5);
  for(int e=0; e< nup; e++) {
    for(int e2=e+1; e2 < nup; e2++) {
      sample->getEEDist(e,e2,dist);
      int place=int(fabs(dist(direction))/resolution);
      if(place < npoints) { 
        avg.vals(place)+=1;
      }
    }
    for(int e2=nup; e2 < nelectrons; e2++) { 
      sample->getEEDist(e,e2,dist);
      int place=int(fabs(dist(direction))/resolution);
      if(place < npoints) { 
        avg.vals(npoints+place)+=1;
      }
    }
  }
  
  for(int e=nup; e < nelectrons; e++) { 
    for(int e2=e+1; e2< nelectrons; e2++) { 
      sample->getEEDist(e,e2,dist);
      int place=int(fabs(dist(direction))/resolution);
      if(place < npoints) { 
        avg.vals(place) += 1.0;
      }
    }
  }
  
  
  
}

//-----------------------------------------------------------------------------

void Average_twobody_correlation::write_init(string & indent, ostream & os) { 
  os << indent << "gr\n";
  os << indent << "npoints " << npoints << endl;
  os << indent << "resolution " << resolution << endl;
  os << indent << "direction " << direction << endl;
}
//-----------------------------------------------------------------------------

void Average_twobody_correlation::read(vector <string> & words) { 
  unsigned int pos=0;
  readvalue(words, pos=0, resolution, "resolution");
  readvalue(words, pos=0, npoints, "npoints");
  readvalue(words, pos=0, direction, "direction");
}
//-----------------------------------------------------------------------------

void Average_twobody_correlation::write_summary(Average_return & avg, Average_return & err, ostream & os) { 
  os << "Electron correlation for like and unlike spins\n";
  os << "    r  g(r) sigma(g(r))   g(r)  sigma(g(r))" << endl;
  assert(avg.vals.GetDim(0) >=2*npoints);
  assert(err.vals.GetDim(0) >=2*npoints);
  
  
  for(int i=0; i< npoints; i++) {
    os << "gr_out " << i*resolution << " " << avg.vals(i) << " " << err.vals(i) 
    << "  " << avg.vals(i+npoints) << "  " << err.vals(i+npoints) << endl;
  }
  
}
//-----------------------------------------------------------------------------

void Average_twobody_correlation::jsonOutput(Average_return & avg,Average_return & err, ostream & os) {
  assert(avg.vals.GetDim(0) >=2*npoints);
  assert(err.vals.GetDim(0) >=2*npoints);

  os << "\"" << avg.type << "\":{" << endl;
  os << "\"r\":[" << endl;
  for(int i=0; i< npoints; i++) {
    os << i*resolution;
    if(i<npoints-1) os << ",";
  }
  os << endl;
  os << "]," << endl;
  
  os << "\"like\":[" << endl;
  for(int i=0; i< npoints; i++) {
    os << avg.vals(i);
    if( i<npoints-1) os << ",";
  }
  os << endl;
  os << "]," << endl;

  os << "\"unlike\":[" << endl;
  for(int i=0; i< npoints; i++) {
    os << avg.vals(i+npoints);
    if( i<npoints-1) os << ",";
  }
  os << endl;
  os << "]," << endl;

  
  os << "\"like_err\":[" << endl;
  for(int i=0; i< npoints; i++) {
    os << err.vals(i);
    if( i<npoints-1) os << ",";
  }
  os << endl;
  os << "]," << endl;
  
  os << "\"unlike_err\":[" << endl;
  for(int i=0; i< npoints; i++) {
    os << err.vals(i+npoints);
    if( i<npoints-1) os << ",";
  }
  os << endl;
  os << "]" << endl;

  os << "}" << endl;

}



//############################################################################


void Average_manybody_polarization::read(System *sys, Wavefunction_data * wfdata, vector <string> & words) { 
  gvec.Resize(3,3);
  if(!sys->getRecipLattice(gvec) ) { 
    error("The manybody polarization operator works only for periodic systems!");
  }
}

//-----------------------------------------------------------------------------
void Average_manybody_polarization::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                                             System * sys, Sample_point * sample, Average_return & avg) {
  avg.type="manybody_pol";
  int nelectrons=sample->electronSize();
  avg.vals.Resize(6);
  Array1 <doublevar> sum(3,0.0);
  Array1 <doublevar> pos(3);
  for(int e=0; e< nelectrons; e++) { 
    sample->getElectronPos(e,pos);
    for(int i=0; i< 3; i++) { 
      for(int d=0; d< 3; d++) { 
        sum(i)-=gvec(i,d)*pos(d);
      }
    }
  }
  int nions=sample->ionSize();  
  for(int at=0; at < nions; at++) {
    sample->getIonPos(at,pos);
    doublevar charge=sample->getIonCharge(at);
    for(int i=0; i< 3; i++) 
      for(int d=0; d< 3; d++) 
        sum(i)+=charge*gvec(i,d)*pos(d);
  }
  
  for(int i=0; i< 3; i++) { 
    avg.vals(2*i)=cos(2*pi*sum(i));
    avg.vals(2*i+1)=sin(2*pi*sum(i));
  }
}
//-----------------------------------------------------------------------------

void Average_manybody_polarization::write_init(string & indent, ostream & os) { 
  os << indent << "manybody_pol" << endl;
  os << indent << "gvec { ";
  for(int i=0; i< 3; i++) 
    for(int j=0; j< 3; j++) os << gvec(i,j) << "  ";
  os << "}\n";
}
//-----------------------------------------------------------------------------
void Average_manybody_polarization::read(vector <string> & words) { 
  unsigned int pos=0;
  vector <string> gvec_sec;
  readsection(words, pos=0, gvec_sec, "gvec");
  int count=0;
  gvec.Resize(3,3);
  for(int i=0; i< 3; i++)
    for(int j=0; j< 3; j++) gvec(i,j)=atof(gvec_sec[count++].c_str());
  
}
//-----------------------------------------------------------------------------


void Average_manybody_polarization::write_summary(Average_return & avg, Average_return & err, ostream & os) { 
  os << "Manybody polarization operator " << endl;
  for(int i=0; i< 3; i++) { 
    os << avg.vals(2*i) << " + i " << avg.vals(2*i+1) << " +/- " << err.vals(2*i) << " + i " << err.vals(2*i+1) << endl;
  }
}

//--------------------------------------------------------------------------------
void Average_manybody_polarization::jsonOutput(Average_return & avg, Average_return & err, ostream & os) {
  os << "\"" << avg.type << "\":{" << endl;
  os <<"\"value\":[";
  for(int i=0; i< 3; i++) {
    os << "[" << avg.vals(2*i) << "," << avg.vals(2*i+1) <<"]";
    if(i<2) os <<",";
  }
  os <<"]," << endl;
    
  os <<"\"error\":[";
  for(int i=0; i< 3; i++) {
    os << "[" << err.vals(2*i) << "," << err.vals(2*i+1) <<"]";
    if(i<2) os <<",";
  }
  os <<"]" << endl;
  os << "}" << endl;
}

//############################################################################

/*!
Auxiliary function returning smallest "height" of the simulation cell
in a periodic system. Used in one- and two-body density matrices to get
a sensible upper cut-off of distances for which elements of the density
matrices are evaluated.
*/
doublevar smallest_height(Array2 <doublevar> & latVec) {

  int ndim=3;
  Array2 <doublevar> crossProduct(ndim,ndim);
  crossProduct(0,0)=(latVec(1,1)*latVec(2,2)-latVec(1,2)*latVec(2,1));
  crossProduct(0,1)=(latVec(1,2)*latVec(2,0)-latVec(1,0)*latVec(2,2));
  crossProduct(0,2)=(latVec(1,0)*latVec(2,1)-latVec(1,1)*latVec(2,0));

  crossProduct(1,0)=(latVec(2,1)*latVec(0,2)-latVec(2,2)*latVec(0,1));
  crossProduct(1,1)=(latVec(2,2)*latVec(0,0)-latVec(2,0)*latVec(0,2));
  crossProduct(1,2)=(latVec(2,0)*latVec(0,1)-latVec(2,1)*latVec(0,0));

  crossProduct(2,0)=(latVec(0,1)*latVec(1,2)-latVec(0,2)*latVec(1,1));
  crossProduct(2,1)=(latVec(0,2)*latVec(1,0)-latVec(0,0)*latVec(1,2));
  crossProduct(2,2)=(latVec(0,0)*latVec(1,1)-latVec(0,1)*latVec(1,0));
  
  doublevar smallestheight=1e99;
  for(int i=0; i< ndim; i++) {
    doublevar tempheight=0;
    doublevar length=0;
    for(int j=0; j< ndim; j++) {
      tempheight+=crossProduct(i,j)*latVec(i,j);
      length+=crossProduct(i,j)*crossProduct(i,j);
    }
    tempheight=fabs(tempheight)/sqrt(length);
    if(tempheight < smallestheight ) smallestheight=tempheight;
  }

  return smallestheight;

}

//-----------------------------------------------------------------------------

void Average_obdm::read(System * sys, Wavefunction_data * wfdata,
			vector <string> & words) {

  single_write(cout, "One-body density matrix will be calculated.\n");

  unsigned int pos=0;
  if(!readvalue(words, pos=0, npoints, "NGRID")) 
    npoints=5;          
  
  // maximum distance for OBDM evaluation

  pos=0;
  doublevar cutoff;
  if(!readvalue(words, pos=0, cutoff, "CUTOFF")) {
    Array2 <doublevar> latVec(3,3);
    Array1 <doublevar> origin(3);
    if(!sys->getBounds(latVec,origin))
      error("In non-periodic systems, CUTOFF has to be given in the OBDM section.");
    cutoff=smallest_height(latVec)/2;
  }
  dR=cutoff/npoints;
  
  // angles and weights for Gaussian quadrature (spherical averaging)

  pos=0;
  if(!readvalue(words, pos=0, np_aver, "AIP")) np_aver=1;
  wt.Resize(np_aver);
  ptc.Resize(np_aver,3);             //!< cartesian coordinates of int. points
  Array1 <doublevar> x(np_aver), y(np_aver), z(np_aver);

  switch (np_aver) {
  case 1:
    // no spherical averaging
    wt=1;
    ptc(0,0)=1.0;
    ptc(0,1)=0.0;
    ptc(0,2)=0.0;
    break;
  default:
    gesqua(np_aver,x,y,z,wt);
    for (int i=0; i<np_aver; i++) {
      ptc(i,0)=x(i);
      ptc(i,1)=y(i);
      ptc(i,2)=z(i);
    }
  }
}  

//-----------------------------------------------------------------------------

void Average_obdm::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
			    System * sys, Sample_point * sample,
			    Average_return & avg) {
  avg.type="obdm";
  avg.vals.Resize(npoints);
  avg.vals=0;

  int nwf=wf->nfunc();
  Wf_return wfval_new(nwf,2);   // this structure I just don't understand
  Wf_return wfval_old(nwf,2);
  Storage_container wfStore;
  Array1 <doublevar> oldpos(3), newpos(3), transl(3);

  Array1 <doublevar> x(3), y(3), z(3);
  Array2 <doublevar> pt(np_aver,3);
  generate_random_rotation(x,y,z);
  for (int i=0; i<np_aver; i++) {
    pt(i,0)=ptc(i,0)*x(0)+ptc(i,1)*y(0)+ptc(i,2)*z(0);
    pt(i,1)=ptc(i,0)*x(1)+ptc(i,1)*y(1)+ptc(i,2)*z(1);
    pt(i,2)=ptc(i,0)*x(2)+ptc(i,1)*y(2)+ptc(i,2)*z(2);
  }

  wf->updateVal(wfdata, sample);
  wfStore.initialize(sample, wf);
  wfStore.saveUpdate(sample,wf,0);
  wf->getVal(wfdata, 0, wfval_old);
  sample->getElectronPos(0,oldpos);

  for( int i=0; i<npoints; i++) {
    for ( int k=0; k<np_aver; k++) {
      // shift makes sense only up to half the lattice vector, we are hitting
      // periodicity for larger distances
      transl(0)=(i+1)*dR*pt(k,0);
      transl(1)=(i+1)*dR*pt(k,1);
      transl(2)=(i+1)*dR*pt(k,2);
      sample->translateElectron(0,transl);
      // TODO: check k-points /= Gamma
      wf->updateVal(wfdata, sample);
      wf->getVal(wfdata, 0, wfval_new);
      // how to handle the possibility of wf_old=0? Can it happen?
      avg.vals(i)+=wt(k)*wfval_new.sign(0)*wfval_old.sign(0)
        *exp(wfval_new.amp(0,0)-wfval_old.amp(0,0));
      sample->setElectronPos(0,oldpos);
      wfStore.restoreUpdate(sample,0);
    }
  }
  
  wfStore.restoreUpdate(sample,wf,0);

}

//-----------------------------------------------------------------------------

void Average_obdm::write_init(string & indent, ostream & os) { 
  os << indent << "obdm" << endl;
  os << indent << "ngrid " << npoints << endl;
  os << indent << "cutoff " << dR*npoints << endl;
  os << indent << "aip " << np_aver << endl;
}

//-----------------------------------------------------------------------------

void Average_obdm::read(vector <string> & words) { 
  unsigned int pos=0;
  readvalue(words, pos, npoints, "ngrid");
  doublevar cutoff;
  readvalue(words, pos=0, cutoff, "cutoff");
  dR=cutoff/npoints;
  readvalue(words, pos=0, np_aver, "aip");
}

//-----------------------------------------------------------------------------

void Average_obdm::write_summary(Average_return & avg, Average_return & err,
				 ostream & os) { 
  os << "One-body density matrix, spherically averaged with AIP = "
     << np_aver << endl;
  os << "    r  rho(r) rho(r)err" << endl;
  assert(avg.vals.GetDim(0) >=npoints);
  assert(err.vals.GetDim(0) >=npoints);
  
  for(int i=0; i< npoints; i++) {
    doublevar rho=avg.vals(i);
    doublevar rhoerr=err.vals(i);
    doublevar r=(i+1)*dR;
    os << "obdm_out " <<  r << "   " << rho << "   " << rhoerr << endl;
  }
  os << endl;

}

//############################################################################


void Average_tbdm::read(System * sys, Wavefunction_data * wfdata,
			vector <string> & words) {

  single_write(cout, "Tne-body density matrix will be calculated.\n");

  nelectrons=sys->nelectrons(0)+sys->nelectrons(1);

  unsigned int pos=0;
  if(!readvalue(words, pos=0, npoints, "NGRID")) 
    npoints=5;          
  
  // maximum distance for TBDM evaluation

  pos=0;
  doublevar cutoff;
  if(!readvalue(words, pos=0, cutoff, "CUTOFF")) {
    Array2 <doublevar> latVec(3,3);
    Array1 <doublevar> origin(3);
    if(!sys->getBounds(latVec,origin))
      error("In non-periodic systems, CUTOFF has to be given in the TBDM section.");
    cutoff=smallest_height(latVec)/2;
  }
  dR=cutoff/npoints;
  
  // angles and weights for Gaussian quadrature (spherical averaging)

  pos=0;
  if(!readvalue(words, pos=0, np_aver, "AIP")) np_aver=1;
  wt.Resize(np_aver);
  ptc.Resize(np_aver,3);             //!< cartesian coordinates of int. points
  Array1 <doublevar> x(np_aver), y(np_aver), z(np_aver);

  switch (np_aver) {
  case 1:
    // no spherical averaging
    wt=1;
    ptc(0,0)=1.0;
    ptc(0,1)=0.0;
    ptc(0,2)=0.0;
    break;
  default:
    gesqua(np_aver,x,y,z,wt);
    for (int i=0; i<np_aver; i++) {
      ptc(i,0)=x(i);
      ptc(i,1)=y(i);
      ptc(i,2)=z(i);
    }
  }
}  

//-----------------------------------------------------------------------------

void Average_tbdm::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
			    System * sys, Sample_point * sample,
			    Average_return & avg) {
  avg.type="tbdm";
  avg.vals.Resize(npoints);
  avg.vals=0;

  int nwf=wf->nfunc();
  Wf_return wfval_new(nwf,2);   // this structure I just don't understand
  Wf_return wfval_old(nwf,2);
  
  Wavefunction_storage * wfStore=NULL;
  Sample_storage * sampleStorage0=NULL, * sampleStorage1 = NULL;
 
  Array1 <doublevar> oldpos(3), newpos(3), transl(3);
  Array1 <doublevar> oldpos2(3); 

  Array1 <doublevar> x(3), y(3), z(3);
  Array2 <doublevar> pt(np_aver,3);
  generate_random_rotation(x,y,z);
  for (int i=0; i<np_aver; i++) {
    pt(i,0)=ptc(i,0)*x(0)+ptc(i,1)*y(0)+ptc(i,2)*z(0);
    pt(i,1)=ptc(i,0)*x(1)+ptc(i,1)*y(1)+ptc(i,2)*z(1);
    pt(i,2)=ptc(i,0)*x(2)+ptc(i,1)*y(2)+ptc(i,2)*z(2);
  }

  wf->updateVal(wfdata, sample);
	
  //We need to generate storage capable to store two electron moves
  wf->generateStorage( wfStore );
  sample->generateStorage( sampleStorage0 ); 
  sample->generateStorage( sampleStorage1 ); 

  //The two elecrons to be moved should have opposite spin
  int e0=0, e1=nelectrons-1;
  //assert( spin(e0) != spin(e1) );

  wf->saveUpdate(sample,e0,e1,wfStore);
  sample->saveUpdate( e0, sampleStorage0 );
  sample->saveUpdate( e1, sampleStorage1 );

  wf->getVal(wfdata, e0, wfval_old);
  sample->getElectronPos(e0,oldpos);
  sample->getElectronPos(e1,oldpos2);

  for ( int k=0; k<np_aver; k++) {
    
    transl(0)=dR*pt(k,0);
    transl(1)=dR*pt(k,1);
    transl(2)=dR*pt(k,2);

    for( int i=0; i<npoints; i++) {
      
      sample->translateElectron(e0,transl);
      sample->translateElectron(e1,transl);
      wf->updateVal(wfdata, sample);
      
      wf->getVal(wfdata, e1, wfval_new);
      
      // how to handle the possibility of wf_old=0? Can it happen?
      avg.vals(i)+=wt(k)*wfval_new.sign(0)*wfval_old.sign(0)
	*exp(wfval_new.amp(0,0)-wfval_old.amp(0,0));
    }

    //The samples need to be restore because we're using the 
    //translateElectron function (important for k-points)
    sample->setElectronPos(e0,oldpos);
    sample->setElectronPos(e1,oldpos2);

    sample->restoreUpdate(e1, sampleStorage1 );
    sample->restoreUpdate(e0, sampleStorage0 );
    wf->restoreUpdate(sample,e0,e1,wfStore);
  }

  if(wfStore != NULL) delete wfStore;
  if(sampleStorage0 != NULL) delete sampleStorage0;
  if(sampleStorage1 != NULL) delete sampleStorage1;

}

//-----------------------------------------------------------------------------

void Average_tbdm::write_init(string & indent, ostream & os) { 
  os << indent << "tbdm" << endl;
  os << indent << "ngrid " << npoints << endl;
  os << indent << "cutoff " << dR*npoints << endl;
  os << indent << "aip " << np_aver << endl;
}

//-----------------------------------------------------------------------------

void Average_tbdm::read(vector <string> & words) { 
  unsigned int pos=0;
  readvalue(words, pos, npoints, "ngrid");
  doublevar cutoff;
  readvalue(words, pos=0, cutoff, "cutoff");
  dR=cutoff/npoints;
  readvalue(words, pos=0, np_aver, "aip");
}

//-----------------------------------------------------------------------------

void Average_tbdm::write_summary(Average_return & avg, Average_return & err,
				 ostream & os) { 
  os << "Projected two-body density matrix, spherically averaged with AIP = "
     << np_aver << endl;
  os << "    r  rho(r) rho(r)err" << endl;
  assert(avg.vals.GetDim(0) >=npoints);
  assert(err.vals.GetDim(0) >=npoints);
  
  for(int i=0; i< npoints; i++) {
    doublevar rho=avg.vals(i);
    doublevar rhoerr=err.vals(i);
    doublevar r=(i+1)*dR;
    os << "tbdm_out " <<  r << "   " << rho << "   " << rhoerr << endl;
  }
  os << endl;

}


//############################################################################


void Average_local_moments::read(System * sys, Wavefunction_data * wfdata,
				 vector <string> & words) {

  single_write(cout, "Local moments and charges will be calculated.\n");

  nup=sys->nelectrons(0);
  nelectrons=nup+sys->nelectrons(1);
  
  sys->getAtomicLabels(atomnames);
  natoms=atomnames.size();
  rMT.Resize(natoms+1);
  rMT=-1;
  for ( int i=0; i < natoms; i++) {
    doublevar r;
    unsigned int pos;
    if(readvalue(words, pos=0, r, atomnames[i].c_str())) rMT(i)=r;
  }
  single_write(cout,"Muffin-tin radii:\n");
  for (int i=0; i < natoms; i++) {
    if ( rMT(i) > 0 ) {
      single_write(cout,atomnames[i],"  ",rMT(i),"\n");
    } else {
      rMT(i)=2;
      single_write(cout,atomnames[i],"  ",rMT(i)," (default used)\n");
    }
  }
  atomnames.push_back("interst.");

}  

//-----------------------------------------------------------------------------

void Average_local_moments::evaluate(Wavefunction_data * wfdata,
				     Wavefunction * wf,
				     System * sys, Sample_point * sample,
				     Average_return & avg) {

  // the content of avg.vals is:
  //    0          ... natoms-1   spin moments on atoms
  //    natoms                    spin moments in interstitial       
  //    natoms+1   ... 2*natoms   charges on atoms
  //    2*natoms+1                charge in interstitial

  avg.type="lm";
  avg.vals.Resize(2*natoms+2);
  avg.vals=0;

  Array1 <doublevar> dist(5);
  sample->updateEIDist();
  
  for(int e=0; e<nup; e++) {
    bool interstitial=true;
    for(int i=0; i<natoms; i++) {
      sample->getEIDist(e,i,dist);
      // zero element of dist is norm, second is norm^2
      if ( dist(0) < rMT(i) ) {
	avg.vals(i)+=1;
	avg.vals(i+natoms+1)+=1;
	interstitial=false;
      }
    }
    if ( interstitial ) {
      avg.vals(natoms)+=1;
      avg.vals(2*natoms+1)+=1;
    }
  }
  
  for(int e=nup; e<nelectrons; e++) {
    bool interstitial=true;
    for(int i=0; i<natoms; i++) {
      sample->getEIDist(e,i,dist);
      if ( dist(0) < rMT(i) ) {
	avg.vals(i)-=1;
	avg.vals(i+natoms+1)+=1;
	interstitial=false;
      }
    }
    if ( interstitial ) {
      avg.vals(natoms)-=1;
      avg.vals(2*natoms+1)+=1;
    }
  }

}

//-----------------------------------------------------------------------------

void Average_local_moments::write_init(string & indent, ostream & os) { 
  os << indent << "lm" << endl;
  for (int i=0; i<natoms+1; i++) {
    os << indent << "atom {"
       << " name " << atomnames[i]
       << " rMT " << rMT(i) << " }" << endl;
  }
}

//-----------------------------------------------------------------------------

void Average_local_moments::read(vector <string> & words) { 
  vector <string> atom_entry;
  vector < vector <string> > atom_entries;
  unsigned int pos=0;
  while(readsection(words, pos, atom_entry, "atom"))
    atom_entries.push_back(atom_entry);
  if ( atomnames.size() > 0 ) atomnames.clear();
  natoms=atom_entries.size()-1;
  rMT.Resize(natoms+1);
  for (int i=0; i<natoms+1; i++) {
    string atomname;
    readvalue(atom_entries[i], pos=0, atomname, "name");
    atomnames.push_back(atomname);
    readvalue(atom_entries[i], pos=0, rMT(i), "rMT");
  }
}

//-----------------------------------------------------------------------------

void Average_local_moments::write_summary(Average_return & avg,
					  Average_return & err,
					  ostream & os) { 
  os << "Local spin moments and charges integrated over muffin-tin spheres"
     << endl;
  os << "       atom name    rMT    moment   moment_err  charge  charge_err"
     << endl;
  ios::fmtflags saved_flags=os.flags(ios::fixed);
  streamsize saved_precision=os.precision(2);

  assert(avg.vals.GetDim(0) >=2*natoms+2);
  assert(err.vals.GetDim(0) >=2*natoms+2);
  
  doublevar totm=0.0;
  doublevar totc=0.0;
  for(int i=0; i< natoms+1; i++) {
    os << "lm_out  " << setw(8) << atomnames[i] << "  " << setw(5);
    if ( rMT(i) < 0.0 ) {
      os << "    ";
    } else {
      os << rMT(i);
    }
    os << setw(10) << avg.vals(i) << setw(10) << err.vals(i) << setw(10)
       << avg.vals(i+natoms+1) << setw(10) << err.vals(i+natoms+1) << endl;
    totm+=avg.vals(i);
    totc+=avg.vals(i+natoms+1);
  }
  os << "total moment in the cell" << setw(7) << totm << endl;
  os << "total charge in the cell" << setw(7) << totc << endl;
  os << endl;

  // restore output stream formating to default
  os.flags(saved_flags);
  os.precision(saved_precision);

}


//############################################################################
//Average_density_moments

void Average_density_moments::read(System * sys, Wavefunction_data * wfdata, vector <string> & words) { 
}

void Average_density_moments::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                              System * sys, Sample_point * sample, Average_return & avg ) { 
  avg.type="density_moments";
  int ndim=3;
  avg.vals.Resize(ndim);
  avg.vals=0.0;
  int nelectrons=sample->electronSize();
  Array1 <doublevar> pos(ndim);
  for(int e=0; e< nelectrons; e++) {
    sample->getElectronPos(e,pos);
    doublevar r=sqrt( pos(0)*pos(0) + pos(1)*pos(1) + pos(2)*pos(2) );
    doublevar r2=r*r;
    doublevar r3=r2*r;
    avg.vals(0)+=r;
    avg.vals(1)+=r2;
    avg.vals(2)+=r3;
  }  
}

void Average_density_moments::write_init(string & indent, ostream & os) { 
  os << indent << "density_moments \n";
}

void Average_density_moments::read(vector <string> & words) { 
  
}

void Average_density_moments::write_summary(Average_return & avg, Average_return & err, ostream & os) {
  int ndim=avg.vals.GetDim(0);
  assert(ndim <= err.vals.GetDim(0));
  //Could put this in Debye if we want to be nice.
  os << "Density moments (a.u.) \n";
  for(int d=0; d< ndim; d++) { 
    if(d==0) os << "|r| ";
    else if(d==1) os << "|r**2| ";
    else if(d==2) os << "|r**3| ";
    os << avg.vals(d) << " +/- " << err.vals(d) << endl;
  }
}

//############################################################################
//derivatives of multideterminant/pfaffian wf without jastrow(stored in symvals)
//needed un SH_DMC
void Average_linear_derivative::read(System * sys, Wavefunction_data * wfdata, vector <string> & words) {
  unsigned int pos=0;
  if(!readvalue(words, pos=0, tau, "TIMESTEP")) 
    tau=0.01;
  if(haskeyword(words, pos=0, "UNR"))
    unr=1;
  else
    unr=0;
}

void Average_linear_derivative::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                              System * sys, Sample_point * sample, Average_return & avg ) { 
  Parm_deriv_return derivatives;
  int nfunctions=wf->nfunc();
  Wf_return vals;
  Wf_return symvals;

  vals.Resize(nfunctions,5);
  symvals.Resize(nfunctions,2);
  avg.type="linear_der";
  if(!wfdata->supports(parameter_derivatives))
    error("Wavefunction needs to supports analytic parameter derivatives");
  
  derivatives.need_hessian=0;

  int ndim=wfdata->nparms();
  avg.vals.Resize(ndim+1);
  avg.vals=0.0;

  if(!unr){
    wf->updateVal(wfdata, sample);
    wf->getVal(wfdata, 0, vals);
  }
  else{
    wf->updateLap(wfdata, sample);
    wf->getLap(wfdata, 0, vals);
  }
    
  
  wf->getSymmetricVal(wfdata, 0, symvals);
  wf->getParmDeriv(wfdata, sample, derivatives);
  //cout <<symvals.amp(0,0)<<"  "<<vals.amp(0,0)<< endl;

  doublevar Jastrow_w2_inverse=exp(-2*symvals.amp(0,0));
  doublevar gamma=1.0;
  if(unr){
    //possibly drift in like UNR sampling
    doublevar drift=vals.amp(0,1)*vals.amp(0,1)+vals.amp(0,2)*vals.amp(0,2)+vals.amp(0,3)*vals.amp(0,3);
    doublevar teff=tau;
    teff*=drift;
  
    gamma=(-1+sqrt(1+2*teff))/teff;
    //if(gamma<0.5)
    // cout <<" gamma "<<gamma<<endl;
  }
  for(int d=0; d< ndim; d++){
    //cout <<" derivatives.gradient("<<d<<")= "<<derivatives.gradient(d)<<endl;
    avg.vals(d)=Jastrow_w2_inverse*derivatives.gradient(d)*gamma;
    //avg.vals(d)=derivatives.gradient(d);
  }
  avg.vals(ndim)=Jastrow_w2_inverse;
}

void Average_linear_derivative::write_init(string & indent, ostream & os) { 
  os << indent << "linear_der\n";
  os << indent << "TIMESTEP "<<tau<<"\n";
  if(unr)
    os << indent << "UNR\n";
}

void Average_linear_derivative::read(vector <string> & words) { 
  
}

void Average_linear_derivative::write_summary(Average_return & avg, Average_return & err, ostream & os) {
  int ndim=avg.vals.GetDim(0)-1;
  assert(ndim <= err.vals.GetDim(0));
  os << "Linear derivatives \n";
  for(int d=0; d< ndim; d++) {
    os << avg.vals(d)/avg.vals(ndim) << " +/- " << err.vals(d)/avg.vals(ndim) << endl;
  }
  os <<" normalization "<<avg.vals(ndim)<<" +/- " << err.vals(ndim) << endl;
      
}
//############################################################################

//############################################################################
//derivatives of multideterminant/pfaffian wf without jastrow(stored in symvals)
//needed un SH_DMC only differences
void Average_linear_delta_derivative::read(System * sys, Wavefunction_data * wfdata, vector <string> & words) { 
  unsigned int pos=0;
  if(!readvalue(words, pos=0, tau, "TIMESTEP")) 
    tau=0.01;
  if(haskeyword(words, pos=0, "UNR"))
    unr=1;
  else
    unr=0;
}

void Average_linear_delta_derivative::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                              System * sys, Sample_point * sample, Average_return & avg ) { 
  Parm_deriv_return derivatives;
  int nfunctions=wf->nfunc();
  Wf_return vals;
  Wf_return symvals;

  vals.Resize(nfunctions,5);
  symvals.Resize(nfunctions,2);
  
  avg.type="linear_delta_der";
  if(!wfdata->supports(parameter_derivatives))
    error("Wavefunction needs to supports analytic parameter derivatives");
  
  derivatives.need_hessian=0;

  int ndim=wfdata->nparms();
  avg.vals.Resize(ndim+1);
  avg.vals=0.0;

  if(!unr){
    wf->updateVal(wfdata, sample);
    wf->getVal(wfdata, 0, vals);
  }
  else{
    wf->updateLap(wfdata, sample);
    wf->getLap(wfdata, 0, vals);
  }
    
  wf->getSymmetricVal(wfdata, 0, symvals);
  wf->getParmDeriv(wfdata, sample, derivatives);
  //cout <<symvals.amp(0,0)<<"  "<<vals.amp(0,0)<< endl;

  doublevar Jastrow_w2_inverse=exp(-2*symvals.amp(0,0));
  doublevar gamma=1.0;
  if(unr){
    //possibly drift in like UNR sampling
    doublevar drift=vals.amp(0,1)*vals.amp(0,1)+vals.amp(0,2)*vals.amp(0,2)+vals.amp(0,3)*vals.amp(0,3);
    doublevar teff=tau;
    teff*=drift;
  
    gamma=(-1+sqrt(1+2*teff))/teff;
    //if(gamma<0.5)
    // cout <<" gamma "<<gamma<<endl;
  }
  
  
  for(int d=0; d< ndim; d++){
    //cout <<" derivatives.gradient("<<d<<")= "<<derivatives.gradient(d)<<endl;
    avg.vals(d)=Jastrow_w2_inverse*derivatives.gradient(d)*gamma;
    //avg.vals(d)=derivatives.gradient(d);
  }
  avg.vals(ndim)=Jastrow_w2_inverse;
}

void Average_linear_delta_derivative::write_init(string & indent, ostream & os) { 
  os << indent << "linear_delta_der\n";
  os << indent << "TIMESTEP "<<tau<<"\n";
  if(unr)
    os << indent << "UNR\n";
}

void Average_linear_delta_derivative::read(vector <string> & words) { 
  
}

void Average_linear_delta_derivative::write_summary(Average_return & avg, Average_return & err, ostream & os) {
  int ndim=avg.vals.GetDim(0)-1;
  assert(ndim <= err.vals.GetDim(0));
  os << "Linear derivatives \n";
  for(int d=0; d< ndim; d++) {
    os << avg.vals(d)/avg.vals(ndim) << " +/- " << err.vals(d)/avg.vals(ndim) << endl;
  }
  os <<" normalization "<<avg.vals(ndim)<<" +/- " << err.vals(ndim) << endl;
      
}

//############################################################################
void Average_spherical_density::read(System * sys, Wavefunction_data * wfdata,
			vector <string> & words) {

  single_write(cout, "Spherically averaged density  will be calculated.\n");

  unsigned int pos=0;
  if(!readvalue(words, pos=0, npoints, "NGRID")) 
    npoints=5;     

  pos=0;
  //doublevar cutoff;
  if(!readvalue(words, pos=0, cutoff, "CUTOFF")) {
    Array2 <doublevar> latVec(3,3);
    Array1 <doublevar> origin(3);
    if(!sys->getBounds(latVec,origin))
      error("In non-periodic systems, CUTOFF has to be given in the OBDM section.");
    cutoff=smallest_height(latVec)/2;
  }
  dR=cutoff/npoints;

  nup=sys->nelectrons(0);
  ndown=sys->nelectrons(1);

  if(!readvalue(words, pos=0, nfunc, "NFUNC")) 
    error ("need nfunc in the the Average_spherical_density");
    
  //vector <string> basistext;
  
  if(!readsection(words, pos=0, basistext, "BASIS"))
     error ("need BASIS in the the Average_spherical_density");

  allocate(basistext, basis);
  string indent="  ";
  basis->showinfo(indent, cout);
  single_write(cout, "Spherically averaged density: done read \n");

}  

//-----------------------------------------------------------------------------

void Average_spherical_density::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
			    System * sys, Sample_point * sample,
			    Average_return & avg) {

  //single_write(cout, "Spherically averaged density: start evaluate \n");
  avg.type="spherical_density";
  //avg.vals.Resize(npoints);
  avg.vals.Resize(2*nfunc);
  avg.vals=0;

  int nwf=wf->nfunc();
  Wf_return wfval_new(nwf,2);   // this structure I just don't understand
  Wf_return wfval_old(nwf,2);

  //  Storage_container wfStore;
  Array1 <doublevar> oldpos(3), newpos(3), transl(3);
  Array1 <doublevar> x(3), y(3), z(3);
  for(int e=0; e < nup+ndown; e++){
    sample->getElectronPos(e,oldpos);
    doublevar distance2=oldpos(0)*oldpos(0)+oldpos(1)*oldpos(1)+oldpos(2)*oldpos(2);
    doublevar distance=sqrt(distance2);
    //for( int i=0; i<npoints; i++) {
    //  if(distance > i*dR && distance< (i+1)*dR)
    //	avg.vals(i)+= 1/distance2;
    //}
    Array1 <doublevar> R(5);
    R(0)=distance;
    Array1 <doublevar> basisvals(nfunc);
    int currfunc=0;
    basis->calcVal(R, basisvals, currfunc);
    for( int i=0; i<nfunc; i++) {
      if(e<nup)
	avg.vals(i)+=basisvals(i)/(4.0*pi);
      //avg.vals(i)+=distance2*basisvals(i);
      else
	avg.vals(i+nfunc)+=basisvals(i)/(4.0*pi);
      //avg.vals(i+nfunc)+=distance2*basisvals(i);
    }
  }
}

//-----------------------------------------------------------------------------

void Average_spherical_density::write_init(string & indent, ostream & os) { 
  os << indent << "spherical_density" << endl;
  os << indent << "ngrid " << npoints << endl;
  os << indent << "cutoff " << dR*npoints << endl;
  os << indent << "nfunc " <<  nfunc << endl;
  os << indent << "BASIS { ";
  for(unsigned int i=0; i<basistext.size();i++)
    os << basistext[i]<<" ";
  os <<" } "<<endl;
}

//-----------------------------------------------------------------------------

void Average_spherical_density::read(vector <string> & words) { 
  unsigned int pos=0;
  readvalue(words, pos, npoints, "ngrid");
  //doublevar cutoff;
  readvalue(words, pos=0, cutoff, "cutoff");
  dR=cutoff/npoints;
  readvalue(words, pos=0, nfunc, "nfunc");
  //vector <string> basistext;
  readsection(words, pos=0, basistext, "BASIS");
  allocate(basistext, basis);
}

//-----------------------------------------------------------------------------

void Average_spherical_density::write_summary(Average_return & avg, Average_return & err,
				 ostream & os) { 
  os << "Spherically averaged density "<< endl;
  
  /*
  assert(avg.vals.GetDim(0) >=npoints);
  assert(err.vals.GetDim(0) >=npoints);

  doublevar integral=0.0;
  for(int i=0; i< npoints; i++) {		
    doublevar r=(i+0.5)*dR;
    doublevar rho=avg.vals(i);
    integral+=rho*r*r*dR;
  }
  cout <<" integral "<<integral<<" nelectrons "<<nelectrons<<endl;
  
  for(int i=0; i< npoints; i++) {
    avg.vals(i)*=1/integral;
    err.vals(i)*=1/integral;
    doublevar rho=avg.vals(i);
    doublevar rhoerr=err.vals(i);
    doublevar r=(i+0.5)*dR;
    os << "density_out " <<  r << "   " << rho << "   " << rhoerr << endl;
  }
  
  
  doublevar summ1=0;
  doublevar summ2=0;
  for(int i=0; i< npoints; i++) {
    doublevar rho=0.0;
    avg.vals(i);
    
    doublevar r=(i+0.5)*dR;
    os << "density_out " <<  r << "   " << rho << "   " << rhoerr << endl;
  }
  */
  
  
  os.setf(ios::scientific);
  os <<"Used basis info :"<<endl;
  string indent="  ";
  basis->showinfo(indent, os);
  os <<"Spin-up and spin-down density in above basis"<<endl;
  //doublevar summ1=0;
  //doublevar summ2=0;
  for(int i=0; i< nfunc; i++) {
    os << i+1 <<setw(20)<<setprecision(10)<< avg.vals(i) <<setw(20)<<setprecision(10)<<  err.vals(i);
    os <<setw(20)<<setprecision(10)<< avg.vals(i+nfunc) <<setw(20)<<setprecision(10)<<  err.vals(i+nfunc) << endl;
  }
  
  
  Array1 <doublevar> R(5);
  Array1 <doublevar> basisvals(nfunc); 
  os << "densities on the grid"<<endl;
  os << "r  rho(r) rho(r)err  rho(r) rho(r)err"  << endl;
  doublevar  norm1=0.0;
  doublevar  norm2=0.0;
  //for(int i=0;i<npoints;i++){
  //R(0)=(i+0.5)*dR;
  R(0)=1.000000000000e-06;
  doublevar lastR=0.0;
  doublevar rf=100.0;
  doublevar ri=1e-6;
  int npts=5001;
  doublevar ratio=exp(log(rf/ri)/(npts-1)); //1.003690930920;
  
  while (R(0)<cutoff) {
    int currfunc=0;
    basis->calcVal(R, basisvals, currfunc);
    Array1 <doublevar> sum(6);
    sum=0;				
    Array1 <doublevar> sum_error(6);
    sum_error=0;
    for(int k=0;k<nfunc;k++){
      sum(0)+=basisvals(k)*avg.vals(k);
      sum_error(0)+=basisvals(k)*err.vals(k)*basisvals(k)*err.vals(k);
      /*
      if(abs(avg.vals(k))>2.0*err.vals(k)){
	sum(1)+=basisvals(k)*avg.vals(k);
	sum_error(1)+=basisvals(k)*err.vals(k)*basisvals(k)*err.vals(k);
      }
      if(abs(avg.vals(k))>4.0*err.vals(k)){
	sum(2)+=basisvals(k)*avg.vals(k);
	sum_error(2)+=basisvals(k)*err.vals(k)*basisvals(k)*err.vals(k);
      }
      */

      sum(3)+=basisvals(k)*avg.vals(k+nfunc);
      sum_error(3)+=basisvals(k)*err.vals(k+nfunc)*basisvals(k)*err.vals(k+nfunc);
	
      /*
      if(abs(avg.vals(k+nfunc))>2.0*err.vals(k+nfunc)){
	sum(4)+=basisvals(k)*avg.vals(k+nfunc);
	sum_error(4)+=basisvals(k)*err.vals(k+nfunc)*basisvals(k)*err.vals(k+nfunc);
      }
      if(abs(avg.vals(k+nfunc))>4.0*err.vals(k+nfunc)){
	sum(5)+=basisvals(k)*avg.vals(k+nfunc);
	sum_error(5)+=basisvals(k)*err.vals(k+nfunc)*basisvals(k)*err.vals(k+nfunc);
      }
      */
      
    }//k
    
    norm1+=4.0*pi*R(0)*R(0)*sum(0)*(R(0)-lastR);
    //norm1+=sum(0)*dR;
    norm2+=4.0*pi*R(0)*R(0)*sum(3)*(R(0)-lastR);
    //norm2+=sum(3)*dR;
    os <<R(0)<<setw(20)<<setprecision(10)<<sum(0)<<setw(20)<<setprecision(10)<<sqrt(sum_error(0))
      //<<setw(20)<<setprecision(10)<<sum(1)<<setw(20)<<setprecision(10)<<sqrt(sum_error(1))
      //<<setw(20)<<setprecision(10)<<sum(2)<<setw(20)<<setprecision(10)<<sqrt(sum_error(2))
      <<setw(20)<<setprecision(10)<<sum(3)<<setw(20)<<setprecision(10)<<sqrt(sum_error(3))
      //<<setw(20)<<setprecision(10)<<sum(4)<<setw(20)<<setprecision(10)<<sqrt(sum_error(4))
      //<<setw(20)<<setprecision(10)<<sum(5)<<setw(20)<<setprecision(10)<<sqrt(sum_error(5))
       <<endl;
    lastR=R(0);
    R(0)*=ratio;
    
  }//i
  os <<"integrated  spin-up density "<<norm1<<" and spin-up density "<<norm2<<endl;
  os << endl;
  os.unsetf(ios::scientific);
  os<<setprecision(6);
}







//-----------------------------------------------------------------------------

void Average_spherical_density_grid::read(System * sys, Wavefunction_data * wfdata,
			vector <string> & words) {

  single_write(cout, "Spherically averaged density  on the grid will be calculated.\n");

  unsigned int pos=0;
  if(!readvalue(words, pos=0, npoints, "NGRID")) 
    npoints=5;     

  pos=0;
  doublevar cutoff;
  if(!readvalue(words, pos=0, cutoff, "CUTOFF")) {
    Array2 <doublevar> latVec(3,3);
    Array1 <doublevar> origin(3);
    if(!sys->getBounds(latVec,origin))
      error("In non-periodic systems, CUTOFF has to be given in the OBDM section.");
    cutoff=smallest_height(latVec)/2;
  }
  dR=cutoff/npoints;

  nup=sys->nelectrons(0);
  ndown=sys->nelectrons(1);
}  

//-----------------------------------------------------------------------------

void Average_spherical_density_grid::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
			    System * sys, Sample_point * sample,
			    Average_return & avg) {

  //single_write(cout, "Spherically averaged density on grid : start evaluate \n");
  avg.type="spherical_density_grid";
  avg.vals.Resize(2*npoints);
  avg.vals=0;

  //  Storage_container wfStore;
  Array1 <doublevar> oldpos(3), newpos(3), transl(3);
  Array1 <doublevar> x(3), y(3), z(3);
  for(int e=0; e < nup+ndown; e++){
    sample->getElectronPos(e,oldpos);
    doublevar distance2=oldpos(0)*oldpos(0)+oldpos(1)*oldpos(1)+oldpos(2)*oldpos(2);
    doublevar distance=sqrt(distance2);
    if(distance2 < 1.0e-6){
      cout <<"electron # "<<e+1<< "very close to origin: R= "<<distance<<endl;
    }
    for( int i=0; i<npoints; i++) {
      if(distance > i*dR && distance< (i+1)*dR){
	//doublevar r=(i+0.5)*dR;
	if(e<nup){
	  // if(i==0){
	  //  cout <<"electron # "<<e+1<< "very close to origin: R= "<<distance<<" value "<<1.0/(distance2*4.0*pi*dR)<<endl;
	  //  if(distance > 0.5*dR)
	  //    avg.vals(i)+= 1.0/(distance2*4.0*pi*dR*0.5);
	  // }
	  //else{
	  //avg.vals(i)+= 1.0/(distance2*4.0*pi*dR);
	  avg.vals(i)+= 1.0/(4.0*pi*dR);
	  
	  //}
	}
	else{
	  //if(i==0){
	  // cout <<"electron # "<<e+1<< "very close to origin: R= "<<distance<<" value "<<1.0/(distance2*4.0*pi*dR)<<endl;
	  //  if(distance > 0.5*dR)
	  //    avg.vals(i+npoints)+= 1.0/(distance2*4.0*pi*dR*0.5);
	  //}
	  //else{
	  //avg.vals(i+npoints)+= 1.0/(distance2*4.0*pi*dR);
	  avg.vals(i+npoints)+= 1.0/(4.0*pi*dR);
	    //}
	}
      }
    }
  }
  
  // for( int i=0; i<npoints; i++)
  //  cout <<avg.vals(i)<<" "<<avg.vals(i+npoints)<<endl;


}

//-----------------------------------------------------------------------------

void Average_spherical_density_grid::write_init(string & indent, ostream & os) { 
  os << indent << "spherical_density_grid" << endl;
  os << indent << "ngrid " << npoints << endl;
  os << indent << "cutoff " << dR*npoints << endl;
  os << indent << "nup " << nup << endl;
  os << indent << "ndown " << ndown << endl;
}

//-----------------------------------------------------------------------------

void Average_spherical_density_grid::read(vector <string> & words) { 
  unsigned int pos=0;
  readvalue(words, pos, npoints, "ngrid");
  doublevar cutoff;
  readvalue(words, pos=0, cutoff, "cutoff");
  dR=cutoff/npoints;
  readvalue(words, pos=0, nup, "nup");
  readvalue(words, pos=0, ndown, "ndown");
}

//-----------------------------------------------------------------------------

void Average_spherical_density_grid::write_summary(Average_return & avg, Average_return & err,
				 ostream & os) { 
  os << "Spherically averaged density on grid"<< endl;
  
  
  assert(avg.vals.GetDim(0) >=2*npoints);
  assert(err.vals.GetDim(0) >=2*npoints);

  doublevar integral1=0.0;
  doublevar integral2=0.0;
  doublevar alpha=10.0/(npoints*dR*npoints*dR);
  doublevar check_intg=0.0;
  doublevar C=sqrt(pi)/(4.0*pow(alpha,1.5));
  
  for(int i=0; i< npoints; i++) {		
    doublevar r=(i+0.5)*dR;
    doublevar rho1=avg.vals(i);
    doublevar rho2=avg.vals(i+npoints);
    //integral1+=4*pi*rho1*r*r*dR;
    integral1+=4*pi*rho1*dR;
    //integral2+=4*pi*rho2*r*r*dR;
    integral2+=4*pi*rho2*dR;
    check_intg+=dR*r*r*exp(-alpha*r*r)/C;
    
  }
  //os.setf(ios::scientific);

  //os <<" alpha "<<alpha<<" C "<<C<<" integration error "<<setw(20)<<abs(1.0-check_intg)<<endl;
  os <<"spin-up integtated density "<<integral1<<" # electrons "<<nup<<"  spin-down integtated density "<<integral2<<" # electrons "<<ndown
     <<" integration error "<<setw(20)<<abs(1.0-check_intg)<<endl;
  
  os.setf(ios::scientific);
  for(int i=0; i< npoints; i++) {
    //avg.vals(i)*=nup/integral1;
    //err.vals(i)*=nup/integral1;
    
    //avg.vals(i+npoints)*=ndown/integral2;
    //err.vals(i+npoints)*=ndown/integral2;
    doublevar r=(i+0.5)*dR; //plot in the middle of the interval
    //doublevar r2=r*r;
    doublevar rho1=avg.vals(i);
    doublevar rhoerr1=err.vals(i);
    doublevar rho2=avg.vals(i+npoints);
    doublevar rhoerr2=err.vals(i+npoints);

    
    os << "density_out " <<setw(20)<<  r<< setw(20)<<"  "<<setw(20)<< rho1 <<"  "<<setw(20)<<rhoerr1<<setw(20)<< rho2 <<"  "<<setw(20)<<rhoerr2<< endl;
  }
  os.unsetf(ios::scientific);
  os<<setprecision(6);
}

//############################################################################
void Average_line_density::read(System * sys, Wavefunction_data * wfdata, 
                                       vector <string> & words) { 
  unsigned int pos=0;
  doublevar length=10.0;
  vec.Resize(3);
  origin.Resize(3);
  vec=0; vec(2)=1.0;
  origin=0;
  if(!readvalue(words,pos=0,resolution, "RESOLUTION"))
    resolution=0.1;
  readvalue(words,pos=0,length,"LENGTH");
  vector<string> dirwords;
  if(readsection(words, pos=0, dirwords, "DIRECTION")) { 
    if(dirwords.size()==3) { 
      for(int i=0; i< 3; i++) { 
        vec(i)=atof(dirwords[i].c_str());
      }
    }
    else { error("Wrong number of elements in DIRECTION"); } 
  }
  vector <string> originwords;
  if(readsection(words,pos=0,originwords, "ORIGIN" )) { 
    if(originwords.size()!=3) error("Wrong number of elements in ORIGIN");
    for(int i=0; i<3; i++) { 
      origin(i)=atof(originwords[i].c_str());
    }
  }

  double norm=0;
  for(int d=0; d< 3; d++) { 
    norm+=vec(d);
  }
  for(int d=0; d< 3; d++) vec(d)/=sqrt(norm);
  
  npoints=int(length/resolution)+1;  
  
}

//-----------------------------------------------------------------------------

void Average_line_density::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                                        System * sys, Sample_point * sample, Average_return & avg) {
  avg.type="line_density";
  int nelectrons=sample->electronSize();
  avg.vals.Resize(2*npoints); //for up and down
  avg.vals=0;
  int nup=sys->nelectrons(0);
  Array1 <doublevar> epos(3);
  for(int e=0; e< nelectrons; e++) {
    sample->getElectronPos(e,epos);
    for(int d=0; d< 3; d++) 
      epos(d)-=origin(d);
    double scalar=0;
    for(int d=0; d< 3; d++) 
      scalar+=epos(d)*vec(d);

    int place=int(fabs(scalar)/resolution);
    if(place < npoints && place > 0) { 
      if(e< nup)
        avg.vals(place)+=1;
      else
        avg.vals(npoints+place)+=1;
    }
  } 
  
}

//-----------------------------------------------------------------------------

void Average_line_density::write_init(string & indent, ostream & os) { 
  os << indent << "line_density\n";
  os << indent << "npoints " << npoints << endl;
  os << indent << "resolution " << resolution << endl;
  os << indent << "direction { ";
  for(int d=0; d< 3; d++) os << vec(d) << " ";
  os << " } " << endl;
  os << indent << "origin { ";
  for(int d=0; d< 3; d++) os << vec(d) << " ";
  os << " } " <<  endl;
}
//-----------------------------------------------------------------------------

void Average_line_density::read(vector <string> & words) { 
  unsigned int pos=0;
  readvalue(words, pos=0, resolution, "resolution");
  readvalue(words, pos=0, npoints, "npoints");
  vector <string> sec;
  readsection(words, pos=0,sec, "direction");
  vec.Resize(sec.size());
  for(unsigned int i=0; i< sec.size(); i++) vec(i)=atof(sec[i].c_str());
  sec.clear();
  readsection(words,pos=0, sec, "origin");
  origin.Resize(sec.size());
  for(unsigned int i=0; i< sec.size(); i++) origin(i)=atof(sec[i].c_str());

}
//-----------------------------------------------------------------------------

void Average_line_density::write_summary(Average_return & avg, Average_return & err, ostream & os) { 
  os << "Electron Density along a line. \n";
  os << "    r  p(r) sigma(p(r))   p(r)  sigma(p(r))" << endl;
  assert(avg.vals.GetDim(0) >=2*npoints);
  assert(err.vals.GetDim(0) >=2*npoints);
  
  for(int i=0; i< npoints; i++) {
    os << "line_out " << i*resolution << " " << avg.vals(i) << " " << err.vals(i) 
    << "  " << avg.vals(i+npoints) << "  " << err.vals(i+npoints) << endl;
  }
  
}

//######################################################################

void Average_wf_parmderivs::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                       System * sys, Sample_point * sample, Average_return & avg) { 
  error("Must use the properties_point version of evaluate with wf_parmderivs");
}


void Average_wf_parmderivs::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
                       System * sys, Sample_point * sample, Properties_point & pt,
                       Average_return & avg) { 
  error("Must use pseudopotential version of evaluate() with wf_parmderivs");
}
//-----------------------------------------------------------------------------


void Average_wf_parmderivs::evaluate(Wavefunction_data * wfdata, Wavefunction * wf,
			System * sys, Pseudopotential * psp, Sample_point * sample,  Properties_point & pt,  Average_return &avg ){ 
  
  Parm_deriv_return deriv;
  deriv.need_hessian=0;
  int nparms=wfdata->nparms();
  avg.type="wf_parmderivs";
  if(!wf->getParmDeriv(wfdata, sample,deriv)) { 
    error("WF needs to support parmderivs for now.");
  }
  //Calculate distance squared from node, per electron
  wf->updateLap(wfdata,sample);
  doublevar d2=0.0;
  assert(wf->nfunc()==1);
  Wf_return lap(1,5);  
  for(int e=0;e<sample->electronSize();e++){
    wf->getLap(wfdata,e,lap);
    for(int i=1;i<4;i++) d2+=lap.amp(0,i)*lap.amp(0,i);
  }
  d2=(1./d2);
  d2/=(sample->electronSize()*sample->electronSize());
  if(d2 < nodal_cutoff*nodal_cutoff) { 
    cout << d2 << " " << nodal_cutoff << endl;
    deriv.gradient=0.0;
    deriv.gradderiv=0.0;
  }

  Array1<doublevar> psp_der(nparms),test(psp->nTest()),totalv(wf->nfunc());
  for(int i=0; i< psp->nTest(); i++) test(i)=rng.ulec();
  psp_der=0;
  if(evaluate_pseudopotential)
    psp->calcNonlocParmDeriv(wfdata, sys,sample,wf,test,totalv,psp_der);
  
  avg.vals.Resize(3*nparms+3*nparms*nparms+1);
  avg.vals=0.0;
  
  for(int i=0; i< nparms; i++) { 
    avg.vals(i)=deriv.gradient(i)*pt.energy(0);
    avg.vals(nparms+i)=deriv.gradient(i);
    //cout << "grad " << deriv.gradient(i) << " en " << pt.energy(0) << endl;
  }
  for(int i=0;i< nparms; i++) { 
    for(int j=0; j< nparms; j++) { 
      avg.vals(3*nparms+i*nparms+j)=deriv.gradient(i)*deriv.gradient(j);
    }
  }
  int offset=3*nparms+nparms*nparms;
  for(int i=0; i< nparms; i++) { 
    for(int j=0; j< nparms; j++) { 
      avg.vals(offset+i*nparms+j)=deriv.gradient(i)*deriv.gradient(j)*pt.energy(0);
    }
  }
  Array1 <doublevar> el(nparms),kin(1);
  int nelectrons=sample->electronSize();
  for(int i=0; i< nparms; i++) { 
     el(i)=0;
     for(int e=0; e< nelectrons; e++) { 
       el(i)+=-0.5*deriv.gradderiv(i,e,3);
     }
     el(i)+=psp_der(i);
  }

  for(int i=0; i< nparms; i++) { 
    avg.vals(2*nparms+i)=el(i);
  }
  offset+=nparms*nparms;
  for(int i=0; i< nparms; i++) {
    for(int j=0; j< nparms; j++) { 
      avg.vals(offset+i*nparms+j)=deriv.gradient(i)*el(j);
    }
  }
  
  avg.vals(offset+nparms*nparms)=exp(2*pt.wf_val.amp(0,0));
}
//-----------------------------------------------------------------------------
void Average_wf_parmderivs::read(System * sys, Wavefunction_data * wfdata, vector
                   <string> & words) { 
  unsigned int pos=0;
  evaluate_pseudopotential=haskeyword(words,pos=0,"EVALUATE_PSEUDOPOTENTIAL");
  if(!readvalue(words,pos=0,nodal_cutoff,"NODAL_CUTOFF"))
    nodal_cutoff=1e-6;

}
//-----------------------------------------------------------------------------
void Average_wf_parmderivs::write_init(string & indent, ostream & os) { 
  os << indent << "WF_PARMDERIV" << endl;
  if(evaluate_pseudopotential) os <<indent << "EVALUATE_PSEUDOPOTENTIAL" << endl;
  os << indent << "NODAL_CUTOFF " << nodal_cutoff << endl;
}
//-----------------------------------------------------------------------------
void Average_wf_parmderivs::read(vector <string> & words) { 
  unsigned int pos=0;
  evaluate_pseudopotential=haskeyword(words,pos=0,"EVALUATE_PSEUDOPOTENTIAL");
  nodal_cutoff=0;
  readvalue(words,pos=0,nodal_cutoff,"NODAL_CUTOFF");
  
}
//-----------------------------------------------------------------------------
void Average_wf_parmderivs::write_summary(Average_return &avg ,Average_return & err, ostream & os) { 
  os << "Wavefunction parameter derivatives. ";
  if(evaluate_pseudopotential) os << "Pseudopotential derivatives evaluated.";
  os << endl;

  int n=sqrt(1.0+avg.vals.GetDim(0))-1;

  os << "energy " << endl;
  for(int i=0; i < n; i++) { 
    os << i << " " << avg.vals(i) << " +/- " << err.vals(i) 
      << " " << avg.vals(n+i) << " +/- " << err.vals(n+i) <<  endl;
  }
   
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void Average_wf_parmderivs::jsonOutput(Average_return &avg ,Average_return & err, ostream & os) { 
  os << "\"" << avg.type << "\":{" << endl;
   int M=avg.vals.GetDim(0);
   int n=(sqrt(1.0-4*(1-M)/3.)-1)/2;
   os << "\"n\":"<< n<<",\n";
   
   Array1<doublevar> dpE(n),dpE_err(n),dp(n),dp_err(n),dE(n),dE_err(n);


   for(int i=0; i < n; i++) { 
     dpE(i)=avg.vals(i);
     dpE_err(i)=err.vals(i);
     dp(i)=avg.vals(n+i);
     dp_err(i)=err.vals(n+i);
     dE(i)=avg.vals(2*n+i);
     dE_err(i)=err.vals(2*n+i);
   }


   Array2<doublevar> dpij(n,n),dpij_err(n,n);
   Array2<doublevar> dpijE(n,n),dpijE_err(n,n),dpidE(n,n),dpidE_err(n,n);
   for(int i=0; i< n; i++) { 
     for(int j=0; j< n; j++) {
       dpij(i,j)=avg.vals(3*n+i*n+j);
       dpij_err(i,j)=err.vals(3*n+i*n+j);
     }
   }
   int offset=3*n+n*n;
   for(int i=0; i< n; i++) { 
     for(int j=0; j< n; j++) {
       dpijE(i,j)=avg.vals(offset+i*n+j);
       dpijE_err(i,j)=err.vals(offset+i*n+j);
     }
   }
   offset+=n*n;
   for(int i=0; i< n; i++) { 
     for(int j=0; j< n; j++) {
       dpidE(i,j)=avg.vals(offset+i*n+j);
       dpidE_err(i,j)=err.vals(offset+i*n+j);
     }
   }
   
   
   os << "\"dpE\":";
   jsonarray(os,dpE);
   os << ",\n";

   os << "\"dpE_err\":";
   jsonarray(os,dpE_err);
   os << ",\n";

   os << "\"dp\":";
   jsonarray(os,dp);
   os << ",\n";

   os << "\"dp_err\":";
   jsonarray(os,dp_err);
   os << ",\n";


   os << "\"dE\":";
   jsonarray(os,dE);
   os << ",\n";

   os << "\"dE_err\":";
   jsonarray(os,dE_err);
   os << ",\n";
   

   os << "\"dpij\":";
   jsonarray(os,dpij);
   os << ",\n";

   os << "\"dpij_err\":";
   jsonarray(os,dpij_err);
   os << ",\n";

   os << "\"dpijE\":";
   jsonarray(os,dpijE);
   os << ",\n";

   os << "\"dpijE_err\":";
   jsonarray(os,dpijE_err);
   os << ",\n";
   
   os << "\"dpidE\":";
   jsonarray(os,dpidE);
   os << ",\n";

   os << "\"dpidE_err\":";
   jsonarray(os,dpidE_err);
   os << "\n";
   
   
   os << "}\n";
   
}
//-----------------------------------------------------------------------------



