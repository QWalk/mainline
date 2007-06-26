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
//------------------------------------------------------------------------
//src/Onebody_int_method.cpp
#include "Onebody_int_method.h"
#include "qmc_io.h"
#include "System.h"
#include "Program_options.h"
/*!
Read the "words" from the method section in the input file
via doinput() parsing, and store section information in private
variables orbs, resolution, and minmax.
Set up MO_matrix and Sample_point objects for wavefunction from input
*/
void Onebody_int_method::read(vector <string> words,
                       unsigned int & pos,
                       Program_options & options)
{
  pos=0;	//always start from first word
  vector <string> Torbs;
  vector <string> Tminmax;
  vector <string> orbtext;


  if(!readsection(words, pos=0, orbtext, "ORBITALS"))
    error("Need ORBITALS in ONEBODY_INT section");

  if(! readvalue(words,pos=0,resolution,"RESOLUTION"))
    resolution=.1;

  sys=NULL;
  allocate(options.systemtext[0],  sys);
  allocate(orbtext, sys, mymomat);


  vector <string> statestxt;
  if(!readsection(words, pos=0, statestxt,"STATES"))
    error("Need STATES in ONEBODY_INT section");

  orbs.Resize(2);
  int nelec=sys->nelectrons(0)+sys->nelectrons(1);
  if(statestxt.size()!=nelec)
    error("states size doesn't match nelectrons");
  
  int counter=0;
  for(int s=0; s< 2; s++) { 
    orbs(s).Resize(sys->nelectrons(s));
    for(int e=0; e< sys->nelectrons(s); e++) {
      orbs(s)(e)=atoi(statestxt[counter].c_str())-1;
      counter++;
    }
  }

  mymomat->buildLists(orbs);

  mywalker=NULL;
  sys->generateSample(mywalker);
}

#include "MatrixAlgebra.h"

void divide_job_1b(int node, int nprocs, int npoints, int & start, int & end) {
  start=node*npoints/nprocs;
  end=(node+1)*npoints/nprocs;
}


/*!
 
*/
void Onebody_int_method::run(Program_options & options, ostream & output) {
  ofstream os; //for writing to *.plt files
  unsigned int electron=0; //# of the electron that will move through grid
  string pltfile; //name of plotfile being written
  Array1 <doublevar> xyz(3),resolution_array(3); //position of electron "in" MO
  Array1 <int> D_array1(3); //dummy array1
  D_array1=0; //sets all 3 components to 0. use as counter for gridpoints
  Array1 <doublevar> minmax(6);
  
  Array2 <doublevar> latvec(3,3);
  if(!sys->getBounds(latvec)) error("sys must support getBounds");
  minmax(0)=minmax(2)=minmax(3)=0;
  for(int i=0; i< 3; i++) { 
    doublevar len=0;
    for(int j=0; j< 3; j++) len+=latvec(i,j)*latvec(i,j);
    len=sqrt(len);
    minmax(2*i+1)=len;
  }


  D_array1(0)=roundoff((minmax(1)-minmax(0))/resolution);
  D_array1(1)=roundoff((minmax(3)-minmax(2))/resolution);
  D_array1(2)=roundoff((minmax(5)-minmax(4))/resolution);

  //resolution_array(0)=(minmax(1)-minmax(0))/(D_array1(0)-1);  
  //resolution_array(1)=(minmax(3)-minmax(2))/(D_array1(1)-1);
  //resolution_array(2)=(minmax(5)-minmax(4))/(D_array1(2)-1);
  resolution_array(0)=1.0/(D_array1(0)-1);  
  resolution_array(1)=1.0/(D_array1(1)-1);
  resolution_array(2)=1.0/(D_array1(2)-1);

  int npts=D_array1(0)*D_array1(1)*D_array1(2);
  //Array2 <doublevar> grid(orbs.GetSize(),npts);
  //Array1 <doublevar> density(npts);
  //generate .xyz file for gOpenMol to view coordinates
  pltfile=options.runid + ".xyz";
  os.open(pltfile.c_str());
  write_xyz(sys,os);
  os.close();

  Array2 <doublevar> mymovals(max(orbs(0).GetDim(0), orbs(1).GetDim(0)),1);

  Array1 <Array2 <dcomplex> > overlap_matrix(2);
  Array1 <Array1 <doublevar> > norm(2);
  Array1 <Array2 <dcomplex> > overlap_com(2);

  Array2 <doublevar> density(sys->nelectrons(0),D_array1(2));
  density=0;


  for(int s=0; s< 2; s++) {
    overlap_matrix(s).Resize(sys->nelectrons(s),sys->nelectrons(s));
    overlap_matrix(s)=dcomplex(0.0,0.0);
    overlap_com(s).Resize(sys->nelectrons(s),sys->nelectrons(s));
    overlap_com(s)=dcomplex(0.0,0.0);
    norm(s).Resize(sys->nelectrons(s));
    norm(s)=0.0;
  }
  

  vector <doublevar> cvg_tmp;
  cvg_tmp.push_back(1.0);
  for(double i=1; i < 20; i++) cvg_tmp.push_back((sys->nelectrons(0)+sys->nelectrons(1))/i);
  //cvg_tmp.push_back(1);
  //cvg_tmp.push_back(sys->nelectrons(0)+sys->nelectrons(1));
  //for(int i=1; i <= sys->nelectrons(0)+sys->nelectrons(1); i++) {
    //if(fabs(sys->nelectrons(0)/i-double(sys->nelectrons(0))/double(i))< 1e-8) {
    //  cvg_tmp.push_back(i);
    //  output << "adding " << i << endl;
    //}
  //}
  Array1 <doublevar> cvg;
  cvg.Resize(cvg_tmp.size());
  for(int i=0; i< cvg.GetDim(0); i++) cvg(i)=cvg_tmp[i];
  int ncvg=cvg.GetDim(0);
  
  Array1 < Array2 <dcomplex> > overlap_com_cvg(ncvg);
  for(int i=0; i< ncvg; i++) {
    overlap_com_cvg(i).Resize(sys->nelectrons(0),sys->nelectrons(0));
    overlap_com_cvg(i)=dcomplex(0.0,0.0);
  }
  


  Array2 <doublevar> gvec;
  if(!sys->getRecipLattice(gvec)) error("need reciplattice");

  int nions=mywalker->ionSize();
  Array2 <doublevar> ioncharge(sys->nelectrons(0),nions,0.0);
  Array1 <doublevar> eidist(5);

  doublevar expansion=1.0;

  output <<"calculating "<<D_array1(0)*D_array1(1)*D_array1(2) <<" grid points"<<endl;
  output <<"for "<< orbs.GetDim(0) <<" molecular orbitals"<<endl;
  xyz(0)=minmax(0);
  xyz(1)=minmax(2);
  xyz(2)=minmax(4); //init elec probe to xmin ymin zmin

  Array2 <doublevar> xavg(sys->nelectrons(0),3,0.0);

 
  int start; int end;

  divide_job_1b(mpi_info.node, mpi_info.nprocs, D_array1(0),
	     start, end);

  //for(int xx=0;xx<D_array1(0);xx++){
  for(int xx=start; xx < end; xx++) {
    double perc=double(xx-start)/double(end-start);
    output << perc*100 << " percent done " << endl;

    for(int yy=0; yy<D_array1(1);yy++){
      for(int zz=0; zz<D_array1(2);zz++){
	xyz=0;
	for(int i=0; i< 3; i++) {
	  xyz(i)+=xx*resolution_array(0)*latvec(0,i)
	    +yy*resolution_array(1)*latvec(1,i)
	    +zz*resolution_array(2)*latvec(2,i);
	}
	mywalker->setElectronPos(electron,xyz); 
	mymomat->updateVal(mywalker,electron,0,mymovals); 
	mywalker->updateEIDist();
	//---Find the closest atom
	doublevar smalldist=1e99;
	int closestatom=-1;
	for(int at=0; at< nions; at++) {
	  mywalker->getEIDist(0,at,eidist);
	  //ut << "at " << at << " eidist0 " << eidist(0) << endl;
	  if(eidist(0)<smalldist) {
	    smalldist=eidist(0);
	    closestatom=at;
	  }
	}
	
	

	doublevar sum=0;
	for(int d=0; d< 3; d++) { 
	  //cout <<"gvec " << gvec(2,d) << "  xyz " << xyz(d) << endl;
	  sum+=expansion*gvec(2,d)*xyz(d);
	}
	
	dcomplex comp_part(cos(2*pi*sum),sin(2*pi*sum));
	
	doublevar sum_com=sum/sys->nelectrons(0);
	dcomplex comp_com(cos(2*pi*sum_com),sin(2*pi*sum_com));
	
	Array1 <dcomplex> comp_cvg(ncvg);
	for(int i=0; i< ncvg; i++) {
	  doublevar sum_tmp=sum/cvg(i);
	  comp_cvg(i)=dcomplex(cos(2*pi*sum_tmp),sin(2*pi*sum_tmp));
	}


	for(int j=0; j< sys->nelectrons(0);j++) { 
	  norm(0)(j)+=mymovals(j,0)*mymovals(j,0)/npts;
	  density(j,zz)+=mymovals(j,0)*mymovals(j,0)/npts;
	  ioncharge(j,closestatom)+=mymovals(j,0)*mymovals(j,0)/npts;
	  for(int d=0; d< 3; d++) 
	    xavg(j,d)+=mymovals(j,0)*mymovals(j,0)*xyz(d)/npts;
	  for(int jp=0; jp< sys->nelectrons(0); jp++) { 
	    doublevar fac=mymovals(j,0)*mymovals(jp,0);
	    overlap_matrix(0)(j,jp)+=fac*comp_part/doublevar(npts);
	    overlap_com(0)(j,jp)+=fac*comp_com/doublevar(npts);
	    for(int i=0; i< ncvg; i++) {
	      overlap_com_cvg(i)(j,jp)+=fac*comp_cvg(i)/doublevar(npts);
	    }
	  }
	}
	
	
      }
    }
  }

  

  for(int j=0; j< sys->nelectrons(0); j++) {
    norm(0)(j)=parallel_sum(norm(0)(j));
    for(int jp=0; jp< sys->nelectrons(0);jp++) {
      overlap_matrix(0)(j,jp)=parallel_sum(overlap_matrix(0)(j,jp));
      overlap_com(0)(j,jp)=parallel_sum(overlap_com(0)(j,jp));
      for(int i=0; i < ncvg; i++) 
	overlap_com_cvg(i)(j,jp)=parallel_sum(overlap_com_cvg(i)(j,jp));
    }
  }



  for(int j=0; j< sys->nelectrons(0); j++) {
    for(int jp=0; jp< sys->nelectrons(0);jp++) {
      overlap_matrix(0)(j,jp)/=sqrt(norm(0)(j)*norm(0)(jp));
      overlap_com(0)(j,jp)/=sqrt(norm(0)(j)*norm(0)(jp));
      for(int i=0; i < ncvg; i++) 
	overlap_com_cvg(i)(j,jp)/=sqrt(norm(0)(j)*norm(0)(jp));
    }
  }


  Array1 <doublevar> totxavg(3,0.0);
  Array1 <doublevar> totcharge(nions,0.0);
  for(int j=0; j< sys->nelectrons(0); j++) {
    for(int d=0; d< 3; d++) xavg(j,d)/=norm(0)(j);
    output << "norm " << norm(0)(j) 
	 << " overlap for this orbital " << overlap_matrix(0)(j,j)<<endl;
    output << "xavg ";
    for(int d=0; d< 3; d++) {
      output << xavg(j,d) <<" ";
      totxavg(d)+=xavg(j,d)/sys->nelectrons(0);
    }
    output << endl;
    output << "ion charges ";
    for(int at=0; at < nions; at++) {
      ioncharge(j,at)/=norm(0)(j);
      ioncharge(j,at)*=2.0; //account for double occupation
      output << ioncharge(j,at) << "   ";
      totcharge(at)+=ioncharge(j,at);
    }
    output << endl;
    
  }

  output << "totcharge " << endl;
  for(int at=0; at  < nions; at++) output << totcharge(at) << " ";
  output << endl;


  output << "totxavg ";
  for(int d=0; d< 3; d++) output << totxavg(d) <<" ";
  output << endl;

  doublevar cellVolume=Determinant(latvec,3);
  doublevar au2Cm=57.216;
  
  Array1 <doublevar> epol(3,0.0);
  for(int d=0; d< 3; d++) 
    epol(d)=(sys->nelectrons(0)+sys->nelectrons(1))*totxavg(d)/cellVolume;


  Array1 <doublevar> estpol(3,0.0);
  Array1 <doublevar> ipol(3,0.0);
  Array1 <doublevar> pos(3);


  for(int at=0; at < nions; at++) {
    mywalker->getIonPos(at,pos);
    doublevar charge=mywalker->getIonCharge(at);
    for(int d=0; d< 3; d++) {
      ipol(d)+=charge*pos(d)/cellVolume;
      estpol(d)+=(charge-totcharge(at))*pos(d)/cellVolume;
    }
  }
  for(int d=0; d< 3; d++) {
    output << "From xavg: electronic pol" << d << "  " << au2Cm*epol(d)
	 << " ionic " << au2Cm*ipol(d) << "  total: " 
	 << au2Cm*(ipol(d)-epol(d)) << endl;
  }

  for(int d=0; d< 3; d++) {
    output << "pol from estimated charges" << d << "   " << au2Cm*estpol(d)
	 << endl;
  }


  
  string outname=options.runid+".polgraph";
  ofstream polgraph(outname.c_str());
  for(doublevar ii=0; ii< 1.0; ii+=.01) {
    dcomplex zpolep(0.0,0.0);
    dcomplex zpolen(0.0,0.0);
    doublevar pol0=0.0;
    doublevar atshift=ii*latvec(2,2);
    doublevar chargeaccp=0.0;
    doublevar chargeaccn=0.0;
    polgraph << atshift << "  ";
    for(int at=0; at < nions; at++) {
      mywalker->getIonPos(at,pos);
      doublevar charge=mywalker->getIonCharge(at)-totcharge(at);
      //output << "charge " << charge << endl;
      chargeaccp+=mywalker->getIonCharge(at);
      chargeaccn+=totcharge(at);
      pos(2)+=atshift;
      if(pos(2) > latvec(2,2)) pos(2)-=latvec(2,2);
      double atpol=charge*(pos(2))/cellVolume;
      pol0+=atpol;
      polgraph << atpol <<  "  ";
      zpolep+=mywalker->getIonCharge(at)*dcomplex(cos(2*pi*gvec(2,2)*pos(2)),
						  sin(2*pi*gvec(2,2)*pos(2)));
      zpolen+= totcharge(at)*dcomplex(cos(2*pi*gvec(2,2)*pos(2)),
						  sin(2*pi*gvec(2,2)*pos(2)));
      
    }
    zpolep/=chargeaccp;
    zpolen/=chargeaccn;
    //out << "chargeacc " << chargeaccp << "  " << chargeaccn << endl;
    polgraph << pol0 << "  " << zpolep.real() << "  " << zpolep.imag() 
	     << "  " << zpolen.real() << "   " << zpolep.imag() << " \n";
    if(ii==0) { 
      output << "zpolep " << zpolep << " zpolen  " << zpolen << endl;
      doublevar phasep=atan2(zpolep.imag(), zpolep.real());
      doublevar phasen=atan2(zpolen.imag(), zpolen.real());
      output << "phases: p " << phasep << "  n " << phasen << endl;
      output << "poldiff " << 40*au2Cm*(phasep-phasen)/(2*pi*latvec(0,0)*latvec(1,1))
	   << endl;
      
    }
    
  }
  polgraph.close();


  dcomplex zpol=Determinant(overlap_matrix(0),sys->nelectrons(0));
  output << "zpol " << zpol << endl;
  zpol=zpol*zpol;
  output << " squared " << zpol << endl;

  
  dcomplex zpol_com=Determinant(overlap_com(0),sys->nelectrons(0));
  output << "zpol_com " << zpol_com << endl;
  zpol_com*=zpol_com;
  output << " com squared " << zpol_com << endl;


  for(int i=0; i< ncvg; i++) {
    dcomplex zpol_cvg=Determinant(overlap_com_cvg(i), sys->nelectrons(0));
    zpol_cvg*=zpol_cvg;
    doublevar angle=atan2(zpol_cvg.imag(), zpol_cvg.real());
    doublevar angscale=cvg(i)*angle;
    while(angscale> 2*pi) angscale-=2*pi;
    while(angscale< 0) angscale+=2*pi;
    output << "zpol_cvg " << cvg(i) << " com " << zpol_cvg  
	 << " angle " << angle
	 << " angscale " << angscale
	 << " mag2 " 
	 << zpol_cvg.real()*zpol_cvg.real()+zpol_cvg.imag()*zpol_cvg.imag() 
	 <<endl;
  }

  /*
  int npartoverlap=sys->nelectrons(0);
  Array1 <Array2 <dcomplex> > overlap_partial(npartoverlap);
  Array1 <int> npartial(npartoverlap);
  for(int i=0;i< npartoverlap; i++) {
    overlap_partial(i).Resize(sys->nelectrons(0),sys->nelectrons(0));
    overlap_partial(i)=dcomplex(0.0,0.0);
  }
  for(int e=0; e< sys->nelectrons(0); e++) {
    npartial(e)=e+1;
  }


  for(int i=0; i< npartoverlap; i++) {
    for(int j=0; j< sys->nelectrons(0); j++) {
      for(int jp=0; jp < npartial(i); jp++) 
	overlap_partial(i)(j,jp)=overlap_matrix(0)(j,jp);
    }
    for(int j=npartial(i); j< sys->nelectrons(0); j++) 
      overlap_partial(i)(j,j)=1.0;
    dcomplex zpol_partial=Determinant(overlap_partial(i),sys->nelectrons(0));
    cout << "**for " << npartial(i) << " partial sums  " << endl;
    cout << "partial zpol(n= " << npartial(i) << " : " << zpol_partial 
     << endl;
    doublevar angle=(sys->nelectrons(0)+sys->nelectrons(1))
      *atan2(zpol_partial.imag(), zpol_partial.real())/npartial(i);
    while(angle > pi) { angle-=2*pi;}
    while(angle < -pi) angle+=2*pi;
    cout << "angle " << angle
	 << endl;

    doublevar area=latvec(0,0)*latvec(1,1);
    doublevar epol=au2Cm*angle/(2*pi*area);
    cout << "epol  " << epol << " totalp " << ipol(2)-epol << endl;
    

  }
  */

  
  dcomplex trace(0.0,0.0);
  for(int e=0; e< sys->nelectrons(0); e++) {
    trace+=overlap_matrix(0)(e,e);
  }

  output << "trace " << trace/(double(sys->nelectrons(0))) << endl;
  
  dcomplex twobody(0.0,0.0);
  double npoints=0;
  for(int e=0; e< sys->nelectrons(0); e++ ) {
    for(int j=e+1; j< sys->nelectrons(0); j++) {
      twobody+=overlap_matrix(0)(e,e)*overlap_matrix(0)(j,j)
	-overlap_matrix(0)(e,j)*overlap_matrix(0)(e,j);
      npoints++;
    }
  }
  twobody/=double(npoints);
  twobody*=twobody;
  output << "twobody " << twobody << endl;

  /*
  for(int e=0; e< sys->nelectrons(0); e++) {
    string densname=options.runid+".dens";
    append_number(densname,e);
    ofstream densout(densname.c_str());
    for(int zz=0; zz<D_array1(2);zz++){
      doublevar zpos=zz*resolution_array(2);
      densout << zpos << "   " << density(e,zz) << endl;
    }
  }
  */

}


/*!
Print information about private variables {orbs,resolution,minmax}
*/
int Onebody_int_method::showinfo(ostream & os)
{
  os<<"#############Onebody_int_method#################\n";
  sys->showinfo(os);
  os<<"resolution="<<resolution<<endl;
  os<<"done"<<endl;
  return 1;
}
