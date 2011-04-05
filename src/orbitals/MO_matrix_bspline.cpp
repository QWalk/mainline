/*
 
Copyright (C) 2007 Jindrich Kolorenc, Michal Bajdich

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



#include "Qmc_std.h"
#include "MO_matrix_bspline.h"
#include "Sample_point.h"
#include "qmc_io.h"

//--------------------------------------------------------------------------

void MO_matrix_bspline::init() {

  single_write(cout, "Bspline MO\n");
  MultiEinsplineComplex=new MultiEinsplineOrbComplex;
  MultiEinspline=MultiEinsplineComplex;
  MultiEinspline->kpoint(kvectors_linear);
  MultiEinspline->read(bandfiles,nsplines);
}


//----------------------------------------------------------------------


void MO_matrix_bspline::buildLists(Array1 < Array1 <int> > & occupations) {
  int numlists=occupations.GetDim(0);
  moLists.Resize(numlists);
  for(int lis=0; lis < numlists; lis++) {
    //cout <<"lis # "<<lis<<endl;
    int nmo_list=occupations(lis).GetDim(0);
    moLists(lis).Resize(nmo_list);
    for(int i=0; i < occupations(lis).GetDim(0); i++) {
      int mo=occupations(lis)(i);
      moLists(lis)(i)=mo;
      //cout <<moLists(lis)(i)<<"  ";
    }
    //cout <<endl;
  }
  
}

//----------------------------------------------------------------------

int MO_matrix_bspline::showinfo(ostream & os)
{
  os << "Bspline Molecular Orbital\n";
  string indent="  ";
  os << "Using orbital value plotfiles : \n";
  //Array1 <doublevar> origin, box_size, spacing;

  for(int i=0;i<bandfiles.size();i++){
    os << indent << bandfiles[i] <<endl;
  }

  os << "Number of molecular orbitals: " << nmo << endl;
  os << indent << "kvector for each spline in reduced (qwalk) coordinates" <<endl;
  for(int i=0;i<nsplines;i++)
    os << indent <<i+1<<" : "<< kvectors_linear(i)(0)/pi << "  "
       << kvectors_linear(i)(1)/pi << "  "
       << kvectors_linear(i)(2)/pi<<endl; 
  os <<endl;
  return 1;
}

int MO_matrix_bspline::writeinput(string & indent, ostream & os)
{
  os << indent << "BSPLINE_MO" << endl;
  os << indent << "NMO " << nmo << endl;
  os << indent << "MAGNIFY " << magnification_factor << endl;
  string indent2=indent+"  ";
  string indent3=indent2+"  ";
  for(int k=0;k<kvectors.GetSize();k++){
    os << indent2 << "KVECTOR {" <<endl;
    //Array1 <doublevar> tmpkpoint(3);
    //tmpkpoint=0;
    //for(int d1=0;d1<3;d1++)
    //for(int d2=0;d2<3;d2++)
    //tmpkpoint(d1)+=PrimRecipLatVec(d1,d2)*kvectors(k)(d2);
    //os << indent2 << tmpkpoint(0) << "  "
    //  << tmpkpoint(1) << "  "
    //  << tmpkpoint(2)<<endl; 

    os << indent2 << kvectors(k)(0)/pi << "  "
       << kvectors(k)(1)/pi << "  "
       << kvectors(k)(2)/pi<<endl; 

    for(int b=0;b<bands_per_kvectors[k].size();b++){
      os << indent3 << "BAND {" <<endl;
      os << indent3 <<bands_per_kvectors[k][b][0]<<endl;
      os << indent3 <<bands_per_kvectors[k][b][1]<<endl;
      os << indent3 << "}"<<endl;
    }//band
    os << indent2 << "}"<<endl;
  }//k
  return 1;
}

void MO_matrix_bspline::read(vector <string> & words, unsigned int & startpos, System * sys){
   unsigned int pos=startpos;

   if(!readvalue(words, pos, nmo, "NMO"))
     {
       error("Need NMO in molecular orbital section");
     }
   if(nmo%2){
     single_write(cout," need even number of molecular orbitals, adding one!  \n");
     nmo++;
   }
   nsplines=int(nmo/2);
   

   if(nmo > 40000) 
     error("You have entered more than 40,000 for NMO.  This seems a bit big; we most likely"
	   " can't handle it.");
   
   
   pos=0;
   if(!readvalue(words, pos, magnification_factor, "MAGNIFY")) {
     magnification_factor=1;
   }
   pos=0;

   origin.Resize(3);
   LatVec.Resize(3,3);
   
   if(sys->getBounds(LatVec)){
     single_write(cout," Using periodic");
     if(!sys->getRecipLattice(RecipLatVec))
       error("Could not read Reciprocal Lattice Vectors from system");
     sys->kpoint(kpoint);
     for(int d=0;d<3;d++){
       kpoint(d)*=pi;
     }
     kpoint_square=dot(kpoint,kpoint);
     origin.Resize(3);
     sys->getorigin(origin);
     single_write(cout," and complex splines\n");

     PrimRecipLatVec.Resize(3,3);
     PrimLatVec.Resize(3,3);
     if(!sys->getPrimLattice(PrimLatVec))
       error("Could not read Prim Lattice Vectors from system");
     else{
       Array1 <doublevar> a(3);
       Array1 <doublevar> b(3);
       Array1 <doublevar> ratio(3);
       Array1 <doublevar> newratio(3);
       Array1 <doublevar> fractpart(3); 
       for(int d1=0;d1<3;d1++){
	 for(int d2=0;d2<3;d2++){
	   a(d2)=LatVec(d1,d2);
	   b(d2)=PrimLatVec(d1,d2);
	 }
	 ratio(d1)=dot(a,b)/dot(b,b);
	 //cout << modf (ratio(d1), &newratio(d1))<<"  "<<newratio(d1) <<endl;
	 if(fabs(modf (ratio(d1), &newratio(d1)))>1e-06)
	   error("Simulation cell lattice has to be integer multiple of the primitive lattice\n");
       }
       single_write(cout,"The simulation cell lattice is ",newratio(0),"x");
       single_write(cout,newratio(1),"x",newratio(2)," of the primitive cell\n");
     }
       
     if(!sys->getPrimRecipLattice(PrimRecipLatVec))
       error("Could not read Prim Reciprocal Lattice Vectors from system");
     
     /*
     single_write(cout,"PrimLatVec: \n");
     for(int d1=0;d1<3;d1++){
       for(int d2=0;d2<3;d2++)
	 single_write(cout,PrimLatVec(d1,d2)," ");
       single_write(cout,"\n");
     }
     single_write(cout,"PrimRecipLatVec: \n");
     for(int d1=0;d1<3;d1++){
       for(int d2=0;d2<3;d2++)
	 single_write(cout,PrimRecipLatVec(d1,d2)," ");
       single_write(cout,"\n");
     }
     */


     
     pos=0;
     vector < vector <string> > strkpoints;
     vector <string> tmpstrkpoints;
     kvectors_linear.Resize(nsplines);

     while(readsection(words, pos, tmpstrkpoints, "KVECTOR"))
       strkpoints.push_back(tmpstrkpoints);


     kvectors.Resize(strkpoints.size());
     int kounter=0;
     for(int k=0;k<kvectors.GetSize();k++){
       vector < vector <string> > strband2;
       Array1 <doublevar> tmp_kpoints(3);
       if(strkpoints[k].size()<3)
	 error("Needs 3 numbers for each KVECTOR \n");
       
       for(int d=0;d<3;d++)
	 tmp_kpoints(d)=atof(strkpoints[k][d].c_str());

       kvectors(k).Resize(3);
       
       /* this when using 1/a.u. units for K vector (e.g. like Abinit)
       Array1 <doublevar> tmp2_kpoints(3);
       tmp2_kpoints=0.0;
       for(int d1=0;d1<3;d1++)
	 for(int d2=0;d2<3;d2++)
	   tmp2_kpoints(d1)+=PrimLatVec(d1,d2)*tmp_kpoints(d2);
       
       kvectors(k)=tmp2_kpoints;
       */
       for(int d1=0;d1<3;d1++){
	 tmp_kpoints(d1)*=pi;
       }
       kvectors(k)=tmp_kpoints;
              
       unsigned int newpos=2;
       vector <string> strband;
       while(readsection(strkpoints[k], newpos, strband, "BAND")){
	 if(strband.size()==2){
	   bandfiles.push_back(strband[0]);
	   bandfiles.push_back(strband[1]);
	 }
	 else{
	   error("Need 2 cube file name for each BAND section");
	 }
	 strband2.push_back(strband);
	 kvectors_linear(kounter).Resize(3);
	 kvectors_linear(kounter)=kvectors(k);
	 kounter++;
       }
       if(!strband2.size()){
	 error("Needs BAND section inside KVECTOR \n");
       }
       bands_per_kvectors.push_back(strband2);
     }//k
   }//end of periodic=1
   else{
     error("Needs to be periodic system");
   }

   init();
}

//------------------------------------------------------------------------

void MO_matrix_bspline::updateVal(Sample_point * sample, int e,
                                   int listnum,
                                   //!< which list to use
                                   Array2 <doublevar> & newvals
                                   //!< The return: in form (MO, val)
) {
  //cout <<"MO_matrix_bspline::updateVal"<<endl;
  newvals=0;
  int nmo_list=moLists(listnum).GetDim(0);
  Array1 <doublevar> xyz(3);
  sample->getElectronPos(e,xyz);
  Array1 <doublevar> Unit_pos_in_prim_latice(3);
  Array1 <doublevar> pos_in_prim_latice(3);
  //cout<<"Unit_pos_in_prim_latice:  "<<endl;
  for(int i=0;i<3;i++){
    Unit_pos_in_prim_latice(i)=0;
    pos_in_prim_latice(i)=0;
    for(int j=0;j<3;j++)
      Unit_pos_in_prim_latice(i)+=(xyz(j)-origin(j))*PrimRecipLatVec(i,j);
    pos_in_prim_latice(i)=Unit_pos_in_prim_latice(i);
    Unit_pos_in_prim_latice(i)-= std::floor (Unit_pos_in_prim_latice(i));
    //cout << Unit_pos_in_prim_latice(i)<<"  ";
  }
  //cout <<endl;
  Array1 <doublevar> multi_val(nmo);
  
  Array1 <dcomplex> cval(nsplines);
  //evaluate spline
  MultiEinspline->evaluate_spline (Unit_pos_in_prim_latice, cval);
  
  //loop over all the splines
  for(int m=0; m < nsplines; m++) {
    //cout <<"cval  "<<cval(m)<<endl;
    //multiply by e(ipikr) to get original values 
    dcomplex eikr=exp(I*dot(pos_in_prim_latice,  kvectors_linear(m)));
    cval(m) *= eikr;
    multi_val(2*m)=cval(m).real();
    multi_val(2*m+1)=cval(m).imag();
  }

  for(int m=0; m < nmo_list; m++) {
    int mo=moLists(listnum)(m);
    newvals(m,0)= magnification_factor*multi_val(mo);
  }//m
}

//------------------------------------------------------------------------


void MO_matrix_bspline::updateLap(Sample_point * sample, int e,
				  int listnum,
				  //!< which list to use
				  Array2 <doublevar> & newvals
				  //!< The return: in form (MO, [val, grad, lap])
) {
  //cout <<"MO_matrix_bspline::updateLap"<<endl;
  newvals=0;
 
  Array1 <doublevar> xyz(3);
  int nmo_list=moLists(listnum).GetDim(0);
  sample->getElectronPos(e,xyz);
  Array1 <doublevar> Unit_pos_in_prim_latice(3);
  Array1 <doublevar> pos_in_prim_latice(3);
  //cout<<"Unit_pos_in_prim_latice:  "<<endl;
  //calculate coordinates in the fractions of Latice vectors 
  for(int i=0;i<3;i++){
    Unit_pos_in_prim_latice(i)=0;
    pos_in_prim_latice(i)=0;
    for(int j=0;j<3;j++){
      Unit_pos_in_prim_latice(i)+=(xyz(j)-origin(j))*PrimRecipLatVec(i,j);
    }
    pos_in_prim_latice(i)=Unit_pos_in_prim_latice(i);
    //cout << Unit_pos_in_prim_latice(i)<<"  ";
    Unit_pos_in_prim_latice(i)-= std::floor (Unit_pos_in_prim_latice(i));
    //cout << Unit_pos_in_prim_latice(i)<<"  ";
  }
  //cout <<endl;
  Array1 <doublevar> multi_val(nmo),multi_lap(nmo);
  Array1 < Array1 <doublevar> > multi_grad(nmo);
  for(int m=0; m < nmo; m++){
    multi_grad(m).Resize(3);
    multi_grad(m)=0;
  }
  multi_val=0;
  multi_lap=0;

  Array1 <dcomplex> cval_spline(nsplines);
  Array1 < Array1 <dcomplex> > cgrad_spline(nsplines);
  Array1 < Array2 <dcomplex> > chess_spline(nsplines);
  for(int m=0; m < nsplines; m++){
    cgrad_spline(m).Resize(3);
    chess_spline(m).Resize(3,3);
  }
  MultiEinspline->evaluate_spline (Unit_pos_in_prim_latice, cval_spline, cgrad_spline, chess_spline);
  //multiply by e(ipikr) to get original values 
  for(int m=0; m < nsplines; m++) {
    dcomplex eikr=exp(I*dot(pos_in_prim_latice,  kvectors_linear(m)));
    //cout <<m<<" cval_spline "<<cval_spline(m)<<" cgrad_spline(0) "<<cgrad_spline(m)(0)<<" chess_spline(0) "<<chess_spline(m)(0,0)<<endl;
    dcomplex cval = eikr*cval_spline(m);
    Array1 <dcomplex> cgrad(3),xgrad(3);
    Array2 <dcomplex> chess(3,3);
    for(int d=0;d<3;d++)
      cgrad(d)= eikr*(I*cval_spline(m)*kvectors_linear(m)(d) + cgrad_spline(m)(d));
    //adjust for coordinate change
    xgrad=dcomplex(0,0);
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        xgrad(i)+=cgrad(j)*PrimRecipLatVec(j,i);
      }
    }
    for(int i=0;i<3;i++)
      for(int j=i;j<3;j++){
        chess(i,j)= eikr*( - kvectors_linear(m)(i)* kvectors_linear(m)(j)*cval_spline(m) + 
			   I*kvectors_linear(m)(i)*cgrad_spline(m)(j)+
			   I*kvectors_linear(m)(j)*cgrad_spline(m)(i)+ 
         chess_spline(m)(i,j));
        if(i!=j)
          chess(j,i)=chess(i,j);
      }
    dcomplex clap=0;
    for(int i=0;i<3;i++)
      for(int k=0;k<3;k++)
        for(int l=0;l<3;l++){
          clap+=PrimRecipLatVec(k,i)*PrimRecipLatVec(l,i)*chess(k,l);
        }
      
    multi_val(2*m)=cval.real();
    multi_val(2*m+1)=cval.imag();
    for(int d=0;d<3;d++){
      multi_grad(2*m)(d)=xgrad(d).real();
      multi_grad(2*m+1)(d)=xgrad(d).imag();
    }
    multi_lap(2*m)=clap.real();
    multi_lap(2*m+1)=clap.imag();
  }//m
 
  for(int m=0; m < nmo_list; m++) {
    int mo=moLists(listnum)(m);
    newvals(m,0)= magnification_factor*multi_val(mo);
      for(int d=0;d<3;d++)
        newvals(m,d+1)= magnification_factor*multi_grad(mo)(d);
      newvals(m,4)=magnification_factor*multi_lap(mo);
  }//m
}

//--------------------------------------------------------------------------


void MO_matrix_bspline::updateHessian(
  Sample_point * sample,
  int e,
  int listnum,
  //const Array1 <int> & occupation,
  //!<A list of the MO's to evaluate
  Array2 <doublevar> & newvals
  //!< The return: in form (MO, [val, grad, dxx,dyy,...])
)
{
  
   //cout <<"MO_matrix_bspline::updateHess"<<endl;
  newvals=0;
  Array1 <doublevar> xyz(3);
  int nmo_list=moLists(listnum).GetDim(0);
  sample->getElectronPos(e,xyz);
  Array1 <doublevar> Unit_pos_in_prim_latice(3);
  Array1 <doublevar> pos_in_prim_latice(3);
  //cout<<"Unit_pos_in_prim_latice:  "<<endl;
  //calculate coordinates in the fractions of Latice vectors 
  for(int i=0;i<3;i++){
    Unit_pos_in_prim_latice(i)=0;
    pos_in_prim_latice(i)=0;
    for(int j=0;j<3;j++){
      Unit_pos_in_prim_latice(i)+=(xyz(j)-origin(j))*PrimRecipLatVec(i,j);
    }
    pos_in_prim_latice(i)=Unit_pos_in_prim_latice(i);
    //cout << Unit_pos_in_prim_latice(i)<<"  ";
    Unit_pos_in_prim_latice(i)-= std::floor (Unit_pos_in_prim_latice(i));
    //cout << Unit_pos_in_prim_latice(i)<<"  ";
  }
  //cout <<endl;
  Array1 <doublevar> multi_val(nmo),multi_lap(nmo);
  Array1 < Array1 <doublevar> > multi_grad(nmo);
  Array1 < Array2 <doublevar> > multi_hess(nmo);
  for(int m=0; m < nmo; m++){
    multi_grad(m).Resize(3);
    multi_hess(m).Resize(3,3);
    multi_grad(m)=0;
    multi_hess(m)=0;
  }
  multi_val=0;

  Array1 <dcomplex> cval_spline(nsplines);
  Array1 < Array1 <dcomplex> > cgrad_spline(nsplines);
  Array1 < Array2 <dcomplex> > chess_spline(nsplines);
  for(int m=0; m < nsplines; m++){
    cgrad_spline(m).Resize(3);
    chess_spline(m).Resize(3,3);
  }
  MultiEinspline->evaluate_spline (Unit_pos_in_prim_latice, cval_spline, cgrad_spline, chess_spline);
  //multiply by e(ipikr) to get original values 
  for(int m=0; m < nsplines; m++) {
    dcomplex eikr=exp(I*dot(pos_in_prim_latice,  kvectors_linear(m)));
    //cout <<m<<" cval_spline "<<cval_spline(m)<<" cgrad_spline(0) "<<cgrad_spline(m)(0)<<" chess_spline(0) "<<chess_spline(m)(0,0)<<endl;
    dcomplex cval = eikr*cval_spline(m);
    Array1 <dcomplex> cgrad(3),xgrad(3);
    Array2 <dcomplex> chess(3,3);
    for(int d=0;d<3;d++)
      cgrad(d)= eikr*(I*cval_spline(m)*kvectors_linear(m)(d) + cgrad_spline(m)(d));
    //adjust for coordinate change
    xgrad=dcomplex(0,0);
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	xgrad(i)+=cgrad(j)*PrimRecipLatVec(j,i);
      }
    }
    for(int i=0;i<3;i++)
      for(int j=i;j<3;j++){
	chess(i,j)= eikr*( - kvectors_linear(m)(i)* kvectors_linear(m)(j)*cval_spline(m) + 
			   I*kvectors_linear(m)(i)*cgrad_spline(m)(j)+
			   I*kvectors_linear(m)(j)*cgrad_spline(m)(i)+ 
			   chess_spline(m)(i,j));
	if(i!=j)
	  chess(j,i)=chess(i,j);
      }
    Array2 < dcomplex > cchess(3,3);
    cchess=dcomplex(0,0);
    for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
	for(int k=0;k<3;k++)
	  for(int l=0;l<3;l++){
	    cchess(i,j)+=PrimRecipLatVec(k,i)*PrimRecipLatVec(l,j)*chess(k,l);
	}
      
    multi_val(2*m)=cval.real();
    multi_val(2*m+1)=cval.imag();
    for(int d=0;d<3;d++){
      multi_grad(2*m)(d)=xgrad(d).real();
      multi_grad(2*m+1)(d)=xgrad(d).imag();
    }

    for(int d1=0;d1<3;d1++){
      for(int d2=0;d2<3;d2++){
	multi_hess(2*m)(d1,d2)=cchess(d1,d2).real();
	multi_hess(2*m+1)(d1,d2)=cchess(d1,d2).imag();
      }
    }
    
  }//m
 
  for(int m=0; m < nmo_list; m++) {
    int mo=moLists(listnum)(m);
    newvals(m,0)= magnification_factor*multi_val(mo);
    for(int d=0;d<3;d++)
      newvals(m,d+1)= magnification_factor*multi_grad(mo)(d);
    int d=4;
    for(int d1=0;d1<3;d1++)
      for(int d2=d1;d2<3;d2++)
	newvals(m,d++)=magnification_factor*multi_hess(mo)(d1,d2);
  }//m
  
}
//--------------------------------------------------------------------------
