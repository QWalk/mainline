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

#include "Qmc_std.h"
#include "qmc_io.h"
#include "Slat_wf_data.h"
#include "Wavefunction_data.h"
#include "Slat_wf.h"
#include "Cslat_wf.h"
#include <algorithm>



/*!
*/
void Slat_wf_data::read(vector <string> & words, unsigned int & pos,
                        System * sys)
{


  vector <string> strdetwt;
  vector <string> strstates;
  vector <vector <string> > statevec;
  unsigned int startpos=pos;
  readsection(words, pos, strdetwt, "DETWT");
  pos=startpos;

  while(readsection(words, pos, strstates, "STATES"))
  {
    statevec.push_back(strstates);
  }


  if(readvalue(words, pos=startpos, mo_place, "OPTIMIZE_MO") ) {
    optimize_mo=1;
  }
  else optimize_mo=0;

  if(haskeyword(words, pos=startpos, "OPTIMIZE_DET")) {
    optimize_det=1;
    if(optimize_mo) error("Can't optimize both mo and det right now");
  }
  else optimize_det=0;

  pos=startpos;
  vector <string> mowords;
  if(readsection(words,pos, mowords, "ORBITALS"))  {
    
    allocate(mowords, sys, molecorb);
    
    nmo=molecorb->getNmo();
    genmolecorb=molecorb;
    use_complexmo=0;
  }
  else if(readsection(words, pos=0, mowords, "CORBITALS")) {
    allocate(mowords, sys, cmolecorb);
    nmo=cmolecorb->getNmo();
    genmolecorb=cmolecorb;
    use_complexmo=1;
    if(optimize_mo) error("Can't optimize MO's with complex MO's");
    if(optimize_det) 
      error("Don't support optimizing determinants with complex MO's yet");
  }
  else {
    error("Need ORBITALS or CORBITALS section in SLATER wave function");
  }

  if(optimize_mo) {
    pos=startpos;
    vector <string> orbitals_for_optimize_mostr;
    if(!readsection(words, pos, orbitals_for_optimize_mostr, "ORBITALS_FOR_OPTIMIZE_MO")){
    cout << "Assuming you want to use all orbitals up to NMO in OPTIMIZE_MO section \n";
     orbitals_for_optimize_mo.Resize(nmo);
     for (int i=0;i<nmo;i++){
        orbitals_for_optimize_mo(i)=i;
     }
    }
    else {
      orbitals_for_optimize_mo.Resize(orbitals_for_optimize_mostr.size());
      for (unsigned int i=0;i<orbitals_for_optimize_mostr.size();i++){
	orbitals_for_optimize_mo(i)=atoi(orbitals_for_optimize_mostr[i].c_str())-1;
	if(orbitals_for_optimize_mo(i)+1 > nmo)
	  error("Suplied incorrect orbital number in  ORBITALS_FOR_OPTIMIZE_MO");
      }
    }
    cout << "orbitals_for_optimize_mo: ";
    for(int i=0;i<orbitals_for_optimize_mo.GetSize();i++)
      cout << orbitals_for_optimize_mo(i) <<"  ";
    cout <<endl;
  }
  



  ndet=strdetwt.size();
  nfunc=statevec.size();


  nelectrons.Resize(2);

  nelectrons(0)=sys->nelectrons(0);
  nelectrons(1)=sys->nelectrons(1);
  pos=startpos;
  vector <string> nspinstr;
  if(readsection(words, pos, nspinstr, "NSPIN"))
  {
    if(nspinstr.size() != 2)
      error("NSPIN must have 2 elements");
    nelectrons(0)=atoi(nspinstr[0].c_str());
    nelectrons(1)=atoi(nspinstr[1].c_str());
    if(nelectrons(0)+nelectrons(1) != sys->nelectrons(0)+sys->nelectrons(1)) {
      error("NSPIN must specify the same number of electrons as the SYSTEM "
            "in SLATER.");
    }
  }

  pos=startpos;
  unsigned int canonstates=ndet*(nelectrons(0)+nelectrons(1));
  for(int i=0; i< nfunc; i++)
  {
    if( canonstates != statevec[i].size())
    {
      error("in STATES section, expecting to find ", canonstates,
            " states(as calculated from NSPIN), but found ",
            statevec[i].size(), " instead.");
    }
  }




  int tote=nelectrons(0)+nelectrons(1);
  ndim=3;


  //Input parameters
  occupation.Resize(nfunc, ndet, 2);
  occupation_orig.Resize(nfunc, ndet, 2);
  detwt.Resize(ndet);


  for(int i=0; i< nfunc; i++)
  {
    for(int det=0; det < ndet; det++)
    {
      for(int s=0; s<2; s++)
      {
        //cout << det << " "  << s << endl;
        occupation(i,det,s).Resize(nelectrons(s));
        occupation_orig(i,det,s).Resize(nelectrons(s));
      }
    }
  }
  //Calculation helpers
  for(int i=0; i< nfunc; i++)
  {
    //cout << "i=" << i << endl;
    int counter=0;
    for(int det=0; det<ndet; det++)
    {
      //cout << "det=" << det << endl;
      for(int s=0; s<2; s++)
      {
        //cout << "s=" << s << endl;
        //cout << "nelectrons " << nelectrons(s) << endl;
        for(int e=0; e<nelectrons(s); e++)
        {
          //cout << "e=" << e << endl;
          occupation_orig(i,det,s)(e)=atoi(statevec[i][counter].c_str())-1;

          counter++;
        }
      }
    }
  }


  spin.Resize(tote);
  opspin.Resize(tote);
  rede.Resize(tote);

  int eup=nelectrons(0);
  for(int e=0; e<eup; e++)
  {
    spin(e)=0;
    opspin(e)=1;
    rede(e)=e;
  }
  for(int e=eup; e<tote; e++)
  {
    spin(e)=1;
    opspin(e)=0;
    rede(e)=e-eup;
  }



  //Find what MO's are necessary for each spin
  totoccupation.Resize(2);
  for(int s=0; s<2; s++)
  {
    vector <int> totocctemp;
    for(int f=0; f< nfunc; f++)
    {
      for(int det=0; det<ndet; det++)
      {
        for(int mo=0; mo < nelectrons(s); mo++)
        {
          int place=-1;
          int ntot=totocctemp.size();
          for(int i=0; i< ntot; i++) {
            //cout << "i " << i << endl;
            if(occupation_orig(f,det,s)(mo)==totocctemp[i]) {
              place=i;
              break;
            }
          }
          //cout << "decide " << endl;
          if(place==-1) { //if we didn't find the MO
            //cout << "didn't find it " << endl;
            occupation(f,det,s)(mo)=totocctemp.size();
            totocctemp.push_back(occupation_orig(f,det,s)(mo));
          }
          else {
            //cout << "found it" << endl;
            occupation(f,det,s)(mo)=place;
          }


        }
      }
    }

    //cout << "done assignment " << endl;

    totoccupation(s).Resize(totocctemp.size());
    for(int i=0; i<totoccupation(s).GetDim(0); i++)
    {
      totoccupation(s)(i) = totocctemp[i];
      //cout << "total occupation for " << s<< " : "
      // << totoccupation(s)(i) << endl;
    }
  }

  //molecorb->buildLists(totoccupation);
  if(genmolecorb) init_mo();


  for(int det=0; det < ndet; det++)
  {
    detwt(det)=atof(strdetwt[det].c_str());
  }

}

//----------------------------------------------------------------------

int Slat_wf_data::supports(wf_support_type support) {
  switch(support) {
  case laplacian_update:
    return 1;
  case density:
    return 1;
  case parameter_derivatives:
    return 1;
  default:
    return 0;
  }
}

//----------------------------------------------------------------------

void Slat_wf_data::init_mo() {
  Array1 <int> nmospin(2);
  nmospin=0;

  for(int i=0; i< nfunc; i++)
  {
    //cout << "i=" << i << endl;
    int counter=0;
    for(int det=0; det<ndet; det++)
    {
      //cout << "det=" << det << endl;
      for(int s=0; s<2; s++)
      {
        //cout << "s=" << s << endl;
        //cout << "nelectrons " << nelectrons(s) << endl;
        for(int e=0; e<nelectrons(s); e++)
        {
          if(nmospin(s) < occupation_orig(i,det,s)(e)+1)
            nmospin(s)=occupation_orig(i,det,s)(e)+1;

          counter++;
        }
      }
    }
  }

  if(nmospin(0) > nmo)
    error("First spin channel contains an orbital higher"
          " than requested NMO's.");
  if(nmospin(1) > nmo)
    error("Second spin channel contains an orbital higher"
          " than requested NMO's.");


  genmolecorb->buildLists(totoccupation);

}

//----------------------------------------------------------------------


void Slat_wf_data::generateWavefunction(Wavefunction *& wf)
{
  assert(wf==NULL);
  if(!genmolecorb)
    error("Slat_wf_data::Need to allocate molecular orbital before generating any Wavefunction objects");

  if(use_complexmo) {
    wf=new Cslat_wf;
    Cslat_wf * slatwf;
    recast(wf,slatwf);
    slatwf->init(this);
    attachObserver(slatwf);
  }
  else { 
    wf=new Slat_wf;
    Slat_wf * slatwf;
    recast(wf, slatwf);
    slatwf->init(this);
    attachObserver(slatwf);
  }
}

int Slat_wf_data::showinfo(ostream & os)
{
  if(!genmolecorb)
    error("Slat_wf_data::showinfo() : Molecular orbital not allocated");

  os << "Slater Determinant" << endl;


  if(optimize_mo) 
    os << "Optimizing molecular orbital coefficients" << endl;

  if(ndet > 1)
    os << ndet << " determinants\n";
  else
    os << "1 determinant" << endl;

  for(int f=0; f< nfunc; f++)
  {
    if(nfunc > 1)
      os << "For function " << f << endl;
    for(int det=0; det<ndet; det++)
    {
      if(ndet > 1) {
        os << "Determinant " << det << ":\n";
        os << "Weight: " << detwt(det) << endl;
      }

      os << "State: \n";
      for(int s=0; s<2; s++)
      {
        if(s==0)
          os << "spin up:\n";
        if(s==1)
          os << "spin down: \n";

        os << "  ";
        for(int e=0; e<nelectrons(s); e++)
        {
          os << occupation_orig(f,det,s)(e)+1 << " ";
          if((e+1)%10 == 0)
            os << endl << "  ";
        }
        os << endl;
      }
    }
  }

  os << "Molecular Orbital object : ";
  genmolecorb->showinfo(os);

  return 1;
}

//----------------------------------------------------------------------

int Slat_wf_data::writeinput(string & indent, ostream & os)
{

  if(!genmolecorb)
    error("Slat_wf_data::writeinput() : Molecular orbital not allocated");

  os << indent << "SLATER" << endl;

  if(optimize_mo) {
    os << indent << "OPTIMIZE_MO" << endl;
    os << indent << "ORBITALS_FOR_OPTIMIZE_MO {";
    for(int i=0;i<orbitals_for_optimize_mo.GetSize();i++)
      os << orbitals_for_optimize_mo(i)+1 <<" ";
    os << "}" << endl;	
  }
  if(optimize_det)
    os << indent << "OPTIMIZE_DET" << endl;

  os << indent << "DETWT { ";
  for(int det=0; det < ndet; det++)
  {
    os << detwt(det) << "  ";
  }
  os << "}" << endl;

  for(int f=0; f< nfunc; f++)
  {
    os << indent << "STATES { " << endl << indent <<"  ";
    for(int det=0; det < ndet; det++)
    {
      for(int s=0; s<2; s++)
      {
        for(int e=0; e< nelectrons(s); e++)
        {
          os << occupation_orig(f,det,s)(e)+1 << " ";
          if((e+1)%20 ==0)
            os << endl << indent << "  ";
        }
        os << endl << indent << "  ";
      }
    }
    os << "}" << endl;
  }

  if(optimize_mo) {
    molecorb->setOrbfile(mo_place);

    ofstream tmp(mo_place.c_str());
    Array2 <doublevar> rotation(nmo, nmo);
    Array1 <int> moList(nmo);
    rotation=0;
    for(int i=0; i< nmo; i++) {
      rotation(i,i)=1;
      moList(i)=i;
    }
    
    molecorb->writeorb(tmp, rotation, moList);
    tmp.close();
  }

  if(use_complexmo) os << indent << "CORBITALS { \n";
  else os << indent << "ORBITALS {\n";
  string indent2=indent+"  ";
  genmolecorb->writeinput(indent2, os);
  os << indent << "}\n";

  return 1;
}

//------------------------------------------------------------------------
void Slat_wf_data::getVarParms(Array1 <doublevar> & parms)
{
  //cout <<"start getVarParms"<<endl;
  if(optimize_mo) {
    Array2 <doublevar> mocoeff;
    molecorb->getMoCoeff(mocoeff);
    int totmo=molecorb->getNmo();
    int totcoeff=molecorb->nMoCoeff();
    int coeff_per_orbital=totcoeff/totmo;
    int nmo_to_optimize=orbitals_for_optimize_mo.GetSize();
    parms.Resize(nmo_to_optimize*coeff_per_orbital);
    int counter=0;
    for(int i=0; i< nmo_to_optimize; i++) {
      for(int j=0; j < mocoeff.GetDim(1); j++) {
        parms(counter)=mocoeff(orbitals_for_optimize_mo(i),j);
        counter++;
      }
    }
  }
  else if(optimize_det) {
    parms.Resize(detwt.GetDim(0)-1);
    for(int i=1; i< detwt.GetDim(0); i++) {
      parms(i-1)=detwt(i);
    }
  }
  else {
    parms.Resize(0);
  }
  //cout <<"done getVarParms"<<endl;
}

void Slat_wf_data::setVarParms(Array1 <doublevar> & parms)
{
  //cout <<"start setVarParms"<<endl;
  if(optimize_mo) {
    int totmo=molecorb->getNmo();
    int totcoeff=molecorb->nMoCoeff();
    int coeff_per_orbital=totcoeff/totmo;
    int nmo_to_optimize=orbitals_for_optimize_mo.GetSize();

    assert(parms.GetDim(0)==nmo_to_optimize*coeff_per_orbital);

    Array2 <doublevar> mocoeff;
    molecorb->getMoCoeff(mocoeff);
    int counter=0;
    for(int i=0; i< nmo_to_optimize; i++) {
      for(int j=0; j < mocoeff.GetDim(1); j++) {
        mocoeff(orbitals_for_optimize_mo(i),j)=parms(counter);
        counter++;
      }
    }
    molecorb->setMoCoeff(mocoeff);
  }
  else if(optimize_det) {
    for(int i=1; i< detwt.GetDim(0); i++) {
      detwt(i)=parms(i-1);
    }
  }
  else {
    parms.Resize(0);
  }
  
  int max=wfObserver.size();
  //cout << "slatmax " << max << endl;
  for(int i=0; i< max; i++) {
    wfObserver[i]->notify(all_wf_parms_change, 0);
  }
  //cout <<"done setVarParms"<<endl;
}
