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
//src/Jastrow_wf_data.cpp

#include "Qmc_std.h"
#include "qmc_io.h"
#include "Jastrow_wf_data.h"
#include "Wavefunction_data.h"
#include "Jastrow_wf.h"
#include "System.h"
/*!
 
*/
void Jastrow_wf_data::read(vector <string> & words, unsigned int & pos,
                           System * sys)
{


  //cout << "jast wf data " << endl;
  unsigned int startpos=pos;
  Array1 <int> espin(2);
  vector <string> nspinstr;
  espin(0)=sys->nelectrons(0);
  espin(1)=sys->nelectrons(1);

  if(readsection(words, pos, nspinstr, "NSPIN"))
  {
    if(nspinstr.size() != 2)
      error("NSPIN must have 2 elements");
    espin(0)=atoi(nspinstr[0].c_str());
    espin(1)=atoi(nspinstr[1].c_str());
    if(espin(0)+espin(1) != sys->nelectrons(0)+sys->nelectrons(1)) {
      error("NSPIN must specify the same number of electrons as the SYSTEM "
            "in JASTROW.");
    }
  }

  nelectrons=espin(0)+espin(1);

  vector <string> atomlabels;
  sys->getAtomicLabels(atomlabels);
  nions=atomlabels.size();


  //faccp and facco are for compatability with the f90 code.
  //They also help with overflow problems in that nice exponential,
  //although since we never actually exponentiate it, it shouldn't
  //be such a big problem.
  pos=startpos;
  if(!haskeyword(words, pos, "NOFACCO")) {
    if(nions==1)
    {
      facco=nelectrons-1;
      faccp=1.0;
    }
    else if(nions <= 10)
    {
      faccp=1.0/((double(nions-1)/double(nions))*double(nelectrons-1));
      facco=(nelectrons-1)*faccp;
    }
    else
    {
      faccp=1.0/nelectrons;
      facco=1.0;
    }
  }
  else { //if there is the keyword NOFACCO
    faccp=1.0;
    facco=1.0;
  }


  firste.Resize(2);
  laste.Resize(2);
  firste(0)=0;
  firste(1)=espin(0);
  laste(0)=espin(0);
  laste(1)=laste(0)+espin(1);

  spin.Resize(nelectrons);
  for(int s=0; s< 2; s++)
  {
    //cout << firste(s) << "   " << laste(s) << endl;
    for(int e=firste(s); e<laste(s); e++)
    {
      spin(e)=s;
    }
  }



  //----------------------------------------------------------------------
  //Normalization
  pos=startpos;
  if(!readvalue(words,pos, normalization, "NORMALIZATION"))
  {
    normalization=0;
  }


  //---------------------------------------------------------------
  //Basis functions

  pos=startpos;
  vector <string> elecIonWords;
  readsection(words, pos, elecIonWords, "EIBASIS");
  allocate(elecIonWords, elecIonBasis);

  pos=startpos;
  vector <string> elecElecWords;
  readsection(words, pos, elecElecWords, "EEBASIS");
  allocate(elecElecWords, elecElecBasis);


  //Find the maximum cutoff radii
  maxeicutoff=0;
  for(int i=0; i< elecIonBasis->nfunc(); i++ )
  {
    if(elecIonBasis->cutoff(i) > maxeicutoff)
    {
      maxeicutoff=elecIonBasis->cutoff(i);
    }
  }

  maxeecutoff=0;
  for(int i=0; i< elecElecBasis->nfunc(); i++)
  {
    if(elecElecBasis->cutoff(i) > maxeecutoff)
    {
      maxeecutoff=elecElecBasis->cutoff(i);
    }
  }

  //The cutoffs only work if the electron-ion cutoff is
  //less than half of the electron-electron
  //(because b_0=1, it's infinite range, and we have to
  //be able to know that when the electrons are out of
  //range of each other, they're also out of range
  //of the ions.  Otherwise, we have to treat the parameters
  //like 2 2 0 and such separately, which is annoying.)
  if(maxeicutoff*2 > maxeecutoff)
  {
    debug_write(cout, "Forcing Jastrow electron-electron"
                 " cutoff to be twice the electron-ion cutoff. May not"
                 " be as efficient as it could be.\n");
    maxeecutoff=maxeicutoff*2;
  }

  debug_write(cout, "Electron-electron cutoff for Jastrow factor: ",
               maxeecutoff, "\n");
  debug_write(cout, "Electron-ion cutoff for Jastrow factor: ",
               maxeicutoff, "\n");

  //----------------------------------------------------------------------

  pos=startpos;
  if(haskeyword(words,pos, "OPTIMIZEBASIS")) {
    single_write(cout, "Turning on basis set optimization in Jastrow factor. "
                 " Don't be surprised if it takes a long time. \n");
    optimize_basis=1;
    nbasisparms=elecElecBasis->nparms()+elecIonBasis->nparms();
  }
  else {
    optimize_basis=0;
    nbasisparms=0;
  }



  //----------------------------------------------------------------
  //Cusp part of the factor


  pos=startpos;
  vector <string> likeCuspWords, unlikeCuspWords;
  if(!readsection(words, pos, likeCuspWords, "LIKECUSP")) {

    error("Need LIKECUSP section in Jastrow factor.  Try: \n"
          "LIKECUSP {\n"
          "  NULL \n"
          "  EXPONENTIAL_CUSP\n"
          "  GAMMA <first gamma> "
          "  CUSP  0.25\n"
          "}");
  }
  pos=startpos;
  if(!readsection(words, pos, unlikeCuspWords, "UNLIKECUSP") ) {
    error("Need UNLIKECUSP section in Jastrow factor.  Try: \n"
          "UNLIKECUSP {\n"
          "  NULL \n"
          "  EXPONENTIAL_CUSP\n"
          "  GAMMA <second gamma> "
          "  CUSP  0.5\n"
          "}");
  }


  cuspBasis.Resize(2);
  cuspBasis=NULL;
  allocate(likeCuspWords, cuspBasis(0));
  allocate(unlikeCuspWords, cuspBasis(1));



  //----------------------------------------------------------------
  //Electron-electron correlation

  pos=startpos;
  vector <string> eeWords;
  readsection(words, pos, eeWords, "EECORRELATION");
  eeCorrelation.Resize(eeWords.size());
  for(unsigned int i=0; i< eeWords.size(); i++)
  {
    eeCorrelation(i)=atof(eeWords[i].c_str());
  }

  //----------------------------------------------------------------
  //Electron-ion correlation

  pos=startpos;
  vector <vector <string> > eiWords;
  vector <string> eiWords_temp;

  int numparms=0;
  while(readsection(words, pos, eiWords_temp, "EICORRELATION"))
  {
    eiWords.push_back(eiWords_temp);
    numparms+=eiWords_temp.size()-1;
    eiWords_temp.erase(eiWords_temp.begin(), eiWords_temp.end());

  }

  eiCorrelation.Resize(numparms);
  eiCorr_start.Resize(nions);
  eiCorr_end.Resize(nions);
  eiCorr_start=0;
  eiCorr_end=0;
  int totnum=0;
  for(unsigned int i=0; i< eiWords.size(); i++)
  {
    int foundone=0;
    for(int k=0; k< nions; k++)
    {
      if(eiWords[i][0] == atomlabels[k])
      {
        if(eiCorr_end(k) != 0) {
          error("There seem to be several EICORRELATION sections for ",
                atomlabels[k]);
        }
        foundone=1;
        eiCorr_start(k)=totnum;
        eiCorr_end(k)=totnum+eiWords[i].size()-1;
      }
    }
    if(foundone==0)
      error("Couldn't find a matching atom for ", eiWords[i][0],
            " in EICORRELATION.");


    for(unsigned int j=1; j< eiWords[i].size(); j++)
    {
      eiCorrelation(totnum)=atof(eiWords[i][j].c_str());
      totnum++;
    }
  }

  //------------------------------------------------------------------
  //Electron-electron-ion correlation

  pos=startpos;
  vector <vector <string> > eeiWords;
  vector <string> eeiWords_temp;
  numparms=0;
  while(readsection(words, pos, eeiWords_temp, "EEICORRELATION"))
  {
    eeiWords.push_back(eeiWords_temp);
    numparms+=eeiWords_temp.size()-1;
    eeiWords_temp.erase(eeiWords_temp.begin(), eeiWords_temp.end());
  }
  if(numparms%4 !=0)
  {
    error("In EEICORRELATION, the number of parameters should"
          "be divisible by four, but is instead ", numparms);
  }


  nUniqueAtoms=eeiWords.size();

  numparms/=4;

  eeiCorrelation.Resize(numparms);
  c_k.Resize(numparms);
  c_l.Resize(numparms);
  c_m.Resize(numparms);
  eeiCorr_start.Resize(nions);
  eeiCorr_end.Resize(nions);
  eeiCorr_start=0;
  eeiCorr_end=0;
  totnum=0;
  for(unsigned int i=0; i< eeiWords.size(); i++)
  {
    for(int k=0; k< nions; k++)
    {
      //cout << eeiWords[i][0] << "   " << options.atoms[k].name << endl;
      if(eeiWords[i][0] == atomlabels[k])
      {
        if(eeiCorr_end(k) !=0) {
          error("There seem to be several EEICORRELATION sections for ",
                atomlabels[k]);
        }
        eeiCorr_start(k)=totnum;
        eeiCorr_end(k)=totnum+(eeiWords[i].size()-1)/4;
        if((eeiWords[i].size()-1)%4 != 0)
          error("need divisable by four in ", atomlabels[k]);
      }
    }

    for(unsigned int j=1; j< eeiWords[i].size(); j++)
    {
      c_k(totnum)=atoi(eeiWords[i][j].c_str());
      j++;
      c_l(totnum)=atoi(eeiWords[i][j].c_str());
      j++;
      c_m(totnum)=atoi(eeiWords[i][j].c_str());
      j++;
      eeiCorrelation(totnum)=atof(eeiWords[i][j].c_str());

      //When k and l are the same, the value a_k*a_l+a_l*a_k
      //should just be a_k*a_k, so we need to prevent double-
      //counting..  Better to do it here than in the inner
      //loop that gets evaluated a bunch of times.
      if(c_k(totnum)==c_l(totnum))
        eeiCorrelation(totnum)*=0.5;

      totnum++;
    }
  }

  //-------------------------------------
  //Temporary variables

  eecusp.Resize(5, nelectrons, nelectrons);

  basisoffset=1;
  bupdate.Resize(nelectrons);

  a_k.Resize(nions, nelectrons);
  //a_kval.Resize(nions, nelectrons);
  for(int i=0; i< nions; i++)
  {
    for(int j=0; j< nelectrons; j++)
    {
      a_k(i,j).Resize(elecIonBasis->nfunc()+1,5);
      //a_kval(i,j).Resize(elecIonBasis->nfunc()+1);

      //The zeroth a_k is defined to be one.
      a_k(i,j)(0,0)=1;
      //a_kval(i,j)(0)=1;
      for(int d=1; d< 5; d++)
      {
        a_k(i,j)(0,d)=0;
      }
    }
  }

  b_m.Resize(nelectrons, nelectrons);
  for(int i=0; i< nelectrons; i++)
  {
    bupdate(i).Resize(elecElecBasis->nfunc()+1);
    bupdate(i)(0)=1;
    for(int j=0; j< nelectrons; j++)
    {
      b_m(i,j).Resize(elecElecBasis->nfunc()+1,5);

      //the zeroth b_m is always 1.
      b_m(i,j)(0,0)=1;
      for(int d=1; d< 5; d++)
      {
        b_m(i,j)(0,d)=0;
      }
    }
  }


  //Keeping track of the names so that we can print out things nicely.
  for(int i=0; i< nions; i++)
  {
    atomname.push_back(atomlabels[i]);
  }

  elecPartial_temp.Resize(nelectrons);

}

int Jastrow_wf_data::supports(wf_support_type support) {
  switch(support) {
  case laplacian_update:
    return 0;
  case density:
    return 1;
  default:
    return 0;
  }
}

//----------------------------------------------------------------------


void Jastrow_wf_data::generateWavefunction(Wavefunction *& wf)
{
  assert(wf==NULL);
  wf=new Jastrow_wf;
  //Jastrow_wf * jaswf;
  //recast(wf, jaswf);
  //jaswf->init(this);
  //attachObserver(jaswf);
  wf->init(this);
  attachObserver(wf);
}

//----------------------------------------------------------------------

void Jastrow_wf_data::renormalize()
{
  Array2 <doublevar> vals(1,2);
  doublevar total=0;
  for(unsigned int i=0; i < wfObserver.size(); i++)
  {
    wfObserver[i]->getVal(this, 0, vals);
    total+=vals(0,1);
  }
  normalization=total/wfObserver.size();
}

void Jastrow_wf_data::resetNormalization()
{
  normalization=0;
}


//----------------------------------------

int Jastrow_wf_data::valSize()
{
  if(optimize_basis) {
    return 0;
  }
  else {
    return eeCorrelation.GetDim(0)+eiCorrelation.GetDim(0)
           +eeiCorrelation.GetDim(0);
  }
}

//----------------------------

int Jastrow_wf_data::showinfo(ostream & os)
{
  os << "Jastrow function\n";

  os << "Normalization: " << normalization << endl;

  os << "\n-------Basis functions----------------\n";
  string indent="  ";
  os << "EIBASIS { \n";
  elecIonBasis->showinfo(indent, os);
  os << " } \n";

  os << "EEBASIS { \n";
  elecElecBasis->showinfo(indent, os);
  os << " } \n";

  os << "--------------Basis Done-------------\n\n";

  //os << "Gamma factor for:\n";
  //os << "Like: " << gamma(0) << "     " << "Unlike: " << gamma(1)
  //   << endl;
  //os << " GAMMA { \n";
  //os << gamma(0) << "     #Like spins \n";
  // os << gamma(1) << "     #Unlike spins \n";
  //os << "} \n";


  os << "\n ---------Cusp basis------------\n";
  os << "Like spins:\n";
  cuspBasis(0)->showinfo(indent, os);
  os << "Unlike spins: \n";
  cuspBasis(1)->showinfo(indent, os);


  //os << "Electron-electron correlation: \n";
  os << "EECORRELATION { \n";
  for(int i=0; i< eeCorrelation.GetDim(0); i++)
  {
    os  << eeCorrelation(i) << endl;
  }
  os << "}\n";

  vector <string> uniquenames;
  for(unsigned int i=0; i< atomname.size(); i++)
  {
    int unique=1;
    for(unsigned int j=0; j< uniquenames.size(); j++)
    {
      if(uniquenames[j]==atomname[i])
      {
        unique=0;
        break;
      }
    }
    if(unique)
    {
      uniquenames.push_back(atomname[i]);
      //os << "Correlation for atom " << atomname[i] << " : \n";
      os << "EICORRELATION { \n";
      os << atomname[i] << endl;
      //cout << "eiCorr_start " << i << "   " << eiCorr_start(i)
      //   <<"   eiCorr_end " << eiCorr_end(i) << endl;
      for(int j=eiCorr_start(i); j< eiCorr_end(i); j++)
      {
        os << eiCorrelation(j) << endl;
      }
      os << "} \n";

      //os << "Electron-electron-ion correlation\n";
      os << "EEICORRELATION { \n";
      os << atomname[i] << endl;
      //cout << i << "  eeiCorr_start  " << eeiCorr_start(i)
      //   << "   eeiCorr_end " << eeiCorr_end(i) << endl;
      for(int j=eeiCorr_start(i); j<eeiCorr_end(i); j++)
      {
        doublevar fac=1;
        if(c_k(j)==c_l(j))
          fac=2.0;
        os << c_k(j) << " " << c_l(j) << " " << c_m(j)
        << "     " << fac*eeiCorrelation(j) << endl;
      }
      os << "}\n";
    }
  }
  os << "########################\n";
  return 1;
}

int Jastrow_wf_data::writeinput(string & indent, ostream & os)
{

  string indent2=indent+"  ";
  os << indent << "JASTROW" << endl;
  os << indent << "NORMALIZATION  " << normalization << endl;

  if(facco==1 && faccp==1) {
    os << indent << "NOFACCO\n";
  }

  os << indent << "EIBASIS { \n";
  elecIonBasis->writeinput(indent2, os);
  os << indent << "} \n";

  os << indent << "EEBASIS { \n";
  elecElecBasis->writeinput(indent2,os);
  os << indent << "} \n";

  os << indent << "LIKECUSP { \n";
  cuspBasis(0)->writeinput(indent2, os);
  os << indent << "} \n";

  os << indent << "UNLIKECUSP { \n";
  cuspBasis(1)->writeinput(indent2, os);
  os << indent << "} \n";
  /*
  os << indent << "CUSP { "
     << cusp(0) << "  " 
     << cusp(1) << " } \n";

  os << indent << "GAMMA { ";
  os << gamma(0) << "  ";
  os << gamma(1);
  os << "} \n";
  */

  //os << "Electron-electron correlation: \n";
  os << indent << "EECORRELATION { ";
  for(int i=0; i< eeCorrelation.GetDim(0); i++)
  {
    os  << eeCorrelation(i) << "  ";
  }
  os << "}\n";

  vector <string> uniquenames;
  for(unsigned int i=0; i< atomname.size(); i++)
  {
    int unique=1;
    for(unsigned int j=0; j< uniquenames.size(); j++)
    {
      if(uniquenames[j]==atomname[i])
      {
        unique=0;
        break;
      }
    }
    if(unique)
    {
      uniquenames.push_back(atomname[i]);
      if(eiCorr_end(i)-eiCorr_start(i) > 0 ) {
        os <<indent <<  "EICORRELATION { \n";
        os <<indent <<  atomname[i] << endl;
        //cout << "eiCorr_start " << i << "   " << eiCorr_start(i)
        //   <<"   eiCorr_end " << eiCorr_end(i) << endl;
        for(int j=eiCorr_start(i); j< eiCorr_end(i); j++)
        {
          os <<indent  << "  " <<  eiCorrelation(j) << endl;
        }
        os << indent << "} \n";
      }

      //os << "Electron-electron-ion correlation\n";
      if(eeiCorr_end(i)-eeiCorr_start(i) > 0 ) {
        os << indent << "EEICORRELATION { \n";
        os << indent <<atomname[i] << endl;
        //cout << i << "  eeiCorr_start  " << eeiCorr_start(i)
        //   << "   eeiCorr_end " << eeiCorr_end(i) << endl;
        for(int j=eeiCorr_start(i); j<eeiCorr_end(i); j++)
        {
          doublevar fac=1;
          if(c_k(j)==c_l(j))
            fac=2.0;
          os << indent << "  "
          << c_k(j) << " " << c_l(j) << " " << c_m(j)
          << "     " << fac*eeiCorrelation(j) << endl;
        }
        os << indent << "}\n";
      }
    }
  }
  return 1;
}


/*
We pass the log of gamma, to put a constraint on the optimization.  
We don't want to it be negative!
*/
void Jastrow_wf_data::getVarParms(Array1 <doublevar> & parms)
{
  parms.Resize(nparms());
  int k=0;
  for(int i=0; i< cuspBasis.GetDim(0); i++) {
    Array1 <doublevar> cuspbasisparms;
    cuspBasis(i)->getVarParms(cuspbasisparms);
    for(int j=0; j< cuspbasisparms.GetDim(0); j++) {
      parms(k)=cuspbasisparms(j);
      k++;
    }
  }
  /*
  for(int i=0; i< gamma.GetDim(0); i++)
  {
    parms(k)=log(gamma(i));
    k++;
  }
  */
  for(int i=0;
      i< eeCorrelation.GetDim(0);
      i++)
  {
    parms(k)=eeCorrelation(i);
    k++;
  }

  for(int i=0; i< eiCorrelation.GetDim(0); i++)
  {
    parms(k)=eiCorrelation(i);
    k++;
  }

  for(int i=0; i< eeiCorrelation.GetDim(0); i++)
  {
    parms(k)=eeiCorrelation(i);
    if(c_k(i)== c_l(i))
      parms(k)*=2.0; //Shouldn't matter, but f90
    //code does it this way.
    //(same in set..)
    k++;
  }

  //Basis set variational parameters
  if(optimize_basis) {
    Array1 <doublevar> eebasisparms;
    elecElecBasis->getVarParms(eebasisparms);
    for(int i=0; i< eebasisparms.GetDim(0); i++) {
      parms(k)=eebasisparms(i);
      k++;
    }

    Array1 <doublevar> eibasisparms;
    elecIonBasis->getVarParms(eibasisparms);
    for(int i=0; i< eibasisparms.GetDim(0); i++) {
      parms(k)=eibasisparms(i);
      k++;
    }
  }

  assert(k==parms.GetDim(0));
}



void Jastrow_wf_data::setVarParms(Array1 <doublevar> & parms)
{
  assert(parms.GetDim(0)==nparms());

  int k=0;


  /*
  for(int i=0; i< gamma.GetDim(0); i++)
  {
    gamma(i)=exp(parms(k));
    k++;
  }
  */

  for(int i=0; i< cuspBasis.GetDim(0); i++) {
    int cuspmax=cuspBasis(i)->nparms();
    Array1 <doublevar> cuspparms(cuspmax);
    for(int j=0; j< cuspmax; j++) {
      cuspparms(j)=parms(k);
      k++;
    }
    cuspBasis(i)->setVarParms(cuspparms);
  }


  for(int i=0;
      i< eeCorrelation.GetDim(0);
      i++)
  {
    eeCorrelation(i)=parms(k);
    k++;
  }

  for(int i=0; i< eiCorrelation.GetDim(0); i++)
  {
    eiCorrelation(i)=parms(k);
    //cout << "eiCorrelation " << i << "   " << eiCorrelation(i) << endl;
    k++;
  }

  for(int i=0; i< eeiCorrelation.GetDim(0); i++)
  {
    eeiCorrelation(i)=parms(k);
    if(c_k(i)==c_l(i))
      eeiCorrelation(i)*=0.5;
    //cout << "eeiCorrelation " << i << "   " << eeiCorrelation(i) << endl;
    k++;
  }


  //Basis set variational parameters
  if(optimize_basis) {
    int eemax=elecElecBasis->nparms();
    Array1 <doublevar> eeparms(eemax);
    for(int i=0; i< eemax; i++) {
      eeparms(i)=parms(k);
      k++;
    }
    elecElecBasis->setVarParms(eeparms);


    int eimax=elecIonBasis->nparms();
    Array1 <doublevar> eiparms(eimax);
    for(int i=0; i< eimax; i++) {
      eiparms(i)=parms(k);
      k++;
    }
    elecIonBasis->setVarParms(eiparms);
  }


  assert(k==parms.GetDim(0));


  int max=wfObserver.size();
  for(int i=0; i< max; i++)
  {
    wfObserver[i]->notify(all_wf_parms_change,0);
  }

}


//------------------------------------------------------------------------
