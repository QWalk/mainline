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

#include "Jastrow2_wf.h"
#include "qmc_io.h"
#include "System.h"
#include "Sample_point.h"

//######################################################################


//----------------------------------------------------------------------
//functions to change the input so that the old style Jastrow 3-body 
//factor will work
//This consists of:
// -Adding a step function with range 1e99 to the beginning of all sections
// -Adding a 0 in the beginning of the two-body and one-body coefficient
//  sections

void fix_onebody_txt(vector <string> & onebodywords) {
  
  for(vector<string>::iterator i=onebodywords.begin();
      i!=onebodywords.end(); i++) {
    if(caseless_eq(*i,"COEFFICIENTS")) { 
      i++; //{
      i++; //label
      i++;
      i=onebodywords.insert(i,"0.0");
    }
  }
}
		     
void fix_onebody_basis(vector <vector <string> > & basistxt,
		       vector <string> atomnames) { 
  vector <string> uniqueat;
  int natoms=atomnames.size();
  for(int at=0; at< natoms; at++) { 
    int unique=1;
    for(vector<string>::iterator i=uniqueat.begin();
	i!=uniqueat.end(); i++) {
      if(*i==atomnames[at]) {
	unique=0;
	break;
      }
    }
    if(unique) { 
      vector <string> tmptxt;
      tmptxt.push_back(atomnames[at]);
      tmptxt.push_back("STEP");
      tmptxt.push_back("CUTOFF");
      tmptxt.push_back("1e99");
      basistxt.push_back(tmptxt);
      uniqueat.push_back(atomnames[at]);
    }
  }

}

void fix_twobody_basis(vector < vector <string> > & basistxt) { 
  vector <string> tmptxt;
  tmptxt.push_back("EE");
  tmptxt.push_back("STEP");
  tmptxt.push_back("CUTOFF");
  tmptxt.push_back("1e99");
  basistxt.push_back(tmptxt);
  
}

void fix_twobody_section(vector <string> & twobodysec) {
  
  for(vector<string>::iterator i=twobodysec.begin();
      i!=twobodysec.end(); i++) {
    if(caseless_eq(*i,"COEFFICIENTS")) { 
      i++; //{
      i++;
      i=twobodysec.insert(i,"0.0");
    }
  }
}

//----------------------------------------------------------------------

void Jastrow_group::set_up(vector <string> & words, System * sys) {
  //cout <<"Start set up"<<endl;
  unsigned int pos=0;
  vector <string> onebodywords;
  //vector <string> atomnames;
  sys->getAtomicLabels(atomnames);
  if(haskeyword(words, pos=0, "OPTIMIZEBASIS"))
    optimize_basis=1;
  else optimize_basis=0;


  int fix_old=haskeyword(words, pos=0,"OLD_FORMAT");
  int natoms=atomnames.size();

  has_one_body=has_two_body=0;
  if(readsection(words, pos=0, onebodywords, "ONEBODY")) {
    has_one_body=1;
    if(fix_old) fix_onebody_txt(onebodywords);
    one_body.set_up(onebodywords, atomnames);
  }
  pos=0;
  vector < vector < string> > basistxt;
  vector <string> tmpbasis;

  if(fix_old) fix_onebody_basis(basistxt,atomnames);

  while(readsection(words, pos, tmpbasis, "EIBASIS"))
    basistxt.push_back(tmpbasis);
  int nbasis=basistxt.size();
  eibasis.Resize(nbasis);
  eibasis=NULL;
  for(int b=0; b< nbasis; b++) {
    allocate(basistxt[b], eibasis(b));
  }

  atom2basis.Resize(natoms, nbasis);
  nbasis_at.Resize(natoms);
  nbasis_at=0;

  for(int b=0; b< nbasis; b++) {
    int found=0;
    for(int at=0; at < natoms; at++) {
      if(eibasis(b)->label() == atomnames[at]) {
        atom2basis(at, nbasis_at(at)++)=b;
        //cout << "atom " << at << " -> basis " << b << endl;
        found=1;
      }
    }
    if(!found)
      error("Couldn't find center for basis label ", eibasis(b)->label());
  }


  nfunc_eib.Resize(nbasis);
  for(int b=0; b< nbasis; b++) {
    nfunc_eib(b)=eibasis(b)->nfunc();
  }
  toteibasis=0;
  for(int b=0; b< nbasis; b++) {
    toteibasis+=nfunc_eib(b);
  }
  maxeibasis=0;
  for(int b=0; b< nbasis; b++) {
    if(maxeibasis < nfunc_eib(b)) maxeibasis=nfunc_eib(b);
  }

  maxbasis_on_center=0;
  for(int at=0; at < natoms; at++) {
    int totbas=0;
    for(int n=0; n< nbasis_at(at); n++) {
      int b=atom2basis(at,n);
      totbas+=nfunc_eib(b);
    }
    if(maxbasis_on_center < totbas) maxbasis_on_center=totbas;
  }


  //Two-body stuff
  vector < vector < string> > eebasistxt;
  vector <string> eebasistmp;
  if(fix_old) fix_twobody_basis(eebasistxt);
  pos=0;
  while(readsection(words, pos, eebasistmp, "EEBASIS"))
    eebasistxt.push_back(eebasistmp);

  //cout << eebasistxt.size() << "  e-e basis objects " << endl;
  eebasis.Resize(eebasistxt.size());
  nfunc_eeb.Resize(eebasistxt.size());
  eebasis=NULL;
  for(int b=0; b< eebasis.GetDim(0); b++) {
    allocate(eebasistxt[b],eebasis(b));
    nfunc_eeb(b)=eebasis(b)->nfunc();
  }
  maxeebasis=0;
  n_eebasis=0;
  for(int b=0; b< eebasis.GetDim(0); b++) {
    n_eebasis+=nfunc_eeb(b);
    if(maxeebasis < nfunc_eeb(b)) maxeebasis=nfunc_eeb(b);
  }

  vector <string> twobodysec;
  if(readsection(words, pos=0, twobodysec, "TWOBODY")) {
    has_two_body=1;
    have_diffspin=0;
    if(fix_old) fix_twobody_section(twobodysec);

    two_body=new Jastrow_twobody_piece;
    two_body->set_up(twobodysec, sys);
  }

  if(readsection(words, pos=0, twobodysec, "TWOBODY_SPIN")) {
    if(has_two_body)
      error("Can only have one of TWOBODY_SPIN or TWOBODY");
    if(fix_old) 
      error("Should not use OLD_FORMAT for spin-dependent Jastrows");
    have_diffspin=1;
    has_two_body=1;
    two_body=new Jastrow_twobody_piece_diffspin;
    two_body->set_up(twobodysec, sys);
  }

  nelectrons=sys->nelectrons(0)+sys->nelectrons(1);
  n_spin_up=sys->nelectrons(0);
  has_three_body=0;
  vector<string> threebodysec;
  if(readsection(words, pos=0,threebodysec, "THREEBODY")) { 
    has_three_body=1;
    three_body.set_up(threebodysec, sys);
  } 
  has_three_body_diffspin=0;
  vector<string> threebodysec_diffspin;
  if(readsection(words, pos=0,threebodysec_diffspin, "THREEBODY_SPIN")) { 
    has_three_body_diffspin=1;
    if(has_three_body)
      error("Can only have one of THREEBODY_SPIN or THREEBODY");
    three_body_diffspin.set_up(threebodysec_diffspin, sys);
  } 

  check_consistency();
}

//----------------------------------------------------------------------

int Jastrow_group::check_consistency() {

  //Check to make sure that the parameters make sense
  if(has_one_body) {
    int natoms=atomnames.size();
    vector <string> uniquenames;
    for(int at=0; at < natoms; at++) {
      int unique=1;
      for(unsigned int i=0; i< uniquenames.size(); i++) {
        if(atomnames[at]==uniquenames[i]) {
          unique=0;
          break;
        }
      }
      if(unique) {

        uniquenames.push_back(atomnames[at]);
        int nbas=0;
        for(int n=0; n< nbasis_at(at); n++)
          nbas+=nfunc_eib(atom2basis(at,n));

        if(one_body.nparms(at) > nbas)
          error("Jastrow_group::the number of electron-ion functions is"
              " less than the number of parameters given for atom ", 
              atomnames[at]);

        if(has_three_body) { 
          if(three_body.eibasis_needed(at) > nbas) 
            error("Threebody needs at least ", three_body.eibasis_needed(at),
                " basis functions for atom ", atomnames[at]);
        }
        if(has_three_body_diffspin) { 
          if(three_body_diffspin.eibasis_needed(at) > nbas) 
            error("Threebody needs at least ", three_body_diffspin.eibasis_needed(at),
                " basis functions for atom ", atomnames[at]);
        }
      }
    }
  }


  if(has_two_body) {

    int neebasis_check=two_body->nbasis_needed();
    if(neebasis_check > n_eebasis)
      error("Jastrow_group::Not enough EE basis functions for the number of parameters");

    if(has_three_body) { 
      if(three_body.eebasis_needed() > n_eebasis) 
        error("Threebody needs at least ", three_body.eebasis_needed(),
            " basis functions for the electron-electron interaction");
    }
    if(has_three_body_diffspin) { 
      if(three_body_diffspin.eebasis_needed() > n_eebasis) 
        error("Threebody needs at least ", three_body_diffspin.eebasis_needed(),
            " basis functions for the electron-electron interaction");
    }
  }

  return 1;
}


//----------------------------------------------------------------------

void Jastrow_group::updateEIBasis(int e, Sample_point * sample,
                                  Array3 <doublevar> & eisave) {
  //cout << "updateEIBasis" << endl;
  int natoms=nbasis_at.GetDim(0);

  assert(sample->ionSize() == natoms);
  assert(atom2basis.GetDim(0)==natoms);

  Array1 <doublevar> R(5);
  Array2 <doublevar> lap(maxeibasis, 5);
  sample->updateEIDist();

  int b; //basis

  for(int at=0; at < natoms; at++) {
     sample->getEIDist(e,at, R);
     int counter=0;
     //eisave(at,counter,0)=1;
     //for(int d=1; d< 5; d++) 
     //  eisave(at,counter,d)=0;
     //counter++;

     for(int n=0; n< nbasis_at(at); n++) {
       b=atom2basis(at, n);
       eibasis(b)->calcLap(R, lap);
       for(int i=0; i< nfunc_eib(b); i++) {
        for(int d=0; d< 5; d++)
          eisave(at,counter, d)=lap(i,d);
        counter++;
       }
     }
  }
  //cout << "done" << endl;
}

//----------------------------------------------------------------------

void Jastrow_group::updateEEBasis(int e, Sample_point * sample,
                                  Array3 <doublevar> & eesave) {
  //cout << "updateEEBasis" << endl;
  Array1 <doublevar> R(5);
  Array2 <doublevar> lap(maxeebasis, 5);

  int neebasis=eebasis.GetDim(0);
  int counter=0;
  sample->updateEEDist();

  //for(int i=0; i< nelectrons; i++) { 
  //  eesave(i,counter,0)=1;
  //  for(int d=1; d< 5; d++)
  //    eesave(i,counter,d)=0;
  //}
  //counter++;

  for(int b=0; b< neebasis; b++) {
    doublevar cutoff=0;
    for(int n=0; n< nfunc_eeb(b); n++) {
      if(cutoff < eebasis(b)->cutoff(n)) cutoff=eebasis(b)->cutoff(n);
    }

    for(int i=0; i< e; i++) {
      sample->getEEDist(i,e,R);
      if(R(0) < cutoff) {
        eebasis(b)->calcLap(R,lap);
        for(int n=0; n< nfunc_eeb(b); n++) {
          for(int d=0; d< 5; d++) {
            eesave(i,counter+n, d)=lap(n,d);
            //if(R(0) > cutoff) 
            //  cout << "beyond " << eesave(i,counter+n, d) << endl;
          }
        }
      }
      else {
        for(int n=0; n< nfunc_eeb(b); n++) {
          for(int d=0; d< 5; d++)
            eesave(i,counter+n, d)=0;
        }
      }
    }
    for(int j=e+1; j< nelectrons; j++) {
      sample->getEEDist(e,j,R);
      if(R(0) < cutoff) { 
        eebasis(b)->calcLap(R,lap);
        for(int n=0; n< nfunc_eeb(b); n++) {
          for(int d=0; d< 5; d++) { 
            eesave(j,counter+n, d)=lap(n,d);
          }
        }
      }
      else { 
        for(int n=0; n< nfunc_eeb(b); n++) {
          for(int d=0; d< 5; d++)
            eesave(j,counter+n, d)=0;
        }        
      }
    }
    counter+=nfunc_eeb(b);
  }
  //cout << "done" << endl;


}

//----------------------------------------------------------------------


int Jastrow_group::writeinput(string & indent, ostream & os) {
  string indent2=indent+"  ";

  if(optimize_basis) {
    os << indent << "OPTIMIZEBASIS" << endl;
  }

  if(has_one_body) {
    os << indent << "ONEBODY { " << endl;
    one_body.writeinput(indent2, os);
    os << indent << "} " << endl;
  }
  if(has_two_body) {
    if(have_diffspin)
      os << indent << "TWOBODY_SPIN { " << endl;
    else
      os << indent << "TWOBODY { " << endl;
    two_body->writeinput(indent2, os);
    os << indent << "} " << endl;
  }

  if(has_three_body) { 
    os << indent << "THREEBODY  { " << endl;
    three_body.writeinput(indent2, os);
    os << indent << "}\n";
  }
  if(has_three_body_diffspin) { 
    os << indent << "THREEBODY_SPIN  { " << endl;
    three_body_diffspin.writeinput(indent2, os);
    os << indent << "}\n";
  }
  
  for(int b=0; b< eibasis.GetDim(0); b++) {
    os << indent << "EIBASIS { " << endl;
    eibasis(b)->writeinput(indent2, os);
    os << indent << "}" << endl;
  }

  for(int b=0; b< eebasis.GetDim(0); b++) {
    os << indent << "EEBASIS { " << endl;
    eebasis(b)->writeinput(indent2, os);
    os << indent << "}" << endl;
  }
  return 1;
}

//----------------------------------------------------------------------

int Jastrow_group::showinfo(string & indent, ostream & os) {
  string indent2=indent+"  ";

  if(optimize_basis) {
    os << indent << "Basis optimization turned on" << endl;
  }

  if(has_one_body) {
    os << indent << "One-body terms " << endl;
    one_body.showinfo(indent2, os);
    os << endl;
  }


  if(has_two_body) {
    if(have_diffspin)
      os << indent << "Separate spin two-body" << endl;
    else
      os << indent << "Two-body terms " << endl;
    two_body->showinfo(indent2, os);
  }
  
  if(has_three_body) { 
    os << indent << " Three-body terms " << endl;
    three_body.showinfo(indent2, os);
    os << endl;
  }
  if(has_three_body_diffspin) { 
    os << indent << " Separate spin three-body terms " << endl;
    three_body_diffspin.showinfo(indent2, os);
    os << endl;
  }
  
  if(has_one_body || has_three_body || has_three_body_diffspin)
    os << indent << "Electron-ion basis" << endl;
  for(int b=0; b< eibasis.GetDim(0); b++) {
    eibasis(b)->showinfo(indent2, os);
    os << endl;
  }
  os << endl;

  
  
  os << indent << "Electron-electron basis" << endl;
  for(int b=0; b< eebasis.GetDim(0); b++) {
    eebasis(b)->showinfo(indent2, os);
  }
  os << endl;
  return 1;
}

//----------------------------------------------------------------------

int Jastrow_group::nparms() {
  int tot=0;
  if(has_one_body)
    tot+=one_body.nparms();

  if(has_two_body)
    tot+=two_body->nparms();

  if(has_three_body)
    tot+=three_body.nparms();
  if(has_three_body_diffspin)
    tot+=three_body_diffspin.nparms();
    
  if(optimize_basis) {
    for(int b=0; b< eibasis.GetDim(0); b++)
      tot+=eibasis(b)->nparms();
    for(int b=0; b< eebasis.GetDim(0); b++)
      tot+=eebasis(b)->nparms();
  }
  return tot;
}

//----------------------------------------------------------------------

int Jastrow_group::valSize() {
  int tot=0;
  return 0; //Temporarily turn off the saving feature..
  //if we're optimizing the basis, we need to recalculate everything
  if(!optimize_basis) {
    if(has_one_body)
      tot+=one_body.nparms();

    if(has_two_body)
      tot+=two_body->nparms();
  }


  return tot;
}

//----------------------------------------------------------------------
void Jastrow_group::getVarParms(Array1 <doublevar> & parms) {

  parms.Resize(nparms());
  Array1 <doublevar> one_parms, two_parms, three_parms, three_parms_diffspin;
  one_body.getParms(one_parms);
  if(has_two_body)
    two_body->getParms(two_parms);

  three_body.getParms(three_parms);
  three_body_diffspin.getParms(three_parms_diffspin);

  int counter=0;
  for(int i=0; i< one_parms.GetDim(0); i++)
    parms(counter++)=one_parms(i);

  for(int i=0; i< two_parms.GetDim(0); i++)
    parms(counter++)=two_parms(i);
  
  for(int i=0; i< three_parms.GetDim(0); i++)
    parms(counter++)=three_parms(i);

  for(int i=0; i< three_parms_diffspin.GetDim(0); i++)
    parms(counter++)=three_parms_diffspin(i);

  if(optimize_basis) {
    Array1 <doublevar> basisparms;
    for(int b=0; b<eibasis.GetDim(0); b++) {
      eibasis(b)->getVarParms(basisparms);
      for(int i=0; i< basisparms.GetDim(0); i++) {
        parms(counter++)=basisparms(i);
      }
    }

    for(int b=0; b<eebasis.GetDim(0); b++) {
      eebasis(b)->getVarParms(basisparms);
      for(int i=0; i< basisparms.GetDim(0); i++) {
        parms(counter++)=basisparms(i);
      }
    }
  }
}

//----------------------------------------------------------------------

void Jastrow_group::setVarParms(Array1 <doublevar> & parms) {
  assert(parms.GetDim(0)==nparms());
  Array1 <doublevar> one_parms, two_parms, three_parms, three_parms_diffspin;
  one_parms.Resize(one_body.nparms());
  if(has_two_body)
    two_parms.Resize(two_body->nparms());
  three_parms.Resize(three_body.nparms());
  three_parms_diffspin.Resize(three_body_diffspin.nparms());

  int counter=0;
  for(int i=0; i< one_parms.GetDim(0); i++)
    one_parms(i)=parms(counter++);

  one_body.setParms(one_parms);

  for(int i=0; i< two_parms.GetDim(0); i++)
    two_parms(i)=parms(counter++);

  if(has_two_body)
    two_body->setParms(two_parms);

  for(int i=0; i< three_parms.GetDim(0); i++)
    three_parms(i)=parms(counter++);
  three_body.setParms(three_parms);  

  for(int i=0; i< three_parms_diffspin.GetDim(0); i++)
    three_parms_diffspin(i)=parms(counter++);
  three_body_diffspin.setParms(three_parms_diffspin);  

  
  if(optimize_basis) {
    Array1 <doublevar> basisparms;
    for(int b=0; b<eibasis.GetDim(0); b++) {
      basisparms.Resize(eibasis(b)->nparms());
      for(int i=0; i< basisparms.GetDim(0); i++) {
        basisparms(i)=parms(counter++);
      }
      eibasis(b)->setVarParms(basisparms);
    }

    for(int b=0; b<eebasis.GetDim(0); b++) {
      basisparms.Resize(eebasis(b)->nparms());

      for(int i=0; i< basisparms.GetDim(0); i++) {
        basisparms(i)=parms(counter++);
      }
      eebasis(b)->setVarParms(basisparms);
    }
  }


}
//----------------------------------------------------------------------

void Jastrow_group::getEEbasisPlot(Array1 <doublevar> &r,
				   Array1 <doublevar> &basisvals) {

  assert(r.GetDim(0) >= 5);
  int nbasis=eebasis.Size();
  int nbasisvals=0;
  for ( int b=0; b<nbasis; b++) nbasisvals+=eebasis(b)->nfunc();
  basisvals.Resize(nbasisvals);
  int startfill=0;
  for ( int b=0; b<nbasis; b++) {
    eebasis(b)->calcVal(r,basisvals,startfill);
    startfill+=eebasis(b)->nfunc();
  }

}

void Jastrow_group::getEIbasisPlotInfo(vector <string> &atnames,
				       Array1 <int> &basissize) {

  int natoms=atomnames.size();
  basissize.Resize(natoms);
  basissize=0;
  for ( int at=0; at<natoms; at++ ) {
    atnames.push_back(atomnames[at]);
    for ( int b=0; b<nbasis_at(at); b++ ) {
      int bb=atom2basis(at,b);
      basissize(at)+=eibasis(bb)->nfunc();
    }
  }

}

void Jastrow_group::getEIbasisPlot(int at, Array1 <doublevar> &r,
				   Array1 <doublevar> &basisvals) {
  assert(r.GetDim(0) >= 5);
  int nbasis=eibasis.Size();
  int nbasisvals=0;
  for ( int b=0; b<nbasis; b++) nbasisvals+=eibasis(b)->nfunc();
  basisvals.Resize(nbasisvals);    // the dimension can be smaller
  int startfill=0;
  for ( int b=0; b<nbasis_at(at); b++) {
    int bb=atom2basis(at,b);
    eibasis(bb)->calcVal(r,basisvals,startfill);
    startfill+=eibasis(bb)->nfunc();
  }
  
}

//----------------------------------------------------------------------



//######################################################################

//----------------------------------------------------------------------
void Jastrow2_wf_data::read(vector <string> & words,
                    unsigned int & pos,
                    System * sys
                  ){
  vector < vector <string> > grouptext;
  vector <string> grouptmp;
  pos=0;
  while(readsection(words, pos, grouptmp, "GROUP"))
    grouptext.push_back(grouptmp);
  group.Resize(grouptext.size());
//  cout << grouptext.size() << " groups found " << endl;
  for(int i=0; i< group.GetDim(0); i++)
    group(i).set_up(grouptext[i], sys);


  nelectrons=sys->nelectrons(0)+sys->nelectrons(1);
  vector <string> atomlabels;
  sys->getAtomicLabels(atomlabels);
  natoms=atomlabels.size();
  nup=sys->nelectrons(0);
  
  
}

//----------------------------------------------------------------------

int Jastrow2_wf_data::supports(wf_support_type support) {
  switch(support) {
  case laplacian_update:
    return 1;
  case density:
    return 1;
  case parameter_derivatives:
    if(nparms()==0) return 1;
    for(int g=0; g< group.GetDim(0); g++) { 
      if(group(g).nparms()>0) { 
        if(group(g).optimizeBasis()) 
          return 0;
      }
    }
    return 1;
  default:
    return 0;
  }
}

//----------------------------------------------------------------------

void Jastrow2_wf_data::getVarParms(Array1 <doublevar> & parms) {
  parms.Resize(nparms());
  Array1 <doublevar> gparms;
  int counter=0;
  for(int g=0; g< group.GetDim(0); g++) {
    group(g).getVarParms(gparms);
    for(int i=0; i< gparms.GetDim(0); i++)
      parms(counter++)=gparms(i);
  }
}

//----------------------------------------------------------------------

void Jastrow2_wf_data::setVarParms(Array1 <doublevar> & parms) {
  assert(parms.GetDim(0)==nparms());
  if(nparms()==0) return;

  Array1 <doublevar> gparms;
  int counter=0;
  for(int g=0; g< group.GetDim(0); g++) {
    gparms.Resize(group(g).nparms());
    for(int i=0; i< gparms.GetDim(0); i++)
      gparms(i)=parms(counter++);
    group(g).setVarParms(gparms);
  }

 int max=wfObserver.size();
 //cout << "max " << max << endl;
  for(int i=0; i< max; i++)
  {
    wfObserver[i]->notify(all_wf_parms_change,0);
  }
}

//----------------------------------------------------------------------

int Jastrow2_wf_data::nparms() {
  int tot=0;
  for(int g=0; g< group.GetDim(0); g++) {
    tot+=group(g).nparms();
  }
  return tot;
}

//----------------------------------------------------------------------

int Jastrow2_wf_data::valSize() {
  if(nparms()==0) { 
    return 1; 
  }
  else { 
    int tot=0;
    for(int g=0; g< group.GetDim(0); g++) {
      tot+=group(g).valSize();
    }
    return tot;
  }
}

//----------------------------------------------------------------------

int Jastrow2_wf_data::writeinput(string & indent, ostream & os) {
  os << indent << "JASTROW2" << endl;
  string indent2=indent+"  ";
  for(int g=0; g< group.GetDim(0); g++) {
    os << indent << "GROUP { " << endl;
    group(g).writeinput(indent2, os);
    os << indent << "}\n";
  }
  //os << "Jastrow2 writeinput" << endl;
  return 1;
}

//----------------------------------------------------------------------

int Jastrow2_wf_data::showinfo( ostream & os) {
  os <<  "Jastrow2 function" << endl;
  string indent="  ";
  for(int g=0; g< group.GetDim(0); g++) {
    os << "Group " << g << endl;
    group(g).showinfo(indent, os);
  }

  return 1;
}

//----------------------------------------------------------------------

void Jastrow2_wf_data::generateWavefunction(Wavefunction * & wf) {
  assert(wf==NULL);
  Jastrow2_wf * jwf=new Jastrow2_wf;
  jwf->init(this);
  attachObserver(jwf);
  wfObserver.push_back(jwf);
  wf=jwf;
}



//--------------------------------------------------------------------------
//##########################################################################
void Jastrow2_wf::notify(change_type change, int num)
{
  
  switch(change)
  {
  case electron_move:
    electronIsStaleVal(num)=1;
    electronIsStaleLap(num)=1;
    //updateEverythingLap=1;
    break;
  case all_electrons_move:
    updateEverythingVal=1;
    updateEverythingLap=1;
    electronIsStaleVal=1;
    electronIsStaleLap=1;
    break;
  case wf_parm_change:
  case all_wf_parms_change:
    parmChanged=1;
    updateEverythingVal=1;
    updateEverythingLap=1;
    electronIsStaleVal=1;
    electronIsStaleLap=1;
    break;
  case sample_attach:
    sampleAttached=1;
    updateEverythingVal=1;
    updateEverythingLap=1;
    electronIsStaleVal=1;
    electronIsStaleLap=1;
    break;
  case data_attach:
    dataAttached=1;
    updateEverythingVal=1;
    updateEverythingLap=1;
    electronIsStaleVal=1;
    electronIsStaleLap=1;
    break;
  default:
    updateEverythingVal=1;
    updateEverythingLap=1;
  }
}

//----------------------------------------------------------

void Jastrow2_wf::init(Wavefunction_data * wfdata) {
  recast(wfdata, parent);
  nelectrons=parent->nelectrons;

  electronIsStaleVal.Resize(nelectrons);
  electronIsStaleLap.Resize(nelectrons);
  electronIsStaleVal=1;
  electronIsStaleLap=1;
  updateEverythingVal=updateEverythingLap=1;
  parmChanged=0;
  two_body_save.Resize(nelectrons, nelectrons, 5);
  two_body_save=0;
  u_twobody=0;
  one_body_save.Resize(nelectrons,5);

  keep_ion_dependent=0;
  int ngroups=parent->group.GetDim(0);
  maxeibasis=0;
  maxeebasis=0;
  for(int g=0; g< ngroups; g++) {
    if(maxeibasis < parent->group(g).maxEIBasis())
      maxeibasis=parent->group(g).maxEIBasis();
    if(maxeebasis < parent->group(g).nEEBasis())
      maxeebasis=parent->group(g).nEEBasis();
  }


  
  eibasis_save.Resize(ngroups);
  for(int g=0; g< ngroups; g++) {
    //if(parent->group(g).hasThreeBody()||parent->group(g).hasThreeBodySpin())
    eibasis_save(g).Resize(nelectrons, parent->natoms, maxeibasis, 5);
  }



}

//----------------------------------------------------------

void Jastrow2_wf::updateVal(Wavefunction_data * wfdata, Sample_point * sample){

  //cout << "updateVal " << endl;
  //cout << "everything " << updateEverythingVal << endl;
  //cout << "stale "; write_array(cout,electronIsStaleVal); cout << endl;
  if(updateEverythingVal) {
    electronIsStaleVal=1;
    updateEverythingVal=0;
  }

  int ngroups=parent->group.GetDim(0);
  Array3 <doublevar> eibasis(parent->natoms, maxeibasis ,5);
  
  bool has3b=false;
  for(int g=0; g< ngroups; g++) 
    has3b=parent->group(g).hasThreeBody() or has3b;
  if(has3b) 
    update_eibasis_save(wfdata,sample);

  for(int e=0; e < nelectrons; e++) {
    if(electronIsStaleVal(e)) {


      Array1 <doublevar> newval_ee(nelectrons);
      doublevar newval_ei;
      newval_ee=0;
      newval_ei=0;

      Array3 <doublevar> eebasis(nelectrons, maxeebasis, 5);

      doublevar old_eval=0;
      for(int i=0; i< e; i++)
        old_eval+=two_body_save(i,e,0);
      for(int j=e+1; j< nelectrons; j++)
        old_eval+=two_body_save(e,j,0);



      if(keep_ion_dependent) { 
        for(int a=0; a< parent->natoms; a++) { 
          for(int d=0; d< 5; d++)
            one_body_ion(e,a,d)=0;
        }
      }


      for(int g=0; g< ngroups; g++) {
        if(parent->group(g).hasOneBody()) 
          parent->group(g).updateEIBasis(e,sample,eibasis);

        if(parent->group(g).hasOneBody()) { 
          parent->group(g).one_body.updateVal(e, eibasis,newval_ei);
          if(keep_ion_dependent) { 
            //cout << "updating " << endl;
            parent->group(g).one_body.updateLap_ion(e, eibasis,one_body_ion);
          }
        }
        

        if(parent->group(g).hasTwoBody() || parent->group(g).hasThreeBody() || parent->group(g).hasThreeBodySpin())
          parent->group(g).updateEEBasis(e,sample, eebasis);
        
        if(parent->group(g).hasTwoBody())
          parent->group(g).two_body->updateVal(e,eebasis, newval_ee);


        if(parent->group(g).hasThreeBody()) 
          parent->group(g).three_body.updateVal(e,eibasis_save(g), 
              eebasis, newval_ee);
        if(parent->group(g).hasThreeBodySpin()) 
          parent->group(g).three_body_diffspin.updateVal(e,eibasis_save(g), 
              eebasis, newval_ee);


      }

      one_body_save(e,0)=newval_ei;

      for(int i=0; i< e; i++)
        two_body_save(i,e,0)=newval_ee(i);
      for(int j=e+1; j< nelectrons; j++)
        two_body_save(e,j,0)=newval_ee(j);

      doublevar new_eval=0;
      for(int i=0; i< e; i++)
        new_eval+=two_body_save(i,e,0);
      for(int j=e+1; j< nelectrons; j++)
        new_eval+=two_body_save(e,j,0);

      u_twobody+=new_eval-old_eval;


      electronIsStaleVal(e)=0;
    }
  }


}
//----------------------------------------------------------

void Jastrow2_wf::update_eibasis_save(Wavefunction_data * wfdata, Sample_point * sample) { 
  int ngroups=parent->group.GetDim(0);
  Array3 <doublevar> eibasis(parent->natoms, maxeibasis ,5);

  for(int g=0; g< ngroups; g++) { 
   // if(parent->group(g).hasThreeBody()|| parent->group(g).hasThreeBodySpin()) { 
      for(int e=0; e< nelectrons; e++) {
        if(electronIsStaleLap(e)) { 
          parent->group(g).updateEIBasis(e,sample,eibasis);
          for(int i=0; i< parent->natoms; i++) {
            for(int j=0; j< maxeibasis; j++) {
              for(int d=0; d< 5; d++) {
                eibasis_save(g)(e,i,j,d)=eibasis(i,j,d);
              }
            }
          }
        }
   //   }
    }
  }
  
}

void Jastrow2_wf::updateLap(Wavefunction_data * wfdata, Sample_point * sample){
  //cout << "Jastrow2_wf::updateLap " << endl;
  if(updateEverythingLap) {
    electronIsStaleLap=1;
    updateEverythingLap=0;
  }

  int ngroups=parent->group.GetDim(0);
  Array3 <doublevar> eibasis(parent->natoms, maxeibasis ,5);
  update_eibasis_save(wfdata,sample);

  for(int e=0; e < nelectrons; e++) {
    if(electronIsStaleLap(e)) {
      //Array3 <doublevar> eibasis(parent->natoms, maxeibasis ,5);
      Array3 <doublevar> eebasis(nelectrons, maxeebasis, 5);

      Array1 <doublevar> newlap_ei(5);
      newlap_ei=0;
      Array2 <doublevar> newlap_ee(nelectrons, 5);
      newlap_ee=0;
      
      Array3 <doublevar> newlap_eei(2,nelectrons,5,0.0);

      newlap_eei=0;

      if(keep_ion_dependent) { 
        for(int a=0; a< parent->natoms; a++) { 
          for(int d=0; d< 5; d++)
            one_body_ion(e,a,d)=0;
        }
      }


      doublevar old_eval=0;
      for(int i=0; i< e; i++)
        old_eval+=two_body_save(i,e,0);
      for(int j=e+1; j< nelectrons; j++)
        old_eval+=two_body_save(e,j,0);

      for(int g=0; g< ngroups; g++) {
        if(parent->group(g).hasOneBody() || parent->group(g).hasThreeBody()
            || parent->group(g).hasThreeBodySpin() )
          parent->group(g).updateEIBasis(e,sample,eibasis);

        if(parent->group(g).hasOneBody()) {
          parent->group(g).one_body.updateLap(e, eibasis,newlap_ei);
          if(keep_ion_dependent) { 
            parent->group(g).one_body.updateLap_ion(e, eibasis,one_body_ion);
          }
        }


        if(parent->group(g).hasTwoBody() || parent->group(g).hasThreeBody() 
            || parent->group(g).hasThreeBodySpin())
          parent->group(g).updateEEBasis(e,sample, eebasis);
        
        if(parent->group(g).hasTwoBody())
          parent->group(g).two_body->updateLap(e,eebasis, newlap_ee);

          
        if(parent->group(g).hasThreeBody()) { 
          for(int i=0; i< parent->natoms; i++) {
            for(int j=0; j< maxeibasis; j++) {
              for(int d=0; d< 5; d++) {
                eibasis_save(g)(e,i,j,d)=eibasis(i,j,d);
              }
            }
          }
          parent->group(g).three_body.updateLap(e,eibasis_save(g),eebasis,newlap_eei);
        }

        if(parent->group(g).hasThreeBodySpin()) { 
          for(int i=0; i< parent->natoms; i++) {
            for(int j=0; j< maxeibasis; j++) {
              for(int d=0; d< 5; d++) {
                eibasis_save(g)(e,i,j,d)=eibasis(i,j,d);
              }
            }
          }
          parent->group(g).three_body_diffspin.updateLap(e,eibasis_save(g),eebasis,newlap_eei);
        }

      }


      for(int d=0; d< 5; d++)
        one_body_save(e,d)=newlap_ei(d);
      

      for(int i=0; i< e; i++) {
        for(int d=0; d< 5; d++) 
          two_body_save(i,e,d)=newlap_ee(i,d);

        two_body_save(i,e,0)+=newlap_eei(0,i,0);
        for(int d=1; d < 5; d++) 
          two_body_save(i,e,d)+=newlap_eei(1,i,d);
        
      }
      for(int i=0; i< e; i++) {
       for(int d=1; d< 4; d++)
          two_body_save(e,i,d)= -newlap_ee(i,d);
        two_body_save(e,i,4)=newlap_ee(i,4);
        
        for(int d=1; d< 5; d++) 
          two_body_save(e,i,d)+=newlap_eei(0,i,d);
      }

      for(int j=e+1; j< nelectrons; j++) {
        for(int d=0; d< 5; d++) 
          two_body_save(e,j,d)=newlap_ee(j,d);
          
        two_body_save(e,j,0)+=newlap_eei(0,j,0);
        
        for(int d=1; d< 5; d++) 
          two_body_save(e,j,d)+=newlap_eei(0,j,d);
        
        
      }
      for(int j=e+1; j< nelectrons; j++) {
        for(int d=1; d< 4; d++)
          two_body_save(j,e,d)= -newlap_ee(j,d);
        two_body_save(j,e,4)=newlap_ee(j,4);
        
        for(int d=1; d< 5; d++) 
          two_body_save(j,e,d)+=newlap_eei(1,j,d);

      }
      

      doublevar new_eval=0;
      for(int i=0; i< e; i++)
        new_eval+=two_body_save(i,e,0);
      for(int j=e+1; j< nelectrons; j++)
        new_eval+=two_body_save(e,j,0);

      u_twobody+=new_eval-old_eval;

      
      /*
      for(int d=0; d< 1; d++) { 
        cout << "updated two_body_save: " << d << endl;
        for(int i=0; i < nelectrons; i++) {
          for(int j=0; j< nelectrons; j++) {
            cout << two_body_save(i,j,d) << "  ";
          }
          cout << endl;
        }
      }
      
      

      cout << "u " << u_twobody << "  new_eval " << new_eval << "  old_eval " << old_eval << endl;
      */
      
      electronIsStaleLap(e)=0;
    }
  }

  //for(int i=0; i< nelectrons; i++) {
  //  for(int j=i+1; j< nelectrons; j++) {
  //    cout << "pair " << i << "   " << j
  //         << " driftion " << two_body_save(i,j,4) << endl;
  //  }
  //}

  electronIsStaleVal=0;
  updateEverythingVal=0;
  //cout << "Jastrow2_wf::updateLap done" << endl;

}
//----------------------------------------------------------
void Jastrow2_wf::updateForceBias(Wavefunction_data * wfdata,
                                  Sample_point * sample){
  updateVal(wfdata, sample);
}
//----------------------------------------------------------

void Jastrow2_wf::getVal(Wavefunction_data * wfdata,
                         int e, Wf_return & val) {

  //assert(val.amp.GetDim(1) >= 2);
  val.phase(0, 0)=0;
  doublevar u=0;
  for(int i=0; i< nelectrons; i++) {
    u+=one_body_save(i,0);
  }

  u+=u_twobody;

  val.amp(0, 0)=u;
  val.cvals(0,0)=u;


}


//----------------------------------------------------------
void Jastrow2_wf::getSymmetricVal(Wavefunction_data * wfdata,
		     int e, Wf_return & val){
  getVal(wfdata,e,val);
} 


//----------------------------------------------------------

void Jastrow2_wf::getLap(Wavefunction_data * wfdata, int e, 
                         Wf_return & lap) {

  lap.amp=0;
  lap.phase=0;

  getVal(wfdata, e, lap);
  doublevar dotproduct=0;


  for(int d=1; d< 4; d++) {
    lap.amp(0, d)+=one_body_save(e,d);
    for(int i=0; i< nelectrons; i++) {
      lap.amp(0,d)+=two_body_save(e,i,d);
      
    }
    dotproduct+=lap.amp(0,d)*lap.amp(0,d);
  }

  for(int i=0; i< nelectrons; i++) {
    lap.amp(0,4)+=two_body_save(e,i,4);
  }
  lap.amp(0,4)+=one_body_save(e,4)+dotproduct;

  for(int i=1; i< 5; i++) 
    lap.cvals(0,i)=lap.amp(0,i);


}

//----------------------------------------------------------

void Jastrow2_wf::getForceBias(Wavefunction_data *wfdata, int e,
                               Wf_return & force){
  assert(force.amp.GetDim(1) >= 5);
  force.amp=0;
  force.phase=0;
  getVal(wfdata, e, force);
}

//----------------------------------------------------------


void Jastrow2_wf::generateStorage(Wavefunction_storage *& wfstore) {
  Jastrow2_storage * store=new Jastrow2_storage(nelectrons);
  wfstore=store;
  int ngroups=eibasis_save.GetDim(0);
  store->eibasis.Resize(ngroups);
  store->eibasis_2.Resize(ngroups);


  for(int g=0; g< ngroups; g++) {
    store->eibasis(g).Resize(eibasis_save(g).GetDim(1), eibasis_save(g).GetDim(2), 
                             eibasis_save(g).GetDim(3));
    store->eibasis_2(g).Resize(eibasis_save(g).GetDim(1), eibasis_save(g).GetDim(2), 
                             eibasis_save(g).GetDim(3));
  }
  
  
}

void Jastrow2_wf::saveUpdate(Sample_point * sample, int e,
                             Wavefunction_storage * wfstore){
  Jastrow2_storage * j2store;
  recast(wfstore, j2store);
 
  for(int d=0; d< 5; d++) {
    j2store->one_body_part(d)=one_body_save(e,d);
  }

  for(int i=0; i< nelectrons; i++) {
    for(int d=0; d< 5; d++) {
      j2store->two_body_part_e(i,d)=two_body_save(e,i,d);
    }
  }
  for(int j=0; j< nelectrons; j++) {
    for(int d=0; d< 5; d++) {
      j2store->two_body_part_others(j,d)=two_body_save(j,e,d);
    }
  }
  
  for(int g=0; g< eibasis_save.GetDim(0); g++) {
    for(int at=0; at < eibasis_save(g).GetDim(1); at++) {
      for(int i=0; i < eibasis_save(g).GetDim(2); i++) {
        for(int d=0; d < eibasis_save(g).GetDim(3); d++) {
          j2store->eibasis(g)(at,i,d)=eibasis_save(g)(e,at,i,d);
        }
      }
    }
  }
  
}
//----------------------------------------------------------
void Jastrow2_wf::restoreUpdate(Sample_point * sample, int e, Wavefunction_storage * wfstore){
  Jastrow2_storage * j2store;
  recast(wfstore, j2store);


  for(int d=0; d< 5; d++) {
    one_body_save(e,d)=j2store->one_body_part(d);
  }

  doublevar old_eval=0;
  for(int i=0; i< e; i++)
    old_eval+=two_body_save(i,e,0);
  for(int j=e+1; j< nelectrons; j++)
    old_eval+=two_body_save(e,j,0);


  for(int i=0; i< nelectrons; i++) {
    for(int d=0; d< 5; d++) {
      two_body_save(e,i,d)=j2store->two_body_part_e(i,d);
    }
  }
  for(int j=0; j< nelectrons; j++) {
    for(int d=0; d< 5; d++) {
      two_body_save(j,e,d)=j2store->two_body_part_others(j,d);
    }
  }
  
  for(int g=0; g< eibasis_save.GetDim(0); g++) {
    for(int at=0; at < eibasis_save(g).GetDim(1); at++) {
      for(int i=0; i < eibasis_save(g).GetDim(2); i++) {
        for(int d=0; d < eibasis_save(g).GetDim(3); d++) {
          eibasis_save(g)(e,at,i,d)=j2store->eibasis(g)(at,i,d);
        }
      }
    }
  }


  doublevar new_eval=0;
  for(int i=0; i< e; i++)
    new_eval+=two_body_save(i,e,0);
  for(int j=e+1; j< nelectrons; j++)
    new_eval+=two_body_save(e,j,0);

   u_twobody+=new_eval-old_eval;

}

//----------------------------------------------------------
void Jastrow2_wf::saveUpdate(Sample_point * sample, int e1, int e2,
                             Wavefunction_storage * wfstore){
  Jastrow2_storage * j2store;
  recast(wfstore, j2store);

  for(int d=0; d< 5; d++) {
    j2store->one_body_part(d)=one_body_save(e1,d);
  }

  for(int i=0; i< nelectrons; i++) {
    for(int d=0; d< 5; d++) {
      j2store->two_body_part_e(i,d)=two_body_save(e1,i,d);
    }
  }
  for(int j=0; j< nelectrons; j++) {
    for(int d=0; d< 5; d++) {
      j2store->two_body_part_others(j,d)=two_body_save(j,e1,d);
    }
  }
  
  for(int g=0; g< eibasis_save.GetDim(0); g++) {
    for(int at=0; at < eibasis_save(g).GetDim(1); at++) {
      for(int i=0; i < eibasis_save(g).GetDim(2); i++) {
        for(int d=0; d < eibasis_save(g).GetDim(3); d++) {
          j2store->eibasis(g)(at,i,d)=eibasis_save(g)(e1,at,i,d);
        }
      }
    }
  }
  
  for(int d=0; d< 5; d++) {
    j2store->one_body_part_2(d)=one_body_save(e2,d);
  }
  
  for(int i=0; i< nelectrons; i++) {
    for(int d=0; d< 5; d++) {
      j2store->two_body_part_e_2(i,d)=two_body_save(e2,i,d);
    }
  }
  
  for(int j=0; j< nelectrons; j++) {
    for(int d=0; d< 5; d++) {
      j2store->two_body_part_others_2(j,d)=two_body_save(j,e2,d);
    }
  }
  
  for(int g=0; g< eibasis_save.GetDim(0); g++) {
    for(int at=0; at < eibasis_save(g).GetDim(1); at++) {
      for(int i=0; i < eibasis_save(g).GetDim(2); i++) {
        for(int d=0; d < eibasis_save(g).GetDim(3); d++) {
          j2store->eibasis_2(g)(at,i,d)=eibasis_save(g)(e2,at,i,d);
        }
      }
    }
  }
}

void Jastrow2_wf::restoreUpdate(Sample_point * sample, int e1, int e2, Wavefunction_storage * wfstore){
  Jastrow2_storage * j2store;
  recast(wfstore, j2store);
  
  //  assert( spin(e1) != spin(e2) );

  doublevar old_eval=0;
  for(int i=0; i< e2; i++)
    old_eval+=two_body_save(i,e2,0);
  for(int j=e2+1; j< nelectrons; j++)
    old_eval+=two_body_save(e2,j,0);
  
  for(int g=0; g< eibasis_save.GetDim(0); g++) {
    for(int at=0; at < eibasis_save(g).GetDim(1); at++) {
      for(int i=0; i < eibasis_save(g).GetDim(2); i++) {
        for(int d=0; d < eibasis_save(g).GetDim(3); d++) {
          eibasis_save(g)(e2,at,i,d)=j2store->eibasis_2(g)(at,i,d);
        }
      }
    }
  }

  for(int d=0; d< 5; d++) 
    one_body_save(e2,d)=j2store->one_body_part_2(d);

  for(int i=0; i< nelectrons; i++) {
    for(int d=0; d< 5; d++) {
      two_body_save(e2,i,d)=j2store->two_body_part_e_2(i,d);
    }
  }
  for(int j=0; j< nelectrons; j++) {
    for(int d=0; d< 5; d++) {
      two_body_save(j,e2,d)=j2store->two_body_part_others_2(j,d);
    }
  }
  


  doublevar new_eval=0;
  for(int i=0; i< e2; i++)
    new_eval+=two_body_save(i,e2,0);
  for(int j=e2+1; j< nelectrons; j++)
    new_eval+=two_body_save(e2,j,0);

  u_twobody += new_eval-old_eval;


  for(int g=0; g< eibasis_save.GetDim(0); g++) {
    for(int at=0; at < eibasis_save(g).GetDim(1); at++) {
      for(int i=0; i < eibasis_save(g).GetDim(2); i++) {
        for(int d=0; d < eibasis_save(g).GetDim(3); d++) {
          eibasis_save(g)(e1,at,i,d)=j2store->eibasis(g)(at,i,d);
        }
      }
    }
  }

  old_eval=0;
  for(int i=0; i< e1; i++)
    old_eval+=two_body_save(i,e1,0);
  for(int j=e1+1; j< nelectrons; j++)
    old_eval+=two_body_save(e1,j,0);

  
  for(int d=0; d< 5; d++) {
    one_body_save(e1,d)=j2store->one_body_part(d);
  }

  for(int i=0; i< nelectrons; i++) {
    for(int d=0; d< 5; d++) {
      two_body_save(e1,i,d)=j2store->two_body_part_e(i,d);
    }
  }
  for(int j=0; j< nelectrons; j++) {
    for(int d=0; d< 5; d++) {
      two_body_save(j,e1,d)=j2store->two_body_part_others(j,d);
    }
  }
  


  new_eval=0;
  for(int i=0; i< e1; i++)
    new_eval+=two_body_save(i,e1,0);
  for(int j=e1+1; j< nelectrons; j++)
    new_eval+=two_body_save(e1,j,0);

  u_twobody+=new_eval-old_eval;
	
}
//----------------------------------------------------------


//Note that one could implement this as updates to improve speed further.
int Jastrow2_wf::getParmDeriv(Wavefunction_data *wfdata , Sample_point * sample,
                               Parm_deriv_return & parm_deriv) {
  sample->updateEEDist();
  sample->updateEIDist();
  //we're not going to support nonlinear parameters for the first iteration of this.
  //hopefully it won't be overly inefficient to do that.
  int ng=parent->group.GetDim(0);
  for(int g=0; g< ng; g++) { 
    if(parent->group(g).optimizeBasis()) return 0;
  }
  update_eibasis_save(wfdata,sample);

  parm_deriv.gradient.Resize(0);
  parm_deriv.hessian.Resize(0,0);
  parm_deriv.gradderiv.Resize(0,nelectrons,4);
  parm_deriv.val_gradient.Resize(nelectrons,3);
  parm_deriv.val_gradient=0.0;
  
  for(int g=0; g< ng; g++) {
    
    if(parent->group(g).hasOneBody() and parent->group(g).one_body.nparms() > 0 ) {
      Parm_deriv_return tmp_parm;
      
      parent->group(g).one_body.getParmDeriv(eibasis_save(g),tmp_parm);
      extend_parm_deriv(parm_deriv,tmp_parm);
    }

    Array4 <doublevar> eetotal(nelectrons, nelectrons, maxeebasis,5);
    eetotal=-1;
    Array3 <doublevar> eebasis(nelectrons, maxeebasis, 5);
    for(int e=0; e< nelectrons; e++) { 
      parent->group(g).updateEEBasis(e,sample, eebasis);
      for(int j=0; j< e; j++) { 
        for(int b=0; b< maxeebasis; b++) {
          for(int d=0; d< 5; d++) 
            eetotal(j,e,b,d)=eebasis(j,b,d);
        }
      }
    }
    if(parent->group(g).hasTwoBody() and parent->group(g).two_body->nparms() > 0 ) {
      Parm_deriv_return tmp_parm;
      parent->group(g).two_body->getParmDeriv(eetotal,tmp_parm);
      extend_parm_deriv(parm_deriv,tmp_parm);
    }

    if(parent->group(g).hasThreeBody() and parent->group(g).three_body.nparms() > 0 ) {
      Parm_deriv_return tmp_parm;
      parent->group(g).three_body.getParmDeriv(eibasis_save(g),eetotal,tmp_parm);
      extend_parm_deriv(parm_deriv,tmp_parm);
    }
    
    //retparm=tmp_parm;
    

    
  }
  
  int np=parent->nparms();
  for(int i=0; i< np; i++) { 
    for(int j=0; j< np; j++) {
      parm_deriv.hessian(i,j)=parm_deriv.gradient(i)*parm_deriv.gradient(j);
    }
  }
  return 1;
}

//--------------------------------------------------------------------------

// JK: to implement analytical derivatives in BCS_wf (which uses two-body
// Jastrow piece as a pair orbital), "electron-resolved" ParamDeriv is needed
int Jastrow2_wf::get_twobody_ParmDeriv(Sample_point * sample,
				       Array3<doublevar>& twobody_parm_deriv ) {
  
  sample->updateEEDist();

  int parm_count=-1;
  int nparm=0;
  
  int ng=parent->group.GetDim(0);
  for(int g=0; g< ng; g++) {
    if(parent->group(g).hasTwoBody()&& parent->group(g).two_body->nparms()) {
      nparm+=parent->group(g).two_body->nparms();
    }
  }

  twobody_parm_deriv.Resize(nelectrons, nelectrons,nparm);

  for(int g=0; g< ng; g++) {
    if(parent->group(g).hasTwoBody()&& parent->group(g).two_body->nparms()) {

      Array3 <doublevar> eetotal(nelectrons, nelectrons, maxeebasis);
      eetotal=-1;
      Array3 <doublevar> eebasis(nelectrons, maxeebasis, 5);
      for(int e=0; e< nelectrons; e++) { 
        parent->group(g).updateEEBasis(e,sample,eebasis);
        for(int j=0; j< e; j++) { 
          for(int b=0; b< maxeebasis; b++) {
            eetotal(j,e,b)=eebasis(j,b,0);
          }
        }
      }
      

      // JK: how do I know that the only parameters are those for linear
      // combination of basis functions?
      // how do I know the ordering is consistent with the rest of the code?


      int np=parent->group(g).two_body->nparms();
      for (int p=0; p<np; p++) {
        parm_count++;
        for (int i=0; i<nelectrons; i++) {
          for (int j=0; j<nelectrons; j++) {
            twobody_parm_deriv(i,j,parm_count)=eetotal(i,j,p);
          }
        }
      }

    }
  }
     
  return 1;

}

//--------------------------------------------------------------------------

void Jastrow2_wf::plot1DInternals(Array1 <doublevar> & xdata,
				  vector <Array1 <doublevar> > & data,
				  vector <string> & desc,
				  string desc0 ) {
  
  Array1 <doublevar> column, column2;
  column.Resize(xdata.Size());

  Array1 <doublevar> r;
  r.Resize(5);

  int ng=parent->group.GetDim(0);
  for(int g=0; g< ng; g++) {
    
    if ( parent->group(g).hasTwoBody() ) {
      
      ostringstream dsc;
      dsc << desc0 << "group " << g << ", two-body";

      Array1 <doublevar> parms;
      // getParms will resize the parms array for us, 
      // BUT: if the parameters are frozen, we don't get them !!
      parent->group(g).two_body->getParms(parms);
      if ( parms.Size()==0 ) {
        parent->group(g).two_body->unfreeze();
        parent->group(g).two_body->getParms(parms);
        // return wave function to the original state (we might want to
        // continue in optimization after plotting)
        parent->group(g).two_body->refreeze();
      }

      int nbasisvals=parent->group(g).two_body->nbasis_needed();
      Array1 <doublevar> basisvals;
      //basisvals.Resize(nbasisvals);
      if ( ( parms.Size()!=nbasisvals ) && ( parms.Size()!=2*nbasisvals ) ) {
        error("Something is wrong in Jastrow2_wf::plot1DInternals, inconsistent number of parameters in the two-body part.");
      }

      int diff_spin=0;
      if ( parms.Size()==2*nbasisvals) {
        diff_spin=1;
        column2.Resize(xdata.Size());
        column2=0.0;
      }

      for ( int i=0; i<xdata.Size(); i++) {
        r=0.0;
        r(2)=xdata(i);    // x-direction
        r(0)=r(2);        // norm
        r(1)=r(0)*r(0);   // norm squared
        parent->group(g).getEEbasisPlot(r,basisvals);
        column(i)=0.0;
        for ( int j=0; j<nbasisvals; j++ ) {
          column(i)+=parms(j)*basisvals(j);
          if ( diff_spin ) column2(i)+=parms(j+nbasisvals)*basisvals(j);
        }
      }

      data.push_back(column);
      if (!diff_spin) {
        desc.push_back(dsc.str());
      } else {
        data.push_back(column2);
        dsc << ", like spins";
        desc.push_back(dsc.str());
        dsc.str("");
        dsc << desc0 << "group " << g << ", two-body, unlike spins";
        desc.push_back(dsc.str());
      }

    } // if ( parent->group(g).hasTwoBody() ) {

    // JK: The one-body implementation is not very elegant. There is
    // probably a better way how to match coefficients (params) with
    // basis functions on atoms.
    
    if ( parent->group(g).hasOneBody() ) {


      // get number of atoms, names and sizes of basis-sets there centered
      vector <string> atomnames;
      Array1 <int> basissize;
      parent->group(g).getEIbasisPlotInfo(atomnames,basissize);
      int natoms=atomnames.size();     

      // get number of different kinds of atoms (or different one-body
      // parts, eventually)
      int natomkinds=0;
      for ( int at=0; at<natoms; at++ ) {
        if ( parent->group(g).one_body.atom_kind(at) > natomkinds ) {
          natomkinds=parent->group(g).one_body.atom_kind(at);
        }
      }
      natomkinds++;

      // find a representative atom to each kind
      Array1 <int> unique_atoms;
      unique_atoms.Resize(natomkinds);
      for ( int at=0; at<natoms; at++ ) {
        unique_atoms(parent->group(g).one_body.atom_kind(at))=at;
      }

      // where is the start (and the end) of coefficients for given
      // kind in the params array
      Array1 <int> basisstart;
      basisstart.Resize(natomkinds+1);
      basisstart(0)=0;
      for ( int kind=0; kind < natomkinds; kind++ ) {
        basisstart(kind+1)=basisstart(kind)+basissize(unique_atoms(kind));
      }

      // getParms will resize the parms array for us, 
      // BUT: if the parameters are frozen, we don't get them !!
      Array1 <doublevar> parms;
      parent->group(g).one_body.getParms(parms);
      if ( parms.Size()== 0 ) {
        parent->group(g).one_body.unfreeze();
        parent->group(g).one_body.getParms(parms);
        // return wave function to the original state (we might want to
        // continue in optimization after plotting)
        parent->group(g).one_body.refreeze();
      }

      Array1 <doublevar> basisvals;
      for ( int kind=0; kind<natomkinds; kind++ ) {

        int at=unique_atoms(kind);

        ostringstream dsc;
        dsc << desc0 << "group " << g << ", one-body, center "
          << atomnames[at];

        for ( int i=0; i<xdata.Size(); i++) {
          r=0.0;
          r(2)=xdata(i);    // x-direction
          r(0)=r(2);        // norm
          r(1)=r(0)*r(0);   // norm squared
          parent->group(g).getEIbasisPlot(at,r,basisvals);
          column(i)=0.0;
          for ( int j=basisstart(kind); j<basisstart(kind+1); j++ ) {
            column(i)+=parms(j)*basisvals(j-basisstart(kind));
          }
        }

        data.push_back(column);
        desc.push_back(dsc.str());

      }

    } // // if ( parent->group(g).hasOneBody() ) {

  } // for(int g=0; g< ng; g++) {

}

//--------------------------------------------------------------------------
void Jastrow2_wf::evalTestPos(Array1 <doublevar> & pos, Sample_point * sample,
    Array1 <Wf_return> & wf) { 

  wf.Resize(nelectrons);
  /*
  Wavefunction_storage * store=NULL;
  generateStorage(store);
  Array1 <doublevar> oldpos(3);

  for(int e=0; e< nelectrons; e++) { 
    wf(e).Resize(1,2);
    saveUpdate(sample,e,store);
    sample->getElectronPos(e,oldpos);
    sample->setElectronPosNoNotify(e,pos);
    electronIsStaleLap(e)=1; electronIsStaleVal(e)=1;
    updateVal(parent,sample);
    getVal(parent,e,wf(e));
    restoreUpdate(sample,e,store);
    sample->setElectronPosNoNotify(e,oldpos);
  }

  delete  store;
  */

  Array1 <doublevar> oldpos(3);
  sample->getElectronPos(0,oldpos);
  sample->setElectronPosNoNotify(0,pos);
  int ngroups=parent->group.GetDim(0);
  Array1 <Array3 <doublevar> > eibasis(ngroups);
  Array1 <Array3 <doublevar> > eebasis(ngroups);
  for(int g=0; g< ngroups; g++) { 
    eibasis(g).Resize(parent->natoms, maxeibasis ,5);
    eebasis(g).Resize(nelectrons,maxeebasis,5);
    parent->group(g).updateEIBasis(0,sample,eibasis(g));
    parent->group(g).updateEEBasis(0,sample,eebasis(g));
  }
  sample->setElectronPosNoNotify(0,oldpos);

  //We have to also get the eebasis for the test with electron 0.
  //Doing this in a somewhat inefficient way
  sample->getElectronPos(1,oldpos);
  sample->setElectronPosNoNotify(1,pos);
  for(int g=0; g< ngroups;g++) { 
    //Using the fact that updateEEBasis(e,..) doesn't touch element e
    parent->group(g).updateEEBasis(1,sample,eebasis(g));
  }
  sample->setElectronPosNoNotify(1,oldpos);
    

  //We've now calculated the basis functions for all the electrons and can 
  //get the wave function values
  //Again doing two electrons per spin channel, slightly inefficient, but the scaling is 
  //still fine.
  //One could save some time by separating out the spin-dependent and spin-independent
  //parts at the cost of code complexity, but so far it doesn't seem to be a 
  //limiting factor for the RDM's

  Array1 <Array1 <doublevar> > newval_ee(2); //one for each spin channel
  Array3 <doublevar> eibasis_tmp(parent->natoms,maxeibasis,5);
  doublevar u_one=0; 
  for(int i=0; i< nelectrons; i++) u_one+=one_body_save(i,0);
  doublevar newval_ei;
  int nup=parent->nup;
  for(int s=0; s< 2; s++) { 
    newval_ee(s).Resize(nelectrons);
    newval_ee(s)=0.0;
    int ne_spin;
    if(s==0) ne_spin=min(2,nup);
    else ne_spin=min(nup+2,nelectrons);
    for(int e=s*nup; e< ne_spin;e++) { 
      newval_ei=0;
      for(int e1=0; e1 < nelectrons; e1++) if(e1!=e) newval_ee(s)(e1)=0;

      for(int g=0; g< ngroups;g++) { 
        if(parent->group(g).hasOneBody()) { 
          parent->group(g).one_body.updateVal(e,eibasis(g),newval_ei);
        }
        if(parent->group(g).hasTwoBody()) 
          parent->group(g).two_body->updateVal(e,eebasis(g),newval_ee(s));

        //Here we have to do some shifting around of values
        if(parent->group(g).hasThreeBody() || parent->group(g).hasThreeBodySpin()) {  
          for(int i=0; i< parent->natoms; i++) { 
            for(int j=0; j< maxeibasis; j++) { 
              for(int d=0; d< 5; d++) { 
                eibasis_tmp(i,j,d)=eibasis_save(g)(e,i,j,d);
                eibasis_save(g)(e,i,j,d)=eibasis(g)(i,j,d);
              }
            }
          }
          if(parent->group(g).hasThreeBody()) 
            parent->group(g).three_body.updateVal(e,eibasis_save(g), 
                eebasis(g), newval_ee(s));
          if(parent->group(g).hasThreeBodySpin()) 
            parent->group(g).three_body_diffspin.updateVal(e,eibasis_save(g), 
                eebasis(g), newval_ee(s));
          for(int i=0; i< parent->natoms; i++) { 
            for(int j=0; j< maxeibasis; j++) { 
              for(int d=0; d< 5; d++) { 
                eibasis_save(g)(e,i,j,d)=eibasis_tmp(i,j,d);
              }
            }
          }

        }
        //----

      }
    }
  }
  //Here we calculate the wave function for each of the electrons going to the 
  //test position
  for(int e=0; e< nelectrons; e++) { 
    int s=0;
    if( e >= nup ) s=1;
    wf(e).Resize(1,1);
    
    doublevar old_eval=0,new_eval=0;
    for(int i=0; i< e; i++) old_eval+=two_body_save(i,e,0);
    for(int j=e+1; j< nelectrons; j++) old_eval+=two_body_save(e,j,0);
    for(int i=0; i< e; i++) new_eval+=newval_ee(s)(i);
    for(int j=e+1; j< nelectrons; j++) new_eval+=newval_ee(s)(j);
    doublevar u=u_twobody+u_one //original
      +new_eval-old_eval+newval_ei-one_body_save(e,0);//updates
    wf(e).amp(0,0)=u;
    wf(e).phase(0,0)=0;
    wf(e).cvals(0,0)=u;
  }
  
}
//--------------------------------------------------------------------------
void create_parm_deriv(const Array3 <doublevar> & func,
    const Array1 <doublevar> & coeff,
    Parm_deriv_return & parm_deriv) { 
  int np=func.GetDim(0);
  int nelectrons=func.GetDim(1);
  assert(coeff.GetDim(0)==np);
  parm_deriv.gradient.Resize(np);
  parm_deriv.gradient=0;
  parm_deriv.hessian.Resize(np,np);
  parm_deriv.hessian=0;
  
  //Create the derivatives and hessian
  for(int e=0; e < nelectrons; e++) {
    for(int p=0; p < np; p++) { 
      parm_deriv.gradient(p)+=func(p,e,0);

    }
  }

  for(int i=0; i< np; i++) 
    parm_deriv.hessian(i,i)=parm_deriv.gradient(i)*parm_deriv.gradient(i);
  
  //Form the derivatives of the electronic gradients
  parm_deriv.gradderiv.Resize(np,nelectrons,4);
  for(int p=0; p < np; p++) { 
    for(int e=0; e< nelectrons; e++) { 
      for(int d=0; d< 4; d++) { 
        parm_deriv.gradderiv(p,e,d)=func(p,e,d+1);
      }
    }
  }
  //sum up the parameter derivative.  This is particularly
  //simple for the linear parameters.

  parm_deriv.val_gradient.Resize(nelectrons,3);
  parm_deriv.val_gradient=0.0;
  for(int p=0; p < np; p++) { 
    for(int e=0; e< nelectrons; e++) { 
      for(int d=0; d< 3; d++) { 
        parm_deriv.val_gradient(e,d)+=coeff(p)*parm_deriv.gradderiv(p,e,d);
      }
    }
  }
  //Now calculate the cross terms and add them to the laplacian
  //derivatives
  for(int p=0; p < np; p++) { 
    for(int e=0; e< nelectrons; e++) { 
      doublevar dot=0;
      for(int d=0;d < 3; d++) dot+=parm_deriv.val_gradient(e,d)*parm_deriv.gradderiv(p,e,d);
      parm_deriv.gradderiv(p,e,3)+=2*dot;
    }
  }

}

//----------------------------------------------------------------------
void create_parm_deriv_frozen(const Array3 <doublevar> & func,
    const Array1 <doublevar> & coeff,
    Parm_deriv_return & parm_deriv) { 
  int np=func.GetDim(0);
  int nelectrons=func.GetDim(1);
  assert(coeff.GetDim(0)==np);
  parm_deriv.gradient.Resize(0);
  parm_deriv.gradient=0;
  parm_deriv.hessian.Resize(0,0);
  parm_deriv.hessian=0;
  parm_deriv.gradderiv.Resize(0,nelectrons,5);
  
  
  //sum up the parameter derivative.  This is particularly
  //simple for the linear parameters.

  parm_deriv.val_gradient.Resize(nelectrons,3);
  parm_deriv.val_gradient=0.0;
  for(int p=0; p < np; p++) { 
    for(int e=0; e< nelectrons; e++) { 
      for(int d=0; d< 3; d++) { 
        parm_deriv.val_gradient(e,d)+=coeff(p)*func(p,e,d+1);
      }
    }
  }
}

