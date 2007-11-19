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



//######################################################################

//----------------------------------------------------------------------
void Jastrow2_wf_data::read(vector <string> & words,
                    unsigned int & pos,
                    System * sys
                  ){

  //cout << "Jastrow2 initialization " << endl;

  vector < vector <string> > grouptext;
  vector <string> grouptmp;
  pos=0;
  while(readsection(words, pos, grouptmp, "GROUP"))
    grouptext.push_back(grouptmp);
  group.Resize(grouptext.size());
  //cout << grouptext.size() << " groups found " << endl;
  for(int i=0; i< group.GetDim(0); i++)
    group(i).set_up(grouptext[i], sys);


  nelectrons=sys->nelectrons(0)+sys->nelectrons(1);
  vector <string> atomlabels;
  sys->getAtomicLabels(atomlabels);
  natoms=atomlabels.size();
  //cout << "jastrow2 done" << endl;
  
  
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
    int ng=group.GetDim(0);
    for(int g=0; g< ng; g++) { 
      if(group(g).nparms()>0) { 
        if(group(g).hasThreeBodySpin() || group(g).optimizeBasis()) 
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
  wf=new Jastrow2_wf;
  wf->init(this);
  attachObserver(wf);
  wfObserver.push_back(wf);
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
  case sample_static:
    staticSample=1;
    break;
  case sample_dynamic:
    staticSample=0;
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
  staticSample=0;
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
    if(parent->group(g).hasThreeBody()||parent->group(g).hasThreeBodySpin())
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

  Array3 <doublevar> eibasis(parent->natoms, maxeibasis ,5);
  int ngroups=parent->group.GetDim(0);
  for(int g=0; g< ngroups; g++) { 
    if(parent->group(g).hasThreeBody()||parent->group(g).hasThreeBodySpin()) { 
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
      }
    }
  }


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
void Jastrow2_wf::updateLap(Wavefunction_data * wfdata, Sample_point * sample){
  //cout << "Jastrow2_wf::updateLap " << endl;
  if(updateEverythingLap) {
    electronIsStaleLap=1;
    updateEverythingLap=0;
  }

  int ngroups=parent->group.GetDim(0);
  Array3 <doublevar> eibasis(parent->natoms, maxeibasis ,5);

  for(int g=0; g< ngroups; g++) { 
    if(parent->group(g).hasThreeBody()|| parent->group(g).hasThreeBodySpin()) { 
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
      }
    }
  }

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
      
        if(parent->group(g).hasOneBody() || parent->group(g).hasThreeBody()|| parent->group(g).hasThreeBodySpin() )
	  parent->group(g).updateEIBasis(e,sample,eibasis);
        
        if(parent->group(g).hasOneBody()) {
          parent->group(g).one_body.updateLap(e, eibasis,newlap_ei);
	  if(keep_ion_dependent) { 
	    parent->group(g).one_body.updateLap_ion(e, eibasis,one_body_ion);
	  }
	}
        

        if(parent->group(g).hasTwoBody() || parent->group(g).hasThreeBody() || parent->group(g).hasThreeBodySpin())
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

void Jastrow2_wf::getDensity(Wavefunction_data *,int,  Array2 <doublevar> &){
  error("Jastrow2 doesn't have density yet");
}

//----------------------------------------------------------

void Jastrow2_wf::generateStorage(Wavefunction_storage *& wfstore) {
  Jastrow2_storage * store=new Jastrow2_storage(nelectrons);
  wfstore=store;
  int ngroups=eibasis_save.GetDim(0);
  store->eibasis.Resize(ngroups);
  for(int g=0; g< ngroups; g++) {
    store->eibasis(g).Resize(eibasis_save(g).GetDim(1), eibasis_save(g).GetDim(2), 
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

void Jastrow2_wf::storeParmIndVal(Wavefunction_data * dataptr, Sample_point *
				  sample, 
                               int e, Array1 <doublevar> & vals){
  //do nothing for now..
  if(parent->nparms()==0 ) { 
    //cout << " saving " << " valsize " << vals.GetDim(0) << endl;
    Wf_return newval(1,1);
    updateVal(parent, sample);
    getVal(parent, e,newval);
    vals(0)=newval.amp(0,0);
  }
    
}

//----------------------------------------------------------
void Jastrow2_wf::getParmDepVal(Wavefunction_data * dataptr,
                             Sample_point * sample,
                             int e,
                             Array1 <doublevar> & oldvals,
                             Wf_return & newval) {
  if(parent->nparms()==0 ) { 
    assert(oldvals.GetDim(0)>=1);
    assert(newval.amp.GetDim(1)>=1);
    assert(newval.amp.GetDim(0)>= 1);
    newval.amp(0,0)=oldvals(0);
    newval.phase(0,0)=0;
  }
  else { 
    updateVal(dataptr, sample);
    getVal(dataptr, e, newval);
  }

}

//--------------------------------------------------------------------------

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
    if(parent->group(g).hasThreeBodySpin())
      return 0;
  }
  Parm_deriv_return retparm;
  for(int g=0; g< ng; g++) {
    Array3 <doublevar> eionbasis(nelectrons,parent->natoms, maxeibasis); //for 3-body terms
    Array3 <doublevar> eibasis(parent->natoms, maxeibasis ,5);
    if(parent->group(g).hasOneBody() || parent->group(g).hasThreeBody()) { 
      Parm_deriv_return tmp_parm;
      int np=parent->group(g).one_body.nparms();
      tmp_parm.gradient.Resize(np);
      tmp_parm.gradient=0;
      tmp_parm.hessian.Resize(np,np);
      tmp_parm.hessian=0;
      for(int e=0; e< nelectrons; e++) { 
        parent->group(g).updateEIBasis(e,sample,eibasis);
        for(int at=0; at< parent->natoms; at++) { 
          for(int b=0; b < maxeibasis; b++) { 
            eionbasis(e,at,b)=eibasis(at,b,0);
          }
        }
        parent->group(g).one_body.getParmDeriv(e,eibasis,tmp_parm);
      }
      for(int i=0; i< np; i++) tmp_parm.hessian(i,i)=tmp_parm.gradient(i)*tmp_parm.gradient(i);
      extend_parm_deriv(retparm,tmp_parm);
    }
    
    Array3 <doublevar> eetotal(nelectrons, nelectrons, maxeebasis);
    eetotal=-1;
    Array3 <doublevar> eebasis(nelectrons, maxeebasis, 5);
    for(int e=0; e< nelectrons; e++) { 
      parent->group(g).updateEEBasis(e,sample, eebasis);
      for(int j=0; j< e; j++) { 
        for(int b=0; b< maxeebasis; b++) {
          eetotal(j,e,b)=eebasis(j,b,0);
        }
      }
    }
    //this part should be redone so that we can do updates instead of 
    //just calculating everything from scratch.
    if(parent->group(g).hasTwoBody()&& parent->group(g).two_body->nparms()) { 
      Parm_deriv_return tmp_parm;
      int np=parent->group(g).two_body->nparms();
      tmp_parm.gradient.Resize(np);
      tmp_parm.hessian.Resize(np,np);
      tmp_parm.gradient=0;
      tmp_parm.hessian=0;

      parent->group(g).two_body->getParmDeriv(eetotal, tmp_parm);
      extend_parm_deriv(retparm,tmp_parm);
    }
    
    
    if(parent->group(g).hasThreeBody() && parent->group(g).three_body.nparms()) { 
      Parm_deriv_return tmp_parm;
      int np=parent->group(g).three_body.nparms();
      tmp_parm.gradient.Resize(np);
      tmp_parm.hessian.Resize(np,np);
      tmp_parm.gradient=0;
      tmp_parm.hessian=0;
      parent->group(g).three_body.getParmDeriv(eionbasis,eetotal, tmp_parm);
      extend_parm_deriv(retparm,tmp_parm);
    }
  }
  parm_deriv=retparm;
  
  
  return 1;
}

//--------------------------------------------------------------------------

