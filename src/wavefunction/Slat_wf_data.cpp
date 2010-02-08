/*
 
Original Copyright (C) 2007 Lucas K. Wagner

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
#include "FSlat_wf.h"
#include <algorithm>
#include "MatrixAlgebra.h"
#include <map>
#include <utility>  // for make_pair


/*!
*/
void Slat_wf_data::read(vector <string> & words, unsigned int & pos,
                        System * sys)
{


  vector <string> strdetwt;
  vector <string> strstates;
  vector <vector <string> > statevec;
  unsigned int startpos=pos;
  
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

  sort=1;
  if(haskeyword(words, pos=startpos, "NOSORT")) {
    sort=0;
  }

  use_iterative_updates=0;
  if(haskeyword(words, pos=startpos, "ITERATIVE_UPDATES")) {
    use_iterative_updates=1;
  }

  all_weights=0;
  if(haskeyword(words, pos=startpos, "ALL_WEIGHTS")) {
    all_weights=1;
  }

  pos=startpos;
  vector <vector <string> > csfstr;
  vector <string> csfsubstr;

  while(readsection(words, pos, csfsubstr , "CSF")){
      csfstr.push_back(csfsubstr);
  }
  ncsf=csfstr.size();
  if(ncsf){
    //cout <<" Using CSF form for determinant weights "<<endl;
    use_csf=1;
    CSF.Resize(ncsf);
    int counter=0;
    for(int csf=0;csf<ncsf;csf++){
      if(csfstr[csf].size()<2)
        error(" Wrong number of elements in the CSF number ",csf+1);
      CSF(csf).Resize(csfstr[csf].size());
      for(int j=0;j<CSF(csf).GetDim(0);j++){
        CSF(csf)(j)=atof(csfstr[csf][j].c_str());
        if(j>0)
          counter++;
      }
    }
    ndet=counter;
    //cout <<" total number of determinants "<<ndet<<endl;

    if(readsection(words, pos, strdetwt, "DETWT")){
      error("Already using suplied CSF's, remove DETWT");
    }
    /*
    for(int csf=0;csf<ncsf;csf++){
      cout <<" CSF { ";
      for(int j=0;j<CSF(csf).GetDim(0);j++){
        cout <<CSF(csf)(j)<<"  ";
      } 
      cout <<"} "<<endl;
    }
    */
    
    counter=0;
    detwt.Resize(ndet);
    for(int csf=0;csf<ncsf;csf++)
      for(int j=1;j<CSF(csf).GetDim(0);j++){
        detwt(counter++)=CSF(csf)(0)*CSF(csf)(j);
      }
  }
  else{
    //cout <<" Using standard  DETWT form for determinant weights "<<endl;
    use_csf=0;
    pos=startpos;
    readsection(words, pos, strdetwt, "DETWT");
    ndet=strdetwt.size();
    detwt.Resize(ndet);
    for(int det=0; det < ndet; det++){
      detwt(det)=atof(strdetwt[det].c_str());
    }
  }

  //no sorting when ndet=1;
  if( ndet==1 && sort)
    sort=0;

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
    if(use_iterative_updates) error("Iterative updates not supported with complex MO's yet (easy code for you to write!)");
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



  //PK Reorder occupations and build work indexes for iterative updates

  // Switch excitation ordering to "natural" ordering"
  // Iterative update code must not accidentally construct a zero value determinant with two equal columns
  if (use_iterative_updates) {
    for (int s=0; s<2; s++)
      {
      for (int f=0; f<nfunc; f++)
	{
	for(int det=1; det<ndet; det++)
	  {

	    int swap=0;

	    // Characterize initial order

	    // Count the number of 'new' (novel) orbitals referenced by determinant det
	    int novel=0,badgs=0,goodgs=0;
	    Array1 <int> isnovel(nelectrons(s)),shouldbe(nelectrons(s)),isgoodgs(nelectrons(s));

	    for (int mo=0; mo<nelectrons(s); mo++) 
            {
	      shouldbe(mo)=0;
	      isgoodgs(mo)=0;
	    }

	    for (int mo=0; mo<nelectrons(s); mo++) 
            {
	      int mo2=0;
	      int found=0;
	      isnovel(mo)=1;
	      while (mo2<nelectrons(s)&&!found)
	      {
		if (occupation_orig(f,0,s)(mo)==occupation_orig(f,det,s)(mo2)) { found=1; isnovel(mo2)=0; }
		mo2++;
	      }
	      if (!found) novel++;
	    }


	    for (int mo=0; mo<nelectrons(s); mo++) 
            {
	      if (!isnovel(mo)) 
              {
		if (occupation_orig(f,0,s)(mo)!=occupation_orig(f,det,s)(mo)) 
		  {
		    badgs++;

		    for (int mo2=0; mo2<nelectrons(s); mo2++) 
		    {
		      if (occupation_orig(f,0,s)(mo2)==occupation_orig(f,det,s)(mo))
		      {
			shouldbe(mo)=mo2;
			break;
                      }
		    }
                  }
		else
		  {
		    isgoodgs(mo)=1;
		    goodgs++;
		  }
	      }
	    }

#ifdef DEBUG_SORT
	    cout << endl << "Det "<<det << " Func " << f << " Spin " << s << " Novel " << novel << " GoodGS " << goodgs << " BadGS " << badgs << endl; 	    
	    cout << "Intial: "; 
	    for (int mot=0; mot<nelectrons(s); mot++) cout << occupation_orig(f,det,s)(mot)  << " ";
	    cout << endl;
	    cout << "GoodGS: " ; 
	    for (int mot=0; mot<nelectrons(s); mot++) cout << isgoodgs(mot)  << " ";
	    cout << endl;
	    cout << "IsNovl: " ; 
	    for (int mot=0; mot<nelectrons(s); mot++) cout << isnovel(mot)  << " ";
	    cout << endl;
#endif

	    // Solve a single "bad position of orbital that occurs in ground state" each iteration
	    // Careful bookeeping to ensure we update records of "good" ground orbitals
	    
	    while (badgs>0) 
	    {

	      // Find 
	      int mo;
	      for (mo=0; mo< nelectrons(s); mo++)
	      {
		if (!isgoodgs(mo)&&!isnovel(mo)) break;
	      }

	      int target=shouldbe(mo);

#ifdef DEBUG_SORT
	      cout << "Swapping indexes " << mo << " " << target << " are " << occupation_orig(f,det,s)(mo) << " " << occupation_orig(f,det,s)(target) << endl;
	      cout << "Currnt: ";
	      for (int mot=0; mot<nelectrons(s); mot++) cout << occupation_orig(f,det,s)(mot)  << " ";
	      cout << endl;
	      cout << "GoodGS: " ; 
	      for (int mot=0; mot<nelectrons(s); mot++) cout << isgoodgs(mot)  << " ";
	      cout << endl;
	      cout << "IsNovl: " ; 
	      for (int mot=0; mot<nelectrons(s); mot++) cout << isnovel(mot)  << " ";
	      cout << endl;
#endif
	      int tmporb=occupation_orig(f,det,s)(target);
	      occupation_orig(f,det,s)(target)=occupation_orig(f,det,s)(mo);
	      occupation_orig(f,det,s)(mo)=tmporb;
	      
	      int tmpshould=shouldbe(mo);
	      shouldbe(mo)=shouldbe(target);
	      shouldbe(target)=tmpshould;
	      
	      int tmpnovel=isnovel(mo);
	      isnovel(mo)=isnovel(target);
	      isnovel(target)=tmpnovel;
	      
	      if (occupation_orig(f,0,s)(mo)==occupation_orig(f,det,s)(mo)) 
	      {
		isgoodgs(mo)=1;
		badgs--;
	      }
	      else
	      {
		isgoodgs(mo)=0;
	      }
	      isgoodgs(target)=1;
	      badgs--;
	      
	      if (swap) { swap=0; } else { swap=1; }

#ifdef DEBUG_SORT
	      cout << "Now   : ";
	      for (int mot=0; mot<nelectrons(s); mot++) cout << occupation_orig(f,det,s)(mot)  << " ";
	      cout << endl;
	      cout << "GoodGS: " ; 
	      for (int mot=0; mot<nelectrons(s); mot++) cout << isgoodgs(mot)  << " ";
	      cout << endl;
	      cout << "IsNovl: " ; 
	      for (int mot=0; mot<nelectrons(s); mot++) cout << isnovel(mot)  << " ";
	      cout << endl;
#endif	    
	    }

	    if (swap) detwt(det)=-detwt(det);

#ifdef DEBUG_SORT
	    cout << "Det "<<det << " Func " << f << " Spin " << s << " Swap " << swap << " Final swapped order " << endl;    
	    for (int mo=0; mo<nelectrons(s); mo++) cout << occupation_orig(f,det,s)(mo)  << " ";
	    cout << endl << endl;
#endif
	  }
	}
      }
  }

  occupation_nchanges.Resize(nfunc,ndet,2);
  evaluation_order.Resize(nfunc,ndet,2);
  occupation_first_diff_change_from_last_det.Resize(nfunc,ndet,2);
  occupation_changes.Resize(nfunc,ndet,2);

  max_occupation_changes=0;
  for (int s=0; s<2; s++)
  {
    for (int f=0; f<nfunc; f++)
    {
      for(int det=0; det<ndet; det++)
      {
	int nchanges=0;
	for (int mo=0; mo<nelectrons(s); mo++)
        {
	  if (occupation_orig(f,det,s)(mo)!=occupation_orig(f,0,s)(mo)) nchanges++;
        }
	occupation_nchanges(f,det,s)=nchanges;
	occupation_changes(f,det,s).Resize(nchanges);
	if (nchanges>max_occupation_changes) max_occupation_changes=nchanges;

	int ichange=0;
	for (int mo=0; mo<nelectrons(s); mo++)
        {
	  if (occupation_orig(f,det,s)(mo)!=occupation_orig(f,0,s)(mo))
	  {
	    occupation_changes(f,det,s)(ichange)=mo;
	    ichange++;
	  }
        }
      }
    }
  }


  // Compute (heuristic) good order for determinant updates
  for (int s=0; s<2; s++)
  {
    for (int f=0; f<nfunc; f++)
    {
      Array1 <int> order(ndet);
      for(int det=0; det<ndet; det++) order(det)=det;
      int dstart=1;
      int dend=ndet-1;
      int dlevel=0;
      calc_determinant_evaluation_order(f,s,occupation_changes,occupation_nchanges,occupation_orig,order,dstart,dend,dlevel);
      for (int idx=0; idx<ndet; idx++) evaluation_order(f,idx,s)=order(idx);
    }
  }

  // Avoid recomputing common updates: note index of last common excitation from previous determinant
  for (int s=0; s<2; s++)
  {
    for (int f=0; f<nfunc; f++)
    {
      int detlast=0;
      for(int idx=0; idx<ndet; idx++)
      {
	int det=evaluation_order(f,idx,s);
	occupation_first_diff_change_from_last_det(f,det,s)=0;
	if (det>0)
	{
	  for (int n=0;n<min(occupation_nchanges(f,detlast,s),occupation_nchanges(f,det,s));n++)
	  {
	    if (occupation_changes(f,det,s)(n)==occupation_changes(f,detlast,s)(n)) 
	    {
	      if (occupation_orig(f,det,s)(occupation_changes(f,det,s)(n))==occupation_orig(f,detlast,s)(occupation_changes(f,detlast,s)(n)))
              {
		occupation_first_diff_change_from_last_det(f,det,s)=n+1;
	      }
	      else
              {
		break;
	      }
	    }
	  }
	}
	detlast=det;
      }
    }
  }

  
  // Check excitation ordering are permissable
  // Iterative update code must not accidentally construct a zero value determinant with two equal columns
  if (use_iterative_updates) {
    for (int s=0; s<2; s++)
      {
      for (int f=0; f<nfunc; f++)
	{
	for(int det=1; det<ndet; det++)
	  {
	  if (occupation_nchanges(f,det,s)>1) // Single excitations out of order are safe
	    {
	    for (int mo=0; mo<nelectrons(s); mo++)
	      {
	      if (occupation_orig(f,det,s)(mo)!=occupation_orig(f,0,s)(mo))
		{
		  // This is an excitation from the zeroth determinant. Must not exist elsewhere is ground state

		for (int mp=0; mp<nelectrons(s); mp++)
		  {
		  if (occupation_orig(f,det,s)(mo)==occupation_orig(f,0,s)(mp))
		    {	    
		      cout << "Orbital " << mo << " in determinant " << det << " function " <<  f << " spin " << s << " is also orbital " << mp << " in first determinant" << endl;
		      cout << "Fast update procedure requires 'natural ordering' with no 'moved' ground state orbitals in excited determinants" << endl;
		      error("Orbital ordering incorrect in determinant ",det);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
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
    if (use_iterative_updates) {
      wf=new FSlat_wf;
      FSlat_wf * slatwf;
      recast(wf, slatwf);
      slatwf->init(this);
      attachObserver(slatwf);
    } else {
      wf=new Slat_wf;
      Slat_wf * slatwf;
      recast(wf, slatwf);
      slatwf->init(this);
      attachObserver(slatwf);
    }
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
  {
    os << ndet << " determinants\n";
  }

  else
    os << "1 determinant" << endl;

  if(use_iterative_updates)
  {
    os << "Using fast/low memory iterative updates for multideterminants" << endl
       << "Ref: Phani K. V. V. Nukala and P. R. C. Kent. J. Chem. Phys. 130 204105 (2009)" << endl;
  }
  else
  {
    // Advice check for iterative updates
    if (ndet>1)
    {
      if (nelectrons(0)>10||nelectrons(1)>10)
      {
	os << "Advice: You have more than 1 determinant and more than 10 electrons in any one spin" << endl
	   << "        ITERATIVE_UPDATES are likely faster and will use less memory" << endl;
      }
    }
  }

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
	// PK Print changes for iterative updates
	os << "Changes: "<< occupation_nchanges(f,det,s) <<" First diff: "<<occupation_first_diff_change_from_last_det(f,det,s) <<endl;
	for (int ichange=0; ichange<occupation_nchanges(f,det,s); ichange++)
	{
	  os << occupation_changes(f,det,s)(ichange)+1 << " ";
	}
	os << endl ;
	// 
      }
    }


    // Output ordering heuristic update costs
    Array1 <int> order(ndet);

    for(int det=0; det<ndet; det++) order(det)=det;

    for(int s=0; s<2; s++)
      {
        if(s==0)
          os << "Update cost, spin up (natural order): ";
        if(s==1)
          os << "Update cost, spin down (natural order): ";
	os << cost_iterative_determinant_order(f,s,occupation_changes,occupation_nchanges,occupation_orig,order) << endl ;
      }

    for(int s=0; s<2; s++)
      {
	for(int det=0; det<ndet; det++) order(det)=evaluation_order(f,det,s);

        if(s==0)
          os << "Update cost, spin up (heuristic sort): ";
        if(s==1)
          os << "Update cost, spin down (heuristic sort): ";
	os << cost_iterative_determinant_order(f,s,occupation_changes,occupation_nchanges,occupation_orig,order) << endl ;
      }
  }

  os << "Maximum occupation changes : " << max_occupation_changes << endl;

  os << "Molecular Orbital object : ";
  genmolecorb->showinfo(os);

  return 1;
}

//----------------------------------------------------------------------
// Compute the heuristic "best" order for evaluating multiderminant expansion
// Added by Paul Kent / CNMS / ORNL 2009
//
// Called recursively to find best order for determinants numbered
// inclusively between dstart and dend. Each level of recursion processes
// one "level" dlevel of excitations (singles, doubles etc.)
//
// Algorithm: at each level, sort by the most common excitation
// for each group of common excitations, sort at the next level
//
// To hash excitation type use orbital index*no. orbitals+replacement orbital
//
// nchanges(f,det,s)
// changes(f,det,s)(ichange)=mo
// occ(f,det,s)(mo) from occ(f,0,s)(mo)

void Slat_wf_data::calc_determinant_evaluation_order(const int & f, const int & s, 
						   const Array3<Array1<int> > & changes,
						   const Array3 <int> & nchanges,
						   const Array3< Array1 <int> > & occ,
						   Array1 <int> & order,
						   int dstart, int dend, int dlevel ) {

  // Ensure no-ops for single entry dstart=dend
  if (dend>dstart) {

    typedef std::multimap<int, int> MultiMap;
    MultiMap excitation_counts;
    for (int i=dstart;i<=dend;i++)
      {
	int det=order(i);
	int nch=nchanges(f,det,s);
	int hash=-1;
	if (nch>dlevel) {
	  int ch=changes(f,det,s)(dlevel);
	  hash=occ(f,0,s)(ch)*nmo+occ(f,det,s)(ch);
	}
	excitation_counts.insert(make_pair(hash,det));
      }
    /*
    cout << "Excitation_counts" << endl;
    for (multimap<int,int>::iterator ii=excitation_counts.begin(); ii!=excitation_counts.end(); ++ii)
      {	cout << "ExcH " << (*ii).first << " det " << (*ii).second << endl; }
    */
    MultiMap excitation_weights;
    MultiMap::iterator iiter=excitation_counts.begin();
    MultiMap::iterator jjter=excitation_counts.end();
    while (iiter!=jjter) {
      int hash=iiter->first;
      int count=excitation_counts.count(hash);
      excitation_weights.insert(make_pair(count,hash));
      std::pair<MultiMap::iterator, MultiMap::iterator> iterpair = excitation_counts.equal_range(hash);
      iiter=iterpair.second;
    }
    /*
    cout << "Excitation_weights" << endl;
    for (multimap<int,int>::reverse_iterator ii=excitation_weights.rbegin(); ii!=excitation_weights.rend(); ++ii)
      {	cout << "Count " << (*ii).first << " ExcH " << (*ii).second << endl; }
    */
    int oout=dstart;
    
    for (MultiMap::reverse_iterator iio=excitation_weights.rbegin(); iio!=excitation_weights.rend(); ++iio)
      {
	int hash=iio->second;
	std::pair<MultiMap::iterator, MultiMap::iterator> iterpair = excitation_counts.equal_range(hash);
	int next_level_start=oout;
	for (multimap<int,int>::iterator jjo=iterpair.first; jjo!=iterpair.second; ++jjo)
	  {
	    order(oout++)=jjo->second;
	  }
	int next_level_end=oout-1;
	int next_level=dlevel+1;
	if (hash!=-1) { // No need to subsort "no change" group
	    if (next_level_end>next_level_start+1) 
	      calc_determinant_evaluation_order(f,s,changes,nchanges,occ,order,next_level_start,next_level_end,next_level);	
	}
      }
  }
}

//----------------------------------------------------------------------
// Simple costing of determinant evaluation order
// Counts total O(kN) of iterative changes to determinants evaluated in order()
// Neglects important prefactors: better to assess actual operations counts based on implementation
int Slat_wf_data::cost_iterative_determinant_order(const int & f, const int & s, 
						   const Array3<Array1<int> > & changes,
						   const Array3 <int> & nchanges,
						   const Array3< Array1 <int> > & occ,
						   const Array1 <int> & order) {
  int cost=0;

  int prev_det=order(0);
  for (int det=1;det<ndet; det++)
    {
      int actual_det=order(det);
      // Count different excitations from previous det      
      int first=0;
      for (int n=0;n<min(nchanges(f,prev_det,s),nchanges(f,actual_det,s));n++)
	{
	  if (changes(f,actual_det,s)(n)==changes(f,prev_det,s)(n)) 
	    {
	      if (occ(f,actual_det,s)(changes(f,actual_det,s)(n))==occ(f,prev_det,s)(changes(f,prev_det,s)(n)))
		{
		  first=n+1;
		}
	      else
		{
		  break;
		}
	    }
	}
      prev_det=actual_det;
      // diffs differences total from the ground state
      // first common differences from last det
      for (int iter=first;iter<nchanges(f,actual_det,s);iter++)
	{
	  cost+=iter+1;
	}
    }
  return cost;
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
  
  if(!sort)
    os << indent << "NOSORT" << endl;

  if(use_iterative_updates)
    os << indent << "ITERATIVE_UPDATES" << endl;

  if(all_weights)
    os << indent << "ALL_WEIGHTS" << endl;

  if(use_csf){
    Array1 <Array1 <doublevar> > CSF_print(ncsf);
    Array1 <doublevar> detwt_print(ndet);
    Array3 < Array1 <int> > occupation_orig_print(nfunc,ndet,2);
    
    if(sort){
      Array1 <doublevar> csf_tmp(ncsf);
      Array1 <int> list;
      for(int csf=0;csf<ncsf;csf++)
	csf_tmp(csf)=CSF(csf)(0);
      sort_abs_values_descending(csf_tmp,csf_tmp,list);
      
      /*
      cout <<"CSF list"<<endl;
      for(int csf=0;csf<ncsf;csf++){
	cout <<csf<<"  "<<list(csf)<<endl;
      }
      */

      Array1 < Array1 <int> > det_pos(ncsf);
      int counterr=0;
      for(int csf=0;csf<ncsf;csf++){
	det_pos(csf).Resize(CSF(csf).GetDim(0)-1);
	for(int j=1;j<CSF(csf).GetDim(0);j++){
	  det_pos(csf)(j-1)=counterr++;
	}
      }

      
      //preparing CSF_print
      int counter_new=0;
      for(int csf=0;csf<ncsf;csf++){
	CSF_print(csf).Resize(CSF(list(csf)).GetDim(0));
	for(int j=0;j<CSF(list(csf)).GetDim(0);j++){   
	  CSF_print(csf)(j)=CSF(list(csf))(j);
	  if(j>0){
	    detwt_print(counter_new++)=CSF_print(csf)(0)*CSF_print(csf)(j);
	  }
	}
      }
      
      Array1 <int> detlist(ndet);
      //cout <<"CSF list"<<endl;
      int counter_old=0;
      for(int csf=0;csf<ncsf;csf++){
	//cout << csf<<"  ";
	for(int j=0;j<CSF(list(csf)).GetDim(0)-1;j++){
	  //cout <<det_pos(list(csf))(j) <<"  ";
	  detlist(counter_old++)=det_pos(list(csf))(j);
	}
	//cout <<endl;
      }

      /*
      cout <<"Determinant list"<<endl;
      for(int det=0; det < ndet; det++){
	cout <<det<<"  "<<detlist(det)<<endl;
      }
      */
      
      //preparing occupation_orig_print
      for(int f=0; f< nfunc; f++)
	for(int det=0; det < ndet; det++)
	  for(int s=0; s<2; s++){
	    occupation_orig_print(f,det,s).Resize(nelectrons(s));
	    occupation_orig_print(f,det,s)=occupation_orig(f,detlist(det),s);
	  }
    }
    else{ //no sorting 
      int counter=0;
      for(int csf=0;csf<ncsf;csf++){
	CSF_print(csf).Resize(CSF(csf).GetDim(0));
	for(int j=0;j<CSF(csf).GetDim(0);j++){   
	  CSF_print(csf)(j)=CSF(csf)(j);
	  if(j>0)
	    detwt_print(counter++)=CSF(csf)(0)*CSF(csf)(j);
	}
      }
      for(int f=0; f< nfunc; f++)
	for(int det=0; det < ndet; det++)
	  for(int s=0; s<2; s++){
	    occupation_orig_print(f,det,s).Resize(nelectrons(s));
	    occupation_orig_print(f,det,s)=occupation_orig(f,det,s);
	  }
    }
    // do printout
    
    for(int csf=0;csf<ncsf;csf++){
      os << indent<<" CSF { ";
      for(int j=0;j<CSF_print(csf).GetDim(0);j++){
        os <<CSF_print(csf)(j)<<"  ";
      } 
      os <<"} "<<endl;
    }
    for(int f=0; f< nfunc; f++){
    os << indent << "STATES { " << endl << indent <<"  ";
    for(int det=0; det < ndet; det++){
      if(ndet>1)
	os <<"#  Determinant "<<det+1<<": weight: "<<detwt_print(det)<<endl<< indent <<"  ";
      for(int s=0; s<2; s++){
        for(int e=0; e< nelectrons(s); e++){
          os << occupation_orig_print(f,det,s)(e)+1 << " ";
          if((e+1)%20 ==0)
            os << endl << indent << "  ";
        }
        os << endl << indent << "  ";
      }
    }
    os << "}" << endl;
    }
    
  }
  else {
    Array1 <doublevar> detwt_print(ndet);
    Array3 < Array1 <int> > occupation_orig_print(nfunc,ndet,2);
    
    if(sort){
      Array1 <int> list;
      sort_abs_values_descending(detwt,detwt_print,list);
     
      for(int f=0; f< nfunc; f++)
	for(int det=0; det < ndet; det++)
	  for(int s=0; s<2; s++){
	    occupation_orig_print(f,det,s).Resize(nelectrons(s));
	    occupation_orig_print(f,det,s)=occupation_orig(f,list(det),s);
	  }
    }
    else{ //no sorting 
      detwt_print=detwt;
      for(int f=0; f< nfunc; f++)
	for(int det=0; det < ndet; det++)
	  for(int s=0; s<2; s++){
	    occupation_orig_print(f,det,s).Resize(nelectrons(s));
	    occupation_orig_print(f,det,s)=occupation_orig(f,det,s);
	  }
    }
    
    //do printout
    os << indent << "DETWT { ";
    for(int det=0; det < ndet; det++)
      {
	os << detwt_print(det) << "  ";
	if ((det+1)%6==0)
	  os <<endl<<indent;
      }
    os << "}" << endl;
    
    for(int f=0; f< nfunc; f++){
    os << indent << "STATES { " << endl << indent <<"  ";
    for(int det=0; det < ndet; det++){
      if(ndet>1)
	os <<"#  Determinant "<<det+1<<": weight: "<<detwt_print(det)<<endl<< indent <<"  ";
      for(int s=0; s<2; s++){
        for(int e=0; e< nelectrons(s); e++){
          os << occupation_orig_print(f,det,s)(e)+1 << " ";
          if((e+1)%20 ==0)
            os << endl << indent << "  ";
        }
        os << endl << indent << "  ";
      }
    }
    os << "}" << endl;
    }
  }//if using determinants
      

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
    if(use_csf){
      if(all_weights){
	parms.Resize(ncsf);
	for(int i=0; i< ncsf; i++) {
	  parms(i)=CSF(i)(0);
	}
      }
      else{
	parms.Resize(ncsf-1);
	for(int i=1; i< ncsf; i++) {
	  parms(i-1)=CSF(i)(0);
	}
      }
    }
    else{//just independent weights
      if(all_weights){
	parms.Resize(detwt.GetDim(0));
	for(int i=0; i< detwt.GetDim(0); i++) {
	  parms(i)=detwt(i);
	}
      }
      else{
	parms.Resize(detwt.GetDim(0)-1);
	for(int i=1; i< detwt.GetDim(0); i++) {
	  parms(i-1)=detwt(i);
	}
      }
    }
  }
  else {
    parms.Resize(0);
  }

  //for(int i=0;i<parms.GetDim(0);i++)
  // cout <<parms(i)<<endl;
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
    if(use_csf){
      if(all_weights)
	for(int csf=0; csf< ncsf; csf++) 
	  CSF(csf)(0)=parms(csf);
      else
	for(int csf=1; csf< ncsf; csf++) 
	  CSF(csf)(0)=parms(csf-1);
      //cout << CSF(csf)(0)<<endl;
      int counter=0;
      for(int csf=0; csf< ncsf; csf++) {
	for(int j=1;j<CSF(csf).GetDim(0);j++){
	  detwt(counter++)=CSF(csf)(0)*CSF(csf)(j);
	}
      }
      assert(counter==ndet);
      // for(int i=0;i<ndet;i++)
      // cout <<detwt(i)<<" ";
      //cout <<endl;
    }
    else{
      if(all_weights)
	for(int i=0; i< detwt.GetDim(0); i++) 
	  detwt(0)=parms(i);
      else
	for(int i=1; i< detwt.GetDim(0); i++) 
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

void Slat_wf_data::renormalize(){

  if(optimize_det && all_weights) {
    doublevar norm=0.0;
    if(use_csf){
      for(int csf=0; csf< ncsf; csf++) 
	norm+=CSF(csf)(0)*CSF(csf)(0);
    }
    else{
      for(int i=0; i< detwt.GetDim(0); i++) 
	norm+=detwt(i)*detwt(i);
    }
    doublevar factor=1.0/sqrt(norm);
    if(use_csf){
      for(int csf=0; csf< ncsf; csf++) 
	CSF(csf)(0)*=factor;
    }
    else{
      for(int i=0; i< detwt.GetDim(0); i++) 
	detwt(i)*=factor;
    }
  }
}
