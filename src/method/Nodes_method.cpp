/*
 
Copyright (C) 2007 Michal Bajdich
Copyright (C) 2012 Shuming Hu

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

#include "Nodes_method.h"
#include "qmc_io.h"
#include "System.h"
#include "Program_options.h"
/*!
Read the "words" from the method section in the input file
via doinput() parsing, and store section information in private
variables isoses, resolution, and minmax.
Set up MO_matrix and Sample_point objects for wavefunction from input
*/
void Nodes_method::read(vector <string> words,
                       unsigned int & pos,
                       Program_options & options)
{
  pos=0;	//always start from first word
  doublevar Tres;
  vector <string> Torbs;
  vector <string> Tminmax;
  vector <string> orbtext;
  vector <string> Tdxyz;


  allocate(options.systemtext[0], sysprop);
  sysprop->generateSample(mywalker);

  allocate(options.twftext[0], sysprop, wfdata);
  wfdata->generateWavefunction(wf);

  mywalker->attachObserver(wf);
 
  pos=0;
  if(readsection(words,pos=0,Torbs,"CONTOURS")){
    //   error("Need CONTOURS in METHOD section");
  plots.Resize(Torbs.size());
  for(unsigned int i=0; i<Torbs.size(); i++){
        plots(i)=atoi(Torbs[i].c_str());
      }
  }
  else {
    cout <<"WARNING: All electrons are used for scanning!"<<endl;
    plots.Resize(mywalker->electronSize());
    for(int i=0; i<plots.GetSize(); i++){
      plots(i)=i+1;
      } 
  }
 
  pos=0;
  if(! readvalue(words,pos,Tres,"RESOLUTION"))
    error("Need RESOLUTION in METHOD section");
  resolution=Tres;


  pos=0;
  minmax.Resize(6);
  if(readsection(words,pos,Tminmax,"MINMAX")) { 
    if(Tminmax.size() != 6)
      error("MINMAX needs 6 values");
    for(unsigned int i=0; i<Tminmax.size(); i++)
    {
      minmax(i)=atof(Tminmax[i].c_str());
    }
  }
  else  { 
    minmax=0.0;
  
    int nions=sysprop->nIons();
    Array1 <doublevar> ionpos(3);
    for(int i=0; i< nions; i++) {
      sysprop->getIonPos(i,ionpos);
      for(int d=0; d< 3; d++) {
        if(ionpos(d) < minmax(2*d)) minmax(2*d) = ionpos(d);
        if(ionpos(d) > minmax(2*d+1)) minmax(2*d+1)=ionpos(d);
      }
    }

    for(int d=0; d< 3; d++) {
      minmax(2*d)-=4.0;
      minmax(2*d+1)+=4.0;
      cout << "minmax " << minmax(2*d) << " " << minmax(2*d+1) << endl;
    }
  } 
  
  
 

  // Makes the Nodes methods use the current implementation of Nodes
  if(readvalue(words, pos=0, readconfig, "READCONFIG")){
    Array1 <Config_save_point> config_pos;
    config_pos.Resize(0);
    if(readconfig!="") { 
      read_configurations(readconfig, config_pos);
    }
    if(config_pos.GetDim(0)<1)
      error("Could not read a single walker from config file");
    //take a first one
    config_pos(0).restorePos(mywalker);
  }
  else {
    mywalker->randomGuess();
    debug_write(cout, 0, " configs read ", 1,
                " configs randomly generated \n");
  }
  
  pos=0;
  doublemove=false;
  if( readsection(words,pos,Tdxyz,"DOUBLEMOVES")){
    doublemove=true;
    if(plots.GetSize()%2)
      error("Needs pairs of countours for doublemoves");
    unsigned int dummy=3*plots.GetSize()/2;
    if(Tdxyz.size()!=dummy)
      error("DOUBLEMOVES needs 3 x # of vectors values");
    dxyz.Resize(plots.GetSize()/2,3);
    for(int i=0; i<plots.GetSize()/2; i++){
      for (int j=0; j<3;j++)
        dxyz(i,j)=atof(Tdxyz[i*3+j].c_str());
    }
  }
}


void Nodes_method::run(Program_options & options, ostream & output)
{
  ofstream os; //for writing to *.plt files
  string pltfile; //name of plotfile being written
  string confile; //name of configuration of electrons to be  being written
  double max_value,min_value;
  int count;
  Array1 <doublevar> xyz(3), xyz2(3),
    resolution_array(3); //position of electron "in" MO
  Array1 <int> D_array1(3); //dummy array1
  Array2 <doublevar> oldpos(mywalker->electronSize(),3);
  Array1 <doublevar> tmp(3);
  D_array1=0; //sets all 3 components to 0. use as counter for gridpoints
  
  
  D_array1(0)=roundoff((minmax(1)-minmax(0))/resolution);
  D_array1(1)=roundoff((minmax(3)-minmax(2))/resolution);
  D_array1(2)=roundoff((minmax(5)-minmax(4))/resolution);

  resolution_array(0)=(minmax(1)-minmax(0))/(D_array1(0)-1);  
  resolution_array(1)=(minmax(3)-minmax(2))/(D_array1(1)-1);
  resolution_array(2)=(minmax(5)-minmax(4))/(D_array1(2)-1);

  Array2 <doublevar> grid(plots.GetSize(),D_array1(0)*D_array1(1)*D_array1(2));
  Wf_return wfvals(wf->nfunc(), 2); //where wfval fist index labels 
    //wf number ie. ground state
    // excited state and so on, and second label is for two values, sign and log(wf).

 //scan electron positions
  cout << "Using these electron positions:"<<endl;
  for(int i=0; i<mywalker->electronSize(); i++){
    mywalker->getElectronPos(i, tmp);
    cout.precision(5);
    cout.width(8);
    cout <<i+1<<" "<<tmp(0)<<" "<<tmp(1)<<" "<<tmp(2)<<endl;
    for (int d=0;d<3;d++)
      oldpos(i,d)=tmp(d);
    //    cout << "pos " << oldpos(i,0) << "   " << oldpos(i,1) << "   " << oldpos(i,2)
    //     << endl;
  }

  //generate .xyz file for gOpenMol to view coordinates
  pltfile=options.runid + ".xyz";
  os.open(pltfile.c_str());
  cout<<"writing to "<<pltfile<<endl;
  vector <string> atomlabels;
  sysprop->getAtomicLabels(atomlabels);
  os<<atomlabels.size()+mywalker->electronSize()<<endl;
  os<<endl;
  
  int spin_up=sysprop->nelectrons(0);
  //cout << "Up electrons "<<spin_up<<endl;

  
  for(unsigned int i=0; i<atomlabels.size(); i++){
    Array1 <doublevar> pos(3);
    mywalker->getIonPos(i, pos);
    os<<atomlabels[i] <<" "<< pos(0)
    <<" "<<pos(1)
    <<" "<< pos(2)<<endl;
  }
  for(int j=0; j<mywalker->electronSize(); j++){
    if (j<spin_up) os<<"Eu";
    else os<<"Ed";
        // mywalker->getElectronPos(j, pos);
       os<<" "<< oldpos(j,0)<<" "<<oldpos(j,1) <<" "<< oldpos(j,2)<<endl;
  }
  os.close();
  

  if (!doublemove) {
    //calculate value of each molecular orbital at each grid point and store in an Array1
    // grid values with x=fastest running variable, and z=slowest

    cout<<"calculating "<<D_array1(0)*D_array1(1)*D_array1(2)
        <<" grid points"<<endl;
    cout<<"for "<< plots.GetDim(0) <<" 3D projections of wavefunction"<<endl;
    count=0;
    xyz(0)=minmax(0);
    xyz(1)=minmax(2);
    xyz(2)=minmax(4); //init elec probe to xmin ymin zmin
    //Array1 <doublevar> oldpos(3)

  
  
    for(int i=0; i<plots.GetSize(); i++){
      cout.width(3);
      cout <<"============================="<<plots(i)
           <<"=============================="
           <<endl;
      count=0;
      for(int xx=0; xx<D_array1(0);xx++){
        xyz(0)=minmax(0)+xx*resolution_array(0);  //move forward on x axis one resolution unit
        max_value=min_value=0.0;
        cout << "x ";
        cout.precision(5);
        cout.width(8);
        cout<< xyz(0);
        for(int yy=0; yy<D_array1(1);yy++){
          xyz(1)=minmax(2)+yy*resolution_array(1);  //move forward on y axis one resolution unit
          for(int zz=0;zz<D_array1(2);zz++){
            xyz(2)=minmax(4)+zz*resolution_array(2); //move forward on z axis one resolution unit
            // mywalker->getElectronPos(plots(i), oldpos);
            mywalker->setElectronPos(plots(i)-1,xyz); //move elec#plots(i) to point specified by xyz
            wf->updateVal(wfdata, mywalker); //update wfdata
            wf->getVal(wfdata, 0, wfvals); //get wf value
            const doublevar cutoff=15;
            if(wfvals.amp(0,0)<cutoff) { 
              grid(i,count)=wfvals.sign(0)*exp(wfvals.amp(0,0));
            }
            else { //cut off the maximum value output so that there aren't overflow errors
              grid(i,count)=wfvals.sign(0)*exp(cutoff);
            }
            //grid(i,count)=exp(2.0*wfvals(0,1));//!square of wavefunction
            if (grid(i,count)>max_value) max_value=grid(i,count);
            if (grid(i,count)<min_value) min_value=grid(i,count);
            //  cout << "grid " << grid(i,count)<< " " << xyz(0) << endl;
            count++; //index for cycling through grid points
          }
        }
        cout.setf(ios::scientific| ios:: showpos);
        cout <<", max. value "<<max_value<<", min. value "<<min_value<< endl;
        cout.unsetf(ios::scientific| ios:: showpos);

      }
      mywalker->setElectronPos(plots(i)-1, oldpos(plots(i)-1));
    }


    //Loop through and generate plot files for each orbital requested
    if(plots.GetSize()<=0)
      error("Number of requested plots is not a positive number");
    cout<<"saving data for "<<plots.GetSize()<<" 3D projections of wavefunction"<<endl;
    for(int i=0; i<plots.GetSize(); i++)
      {
        //output to file with orbital number in it
        char strbuff[40];
        sprintf(strbuff, "%d", plots(i));
        confile=pltfile = options.runid;
        confile += ".orb.";
        pltfile += ".orb.";
        pltfile += strbuff;
        confile += strbuff;
        pltfile += ".cube"; /*FIGURE OUT HOW TO CONVERT INT TO STRING*/
        confile += ".xyz";
        
        os.open(pltfile.c_str());
        cout<<"writing to "<<pltfile<<endl;

        os << "QWalk nodes output\n";
        os << "Wavefunction single scan with " <<   plots(i) <<" electron"<<endl;
        int natoms=sysprop->nIons();
        os << "  " << natoms+ mywalker->electronSize()-1 << "   " << minmax(0) << "   "
          << minmax(2) << "   " << minmax(4) << endl;
        os << D_array1(0) << "   " << resolution_array(0) << "  0.0000   0.0000" << endl;
        os << D_array1(1) << "   0.0000   " << resolution_array(1) << "  0.0000" << endl;
        os << D_array1(2) << "   0.0000    0.0000    " << resolution_array(2) << endl;
        Array1 <doublevar> pos(3);
        for(int at=0; at< natoms; at++) {
          mywalker->getIonPos(at,pos);
          os << "   " << mywalker->getIonCharge(at) << "   0.0000    " << pos(0) 
            <<"    " << pos(1) << "   " << pos(2) << endl;
        }
        for(int j=0; j<mywalker->electronSize(); j++){
          if (j!=plots(i)-1){
            if (j<spin_up) os<<"     1   0.0000    ";
            else os<<"     2     0.0000    ";
            os<< oldpos(j,0)<<"   "<<oldpos(j,1)<<"   "<< oldpos(j,2)<<endl;
          }
        }
        os.setf(ios::scientific);
        for(int j=0; j< D_array1(0)*D_array1(1)*D_array1(2); j++) {
          os <<setw(16)<<setprecision(8)<<grid(i,j);
          if(j%6 ==5) os << endl;
        }
        os << endl;
        os.unsetf(ios::scientific);
        os<<setprecision(6);
        os.close();

	
        /*
	  old plt file plots
        // http://www.csc.fi/gopenmol/developers/plt_format.phtml
        os<<"3 "; //rank=3 always
        os<<"2\n"; //dummy variable => "Orbital/density surface"
        //number of grid points for x, y, & z direction
        os <<D_array1(2)<<" "<<D_array1(1)<<" "<<D_array1(0)<<endl;
        os <<minmax(4)<<" "<<minmax(5)<<" "<<minmax(2)<<" "<<minmax(3)
           <<" "<<minmax(0)<<" "<<minmax(1)<<endl;
        
        for(int j=0; j<(D_array1(0)*D_array1(1)*D_array1(2)); j++)
          os<<grid(i,j)<<endl;
        os.close();
        
        os.open(confile.c_str());
        cout<<"writing to "<<confile<<endl;
        Array1 <doublevar> pos(3);
        sysprop->getAtomicLabels(atomlabels);
        os<<atomlabels.size()+ mywalker->electronSize()-1 <<endl;
        os << endl;
        for(unsigned int j=0; j<atomlabels.size();j++){
          mywalker->getIonPos(j, pos);
          os<<atomlabels[j] <<" "<< pos(0)
            <<" "<<pos(1)
            <<" "<< pos(2)<<endl;
        }
        for(int j=0; j<mywalker->electronSize(); j++){
          if (j!=plots(i)-1){
            if (j<spin_up) os<<"Eu";
            else os<<"Ed";
            // mywalker->getElectronPos(j, pos);
            os<<" "<< oldpos(j,0)<<" "<<oldpos(j,1)
              <<" "<< oldpos(j,2)<<endl;
          }
        }
        os.close();   
	*/
	
      }
    
  }
  else {
    cout << "Using translation vector(s): "<<endl;
    for(int i=0;i<dxyz.GetDim(0);i++)    
      cout<<"( "<<dxyz(i,0)<<" , "<<dxyz(i,1)<<" , "<<dxyz(i,2)<<" )"<<endl;

    for(int i=0; i<plots.GetSize(); i=i+2){
      count=0;
      xyz(0)=minmax(0);
      xyz(1)=minmax(2);
      xyz(2)=minmax(4); //init elec probe to xmin ymin zmin
      //Array1 <doublevar> oldpos(3)
      cout.width(3);
      cout <<"============================="<<plots(i)<<" and "<<plots(i+1)
           <<"=============================="
           <<endl;
      for(int xx=0; xx<D_array1(0);xx++){
        max_value=min_value=0.0;
        xyz(0)=minmax(0)+xx*resolution_array(0); 
        xyz2(0)=xyz(0)+dxyz(i/2,0);  
        cout << "x ";
        cout.precision(5);
        cout.width(8);
        cout<< xyz(0);
        for(int yy=0; yy<D_array1(1);yy++){
          xyz(1)=minmax(2)+yy*resolution_array(1); 
          xyz2(1)=xyz(1)+dxyz(i/2,1);
          for(int zz=0;zz<D_array1(2);zz++){
            xyz(2)=minmax(4)+zz*resolution_array(2); 
            xyz2(2)=xyz(2)+dxyz(i/2,2);
            // mywalker->getElectronPos(plots(i), oldpos);
            mywalker->setElectronPos(plots(i)-1,xyz);//move elec#plots(i) 
            //to point specified by xyz
            mywalker->setElectronPos(plots(i+1)-1,xyz2);//move elec#plots(i) 
            //to point specified by xyz
            wf->updateVal(wfdata, mywalker); //update wfdata
            wf->getVal(wfdata, 0, wfvals); //get wf value
            grid(i,count)=wfvals.sign(0)*exp(wfvals.amp(0,0));
            //grid(i,count)=exp(2.0*wfvals(0,1));//!square of wavefunction
            //  cout << "grid " << grid(i,count)<< " " << xyz(0) << endl;
            if (grid(i,count)>max_value) max_value=grid(i,count);
            if (grid(i,count)<min_value) min_value=grid(i,count);
            
            count++; //index for cycling through grid points
          }
        }
        cout.setf(ios::scientific| ios:: showpos);
        cout <<", max. value "<<max_value<<", min. value "<<min_value<< endl;
        cout.unsetf(ios::scientific| ios:: showpos);
      }
      
      mywalker->setElectronPos(plots(i)-1, oldpos(plots(i)-1));
      mywalker->setElectronPos(plots(i+1)-1, oldpos(plots(i+1)-1)); 
    }
    
    //Loop through and generate plot files for each orbital requested
    if(plots.GetSize()<=0)
      error("Number of requested plots is not a positive number");
    cout<<"saving data for "<<plots.GetSize()<<" 3D projections of wavefunction"<<endl;
    for(int i=0; i<plots.GetSize(); i=i+2)
      {
        char strbuff2[40],strbuff[40];
        //output to file with orbital number in it
        sprintf(strbuff, "%d", plots(i));
        sprintf(strbuff2, "%d", plots(i+1));
        confile=pltfile = options.runid;
        confile += ".orb.";
        pltfile += ".orb.";
        pltfile += strbuff;
        confile += strbuff;
        pltfile += "with";
        confile += "with";
        pltfile += strbuff2;
        confile += strbuff2;
        
        pltfile += ".cube"; /*FIGURE OUT HOW TO CONVERT INT TO STRING*/
        confile += ".xyz";
        
        os.open(pltfile.c_str());
        cout<<"writing to "<<pltfile<<endl;
	os << "GOS nodes output\n";
	os << "Wave function double scan with "<<   plots(i) <<" & "<<plots(i+1)<<" electrons"<< endl;
        int natoms=sysprop->nIons();
        os << "  " << natoms + mywalker->electronSize()-2 << "   " << minmax(0) << "   "
           << minmax(2) << "   " << minmax(4) << endl;
        os << D_array1(0) << "   " << resolution_array(0) << "  0.0000   0.0000" << endl;
        os << D_array1(1) << "   0.0000   " << resolution_array(1) << "  0.0000" << endl;
        os << D_array1(2) << "   0.0000    0.0000    " << resolution_array(2) << endl;
        Array1 <doublevar> pos(3);
        for(int at=0; at< natoms; at++) {
          mywalker->getIonPos(at,pos);
          os << "   " << mywalker->getIonCharge(at) << "   0.0000    " << pos(0) 
             <<"    " << pos(1) << "   " << pos(2) << endl;
        }
        for(int j=0; j<mywalker->electronSize(); j++){
          if ((j!=plots(i)-1)&&(j!=plots(i+1)-1)){
            if (j<spin_up) os<<"    3    0.0000    ";
            else os<<"    3    0.0000   ";
            os<< oldpos(j,0)<<"    "<<oldpos(j,1)<<"    "<< oldpos(j,2)<<endl;
          }
        }
	os.setf(ios::scientific);
        for(int j=0; j< D_array1(0)*D_array1(1)*D_array1(2); j++) {
          os << setw(16) << setprecision(8) << grid(i,j);
          if(j%6 ==5) os << endl;
        }
        os << endl;
	os.unsetf(ios::scientific);
	os<<setprecision(6);
        os.close();

        /*
          old plot to plt file version
        // http://www.csc.fi/gopenmol/developers/plt_format.phtml
        os<<"3 "; //rank=3 always
        os<<"2\n"; //dummy variable => "Orbital/density surface"
        //number of grid points for x, y, & z direction
        os <<D_array1(2)<<" "<<D_array1(1)<<" "<<D_array1(0)<<endl;
        os <<minmax(4)<<" "<<minmax(5)<<" "<<minmax(2)<<" "<<minmax(3)
           <<" "<<minmax(0)<<" "<<minmax(1)<<endl;
        
        for(int j=0; j<(D_array1(0)*D_array1(1)*D_array1(2)); j++)
          os<<grid(i,j)<<endl;
        os.close();
        
        os.open(confile.c_str());
        cout<<"writing to "<<confile<<endl;
        Array1 <doublevar> pos(3);
        sysprop->getAtomicLabels(atomlabels);
        os<<atomlabels.size()+ mywalker->electronSize()-2 <<endl;
        os << endl;
        for(unsigned int j=0; j<atomlabels.size();j++){
          mywalker->getIonPos(j, pos);
          os<<atomlabels[j] <<" "<< pos(0)
            <<" "<<pos(1)
            <<" "<< pos(2)<<endl;
        }
        for(int j=0; j<mywalker->electronSize(); j++){
          if ((j!=plots(i)-1)&&(j!=plots(i+1)-1)){
            if (j<spin_up) os<<"Eu";
            else os<<"Ed";
            // mywalker->getElectronPos(j, pos);
            os<<" "<< oldpos(j,0)<<" "<<oldpos(j,1)
              <<" "<< oldpos(j,2)<<endl;
          }
        }
        os.close();    
        */
    
      } 
  }
  cout <<"End of Nodes Method"<<endl;  
}


/*!
Print information about private variables {plots,resolution,minmax}
*/
int Nodes_method::showinfo(ostream & os)
{
   sysprop->showinfo(os);
    os << endl << endl;
    os << "                               Wavefunction  "
       << endl << endl;
    wfdata->showinfo(os);
  os<<"#############Nodes_method#################\n";
  if (doublemove)
    os <<"Doing doublemoves with"<<endl;
  os<<"Plots= "<<plots(0);
  for(int i=1;i<plots.GetSize();i++)
    os<<", "<<plots(i);
  os<<endl;

  if (doublemove){
    os << "Translation vector(s): "<<endl;
    for(int i=0;i<dxyz.GetDim(0);i++)    
      os<<"( "<<dxyz(i,0)<<" , "<<dxyz(i,1)<<" , "<<dxyz(i,2)<<" )"<<endl;
    }
  

  os<<"resolution "<<resolution<<endl;
  os<<"xmin="<<minmax(0)<<" xmax="<<minmax(1)<<endl;
  os<<"ymin="<<minmax(2)<<" ymax="<<minmax(3)<<endl;
  os<<"zmin="<<minmax(4)<<" zmax="<<minmax(5)<<endl;
  os<<"done"<<endl;
  return 1;
}
