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
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include "dataVisF.h"
using namespace std;

/*!
<h2>USE: dataVis -h 20 -ooqmc inputfile</h2>
 <h3>Output Data Format Options:</h3>
  "-h 20" will cause dataVis to bin data from inputfile
  into 20 bins for a "histogram". Default # of bins = 20<br>
  "-o file.out" will send output data to file.out<br>
  "-sc 1" will set serial correlation index to 1<br>
  "-sc -1" will look at correlation between weights
   psi[1]/(psi[1]+psi[2]) and psi[2]/(psi[1]+psi[2])
   for two wavefunctions with GUIDETYPE SUM<br>
 <h3>Input Data Format Options:</h3>
  "-ooqmc li2.path" tells dataVis to read ooqmc formatted
  .path file data. Default is ooqmc format, and if there
  is no "format" option, then the last argument will be
  read as the input data file.<br>
  "-wf 1" tells the code that you want to plot data for
  the first wavefunction.<br>
 <h3>Assumptions</h3>
  expects consistent file format<br>
  .path file with 8 items per line<br>
*/
int main(int argc, char * argv[]){
	//VARIABLE DECLARATIONS AND DEFAULT INITIALIZATIONS
	int nbins=0; //number of bins for histogram
	int scvar=0; //serial correlation index
	int wf=1; //wave function you want data for {1,2,3,...}
	int nwalk=0; //total number of walkers {1,2,3,...}
	int nwf=0; //total number of wavefunctions {1,2,3,...}
	int startindex=0; //for parsing
	int skipindex=0; //for parsing
	double runningsum=0;
	double sercorr;
	string idformat="OOQMC"; //input data format
	string ifilename=argv[argc-1]; //input data file name
	string ofilename; //default output filename
	string temp;
	vector <double> datapoints;
	vector <double> weights;

	//READ COMMAND LINE OPTIONS
	if(argc==1) error("USE: dataVis -h 20 -wf 1 -sc 1 -o file.out -ooqmc file.path");

	if(readoption(argv, argc, "-h", temp)) nbins=atoi(temp.c_str());

	//user specifies type of input file and name
	if(readoption(argv, argc, "-ooqmc", ifilename)){
		idformat="OOQMC";
		if(ifilename==""){ error("must enter a filename after -ooqmc"); }
	}
	cout<<"Input type = "<<idformat<<"\nRead data from "<<ifilename<<endl;

	//get serial correlation index (default=1)
	if(readoption(argv, argc, "-sc", temp)){
		scvar=atoi(temp.c_str());
	}

	//find out which wavefunction user wants data for
	if(readoption(argv, argc, "-wf", temp)){
		wf=atoi(temp.c_str());
		if(wf<1) error("-wf number must be larger than 0");
		if(wf>special_getwfmax(ifilename))
			error("err: you asked for a wf that is not in the input file");
		if(idformat=="OOQMC") cout<<"looking at data for wavefunction # "<<wf<<endl;
	}

	//Generate Histogram with nbins # of bins
	if(nbins>0){
		if(readoption(argv, argc, "-o", ofilename)) ofilename+=ifilename+".histogram";
		else ofilename=ifilename+".histogram";

		//get datapoints from column 5 and weights from col 8
		startindex=wf-1;
		nwf=special_getwfmax(ifilename);
		skipindex=nwf-1;
		datapoints=parsefile(ifilename,idformat,startindex,skipindex,5);
		weights=parsefile(ifilename,idformat,startindex,skipindex,8);
		generateHistogram(datapoints,weights,nbins,ofilename);
		
		//print out histogram of the weights
		vector <double> nullweight;
		nullweight.reserve(weights.size());
		for(uint i=0; i< weights.size(); i++) nullweight.push_back(1);
		string weightfilename=ofilename+"w";
		generateHistogram(weights, nullweight, nbins, weightfilename);
	}

	//Calculate Serial Correlation with index scvar
	if(scvar>0){
		if(readoption(argv, argc, "-o", ofilename)) ofilename+=ifilename+".serCorr";
		else ofilename=ifilename+".serCorr";
		//wipe ofilename clean of previous data since serial_correlation() appends
		ofstream zout453;
		zout453.open(ofilename.c_str());
		zout453.close();

		startindex=wf-1;
		nwf=special_getwfmax(ifilename);
		nwalk=special_getwalkermax(ifilename);
		skipindex=nwf*nwalk-1;
		runningsum=0;
		for(int i=0; i<nwalk; i++){
			datapoints=parsefile(ifilename,idformat,startindex+i*nwf,skipindex,5);
			weights=parsefile(ifilename,idformat,startindex+i*nwf,skipindex,8);
			//verify that serial correlation is between -1 and 1
			sercorr=serial_correlation(datapoints,weights,scvar,ofilename);
			runningsum+=sercorr;
			cout<<"walker # "<<i+1<<" serial correlation = "
				<<sercorr<<endl;
		}
		runningsum/=nwalk;
		cout<<"average serial correlation = "<<runningsum<<endl;

	}

	//analyze correlation between weights psi[1]/(psi[1]+psi[2]) and psi[2]/(psi[1]+psi[2])
	if(scvar==-1){
		nwf=special_getwfmax(ifilename);
		if(nwf!=2) error("'-sc -1' requires a pathfile with exactly 2 wavefunctions");

		vector <double> data1,data2,data3,weight1,weight2,weight3;
		double psi1,psi2;

		if(readoption(argv, argc, "-o", ofilename)) ofilename+=ifilename+".wfWeight";
		else ofilename=ifilename+".wfWeight";

		nwalk=special_getwalkermax(ifilename);
		startindex=0; //start from first line of input file
		skipindex=0; //

		//datapoints=ln(wfsample) weights=sign(wfsample)
		datapoints=parsefile(ifilename,idformat,startindex,skipindex,7);
		weights=parsefile(ifilename,idformat,startindex,skipindex,6);
		for(uint i=0; i<weights.size(); i++){
			psi1=weights[i]*exp(datapoints[i]);
			//cout<<" weights["<<i<<"] "<<weights[i]<<" datapoints[i] "<<datapoints[i]<<endl;
			i++;
			psi2=weights[i]*exp(datapoints[i]);
			//cout<<" weights["<<i<<"] "<<weights[i]<<" datapoints[i] "<<datapoints[i]<<endl;

			data1.push_back(psi1/(psi1+psi2)); //wf1 sample
			data2.push_back(psi2/(psi1+psi2)); //wf2 sample
			weight1.push_back(1);
			weight2.push_back(1);

			data3.push_back(psi1/(psi1+psi2)); //wf1 and wf2 in same vector
			data3.push_back(psi2/(psi1+psi2));
			weight3.push_back(1);
			weight3.push_back(1);
		}

		cout<<"covariance of weights psi[1]/(psi[1]+psi[2]) and psi[2]/(psi[1]+psi[2])\n"
			<<" for all blocks, steps, and walkers = "
			<<covariance(data1,data2,weight1,weight2,ofilename)<<endl;

		if(nbins<=0) nbins=100;
		ofilename=ifilename+".wf1.hist";
		generateHistogram(data1,weight1,nbins,ofilename);
		ofilename=ifilename+".wf2.hist";
		generateHistogram(data2,weight2,nbins,ofilename);
		ofilename=ifilename+".wf0.hist";
		generateHistogram(data3,weight3,nbins,ofilename);
	}

	/*/INCOMPLETE: statistics on ooqmc.path walker/step/block energies
	if(nbins==-1){
		ifstream zin863;
		double dud;
		const int nblocks=10;
		int temp_block;
		int temp_wf;
		double temp_energy;
		double temp_weight;
		double data[nblocks];
		double weightsum[nblocks];

		for(int i=0; i<nblocks; i++) data[i]=0;
		for(int i=0; i<nblocks; i++) weightsum[i]=0;

		zin863.open(ifilename.c_str());
		while(zin863>>dud){
			temp_block=(int)dud;
			zin863>>dud; zin863>>dud; zin863>>dud;
			temp_wf=(int)dud;
			zin863>>dud;
			temp_energy=dud;
			zin863>>dud; zin863>>dud; zin863>>dud;
			temp_weight=dud;

			if(temp_wf==(wf-1)) {
				data[temp_block]+=temp_energy*temp_weight;
				weightsum[temp_block]+=temp_weight;
			}
		}

		for(int i=0; i<nblocks; i++) 
			cout<<"block "<<i<<" avg = "<<(data[i]/weightsum[i])<<endl;

		zin863.close();
	}*/

/*////////////////////////////////////////////////
// check that the set of integers from 0 to 5000 has serial correlation of ~1
vector <double> v1,v2;
for(double i=0; i<5000; i++){
	v1.push_back(i);
	v2.push_back(1);
}
cout<<"TEST: serial correlation = "
	<<serial_correlation(v1,v2,1,"testfile")<<endl;

////////////////////////////////////////////////*/

	cout<<endl;
	return 0;
}







