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

bool readoption(char* elem[], int & size, const string val1, string & val2){
	for(int i=0; i<size; i++){
		if(0==strcmp(val1.c_str(), elem[i])){
			val2=elem[i+1];
			return true;
		}
	}
	return false;
}

void error(char * str) {
	cout << str << endl << "Exiting \n";
	exit(1);
}

double wavg(vector <double> & v, vector <double> & w){
	double runningsum=0;
	double weightsum=0;
	for(uint i=0; i<v.size(); i++){
		runningsum+=v[i]*w[i];
		weightsum+=w[i];
	}
	return (runningsum/weightsum);
}

double wstdev(vector <double> & v, vector <double> & w){
	double m=wavg(v,w);
	double runningsum=0;
	double weightsum=0;
	for(uint i=0; i<v.size(); i++){
		runningsum+=(v[i]-m)*(v[i]-m)*w[i];
		weightsum+=w[i];
	}
	return (sqrt(runningsum/weightsum));
}

double sum(vector <double> & v){
	double s=0;
	for(uint i=0; i<v.size(); i++) s+=v[i];
	return s;
}

double getmin(vector <double> & v){
	double m=v[0];
	for(uint i=0; i<v.size(); i++) if(m>v[i]) m=v[i];
	return m;
}

double getmax(vector <double> & v){
	double m=v[0];
	for(uint i=0; i<v.size(); i++) if(m<v[i]) m=v[i];
	return m;
}

int special_getwfmax(string ifilename){
	int maxwf=-1;
	double dud=0;
	ifstream zin;

	zin.open(ifilename.c_str());
	if(!zin) error("Unable to open input file");

	while(zin>>dud){//ignore first item on line, also ensures not at eof
		zin>>dud;zin>>dud; //ignore first two items
		zin>>dud; //read in wf#
		if(dud<=maxwf){//wf number should increase or stay at 0
			zin.close();
			return (maxwf+1);
		}
		if(dud>maxwf) maxwf=(int)dud;
		zin>>dud;zin>>dud;zin>>dud;zin>>dud; //ignore energy and last 3 items
	}

	//should not get here
	zin.close();
	error("can't tell how many wf's there are");
	return 0;
}

int special_getwalkermax(string ifilename){
	int maxwf=-1;
	double dud=0;
	ifstream zin;

	zin.open(ifilename.c_str());
	if(!zin) error("Unable to open input file");

	while(zin>>dud){//ignore first item on line, also ensures not at eof
		zin>>dud; //ignore first item
		zin>>dud; //read in wf#
		if(dud<maxwf){//walker number should increase
			zin.close();
			return (maxwf+1);
		}
		if(dud>maxwf) maxwf=(int)dud;
		zin>>dud;zin>>dud;zin>>dud;zin>>dud;zin>>dud; //ignore last 5 items
	}

	//should not get here
	zin.close();
	error("can't tell how many walker's there are");
	return 0;
}

vector <double> parsefile(string ifilename, string idformat, int opt1, int opt2, int opt3){
	vector <double> data;
	double dud;
	string ofilename;
	ifstream zin;
	int num_items=8; //#items per line in input file

	if(idformat!="OOQMC") error("ERROR: OOQMC path data is the only supported format");
	if(opt3>num_items) error("within the code, you are asking parsefile() for a column# > num_items");

	//OPEN INPUT FILE FOR READING
	zin.open(ifilename.c_str());
	if(!zin) error("Unable to open input file");

	//skip to appropriate wf, opt1=0 => get data for first wf
	for(int i=0; i<opt1; i++) for(int j=0; j<num_items; j++) zin>>dud;

	//read in data points
	while(zin>>dud){//ignore first item on line, also ensures not at eof
		for(int i=1; i<opt3; i++) zin>>dud; //ignore preceding data items
		data.push_back(dud); //store energy in vector
		for(int i=0; i<(num_items-opt3); i++) zin>>dud; //ignore following data items
		//ignore opt2 # of lines between data points we want
		for(int i=0; i<opt2; i++) for(int j=0; j<num_items; j++) zin>>dud;
	}

	//CLOSE INPUT FILE UPON COMPLETION
	zin.close();

	return data;
}

vector <double> generateHistogram(vector <double> & data, vector <double> & weight, int nbins, string ofilename){
	cout<<"-----GENERATING HISTOGRAM\n";
	double myavg=wavg(data,weight);
	double mystdev=wstdev(data,weight);
	double min=getmin(data);
	double max=getmax(data);
	vector <double> bincount;
	double total_weight=0; //for normalizing bincounts
	double bin_width=0; //for normalizing
	int bin=0;
	int stdev_width=10;
	ofstream zout;

	//initialize nbin elements of bincount to 0
	for(int i=0; i<nbins; i++) bincount.push_back(0);

	//print info about data
	cout<<"number of data points = "<<data.size()<<endl;
	cout<<"min = "<<min<<" max = "<<max<<endl;
	cout<<"weighted average = "<<myavg<<endl;
	cout<<"standard deviation = "<<mystdev<<endl;

	//to exclude outliers, plot datapoints only within "stdev_width" stdev of avg
	cout<<"Generating histogram using data within\n "<<stdev_width
		<<" standard deviations of the sample average."<<endl;
	if((myavg-stdev_width*mystdev)>min) min=myavg-stdev_width*mystdev;
	if((myavg+stdev_width*mystdev)<max) max=myavg+stdev_width*mystdev;
	bin_width=(max-min)/nbins; //for normalizing

	//for each datapoint, increment the bin which it belongs in
	for(uint i=0; i<=data.size(); i++){
		if((min<=data[i]) && (data[i]<max)){
			bin=(int)(nbins*(data[i]-min)/(max-min));
			bincount[bin]+=weight[i]; //weighted increment
			total_weight+=weight[i]; //determine total weight
		}
		else cout<<"OUTLIER: energy="<<data[i]<<endl;
	}

	//OPEN OUTPUT FILE FOR WRITING
	zout.open(ofilename.c_str());
	for(int i=0; i<nbins; i++){
		//print in "x y" format for xmgr, where the x points
		// are bin-min values, & y points are normalized bin counts
		zout<<(i*(max-min)/(nbins)+min)<<" ";
		zout<<(bincount[i]/(total_weight*bin_width))<<endl;
		//now output bin-max value to make it look like a histogram
		zout<<((i+1)*(max-min)/(nbins)+min)<<" ";
		zout<<(bincount[i]/(total_weight*bin_width))<<endl;
	}

	//CLOSE OUTPUT FILE UPON COMPLETION
	zout.close();

	for(uint i=0; i<bincount.size(); i++) bincount[i]/=(total_weight*bin_width);
	return bincount;
}

double serial_correlation(vector <double> & X, vector <double> & Xw, int index, string ofilename){
	cout<<"-----CALCULATING SERIAL CORRELATION\n";
	vector <double> U;
	vector <double> Uw;
	vector <double> V;
	vector <double> Vw;
	double cov_of_UV=0;
	double var_of_X=wstdev(X,Xw);var_of_X*=var_of_X;

	if((index<1) || (index>=(int)X.size()))
		error("err: in function 'serial_correlation', the index must be within bounds of data vector size");

	//this step is probably not necessary if you want to be efficient
	for(uint i=0; i<(X.size()-index); i++){
		U.push_back(X[i]);
		Uw.push_back(Xw[i]);
		V.push_back(X[i+index]);
		Vw.push_back(Xw[i+index]);
	}

	cov_of_UV=covariance(U,V,Uw,Vw,ofilename);

	return (cov_of_UV/var_of_X);
}

double covariance(vector <double> & v1, vector <double> & v2, vector <double> & w1, vector <double> & w2, string ofilename){
	//assertions
	if(v1.size()!=v2.size())
		error("unable to calculate covariance of datasets with different size");
	cout<<"-----CALCULATING COVARIANCE\n";

	double cov=0;
	double mean_of_v1=wavg(v1,w1);
	double mean_of_v2=wavg(v2,w2);

	//OUTPUT DATA PAIRS TO FILE
	ofstream zout;
	zout.open(ofilename.c_str());
	for(uint i=0; i<v1.size(); i++)
		zout<<v1[i]<<" "<<v2[i]<<endl;
	zout<<endl;
	zout.close();

	//CALCULATE COVARIANCE
	for(uint i=0; i<v1.size(); i++)
		for(uint j=0; j<v2.size(); j++)
			cov+=(v1[i]-mean_of_v1)*(v2[j]-mean_of_v2);
	cov/=v1.size()*v2.size();

	return cov;
}



