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
#include <string>
#include <fstream>
#include <iomanip>	//setw()
#include <cstdio>	//sprintf()
#include <cmath>	//floor()
#include <ctime>
using namespace std;

///////////////////////////////////////////////////////////////////////
//function declarations//function declarations//function declarations//
///////////////////////////////////////////////////////////////////////
/*! Print an error message and exit the program */
void error(char * arg1);

////////////////////////////////////////////////////////////////////
//main//main//main//main//main//main//main//main//main//main//main//
////////////////////////////////////////////////////////////////////
/*!
<h2>USE: walkerReDist c2h6.hf.config 7 c2h6.hf.newconfig 10</h2>
 detail: arg1 - utility name (walkerReDist) <br>
 arg2 - base name for new configuration files<br>
 arg3 - number of old config files to redistribute walkers from<br>
 arg4 - base name for new configuration files<br>
 arg5 - number of new config files to redistribute walkers to<br>
 <h3>Assumptions</h3>
  if file begins with "Array { nconfigs" where nconfigs is the # of
   configurations in the file, then the rest of the file is formatted
   like an ooqmc .config file with number of ions/electrons in each 
   section given by the numbers in the file.<br>
*/
int main(int argc, char * argv[]){
	time_t starttime;
	time(&starttime);

	cout<<endl;
	if(argc!=5) error("USE: walkerReDist c2h6.hf.config 'numConfigFiles you have' c2h6.opt.config 'numConfigFiles you want'");

	string ifilebase = argv[1];	//base name for old configuration files
	string ifilename;		//whole name for old configuration files
	string ofilebase = argv[3];	//base name for new configuration files
	string ofilename;		//whole name for new configuration files
	int n1 = atoi(argv[2]);		//number of old config files to redistribute walkers from
	int n2 = atoi(argv[4]);		//number of new config files to redistribute walkers to
	ofstream zout;			//re-usable output stream
	ifstream zin;			//re-usable input stream
	string text;			//temporary string holder for zin
	double text_d;			//temporary double holder for zin
	int text_i;			//temporary double holder for zin
	int nconfigs=0;			//total number of configurations
	char strbuff[40];		//character buffer to hold the file# as a character
	int numleftover;		//number of files that get maxConfigsPerFile
	int maxConfigsPerFile;		//count for the new config files
	int minConfigsPerFile;		//count for the new config files
	string sampleType;		//Periodic_sample  OR  Molecular_sample
	int numIons;			//number of 3d ions per configuration
	int numElectrons;		//number of 3d electrons per configuration
	int numWordsPerConfig;		//number of words per configuration
	//int numWordsPerFile;		//number of words in oldconfigfile minus first 3 and last bracket
	int numConfigsInNewfile=0;	//count of the # of configurations written in a newconfigfile
	int newFileCount=0;		//the number of the new file to write to



	/////////////////////////////////////////////////////////////////////
	//check that all input config files exist and have correct formatting
	/////////////////////////////////////////////////////////////////////
	if (n1<1) error("WARNING: not a valid number of old configuration files");
	if (n2<1) error("WARNING: not a valid number of new configuration files");
	for(int i=0;i<n1;i++){
		sprintf(strbuff, "%d", i);
		ifilename = ifilebase;
		ifilename += "_";
		ifilename += strbuff;
		zin.open(ifilename.c_str());
		//cout<<"Checking "<<ifilename<<" for consistency\n";
		if(!zin){
			cout<<"WARNING: Unable to find "<<ifilename<<endl;
			error("");
		}

		zin>>text;
		if(0!=strcmp(text.c_str(),"Array")){
			cout<<"BAD FORMAT: first word != Array in "<<ifilename<<endl;
			error("");
		}

		zin>>text;
		if(0!=strcmp(text.c_str(),"{")){
			cout<<"BAD FORMAT: second word != { in "<<ifilename<<endl;
			error("");
		}

		zin>>text_i;
		nconfigs += text_i;
		cout<<"     # of configurations in "<<ifilename<<" = "<<text_i<<endl;

		zin>>text;
		sampleType=text;
		zin>>text; zin>>text;
		numIons=atoi(text.c_str());
		for(int j=0;j<(numIons*3);j++) zin>>text;
		zin>>text; zin>>text;
		numElectrons=atoi(text.c_str());
		numWordsPerConfig=5+3*(numIons+numElectrons);
		//cout<<"numIons="<<numIons<<" numElectrons="<<numElectrons<<" sampleType="<<sampleType<<endl;

		zin.close();
	}



	///////////////////////////////////////////////////
	//output information about what code is going to do
	///////////////////////////////////////////////////
	numleftover = (int)(nconfigs-(n2*floor((double)nconfigs/(double)n2)));
	maxConfigsPerFile = (int)ceil((double)nconfigs/(double)n2);
	minConfigsPerFile = (int)floor((double)nconfigs/(double)n2);
	cout<<nconfigs<<" configurations stored in "<<n1<<" old files with base name \n";
	cout<<"     "<<ifilebase<<endl<<endl;
	cout<<nconfigs<<" configurations will be redistributed to ";
	cout<<n2<<" new files with base name\n";
	cout<<"     "<<ofilebase<<endl;

	


	/////////////////////////////////////////////////////////////
	//re-distribute walkers from oldconfigfiles to newconfigfiles
	/////////////////////////////////////////////////////////////
	//open first zout ouput stream here (I know there will be at least one)
	newFileCount=0;
	sprintf(strbuff, "%d", newFileCount);
	ofilename = ofilebase;
	ofilename += "_";
	ofilename += strbuff;
	zout.open(ofilename.c_str());
	cout<<"     Writing "<<maxConfigsPerFile<<" configurations to "<<ofilename<<endl;
	zout<<"Array { "<< maxConfigsPerFile << endl;
	//for each old config file
	for(int i=0;i<n1;i++){
		sprintf(strbuff, "%d", i);
		ifilename = ifilebase;
		ifilename += "_";
		ifilename += strbuff;
		zin.open(ifilename.c_str());
		zin>>text;zin>>text;zin>>text_i; //next zin will get 1st element of config

		nconfigs=text_i;
		//for each config in old config file
		for(int j=0; j<nconfigs; j++){
			//is it time to prepare a new zout output stream?
			if((numConfigsInNewfile >= minConfigsPerFile)){ //if I have min#configs in new file
				// open new zout stream if 'writing to minConfigFile' or if I have max#configs in new file
				if( (newFileCount>=numleftover) || (numConfigsInNewfile >= maxConfigsPerFile)){
					//zout '}', close zout, open zout on next output file, put first 3 words in
					zout<<" } "<<endl;
					zout.close();
					newFileCount++;		//just closed one more new file
					sprintf(strbuff, "%d", newFileCount);
					ofilename = ofilebase;
					ofilename += "_";
					ofilename += strbuff;
					zout.open(ofilename.c_str());
					numConfigsInNewfile=0;	//there are no configs in this new file
					if(newFileCount<numleftover){
						zout<<"Array { "<<maxConfigsPerFile<< endl;
						cout<<"     Writing "<<maxConfigsPerFile;
						cout<<" configurations to "<<ofilename<<endl;
					}
					else {
						zout<<"Array { "<<minConfigsPerFile<< endl;
						cout<<"     Writing "<<minConfigsPerFile;
						cout<<" configurations to "<<ofilename<<endl;
					}
				}
			}

			//transfer next config to new config file
			zin>>text;
			zout<<text<<endl;	//sample type
			zin>>text;
			zout<<text<<"  ";	//'numIons'
			zin>>text;
			zout<<text<<endl;	//numIons
			for(int k=0; k<numIons; k++){
				for(int m=0; m< 3; m++)
				{
					zin >> text_d;
					zout << setw(20) << setprecision(12) << text_d;
				}
				zout<<endl;
			}
			zin>>text;		//'numElectrons'
			zout<<"   "<<text<<"   ";
			zin>>text;		//numElectrons
			zout<<text<<endl;
			for(int k=0; k<numElectrons; k++){
				for(int m=0; m< 3; m++)
				{
					zin >> text_d;
					zout << setw(20) << setprecision(12) << text_d;
				}
				zout<<endl;
			}


			numConfigsInNewfile++;	//increment this counter
		}

		zin.close();
	}

	//if the last file was not closed, go ahead and close it properly
	zout<<" } "<<endl;	//if its already been closed, this wont get written to file
	zout.close();		//cant hurt to close stream if its already been closed
	newFileCount++;		//just closed one more new file

	cout<<endl<<"NOTE: after a VMC run each new config file will have ";
	cout<<maxConfigsPerFile;
	cout<<" configurations per file for a total of ";
	cout<<(n2*maxConfigsPerFile);
	cout<<" configurations (because VMC generates random walkers in any";
	cout<<" file that dosent have the number asked for in the input.\n";
	cout<<endl;
	
	time_t endtime;
	time(&endtime);
	cout<<"finished run. Total time: "<<difftime(endtime, starttime)<<" seconds. "<<endl<<endl;

}





////////////////////////////////////////////////////////////////////
//function definitions//function definitions//function definitions//
////////////////////////////////////////////////////////////////////
void error(char * str) {
	cout << str << endl << "Exiting \n";
	exit(1);
}







