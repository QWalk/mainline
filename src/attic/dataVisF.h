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
#ifndef DATAVISF_H
#define DATAVISF_H

#include <string>
#include <vector>

/*! for reading command line options. <br>
arg1 is the array of std::strings to look through<br>
arg2 is the number of std::strings to look through<br>
arg3 is the option tag like -h or -ooqmc<br>
arg4 is the <i>return path</i> for the std::string after the option tag<br>
function returns 'true' if it successfully reads an option value */
bool readoption(char* arg1[], int & arg2, const std::string arg3, std::string & arg4);

/*! Print an error message and exit the program */
void error(char * arg1);

/*! returns the weighted average of doubles stored in a std::vector<br>
arg1 is data<br>
arg2 is weights*/
double wavg(std::vector <double> & arg1, std::vector <double> & arg2);

/*! returns the weighted standard deviation of doubles stored in a std::vector<br>
arg1 is data<br>
arg2 is weights*/
double wstdev(std::vector <double> & arg1, std::vector <double> & arg2);

/*! returns the sum of items in the std::vector*/
double sum(std::vector <double> & arg1);

/*! returns the maximum value stored in the std::vector */
double getmax(std::vector <double> & arg1);

/*! returns the minimum value stored in the std::vector */
double getmin(std::vector <double> & arg1);

/*! a quick routine that opens the OOQMC .path input file
 and determines what the max number of wavefunctions is.<br>
arg1 is the input file name<br>
returns the max number of wavefunctions {1,2,3,...}*/
int special_getwfmax(std::string arg1);

/*! opens the OOQMC .path input file and determines how 
many walkers there are.<br>
arg1 is the input file name<br>
returns the number of walkers {1,2,3,...}*/
int special_getwalkermax(std::string arg1);

/*! Opens an input stream from the input file, extracts
 appropriate data for the given file type and options.
arg1 is the input data file name<br>
arg2 is the input data file format {OOQMC, etc...}<br>
arg3 is an option for which line to start on {0,1,2,...}<br>
arg4 is an option for how many lines to skip between
arg5 is the column number of the data item to extract
 each data point.<br>
returns std::vector containing data*/
std::vector <double> parsefile(std::string arg1, std::string arg2, int arg3, int arg4, int arg5);

/*! This function takes the data in a std::vector and outputs
data in a format ready for xmgr to plot. Data points
outside of 4 standard deviations of the sample average
are not included in the histogram data, but instead
the data point is output to the screen with a note.<br>
arg1 is a std::vector of data<br>
arg2 is a std::vector of weights<br>
arg3 is the number of bins<br>
arg4 is the output file name<br>
returns the std::vector of normalized bin counts*/
std::vector <double> generateHistogram(std::vector <double> & arg1, std::vector <double> & arg2, int arg3, std::string arg4);

/*! A measure of the correlation between U[i] and U[i+j]
 where j in the index of the j'th data element that comes
 after the i'th element. the equation reads: r = Cov(U,V)/var(X)
 where X is the set of data, and U and V are overlapping
 subsequences of X with U={x[0],x[1],...,x[N-j-1]} and
 V={x[j],x[1+j],...,x[N-1]<br>
arg1 is a std::vector containing the data<br>
arg2 is a std::vector containing the weights<br>
arg3 is the index 'j'<br>
arg4 is the output file name<br>
returns the serial correlation which ranges from -1 to 1
 with 0 indicating independant data points*/
double serial_correlation(std::vector <double> & arg1, std::vector <double> & arg2, int arg3, std::string arg4);

/*! calculate and return cov(arg1,arg2) assuming the joint probability mass function
is given by p(arg1,arg2)=1/(arg1.size()*arg2.size()) (i.e. equally weighted data points)<br>
quits with error message if arg1.size() != arg2.size()<br>
arg1 is a std::vector of data1<br>
arg2 is a std::vector of data2<br>
arg3 is a std::vector of weights1<br>
arg4 is a std::vector of weights2<br>
arg5 is the output file name for ordered pairs (arg1[i],arg2[i])
*/
double covariance(std::vector <double> & arg1, std::vector <double> & arg2, std::vector <double> & arg3, std::vector <double> & arg4, std::string arg5);


typedef unsigned int uint;

#endif //DATAVISF_H
