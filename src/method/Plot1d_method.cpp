/*
 
Copyright (C) 2008 Jindrich Kolorenc

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

#include "Plot1d_method.h"
#include "qmc_io.h"
#include "System.h"
#include "Program_options.h"

/*!

*/
void Plot1D_method::read(vector <string> words,
                       unsigned int & pos,
                       Program_options & options)
{

  if(!readvalue(words, pos=0, ngrid, "NGRID")) 
    ngrid=100;

  if(!readvalue(words, pos=0, cutoff, "CUTOFF")) 
    error("Need CUTOFF in PLOT1D section.");

  if(!readvalue(words, pos=0, suffix, "SUFFIX"))
    suffix="plt1d";

}

/*!
  mywfdata has to be allocated already in generateVariables, otherwise
  we could not use it in showinfo.
*/
int Plot1D_method::generateVariables(Program_options & options) {

  if(options.twftext.size() < 1)
    error("Need TRIALFUNC section for PLOT1D.");
  allocate(options.systemtext[0], sysprop );
  allocate(options.twftext[0], sysprop, mywfdata);

  return 1;
}

/*!
 
*/
void Plot1D_method::run(Program_options & options, ostream & output) {

  mywfdata->generateWavefunction(wf);

  string plotfilename=options.runid+"."+suffix;
  ofstream plotfile(plotfilename.c_str());

  Array1 <doublevar> xdata;
  xdata.Resize(ngrid);
  for (int i=0; i<ngrid; i++) xdata(i)=i*cutoff/ngrid;
  
  vector <Array1 <doublevar> > data;
  vector <string> desc;
  
  wf->plot1DInternals(xdata,data,desc,"");

  plotfile << "# column description:" << endl;
  plotfile << "#   distance" << endl;
  for ( unsigned int i=0; i<desc.size(); i++) {
    plotfile << "#   " << desc[i] << endl;
  }
  plotfile << endl;

  for ( int j=0; j<ngrid; j++) {
    plotfile << xdata[j] << " ";
    for ( unsigned int i=0; i<data.size(); i++) {
      plotfile << data[i](j) << " ";
    }
    plotfile << endl;
  }

  output << desc.size() << " functions found and \"plotted\" to "
	 << plotfilename << " :" << endl;
  for ( unsigned int i=0; i<desc.size(); i++) {
    output << "  " << desc[i] << endl;
  }
  output << endl;

  plotfile.close();
  mywfdata->clearObserver();
  if (wf) delete wf;
  wf=NULL;

}


/*!

*/
int Plot1D_method::showinfo(ostream & os)
{
  
  os << "Wavefunction" << endl << endl;
  mywfdata->showinfo(os);
  os << endl << "Plot1D method" << endl;
  os << "  " << ngrid << " (NGRID) mesh points between 0 and "
     << cutoff << " (CUTOFF)" << endl << endl;
  return 1;
}
