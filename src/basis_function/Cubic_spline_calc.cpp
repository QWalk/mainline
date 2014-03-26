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

#include "Qmc_std.h"
#include "Cubic_spline.h"


void Cubic_spline::calcVal(const Array1 <doublevar> & r,
                           Array1 <doublevar> & symvals,
                           const int startfill)
{

  if(r(0) >= threshold)
  {
    int end=startfill+nfunctions;
    for(int i=startfill; i< end; i++)
    {
      symvals(i)=0;
    }
  }
  else {

    doublevar v, v1, v2, v3, v4, v5;
    doublevar temp;

    int totf=startfill;
    for(int i=0; i<nsplines; i++)  {

      int interval=splines(i).getInterval(r(0));

      //temp=coeff(i,interval,0)+height*(coeff(i,interval,1)
      //                                 +height*(coeff(i,interval,2)
      //                                          +height*coeff(i,interval,3)));
      temp=splines(i).getVal(r(0), interval);
      //cout << "phir " << temp << endl;
      int end=totf+nfuncspline(i);
      for(; totf <end; totf++) {
        switch(indiv_symmetry(totf-startfill)) {
        case isym_S: // s
          symvals(totf)=temp;
          break;

        case isym_Px: //p
          symvals(totf)=temp*r(2);
          break;
        case isym_Py:
          symvals(totf)=temp*r(3);
          break;
        case isym_Pz:
          symvals(totf)=temp*r(4);
          break;

        case isym_Dxx: //d
          v=r(2)*r(2);  //xx
          symvals(totf)=v*temp;
          break;
        case isym_Dyy:
          v=r(3)*r(3); //yy
          symvals(totf)=v*temp;
          break;
        case isym_Dzz:
          v=r(4)*r(4); //zz
          symvals(totf)=v*temp;
          break;
        case isym_Dxy:
          v=r(2)*r(3); //xy
          symvals(totf)=v*temp;
          break;
        case isym_Dxz:
          v=r(2)*r(4);  //xz
          symvals(totf)=v*temp;
          break;
        case isym_Dyz:
          v=r(3)*r(4);  //yz
          symvals(totf)=v*temp;
          break;

        case isym_Dx2y2:
          v=r(2)*r(2)-r(3)*r(3);
          symvals(totf)=v*temp;
          break;

        case isym_Dz2r2:
          v=2.*r(4)*r(4)-r(2)*r(2)-r(3)*r(3);
          symvals(totf)=v*temp;
          break;


        case isym_Fxxx: //f
          //xxx
          v2=r(2)*r(2);
          v3=v2*r(2);
          symvals(totf)=v3*temp;
          break;

        case isym_Fyyy: //yyy
          v2=r(3)*r(3);
          v3=v2*r(3);
          symvals(totf)=v3*temp;
          break;

        case isym_Fzzz://zzz
          v2=r(4)*r(4);
          v3=v2*r(4);
          symvals(totf)=v3*temp;
          break;

        case isym_Fxxy://xxy
          v2=r(2)*r(2);
          v3=v2*r(3);
          symvals(totf)=v3*temp;
          break;

        case isym_Fxxz://xxz
          v2=r(2)*r(2);
          v3=v2*r(4);
          symvals(totf)=v3*temp;
          break;

        case isym_Fyyx://yyx
          v2=r(3)*r(3);
          v3=v2*r(2);
          symvals(totf)=v3*temp;
          break;

        case isym_Fyyz://yyz
          v2=r(3)*r(3);
          v3=v2*r(4);
          symvals(totf)=v3*temp;
          break;

        case isym_Fzzx://zzx
          v2=r(4)*r(4);
          v3=r(2)*v2;
          symvals(totf)=v3*temp;
          break;

        case isym_Fzzy://zzy
          v2=r(4)*r(4);
          v3=v2*r(3);
          symvals(totf)=v3*temp;
          break;

        case isym_Fxyz://xyz     
          v3=r(2)*r(3)*r(4);
          symvals(totf)=v3*temp;
          break;
          
        case isym_Fm3: //y(3*x*x-y*y)
          v3=r(3)*(3.0*r(2)*r(2)-r(3)*r(3));
          symvals(totf)=v3*temp;
          break;
          
        case isym_Fm1://(4z2 - x2 - y2)*y
          v3=r(3)*(4.0*r(4)*r(4)-r(2)*r(2)-r(3)*r(3));
          symvals(totf)=v3*temp;
          break;
        
        case isym_F0:
          v3=r(4)*(2.0*r(4)*r(4)-3.0*(r(2)*r(2)+r(3)*r(3)));
          symvals(totf)=v3*temp;
          break;
          
        case isym_Fp1://(4z2 - x2 - y2)x
          v3=r(2)*(4.0*r(4)*r(4)-r(2)*r(2)-r(3)*r(3));
          symvals(totf)=v3*temp;
          break;
          
        case isym_Fp2:
          v3=r(4)*(r(2)*r(2)-r(3)*r(3));
          symvals(totf)=v3*temp;
          break;
        
        case isym_Fp3:
          v3=r(2)*(r(2)*r(2)+r(3)*r(3));
          symvals(totf)=v3*temp;
          break;

        case isym_Fp3mod:
          v3=r(2)*(r(2)*r(2)-3.0*r(3)*r(3));
          symvals(totf)=v3*temp;
          break;  

                  
        case isym_Gxxxx: //g xxxx (15)
          v4=r(2)*r(2)*r(2)*r(2);
          symvals(totf)=v4*temp;
          break;

        case isym_Gyyyy: //yyyy
          v4=r(3)*r(3)*r(3)*r(3);
          symvals(totf)=v4*temp;
          break;

        case isym_Gzzzz://zzzz
          v4=r(4)*r(4)*r(4)*r(4);
          symvals(totf)=v4*temp;
          break;

        case isym_Gxxxy://xxxy
          v4=r(2)*r(2)*r(2)*r(3);
          symvals(totf)=v4*temp;
          break;

        case isym_Gxxxz://xxxz
          v4=r(2)*r(2)*r(2)*r(4);
          symvals(totf)=v4*temp;
          break;

        case isym_Gyyyx: //g yyyx 
          v4=r(3)*r(3)*r(3)*r(2);
          symvals(totf)=v4*temp;
          break;

        case isym_Gyyyz: //yyyz
          v4=r(3)*r(3)*r(3)*r(4);
          symvals(totf)=v4*temp;
          break;

        case isym_Gzzzx://zzzx
          v4=r(4)*r(4)*r(4)*r(2);
          symvals(totf)=v4*temp;
          break;

        case isym_Gzzzy://zzzy
          v4=r(4)*r(4)*r(4)*r(3);
          symvals(totf)=v4*temp;
          break;

        case isym_Gxxyy://xxyy
          v4=r(2)*r(2)*r(3)*r(3);
          symvals(totf)=v4*temp;
          break;

        case isym_Gxxzz: //g xxzz 
          v4=r(2)*r(2)*r(4)*r(4);
          symvals(totf)=v4*temp;
          break;

        case isym_Gyyzz: //yyzz
          v4=r(3)*r(3)*r(4)*r(4);
          symvals(totf)=v4*temp;
          break;

        case isym_Gxxyz://xxyz
          v4=r(2)*r(2)*r(3)*r(4);
          symvals(totf)=v4*temp;
          break;

        case isym_Gyyxz://yyxz
          v4=r(3)*r(3)*r(2)*r(4);
          symvals(totf)=v4*temp;
          break;

        case isym_Gzzxy://zzxy  end of g
          v4=r(4)*r(4)*r(2)*r(3);
          symvals(totf)=v4*temp;
          break;


	case isym_G0: //g0
	  v1=r(2)*r(2);
	  v2=r(3)*r(3);
	  v3=r(4)*r(4);
	  v5=v1+v2+v3;
          v4=35.*v3*v3-30.*v3*v5+3.*v5*v5;
          symvals(totf)=v4*temp;
          break;
	case isym_G1: //g1
	  v1=r(2)*r(2);
	  v2=r(3)*r(3);
	  v3=r(4)*r(4);
	  v5=v1+v2+v3;
          v4=r(2)*r(4)*(7.*v3-3.*v5);
          symvals(totf)=v4*temp;
          break;
	case isym_G2: //g2
	  v1=r(2)*r(2);
	  v2=r(3)*r(3);
	  v3=r(4)*r(4);
	  v5=v1+v2+v3;
          v4=r(3)*r(4)*(7.*v3-3.*v5);
          symvals(totf)=v4*temp;
          break;
	case isym_G3: //g3
	  v1=r(2)*r(2);
	  v2=r(3)*r(3);
	  v3=r(4)*r(4);
	  v5=v1+v2+v3;
          v4=(v1-v2)*(7.*v3-v5);
          symvals(totf)=v4*temp;
          break;
        case isym_G4: //g4
	  v1=r(2)*r(2);
	  v2=r(3)*r(3);
	  v3=r(4)*r(4);
	  v5=v1+v2+v3;
          v4=r(2)*r(3)*(7.*v3-v5);
          symvals(totf)=v4*temp;
          break;
	case isym_G5: //g5
	  v1=r(2)*r(2);
	  v2=r(3)*r(3);
          v4=(v1-3.*v2)*r(2)*r(4);
          symvals(totf)=v4*temp;
          break;  
	case isym_G6: //g6
	  v1=r(2)*r(2);
	  v2=r(3)*r(3);
          v4=(3.*v1-v2)*r(3)*r(4);
          symvals(totf)=v4*temp;
          break; 
	case isym_G7: //g7
	  v1=r(2)*r(2);
	  v2=r(3)*r(3);
          v4=(v1-3.*v2)*v1-(3.*v1-v2)*v2;
          symvals(totf)=v4*temp;
          break;
	case isym_G8: //g8
	  v1=r(2)*r(2);
	  v2=r(3)*r(3);
          v4=(v1-v2)*r(2)*r(3);
          symvals(totf)=v4*temp;
          break;

        default:
          error("Bad symmetry in Cubic_spline::calcVal", indiv_symmetry(totf-startfill));
        }

      }
    }
  }

}

void Cubic_spline::calcLap(
  const Array1 <doublevar> & r,
  //!< in form r, r^2, x, y, z
  Array2 <doublevar> & symvals,
  //!< The values of the spline propogated through symmetry.  For example, a p state would be a 3x3 matrix, and an s state a 1x3.
  const int startfill
)
{

  assert(r.GetDim(0) >= 5);
  //cout << "spline interval " << interval << "   r   " << r(0) << endl;
  if(r(0) >= threshold)
  {
    int end=startfill+nfunctions;
    for(int i=startfill; i< end; i++)
    {
      for(int d=0; d< 5; d++)
      {
        symvals(i,d)=0;
      }
    }
  }
  else
  {
    //cout << "calculating spline\n";

    //int interval=splines(0).getInterval(r(0));

    doublevar rdp, v, h, x, y, z, v1, v2, v3, v4, v5, u2, v2u2, dfdx, dfdy, dfdz;
    doublevar func, fdir, f2dir;
    //doublevar funct[3];
    int totf=startfill;

    assert(r(0) < threshold);
    assert(symvals.GetDim(0) >= nfunctions);


    for(int i=0; i<nsplines; i++) {

      int interval=splines(i).getInterval(r(0));
      //cout << "interval " << interval << " r " << r(0) << endl;
      //func=coeff(i,interval,0)+height*(coeff(i,interval,1)
      //                                 +height*(coeff(i,interval,2)
      //                                          +height*coeff(i,interval,3)));

      //fdir=(coeff(i,interval,1)+height*(2*coeff(i,interval,2)
      //                                  +height*(3*coeff(i,interval,3))))/r(0);

      //f2dir=2*coeff(i,interval,2)+6*height*coeff(i,interval,3);
      splines(i).getDers(r(0), interval,func, fdir, f2dir);
      
      //cout << "r(0) " << r(0) << " interval " << interval << " func " << func << endl;

      //cout << "height " << height << endl;
      //for(int j=0; j<4; j++) {
      //  cout << "coeff " << j << "  " << coeff(i,interval,j) << endl;
      //}
      //cout << "rad fn " << func << "   " << fdir << "   "
      //    << f2dir << endl;

      int end=totf+nfuncspline(i);
      for(; totf <end; totf++)
      {
        switch(indiv_symmetry(totf-startfill))
        {
        case isym_S: // s
          symvals(totf,0)=func;
          symvals(totf,1)=fdir*r(2);
          symvals(totf,2)=fdir*r(3);
          symvals(totf,3)=fdir*r(4);
          symvals(totf,4)=f2dir+2.*fdir;
          break;

        case isym_Px: //p
          symvals(totf,0)=func*r(2);
          rdp=fdir*r(2);
          symvals(totf,1)=rdp*r(2)+func;
          symvals(totf ,2)=rdp*r(3);
          symvals(totf ,3)=rdp*r(4);
          symvals(totf ,4)=f2dir*r(2)+4.*rdp;
          break;

        case isym_Py:
          symvals(totf ,0)=func*r(3);
          rdp=fdir*r(3);
          symvals(totf ,1)=rdp*r(2);
          symvals(totf ,2)=rdp*r(3)+func;
          symvals(totf ,3)=rdp*r(4);
          symvals(totf ,4)=f2dir*r(3)+4.*rdp;
          break;

        case isym_Pz:
          symvals(totf ,0)=func*r(4);
          rdp=fdir*r(4);
          symvals(totf ,1)=rdp*r(2);
          symvals(totf ,2)=rdp*r(3);
          symvals(totf ,3)=rdp*r(4)+func;
          symvals(totf ,4)=f2dir*r(4)+4.*rdp;
          break;

        case isym_Dxx: //d
          v=r(2)*r(2);  //xx
          symvals(totf,0)=v*func;
          h=v*fdir;
          symvals(totf,1)=h*r(2)+func*2.*r(2);
          symvals(totf,2)=h*r(3);
          symvals(totf,3)=h*r(4);
          symvals(totf,4)=v*f2dir+6.*h+2.*func;
          break;

        case isym_Dyy:
          v=r(3)*r(3); //yy
          symvals(totf,0)=v*func;
          h=v*fdir;
          symvals(totf,1)=h*r(2);
          symvals(totf,2)=h*r(3)+func*2.*r(3);
          symvals(totf,3)=h*r(4);
          symvals(totf,4)=v*f2dir+6.*h+2.*func;
          break;

        case isym_Dzz:
          v=r(4)*r(4); //zz
          symvals(totf,0)=v*func;
          h=v*fdir;
          symvals(totf,1)=h*r(2);
          symvals(totf,2)=h*r(3);
          symvals(totf,3)=h*r(4)+func*2.*r(4);
          symvals(totf,4)=v*f2dir+6.*h+2.*func;
          break;

        case isym_Dxy:
          v=r(2)*r(3); //xy
          symvals(totf,0)=v*func;
          h=v*fdir;
          symvals(totf,1)=h*r(2)+func*r(3);
          symvals(totf,2)=h*r(3)+func*r(2);
          symvals(totf,3)=h*r(4);
          symvals(totf,4)=v*f2dir+6.*h;
          break;

        case isym_Dxz:
          v=r(2)*r(4);  //xz
          symvals(totf,0)=v*func;
          h=v*fdir;
          symvals(totf,1)=h*r(2)+func*r(4);
          symvals(totf,2)=h*r(3);
          symvals(totf,3)=h*r(4)+func*r(2);
          symvals(totf,4)=v*f2dir+6.*h;
          break;

        case isym_Dyz:
          v=r(3)*r(4);  //yz
          symvals(totf,0)=v*func;
          h=v*fdir;
          symvals(totf,1)=h*r(2);
          symvals(totf,2)=h*r(3)+func*r(4);
          symvals(totf,3)=h*r(4)+func*r(3);
          symvals(totf,4)=v*f2dir+6.*h;
          break;

        case isym_Dx2y2:
          v=r(2)*r(2)-r(3)*r(3);
          symvals(totf,0)=v*func;
          h=v*fdir;
          symvals(totf,1)=h*r(2)+2.*func*r(2);
          symvals(totf,2)=h*r(3)-2.*func*r(3);
          symvals(totf,3)=h*r(4);
          symvals(totf,4)=v*f2dir+6.*h;
          break;

        case isym_Dz2r2:
          v=2.*r(4)*r(4)-r(2)*r(2)-r(3)*r(3);
          symvals(totf,0)=v*func;
          h=v*fdir;
          symvals(totf,1)=h*r(2)-2*func*r(2);
          symvals(totf,2)=h*r(3)-2*func*r(3);
          symvals(totf,3)=h*r(4)+4*func*r(4);
          symvals(totf,4)=v*f2dir+6*h;
          break;

          //f
        case isym_Fxxx: //xxx
          v2=r(2)*r(2);
          v3=v2*r(2);
          symvals(totf,0)=v3*func;
          h=v3*fdir;
          symvals(totf, 1)=h*r(2)+func*3.*v2;
          symvals(totf, 2)=h*r(3);
          symvals(totf, 3)=h*r(4);
          symvals(totf, 4)=v3*f2dir+8.*h+6.*func*r(2);
          break;

        case isym_Fyyy://yyy
          v2=r(3)*r(3);
          v3=v2*r(3);
          symvals(totf,0)=v3*func;
          h=v3*fdir;
          symvals(totf, 1)=h*r(2);
          symvals(totf, 2)=h*r(3)+3.*func*v2;
          symvals(totf, 3)=h*r(4);
          symvals(totf, 4)=v3*f2dir+8.*h+6.*func*r(3);
          break;


        case isym_Fzzz: //zzz
          v2=r(4)*r(4);
          v3=v2*r(4);
          symvals(totf, 0)=v3*func;
          h=v3*fdir;
          symvals(totf, 1)=h*r(2);
          symvals(totf, 2)=h*r(3);
          symvals(totf, 3)=h*r(4)+3.*func*v2;
          symvals(totf, 4)=v3*f2dir+8.*h+6.*func*r(4);
          break;

        case isym_Fxxy: //xxy
          v2=r(2)*r(2);
          v3=v2*r(3);
          symvals(totf, 0)=v3*func;
          h=v3*fdir;
          symvals(totf, 1)=h*r(2)+2.*func*r(2)*r(3);
          symvals(totf, 2)=h*r(3)+func*v2;
          symvals(totf, 3)=h*r(4);
          symvals(totf, 4)=v3*f2dir+8.*h+2.*func*r(3);
          break;

        case isym_Fxxz://xxz
          v2=r(2)*r(2);
          v3=v2*r(4);
          symvals(totf, 0)=v3*func;
          h=v3*fdir;
          symvals(totf, 1)=h*r(2)+2.*func*r(2)*r(4);
          symvals(totf, 2)=h*r(3);
          symvals(totf, 3)=h*r(4)+func*v2;
          symvals(totf, 4)=v3*f2dir+8.*h+2.*func*r(4);
          break;

        case isym_Fyyx://yyx
          v2=r(3)*r(3);
          v3=v2*r(2);
          symvals(totf, 0)=v3*func;
          h=v3*fdir;
          symvals(totf, 1)=h*r(2)+v2*func;
          symvals(totf, 2)=h*r(3)+2.*func*r(2)*r(3);
          symvals(totf, 3)=h*r(4);
          symvals(totf, 4)=v3*f2dir+8.*h+2.*func*r(2);
          break;

        case isym_Fyyz://yyz
          v2=r(3)*r(3);
          v3=v2*r(4);
          symvals(totf, 0)=v3*func;
          h=v3*fdir;
          symvals(totf, 1)=h*r(2);
          symvals(totf, 2)=h*r(3)+2.*func*r(3)*r(4);
          symvals(totf, 3)=h*r(4)+func*v2;
          symvals(totf, 4)=v3*f2dir+8.*h+2.*func*r(4);
          break;

        case isym_Fzzx://zzx
          v2=r(4)*r(4);
          v3=r(2)*v2;
          symvals(totf, 0)=v3*func;
          h=v3*fdir;
          symvals(totf, 1)=h*r(2)+func*v2;
          symvals(totf, 2)=h*r(3);
          symvals(totf, 3)=h*r(4)+2.*func*r(2)*r(4);
          symvals(totf, 4)=v3*f2dir+8.*h+2.*func*r(2);
          break;

        case isym_Fzzy://zzy
          v2=r(4)*r(4);
          v3=v2*r(3);
          symvals(totf, 0)=v3*func;
          h=v3*fdir;
          symvals(totf, 1)=h*r(2);
          symvals(totf, 2)=h*r(3)+func*v2;
          symvals(totf, 3)=h*r(4)+2.*func*r(3)*r(4);
          symvals(totf, 4)=v3*f2dir+8.*h+2.*func*r(3);
          break;

        case isym_Fxyz://xyz
          v3=r(2)*r(3)*r(4);
          symvals(totf, 0)=v3*func;
          v2=func*r(4);
          h=v3*fdir;
          symvals(totf, 1)=h*r(2)+v2*r(3);
          symvals(totf, 2)=h*r(3)+v2*r(2);
          symvals(totf, 3)=h*r(4)+func*r(2)*r(3);
          symvals(totf, 4)=v3*f2dir+8.*h;
          break;
          
          
        case isym_Fm3: 
          v3=r(3)*(3.0*r(2)*r(2)-r(3)*r(3));
          symvals(totf,0)=v3*func;
          h=v3*fdir;
          v2=6.0*r(2)*r(3);
          v4=3.0*(r(2)*r(2)-r(3)*r(3));
          symvals(totf,1)=h*r(2)+v2*func;
          symvals(totf,2)=h*r(3)+v4*func;
          symvals(totf,3)=h*r(4);
          symvals(totf,4)=v3*f2dir+8.*h;
          break;
          
        case isym_Fm1:
          v3=r(3)*(4.0*r(4)*r(4)-r(2)*r(2)-r(3)*r(3));
          symvals(totf,0)=v3*func;
          h=v3*fdir;
          symvals(totf,1)=h*r(2)-2.*r(2)*r(3)*func;
          symvals(totf,2)=h*r(3)+(4.*r(4)*r(4)-r(2)*r(2)-3.*r(3)*r(3))*func;
          symvals(totf,3)=h*r(4)+8.*r(3)*r(4)*func;
          symvals(totf,4)=v3*f2dir+8.*h;
          break;
          
        case isym_F0:
          v3=r(4)*(2.0*r(4)*r(4)-3.0*(r(2)*r(2)+r(3)*r(3)));
          symvals(totf,0)=v3*func;
          h=v3*fdir;
          symvals(totf,1)=h*r(2)-6.*r(2)*r(4)*func;
          symvals(totf,2)=h*r(3)-6.*r(3)*r(4)*func;
          symvals(totf,3)=h*r(4)+func*(6.*r(4)*r(4)-3.*(r(2)*r(2)+r(3)*r(3)));
          symvals(totf,4)=v3*f2dir+8.*h;
          break;
          
        case isym_Fp1:
          v3=r(2)*(4.0*r(4)*r(4)-r(2)*r(2)-r(3)*r(3));
          symvals(totf,0)=v3*func;
          h=v3*fdir;
          symvals(totf,1)=h*r(2)+func*(4.*r(4)*r(4)-3.*r(2)*r(2)-r(3)*r(3));
          symvals(totf,2)=h*r(3)-2.*r(2)*r(3)*func;
          symvals(totf,3)=h*r(4)+8*r(2)*r(4)*func;
          symvals(totf,4)=v3*f2dir+8.*h;

          break;
          
        case isym_Fp2:
          v3=r(4)*(r(2)*r(2)-r(3)*r(3));
          symvals(totf,0)=v3*func;
          h=v3*fdir;
          symvals(totf,1)=h*r(2)+2.*r(2)*r(4)*func;
          symvals(totf,2)=h*r(3)-2.*r(4)*r(3)*func;
          symvals(totf,3)=h*r(4)+func*(r(2)*r(2)-r(3)*r(3));
          symvals(totf,4)=v3*f2dir+8.*h;
          break;
          
        case isym_Fp3:
          v3=r(2)*(r(2)*r(2)+r(3)*r(3));
          symvals(totf,0)=v3*func;
          h=v3*fdir;
          symvals(totf,1)=h*r(2)+func*(3.*r(2)*r(2)+r(3)*r(3));
          symvals(totf,2)=h*r(3)+2.*r(2)*r(3)*func;
          symvals(totf,3)=h*r(4);
          symvals(totf,4)=v3*f2dir+8.*h+func*8.*r(2);
          break;
       case isym_Fp3mod:
          v3=r(2)*(r(2)*r(2)-3.*r(3)*r(3));
          symvals(totf,0)=v3*func;
          h=v3*fdir;
          symvals(totf,1)=h*r(2)+func*(3.*r(2)*r(2)-3.0*r(3)*r(3));
          symvals(totf,2)=h*r(3)-6.0*r(2)*r(3)*func;
          symvals(totf,3)=h*r(4);
          symvals(totf,4)=v3*f2dir+8.*h;
          break;   
          
          
          
        case isym_Gxxxx: // g xxxx done by MB
	  v2=r(2)*r(2);
	  v3=v2*r(2);
	  v4=v2*v2;
	  h=v4*fdir;
          symvals(totf,0)=v4*func;
          symvals(totf,1)=4.*v3*func+h*r(2);
          symvals(totf,2)=h*r(3);
          symvals(totf,3)=h*r(4);
          symvals(totf,4)=12.*v2*func+10.*h+v4*f2dir;
          break;

        case isym_Gyyyy: // g yyyy
	  v2=r(3)*r(3);
	  v3=v2*r(3);
	  v4=v2*v2;
	  h=v4*fdir;
          symvals(totf,0)=v4*func;
          symvals(totf,1)=h*r(2);
          symvals(totf,2)=4.*v3*func+h*r(3);
          symvals(totf,3)=h*r(4);
          symvals(totf,4)=12.*v2*func+10.*h+v4*f2dir;
	  break;

        case isym_Gzzzz: // g zzzz
	  v2=r(4)*r(4);
	  v3=v2*r(4);
	  v4=v2*v2;
	  h=v4*fdir;
          symvals(totf,0)=v4*func;
          symvals(totf,1)=h*r(2);
          symvals(totf,2)=h*r(3);
          symvals(totf,3)=4.*v3*func+h*r(4);
          symvals(totf,4)=12.*v2*func+10.*h+v4*f2dir;
          break;

        case isym_Gxxxy: // g xxxy
	  v2=r(2)*r(2);
	  v3=v2*r(2);
	  v4=v2*v2;
	  h=fdir*v3;	
	  symvals(totf,0)=v3*r(3)*func;
          symvals(totf,1)=r(3)*(3.*v2*func+v4*fdir);
          symvals(totf,2)=v3*func+r(3)*r(3)*h;
          symvals(totf,3)=h*r(3)*r(4);
          symvals(totf,4)=r(2)*r(3)*(6.*func+v2*(10.*fdir+f2dir));
	  break;

        case isym_Gxxxz: // g xxxzarXiv:cond-mat/9503150v2
	  v2=r(2)*r(2);
	  v3=v2*r(2);
	  v4=v2*v2;
	  h=fdir*v3;
	  symvals(totf,0)=v3*r(4)*func;
          symvals(totf,1)=r(4)*(3.*v2*func+v4*fdir);
          symvals(totf,2)=h*r(3)*r(4);
          symvals(totf,3)=v3*func+r(4)*r(4)*h;
          symvals(totf,4)=r(2)*r(4)*(6.*func+v2*(10.*fdir+f2dir));
          break;

        case isym_Gyyyx: // g yyyx
	  v2=r(3)*r(3);
	  v3=v2*r(3);
	  v4=v2*v2;
          h=fdir*v3;
	  symvals(totf,0)=v3*r(2)*func;
          symvals(totf,1)=v3*func+r(2)*r(2)*h;
          symvals(totf,2)=r(2)*(3.*v2*func+v4*fdir);
          symvals(totf,3)=h*r(2)*r(4);
          symvals(totf,4)=r(2)*r(3)*(6.*func+v2*(10.*fdir+f2dir));
          break;

        case isym_Gyyyz: // g yyyz
	  v2=r(3)*r(3);
	  v3=v2*r(3);
	  v4=v2*v2;
          h=fdir*v3;
	  symvals(totf,0)=v3*r(4)*func;
          symvals(totf,1)=h*r(2)*r(4);
          symvals(totf,2)=r(4)*(3.*v2*func+v4*fdir);
          symvals(totf,3)=v3*func+r(4)*r(4)*h;
          symvals(totf,4)=r(4)*r(3)*(6.*func+v2*(10.*fdir+f2dir));
          break;

        case isym_Gzzzx: // g zzzx
	  v2=r(4)*r(4);
	  v3=v2*r(4);
	  v4=v2*v2;
	  h=fdir*v3;
	  symvals(totf,0)=v3*r(2)*func;
          symvals(totf,1)=v3*func+r(2)*r(2)*h;
          symvals(totf,2)=h*r(2)*r(3);
          symvals(totf,3)=r(2)*(3.*v2*func+v4*fdir);
          symvals(totf,4)=r(2)*r(4)*(6.*func+v2*(10.*fdir+f2dir));
          break;

        case isym_Gzzzy: // g zzzy
	  v2=r(4)*r(4);
	  v3=v2*r(4);
	  v4=v2*v2;
          h=fdir*v3;
	  symvals(totf,0)=v3*r(3)*func;
          symvals(totf,1)=h*r(3)*r(2);
          symvals(totf,2)=v3*func+r(3)*r(3)*h;
          symvals(totf,3)=r(3)*(3.*v2*func+v4*fdir);
          symvals(totf,4)=r(3)*r(4)*(6.*func+v2*(10.*fdir+f2dir));
          break;

        case isym_Gxxyy: // g xxyy
	  v2=r(2)*r(2);
	  u2=r(3)*r(3);
	  v2u2=v2*u2;
	  symvals(totf,0)=v2u2*func;
	  symvals(totf,1)=r(2)*u2*(2.*func+v2*fdir);
          symvals(totf,2)=r(3)*v2*(2.*func+u2*fdir);
          symvals(totf,3)=v2u2*r(4)*fdir;
	  symvals(totf,4)=2.*(v2+u2)*func+v2u2*(10.*fdir+f2dir);
	  break;

        case isym_Gxxzz: // g xxzz
	  v2=r(2)*r(2);
	  u2=r(4)*r(4);
	  v2u2=v2*u2;
	  symvals(totf,0)=v2u2*func;
	  symvals(totf,1)=r(2)*u2*(2.*func+v2*fdir);
          symvals(totf,2)=v2u2*r(3)*fdir;
          symvals(totf,3)=r(4)*v2*(2.*func+u2*fdir);
	  symvals(totf,4)=2.*(v2+u2)*func+v2u2*(10.*fdir+f2dir);
          break;

        case isym_Gyyzz: // g yyzz
	  v2=r(3)*r(3);
	  u2=r(4)*r(4);
	  v2u2=v2*u2;
	  symvals(totf,0)=v2u2*func;
	  symvals(totf,1)=v2u2*r(2)*fdir;
          symvals(totf,2)=r(3)*u2*(2.*func+v2*fdir);
          symvals(totf,3)=r(4)*v2*(2.*func+u2*fdir);
	  symvals(totf,4)=2.*(v2+u2)*func+v2u2*(10.*fdir+f2dir);
          break;

        case isym_Gxxyz: // g xxyz
	  v2=r(2)*r(2);
	  u2=r(3)*r(4);
	  symvals(totf,0)=v2*u2*func;
	  symvals(totf,1)=r(2)*u2*(2.*func+v2*fdir);
          symvals(totf,2)=v2*r(4)*(func+r(3)*r(3)*fdir);
	  symvals(totf,3)=v2*r(3)*(func+r(4)*r(4)*fdir);
	  symvals(totf,4)=u2*(2.*func+v2*(10.*fdir+f2dir));
          break;

        case isym_Gyyxz: // g yyxz
	  v2=r(3)*r(3);
	  u2=r(2)*r(4);
	  symvals(totf,0)=v2*u2*func;
	  symvals(totf,1)=v2*r(4)*(func+r(2)*r(2)*fdir);
          symvals(totf,2)=r(3)*u2*(2.*func+v2*fdir);
	  symvals(totf,3)=v2*r(2)*(func+r(4)*r(4)*fdir);
	  symvals(totf,4)=u2*(2.*func+v2*(10.*fdir+f2dir));
          break;

        case isym_Gzzxy: // g zzxy
	  v2=r(4)*r(4);
	  u2=r(2)*r(3);
	  symvals(totf,0)=v2*u2*func;
	  symvals(totf,1)=v2*r(3)*(func+r(2)*r(2)*fdir);
          symvals(totf,2)=v2*r(2)*(func+r(3)*r(3)*fdir);
	  symvals(totf,3)=r(4)*u2*(2.*func+v2*fdir);
	  symvals(totf,4)=u2*(2.*func+v2*(10.*fdir+f2dir));
          break;

	case isym_G0: //g0
	  x=r(2);
	  y=r(3);
	  z=r(4);
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;
          v4=35.*v3*v3-30.*v3*v5+3.*v5*v5;
	  dfdx=  12.0*(v1 + v2 - 4*v3)*x;
	  dfdy=  12.0 *(v1 + v2 - 4*v3)*y;
	  dfdz= -48.0*(v1 + v2 + v3)*z + 80*v3*z;
	  h=fdir*v4;
          symvals(totf,0)=v4*func;
	  symvals(totf,1)=func*dfdx+h*x;
	  symvals(totf,2)=func*dfdy+h*y;
	  symvals(totf,3)=func*dfdz+h*z;
	  symvals(totf,4)=2*h+2*fdir*(dfdx*x+dfdy*y+dfdz*z)+v4*f2dir;
          break;

	case isym_G1: //g1
	  x=r(2);
	  y=r(3);
	  z=r(4);
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;
          v4=x*z*(7.*v3-3.*v5);
	  dfdx= (-9.*v1 - 3.*v2 + 4.*v3)*z;
	  dfdy= -6.*x*y*z;
	  dfdz= (-3.*v1 - 3.*v2 + 12.*v3)*x;
	  h=fdir*v4;
          symvals(totf,0)=v4*func;
	  symvals(totf,1)=func*dfdx+h*x;
	  symvals(totf,2)=func*dfdy+h*y;
	  symvals(totf,3)=func*dfdz+h*z;
	  symvals(totf,4)=2*h+2*fdir*(dfdx*x+dfdy*y+dfdz*z)+v4*f2dir;
	  break;

	case isym_G2: //g2
	  x=r(2);
	  y=r(3);
	  z=r(4);
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;
          v4=y*z*(7.*v3-3.*v5);
	  dfdx= -6.*x*y*z;
	  dfdy= (-3.*v1 - 9.*v2 + 4.*v3)*z;
	  dfdz= (-3.*v1 - 3.*v2 + 12.*v3)*y;
	  h=fdir*v4;
          symvals(totf,0)=v4*func;
	  symvals(totf,1)=func*dfdx+h*x;
	  symvals(totf,2)=func*dfdy+h*y;
	  symvals(totf,3)=func*dfdz+h*z;
	  symvals(totf,4)=2*h+2*fdir*(dfdx*x+dfdy*y+dfdz*z)+v4*f2dir;
	  break;

	case isym_G3: //g3
	  x=r(2);
	  y=r(3);
	  z=r(4);
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;
          v4=(v1 - v2)*(-v5 + 7. *v3);
	  dfdx= -4.*v1*x + 12.*v3*x;
	  dfdy= 4.*v2*y - 12.*v3*y;
	  dfdz= 12.*(v1 - 1.*v2)*z;
	  h=fdir*v4;
          symvals(totf,0)=v4*func;
	  symvals(totf,1)=func*dfdx+h*x;
	  symvals(totf,2)=func*dfdy+h*y;
	  symvals(totf,3)=func*dfdz+h*z;
	  symvals(totf,4)=2*h+2*fdir*(dfdx*x+dfdy*y+dfdz*z)+v4*f2dir;
	  break;
	  
	case isym_G4: //g4
	  x=r(2);
	  y=r(3);
	  z=r(4);
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;
          v4=x*y*(-v5 + 7. *v3);
	  dfdx= (-3.*v1 - 1.*v2 + 6.*v3)*y;
	  dfdy= (-1.*v1 - 3.*v2 + 6.*v3)*x;
	  dfdz= 12.*x*y*z;
	  h=fdir*v4;
          symvals(totf,0)=v4*func;
	  symvals(totf,1)=func*dfdx+h*x;
	  symvals(totf,2)=func*dfdy+h*y;
	  symvals(totf,3)=func*dfdz+h*z;
	  symvals(totf,4)=2*h+2*fdir*(dfdx*x+dfdy*y+dfdz*z)+v4*f2dir;
	  break;

       case isym_G5: //g5
	  x=r(2);
	  y=r(3);
	  z=r(4);
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;
          v4=x*z*(v1-3.*v2);
	  dfdx= 3*(v1 - v2)*z;
	  dfdy= -6*x*y*z;
	  dfdz= (v1 - 3*v2)*x;
	  h=fdir*v4;
          symvals(totf,0)=v4*func;
	  symvals(totf,1)=func*dfdx+h*x;
	  symvals(totf,2)=func*dfdy+h*y;
	  symvals(totf,3)=func*dfdz+h*z;
	  symvals(totf,4)=2*h+2*fdir*(dfdx*x+dfdy*y+dfdz*z)+v4*f2dir;
	  break;

	case isym_G6: //g6
	  x=r(2);
	  y=r(3);
	  z=r(4);
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;
          v4=y*z*(3*v1-v2);
	  dfdx= 6*x*y*z;
	  dfdy= 3*(v1 - v2)*z;
	  dfdz= (3*v1 - v2)*y;
	  h=fdir*v4;
          symvals(totf,0)=v4*func;
	  symvals(totf,1)=func*dfdx+h*x;
	  symvals(totf,2)=func*dfdy+h*y;
	  symvals(totf,3)=func*dfdz+h*z;
	  symvals(totf,4)=2*h+2*fdir*(dfdx*x+dfdy*y+dfdz*z)+v4*f2dir;
	  break;

	case isym_G7: //g7
	  x=r(2);
	  y=r(3);
	  z=r(4);
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;
          v4=(v1-3.*v2)*v1-(3.*v1-v2)*v2;
	  dfdx= 2*x*(2.*v1 - 6*v2);
	  dfdy= 2*y*(-6*v1 + 2*v2);
	  h=fdir*v4;
          symvals(totf,0)=v4*func;
	  symvals(totf,1)=func*dfdx+h*x;
	  symvals(totf,2)=func*dfdy+h*y;
	  symvals(totf,3)=h*z;
	  symvals(totf,4)=2*h+2*fdir*(dfdx*x+dfdy*y)+v4*f2dir;
	  break;

	case isym_G8: //g8
	  x=r(2);
	  y=r(3);
	  z=r(4);
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;
	  v4=x*y*(v1-v2);
	  dfdx= (3*v1 - v2)*y;
	  dfdy= (v1 - 3*v2)*x;
	  h=fdir*v4;
	  symvals(totf,0)=v4*func;
	  symvals(totf,1)=func*dfdx+h*x;
	  symvals(totf,2)=func*dfdy+h*y;
	  symvals(totf,3)=h*z;
	  symvals(totf,4)=2*h+2*fdir*(dfdx*x+dfdy*y)+v4*f2dir;
	  break;


        default:
          error("Bad symmetry in Cubic_spline::calcLap.\n");
        }
        //if(fabs(r(4)-3.39275) < 1e-3) { 
        //  cout << "r " << r(0) << " " << r(1) << "  " 
        //  << r(2) << "  " << r(3) << " " << r(4) << endl;
        //  cout << "symmvals " ;
        //  for(int i=0; i< 5; i++) cout << symvals(totf,i) << " ";
        //  cout << endl;
        //  cout << "interval " << interval << " func " << func << endl;
        //}
      }
    }
  }
}

//--------------------------------------------------------------------------
void Cubic_spline::calcHessian(const Array1 <doublevar> & r,
			       Array2 <doublevar> & symvals,
			       const int startfill) {
  assert(r.GetDim(0) >= 5);
  assert(symvals.GetDim(1)==10);
  assert(symvals.GetDim(0)>= startfill+nfunctions);
  //cout << "spline interval " << interval << "   r   " << r(0) << endl;
  if(r(0) >= threshold)
  {
    int end=startfill+nfunctions;
    for(int i=startfill; i< end; i++)
    {
      for(int d=0; d< 10; d++)
      {
        symvals(i,d)=0;
      }
    }
  }
  else
  {
    register doublevar v, h, v1, v2, v3, v4, v5, u2;
    register doublevar dgdx, dgdy, dgdz, d2gdxdx, d2gdydy, d2gdzdz, d2gdxdy, d2gdxdz, d2gdydz;
    register doublevar func, fdir, f2dir;
    int totf=startfill;
    doublevar ovr2=1.0/r(1);
    doublevar x2r2=r(2)*r(2)*ovr2;
    doublevar y2r2=r(3)*r(3)*ovr2;
    doublevar z2r2=r(4)*r(4)*ovr2;
    doublevar xyr2=r(2)*r(3)*ovr2;
    doublevar xzr2=r(2)*r(4)*ovr2;
    doublevar yzr2=r(3)*r(4)*ovr2;
    doublevar x=r(2),y=r(3),z=r(4);

    for(int i=0; i<nsplines; i++) {
      int interval=splines(i).getInterval(r(0));

      splines(i).getDers(r(0), interval,func, fdir, f2dir);
      //Writing everything in terms of the radial function f(r),
      //which makes it easy to derive the derivatives and hessian matrix
      doublevar fp=f2dir-fdir;
      doublevar dfdx2=x2r2*fp+fdir;
      doublevar dfdy2=y2r2*fp+fdir;
      doublevar dfdz2=z2r2*fp+fdir;
      doublevar dfdxy=xyr2*fp;
      doublevar dfdxz=xzr2*fp;
      doublevar dfdyz=yzr2*fp;
      doublevar dfdx=fdir*x;
      doublevar dfdy=fdir*y;
      doublevar dfdz=fdir*z;

      int end=totf+nfuncspline(i);
      for(; totf <end; totf++)
      {
	//assigning the pointer here so we can 
	//gradually change on a per-section basis
	doublevar * ptr=symvals.v+totf*10;
        switch(indiv_symmetry(totf-startfill))
        {

        case isym_S: // s
          *(ptr++)=func;
          *(ptr++)=dfdx;
          *(ptr++)=dfdy;
          *(ptr++)=dfdz;
          *(ptr++)=dfdx2;
	  *(ptr++)=dfdy2;
          *(ptr++)=dfdz2;
	  *(ptr++)=dfdxy;
          *(ptr++)=dfdxz;
	  *(ptr++)=dfdyz;
          break;

        case isym_Px: //p
          *(ptr++)=func*x;
          *(ptr++)=x*dfdx+func;
          *(ptr++)=x*dfdy;
          *(ptr++)=x*dfdz;
          *(ptr++)=x*dfdx2+2*dfdx;
	  *(ptr++)=x*dfdy2;
	  *(ptr++)=x*dfdz2;
	  *(ptr++)=dfdy+x*dfdxy;
	  *(ptr++)=dfdz+x*dfdxz;
	  *(ptr++)=x*dfdyz;
          break;

        case isym_Py:
          *(ptr++)=func*y;
          *(ptr++)=y*dfdx;
          *(ptr++)=y*dfdy+func;
          *(ptr++)=y*dfdz;
          *(ptr++)=y*dfdx2;
	  *(ptr++)=y*dfdy2+2*dfdy;
	  *(ptr++)=y*dfdz2;
	  *(ptr++)=dfdx+y*dfdxy;
	  *(ptr++)=y*dfdxz;
	  *(ptr++)=dfdz+y*dfdyz;
          break;

        case isym_Pz:
          *(ptr++)=func*z;
          *(ptr++)=z*dfdx;
          *(ptr++)=z*dfdy;
          *(ptr++)=z*dfdz+func;
          *(ptr++)=z*dfdx2;
	  *(ptr++)=z*dfdy2;
	  *(ptr++)=z*dfdz2+2*dfdz;
	  *(ptr++)=z*dfdxy;
	  *(ptr++)=dfdx+z*dfdxz;
	  *(ptr++)=dfdy+z*dfdyz;
          break;

        case isym_Dxx: //d
          v=x*x;  //xx
	  *(ptr++)=v*func;
	  *(ptr++)=2*x*func+v*dfdx;
	  *(ptr++)=v*dfdy;
	  *(ptr++)=v*dfdz;
	  *(ptr++)=2*func+4*x*dfdx+v*dfdx2;
	  *(ptr++)=v*dfdy2;
	  *(ptr++)=v*dfdz2;
	  *(ptr++)=2*x*dfdy+v*dfdxy;
	  *(ptr++)=2*x*dfdz+v*dfdxz;
	  *(ptr++)=v*dfdyz;
          break;

        case isym_Dyy:
          v=y*y; //yy
	  *(ptr++)=v*func;
	  *(ptr++)=v*dfdx;
	  *(ptr++)=v*dfdy+2*y*func;
	  *(ptr++)=v*dfdz;
	  *(ptr++)=v*dfdx2;
	  *(ptr++)=v*dfdy2+2*func+4*y*dfdy;
	  *(ptr++)=v*dfdz2;
	  *(ptr++)=2*y*dfdx+v*dfdxy;
	  *(ptr++)=v*dfdxz;
	  *(ptr++)=2*y*dfdz+v*dfdyz;
          break;

        case isym_Dzz:
          v=z*z; //zz
	  *(ptr++)=v*func;
	  *(ptr++)=v*dfdx;
	  *(ptr++)=v*dfdy;
	  *(ptr++)=v*dfdz+2*z*func;
	  *(ptr++)=v*dfdx2;
	  *(ptr++)=v*dfdy2;
	  *(ptr++)=v*dfdz2+2*func+4*z*dfdz;
	  *(ptr++)=v*dfdxy;
	  *(ptr++)=2*z*dfdx+v*dfdxz;
	  *(ptr++)=2*z*dfdy+v*dfdyz;
          break;

        case isym_Dxy:
          v=x*y; 
	  *(ptr++)=v*func;
	  *(ptr++)=v*dfdx+y*func;
	  *(ptr++)=v*dfdy+x*func;
	  *(ptr++)=v*dfdz;
	  *(ptr++)=v*dfdx2+2*y*dfdx;
	  *(ptr++)=v*dfdy2+2*x*dfdy;
	  *(ptr++)=v*dfdz2;
	  *(ptr++)=v*dfdxy+x*dfdx+y*dfdy+func;
	  *(ptr++)=v*dfdxz+y*dfdz;
	  *(ptr++)=v*dfdyz+x*dfdz;
          break;

        case isym_Dxz:
          v=x*z;  //xz
	  *(ptr++)=v*func;
	  *(ptr++)=v*dfdx+z*func;
	  *(ptr++)=v*dfdy;
	  *(ptr++)=v*dfdz+x*func;
	  *(ptr++)=v*dfdx2+2*z*dfdx;
	  *(ptr++)=v*dfdy2;
	  *(ptr++)=v*dfdz2+2*x*dfdz;
	  *(ptr++)=v*dfdxy+z*dfdy;
	  *(ptr++)=v*dfdxz+x*dfdx+z*dfdz+func;
	  *(ptr++)=v*dfdyz+x*dfdy;
          break;

        case isym_Dyz:
          v=y*z;  //yz
	  *(ptr++)=v*func;
	  *(ptr++)=v*dfdx;
	  *(ptr++)=v*dfdy+z*func;
	  *(ptr++)=v*dfdz+y*func;
	  *(ptr++)=v*dfdx2;
	  *(ptr++)=v*dfdy2+2*z*dfdy;
	  *(ptr++)=v*dfdz2+2*y*dfdz;
	  *(ptr++)=v*dfdxy+z*dfdx;
	  *(ptr++)=v*dfdxz+y*dfdx;
	  *(ptr++)=v*dfdyz+y*dfdy+z*dfdz+func;	  
          break;


        case isym_Dx2y2:
          v=x*x-y*y;
	  *(ptr++)=v*func;
	  *(ptr++)=v*dfdx+2*x*func;
	  *(ptr++)=v*dfdy-2*y*func;
	  *(ptr++)=v*dfdz;
	  *(ptr++)=v*dfdx2+4*x*dfdx+2*func;
	  *(ptr++)=v*dfdy2-4*y*dfdy-2*func;
	  *(ptr++)=v*dfdz2;
	  *(ptr++)=v*dfdxy-2*y*dfdx+2*x*dfdy;
	  *(ptr++)=v*dfdxz+2*x*dfdz;
	  *(ptr++)=v*dfdyz-2*y*dfdz;	
          break;

        case isym_Dz2r2:
          v=2.*z*z-x*x-y*y;
          *(ptr++)=v*func;
          *(ptr++)=v*dfdx-2*x*func;
          *(ptr++)=v*dfdy-2*y*func;
          *(ptr++)=v*dfdz+4*z*func;
          *(ptr++)=v*dfdx2-4*x*dfdx-2*func;
	  *(ptr++)=v*dfdy2-4*y*dfdy-2*func;
          *(ptr++)=v*dfdz2+8*z*dfdz+4*func;
	  *(ptr++)=v*dfdxy-2*x*dfdy-2*y*dfdx;
          *(ptr++)=v*dfdxz-2*x*dfdz+4*z*dfdx;
	  *(ptr++)=v*dfdyz-2*y*dfdz+4*z*dfdy;
          break;

          //f
        case isym_Fxxx: //xxx
          v2=x*x;
          v3=v2*x;
	  *(ptr++)=v3*func;
          *(ptr++)=3*v2*func+v3*dfdx;
          *(ptr++)=v3*dfdy;
          *(ptr++)=v3*dfdz;
          *(ptr++)=v3*dfdx2+6*v2*dfdx+6*x*func;
	  *(ptr++)=v3*dfdy2;
          *(ptr++)=v3*dfdz2;
	  *(ptr++)=v3*dfdxy+3*v2*dfdy;
          *(ptr++)=v3*dfdxz+3*v2*dfdz;
	  *(ptr++)=v3*dfdyz;       
          break;

        case isym_Fyyy://yyy
          v2=y*y;
          v3=v2*y;
   	  *(ptr++)=v3*func;
	  *(ptr++)=v3*dfdx;
          *(ptr++)=3*v2*func+v3*dfdy;
          *(ptr++)=v3*dfdz;
          *(ptr++)=v3*dfdx2;
	  *(ptr++)=v3*dfdy2+6*v2*dfdy+6*y*func;
	  *(ptr++)=v3*dfdz2;
	  *(ptr++)=3*v2*dfdx+dfdxy*v3;
	  *(ptr++)=v3*dfdxz;
	  *(ptr++)=v3*dfdyz+3*v2*dfdz;

          break;


        case isym_Fzzz: //zzz
          v2=z*z;
          v3=v2*z;
    	  *(ptr++)=v3*func;
	  *(ptr++)=v3*dfdx;
          *(ptr++)=v3*dfdy;
          *(ptr++)=3*v2*func+v3*dfdz;
          *(ptr++)=v3*dfdx2;
	  *(ptr++)=v3*dfdy2;
	  *(ptr++)=v3*dfdz2+6*v2*dfdz+6*z*func;
	  *(ptr++)=v3*dfdxy;
	  *(ptr++)=3*v2*dfdx+dfdxz*v3;
	  *(ptr++)=3*v2*dfdy+dfdyz*v3;

          break;

        case isym_Fxxy: //xxy
          v2=x*x;
          v3=v2*y;
          *(ptr++)=v3*func;
          *(ptr++)=v3*dfdx+2*x*y*func;
          *(ptr++)=v3*dfdy+v2*func;
          *(ptr++)=v3*dfdz;
          *(ptr++)=v3*dfdx2+4*x*y*dfdx+2*y*func;
	  *(ptr++)=v3*dfdy2+2*v2*dfdy;
          *(ptr++)=v3*dfdz2;
	  *(ptr++)=v3*dfdxy+v2*dfdx+2*x*y*dfdy+2*x*func;
          *(ptr++)=v3*dfdxz+2*x*y*dfdz;
	  *(ptr++)=v3*dfdyz+v2*dfdz;
          break;

        case isym_Fxxz://xxz
          v2=x*x;
          v3=v2*z;
          *(ptr++)=v3*func;
          *(ptr++)=v3*dfdx+2*x*z*func; 
          *(ptr++)=v3*dfdy;
          *(ptr++)=v3*dfdz+v2*func;
          *(ptr++)=v3*dfdx2+4*x*z*dfdx+2*z*func;
	  *(ptr++)=v3*dfdy2;
          *(ptr++)=v3*dfdz2+2*v2*dfdz;
	  *(ptr++)=v3*dfdxy+2*x*z*dfdy;
	  *(ptr++)=v3*dfdxz+v2*dfdx+2*x*z*dfdz+2*x*func;
          *(ptr++)=v3*dfdyz+v2*dfdy;
          break;

        case isym_Fyyx://yyx
          v2=y*y;
          v3=v2*x;
          *(ptr++)=v3*func;   
          *(ptr++)=v3*dfdx+v2*func;
          *(ptr++)=v3*dfdy+2*y*x*func; 
          *(ptr++)=v3*dfdz;
	  *(ptr++)=v3*dfdx2+2*v2*dfdx;
          *(ptr++)=v3*dfdy2+4*y*x*dfdy+2*x*func;
          *(ptr++)=v3*dfdz2;
	  *(ptr++)=v3*dfdxy+2*x*y*dfdx+v2*dfdy+2*y*func;
	  *(ptr++)=v3*dfdxz+v2*dfdz;
          *(ptr++)=v3*dfdyz+2*x*y*dfdz;


          break;

        case isym_Fyyz://yyz
          v2=y*y;
          v3=v2*z;
          *(ptr++)=v3*func;   
          *(ptr++)=v3*dfdx; 
          *(ptr++)=v3*dfdy+2*y*z*func;
          *(ptr++)=v3*dfdz+v2*func;
          *(ptr++)=v3*dfdx2; 
          *(ptr++)=v3*dfdy2+4*y*z*dfdy+2*z*func;
	  *(ptr++)=v3*dfdz2+2*v2*dfdz;
	  *(ptr++)=v3*dfdxy+2*y*z*dfdx;
	  *(ptr++)=v3*dfdxz+v2*dfdx;
          *(ptr++)=v3*dfdyz+v2*dfdy+2*y*z*dfdz+2*y*func;

          break;

        case isym_Fzzx://zzx
          v2=z*z;
          v3=x*v2;
          *(ptr++)=v3*func;
          *(ptr++)=v3*dfdx+v2*func;
          *(ptr++)=v3*dfdy;
          *(ptr++)=v3*dfdz+2*z*x*func; 
          *(ptr++)=v3*dfdx2+2*v2*dfdx;
	  *(ptr++)=v3*dfdy2;
          *(ptr++)=v3*dfdz2+4*z*x*dfdz+2*x*func;
	  *(ptr++)=v3*dfdxy+v2*dfdy;
	  *(ptr++)=v3*dfdxz+2*z*x*dfdx+v2*dfdz+2*z*func;
          *(ptr++)=v3*dfdyz+2*z*x*dfdy;

          break;

        case isym_Fzzy://zzy
          v2=z*z;
          v3=v2*y;
          *(ptr++)=v3*func;
          *(ptr++)=v3*dfdx;
          *(ptr++)=v3*dfdy+v2*func; 
          *(ptr++)=v3*dfdz+2*z*y*func; 
	  *(ptr++)=v3*dfdx2;
          *(ptr++)=v3*dfdy2+2*v2*dfdy;
          *(ptr++)=v3*dfdz2+4*z*y*dfdz+2*y*func;
	  *(ptr++)=v3*dfdxy+v2*dfdx;
	  *(ptr++)=v3*dfdxz+2*z*y*dfdx;
          *(ptr++)=v3*dfdyz+2*z*y*dfdy+v2*dfdz+2*z*func;
          break;

        case isym_Fxyz://xyz
	  
          v3=x*y*z;
          *(ptr++)=v3*func;
          *(ptr++)=v3*dfdx+y*z*func;
          *(ptr++)=v3*dfdy+x*z*func; 
          *(ptr++)=v3*dfdz+x*y*func; 
	  *(ptr++)=v3*dfdx2+2*y*z*dfdx;
          *(ptr++)=v3*dfdy2+2*x*z*dfdy;
          *(ptr++)=v3*dfdz2+2*x*y*dfdz;
	  *(ptr++)=v3*dfdxy+y*z*dfdy+z*x*dfdx+z*func;
	  *(ptr++)=v3*dfdxz+dfdz*y*z+x*y*dfdx+y*func;
          *(ptr++)=v3*dfdyz+x*y*dfdy+x*z*dfdz+x*func;
	  
	  /*
	  //new way does the same thing 
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;

          v4=x*y*z;
	  dgdx=y*z;
	  dgdy=x*z;
	  dgdz=x*y;
	  d2gdxdx=0;
	  d2gdydy=0;
	  d2gdzdz=0;
	  d2gdxdy=z;
	  d2gdxdz=y;
	  d2gdydz=x;
	  
	  h=fdir*v4;
          *(ptr++)=v4*func;
	  *(ptr++)=func*dgdx+h*x;
	  *(ptr++)=func*dgdy+h*y;
	  *(ptr++)=func*dgdz+h*z;
	  *(ptr++)=func*d2gdxdx+2*fdir*dgdx*x+v4*fdir+v4*(v1/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydy+2*fdir*dgdy*y+v4*fdir+v4*(v2/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdzdz+2*fdir*dgdz*z+v4*fdir+v4*(v3/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdy+y*fdir*dgdx+x*fdir*dgdy+v4*(x*y/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdz+z*fdir*dgdx+x*fdir*dgdz+v4*(x*z/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydz+z*fdir*dgdy+y*fdir*dgdz+v4*(y*z/v5)*(f2dir-fdir);
	  */
          break;

        case isym_Fm3:
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;

          v4=(3*v1-v2)*y;
	  dgdx=6*x*y;
	  dgdy=3*(v1 - v2);
	  dgdz=0;
	  d2gdxdx=6*y;
	  d2gdydy=-6*y;
	  d2gdzdz=0;
	  d2gdxdy=6*x;
	  d2gdxdz=0;
	  d2gdydz=0;
	  
	  h=fdir*v4;
          *(ptr++)=v4*func;
	  *(ptr++)=func*dgdx+h*x;
	  *(ptr++)=func*dgdy+h*y;
	  *(ptr++)=h*z;
	  *(ptr++)=func*d2gdxdx+2*fdir*dgdx*x+v4*fdir+v4*(v1/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydy+2*fdir*dgdy*y+v4*fdir+v4*(v2/v5)*(f2dir-fdir);
	  *(ptr++)=                           v4*fdir+v4*(v3/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdy+y*fdir*dgdx+x*fdir*dgdy+v4*(x*y/v5)*(f2dir-fdir);
	  *(ptr++)=             z*fdir*dgdx+v4*(x*z/v5)*(f2dir-fdir);
	  *(ptr++)=             z*fdir*dgdy+v4*(y*z/v5)*(f2dir-fdir);
	  break;
        case isym_Fm1:
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;

          v4=y*(4*v3-v1-v2);
	  dgdx=-2*x*y;
	  dgdy=-v1 - 3*v2 + 4*v3;
	  dgdz=8*y*z;
	  d2gdxdx=-2*y;
	  d2gdydy=-6*y;
	  d2gdzdz=8*y;
	  d2gdxdy=-2*x;
	  d2gdxdz=0.0;
	  d2gdydz=8*z;
	  
	  h=fdir*v4;
	  *(ptr++)=v4*func;
	  *(ptr++)=func*dgdx+h*x;
	  *(ptr++)=func*dgdy+h*y;
	  *(ptr++)=func*dgdz+h*z;
	  *(ptr++)=func*d2gdxdx+2*fdir*dgdx*x+v4*fdir+v4*(v1/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydy+2*fdir*dgdy*y+v4*fdir+v4*(v2/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdzdz+2*fdir*dgdz*z+v4*fdir+v4*(v3/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdy+y*fdir*dgdx+x*fdir*dgdy+v4*(x*y/v5)*(f2dir-fdir);
	  *(ptr++)=             z*fdir*dgdx+x*fdir*dgdz+v4*(x*z/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydz+z*fdir*dgdy+y*fdir*dgdz+v4*(y*z/v5)*(f2dir-fdir);
	  break;
        case isym_F0:
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;

          v4=z*(2*v3-3*v1-3*v2);
	  dgdx=-6*x*z;
	  dgdy=-6*y*z;
	  dgdz=-3*(v1 + v2 - 2*v3);
	  d2gdxdx=-6*z;
	  d2gdydy=-6*z;
	  d2gdzdz=12*z;
	  d2gdxdy=0;
	  d2gdxdz=-6*x;
	  d2gdydz=-6*y;
	  
	  h=fdir*v4;
          *(ptr++)=v4*func;
	  *(ptr++)=func*dgdx+h*x;
	  *(ptr++)=func*dgdy+h*y;
	  *(ptr++)=func*dgdz+h*z;
	  *(ptr++)=func*d2gdxdx+2*fdir*dgdx*x+v4*fdir+v4*(v1/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydy+2*fdir*dgdy*y+v4*fdir+v4*(v2/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdzdz+2*fdir*dgdz*z+v4*fdir+v4*(v3/v5)*(f2dir-fdir);
	  *(ptr++)=             y*fdir*dgdx+x*fdir*dgdy+v4*(x*y/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdz+z*fdir*dgdx+x*fdir*dgdz+v4*(x*z/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydz+z*fdir*dgdy+y*fdir*dgdz+v4*(y*z/v5)*(f2dir-fdir);
	  break;
        case isym_Fp1:
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;

          v4=x*(4*v3-v1-v2);
	  dgdx=-3*v1 - v2 + 4*v3;
	  dgdy=-2*x*y;
	  dgdz=8*x*z;
	  d2gdxdx=-6*x;
	  d2gdydy=-2*x;
	  d2gdzdz=8*x;
	  d2gdxdy=-2*y;
	  d2gdxdz=8*z;
	  d2gdydz=0.0;
	  
	  h=fdir*v4;
          *(ptr++)=v4*func;
	  *(ptr++)=func*dgdx+h*x;
	  *(ptr++)=func*dgdy+h*y;
	  *(ptr++)=func*dgdz+h*z;
	  *(ptr++)=func*d2gdxdx+2*fdir*dgdx*x+v4*fdir+v4*(v1/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydy+2*fdir*dgdy*y+v4*fdir+v4*(v2/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdzdz+2*fdir*dgdz*z+v4*fdir+v4*(v3/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdy+y*fdir*dgdx+x*fdir*dgdy+v4*(x*y/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdz+z*fdir*dgdx+x*fdir*dgdz+v4*(x*z/v5)*(f2dir-fdir);
	  *(ptr++)=             z*fdir*dgdy+y*fdir*dgdz+v4*(y*z/v5)*(f2dir-fdir);
	  break;
        case isym_Fp2:
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;

          v4=(v1-v2)*z;
	  dgdx=2*x*z;
	  dgdy=-2*y*z;
	  dgdz=v1 - v2;
	  d2gdxdx=2*z;
	  d2gdydy=-2*z;
	  d2gdzdz=0;
	  d2gdxdy=0;
	  d2gdxdz=2*x;
	  d2gdydz=-2*y;
	  
	  h=fdir*v4;
          *(ptr++)=v4*func;
	  *(ptr++)=func*dgdx+h*x;
	  *(ptr++)=func*dgdy+h*y;
	  *(ptr++)=func*dgdz+h*z;
	  *(ptr++)=func*d2gdxdx+2*fdir*dgdx*x+v4*fdir+v4*(v1/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydy+2*fdir*dgdy*y+v4*fdir+v4*(v2/v5)*(f2dir-fdir);
	  *(ptr++)=            2*fdir*dgdz*z+v4*fdir+v4*(v3/v5)*(f2dir-fdir);
	  *(ptr++)=             y*fdir*dgdx+x*fdir*dgdy+v4*(x*y/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdz+z*fdir*dgdx+x*fdir*dgdz+v4*(x*z/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydz+z*fdir*dgdy+y*fdir*dgdz+v4*(y*z/v5)*(f2dir-fdir);
	  break;
	case isym_Fp3mod:
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;

          v4=x*(v1-3.0*v2);
	  dgdx=3*(v1 - v2);
	  dgdy=-6*x*y;
	  dgdz=0;
	  d2gdxdx=6*x;
	  d2gdydy=-6*x;
	  d2gdzdz=0;
	  d2gdxdy=-6*y;
	  d2gdxdz=0;
	  d2gdydz=0;
	  
	  h=fdir*v4;
          *(ptr++)=v4*func;
	  *(ptr++)=func*dgdx+h*x;
	  *(ptr++)=func*dgdy+h*y;
	  *(ptr++)=          h*z;
	  *(ptr++)=func*d2gdxdx+2*fdir*dgdx*x+v4*fdir+v4*(v1/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydy+2*fdir*dgdy*y+v4*fdir+v4*(v2/v5)*(f2dir-fdir);
	  *(ptr++)=                          +v4*fdir+v4*(v3/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdy+y*fdir*dgdx+x*fdir*dgdy+v4*(x*y/v5)*(f2dir-fdir);
	  *(ptr++)=             z*fdir*dgdx            +v4*(x*z/v5)*(f2dir-fdir);
	  *(ptr++)=             z*fdir*dgdy            +v4*(y*z/v5)*(f2dir-fdir);
	  break;
        case isym_Fp3:
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;

          v4=x*(v1+v2);
	  dgdx=3*v1 + v2;
	  dgdy=2*x*y;
	  dgdz=0.0;
	  d2gdxdx=6*x;
	  d2gdydy=2*x;
	  d2gdzdz=0;
	  d2gdxdy=2*y;
	  d2gdxdz=0;
	  d2gdydz=0;
	  
	  h=fdir*v4;
          *(ptr++)=v4*func;
	  *(ptr++)=func*dgdx+h*x;
	  *(ptr++)=func*dgdy+h*y;
	  *(ptr++)=          h*z;
	  *(ptr++)=func*d2gdxdx+2*fdir*dgdx*x+v4*fdir+v4*(v1/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydy+2*fdir*dgdy*y+v4*fdir+v4*(v2/v5)*(f2dir-fdir);
	  *(ptr++)=                          +v4*fdir+v4*(v3/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdy+y*fdir*dgdx+x*fdir*dgdy+v4*(x*y/v5)*(f2dir-fdir);
	  *(ptr++)=             z*fdir*dgdx            +v4*(x*z/v5)*(f2dir-fdir);
	  *(ptr++)=             z*fdir*dgdy            +v4*(y*z/v5)*(f2dir-fdir);
	  break;

        case isym_Gxxxx: // g xxxx
	  v2=x*x;
	  v3=v2*x;
	  v4=v2*v2;
          *(ptr++)=v4*func;
          *(ptr++)=v4*dfdx+4*v3*func;
          *(ptr++)=v4*dfdy; 
          *(ptr++)=v4*dfdz; 
	  *(ptr++)=v4*dfdx2+8*v3*dfdx+12*v2*func;
          *(ptr++)=v4*dfdy2;
          *(ptr++)=v4*dfdz2;
	  *(ptr++)=v4*dfdxy+4*v3*dfdy;
	  *(ptr++)=v4*dfdxz+4*v3*dfdz;
          *(ptr++)=v4*dfdyz;

          break;

        case isym_Gyyyy: // g yyyy
	  v2=y*y;
	  v3=v2*y;
	  v4=v2*v2;
          *(ptr++)=v4*func;
          *(ptr++)=v4*dfdx;
          *(ptr++)=v4*dfdy+4*v3*func; 
          *(ptr++)=v4*dfdz; 
          *(ptr++)=v4*dfdx2;
	  *(ptr++)=v4*dfdy2+8*v3*dfdy+12*v2*func;
          *(ptr++)=v4*dfdz2;
	  *(ptr++)=v4*dfdxy+4*v3*dfdx;
          *(ptr++)=v4*dfdxz;
	  *(ptr++)=v4*dfdyz+4*v3*dfdz;

	  break;

        case isym_Gzzzz: // g zzzz
	  v2=z*z;
	  v3=v2*z;
	  v4=v2*v2;
          *(ptr++)=v4*func;
          *(ptr++)=v4*dfdx;
          *(ptr++)=v4*dfdy;
          *(ptr++)=v4*dfdz+4*v3*func; 
          *(ptr++)=v4*dfdx2;
          *(ptr++)=v4*dfdy2;
	  *(ptr++)=v4*dfdz2+8*v3*dfdz+12*v2*func;
	  *(ptr++)=v4*dfdxy;
	  *(ptr++)=v4*dfdxz+4*v3*dfdx;
          *(ptr++)=v4*dfdyz+4*v3*dfdy;
          break;

        case isym_Gxxxy: // g xxxy
	  v2=x*x;
	  v3=v2*x;
	  v4=v3*y;
          *(ptr++)=v4*func;
          *(ptr++)=v4*dfdx+3*v2*y*func;
          *(ptr++)=v4*dfdy+v3*func;
          *(ptr++)=v4*dfdz; 
          *(ptr++)=v4*dfdx2+6*v2*y*dfdx+6*x*y*func;
          *(ptr++)=v4*dfdy2+2*v3*dfdy;
	  *(ptr++)=v4*dfdz2;
	  *(ptr++)=v4*dfdxy+v3*dfdx+3*v2*y*dfdy+3*v2*func;
	  *(ptr++)=v4*dfdxz+3*v2*y*dfdz;
          *(ptr++)=v4*dfdyz+v3*dfdz;

	  break;

        case isym_Gxxxz: // g xxxz
	  v2=x*x;
	  v3=v2*x;
	  v4=v3*z;
          *(ptr++)=v4*func;
          *(ptr++)=v4*dfdx+3*v2*z*func;
          *(ptr++)=v4*dfdy; 
          *(ptr++)=v4*dfdz+v3*func;
          *(ptr++)=v4*dfdx2+6*v2*z*dfdx+6*x*z*func;
	  *(ptr++)=v4*dfdy2;
          *(ptr++)=v4*dfdz2+2*v3*dfdz;
	  *(ptr++)=v4*dfdxy+3*v2*z*dfdy;
	  *(ptr++)=v4*dfdxz+v3*dfdx+3*v2*z*dfdz+3*v2*func;
          *(ptr++)=v4*dfdyz+v3*dfdy;
          break;

        case isym_Gyyyx: // g yyyx
	  v2=y*y;
	  v3=v2*y;
	  v4=v3*x;
          *(ptr++)=v4*func;
          *(ptr++)=v4*dfdx+v3*func;
          *(ptr++)=v4*dfdy+3*v2*x*func;
          *(ptr++)=v4*dfdz; 
          *(ptr++)=v4*dfdx2+2*v3*dfdx;
          *(ptr++)=v4*dfdy2+6*v2*x*dfdy+6*y*x*func;
	  *(ptr++)=v4*dfdz2;
	  *(ptr++)=v4*dfdxy+3*v2*x*dfdx+v3*dfdy+3*v2*func;
	  *(ptr++)=v4*dfdxz+v3*dfdz;
          *(ptr++)=v4*dfdyz+3*v2*x*dfdz;
 
          break;

        case isym_Gyyyz: // g yyyz
	  v2=y*y;
	  v3=v2*y;
	  v4=v3*z;
          *(ptr++)=v4*func;
          *(ptr++)=v4*dfdx; 
          *(ptr++)=v4*dfdy+3*v2*z*func;
          *(ptr++)=v4*dfdz+v3*func;
	  *(ptr++)=v4*dfdx2;
          *(ptr++)=v4*dfdy2+6*v2*z*dfdy+6*y*z*func;
          *(ptr++)=v4*dfdz2+2*v3*dfdz;
	  *(ptr++)=v4*dfdxy+3*v2*z*dfdx;
	  *(ptr++)=v4*dfdxz+v3*dfdx;
          *(ptr++)=v4*dfdyz+v3*dfdy+3*v2*z*dfdz+3*v2*func;

          break;

        case isym_Gzzzx: // g zzzx
	  v2=z*z;
	  v3=v2*z;
	  v4=v3*x;
          *(ptr++)=v4*func;
          *(ptr++)=v4*dfdx+v3*func;
          *(ptr++)=v4*dfdy;
          *(ptr++)=v4*dfdz+3*v2*x*func;
          *(ptr++)=v4*dfdx2+2*v3*dfdx;
 	  *(ptr++)=v4*dfdy2;
          *(ptr++)=v4*dfdz2+6*v2*x*dfdz+6*z*x*func;

	  *(ptr++)=v4*dfdxy+v3*dfdy;
	  *(ptr++)=v4*dfdxz+3*v2*x*dfdx+3*v2*func+dfdz*v3;
          *(ptr++)=v4*dfdyz+3*v2*x*dfdy;

          break;

        case isym_Gzzzy: // g zzzy
	  v2=z*z;
	  v3=v2*z;
	  v4=v3*y;
          *(ptr++)=v4*func;
          *(ptr++)=v4*dfdx;
          *(ptr++)=v4*dfdy+v3*func;  
          *(ptr++)=v4*dfdz+3*v2*y*func;
 	  *(ptr++)=v4*dfdx2;
          *(ptr++)=v4*dfdy2+2*v3*dfdy; 
          *(ptr++)=v4*dfdz2+6*v2*y*dfdz+6*z*y*func;
	  *(ptr++)=v4*dfdxy+v3*dfdx;
          *(ptr++)=v4*dfdxz+3*v2*y*dfdx;
	  *(ptr++)=v4*dfdyz+3*v2*y*dfdy+3*v2*func+dfdz*v3;
 
          break;

        case isym_Gxxyy: // g xxyy
	  v2=x*x;
	  u2=y*y;
	  v4=v2*u2;
          *(ptr++)=v4*func;
          *(ptr++)=v4*dfdx+2*x*u2*func;
          *(ptr++)=v4*dfdy+2*y*v2*func;  
          *(ptr++)=v4*dfdz;
 	  *(ptr++)=v4*dfdx2+4*x*u2*dfdx+2*u2*func;
 	  *(ptr++)=v4*dfdy2+4*y*v2*dfdy+2*v2*func;
          *(ptr++)=v4*dfdz2;
	  *(ptr++)=v4*dfdxy+4*x*y*func+2*x*u2*dfdy+v2*2*y*dfdx;
          *(ptr++)=v4*dfdxz+2*x*u2*dfdz;
	  *(ptr++)=v4*dfdyz+2*y*v2*dfdz;

	  break;

        case isym_Gxxzz: // g xxzz
	  v2=x*x;
	  u2=z*z;
	  v4=v2*u2;
          *(ptr++)=v4*func;
          *(ptr++)=v4*dfdx+2*x*u2*func;
          *(ptr++)=v4*dfdy;
          *(ptr++)=v4*dfdz+2*z*v2*func; 
 	  *(ptr++)=v4*dfdx2+4*x*u2*dfdx+2*u2*func;
          *(ptr++)=v4*dfdy2;
 	  *(ptr++)=v4*dfdz2+4*z*v2*dfdz+2*v2*func;
          *(ptr++)=v4*dfdxy+2*x*u2*dfdy;
	  *(ptr++)=v4*dfdxz+4*x*z*func+2*x*u2*dfdz+v2*2*z*dfdx;     
	  *(ptr++)=v4*dfdyz+2*v2*z*dfdy;
          break;

        case isym_Gyyzz: // g yyzz
	  v2=y*y;
	  u2=z*z;
	  v4=v2*u2;
          *(ptr++)=v4*func;
          *(ptr++)=v4*dfdx;
          *(ptr++)=v4*dfdy+2*y*u2*func;
          *(ptr++)=v4*dfdz+2*z*v2*func; 
          *(ptr++)=v4*dfdx2;
	  *(ptr++)=v4*dfdy2+4*y*u2*dfdy+2*u2*func; 
 	  *(ptr++)=v4*dfdz2+4*z*v2*dfdz+2*v2*func;
          *(ptr++)=v4*dfdxy+2*y*u2*dfdx;
	  *(ptr++)=v4*dfdxz+2*v2*z*dfdx;
          *(ptr++)=v4*dfdyz+4*y*z*func+2*y*u2*dfdz+v2*2*z*dfdy;     

          break;

        case isym_Gxxyz: // g xxyz
	  v2=x*x;
	  u2=y*z;
	  v4=v2*u2;
          *(ptr++)=v4*func;
          *(ptr++)=v4*dfdx+2*x*u2*func;
          *(ptr++)=v4*dfdy+v2*z*func;
          *(ptr++)=v4*dfdz+v2*y*func; 
          *(ptr++)=v4*dfdx2+2*u2*func+4*dfdx*x*u2;
	  *(ptr++)=v4*dfdy2+2*v2*z*dfdy; 
 	  *(ptr++)=v4*dfdz2+2*v2*y*dfdz;
          *(ptr++)=v4*dfdxy+2*x*z*func+dfdy*2*x*u2+v2*z*dfdx;
	  *(ptr++)=v4*dfdxz+2*x*y*func+dfdz*2*x*u2+v2*y*dfdx;
          *(ptr++)=v4*dfdyz+v2*func+dfdz*v2*z+v2*y*dfdy;     
	
          break;

        case isym_Gyyxz: // g yyxz
	  v2=y*y;
	  u2=x*z;
	  v4=v2*u2;
          *(ptr++)=v4*func;
          *(ptr++)=v4*dfdx+v2*z*func;
          *(ptr++)=v4*dfdy+2*y*u2*func; 
          *(ptr++)=v4*dfdz+v2*x*func; 
     	  *(ptr++)=v4*dfdx2+2*v2*z*dfdx; 
          *(ptr++)=v4*dfdy2+2*u2*func+4*dfdy*y*u2;
 	  *(ptr++)=v4*dfdz2+2*v2*x*dfdz;
          *(ptr++)=v4*dfdxy+dfdx*2*y*u2+2*y*z*func+v2*z*dfdy;
	  *(ptr++)=v4*dfdxz+v2*x*dfdx+v2*func+dfdz*v2*z;
          *(ptr++)=v4*dfdyz+v2*x*dfdy+2*y*x*func+dfdz*2*y*u2;     
	
          break;

        case isym_Gzzxy: // g zzxy
	  v2=z*z;
	  u2=x*y;
	  v4=v2*u2;
          *(ptr++)=v4*func;
          *(ptr++)=v4*dfdx+v2*y*func;
          *(ptr++)=v4*dfdy+v2*x*func;
          *(ptr++)=v4*dfdz+2*z*u2*func;         
     	  *(ptr++)=v4*dfdx2+2*v2*y*dfdx; 
          *(ptr++)=v4*dfdy2+2*v2*x*dfdy;
          *(ptr++)=v4*dfdz2+2*u2*func+4*dfdz*z*u2;
 	  *(ptr++)=v4*dfdxy+v2*x*dfdx+v2*func+dfdy*v2*y;
          *(ptr++)=v4*dfdxz+dfdx*2*z*u2+2*y*z*func+v2*y*dfdz;	
          *(ptr++)=v4*dfdyz+2*z*u2*dfdy+2*z*x*func+dfdz*v2*x;  
          break;


	case isym_G0:
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;

          v4=35.*v3*v3-30.*v3*v5+3.*v5*v5;
	  dgdx=  12.0*(v1 + v2 - 4*v3)*x;
	  dgdy=  12.0*(v1 + v2 - 4*v3)*y;
	  dgdz= -48.0*(v1 + v2 + v3)*z + 80*v3*z;
	  d2gdxdx=12*(3*v1 + v2 - 4*v3);
	  d2gdydy=12*(v1 + 3*v2 - 4*v3);
	  d2gdzdz=-48*(v1 + v2 - 2*v3);
	  d2gdxdy=24*x*y;
	  d2gdxdz=-96*x*z;
	  d2gdydz=-96*y*z;
	  
	  h=fdir*v4;
          *(ptr++)=v4*func;
	  *(ptr++)=func*dgdx+h*x;
	  *(ptr++)=func*dgdy+h*y;
	  *(ptr++)=func*dgdz+h*z;
	  *(ptr++)=func*d2gdxdx+2*fdir*dgdx*x+v4*fdir+v4*(v1/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydy+2*fdir*dgdy*y+v4*fdir+v4*(v2/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdzdz+2*fdir*dgdz*z+v4*fdir+v4*(v3/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdy+y*fdir*dgdx+x*fdir*dgdy+v4*(x*y/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdz+z*fdir*dgdx+x*fdir*dgdz+v4*(x*z/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydz+z*fdir*dgdy+y*fdir*dgdz+v4*(y*z/v5)*(f2dir-fdir);
	  break;
        case isym_G1:
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;

	  v4=x*z*(7.*v3-3.*v5);
	  dgdx= (-9.*v1 - 3.*v2 + 4.*v3)*z;
	  dgdy= -6.*x*y*z;
	  dgdz= (-3.*v1 - 3.*v2 + 12.*v3)*x;
	  d2gdxdx=-18.*x*z;
	  d2gdydy=-6.*x*z;
	  d2gdzdz=24.*x*z;
	  d2gdxdy=-6.*y*z;
	  d2gdxdz=-9.*v1 - 3.*v2 + 12.*v3;
	  d2gdydz=-6.*x*y;
	  
	  h=fdir*v4;
          *(ptr++)=v4*func;
	  *(ptr++)=func*dgdx+h*x;
	  *(ptr++)=func*dgdy+h*y;
	  *(ptr++)=func*dgdz+h*z;
	  *(ptr++)=func*d2gdxdx+2*fdir*dgdx*x+v4*fdir+v4*(v1/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydy+2*fdir*dgdy*y+v4*fdir+v4*(v2/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdzdz+2*fdir*dgdz*z+v4*fdir+v4*(v3/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdy+y*fdir*dgdx+x*fdir*dgdy+v4*(x*y/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdz+z*fdir*dgdx+x*fdir*dgdz+v4*(x*z/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydz+z*fdir*dgdy+y*fdir*dgdz+v4*(y*z/v5)*(f2dir-fdir);
	  break;
        case isym_G2:
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;

	  v4=y*z*(7.*v3-3.*v5);
	  dgdx= -6.*x*y*z;
	  dgdy= (-3.*v1 - 9.*v2 + 4.*v3)*z;
	  dgdz= (-3.*v1 - 3.*v2 + 12.*v3)*y;
	  d2gdxdx=-6.*y*z;
	  d2gdydy=-18.*y*z;
	  d2gdzdz=24.*y*z;
	  d2gdxdy=-6.*x*z;
	  d2gdxdz=-6.*x*y;
	  d2gdydz=-3.*v1 - 9.*v2 + 12.*v3;
	  
	  h=fdir*v4;
          *(ptr++)=v4*func;
	  *(ptr++)=func*dgdx+h*x;
	  *(ptr++)=func*dgdy+h*y;
	  *(ptr++)=func*dgdz+h*z;
	  *(ptr++)=func*d2gdxdx+2*fdir*dgdx*x+v4*fdir+v4*(v1/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydy+2*fdir*dgdy*y+v4*fdir+v4*(v2/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdzdz+2*fdir*dgdz*z+v4*fdir+v4*(v3/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdy+y*fdir*dgdx+x*fdir*dgdy+v4*(x*y/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdz+z*fdir*dgdx+x*fdir*dgdz+v4*(x*z/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydz+z*fdir*dgdy+y*fdir*dgdz+v4*(y*z/v5)*(f2dir-fdir);
	  break;
        case isym_G3:
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;

	  v4=(v1 - v2)*(-v5 + 7. *v3);
	  dgdx= -4.*v1*x + 12.*v3*x;
	  dgdy= 4.*v2*y - 12.*v3*y;
	  dgdz= 12.*(v1 - 1.*v2)*z;
	  d2gdxdx=-12.*v1 + 12.*v3;
	  d2gdydy=12.*v2 - 12.*v3;
	  d2gdzdz=12.*v1 - 12.*v2;
	  d2gdxdy=0;
	  d2gdxdz=24.*x*z;
	  d2gdydz=-24.*y*z;
	  
	  h=fdir*v4;
          *(ptr++)=v4*func;
	  *(ptr++)=func*dgdx+h*x;
	  *(ptr++)=func*dgdy+h*y;
	  *(ptr++)=func*dgdz+h*z;
	  *(ptr++)=func*d2gdxdx+2*fdir*dgdx*x+v4*fdir+v4*(v1/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydy+2*fdir*dgdy*y+v4*fdir+v4*(v2/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdzdz+2*fdir*dgdz*z+v4*fdir+v4*(v3/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdy+y*fdir*dgdx+x*fdir*dgdy+v4*(x*y/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdz+z*fdir*dgdx+x*fdir*dgdz+v4*(x*z/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydz+z*fdir*dgdy+y*fdir*dgdz+v4*(y*z/v5)*(f2dir-fdir);
	  break;
        case isym_G4:
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;
	  
	  v4=x*y*(-v5 + 7. *v3);
	  dgdx= (-3.*v1 - 1.*v2 + 6.*v3)*y;
	  dgdy= (-1.*v1 - 3.*v2 + 6.*v3)*x;
	  dgdz= 12.*x*y*z;
	  d2gdxdx=-6*x*y;
	  d2gdydy=-6*x*y;
	  d2gdzdz=12.*x*y;
	  d2gdxdy=-3*v1 - 3*v2 + 6.*v3;
	  d2gdxdz=12.*y*z;
	  d2gdydz=12.*x*z;
	  
	  h=fdir*v4;
          *(ptr++)=v4*func;
	  *(ptr++)=func*dgdx+h*x;
	  *(ptr++)=func*dgdy+h*y;
	  *(ptr++)=func*dgdz+h*z;
	  *(ptr++)=func*d2gdxdx+2*fdir*dgdx*x+v4*fdir+v4*(v1/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydy+2*fdir*dgdy*y+v4*fdir+v4*(v2/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdzdz+2*fdir*dgdz*z+v4*fdir+v4*(v3/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdy+y*fdir*dgdx+x*fdir*dgdy+v4*(x*y/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdz+z*fdir*dgdx+x*fdir*dgdz+v4*(x*z/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydz+z*fdir*dgdy+y*fdir*dgdz+v4*(y*z/v5)*(f2dir-fdir);
	  break;
	case isym_G5:
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;

	  v4=x*z*(v1-3.*v2);
	  dgdx= 3*(v1 - v2)*z;
	  dgdy= -6*x*y*z;
	  dgdz= (v1 - 3*v2)*x;
	  d2gdxdx=6*x*z;
	  d2gdydy=-6*x*z;
	  d2gdzdz=0;
	  d2gdxdy=-6*y*z;
	  d2gdxdz=3*(v1 - v2);
	  d2gdydz=-6*x*y;
	  
	  h=fdir*v4;
          *(ptr++)=v4*func;
	  *(ptr++)=func*dgdx+h*x;
	  *(ptr++)=func*dgdy+h*y;
	  *(ptr++)=func*dgdz+h*z;
	  *(ptr++)=func*d2gdxdx+2*fdir*dgdx*x+v4*fdir+v4*(v1/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydy+2*fdir*dgdy*y+v4*fdir+v4*(v2/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdzdz+2*fdir*dgdz*z+v4*fdir+v4*(v3/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdy+y*fdir*dgdx+x*fdir*dgdy+v4*(x*y/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdz+z*fdir*dgdx+x*fdir*dgdz+v4*(x*z/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydz+z*fdir*dgdy+y*fdir*dgdz+v4*(y*z/v5)*(f2dir-fdir);
	  break;
        case isym_G6:
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;

	  v4=y*z*(3*v1-v2);
	  dgdx= 6*x*y*z;
	  dgdy= 3*(v1 - v2)*z;
	  dgdz= (3*v1 - v2)*y;
	  d2gdxdx=6*y*z;
	  d2gdydy=-6*y*z;
	  d2gdzdz=0;
	  d2gdxdy=6*x*z;
	  d2gdxdz=6*x*y;
	  d2gdydz=3*(v1 - v2);
	  
	  h=fdir*v4;
          *(ptr++)=v4*func;
	  *(ptr++)=func*dgdx+h*x;
	  *(ptr++)=func*dgdy+h*y;
	  *(ptr++)=func*dgdz+h*z;
	  *(ptr++)=func*d2gdxdx+2*fdir*dgdx*x+v4*fdir+v4*(v1/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydy+2*fdir*dgdy*y+v4*fdir+v4*(v2/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdzdz+2*fdir*dgdz*z+v4*fdir+v4*(v3/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdy+y*fdir*dgdx+x*fdir*dgdy+v4*(x*y/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdz+z*fdir*dgdx+x*fdir*dgdz+v4*(x*z/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydz+z*fdir*dgdy+y*fdir*dgdz+v4*(y*z/v5)*(f2dir-fdir);
	  break;
	case isym_G7:
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;

	  v4=(v1-3.*v2)*v1-(3.*v1-v2)*v2;
	  dgdx= 2*x*(2.*v1 - 6*v2);
	  dgdy= 2*y*(-6*v1 + 2*v2);
	  dgdz=0.0;
	  d2gdxdx=12*(v1 - v2);
	  d2gdydy=-12*(v1 - v2);
	  d2gdzdz=0;
	  d2gdxdy=-24*x*y;
	  d2gdxdz=0;
	  d2gdydz=0;
	  
	  h=fdir*v4;
          *(ptr++)=v4*func;
	  *(ptr++)=func*dgdx+h*x;
	  *(ptr++)=func*dgdy+h*y;
	  *(ptr++)=func*dgdz+h*z;
	  *(ptr++)=func*d2gdxdx+2*fdir*dgdx*x+v4*fdir+v4*(v1/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydy+2*fdir*dgdy*y+v4*fdir+v4*(v2/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdzdz+2*fdir*dgdz*z+v4*fdir+v4*(v3/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdy+y*fdir*dgdx+x*fdir*dgdy+v4*(x*y/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdz+z*fdir*dgdx+x*fdir*dgdz+v4*(x*z/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydz+z*fdir*dgdy+y*fdir*dgdz+v4*(y*z/v5)*(f2dir-fdir);
	  break;
	case isym_G8:
	  v1=x*x;
	  v2=y*y;
	  v3=z*z;
	  v5=v1+v2+v3;

	  v4=x*y*(v1-v2);
	  dgdx= (3*v1 - v2)*y;
	  dgdy= (v1 - 3*v2)*x;
	  dgdz= 0.0;
	  d2gdxdx=6*x*y;
	  d2gdydy=-6*x*y;
	  d2gdzdz=0;
	  d2gdxdy=3*(v1 - v2);
	  d2gdxdz=0;
	  d2gdydz=0;
	  
	  h=fdir*v4;
          *(ptr++)=v4*func;
	  *(ptr++)=func*dgdx+h*x;
	  *(ptr++)=func*dgdy+h*y;
	  *(ptr++)=func*dgdz+h*z;
	  *(ptr++)=func*d2gdxdx+2*fdir*dgdx*x+v4*fdir+v4*(v1/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydy+2*fdir*dgdy*y+v4*fdir+v4*(v2/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdzdz+2*fdir*dgdz*z+v4*fdir+v4*(v3/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdy+y*fdir*dgdx+x*fdir*dgdy+v4*(x*y/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdxdz+z*fdir*dgdx+x*fdir*dgdz+v4*(x*z/v5)*(f2dir-fdir);
	  *(ptr++)=func*d2gdydz+z*fdir*dgdy+y*fdir*dgdz+v4*(y*z/v5)*(f2dir-fdir);
          break;
        default:
          error("Bad symmetry in Cubic_spline::calcHessian.\n");
        }

      }
    }
  }
}
