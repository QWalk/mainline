
/*   macopt library          release 1.1          gradient-based optimizer

     Copyright   (c) 2002   David J.C. MacKay and Steve Waterhouse

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    GNU licenses are here :
    http://www.gnu.org/licenses/licenses.html

    Author contact details are here :
    http://www.inference.phy.cam.ac.uk/mackay/c/macopt.html       mackay@mrao.cam.ac.uk
*/
#ifndef MACOPT_H
#define MACOPT_H

class Macopt 
{
public:

  Macopt(int _n, 
	 int _verbose = 0, 
	 double _tolerance = 0.001, 
	 int _itmax = 100,
	 int _rich = 1);
  virtual ~Macopt();
  void macoptII(double *p, int    dim ); 
  void maccheckgrad(double *p,
		    int    n,
		    double epsilon, 
		    int stopat);	
protected:
  virtual double func(double* _p) = 0;
  virtual double dfunc(double* _p, double* _xi)  = 0; //returns the value..
  virtual void iteration_print(double val, double gg, double tol,  int it);
  int a_n;     /* dimension of parameter space */
private:
  double maclinminII(double *p);
  double macprodII (double * , double * , double ) ;
  void macopt_restart ( int ) ;

private:

  double a_tol ;    /* convergence declared when the gradient vector is smaller
		     in magnitude than this, or when the mean absolute 
		     step is less than this (see above) */
  double a_grad_tol_tiny ; /* if gradient is less than this, we definitely 
			    stop, even if we are not using a gradient 
			    tolerance */
  double a_step_tol_tiny ; /* if step is less than this, we stop, even if 
			    we are not using a step tolerance */
  int a_end_if_small_step ; /* defines the role of tol -- alternative is
			     end_on_small_grad */
  int a_its ;               /* number of its */
  int a_itmax ;             /* max */
  int a_rich ; /* whether to do the extra gradient evaluation at the beginning 
	      of each new line min */
  int a_verbose ; 
  double a_stepmax ;        /* largest step permitted (not used in macopt) */

  int a_linmin_maxits ;     /* in maclinmin */
  double a_linmin_g1 ;      /* factors for growing and shrinking the interval */
  double a_linmin_g2 ;
  double a_linmin_g3 ;
  double a_lastx     ;      /* keeps track of typical step length */
  double a_lastx_default ;  /* if maclinmin is reset, lastx is set to this */

/* These should not be touched by the user. They are handy pointers for macopt
   to use 
*/
  double a_gtyp ; /* stores the rms gradient for linmin */
  double *a_pt , *a_gx , *a_gy , *a_gunused ;
  double *a_xi , *a_g , *a_h ; 

  int a_restart ;           /* whether to restart macopt - fresh cg directions */

/* a_lastx :--- 1.0 might make general sense, (cf N.R.)
				  but the best setting of all is to have 
				  a prior idea of the eigenvalues. If 
				  the objective function is equal to sum of N
				  terms then set this to 1/N, for example 
				  Err on the small side to be conservative. */
};

#endif

