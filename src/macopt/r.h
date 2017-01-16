
/*   r.h         header for assorted minor subroutines library          release 1.1    

     Copyright   (c) 2002   David J.C. MacKay

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
/*	ANSI Version 1 	*/
/* 9 6 92 */

#ifdef __cplusplus
extern "C" {
#endif

#include	<stdio.h>
#include	<math.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include 	"nrutil.h"

#ifndef INT
#define INT( a ) ( ( ( a ) > 0.0 ) ? ( (int) ( (a)+0.5 ) ) : ( (int) ( (a)-0.5 ) ) ) 
#endif
#ifndef MIN
#define MIN( a , b ) ( ( a ) < ( b ) ? ( a ) : ( b ) ) 
#endif
#ifndef MAX
#define MAX( a , b ) ( ( a ) > ( b ) ? ( a ) : ( b ) ) 
#endif
#ifndef PI                    
#define	PI	3.1415926535
#endif
#ifndef STPI
#define	STPI	2.50663 		/* sqrt ( 2*PI ) */
#endif
#ifndef LTPI
#define	LTPI	1.83788 		/* log ( 2*PI ) */
#endif
#ifndef FLUSH
#define	FLUSH	fflush( stdout ) ; 
#endif
#define fnewline fprintf( fp , "\n" )
#define junkstring fscanf( fp , "%s" , junk )
#define FNEWLINE fprintf( fp , "\n" ) 
#define JUNKSTRING fscanf( fp , "%s" , junk )

/* list of functions that should be accepted from maths library */

/*
void   srandom ( long ) ;
int    random ( void ) ; 

double	drand48( ) ;
void 	srand48( long ) ;
int 	abs( int ) ; 
double 	fabs( double ) ; 
double	exp( double ) ;
double	sqrt( double ) ; 
*/

#ifndef BELL
#define BELL '\007'
#endif
#ifndef ALERT
#define ALERT '\a' /* these are integers with the relevant values */
#endif

/* routines in ansi/r.c */

typedef struct {
  unsigned char **m ; /* the matrix (this copy gets munged */
  unsigned char **mo; /* original copy -- unchanged, except by the 
		       modify_cmatrix_row function */
  unsigned char **mi; /* stores the inverse */
  unsigned char **mt; /* used to assist in undoing permutation */
  int *perm ;  /* the permutation */
  int *iperm ;  /* and its inverse */
  int i ;       /* what row we have got up to in the LU decomposition */
  int l , N ;   /* e.g. 1 , N */
} cm_inversion ;


double ***fancy_dmatrix3( int , int , int , int , int * , int * ) ;
double **fancy_dmatrix2( int , int , int * , int * ) ;
void enterdmatrix ( double ** , int * , int * , int * , int * ) ;
double 	random_1( void ) ;
void 	randomise( long ) ;
void randomdmatrix( double ** , int , int , int , int , int , double ) ;
void IJdmatrix( double ** , int , int , int , int , double , double ) ;
void constantdmatrix( double ** , int , int , int , int , double ) ;
double readindmatrix ( double ** , int , int , int , int , char * ) ;
int readinimatrix ( int ** , int , int , int , int , char * ) ;
int fread_imatrix ( int ** , int , int , int , int , FILE * ) ;
void readinlumatrix ( double ** , int * , int , char * ) ;
void readindvector ( double * , int , int , char * ) ;
int readdvector ( double * , int , int , char * ) ;
int writedvector ( double * , int , int , char * ) ;
void readinivector ( int * , int , int , char * ) ;
int fread_ivector ( int * , int , int , FILE * ) ;
int fread_cvector ( unsigned char * , int , int , FILE * ) ;
int fread_dvector ( double * , int , int , FILE * ) ;
void inputf( float * ) ;
void inputd( double * ) ;
void inputi( int * ) ;
void inputc( unsigned char * ) ;
void inputrc( unsigned char * ) ;
void inputrf( float * ) ;
void inputrd( double * ) ;
void inputri( int * ) ;
void clearscan( void ) ;
void typeindvector( double * , int , int ) ;
void set_dvector_const( double * , int , int , double ) ;
void set_dvector_c_dvector ( double * , int , int , double , double * ) ;
void set_ivector_const( int * , int , int , int ) ;
void typeindmatrix ( double ** , int , int , int , int ) ;
void typeincmatrix ( unsigned char ** , int , int , int , int ) ;
int  ***imatrix3( int  , int , int , int , int , int ) ;
long int ***limatrix3( int , int , int , int , int , int ) ;
long int **limatrix( int , int , int , int ) ;
int ipower( int , int ) ;

double ***dmatrix3( int   , int , int , int , int , int ) ;
void printoutimatrix ( int **, int , int , int , int ) ;
void write_imatrix ( FILE * , int **, int , int , int , int ) ;
void write_imatrix2 ( FILE * , int **, int , int , int , int ) ;
void cmatrix2pbm ( unsigned char **, int , int , int , int , FILE *) ;
void print_cm_inversion ( cm_inversion *) ;
void free_cm_inversion ( cm_inversion *) ;
void undo_cm_perm ( cm_inversion *) ;
int invert_cmatrix( cm_inversion *) ;
int modify_cmatrix_row ( cm_inversion * ) ;
int invert_utriangularc ( cm_inversion * ) ;
void allocate_cm_inversion ( unsigned char **, int , int , unsigned char **,
			    cm_inversion *) ;
void printoutcmatrix ( unsigned char **, int ,int , int , int ) ;
void printoutcmatrix1 ( unsigned char **, int ,int , int , int ) ;
void printoutivector ( int * , int , int ) ;
void write_ivector ( FILE * , int * , int , int ) ;
void write_cvector ( FILE * , unsigned char * , int , int ) ;
void printoutcvector ( unsigned char *, int , int ) ;
void printoutcvector1 ( unsigned char *, int , int ) ;
void printoutdmatrix ( double **, int ,int , int , int , int ) ;
int writedmatrix ( double **, int ,int , int , int , char * ) ;
void printoutdmatrix3 ( double ***, int ,int ,int ,int , int , int , int ) ;
void pd( double , int ) ;
void pdv( double * , int, int, int ) ;
void 	pause_for_return( void ) ;
double gammln ( double ) ;
int find_rank ( double , double * , int , int ) ;
int dotprod_mod2 ( int * , int * , int  , int  ) ;
int idotprod ( int * , int * , int  , int  ) ;
unsigned char cdotprod_mod2 ( unsigned char * , unsigned char * , int  , int  ) ;
void mult_cms ( unsigned char ** , unsigned char ** , unsigned char ** , 
	       int , int  ) ;
void mult_cm_cv ( unsigned char ** , unsigned char * , unsigned char * , 
	       int , int  ) ;
int cdotprod ( unsigned char * , unsigned char * , int  , int  ) ;
float ran3(int*) ;

#ifdef __cplusplus
}
#endif
