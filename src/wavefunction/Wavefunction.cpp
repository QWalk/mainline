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

#include "Wavefunction.h"
#include "Program_options.h"
#include "MatrixAlgebra.h"

//----------------------------------------------------------------------

void extend_parm_deriv(Parm_deriv_return & ret1, const Parm_deriv_return & ret2) { 
  int nparms=ret1.gradient.GetDim(0)+ret2.gradient.GetDim(0);
  int nparms1=ret1.gradient.GetDim(0);
  int nparms2=ret2.gradient.GetDim(0);
  int nelectrons=ret1.val_gradient.GetDim(0);
  Parm_deriv_return derivatives;
  if(nparms1==0 and nparms2==0) return;
  if(nparms2==0) { 
    for(int e=0; e < nelectrons; e++) {
      for(int d=0;d < 3; d++)  { 
        ret1.val_gradient(e,d)+=ret2.val_gradient(e,d);
      }
    }
    return;
  }
  derivatives.need_hessian=ret1.need_hessian;
  derivatives.gradient.Resize(nparms);
  derivatives.hessian.Resize(nparms,nparms);
  for(int i=0; i< nparms1; i++) { 
    derivatives.gradient(i)=ret1.gradient(i);
    derivatives.hessian(i,i)=ret1.hessian(i,i);
  }
  for(int i=nparms1; i< nparms1+nparms2; i++) { 
    derivatives.gradient(i)=ret2.gradient(i-nparms1);
    derivatives.hessian(i,i)=ret2.hessian(i-nparms1,i-nparms1);
  }
  for(int i=0; i< nparms1; i++) { 
    for(int j=i+1; j< nparms1; j++) { 
      derivatives.hessian(i,j)=derivatives.hessian(j,i)=ret1.hessian(i,j);
    }
  }//cout << "t " << endl;
  for(int i=0; i< nparms1; i++) { 
    for(int j=nparms1; j< nparms1+nparms2; j++) { 
      derivatives.hessian(i,j)=derivatives.hessian(j,i)=ret1.gradient(i)*ret2.gradient(j-nparms1);
    }
  }//cout << " q " << endl;
  for(int i=nparms1; i< nparms1+nparms2; i++) { 
    for(int j=i+1; j< nparms1+nparms2; j++) { 
      derivatives.hessian(i,j)=derivatives.hessian(j,i)=ret2.hessian(i-nparms1,j-nparms1);
    }
  }


  derivatives.val_gradient.Resize(nelectrons,3);
  for(int e=0; e< nelectrons; e++) { 
    for(int d=0;d < 3; d++) {
      derivatives.val_gradient(e,d)=ret1.val_gradient(e,d)+ret2.val_gradient(e,d);
    }
  }
  
  //  Extending the derivative of the laplacian
  derivatives.gradderiv.Resize(nparms,nelectrons,4);

  for(int p=0; p < nparms1; p++) { 
    for(int e=0; e< nelectrons; e++) { 
      for(int d=0; d< 4; d++) { 
        derivatives.gradderiv(p,e,d)=ret1.gradderiv(p,e,d);
      }
    }
  }
  for(int p=0; p < nparms2; p++) { 
    for(int e=0; e< nelectrons; e++) { 
      for(int d=0; d< 4; d++) { 
        derivatives.gradderiv(p+nparms1,e,d)=ret2.gradderiv(p,e,d);
      }
    }
  }
  //Now do the cross-terms
  for(int p=0; p < nparms1; p++) { 
    for(int e=0; e< nelectrons; e++) { 
      doublevar dot=0;
      for(int d=0; d< 3; d++) { 
        dot+=ret1.gradderiv(p,e,d)*ret2.val_gradient(e,d);
      }
      derivatives.gradderiv(p,e,3)+=2*dot;
    }
  }

  for(int p=0; p < nparms2; p++) { 
    for(int e=0; e< nelectrons; e++) { 
      doublevar dot=0;
      for(int d=0; d< 3; d++) { 
        dot+=ret2.gradderiv(p,e,d)*ret1.val_gradient(e,d);
      }
      derivatives.gradderiv(p+nparms1,e,3)+=2*dot;
    }
  }

  ret1=derivatives;
}
//----------------------------------------------------------------------

int deallocate(Wavefunction * & wfptr)
{
  if(wfptr == NULL)
    return 0;

  delete wfptr;
  wfptr=NULL;
  return 1;
}


//------------------------------------------------------------------------
