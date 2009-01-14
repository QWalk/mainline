#ifndef GESQUA_H_INCLUDED
#define GESQUA_H_INCLUDED

/*!
Computes the integration points for a quadrature for spherical harmonics
*/
void gesqua(int & nq, Array1 <doublevar> & xq,Array1 <doublevar> & yq,
            Array1 <doublevar> & zq, Array1 <doublevar> & wq);

#endif  // GESQUA_H_INCLUDED
