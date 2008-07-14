c     from TETPACK Library
c     http://perso.fundp.ac.be/~jpvigner/homepage/tetpack/tetpack.html
c
c
      subroutine setk01(na,a,ptk,nptk,idef,ntet,nkmax,ntmax)
c     ======================================================================
c     set the k-points in the irreductible wedge of the brillouin zone
c     for a simple cubic lattice with direct lattice parameter a,
c     symmetry is oh
c
c     Sampling volume = 1/48 irreductible wedge of the first Brillouin zone.
c     Cubic mesh of points with step size pi/a/NA along the three coordinate
c     axes.
c
c        CALL SETK01(NA,A,PTK,NPTK,IDEF,NTET,NKMAX,NTMAX)
c
c     with DIMENSION PTK(4,NKMAX),IDEF(5,NTMAX)
c     NA is the k-mesh discretization parameter,
c     A is the direct-space cubic lattice parameter "a".
c
c     Output:
c
c       - the number NPTK of k-points  
c           NPTK = (NA+1)*(NA+2)*(NA+3)/6
c       - the number NTET of tetrahedra
c           NTET = NA**3
c       - the cartesian x, y, z components of the k-points in PTK(1,K),
c         PTK(2,K), PTK(3,K) and the corresponding sampling weigths in
c         PTK(4,K) , K = 1 ... NPTK
c       - a table defining the tetrahedra (list of corner points and the
c         volumes of the tetrahedra) in a two-dimensional array IDEF
c
c     PTK is a REAL*8 array and IDEF is an INTEGER*4 array with dimensions
c     PTK(4,NKMAX) and IDEF(5,NTMAX) in the calling program. Execution is
c     bypassed whenever NPTK or NTET exceed NKMAX or NTMAX. 
c     ======================================================================
      implicit real*8(a-h,o-z)
      real*4 avol
      dimension ptk(4,nkmax),idef(5,ntmax)
      equivalence (ivol,avol)
      pi = 3.141592653589793238d0
      if(na.le.0) goto 97
      if(a.le.0.0d0) goto 98
      nptk = (na+1)*(na+2)*(na+3)/6
      if(nptk.gt.nkmax) stop '***  nptk exceeds nkmax ***'
      ntet = na**3
      if(ntet.gt.ntmax) stop '***  ntet exceeds ntmax ***'
c *** set the k-points
      dk=pi/a/na
      n=na
      write(6,100) nptk,ntet,n*dk,n*dk,n*dk
      w = 6.0d0/na**3
      nptk=0
      i=0
    1 j=0
    2 k=0
    3 nptk=nptk+1
c          nptk = i*(i+1)*(i+2)/6 + j*(j+1)/2 + k + 1
      wk = w
      if(i.eq.j .or. j.eq.k .or. k.eq.i) wk = wk/2.0d0
      if(i.eq.j .and. j.eq.k) wk = wk/3.0d0
      if(k.eq.0) wk = wk/2.0d0
      if(j.eq.0) wk = wk/2.0d0
      if(j.eq.na) wk = wk/2.0d0
      if(i.eq.na) wk = wk/2.0d0
      ptk(1,nptk)=dfloat(i)*dk
      ptk(2,nptk)=dfloat(j)*dk
      ptk(3,nptk)=dfloat(k)*dk
      ptk(4,nptk) = wk
      k=k+1
      if(k.le.j) goto 3
      j=j+1
      if(j.le.i) goto 2
      i=i+1
      if(i.le.n) goto 1
      ptk(4,1) = w/48.0d0
      ptk(4,nptk) = w/48.0d0
c *** define the tetrahedra
      ntet=0
      i7=0
      i=0
    4 ix=(i+1)*(i+2)/2
      j=0
    5 k=0
    6 i7=i7+1
      i6=i7+ix
      i2=i6+j+1
      i1=i2+1
      i8=i7+1
      i5=i6+1
      i3=i7+j+1
      i4=i3+1
      ntet=ntet+1
      idef(1,ntet)=i7
      idef(2,ntet)=i6
      idef(3,ntet)=i2
      idef(4,ntet)=i1
      if(k.eq.j) goto 7
      ntet=ntet+1
      idef(1,ntet)=i7
      idef(2,ntet)=i6
      idef(3,ntet)=i5
      idef(4,ntet)=i1
      ntet=ntet+1
      idef(1,ntet)=i7
      idef(2,ntet)=i8
      idef(3,ntet)=i5
      idef(4,ntet)=i1
    7 if(j.eq.i) goto 8
      ntet=ntet+1
      idef(1,ntet)=i7
      idef(2,ntet)=i3
      idef(3,ntet)=i2
      idef(4,ntet)=i1
      ntet=ntet+1
      idef(1,ntet)=i7
      idef(2,ntet)=i3
      idef(3,ntet)=i4
      idef(4,ntet)=i1
      if(k.eq.j) goto 8
      ntet=ntet+1
      idef(1,ntet)=i7
      idef(2,ntet)=i8
      idef(3,ntet)=i4
      idef(4,ntet)=i1
    8 k=k+1
      if(k.le.j) goto 6
      j=j+1
      if(j.le.i) goto 5
      i=i+1
      if(i.lt.n) goto 4
      avol=1.0d0/dfloat(ntet)
      do 15 it=1,ntet
   15 idef(5,it)=ivol
      return
   97 write(6,101)
      goto 99
   98 write(6,102)
   99 stop
  100 format(' sampling the 48th part of the b.z. of the',
     ,' simple cubic lattice'/1x,i5,' k-points',i7,' tetrahedra'/
     .' kxmax =',d11.4,'  kymax =',d11.4,'  kzmax =',d11.4)
  101 format(' ***  na is not a positive integer ***')
  102 format(' ***  a is not a positive number ***')
      end
