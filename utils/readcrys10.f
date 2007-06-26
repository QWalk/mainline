!    Copyright (C) 2007 Lucas K. Wagner
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.
!    
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!   
!    You should have received a copy of the GNU General Public License along
!    with this program; if not, write to the Free Software Foundation, Inc.,
!    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!    Read the fortran unit 10 file that crystal puts out
!    and write to stdout the MO's.  size should be the number
!    of functions per molecular orbital, and this will not work
!    if you used symmetry adapted functions(NOSYMADA must have
!    been set for Crystal version 98 and higher) 
      program read10
      implicit double precision (A-H, O-Z)
      integer kpoint
      parameter (size=104)
      DIMENSION A(size*size)
      iniz=0
      istart=0
      nspin=1
      read(*,*) kpoint
      do k=1, kpoint-1
         read(10)A
      enddo
      do spin=1,nspin
         write(*,*) "---------Spin ", spin
         read(10)A
         do i=1,size
            write(*,*) "====================="
            write(*,*) "CO ", i
            do j=1,size
               write(*,*) j,A((istart+i-1)*size+j)
            enddo
         enddo
         istart=istart+size
      enddo
      end
