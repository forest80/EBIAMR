
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PROJECTION_F.H>
#include <AMReX_ArrayLim.H>

#define SDIM 3

       subroutine FORT_ACCEL_TO_VEL( lo, hi, 
     &     uold,DIMS(uold),
     &     dt,
     &     unew,DIMS(unew) )
c
c     This function converts unew into a velocity via
c     Unew = Uold + alpha*Unew
c
       implicit none
       integer    lo(SDIM), hi(SDIM)
       REAL_T     dt
       integer    DIMDEC(uold),DIMDEC(unew)
       REAL_T     uold(DIMV(uold),SDIM)
       REAL_T     unew(DIMV(unew),SDIM)

       integer i, j, k, n

       do n = 1, SDIM
          do k = lo(3)-1, hi(3)+1
             do j = lo(2)-1, hi(2)+1
                do i = lo(1)-1, hi(1)+1
                   unew(i,j,k,n) = uold(i,j,k,n) + dt*unew(i,j,k,n)
                end do
             end do
          end do
       end do

       end

      subroutine FORT_VEL_TO_ACCEL( lo, hi, 
     &     unew,DIMS(unew),
     &     uold,DIMS(uold),
     &     dt )
c     
c     This function converts unew into an acceleration
c
      implicit none
      integer    lo(SDIM), hi(SDIM)
      REAL_T     dt
      integer    DIMDEC(unew),DIMDEC(uold)
      REAL_T     uold(DIMV(uold),SDIM)
      REAL_T     unew(DIMV(unew),SDIM)

      integer i, j, k, n

      do n = 1, SDIM
         do k = lo(3)-1, hi(3)+1
            do j = lo(2)-1, hi(2)+1
               do i = lo(1)-1, hi(1)+1
                  unew(i,j,k,n) = (unew(i,j,k,n)-uold(i,j,k,n))/dt
               end do
            end do
         end do
      end do

      end

      subroutine FORT_PROJ_UPDATE(
     &     boxlo, boxhi, nvar, ngrow,
     &     un, DIMS(un),
     &     alpha,
     &     uo, DIMS(uo) )
c     
c     This function updates un via un = un + alpha*uo
c     The loop bounds are determined in the C++
c
      implicit none
      integer    boxlo(SDIM), boxhi(SDIM), nvar, ngrow
      REAL_T     alpha
      integer    DIMDEC(un),DIMDEC(uo)
      REAL_T     un(DIMV(un),nvar)
      REAL_T     uo(DIMV(uo),nvar)

      integer i, j, k, n

!$omp parallel do private(i,j,k,n) collapse(2)
      do n = 1, nvar
         do k = boxlo(3), boxhi(3)
            do j = boxlo(2), boxhi(2)
               do i = boxlo(1), boxhi(1)
                  un(i,j,k,n) = un(i,j,k,n) + alpha*uo(i,j,k,n)
               end do
            end do
         end do
      end do
!$omp end parallel do

      end

      subroutine FORT_ANELCOEFFMPY(a,DIMS(grid),domlo,domhi,ng,
     $                             anel_coeff,nr,bogus_value,mult)
c 
c     multiply A by anel_coeff
c
c 
c     NOTE: THIS ROUTINE HAS BEEN MODIFIED SO THAT ALL VALUES
c           OUTSIDE THE DOMAIN ARE SET TO BOGUS VALUE
c
      implicit none
      integer    ng,nr
      integer    DIMDEC(grid)
      integer    domlo(3), domhi(3)
      REAL_T     a(ARG_L1(grid)-ng:ARG_H1(grid)+ng, 
     $             ARG_L2(grid)-ng:ARG_H2(grid)+ng,
     $             ARG_L3(grid)-ng:ARG_H3(grid)+ng)
      REAL_T     anel_coeff(ARG_L3(grid)-nr:ARG_H3(grid)+nr)
      REAL_T     bogus_value
      integer    mult

      integer i, j, k

      if (mult .eq. 1) then
         do k = ARG_L3(grid)-ng, ARG_H3(grid)+ng
         do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do i = ARG_L1(grid)-ng, ARG_H1(grid)+ng
           a(i,j,k) = a(i,j,k) * anel_coeff(k)
         end do
         end do
         end do
      else if (mult .eq. 0) then
         do k = ARG_L3(grid)-ng, ARG_H3(grid)+ng
         do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do i = ARG_L1(grid)-ng, ARG_H1(grid)+ng
           a(i,j,k) = a(i,j,k) / anel_coeff(k)
         end do
         end do
         end do
      else 
         print *,'BOGUS MULT IN ANELCOEFFMULT ',mult
         stop
      end if

      if (ARG_L3(grid)-ng .lt. domlo(3)) then
         do i = ARG_L1(grid)-ng, ARG_H1(grid)+ng
         do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do k = ARG_L3(grid)-ng, domlo(3)-1
           a(i,j,k) = bogus_value
         end do
         end do
         end do
      end if

      if (ARG_H3(grid)+ng .gt. domhi(3)) then
         do i = ARG_L1(grid)-ng, ARG_H1(grid)+ng
         do j = ARG_L2(grid)-ng, ARG_H2(grid)+ng
         do k = domhi(3)+1, ARG_H3(grid)+ng
           a(i,j,k) = bogus_value
         end do
         end do
         end do
      end if

      end
