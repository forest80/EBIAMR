
#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_ArrayLim.H>
#include <infl_frc.H>
#include <FLUCTFILE.H>
#include <AMReX_BC_TYPES.H>

#define SDIM BL_SPACEDIM

#define FF_UNIT       20


c ::: -----------------------------------------------------------
c ::: This routine is used by INITDATA and the fill routines to 
c ::: extrapolate the perturbations from the flct_file to 
c ::: fill the data required for forcing the inflow.  Mostly this 
c ::: routine manages the reading of the data from the flct_file
c ::: and then passes off to the XTR_DAT routine to actually
c ::: extrapolate the data and fill the arrays.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: fillComp       =>  Component to fill
c ::: DIMS(inflDat)  =>  Dimensions of inflDat
c ::: inflDat       <=   Array to fill
c ::: dx             =>  Grid spacing
c ::: time           =>  Time for the fill
c ::: -----------------------------------------------------------
c
      subroutine INFL_FILL_PERIODIC(fillComp, DIMS(inflDat), inflDat,
     &     xlo, dx, time, bc, probLo, probHi)

      implicit none

      integer fillComp, DIMDEC(inflDat), bc(SDIM,2)

      REAL_T time, xlo(SDIM), dx(SDIM)
      REAL_T inflDat(DIMV(inflDat)), probLo(SDIM), probHi(SDIM)

      integer lo(SDIM), hi(SDIM), nCmpFile
      integer n, npass, filePnt, storePnt, filePntOld, npassOld, proc
      REAL_T dtFile, timeMaxFile, timeOffset, timePnt
c
c     The arrays below are dimensioned as 3-d no matter what the 
c     BL_SPACEDIM is.  This is to allow for the fact that the inflow forcing
c     data arrays are always 3-d.
c
      logical full_file_in_memory
      integer dimFile(3), loStoreDim(3), hiStoreDim(3)
      integer FF_DIMDEC(storDat), ierr
      REAL_T dxFile(3), probSizeFile(3), xloFile(3), xhiFile(3)
      REAL_T, allocatable :: storDat(:,:,:,:)

      save full_file_in_memory, loStoreDim, hiStoreDim, storDat,
     &     dimFile, xloFile, xhiFile, dxFile, dtFile, timeMaxFile, 
     &     storePnt, filePnt, timePnt, timeOffset, npass

#include <INFL_FORCE_F.H>

      call bl_pd_myproc(proc)

      if (.NOT. ALLOCATED(storDat)) then

         ierr = 0

         open(FF_UNIT, file=trim(flct_file)//'/HDR', form='formatted', action='read', status='old', iostat=ierr)

         if (ierr .ne. 0) then
            call bl_abort('Problem opening file: ' // trim(flct_file) // '/HDR')
         end if

         call RD_SCL_FLCTHD(FF_UNIT,nCmpFile,dimFile,probSizeFile,dxFile)

         close(FF_UNIT)

         call SET_LOHI(DIMS(inflDat), lo, hi)

         xloFile(3) = zero
         xhiFile(3) = zero
         do n = 1, SDIM
            xloFile(n) = half * (probHi(n) + probLo(n)) - half * probSizeFile(n)
            xhiFile(n) = xloFile(n) + probSizeFile(n)
         enddo

         do n = 1, 3
            loStoreDim(n) = 1
            hiStoreDim(n) = dimFile(n)
         enddo

         if (numInflPlanesStore .GT. 0 .AND. numInflPlanesStore .LT. hiStoreDim(strmwse_dir)) then
             hiStoreDim(strmwse_dir) = numInflPlanesStore
             full_file_in_memory = .false.
          else
             full_file_in_memory = .true.
          endif

         if (full_file_in_memory) loStoreDim(strmwse_dir) = 0

         call FF_SET_ARGS(FF_DIMS(storDat), loStoreDim, hiStoreDim)
         ALLOCATE(storDat(FF_DIMV(storDat),nCompInflow))
c
c       ::::   Convert the streamwise direction lengths to times   ::::
c       :::: and set up the pointers into the data arrays and file ::::
c
         timeMaxFile = probSizeFile(strmwse_dir) / convVel
         dtFile      = dxFile(strmwse_dir) / convVel
         npass       = INT(time / timeMaxFile)
         timeOffset  = DBLE(npass) * timeMaxFile
         filePnt     = INT((time - timeOffset) / dtFile) + 1
         timePnt     = timeOffset + DBLE(filePnt-1) * dtFile

         if (full_file_in_memory) then
            storePnt = filePnt + 1
         else
            storePnt = 2
         endif

         if (full_file_in_memory) then
            call FILL_FRCARRYS(1, 1, dimFile, nCompInflow,
     &           FF_DIMS(storDat), storDat)
            if (strmwse_dir.eq.2) then
               do n=1,nCompInflow
                  storDat(:,0,:,n) = storDat(:,hiStoreDim(2)-1,:,n)
               enddo        
            else if (strmwse_dir.eq.3) then
               do n=1,nCompInflow
                  storDat(:,:,0,n) = storDat(:,:,hiStoreDim(3)-1,n)
               end do
            else
               call bl_abort('Reflections to X direction not coded')
            endif
         else
            call FILL_FRCARRYS(filePnt - 1, 1, dimFile, nCompInflow,
     &           FF_DIMS(storDat), storDat)
         endif
      endif

      call FF_SET_ARGS(FF_DIMS(storDat), loStoreDim, hiStoreDim)

      npassOld   = npass
      filePntOld = filePnt

      do while ( (time .lt. timePnt - half*dtFile) .or. (time .gt. timePnt + half*dtFile) )
         if (time .gt. timePnt + half*dtFile) then
            filePnt  = filePnt + 1
            storePnt = storePnt + 1
            if (filePnt .gt. dimFile(strmwse_dir) - 1) then
               npass      = npass + 1
               timeOffset = DBLE(npass) * timeMaxFile
               filePnt    = filePnt - dimFile(strmwse_dir) + 1
               if (full_file_in_memory) storePnt = filePnt + 1
            endif
         else if (time .lt. timePnt - half*dtFile) then
            filePnt  = filePnt - 1
            storePnt = storePnt - 1
            if (filePnt .lt. 1) then
               npass      = npass - 1
               timeOffset = DBLE(npass) * timeMaxFile
               filePnt    = filePnt + dimFile(strmwse_dir) - 1
               if (full_file_in_memory) storePnt = filePnt + 1
            endif
         endif
         timePnt = timeOffset + DBLE(filePnt-1) * dtFile
      enddo

      if (full_file_in_memory) then
         if ((filePnt + 1 .gt. hiStoreDim(strmwse_dir)) .or. (filePnt .lt. 1)) then
            print*,"Problem "
            print*,"time = ",time
            print*,"timeMaxFile = ",timeMaxFile
            print*,"dtFile = ",dtFile
            print*,"npass = ",npass
            print*,"timeOffset = ",timeOffset
            print*,"filePnt,timePnt,storePnt = ",filePnt,timePnt,storePnt
            call flush(6)
            call bl_abort('INFL_FILL_PERIODIC: how did we get here?')
         endif
      else
         if ((storePnt + 1 .gt. hiStoreDim(strmwse_dir)) .or. (storePnt - 1 .lt. 1)) then
            call FILL_FRCARRYS(filePnt - 1, 1, dimFile, nCompInflow,
     &           FF_DIMS(storDat), storDat)
            storePnt = 2
         endif
      endif
 
      if (full_file_in_memory) then
         !
         ! Should interpolate about filePnt in 2d not storePnt.
         !
         call INTRP_DATA(time, dx, filePnt, fillComp, storePnt,
     &        timePnt, dtFile, dxFile, xloFile, xhiFile, nCompInflow, 
     &        FF_DIMS(storDat), storDat, DIMS(inflDat), inflDat,
     &        bc, probLo, probHi)
      else
         call INTRP_DATA(time, dx, storePnt, fillComp, filePnt,
     &        timePnt, dtFile, dxFile, xloFile, xhiFile, nCompInflow, 
     &        FF_DIMS(storDat), storDat, DIMS(inflDat), inflDat,
     &        bc, probLo, probHi)
      endif

      end

      subroutine INFL_FILL_SWIRL(fillComp, DIMS(inflDat), inflDat, xlo,
     &     dx, time, bc, probLo, probHi)

      implicit none

      integer fillComp, DIMDEC(inflDat), bc(SDIM,2)
      REAL_T time, xlo(SDIM), dx(SDIM)
      REAL_T inflDat(DIMV(inflDat)), probLo(SDIM), probHi(SDIM)

      integer lo(SDIM), hi(SDIM), nCmpFile, storePnt, shft, proc
      integer n, filePnt, fileBase, FF_DIMDEC(storDat), ierr
      integer dimFile(3), loStoreDim(3), hiStoreDim(3)
      REAL_T dxFile(3), probSizeFile(3), xloFile(3), xhiFile(3), timePnt
      REAL_T storDat(:,:,:,:), dtFile

      allocatable storDat

      save loStoreDim, hiStoreDim, storDat, dimFile, xloFile, xhiFile,
     &     dxFile, fileBase

#include <INFL_FORCE_F.H>

      call bl_pd_myproc(proc)

      if (.NOT. ALLOCATED(storDat)) then

         ierr = 0

         open(FF_UNIT, file=trim(flct_file)//'/HDR', form='formatted', action='read', status='old', iostat=ierr)

         if (ierr .ne. 0) then
            call bl_abort('Problem opening file: ' // trim(flct_file) // '/HDR')
         end if

         call RD_SCL_FLCTHD(FF_UNIT,nCmpFile,dimFile,probSizeFile,dxFile)

         close(FF_UNIT)

         call SET_LOHI(DIMS(inflDat), lo, hi)

         xloFile(3) = zero
         xhiFile(3) = zero

         do n = 1, SDIM
            xloFile(n) = half * (probHi(n) + probLo(n)) - half * probSizeFile(n)
            xhiFile(n) = xloFile(n) + probSizeFile(n)
         enddo

         do n = 1, 3
            loStoreDim(n) = 1
            hiStoreDim(n) = dimFile(n)
         enddo

         if (numInflPlanesStore .GT. 0 .AND.
     &        numInflPlanesStore .LT. hiStoreDim(strmwse_dir)) then
            hiStoreDim(strmwse_dir) = numInflPlanesStore
         endif

         call FF_SET_ARGS(FF_DIMS(storDat), loStoreDim, hiStoreDim)

         ALLOCATE(storDat(FF_DIMV(storDat),nCompInflow))

         if ((time .le. fluct_times(1)) .or.
     &        (time .ge. fluct_times(dimFile(strmwse_dir)))) then
            write(6,101) time, fluct_times(1), fluct_times(dimFile(strmwse_dir))
            call bl_abort('INFL_FILL_SWIRL: time is out of range')
         endif

         filePnt = 2
         do while (fluct_times(filePnt) .lt. time)
            filePnt = filePnt + 1
         enddo

         if (filePnt .eq. dimFile(strmwse_dir)) then
            filePnt = filePnt - 1
         endif

         shft = filePnt + hiStoreDim(strmwse_dir) - 2 - dimFile(strmwse_dir)

         if (shft .gt. 0) then
            fileBase = filePnt - shft - 1
         else
            fileBase = filePnt - 1
         endif

         call FILL_FRCARRYS_SWIRL(fileBase, 1, nCompInflow,
     &        FF_DIMS(storDat), storDat)

      endif

      call FF_SET_ARGS(FF_DIMS(storDat), loStoreDim, hiStoreDim)

      if ((time .le. fluct_times(1)) .or.
     &     (time .ge. fluct_times(dimFile(strmwse_dir)))) then
            write(6,101) time, fluct_times(1), fluct_times(dimFile(strmwse_dir))
         call bl_abort('INFL_FILL_SWIRL: time is out of range')
      endif

101   format('time: ',f10.5,' range: ',f10.5,' ... ',f10.5)

      filePnt = 2
      do while (fluct_times(filePnt) .lt. time)
         filePnt = filePnt + 1
      enddo

      if (filePnt .eq. dimFile(strmwse_dir)) then
         filePnt = filePnt - 1
      endif

      if ((filePnt .gt. fileBase) .and.
     &     (filePnt .lt. (fileBase + hiStoreDim(strmwse_dir) - 1))) then
         !
         ! We've got enough data in storDat to do the interpolation.
         !
         storePnt = 1 + (filePnt - fileBase)
      else
         !
         ! Got to reload storDat.
         ! Always load exactly hiStoreDim(strmwse_dir) platters of data.
         !
         shft = filePnt + hiStoreDim(strmwse_dir) - 2 - dimFile(strmwse_dir)

         if (shft .gt. 0) then
            storePnt = 2 + shft
            fileBase = filePnt - shft - 1
         else
            storePnt = 2
            fileBase = filePnt - 1
         endif

         call FILL_FRCARRYS_SWIRL(fileBase, 1, nCompInflow,
     &        FF_DIMS(storDat), storDat)
      endif

      timePnt = fluct_times(filePnt)
      dtFile  = -1 ! Not used

      call INTRP_DATA(time, dx, storePnt, fillComp, filePnt,
     &     timePnt, dtFile, dxFile, xloFile, xhiFile, nCompInflow, 
     &     FF_DIMS(storDat), storDat, DIMS(inflDat), inflDat,
     &     bc, probLo, probHi)

      end

      subroutine INFL_FILL(fillComp, DIMS(inflDat), inflDat, xlo, dx, time,
     &                     bc, probLo, probHi)

      implicit none

      integer fillComp, DIMDEC(inflDat), bc(SDIM,2)

      REAL_T time, xlo(SDIM), dx(SDIM)
      REAL_T inflDat(DIMV(inflDat)), probLo(SDIM), probHi(SDIM)

      REAL_T offset_time

#include <INFL_FORCE_F.H>

      if (infl_type .eq. -1) then
         call bl_abort('INFL_FILL: infl_type is not set')
      endif

      offset_time = time + tstart_turb

      if (infl_type .eq. infl_periodic_type) then
         call INFL_FILL_PERIODIC(fillComp, DIMS(inflDat), inflDat,
     &        xlo, dx, offset_time, bc, probLo, probHi)
      else if (infl_type .eq. infl_swirl_type) then
         call INFL_FILL_SWIRL(fillComp, DIMS(inflDat), inflDat,
     &        xlo, dx, offset_time, bc, probLo, probHi)
      else
         call bl_abort('INFL_FILL: unknown infl_type')
      endif
      end
c
c ::: -----------------------------------------------------------
c ::: This routine fills the inflow forcing data array from the file.
c ::: A basepoint is specified for the array as well as for the file.
c ::: These are the points in the strmwse_dir at which reading from the
c ::: file is started and at which the array is filled from.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: 
c ::: baseFilePnt   => Basepoint in the file in the strmwse_dir to 
c :::                    start reading from.
c ::: baseArrayPnt  => Basepoint in the array in the strmwse_dir to 
c :::                    start filling from.  The array is filled from
c :::                    this point in strmwse_dir to the end of the 
c :::                    array.
c ::: dimFile       => Dimensions from the header of the fluctuations
c ::: nComp         => Number of components in the array
c ::: FF_DIMS(dat)  => Dimensions of the array dat
c ::: dat          <=  Array to fill
c ::: -----------------------------------------------------------
c
      subroutine FILL_FRCARRYS(baseFilePnt, baseArrayPnt, dimFile, 
     &                         nComp, FF_DIMS(dat), dat)

      implicit none

      integer baseFilePnt, baseArrayPnt, nComp, dimFile(3)
      integer FF_DIMDEC(dat)
      REAL_T dat(FF_DIMV(dat),nComp)

      integer n, proc, lo(3), hi(3), loRd(3), hiRd(3), filLo(3)

#include <INFL_FORCE_F.H>

      call bl_pd_myproc(proc)

      call FF_SET_LOHI(FF_DIMS(dat), lo, hi)

      if (SDIM .eq. 2) then
         do n = 1, 3
            filLo(n) = 1
            loRd(n) = lo(n)
            hiRd(n) = hi(n)
         enddo
         loRd(strmwse_dir) = baseArrayPnt
         filLo(strmwse_dir) = baseFilePnt
         do n = 1, nComp
            call RD_FLCTREC(loRd, hiRd, filLo, FF_DIMS(dat),
     &           dat(lo(1),lo(2),lo(3),n),n)
         enddo
         
      else
c     
c     *** Wrap baseFilePnt into the box if it is outside ***
c     
c     Note: This assumes the first and last point in the file 
c     are the same.  This is the case for data generated with mkInitFlct.
c     
         if (baseFilePnt .lt. 1) then
            baseFilePnt = dimFile(strmwse_dir) - (1 - baseFilePnt)
         else if (baseFilePnt .gt. dimFile(strmwse_dir) - 1) then
            baseFilePnt = baseFilePnt - (dimFile(strmwse_dir) - 1)
         endif
c     
c     *** If the data can be filled in one pass, do so ***
c     
         if (hi(strmwse_dir) - baseArrayPnt + 1 
     &        .LE. dimFile(strmwse_dir) - baseFilePnt + 1) then
            do n = 1, 3
               filLo(n) = 1
               loRd(n) = lo(n)
               hiRd(n) = hi(n)
            enddo
            loRd(strmwse_dir) = baseArrayPnt
            filLo(strmwse_dir) = baseFilePnt

            do n = 1, nComp
               call RD_FLCTREC(loRd, hiRd, filLo, FF_DIMS(dat),
     &              dat(lo(1),lo(2),lo(3),n),n)
            enddo
c     
c     Note: In this case, we are guaranteed that the arrays can be filled 
c     in two passes since the array is guaranteed not to have 
c     dimensions larger than the data in the file.
c     
         else
c     
c     Fill as much data as can be read without reading beyond the end of
c     the file
c     
            do n = 1, 3
               filLo(n) = 1
               loRd(n) = lo(n)
               hiRd(n) = hi(n)
            enddo
            filLo(strmwse_dir) = baseFilePnt
            loRd(strmwse_dir) = baseArrayPnt
            hiRd(strmwse_dir) = loRd(strmwse_dir) + dimFile(strmwse_dir) 
     &           - baseFilePnt

            do n = 1, nComp
               call RD_FLCTREC(loRd, hiRd, filLo, FF_DIMS(dat),
     &              dat(lo(1),lo(2),lo(3),n),n)
            enddo
c     
c     Now fill the rest of the array starting from the beginning of the file
c     
c     Note: The first point in the file in the streamwise direction is 
c     skipped because it is identical to the last point in the file.
c     
            filLo(strmwse_dir) = 2
            loRd(strmwse_dir) = hiRd(strmwse_dir) + 1
            hiRd(strmwse_dir) = hi(strmwse_dir)

            do n = 1, nComp
               call RD_FLCTREC(loRd, hiRd, filLo, FF_DIMS(dat),
     &              dat(lo(1),lo(2),lo(3),n),n)
            enddo

         endif
      endif

      END

      subroutine FILL_FRCARRYS_SWIRL(baseFilePnt, baseArrayPnt,
     &     nComp, FF_DIMS(dat), dat)

      implicit none

      integer baseFilePnt, baseArrayPnt, nComp
      integer FF_DIMDEC(dat)
      REAL_T dat(FF_DIMV(dat),nComp)

      integer n
      integer lo(3), hi(3), loRd(3), hiRd(3), filLo(3)

#include <INFL_FORCE_F.H>

      call FF_SET_LOHI(FF_DIMS(dat), lo, hi)

      do n = 1, 3
         filLo(n) = 1
         loRd(n)  = lo(n)
         hiRd(n)  = hi(n)
      enddo
      loRd(strmwse_dir)  = baseArrayPnt
      filLo(strmwse_dir) = baseFilePnt

      do n = 1, nComp
         call RD_FLCTREC(loRd, hiRd, filLo,
     &        FF_DIMS(dat), dat(lo(1),lo(2),lo(3),n),n)
      enddo

      end
