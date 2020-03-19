Module Scenario_debug_m

!  Debug scenario output
!     1. plot_psddbg : Particle size distribution plot
!                      Only for the fixed aerosol PSD properties
!  Type definitions

   use vlidort_pars, only : fpk, zero

!  Mie distributions (also good for T-matrix)

  use RTSMie_distributions_m

!  Everything public

public

contains

subroutine plot_psddbg & 
    ( filename, do_bimodal, PSDIndex, PSDpars, R1, R2, &
      message_psddbg, fail_psddbg )

!  subroutine for a Plot of the PSDs
!    Divides radius regions into 2 parts for ease of calculation

   implicit none

!  Inputs
!  ------

   character*(*), intent(in) :: filename
   logical      , intent(in) :: do_bimodal
   integer      , intent(in) :: PSDIndex(2)
   real(fpk)    , intent(in) :: PSDpars(3,2)
   real(fpk)    , intent(in) :: R1(2)
   real(fpk)    , intent(in) :: R2(2)

!  Error handling

   character*(*), intent(out) :: message_psddbg
   logical      , intent(out) :: fail_psddbg

!  Local variables

   real(fpk)  :: prad(1),division,diffrad,minrad,maxrad,nwithr1(1),nwithr2(1)
   integer    :: i, nplot_psd_radii(2)

!  two plot regions

   division = 2.0              ! Radius (micron) for division of regions 
   nplot_psd_radii(1) = 1001   ! Number of points in primary region
   nplot_psd_radii(2) = 101    ! Number of points in secondary (tail) region

!  Open file

   open(444,File=trim(filename),status='unknown')

!  One mode
!  --------

   if ( .not. do_bimodal ) then

!  First region

      diffrad = (division-r1(1))/real(nplot_psd_radii(1)-1)
      do i = 1, nplot_psd_radii(1)
         prad(1) = r1(1) + real(i-1)*diffrad
         call sizedis (1,PSDIndex(1),PSDpars(:,1),prad(1),1, &
                       nwithr1(1), message_psddbg, fail_psddbg )
         if (fail_psddbg) go to 56
         write(444,*)i, prad(1), nwithr1(1), zero
      enddo

!  Second region

      if ( r2(1).gt.division )  then
         diffrad = (r2(1)-diffrad)/real(nplot_psd_radii(2)-1)
         do i = 2, nplot_psd_radii(2)
            prad(1) = diffrad + real(i-1)*diffrad
            call sizedis (1,PSDIndex(1),PSDpars(:,1),prad(1),1, &
                          nwithr1(1), message_psddbg, fail_psddbg )
            if (fail_psddbg) go to 56
            write(444,*)i, prad(1), nwithr1(1), zero
         enddo
      endif

   endif

!  Bimodal
!  -------

   if ( do_bimodal ) then

!  First region

      minrad = min(r1(1),r1(2)) ; maxrad = division
      diffrad = ( maxrad - minrad ) / real(nplot_psd_radii(1)-1)
      do i = 1, nplot_psd_radii(1)
         prad(1) = minrad + real(i-1)*diffrad
         nwithr1 = zero ; nwithr2 = zero
         if ( prad(1).ge.r1(1) .and.prad(1).le.r2(1) ) then 
            call sizedis &
              ( 1, PSDIndex(1), PSDpars(:,1) , prad(1), 1, &
                nwithr1(1), message_psddbg, fail_psddbg )
            if (fail_psddbg) go to 56
         endif
         if ( prad(1).ge.r1(2) .and.prad(1).le.r2(2) ) then 
            call sizedis &
              ( 1, PSDIndex(2), PSDpars(:,2) , prad(1), 1, &
                nwithr2(1), message_psddbg, fail_psddbg )
            if (fail_psddbg) go to 56
         endif
         write(444,*)i, prad(1), nwithr1(1), nwithr2(1)
      enddo

!  second region

      if ( (r2(1).gt.division) .or. (r2(2).gt.division) )  then
         minrad = division ; maxrad =  max(r2(1),r2(2))
         diffrad = ( maxrad - minrad ) / real(nplot_psd_radii(2)-1)
         do i = 2, nplot_psd_radii(2)
            prad(1) = division + real(i-1)*diffrad
            nwithr1 = zero ; nwithr2 = zero
            if ( prad(1).ge.r1(1) .and.prad(1).le.r2(1) ) then 
               call sizedis &
               ( 1, PSDIndex(1), PSDpars(:,1) , prad(1), 1, &
                   nwithr1(1), message_psddbg, fail_psddbg )
               if (fail_psddbg) go to 56
            endif
            if ( prad(1).ge.r1(2) .and.prad(1).le.r2(2) ) then 
               call sizedis &
                 ( 1, PSDIndex(2), PSDpars(:,2) , prad(1), 1, &
                   nwithr2(1), message_psddbg, fail_psddbg )
               if (fail_psddbg) go to 56
            endif
            write(444,*)i, prad(1), nwithr1(1), nwithr2(1)
         enddo
      endif

   endif

!  Finish

56 continue
   return
end subroutine plot_psddbg

!  End module

end Module Scenario_debug_m
