module GEMSTOOL_Albedo_m

use GEMSTOOL_Input_Types_m
public

contains

subroutine GEMSTOOL_Albedo &
      ( Maxwav, MaxAlbCoeffs, nwav, wav, do_wavnums, GEMSTOOL_INPUTS, &
        PhysicsDataPath, alb_filename,                                &
        lambertian_albedo, ClosureTerms )
     
   implicit none

!  precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  dimensioning

   integer  , intent(in) :: Maxwav, MaxAlbCoeffs

!  Wavelength input

   integer  , intent(in) :: nwav
   real(fpk), intent(in) :: wav(Maxwav)

!  Flag for using wavenumbers

   logical  , intent(in) :: do_wavnums

! GEOMTOOL inputs, structure. Intent(in) here

   TYPE(GEMSTOOL_Config_Inputs), INTENT(IN) :: GEMSTOOL_INPUTS

!  Albedo output

   real(fpk), intent(out) :: lambertian_albedo(Maxwav)
   real(fpk), intent(out) :: ClosureTerms(Maxwav,MaxAlbCoeffs)

!  Local

   integer   :: w, wc, k
   logical   :: trawl
   real(fpk) :: albedo, poly, fac, closure_ref

!  File paths and file names

   character*(*), intent(in) :: PhysicsDataPath, alb_filename
   character*256             :: file_pth

!  initialize

   lambertian_albedo = 0.0d0
   ClosureTerms      = 0.0d0

!  Simple constant-albedo code for the NSW case.
!    Current default, 15 August 2013

!   if ( do_wavnums ) then
!      if (GEMSTOOL_INPUTS%Closure%n_closure_bands.eq.1 ) then
!         albedo = GEMSTOOL_INPUTS%Closure%closure_coeffs(1,1)
!         lambertian_albedo(1:nwav) = albedo
!         ClosureTerms(1:nwav,1) = 1.0d0
!         write(*,*) lambertian_albedo(1:nwav),'----AAAA'
!      endif
!      return
!   endif

!  Multi-band closure case with wavelengths, UVN application


!  file_pth =  adjustl(trim(PhysicsDataPath))//'Albedo_data/'//adjustl(trim(alb_filename))
!  print*, file_pth
!   if ( do_wavnums ) then
!        open(1, file = trim(file_pth), status = 'old')
!        do w = 1, nwav
!            read(1,*) fac, albedo
!            lambertian_albedo(w) = albedo
!        enddo
!         close(1)
!    endif

!   if ( do_wavnums ) then

   do w = 1, nwav

!  Step 1, find closure band
      wc = 0 ; trawl = .true.
      do while (trawl)
         wc = wc + 1
         if ( wav(w) .ge. GEMSTOOL_INPUTS%Closure%closure_start(wc) .and. &
              wav(w) .le. GEMSTOOL_INPUTS%Closure%closure_finish(wc) ) trawl = .false.
      enddo

!!  Step 2, polynomial factor
     closure_ref = 0.5_fpk * ( GEMSTOOL_INPUTS%Closure%closure_finish(wc) + &
                                GEMSTOOL_INPUTS%Closure%closure_start (wc) )
      poly = ( wav(w) - closure_ref ) ; fac = 1.0_fpk
      ClosureTerms(w,1) = fac

!!  Albedo as polynomial
      albedo = GEMSTOOL_INPUTS%Closure%closure_coeffs(wc,1)
      do k = 2, GEMSTOOL_INPUTS%Closure%closure_ncoeffs(wc)
         fac = fac * poly ; albedo = albedo + fac * GEMSTOOL_INPUTS%Closure%closure_coeffs(wc,k)
        ClosureTerms(w,k) = fac
      enddo
      lambertian_albedo(w) = albedo
   enddo

!   endif

!  Debug
!   open(5, file = 'alb_model_input.out')
!   do w = 1, nwav
!        write(5,*) w, wav(w), lambertian_albedo(w)
!   enddo
!  close(5) 

!  Done

   return
end subroutine GEMSTOOL_Albedo

!  End module

end module GEMSTOOL_Albedo_m
