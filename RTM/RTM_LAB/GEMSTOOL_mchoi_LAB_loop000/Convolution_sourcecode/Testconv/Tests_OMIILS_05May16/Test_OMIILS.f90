program Test_OMIILS

   use OMIILS_Convolution_Tool_m

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Dimensions

   integer, parameter :: E_ndat_c = 200, E_ndat_ILS = 505, E_ndat = 2871

!  Numbers

   integer    :: ndat_c, ndat, nspec
   integer    :: ndat_ILS

!  Fine grid data and wavelengths

   real(fpk)  :: intensity  (E_ndat,3)
   real(fpk)  :: lambdas    (E_ndat)

!  Corase grid Wavelengths, and results

   real(fpk)  :: lambdas_c   (E_ndat_c)
   real(fpk)  :: intensity_c (E_ndat_c,3)

!  OMI control

   integer    :: OMI_PixelIndex, OMI_Channel, OMI_Window

!  ILS wavelengths and Values
!    -- Values will be regridded to the fine wavelengths

   real(fpk) :: ILS_Values  (E_ndat_ILS,E_ndat_c)
   real(fpk) :: ILS_Wavs    (E_ndat_ILS,E_ndat_c)

!  ILS Offsets and Npoints

   integer   :: ILS_offsets (E_ndat_c)
   integer   :: ILS_npoints (E_ndat_c)

!  Exception handling

   logical       :: fail_ILS
   character*100 :: message_ILS

!  Local

   integer   :: w, j, jdum
   real(fpk) :: dum

!  OMI settings

   OMI_PixelIndex = 21
   OMI_Channel    = 2
   OMI_Window     = 3

!  Get coarse

   open(1,file='OMI_MW2_Coarse.dat',status='old')
   read(1,*)ndat_c
   do w = 1, ndat_c
      read(1,*)lambdas_c(w), dum
   enddo
   close(1)

!  Get Fine data

   ndat = 2871 ; nspec = 3
   open(20,file='ExactScalar_D04_Ray_RS_Obsg_S30V00A010_MW2_CF022_Tr106.Out',status='old')
   do j = 1, ndat
       read(20,*)jdum,lambdas(j), dum, intensity(j,1:3), dum,dum,dum,dum,dum,dum,dum,dum,dum
   enddo
   close(20)

!  Extract

   Call Extract_OMIILS &
         ( E_ndat_c, E_ndat_ILS, ndat_c, ndat_ILS, lambdas_c, &
           OMI_PixelIndex, OMI_Channel, OMI_Window,           &
           ILS_wavs, ILS_Values, fail_ILS, message_ILS )

   if ( Fail_ILS ) then
      write(*,*)Trim(message_ILS) ; stop
   endif

   do w = 1, ndat_c
      write(21,44)w,lambdas_c(w),ILS_wavs(1:505,w),ILS_values(1:505,w)
   enddo
   do j = 1, ndat_ILS
      write(31,54)j,ILS_wavs(j,1:ndat_c),ILS_values(j,1:ndat_c)
   enddo

   write(*,*)'Done extraction'

!  Spline

   Call Spline_OMIILS &
         ( E_ndat, E_ndat_c, E_ndat_ILS, ndat, lambdas, ndat_c,      &
           ndat_ILS, ILS_wavs, ILS_Values, ILS_offsets, ILS_npoints, &
           fail_ILS, message_ILS )

   if ( Fail_ILS ) then
      write(*,*)Trim(message_ILS) ; stop
   endif

   do w = 1, ndat_c
      write(11,64)w,lambdas_c(w),ILS_offsets(w),ILS_npoints(w),ILS_values(1:ILS_npoints(w),w)
!      write(*,64)w,lambdas_c(w),ILS_offsets(w),ILS_npoints(w),sum(ILS_values(1:ILS_npoints(w),w))*0.01
   enddo

   write(*,*)'Done Splining'

!  convolve

   Call Convolve_OMIILS ( E_ndat, E_ndat_c, E_ndat_ILS, ndat_c, nspec, &
                          ILS_Values, ILS_offsets, ILS_npoints,        &
                          lambdas, Intensity(:,1:nspec), Intensity_c(:,1:nspec) )

   open(24,file= 'ExactScalar_D04_Ray_RS_Obsg_S30V00A010_MW2_CF022_Tr106.COA', status = 'unknown')
   do w = 1, ndat_c
      write(24,'(f10.5,1p3e17.7)')lambdas_c(w), Intensity_c(w,1:nspec)
   enddo
   close(24)

   write(*,*)'Done Convolution'

!  format

44 format(i4,f10.5,3x,505f10.5,3x,505e11.3)
54 format(i4,3x,152f10.5,3x,152e11.3)
64 format(i4,f10.5,2i5,3x,505e11.3)

!  done

   stop
end program Test_OMIILS
