module FastV2p7_LRRS2p3_Master_m

!  This is an interface to the LIDORT RRS code.
!   -- 01 April 2016
!  Modified by Myungje Choi 04/19/2019
!  share parameters for mono calculation
!  NPOINTS_MONO_read, W_EXCIT_read, LAMBDA_EXCIT_read,
!   -- 07 Julay 2019
!  N_LOUTPUT_read added
   
   USE LRRS_PARS, only : fpk, zero, one, half, LRRS_SUCCESS, max_bins, max_points
   USE LRRS_INPUTS_DEF
   USE LRRS_OUTPUTS_DEF
   USE LRRS_IO_CHECK
   USE LRRS_MAIN_MASTER

public
! private

contains

subroutine FastV2p7_LRRS2p3_Master &
    ( E_ndat_c, E_nlayers, E_nmoments,                    &
      ndat_c, nlayers, nstreams, nmoments_input, Eradius, &    
      sza, vza, azm, heights, temperatures, Aircolumns,   &
      lambdas_c, fluxes_c, Albedos_c,                     &
      RayXsec_c, Depol_c, deltau_c, omegamoms_c,          &
      elastic_c, raman_c, filling_c, Filling_SS_c,        &
      NPOINTS_MONO_read, W_EXCIT_read, LAMBDA_EXCIT_read, N_LOUTPUT_read, & !by mchoi; reading
      fail, n_messages, messages )
   !  
   implicit none

!  inputs
!  ------

!  External Dimensioning

   integer, intent(in)  :: E_ndat_c, E_nlayers, E_nmoments

!  numbers

   integer, intent(in)  :: ndat_c, nlayers, nstreams, nmoments_input

!  Earth radius

   real(kind=fpk), intent(in) :: Eradius

!  Angles (degrees)

   real(kind=fpk), intent(in) :: sza, vza, azm

!  Heights, temperatures, aircolumns

   real(kind=fpk), intent(in) :: heights       ( 0:E_nlayers )
   real(kind=fpk), intent(in) :: Temperatures  ( E_nlayers )
   real(kind=fpk), intent(in) :: Aircolumns    ( E_nlayers )

!  Wavelengths, fluxes, albedos

   real(kind=fpk), intent(in) :: lambdas_c ( E_ndat_c )
   real(kind=fpk), intent(in) :: fluxes_c  ( E_ndat_c )
   real(kind=fpk), intent(in) :: albedos_c ( E_ndat_c )

!  optical information

   real(kind=fpk), intent(in) :: RayXsec_c   ( E_ndat_c )
   real(kind=fpk), intent(in) :: Depol_c     ( E_ndat_c )
   real(kind=fpk), intent(in) :: Deltau_c    ( E_nlayers, E_ndat_c )
   real(kind=fpk), intent(in) :: Omegamoms_c ( E_nlayers, 0:E_nmoments, E_ndat_c )

!  outputs
!  -------

!  RT results

   real(kind=fpk), intent(out) :: Elastic_c   ( E_ndat_c )
   real(kind=fpk), intent(out) :: Raman_c     ( E_ndat_c )
   real(kind=fpk), intent(out) :: Filling_c   ( E_ndat_c )
   real(kind=fpk), intent(out) :: Filling_SS_c ( E_ndat_c )

!  Exception handling

   logical      , intent(out) :: fail
   integer      , intent(out) :: n_messages

   character*(*), intent(out) :: messages(100)

!  Local
!  -----

!  LIDORT-RRS input & output type structures

   TYPE(LRRS_Fixed_Inputs)       :: LRRS_FixIn
   TYPE(LRRS_Modified_Inputs)    :: LRRS_ModIn
   TYPE(LRRS_Outputs)            :: LRRS_Out

!  other variables

   integer   :: w, wc, woff1, woff2, m
   real(fpk) :: num, nup, sig_0, sig_t, lambda_first, lambda_last, bin1, bin2, binw


! !  for momo calculation; mchoi
   integer, intent(in) :: NPOINTS_MONO_read, W_EXCIT_read, n_loutput_read
   real(kind=fpk), intent(in) :: LAMBDA_EXCIT_read


!  INITIALIZE
!  ==========

!  outputs

   elastic_c    = zero
   Raman_c      = zero
   Filling_c    = zero
   Filling_SS_c = zero

!  Exception handling

   fail       = .false.
   n_messages = 0
   messages   = ' '

!  Fixed  Booleans

   LRRS_FixIn%Bool%DO_RRS_OVERALL    = .true.
   LRRS_FixIn%Bool%DO_DIRECTRRS_ONLY = .false.

   LRRS_FixIn%Bool%DO_UPWELLING   = .true.
   LRRS_FixIn%Bool%DO_DNWELLING   = .false.

   LRRS_FixIn%Bool%DO_DIRECT_BEAM      = .true.
   LRRS_FixIn%Bool%DO_SSCORR_OUTGOING  = .true.
   LRRS_FixIn%Bool%DO_SSFULL           = .false.

   if ( nmoments_input.gt.2 ) then
      LRRS_FixIn%Bool%DO_MOLECSCAT_ONLY  = .false.
      LRRS_FixIn%Bool%DO_DELTAM_SCALING  = .true.
      LRRS_FixIn%Bool%DO_DOUBLE_CONVTEST = .true.
   else
      LRRS_FixIn%Bool%DO_MOLECSCAT_ONLY  = .true. ! Rayleigh-only
      LRRS_FixIn%Bool%DO_DELTAM_SCALING  = .false.
      LRRS_FixIn%Bool%DO_DOUBLE_CONVTEST = .false.
   endif

   LRRS_FixIn%Bool%DO_LBOUNDARIES     = .true.  ! Always true

   LRRS_FixIn%Bool%DO_PLANE_PARALLEL  = .false.

   LRRS_FixIn%Bool%DO_MVOUT_ONLY        = .false.
   LRRS_FixIn%Bool%DO_ADDITIONAL_MVOUT  = .false.

   LRRS_FixIn%Bool%DO_LAMBERTIAN_SURFACE  = .true.
!mick fix 11/11/2016 - initialized DO_WATERLEAVING
   LRRS_FixIn%Bool%DO_WATERLEAVING        = .false.
   LRRS_FixIn%Bool%DO_SHADOW_EFFECT       = .false.
   LRRS_FixIn%Bool%DO_GLITTER_DBMS        = .false.
   LRRS_FixIn%Bool%DO_LAMBERTIAN_FRACTION = .false.

   LRRS_FixIn%Bool%DO_LRRS_WRITEINPUT    = .false.
   LRRS_FixIn%Bool%DO_LRRS_WRITESCENARIO = .false.
   LRRS_FixIn%Bool%DO_LRRS_WRITEFOURIER  = .false.
   LRRS_FixIn%Bool%DO_LRRS_WRITERESULTS  = .false.

!  Modified Inputs (All)

   LRRS_ModIn%MBool%DO_BIN_REALIZATION    = .false.  ! Always true
   LRRS_ModIn%MBool%DO_MONO_REALIZATION   = .true.

   LRRS_ModIn%MBool%DO_ENERGY_BALANCING   = .true.   ! Always true
   LRRS_ModIn%MBool%DO_CABANNES_RAMAN     = .false.

   LRRS_ModIn%MBool%DO_USER_STREAMS   = .true.
   LRRS_ModIn%MBool%DO_MSMODE_LRRS    = .false.
   LRRS_ModIn%MBool%DO_ELASTIC_ONLY   = .false.
   LRRS_ModIn%MBool%DO_NO_AZIMUTH     = .false.
   LRRS_ModIn%MAtmos%GEOMETRY_SPECHEIGHT = zero    !  Zeroed here, set later

!  Fixed integers and numbers

   LRRS_FixIn%Cont%NSTREAMS       = nstreams       ! Number of streams
   LRRS_FixIn%Cont%NLAYERS        = nlayers        ! Number of layers
   LRRS_FixIn%Cont%NMOMENTS_INPUT = nmoments_input ! 2 for Rayleigh-only
   LRRS_FixIn%Cont%NFINELAYERS    = 2
   LRRS_FixIn%Cont%PROGRESS       = 20 ! Screen output: Point frequency

   LRRS_FixIn%Cont%FLUX_FACTOR           = one
   LRRS_FixIn%Cont%SOLAR_ANGLE           = sza       ! Solar senith angle
   LRRS_FixIn%Cont%EARTH_RADIUS          = Eradius   ! From Main program
   LRRS_FixIn%Cont%LRRS_ACCURACY         = 1.0d-05
   LRRS_FixIn%Cont%LRRS_FGSMALL          = 1.0d-03

   LRRS_FixIn%Cont%DEBUG_FILENAMES       = ' '

   LRRS_FixIn%UserVal%N_USER_STREAMS     = 1
   LRRS_FixIn%UserVal%USER_ANGLES(1)     = vza ! View zenith angle
   LRRS_FixIn%UserVal%N_USER_RELAZMS     = 1
   LRRS_FixIn%UserVal%USER_RELAZMS(1)    = azm ! Azimuth angle

   LRRS_FixIn%UserVal%N_LOUTPUT          = N_LOUTPUT_read
!mick fix 10/11/2016 - initialized LPARTIALS_OUTPUT
   LRRS_FixIn%UserVal%LPARTIALS_OUTPUT   = 1
   LRRS_FixIn%UserVal%LBOUNDARIES_OUTPUT = 0
!  for mono calculation - mchoi
   LRRS_FixIn%Spect%NPOINTS_MONO  = NPOINTS_MONO_read !fixed for mono
   LRRS_FixIn%Spect%W_EXCIT       = W_EXCIT_read !fixed for mono and wavelength order
   LRRS_FixIn%Spect%LAMBDA_EXCIT  = LAMBDA_EXCIT_read


!    write(*,* ) LRRS_FixIn%Spect%NPOINTS_MONO
!    write(*,* ) LRRS_FixIn%Spect%W_EXCIT
!    write(*,* ) LRRS_FixIn%Spect%LAMBDA_EXCIT
! pause


!mick fix 10/11/2016 - initialized LAMBERTIAN_FRACTION
   LRRS_FixIn%Surf%LAMBERTIAN_FRACTION = one
   LRRS_FixIn%Surf%NSTREAMS_BRDF  = 0
   LRRS_FixIn%Surf%BRDF_NAMES     = ' '
   LRRS_FixIn%Surf%WHICH_BRDF     = 0
   LRRS_FixIn%Surf%BRDF_FACTOR    = zero
   LRRS_FixIn%Surf%BRDF_NPARS     = 0
   LRRS_FixIn%Surf%BRDF_PARS      = zero

!  Following are Zeroed here, will be set in the rest of the routine.

   LRRS_FixIn%Spect%NPOINTS_INNER = 0
   LRRS_FixIn%Spect%OFFSET_INNER  = 0
   LRRS_FixIn%Spect%NPOINTS_OUTER = 0

   LRRS_FixIn%Spect%LAMBDAS_RANKED = zero
   LRRS_FixIn%Spect%FLUXES_RANKED  = zero
   LRRS_FixIn%Spect%BINLOWER       = zero
   LRRS_FixIn%Spect%BINUPPER       = zero

   LRRS_FixIn%Atmos%HEIGHT_GRID         = zero
   LRRS_FixIn%Atmos%LAYER_TEMPERATURES  = zero
   LRRS_FixIn%Atmos%LAYER_AIRCOLUMNS    = zero

   LRRS_FixIn%Atmos%RAYLEIGH_XSEC   = zero
   LRRS_FixIn%Atmos%RAYLEIGH_DEPOL  = zero

   LRRS_FixIn%Atmos%DELTAU_INPUT_UNSCALED       = zero
   LRRS_FixIn%Atmos%OMEGAMOMS_ELASTIC_UNSCALED  = zero

   LRRS_FixIn%Surf%ALBEDOS_RANKED = zero

!  1. Do the Buffering
!  -------------------

!  find inner wavelength limits, determined by the highest Raman
!  Stokes and Anti-Stokes shifts (194.30 and 196.83 wavenumbers)
!                                 =>  

!Rob/mick fix 10/11/2016 - chg sign on wvn delta value of "num"
   !num = 194.3015d0
   num = - 194.3015d0
   sig_0  = 1.0d+7 / lambdas_c(1)
   sig_t  = sig_0 - num
   lambda_first = 1.0d+7 / sig_t

   nup = + 196.8269d0
   sig_0  = 1.0d+7 / lambdas_c(Ndat_c)
!mick fix 10/11/2016 - subtract "nup"
   !sig_t  = sig_0 + nup
   sig_t  = sig_0 - nup
   lambda_last = 1.0d+7 / sig_t

!  Trawl through OMI wavelengths, find number of inner points and offset.

   woff1 = 0 ; woff2 = 0 ; wc = 0
   do w = 1, ndat_c
      if ( lambdas_c(w) .lt. lambda_first ) then
         woff1 = woff1 + 1
      else
         if ( lambdas_c(w) .gt. lambda_last ) then
            woff2 = woff2 + 1
         else
            wc = wc + 1
         endif
      endif
   enddo

!  exception handling

   if ( woff1+woff2.gt.max_bins ) then
      messages(n_messages + 1) = 'increase the Dimensioning for Raman bins'
      n_messages = n_messages + 1
      fail = .true. ; return
   endif
   if ( ndat_c.gt.max_points ) then
      messages(n_messages + 1) = 'increase the Dimensioning for Maxpoints in LRRS'
      n_messages = n_messages + 1
      fail = .true. ; return
   endif

!  Set RRS buffering for binning calculations:

   LRRS_FixIn%Spect%NPOINTS_INNER = wc
   LRRS_FixIn%Spect%OFFSET_INNER  = woff1
   LRRS_FixIn%Spect%NPOINTS_OUTER = ndat_c
   LRRS_FixIn%Spect%LAMBDAS_RANKED(1:ndat_c) = lambdas_c(1:ndat_c)

!  set the Bin limits

   bin1 = half * ( lambdas_c(2) - lambdas_c(1) )
   LRRS_FixIn%Spect%BINLOWER(1) = lambdas_c(1) - bin1
   bin2 = half * ( lambdas_c(ndat_c) - lambdas_c(ndat_c-1) )
   LRRS_FixIn%Spect%BINUPPER(ndat_c) = lambdas_c(ndat_c) +bin2
   do w = 2, ndat_c
      binw = half * ( lambdas_c(w) + lambdas_c(w-1) )
      LRRS_FixIn%Spect%BINLOWER(w)   = binw
      LRRS_FixIn%Spect%BINUPPER(w-1) = binw
   enddo
   bin2 = half * ( lambdas_c(ndat_c) - lambdas_c(ndat_c-1) )
   LRRS_FixIn%Spect%BINUPPER(ndat_c) = lambdas_c(ndat_c) + bin2

!  2. Fixed Atmosphere/Optical Inputs
!  ----------------------------------

!  Fluxes and Albedos

   LRRS_FixIn%Spect%FLUXES_RANKED (1:ndat_c) = fluxes_c(1:ndat_c)
   LRRS_FixIn%Surf%ALBEDOS_RANKED (1:ndat_c) = albedos_c(1:ndat_c)
   ! write(*,*) LRRS_FixIn%Spect%FLUXES_RANKED (1:ndat_c)
   ! write(*,*) LRRS_FixIn%Surf%ALBEDOS_RANKED (1:ndat_c)
   ! pause

!  Atmospheric quantities

   LRRS_FixIn%Atmos%HEIGHT_GRID(0:nlayers)        = heights     (0:nlayers)
   LRRS_FixIn%Atmos%LAYER_TEMPERATURES(1:nlayers) = Temperatures(1:nlayers)
   LRRS_FixIn%Atmos%LAYER_AIRCOLUMNS(1:nlayers)   = Aircolumns  (1:nlayers)
   LRRS_ModIn%MAtmos%GEOMETRY_SPECHEIGHT = heights (nlayers)


   ! write(*,*) LRRS_FixIn%Atmos%HEIGHT_GRID(0:nlayers)
   ! write(*,*)
   ! write(*,*) LRRS_FixIn%Atmos%LAYER_TEMPERATURES(1:nlayers)
   ! write(*,*)
   ! write(*,*) LRRS_FixIn%Atmos%LAYER_AIRCOLUMNS(1:nlayers)
   ! write(*,*)
   ! write(*,*) LRRS_ModIn%MAtmos%GEOMETRY_SPECHEIGHT

   ! pause


!  optical input

   LRRS_FixIn%Atmos%RAYLEIGH_XSEC(1:ndat_c)  = RayXsec_c(1:ndat_c)
   LRRS_FixIn%Atmos%RAYLEIGH_DEPOL(1:ndat_c) = Depol_c(1:ndat_c)

   LRRS_FixIn%Atmos%DELTAU_INPUT_UNSCALED(1:nlayers,1:ndat_c)      = &
                       Deltau_c(1:nlayers,1:ndat_c)
   LRRS_FixIn%Atmos%OMEGAMOMS_ELASTIC_UNSCALED(1:nlayers,0:nmoments_input,1:ndat_c) = &
                       Omegamoms_c(1:nlayers,0:nmoments_input,1:ndat_c)

   ! write(*,*) LRRS_FixIn%Atmos%RAYLEIGH_XSEC(1:ndat_c)
   ! write(*,*)
   ! write(*,*) LRRS_FixIn%Atmos%RAYLEIGH_DEPOL(1:ndat_c)
   ! write(*,*)
   ! ! write(*,*) LRRS_FixIn%Atmos%LAYER_AIRCOLUMNS(1:nlayers)
   ! ! write(*,*)
   ! ! write(*,*) LRRS_ModIn%MAtmos%GEOMETRY_SPECHEIGHT

   ! pause

!  Debug

!   do w = 1, LRRS_FixIn%Spect%NPOINTS_OUTER
!     write(*,*)w,LRRS_FixIn%Spect%LAMBDAS_RANKED(w),LRRS_FixIn%Atmos%RAYLEIGH_DEPOL(w),LRRS_FixIn%Atmos%RAYLEIGH_XSEC(w)
!   enddo
!   do w = 1, nlayers
!      write(*,*)w, LRRS_FixIn%Atmos%HEIGHT_GRID(w),LRRS_FixIn%Atmos%LAYER_TEMPERATURES(w),LRRS_FixIn%Atmos%LAYER_AIRCOLUMNS(w)
!   enddo
!   do w = 1, nlayers
!      write(*,*)w,LRRS_FixIn%Atmos%OMEGAMOMS_ELASTIC_UNSCALED(w,0:2,22),LRRS_FixIn%Atmos%DELTAU_INPUT_UNSCALED(w,22)
!   enddo
!   pause

!  3. Call RRS and set filling
!  ---------------------------

!  Call to the LRRS Master

   CALL RAMAN_MASTER (LRRS_FixIn, LRRS_ModIn, LRRS_Out)

!  Exception handling. If set, collect errors and return

   IF ( LRRS_Out%Status%STATUS_CALCULATION .NE. LRRS_SUCCESS ) THEN
     N_MESSAGES = LRRS_Out%Status%N_MESSAGES
     do m = 1, n_messages
       MESSAGES(m)   = Trim(LRRS_Out%Status%MESSAGES(m))
     enddo  
     fail = .true. ; return
   END IF

!  compute filling and results

   ! do w = 1, LRRS_FixIn%Spect%NPOINTS_INNER
   !    wc = w + LRRS_FixIn%Spect%OFFSET_INNER
   !    elastic_c(wc) = LRRS_Out%Main%ELASTIC_UP(1,1,w)
   !    Raman_c(wc)   = LRRS_Out%Main%RAMAN_UP(1,1,w)
   !    Filling_c(wc)  = one - (elastic_c(wc)/Raman_c(wc))
   ! enddo

   elastic_c(1) = LRRS_Out%Main%ELASTIC_UP(1,1,1)
   Raman_c(1)   = LRRS_Out%Main%RAMAN_UP(1,1,1)
   Filling_c(1)  = one - (elastic_c(1)/Raman_c(1))
   !one - (elastic_c(W_EXCIT_read)/Raman_c(W_EXCIT_read))
   ! write(*,*)LRRS_Out%Main%ELASTIC_UP(1,1,1), LRRS_Out%Main%RAMAN_UP(1,1,1), &
   !          one - (LRRS_Out%Main%ELASTIC_UP(1,1,1))/(LRRS_Out%Main%RAMAN_UP(1,1,1))
   

!  Finish

   return
end subroutine FastV2p7_LRRS2p3_Master

!  End module

end module FastV2p7_LRRS2p3_Master_m

