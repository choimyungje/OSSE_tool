module GEMSTOOL_RTCALC_PLUS_m

!  Module files for VLIDORT

   USE VLIDORT_PARS
   USE VLIDORT_IO_DEFS
   USE VLIDORT_MASTERS

   USE VLIDORT_LIN_IO_DEFS
   USE VLIDORT_LPS_MASTERS
   USE VLIDORT_LCS_MASTERS

public  :: GEMSTOOL_RTCALC_PLUS_Actual, GEMSTOOL_RTCALC_PLUS_Final
private :: makechar3

contains

!   (  do_firstorder_option, FO_do_regular_ps,                        & ! FO control

subroutine GEMSTOOL_RTCALC_PLUS_Actual &
   (  do_Profile_Jacobians, do_AerOpdepProfile_Jacobians,                 & ! Main choice of Jacobian
      VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,                          & ! VLIDORT regular inputs
      VLIDORT_LinFixIn, VLIDORT_LinModIn, VLIDORT_LinSup,                 & !VLIDORT linearized inputs
      nlayers, NGREEK_MOMENTS_INPUT, height_grid, lambertian_albedo,      & ! Proxies for VLIDORT
      DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,         & ! Proxies for VLIDORT
      L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,   & ! Linearized Proxies for VLIDORT
      VLIDORT_Out, VLIDORT_LinOut ) 

   implicit none

!  precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  First order flags

!   logical, INTENT(IN) :: do_firstorder_option  ! turn on, using new FO code
!   logical, INTENT(IN) :: FO_do_regular_ps      ! turn on, using FO code regular PS mode

!  VLIDORT Type structures
!  -----------------------

!  Profile Jacobian flags

   logical, intent(in) ::  do_Profile_Jacobians
   logical, intent(in) ::  do_AerOpdepProfile_Jacobians

!  VLIDORT settings structure, Intent(out) here

   TYPE(VLIDORT_Fixed_Inputs), INTENT(INOUT)       :: VLIDORT_FixIn
   TYPE(VLIDORT_Modified_Inputs), INTENT(INOUT)    :: VLIDORT_ModIn

!  VLIDORT supplements i/o structure

   TYPE(VLIDORT_Sup_InOut), INTENT(INOUT)          :: VLIDORT_Sup

!  VLIDORT linearized input structures

   TYPE(VLIDORT_Fixed_LinInputs), INTENT(INOUT)    :: VLIDORT_LinFixIn
   TYPE(VLIDORT_Modified_LinInputs), INTENT(INOUT) :: VLIDORT_LinModIn

!  VLIDORT linearized supplements i/o structure

   TYPE(VLIDORT_LinSup_InOut), INTENT(INOUT)       :: VLIDORT_LinSup

!  Proxy inputs for VLIDORT
!  ------------------------

   INTEGER, intent(in)      :: NLAYERS
   INTEGER, intent(in)      :: NGREEK_MOMENTS_INPUT

   REAL(fpk), intent(in) :: HEIGHT_GRID           ( 0:MAXLAYERS )
   REAL(fpk), intent(in) :: DELTAU_VERT_INPUT     ( MAXLAYERS )
   REAL(fpk), intent(in) :: OMEGA_TOTAL_INPUT     ( MAXLAYERS )
   REAL(fpk), intent(in) :: GREEKMAT_TOTAL_INPUT  ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
   REAL(fpk), intent(in) :: LAMBERTIAN_ALBEDO

!  Linearized

   REAL(fpk), intent(in) :: L_DELTAU_VERT_INPUT    ( MAX_ATMOSWFS, MAXLAYERS )
   REAL(fpk), intent(in) :: L_OMEGA_TOTAL_INPUT    ( MAX_ATMOSWFS, MAXLAYERS )
   REAL(fpk), intent(in) :: L_GREEKMAT_TOTAL_INPUT ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  Outputs
!  -------

!  VLIDORT output structure

   TYPE(VLIDORT_Outputs), intent(inout)              :: VLIDORT_Out

!  VLIDORT linearized output structure

   TYPE(VLIDORT_LinOutputs), intent(inout)           :: VLIDORT_LinOut

!  Local

   logical, parameter :: skip_vlidort = .false.

!  START Of CODE
!  =============

!  Set Final VLIDORT inputs from the PROXIES
!  -----------------------------------------

!  Regular inputs

   VLIDORT_FixIn%Cont%TS_NLAYERS                  = NLAYERS
   VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID           = HEIGHT_GRID
   VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT  = HEIGHT_GRID(NLAYERS)

   VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT    = NGREEK_MOMENTS_INPUT
   VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT     = DELTAU_VERT_INPUT
   VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT    = OMEGA_TOTAL_INPUT
   VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT  = GREEKMAT_TOTAL_INPUT
   VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO     = LAMBERTIAN_ALBEDO

!  Linearized inputs

   VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT    = L_DELTAU_VERT_INPUT
   VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT    = L_OMEGA_TOTAL_INPUT
   VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT = L_GREEKMAT_TOTAL_INPUT

! SANITY CHECK, 27 July 2013
!   write(*,*)VLIDORT_FixIn%Cont%TS_NSTREAMS
!   write(*,*)VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT
!   write(*,*)VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR
!   write(*,*)VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY
!   write(*,*)VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING
!   write(*,*) VLIDORT_ModIn%MBool%TS_DO_SSCORR_NADIR
!   write(*,*) VLIDORT_ModIn%MBool%TS_DO_SSCORR_OUTGOING
!   write(*,*) VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY
!   write(*,*)VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS
!   write(*,*)VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT
!   write(*,*)VLIDORT_FixIn%Cont%TS_NLAYERS
!   write(*,*)VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT
!   write(*,*)VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID
!   write(*,*)VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT
!   write(*,*)VLIDORT_FixIn%Bool%TS_DO_UPWELLING
!   write(*,*)VLIDORT_FixIn%Bool%TS_DO_DNWELLING
!   write(*,*)VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS
!   write(*,*)VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1)
! pause

!  Perform the VLIDORT calculation

   if (skip_vlidort) then
      VLIDORT_Out%Main%TS_Stokes = 1.0d0
      VLIDORT_Out%Status%TS_STATUS_INPUTCHECK  = 0
      VLIDORT_Out%Status%TS_NCHECKMESSAGES     = 0
      VLIDORT_Out%Status%TS_CHECKMESSAGES      = ' '
      VLIDORT_Out%Status%TS_ACTIONS            = ' '
      VLIDORT_Out%Status%TS_STATUS_CALCULATION = 0
      VLIDORT_Out%Status%TS_MESSAGE            = ' '
      VLIDORT_Out%Status%TS_TRACE_1            = ' '
      VLIDORT_Out%Status%TS_TRACE_2            = ' '
      VLIDORT_Out%Status%TS_TRACE_3            = ' '
   else
      if ( do_Profile_Jacobians .or.do_AerOpdepProfile_Jacobians ) then        !  profile Jacobians, also Surface
        CALL vlidort_LPS_master( &
           VLIDORT_FixIn, &
           VLIDORT_ModIn, &
           VLIDORT_Sup,   &
           VLIDORT_Out,   &
           VLIDORT_LinFixIn, &
           VLIDORT_LinModIn, &
           VLIDORT_LinSup,   &
           VLIDORT_LinOut )
      else                                   !  Column (Tshift/Surfp/Aerosol) Jacobians, also surface 
        CALL vlidort_LCS_master( &
           VLIDORT_FixIn, &
           VLIDORT_ModIn, &
           VLIDORT_Sup,   &
           VLIDORT_Out,   &
           VLIDORT_LinFixIn, &
           VLIDORT_LinModIn, &
           VLIDORT_LinSup,   &
           VLIDORT_LinOut )
       endif
   end if

!  Finish

   return
end subroutine GEMSTOOL_RTCALC_PLUS_Actual

subroutine GEMSTOOL_RTCALC_PLUS_Final &
   (  MAXWAV, MaxAlbCoeffs, maxgases, maxmessages, maxaerwfs, W, do_SphericalAlbedo,       & ! GEMSTOOL Dimensions and W
      do_Profile_Jacobians, do_AerOpdepProfile_Jacobians, do_AerBulk_Jacobians,            & ! GEMSTOOL Jacobian control
      do_Surface_Jacobians, do_Tshift_Jacobian, do_Surfp_Jacobian, do_normalized_wfoutput, & ! GEMSTOOL Jacobian control
      do_gases, do_gas_wfs, ngases, nlayers, nstokes, n_geometries, n_aerosol_wfs, dir,    & ! Control numbers/gases
      Albedo, Tshift, SurfPress, Level_Vmrs, AERTAU_UNSCALED, aerosol_pars,                & ! necessary for normalizing
      VLIDORT_Out, VLIDORT_LinOut,                                               & ! VLIDORT Results (includes linearized)
      LAYER_GASOPDEPS, dGASOPDEPS_dV,                                            & ! Layer-to-Level transformations
      STOKES, ACTINIC, REGFLUX, DOLP, DOCP,                                      & ! Main program, GEMSTOOL Stokes Output
      RADIANCE_GASWFS, RADIANCE_AERPROFWFS, RADIANCE_AERBULKWFS,                 & ! Main program, GEMSTOOL Jacobians Output
      RADIANCE_ALBEDOWFS, RADIANCE_TSHIFTWFS, RADIANCE_SURFPWFS,                 & ! Main program, GEMSTOOL Jacobians Output
      AMFSCATWTS, AMFTOTAL, TOTALWF_SAVED,                                       & ! Main program, GEMSTOOL AMF Output
      Errorstatus, nmessages, messages )                                           ! Main program, exception handling

   implicit none

!  precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  GEMSTOOL wavenumber, ClosureCoeffs, gases and MAX-Messages dimensions

   integer, intent(in) :: MAXWAV, MAXALBCOEFFS, maxgases, maxaerwfs, maxmessages

!  wavelength point

   integer, intent(in) :: w

!  Proxy control inputs

   integer, intent(in) :: nlayers, nstokes, n_geometries

!  Flux control flag (GEMSTOOL control)

   logical, intent(in) :: do_SphericalAlbedo

!  Linearization control

   logical, intent(in) :: do_Profile_Jacobians
   logical, intent(in) :: do_AerOpdepProfile_Jacobians
   logical, intent(in) :: do_AerBulk_Jacobians

   logical, intent(in) :: do_Surface_Jacobians

   logical, intent(in) :: do_Tshift_Jacobian
   logical, intent(in) :: do_Surfp_Jacobian

   logical, intent(in) :: do_normalized_wfoutput

!  Direction (TOA upwelling = 1, BOA downwelling = 2)

   integer, intent(in) :: dir

!  Albedo, Tshift, surface-pressure, Level VMRs, aerosol optical depth profile
!    ( necessary for normalized WFs)

   real(fpk), intent(in)  :: ALBEDO (maxwav)
   real(fpk), intent(in)  :: Tshift 
   real(fpk), intent(in)  :: Surfpress
   REAL(fpk), intent(in)  :: LEVEL_VMRS ( 0:maxlayers, maxgases )
   REAL(fpk), intent(in)  :: AERTAU_UNSCALED ( maxlayers )

!  aerosol Bulk property input

   REAL(fpk), intent(in)  :: aerosol_pars ( maxaerwfs )
   integer  , intent(in)  :: n_aerosol_wfs

!  VLIDORT output structures 

   TYPE(VLIDORT_Outputs)   , intent(in)   :: VLIDORT_Out
   TYPE(VLIDORT_LinOutputs), intent(in)   :: VLIDORT_LinOut

!  Additional inputs to develop LEVEL-VMR weighting functions
!  ----------------------------------------------------------

!  Gas control

   integer    , intent(in) :: ngases
   logical    , dimension ( maxgases), INTENT(IN)  :: do_gases
   logical    , dimension ( maxgases), INTENT(IN)  :: do_gas_wfs

!  Layer gas optical depths

    REAL(fpk), intent(in) :: LAYER_GASOPDEPS       ( maxlayers, maxgases )

!  Chain-Rule Transformation from Layer-GASOPDEP to Level-VMR Jacobians

    REAL(fpk), intent(in) :: dGASOPDEPS_dV       ( maxlayers, maxgases, 2 )

!  Output arrays
!  =============

!  Radiances
      
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,4,MAXWAV)      :: STOKES

!  Degree of linear/circular polarization (DOLP/DOCP)

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV)        :: DOLP
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV)        :: DOCP

!  Actinic and Regular Flux output

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,4,MAXWAV)      :: ACTINIC
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,4,MAXWAV)      :: REGFLUX

!  Jacobian output
!  ===============

!  Gas weighting functions w.r.t. LEVEL VMRS --- Radiance only.

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV,0:MAXLAYERS,MAXGASES)  :: RADIANCE_GASWFS

!  Albedo Jacobians

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV,MAXALBCOEFFS)  :: RADIANCE_ALBEDOWFS

!  T-Shift Jacobians

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV)  :: RADIANCE_TSHIFTWFS

!  Surface pressure Jacobians

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV)  :: RADIANCE_SURFPWFS

!  Aerosol profile Jacobians. LAYER OPTICAL DEPTHS

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV,MAXLAYERS)  :: RADIANCE_AERPROFWFS

!  Aerosol Bulk-property Jacobians. New June 2014

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV,MAXAERWFS)  :: RADIANCE_AERBULKWFS

!  AMF Output
!  ==========

!   AMF Scattering weights w.r.t LAYER GAS ODS, AMF TOTAL for each gas

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV,MAXLAYERS,MAXGASES)  :: AMFSCATWTS
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV,MAXGASES)            :: AMFTOTAL
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV,MAXGASES)            :: TOTALWF_SAVED

!  Exception handling
!  ------------------
!    ( Errorstatus = 0 = SUCCCESS,  Errorstatus = 1 = FAILURE,  Errorstatus = 2 = WARNING )

   integer         , intent(out)     ::  Errorstatus
   integer         , intent(inout)   ::  nmessages
   CHARACTER(LEN=*), intent(inout)   ::  messages(maxmessages)

!  Local

   real(fpk)           :: sq2u2, sv2, KN1_GASABS, KN2_GASABS, KN1, KN2, TOTALWF, TOTALOD
   integer             :: nm, v, o1, n, n1, q, q1, qprof, qbulk, nc, g
   character(LEN=3)    :: c3

!  Initialization

   STOKES (:,:,W) = ZERO
   ACTINIC(:,:,W) = ZERO
   REGFLUX(:,:,W) = ZERO
   DOLP   (:,W)   = ZERO
   DOCP   (:,W)   = ZERO

   RADIANCE_GASWFS   (:,W,:,:) = ZERO
   RADIANCE_ALBEDOWFS(:,W,:)   = ZERO
   RADIANCE_TSHIFTWFS(:,W)     = ZERO
   RADIANCE_SURFPWFS(:,W)      = ZERO
   RADIANCE_AERPROFWFS(:,W,:)  = ZERO
   RADIANCE_AERBULKWFS(:,W,:)  = ZERO      ! New, June 2014

   AMFSCATWTS        (:,W,:,:) = ZERO
   AMFTOTAL          (:,W,:)   = ZERO

!  Exception handling

   NM = Nmessages
   Errorstatus = 0

!  Input check - failure in VLIDORT

   if ( VLIDORT_Out%Status%TS_status_inputcheck.eq.vlidort_serious ) then
      call makechar3(w,c3)
      messages(nm + 1) = 'vlidort input check, failed for wavenumber # '//c3
      messages(nm + 2) = 'vlidort will not execute, here are the Input Check messages and Actions:--'
      nm = nm + 2
      nc = VLIDORT_Out%Status%TS_NCHECKMESSAGES
      DO N = 1, VLIDORT_Out%Status%TS_NCHECKMESSAGES
         call makechar3(n,c3)
         messages(nm + 2*n-1) = ' VLIDORT Message # '//C3//' : '// &
                            adjustl(trim(VLIDORT_Out%Status%TS_CHECKMESSAGES(N)))
         messages(nm + 2*n)   = ' VLIDORT Action  # '//C3//' : '// &
                            adjustl(trim(VLIDORT_Out%Status%TS_ACTIONS(N)))
      ENDDO
      nm = nm + 2*nc
      nmessages = nm ; Errorstatus = 1 ; return
   endif

!  Input check - Warning in VLIDORT, which will go on to Execute

   if ( VLIDORT_Out%Status%TS_status_inputcheck.eq.vlidort_serious ) then
      call makechar3(w,c3)
      messages(nm + 1) = 'vlidort input check, Warning for wavenumber # '//c3
      messages(nm + 2) = 'vlidort will carry on with defaults; here are Input Check messages and Actions:--'
      nm = nm + 2
      nc = VLIDORT_Out%Status%TS_NCHECKMESSAGES
      DO N = 1, VLIDORT_Out%Status%TS_NCHECKMESSAGES
         call makechar3(n,c3)
         messages(nm + 2*n-1) = ' VLIDORT Message # '//C3//' : '// &
                            adjustl(trim(VLIDORT_Out%Status%TS_CHECKMESSAGES(N)))
         messages(nm + 2*n)   = ' VLIDORT Action  # '//C3//' : '// &
                            adjustl(trim(VLIDORT_Out%Status%TS_ACTIONS(N)))
      ENDDO
      nm = nm + 2*nc
      nmessages = nm ; Errorstatus = 2
   endif

!  Execution - failure in VLIDORT

   if ( VLIDORT_Out%Status%TS_status_calculation.eq.vlidort_serious ) then
      call makechar3(w,c3)
      messages(nm + 1) = 'vlidort execution, failed for wavenumber # '//c3
      messages(nm + 2) = 'here are the Execution-failure messages :--'
      nm = nm + 2
      messages(nm + 1) = ' VLIDORT Execution Message : '//adjustl(trim(VLIDORT_Out%Status%TS_MESSAGE))
      messages(nm + 2) = ' VLIDORT Execution Trace 1 : '//adjustl(trim(VLIDORT_Out%Status%TS_TRACE_1))
      messages(nm + 3) = ' VLIDORT Execution Trace 2 : '//adjustl(trim(VLIDORT_Out%Status%TS_TRACE_2))
      messages(nm + 4) = ' VLIDORT Execution Trace 3 : '//adjustl(trim(VLIDORT_Out%Status%TS_TRACE_3))
      nm = nm + 4
      nmessages = nm ; Errorstatus = 1 ; return
   endif

!  CARRY ON with CALCULATION (Errorstatus = 0 or 2)
!  ------------------------------------------------

!  Loop over geometries and Number of stokes parameters

   DO V = 1, N_GEOMETRIES
      DO O1 = 1, NSTOKES
         STOKES(V,O1,W) = VLIDORT_Out%Main%TS_STOKES(1,V,O1,DIR)
      ENDDO
   ENDDO

!  Flux Output: Loop over solar_angles and Number of stokes parameters

   if ( do_SphericalAlbedo ) then
      DO V = 1, N_GEOMETRIES
         DO O1 = 1, NSTOKES
             ACTINIC(V,O1,W) = VLIDORT_Out%Main%TS_MEAN_STOKES(1,V,O1,DIR)
             REGFLUX(V,O1,W) = VLIDORT_Out%Main%TS_FLUX_STOKES(1,V,O1,DIR)
         ENDDO
      ENDDO
   endif

!  Set DOLP to zero if no polarization (Nstokes = 1), and skip DOLP calculation
!         degree of polarization (only for NSTOKES = 3 or 4)

   IF ( NSTOKES.gt.1 ) THEN
      DO V = 1, N_GEOMETRIES
             !------------------------------------- DOLP and DOCP calculation
         SQ2U2 = SQRT(VLIDORT_Out%Main%TS_STOKES(1,V,2,DIR)*VLIDORT_Out%Main%TS_STOKES(1,V,2,DIR) + &
                      VLIDORT_Out%Main%TS_STOKES(1,V,3,DIR)*VLIDORT_Out%Main%TS_STOKES(1,V,3,DIR))
         DOLP(V,W) = SQ2U2 / VLIDORT_Out%Main%TS_STOKES(1,V,1,DIR)
         if ( nstokes .eq. 4 ) then
            SV2 = SQRT(VLIDORT_Out%Main%TS_STOKES(1,V,4,DIR)*VLIDORT_Out%Main%TS_STOKES(1,V,4,DIR))
            DOCP(V,W) = SV2 / VLIDORT_Out%Main%TS_STOKES(1,V,1,DIR)
         endif
      ENDDO
   ENDIF

!  Weighting function output, VLIDORT Result = NORMALIZED
!    1 = upper boundary of layer,  2 = lower boundary of layer

   QPROF = 0
   if ( do_Profile_Jacobians ) then

!  Gas loop

     Q = 0
     DO G = 1, NGASES
       IF ( do_gases(g) .and. do_gas_wfs(g) ) then
          Q = Q + 1

!  Transform VLIDORT Layer Jacobians to Level Profile Jacobians
!      ( Unnormalized results )

          DO V = 1, N_GEOMETRIES
            do n = 0, nlayers
              n1 = n + 1
              if ( n.eq.0 ) then
                 if ( LAYER_GASOPDEPS(N1,G) .gt. zero ) then
                    KN1 = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N1,1,V,1,DIR)
                    KN1_GASABS = KN1 / LAYER_GASOPDEPS(N1,G)
                    RADIANCE_GASWFS(V,W,N,G) = KN1_GASABS * dGASOPDEPS_dV(N1,G,1)
                 endif
              else if ( n.eq.nlayers ) then
                 if ( LAYER_GASOPDEPS(N,G) .gt. zero ) then
                    KN2 = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,1,V,1,DIR)
                    KN2_GASABS = KN2 / LAYER_GASOPDEPS(N,G)
                    RADIANCE_GASWFS(V,W,N,G) = KN2_GASABS * dGASOPDEPS_dV(N,G,2)
                 endif
              else
                 if ( (LAYER_GASOPDEPS(N,G).gt.zero).and.(LAYER_GASOPDEPS(N1,G).gt.zero) ) then
                    KN2 = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,1,V,1,DIR)   
                    KN1 = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N1,1,V,1,DIR)
                    KN2_GASABS = KN2 / LAYER_GASOPDEPS(N,G)
                    KN1_GASABS = KN1 / LAYER_GASOPDEPS(N1,G)
                    RADIANCE_GASWFS(V,W,N,G) = KN2_GASABS * dGASOPDEPS_dV(N,G,2) + KN1_GASABS * dGASOPDEPS_dV(N1,G,1)
                 endif
              endif
            enddo
          enddo

!  Generate Normalized Results

          if ( do_normalized_wfoutput ) then
            DO V = 1, N_GEOMETRIES
              do n = 0, nlayers
                RADIANCE_GASWFS(V,W,N,G) = RADIANCE_GASWFS(V,W,N,G) * Level_vmrs(n,g)
              enddo
            ENDDO
          endif

!  End gas loop

       endif
       QPROF = Q
     enddo
   endif

!  AMF scattering weights. Layer quantities

   if ( do_Profile_Jacobians ) then
     DO V = 1, N_GEOMETRIES
       Q = 0
       DO G = 1, NGASES
         IF ( do_gases(g) .and. do_gas_wfs(g) ) then
           Q = Q + 1
!  ../scattering weight AMFS
           do n = 1, nlayers
              if ( LAYER_GASOPDEPS(N,G).gt.zero ) then
                 KN2 = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,1,V,1,DIR)    
                 AMFSCATWTS(V,W,N,G) = - KN2 / STOKES(V,1,W) / LAYER_GASOPDEPS(N,G)
              endif
           enddo
! ../Total AMFs
           TOTALWF = SUM( VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,1:nlayers,1,V,1,DIR) )
           TOTALWF_SAVED(V,W,G) = TOTALWF
           TOTALOD = SUM( LAYER_GASOPDEPS(1:nlayers,G) )
           AMFTOTAL(V,W,G) = - TOTALWF / TOTALOD / STOKES(V,1,W)
         endif
       enddo
     enddo
   endif

!  Aerosol profile Jacobians. VLIDORT Result = NORMALIZED

   if ( do_AerOpdepProfile_Jacobians ) then
      Q = QPROF + 1
      if ( do_normalized_wfoutput ) then
         DO V = 1, N_GEOMETRIES 
           do n = 1, nlayers
              KN1 = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,1,V,1,DIR)
              IF ( AERTAU_UNSCALED(N) .ne. 0.0d0 ) then
                 RADIANCE_AERPROFWFS(V,W,N) = KN1 
              endif
           enddo
        enddo
      else
         DO V = 1, N_GEOMETRIES 
           do n = 1, nlayers
              KN1 = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,1,V,1,DIR)
              IF ( AERTAU_UNSCALED(N) .ne. 0.0d0 ) then
                 RADIANCE_AERPROFWFS(V,W,N) = KN1 / AERTAU_UNSCALED(N)
!                write(*,*), AERTAU_UNSCALED(N), KN1, RADIANCE_AERPROFWFS(V,W,N)
              endif
           enddo
        enddo
      endif
   endif

!  Surface jacobians. VLIDORT result = UNNORMALIZED value.

   if ( do_Surface_Jacobians ) then
      if ( do_normalized_wfoutput ) then
         DO V = 1, N_GEOMETRIES
            RADIANCE_ALBEDOWFS(V,W,1) = albedo(w) * VLIDORT_LinOut%Surf%TS_SURFACEWF(1,1,V,1,DIR)
         ENDDO
      else
         DO V = 1, N_GEOMETRIES
            RADIANCE_ALBEDOWFS(V,W,1) = VLIDORT_LinOut%Surf%TS_SURFACEWF(1,1,V,1,DIR)
         ENDDO
      endif
   endif

!  Column Jacobians
!  ----------------

! (T-shift and/or Surface pressure). VLIDORT results = NORMALIZED VALUE.

   Q = 0
   if ( do_Tshift_Jacobian ) then
      Q = Q + 1
      if ( do_normalized_wfoutput ) then
         DO V = 1, N_GEOMETRIES
            RADIANCE_TSHIFTWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,1,DIR)
         ENDDO
      else
         DO V = 1, N_GEOMETRIES
            RADIANCE_TSHIFTWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,1,DIR) / TSHIFT
         ENDDO
      endif
   endif
   if ( do_Surfp_Jacobian ) then
      Q = Q + 1
      if ( do_normalized_wfoutput ) then
         DO V = 1, N_GEOMETRIES
            RADIANCE_SURFPWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,1,DIR)
         ENDDO
      else
         DO V = 1, N_GEOMETRIES
            RADIANCE_SURFPWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,1,DIR) / SURFPRESS
         ENDDO
      endif
   endif

!  Aerosol Bulk-property Jacobians. New section, June 2014

   QBULK = Q
   if ( do_AerBulk_Jacobians ) then
      DO Q = 1, n_aerosol_wfs
         Q1 = Q + QBULK
         if ( do_normalized_wfoutput ) then
            DO V = 1, N_GEOMETRIES
               RADIANCE_AERBULKWFS(V,W,Q) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q1,1,V,1,DIR)
            ENDDO
         else
            DO V = 1, N_GEOMETRIES
               RADIANCE_AERBULKWFS(V,W,Q) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q1,1,V,1,DIR) / AEROSOL_PARS(Q)
            ENDDO
         endif
      ENDDO
   endif

!  End subroutine

   return
end subroutine GEMSTOOL_RTCALC_PLUS_Final

subroutine makechar3(w,c3)
  implicit none
  integer, intent(in)          :: w
  character(LEN=3),intent(out) :: c3
  c3 = '000'
  if ( w.lt.10 ) then
     write(c3(3:3),'(I1)')W
  else if ( w.gt.99 ) then
     write(c3(1:3),'(I3)')W
  else
     write(c3(2:3),'(I2)')W
  endif
end subroutine makechar3

!  End module

end module GEMSTOOL_RTCALC_PLUS_m
