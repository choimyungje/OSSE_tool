module GEMSTOOL_NSW_SetJacobians_m

!  Rob Fix, 10/18/16. Introduced SIF
!    GEMSTOOL_RTCALC_PLUS_Actual - SIF inputs added, VLIDORT variables set

!  Rob Fix, 10/25/16. Introduced BRDF, No linearization.
!    GEMSTOOL_RTCALC_PLUS_Actual - BRDF inputs added, VLIDORT variables set

!  Rob Fix 10/25/16
!    GEMSTOOL_RTCALC_PLUS_Final  - New SIF weighting function output.
!    GEMSTOOL_RTCALC_PLUS_Final  - Given UVN/NWS separation  --> SETJACOBIANS

!  Rob Fix, 11/30/16
!    Allowing for linear parameterization of SIF.
!    Renaming of SIF Jacobian variables
       
!  Module files for VLIDORT

   USE VLIDORT_PARS
   USE VLIDORT_IO_DEFS
   USE VLIDORT_LIN_IO_DEFS

public  :: GEMSTOOL_NSW_SetJacobians
private :: makechar3

contains

!

subroutine GEMSTOOL_NSW_SetJacobians &
   (  MAXWAV, MaxAlbCoeffs, maxgases, maxmessages, maxaerwfs, W, dir, do_SphericalAlbedo,    & ! GEMSTOOL Dimensions and W
      do_normalized_wfoutput, do_Profile_Jacobians, do_AerOpdepProfile_Jacobians,            & ! GEMSTOOL Jacobian control
      do_Surface_Jacobians, do_SIF_Jacobians, do_AerBulk_Jacobians,                          & ! GEMSTOOL Jacobian control
      do_Tshift_Jacobian, do_Surfp_Jacobian, do_H2OScaling_Jacobian, do_CH4Scaling_Jacobian, & ! GEMSTOOL Jacobian control
      do_gases, do_gas_wfs, ngases, nlayers, nstokes, n_geometries, n_aerosol_wfs, n_SIFPars,& ! Control numbers/gases
      Albedo, SIFPars, Tshift, SurfPress, H2OScaling, CH4Scaling,                            & ! necessary for normalizing
      Level_Vmrs, AERTAU_UNSCALED, aerosol_pars,                                             & ! necessary for normalizing
      VLIDORT_Out, VLIDORT_LinOut,                                        & ! VLIDORT Results (includes linearized)
      LAYER_GASOPDEPS, dGASOPDEPS_dV,                                     & ! Layer-to-Level transformations
      STOKES, ACTINIC, REGFLUX, DOLP, DOCP,                               & ! Main program, GEMSTOOL Stokes Output
      RADIANCE_GASWFS, RADIANCE_AERPROFWFS, RADIANCE_AERBULKWFS,          & ! Main program, GEMSTOOL Jacobians Output
      RADIANCE_ALBEDOWFS, RADIANCE_SIFPARAMWFS, RADIANCE_TSHIFTWFS,       & ! Main program, GEMSTOOL Jacobians Output
      RADIANCE_SURFPWFS,  RADIANCE_WSCALEWFS,   RADIANCE_MSCALEWFS,       & ! Main program, GEMSTOOL Jacobians Output
      AMFSCATWTS, AMFTOTAL, TOTALWF_SAVED,                                & ! Main program, GEMSTOOL AMF Output
      Errorstatus, nmessages, messages )                                    ! Main program, exception handling

   implicit none

!  precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  GEMSTOOL wavenumber, ClosureCoeffs, gases and MAX-Messages dimensions

   integer, intent(in) :: MAXWAV, MAXALBCOEFFS, maxgases, maxaerwfs, maxmessages

!  wavelength point

   integer, intent(in) :: w

!  Direction (TOA upwelling = 1, BOA downwelling = 2)

   integer, intent(in) :: dir

!  Proxy control inputs

   integer, intent(in) :: nlayers, nstokes, n_geometries

!  Flux control flag (GEMSTOOL control)

   logical, intent(in) :: do_SphericalAlbedo

!  Linearization control

   logical, intent(in) :: do_Profile_Jacobians
   logical, intent(in) :: do_AerOpdepProfile_Jacobians
   logical, intent(in) :: do_AerBulk_Jacobians

   logical, intent(in) :: do_Surface_Jacobians
   logical, intent(in) :: do_SIF_Jacobians ! New 10/18/16, renamed 11/30/16

   logical, intent(in) :: do_Tshift_Jacobian
   logical, intent(in) :: do_Surfp_Jacobian
   logical, intent(in) :: do_H2OScaling_Jacobian
   logical, intent(in) :: do_CH4Scaling_Jacobian

   logical, intent(in) :: do_normalized_wfoutput

!  Number of SIF parameters. 11/30/16

   integer  , intent(in)  :: n_SIFPars

!  Albedo, Tshift, surface-pressure, H2O Scaling, CH4 Scaling, Level VMRs, aerosol optical depth profile
!    ( necessary for normalized WFs ). 
!      SIF parameters added to this list, 10/18/16, 11/30/16 (renamed)

   real(fpk), intent(in)  :: ALBEDO (maxwav)
   real(fpk), intent(in)  :: Tshift 
   real(fpk), intent(in)  :: Surfpress
   real(fpk), intent(in)  :: H2OScaling
   real(fpk), intent(in)  :: CH4Scaling
   REAL(fpk), intent(in)  :: LEVEL_VMRS ( 0:maxlayers, maxgases )
   REAL(fpk), intent(in)  :: AERTAU_UNSCALED ( maxlayers )
   real(fpk), intent(in)  :: SIFPars(2)

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

!  Rob, 10/18/16, Derivatives of Stokes Radiance quantities w.r.t. SIF Parameters.
!       11/30/16, Renamed, added 2 parameters for linear parameterization.

!   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV)  :: RADIANCE_SIFSCALEWFS  ! (old)
   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV,2)  :: RADIANCE_SIFPARAMWFS

!  T-Shift Jacobians

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV)  :: RADIANCE_TSHIFTWFS

!  Surface pressure Jacobians

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV)  :: RADIANCE_SURFPWFS

!  Derivatives of Stokes Radiance quantities w.r.t. H2O Scaling

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV)  :: RADIANCE_WSCALEWFS

!  Derivatives of Stokes Radiance quantities w.r.t. CH4 Scaling        
!  Y.Jung fix 2/1/15

   real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV)  :: RADIANCE_MSCALEWFS

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
   integer             :: nm, v, o1, n, n1, q, q1, qprof, qbulk, nc, g, qs, qoffset
   character(LEN=3)    :: c3

!  Initialization

   STOKES (:,:,W) = ZERO
   ACTINIC(:,:,W) = ZERO
   REGFLUX(:,:,W) = ZERO
   DOLP   (:,W)   = ZERO
   DOCP   (:,W)   = ZERO

   RADIANCE_GASWFS   (:,W,:,:) = ZERO
   RADIANCE_ALBEDOWFS(:,W,:)   = ZERO
   RADIANCE_SIFPARAMWFS(:,W,:) = ZERO      ! New, 11/30/16
   RADIANCE_TSHIFTWFS(:,W)     = ZERO
   RADIANCE_SURFPWFS(:,W)      = ZERO
   RADIANCE_WSCALEWFS(:,W)     = ZERO      ! New, July 2014
   RADIANCE_MSCALEWFS(:,W)     = ZERO      ! New, Feb 2015
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

!  @@ Rob fix 10/18/16, add SIF scaling derivative
!             11/30/16, Add second derivative for linear parameterization
!  SIF-Scaling Jacobians depend on the VLIDORT surface Jacobians results = UNNORMALIZED VALUS
!    --> the first is the albedo Jacobian

   QOFFSET = 0 ; if ( do_Surface_Jacobians ) QOFFSET = 1
   if ( do_SIF_Jacobians ) then
     DO QS = 1, n_SIFPars
       Q = QOFFSET + QS
       if ( do_normalized_wfoutput ) then
         DO V = 1, N_GEOMETRIES
            RADIANCE_SIFPARAMWFS(V,W,QS) = SIFPars(qs) * VLIDORT_LinOut%Surf%TS_SURFACEWF(Q,1,V,1,DIR)
         ENDDO
       else
         DO V = 1, N_GEOMETRIES
            RADIANCE_SIFPARAMWFS(V,W,QS) = VLIDORT_LinOut%Surf%TS_SURFACEWF(Q,1,V,1,DIR)
         ENDDO
       endif
     ENDDO
   endif

!  Column Jacobians
!  ----------------

! (T-shift and/or Surface pressure and/or/ H2O Scaling). VLIDORT results = NORMALIZED VALUE.
!  @@ Rob fix 7/23/14, add H2O scaling derivative

   Q = 0
   if ( do_H2OScaling_Jacobian ) then
      Q = Q + 1
      if ( do_normalized_wfoutput ) then
         DO V = 1, N_GEOMETRIES
            RADIANCE_WSCALEWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,1,DIR)
         ENDDO
      else
         DO V = 1, N_GEOMETRIES
            RADIANCE_WSCALEWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,1,DIR) / H2OSCALING
         ENDDO
      endif
   endif

   if ( do_CH4Scaling_Jacobian ) then
      Q = Q + 1
      if ( do_normalized_wfoutput ) then
         DO V = 1, N_GEOMETRIES
            RADIANCE_MSCALEWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,1,DIR)
         ENDDO
      else
         DO V = 1, N_GEOMETRIES
            RADIANCE_MSCALEWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,1,DIR) / CH4SCALING
         ENDDO
      endif
   endif

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
end subroutine GEMSTOOL_NSW_SetJacobians

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

end module GEMSTOOL_NSW_SetJacobians_m
