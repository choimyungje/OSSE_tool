module GEMSTOOL_UVN_SetJacobians_IQU_m

  !  Rob Fix, 10/25/16. Introduced BRDF, No linearization.
  !    GEMSTOOL_RTCALC_PLUS_Actual - BRDF inputs added, VLIDORT variables set
  
  !  Rob Fix 10/25/16
  !    This used to be called GEMSTOOL_RTCALC_PLUS_Final, Now UVN_SETJACOBIANS
  
  !  Module files for VLIDORT
  
     USE VLIDORT_PARS
     USE VLIDORT_IO_DEFS
     USE VLIDORT_LIN_IO_DEFS
  
  public  :: GEMSTOOL_UVN_SetJacobians_IOU
  private :: makechar3
  
  contains
  
  !
  
  subroutine GEMSTOOL_UVN_SetJacobians_IQU &
     (  MAXWAV, MaxAlbCoeffs, maxgases, maxmessages, maxaerwfs, W, do_SphericalAlbedo,         & ! GEMSTOOL Dimensions and W
        do_normalized_wfoutput, do_Profile_Jacobians, do_AerOpdepProfile_Jacobians,            & ! GEMSTOOL Jacobian control
        do_Surface_Jacobians, do_AerBulk_Jacobians, do_Tshift_Jacobian, do_Surfp_Jacobian,     & ! GEMSTOOL Jacobian control
        do_gases, do_gas_wfs, ngases, nlayers, nstokes, n_geometries, n_aerosol_wfs, dir,      & ! Control numbers/gases
        Albedo,  Tshift, SurfPress, Level_Vmrs, AERTAU_UNSCALED, aerosol_pars,                 & ! necessary for normalizing
        VLIDORT_Out, VLIDORT_LinOut,                                        & ! VLIDORT Results (includes linearized)
        LAYER_GASOPDEPS, dGASOPDEPS_dV,                                     & ! Layer-to-Level transformations
        STOKES, ACTINIC, REGFLUX, DOLP, DOCP,                               & ! Main program, GEMSTOOL Stokes Output
        RADIANCE_GASWFS, RADIANCE_AERPROFWFS, RADIANCE_AERBULKWFS,          & ! Main program, GEMSTOOL Jacobians Output
        RADIANCE_ALBEDOWFS, RADIANCE_TSHIFTWFS, RADIANCE_SURFPWFS,          & ! Main program, GEMSTOOL Jacobians Output
        AMFSCATWTS, AMFTOTAL, TOTALWF_SAVED,                                & ! Main program, GEMSTOOL AMF Output
        STOKES_Q_GASWFS,    STOKES_Q_AERPROFWFS,  STOKES_Q_AERBULKWFS,   & ! Main program, GEMSTOOL Q-Jacobians Output
        STOKES_Q_ALBEDOWFS, STOKES_Q_TSHIFTWFS,   STOKES_Q_SURFPWFS,     & ! Main program, GEMSTOOL Q-Jacobians Output
        STOKES_U_GASWFS,    STOKES_U_AERPROFWFS,  STOKES_U_AERBULKWFS,   & ! Main program, GEMSTOOL U-Jacobians Output
        STOKES_U_ALBEDOWFS, STOKES_U_TSHIFTWFS,   STOKES_U_SURFPWFS,     & ! Main program, GEMSTOOL U-Jacobians Output
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
  !    ( necessary for normalized WFs )
  
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
  
  !  Stokes-Q and Stokes-U Jacobian output. Rob, added 12/19/16
  !  ==========================================================

  !  Gas weighting functions w.r.t. LEVEL VMRS

     real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV,0:MAXLAYERS,MAXGASES)  :: STOKES_Q_GASWFS
     real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV,0:MAXLAYERS,MAXGASES)  :: STOKES_U_GASWFS
  
  !  Albedo Jacobians
  
     real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV,MAXALBCOEFFS)  :: STOKES_Q_ALBEDOWFS
     real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV,MAXALBCOEFFS)  :: STOKES_U_ALBEDOWFS
  
  ! !  Rob, 10/18/16, Derivatives of Stokes quantities w.r.t. SIF Parameters.
  ! !       11/30/16, Renamed, added 2 parameters for linear parameterization.
  
  !    real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV,2)  :: STOKES_Q_SIFPARAMWFS
  !    real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV,2)  :: STOKES_U_SIFPARAMWFS
  
  !  T-Shift Jacobians
  
     real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV)  :: STOKES_Q_TSHIFTWFS
     real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV)  :: STOKES_U_TSHIFTWFS
  
  !  Surface pressure Jacobians
  
     real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV)  :: STOKES_Q_SURFPWFS
     real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV)  :: STOKES_U_SURFPWFS
  
  ! !  Derivatives of Stokes Radiance quantities w.r.t. H2O Scaling
  
  !    real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV)  :: STOKES_Q_WSCALEWFS
  !    real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV)  :: STOKES_U_WSCALEWFS
  
  ! !  Derivatives of Stokes Radiance quantities w.r.t. CH4 Scaling        
  ! !  Y.Jung fix 2/1/15
  
  !    real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV)  :: STOKES_Q_MSCALEWFS
  !    real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV)  :: STOKES_U_MSCALEWFS
  
  !  Aerosol profile Jacobians. LAYER OPTICAL DEPTHS
  
     real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV,MAXLAYERS)  :: STOKES_Q_AERPROFWFS
     real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV,MAXLAYERS)  :: STOKES_U_AERPROFWFS
  
  !  Aerosol Bulk-property Jacobians. New June 2014
  
     real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV,MAXAERWFS)  :: STOKES_Q_AERBULKWFS
     real(fpk)   , DIMENSION (MAX_GEOMETRIES,MAXWAV,MAXAERWFS)  :: STOKES_U_AERBULKWFS
     
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
  
  !  Stokes Q/U Jacobians, initialization added 12/19/16 by Rob

     STOKES_Q_GASWFS   (:,W,:,:) = ZERO ; STOKES_U_GASWFS   (:,W,:,:) = ZERO
     STOKES_Q_ALBEDOWFS(:,W,:)   = ZERO ; STOKES_U_ALBEDOWFS(:,W,:)   = ZERO
    !  STOKES_Q_SIFPARAMWFS(:,W,:) = ZERO ; STOKES_U_SIFPARAMWFS(:,W,:) = ZERO
     STOKES_Q_TSHIFTWFS(:,W)     = ZERO ; STOKES_U_TSHIFTWFS(:,W)     = ZERO
     STOKES_Q_SURFPWFS(:,W)      = ZERO ; STOKES_U_SURFPWFS(:,W)      = ZERO
    !  STOKES_Q_WSCALEWFS(:,W)     = ZERO ; STOKES_U_WSCALEWFS(:,W)     = ZERO
    !  STOKES_Q_MSCALEWFS(:,W)     = ZERO ; STOKES_U_MSCALEWFS(:,W)     = ZERO
     STOKES_Q_AERPROFWFS(:,W,:)  = ZERO ; STOKES_U_AERPROFWFS(:,W,:)  = ZERO
     STOKES_Q_AERBULKWFS(:,W,:)  = ZERO ; STOKES_U_AERBULKWFS(:,W,:)  = ZERO


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
  
!  Stokes Q and U Jacobians, Rob 12/19/16
!    Procedure same as above. Only present if NSTOKES > 1

            if ( nstokes.gt.1 ) then
              DO V = 1, N_GEOMETRIES
                do n = 0, nlayers
                  n1 = n + 1
                  if ( n.eq.0 ) then
                    if ( LAYER_GASOPDEPS(N1,G) .gt. zero ) then
                      KN1 = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N1,1,V,2,DIR) ; KN1_GASABS = KN1 / LAYER_GASOPDEPS(N1,G)
                      STOKES_Q_GASWFS(V,W,N,G) = KN1_GASABS * dGASOPDEPS_dV(N1,G,1)
                      KN1 = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N1,1,V,3,DIR) ; KN1_GASABS = KN1 / LAYER_GASOPDEPS(N1,G)
                      STOKES_U_GASWFS(V,W,N,G) = KN1_GASABS * dGASOPDEPS_dV(N1,G,1)
                    endif
                  else if ( n.eq.nlayers ) then
                    if ( LAYER_GASOPDEPS(N,G) .gt. zero ) then
                      KN2 = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,1,V,2,DIR) ; KN2_GASABS = KN2 / LAYER_GASOPDEPS(N,G)
                      STOKES_Q_GASWFS(V,W,N,G) = KN2_GASABS * dGASOPDEPS_dV(N,G,2)
                      KN2 = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,1,V,3,DIR) ; KN2_GASABS = KN2 / LAYER_GASOPDEPS(N,G)
                      STOKES_U_GASWFS(V,W,N,G) = KN2_GASABS * dGASOPDEPS_dV(N,G,2)
                    endif
                  else
                    if ( (LAYER_GASOPDEPS(N,G).gt.zero).and.(LAYER_GASOPDEPS(N1,G).gt.zero) ) then
                      KN2 = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,1,V,2,DIR)  ; KN2_GASABS = KN2 / LAYER_GASOPDEPS(N,G)  
                      KN1 = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N1,1,V,2,DIR) ; KN1_GASABS = KN1 / LAYER_GASOPDEPS(N1,G)
                      STOKES_Q_GASWFS(V,W,N,G) = KN2_GASABS * dGASOPDEPS_dV(N,G,2) + KN1_GASABS * dGASOPDEPS_dV(N1,G,1)
                      KN2 = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,1,V,3,DIR)  ; KN2_GASABS = KN2 / LAYER_GASOPDEPS(N,G)  
                      KN1 = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N1,1,V,3,DIR) ; KN1_GASABS = KN1 / LAYER_GASOPDEPS(N1,G)
                      STOKES_U_GASWFS(V,W,N,G) = KN2_GASABS * dGASOPDEPS_dV(N,G,2) + KN1_GASABS * dGASOPDEPS_dV(N1,G,1)
                    endif
                  endif
                  if ( do_normalized_wfoutput ) then
                    STOKES_Q_GASWFS(V,W,N,G) = STOKES_Q_GASWFS(V,W,N,G) * Level_vmrs(n,g)
                    STOKES_U_GASWFS(V,W,N,G) = STOKES_U_GASWFS(V,W,N,G) * Level_vmrs(n,g)
                  endif
                enddo
              enddo
            endif
  
  !  End gas loop
  
         endif
         QPROF = Q
       enddo
     endif
  
  !  AMF scattering weights. Layer quantities
  !     Make sure optical depths are non-zero,even for total values, Bug Fix 10/28/16
  
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
             if ( TOTALOD.gt.zero ) then
                AMFTOTAL(V,W,G) = - TOTALWF / TOTALOD / STOKES(V,1,W)
             endif
           endif
         enddo
       enddo
     endif
  
  !  debug. 10/28/16
  !   write(997,*)'Next point',W,ngases,do_gases(1:ngases),do_gas_wfs(1:ngases)
  !   write(998,*)'Next point',W
  !   do n = 1, nlayers
  !      write(997,*)n,LAYER_GASOPDEPS(N,1:ngases)
  !      write(998,*)n,AMFSCATWTS(1,W,N,1:ngases)
  !   enddo
  
 !  Aerosol profile Jacobians.
 !  --------------------------

     if ( do_AerOpdepProfile_Jacobians ) then

      !  Increase profile Jacobian count by 1
      
            Q = QPROF + 1
      
      !  Radiance (Stokes-I Jacobians)
      
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
                    endif
                 enddo
              enddo
            endif
      
      !  Stokes Q and U Jacobians, Rob 12/19/16
      !    Procedure same as above. Only present if NSTOKES > 1
      
            if ( nstokes.gt.1 ) then
              if ( do_normalized_wfoutput ) then
                DO V = 1, N_GEOMETRIES 
                  do n = 1, nlayers
                    KN1 = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,1,V,2,DIR)
                    KN2 = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,1,V,3,DIR)
                    IF ( AERTAU_UNSCALED(N) .ne. 0.0d0 ) then
                      STOKES_Q_AERPROFWFS(V,W,N) = KN1
                      STOKES_U_AERPROFWFS(V,W,N) = KN2
                    endif
                  enddo
                enddo
              else
                DO V = 1, N_GEOMETRIES 
                  do n = 1, nlayers
                    KN1 = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,1,V,2,DIR)
                    KN2 = VLIDORT_LinOut%Prof%TS_PROFILEWF(Q,N,1,V,3,DIR)
                    IF ( AERTAU_UNSCALED(N) .ne. 0.0d0 ) then
                      STOKES_Q_AERPROFWFS(V,W,N) = KN1 / AERTAU_UNSCALED(N)
                      STOKES_U_AERPROFWFS(V,W,N) = KN2 / AERTAU_UNSCALED(N)
                    endif
                  enddo
                enddo
              endif
            endif
      
      !  end clause
      endif
  
    !  Surface jacobians
    !  -----------------
    if ( do_Surface_Jacobians ) then

      !  Stokes-I Jacobians
      
            if ( do_normalized_wfoutput ) then
                DO V = 1, N_GEOMETRIES
                  RADIANCE_ALBEDOWFS(V,W,1) = albedo(w) * VLIDORT_LinOut%Surf%TS_SURFACEWF(1,1,V,1,DIR)
                ENDDO
            else
                DO V = 1, N_GEOMETRIES
                  RADIANCE_ALBEDOWFS(V,W,1) = VLIDORT_LinOut%Surf%TS_SURFACEWF(1,1,V,1,DIR)
                ENDDO
            endif
      
      !  Stokes Q and U Jacobians, Rob 12/19/16
      !    Procedure same as above. Only present if NSTOKES > 1
      
            if ( nstokes.gt.1 ) then
              if ( do_normalized_wfoutput ) then
                DO V = 1, N_GEOMETRIES
                  STOKES_Q_ALBEDOWFS(V,W,1) = albedo(w) * VLIDORT_LinOut%Surf%TS_SURFACEWF(1,1,V,2,DIR)
                  STOKES_U_ALBEDOWFS(V,W,1) = albedo(w) * VLIDORT_LinOut%Surf%TS_SURFACEWF(1,1,V,3,DIR)
                ENDDO
              else
                DO V = 1, N_GEOMETRIES
                  STOKES_Q_ALBEDOWFS(V,W,1) = VLIDORT_LinOut%Surf%TS_SURFACEWF(1,1,V,2,DIR)
                  STOKES_U_ALBEDOWFS(V,W,1) = VLIDORT_LinOut%Surf%TS_SURFACEWF(1,1,V,3,DIR)
                ENDDO
              endif
            endif
      
      endif
  
  !  Column Jacobian for T-shift
!  ---------------------------

      if ( do_Tshift_Jacobian ) then

        !  Increase column Jacobian count by 1
        
              Q = Q + 1
        
        !  Radiance (Stokes-I) Jacobian
        
              if ( do_normalized_wfoutput ) then
                 DO V = 1, N_GEOMETRIES
                    RADIANCE_TSHIFTWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,1,DIR)
                 ENDDO
              else
                 DO V = 1, N_GEOMETRIES
                    RADIANCE_TSHIFTWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,1,DIR) / TSHIFT
                 ENDDO
              endif
        
        !  Stokes Q and U Jacobians, Rob 12/19/16
        !    Procedure same as above. Only present if NSTOKES > 1
        
              if ( nstokes.gt.1 ) then
                if ( do_normalized_wfoutput ) then
                  DO V = 1, N_GEOMETRIES
                    STOKES_Q_TSHIFTWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,2,DIR)
                    STOKES_U_TSHIFTWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,3,DIR)
                  ENDDO
                else
                  DO V = 1, N_GEOMETRIES
                    STOKES_Q_TSHIFTWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,2,DIR) / TSHIFT
                    STOKES_U_TSHIFTWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,3,DIR) / TSHIFT
                  ENDDO
                endif
              endif
        
           endif
        
        !  Column Jacobian for Surface pressure
        !  ------------------------------------
        
           if ( do_Surfp_Jacobian ) then
        
        !  Increase column Jacobian count by 1
        
              Q = Q + 1
        
        !  Radiance (Stokes-I) Jacobian
        
              if ( do_normalized_wfoutput ) then
                 DO V = 1, N_GEOMETRIES
                    RADIANCE_SURFPWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,1,DIR)
                 ENDDO
              else
                 DO V = 1, N_GEOMETRIES
                    RADIANCE_SURFPWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,1,DIR) / SURFPRESS
                 ENDDO
              endif
        
        !  Stokes Q and U Jacobians, Rob 12/19/16
        !    Procedure same as above. Only present if NSTOKES > 1
        
              if ( nstokes.gt.1 ) then
                if ( do_normalized_wfoutput ) then
                  DO V = 1, N_GEOMETRIES
                    STOKES_Q_SURFPWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,2,DIR)
                    STOKES_U_SURFPWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,3,DIR)
                  ENDDO
                else
                  DO V = 1, N_GEOMETRIES
                    STOKES_Q_SURFPWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,2,DIR) / SURFPRESS
                    STOKES_U_SURFPWFS(V,W) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q,1,V,3,DIR) / SURFPRESS
                  ENDDO
                endif
              endif
        
           endif
        
        !  Column Jacobian for Aerosol Bulk-properties
        !  -------------------------------------------
        
        !  Aerosol Bulk-property Jacobians. New section, June 2014
        !   Initialize the bulk count overall
        
           QBULK = Q
        
           if ( do_AerBulk_Jacobians ) then
              DO Q = 1, n_aerosol_wfs
        
        !  Increase bulk properties count by 1
        
                 Q1 = Q + QBULK
        
        !  Radiance (Stokes-I) Jacobian
        
                 if ( do_normalized_wfoutput ) then
                    DO V = 1, N_GEOMETRIES
                       RADIANCE_AERBULKWFS(V,W,Q) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q1,1,V,1,DIR)
                    ENDDO
                 else
                    DO V = 1, N_GEOMETRIES
                       RADIANCE_AERBULKWFS(V,W,Q) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q1,1,V,1,DIR) / AEROSOL_PARS(Q)
                    ENDDO
                 endif
        
        !  Stokes Q and U Jacobians, Rob 12/19/16
        !    Procedure same as above. Only present if NSTOKES > 1
        
                 if ( nstokes.gt.1 ) then
                   if ( do_normalized_wfoutput ) then
                     DO V = 1, N_GEOMETRIES
                       STOKES_Q_AERBULKWFS(V,W,Q) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q1,1,V,2,DIR)
                       STOKES_U_AERBULKWFS(V,W,Q) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q1,1,V,3,DIR)
                     ENDDO
                   else
                     DO V = 1, N_GEOMETRIES
                       STOKES_Q_AERBULKWFS(V,W,Q) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q1,1,V,2,DIR) / AEROSOL_PARS(Q)
                       STOKES_U_AERBULKWFS(V,W,Q) = VLIDORT_LinOut%Col%TS_COLUMNWF(Q1,1,V,3,DIR) / AEROSOL_PARS(Q)
                     ENDDO
                   endif
                 endif
        
              ENDDO
           endif
        
        !  End subroutine
  
     return
  end subroutine GEMSTOOL_UVN_SetJacobians_IQU
  
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
  
  end module GEMSTOOL_UVN_SetJacobians_IQU_m
  