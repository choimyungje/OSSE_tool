module GEMSTOOL_NSW_linearized_iops_m

!  routines public

private
public :: GEMSTOOL_NSW_linearized_iops

contains

subroutine GEMSTOOL_NSW_linearized_iops &
     ( maxlayers, maxmoments_input, maxstokes_sq, max_atmoswfs,                & ! VLIDORT dimensions
       maxwav, maxaermoms, maxcldmoms, maxgases,                               & ! GEMSTOOL dimensions
       IOP_debug_filename, IOP_Debug_flag,                                     & ! Debug IOP check 
       do_Profile_Jacobians, do_AerOpdepProfile_Jacobians,                     & ! Profile Jacobian control flags
       do_Tshift_Jacobian, do_Surfpress_Jacobian,                              & ! Column  Jacobian control flags
       w, wavnum, omega_lim, do_gases, do_gas_wfs, ngases, nstokes, nlayers, nwav,     & ! GEMSTOOL Inputs (Control)
       do_aerosols, aerosol_nscatmoms, do_clouds_V2, cloud_nscatmoms,          & ! GEMSTOOL Inputs (Aer/clouds)
       gascolumns, dGASCOLUMNS_dS, dGASCOLUMN_dP, dGASCOLUMNS_dV,              & ! GEMSTOOL Inputs (GasColumns + derivs)
       aircolumns, dAIRCOLUMNS_dS, dAIRCOLUMN_dP, Tshift, SurfPress,           & ! GEMSTOOL Inputs (AirColumns + derivs)
       gas_xsecs, dgasxsecs_dS, dgasxsecs_dP, Rayleigh_xsecs, Rayleigh_depols, & ! GEMSTOOL Inputs (Cross-sections, depol)
       aerosol_layerflags, aerosol_deltau, aerosol_ssalbs, aerosol_scatmoms,   & ! GEMSTOOL Inputs (Aerosol)
       cloud_layerflags, cloud_deltau, cloud_ssalbs, cloud_scatmoms,           & ! GEMSTOOL Inputs (Clouds)
       DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT,                                   & ! VLIDORT Outputs (Regular iops)
       NGREEK_MOMENTS_INPUT, GREEKMAT_TOTAL_INPUT,                             & ! VLIDORT Outputs (Regular iops)
       L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,       & ! VLIDORT Outputs (linearized)
       GASOPDEPS, dGASOPDEPS_dV )                                                ! Layer-to-Level transformation

   implicit none

!  precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  VLIDORT INPUTS
!  ==============

!  Dimensioning

    integer, INTENT(IN)  :: maxlayers, maxmoments_input, maxstokes_sq, max_atmoswfs

!  GEMSTOOL INPUTS
!  ==============

!  GEMSTOOL Dimensioning

   integer, INTENT(IN)  ::  maxwav, maxaermoms, maxcldmoms, maxgases

!  Debug Iop check

   logical      , intent(in) :: IOP_debug_flag
   character*(*), intent(in) :: IOP_debug_filename
  
!  Control

    integer, INTENT(IN)  ::  nstokes, nlayers

!  Wavenumber index and value

    integer  , INTENT(IN)  :: w, nwav
    real(fpk), INTENT(IN)  :: wavnum

!  Limiting single scatter albedo
    !	This should be pre-set to coincide with the internal limit in
    !	the LIDORT models. Suggested value = 0.999999

    real(fpk),    INTENT(IN) :: omega_lim

!  Linearization control

    logical      , intent(in) :: do_Profile_Jacobians
    logical      , intent(in) :: do_AerOpdepProfile_Jacobians

    logical      , intent(in) :: do_Tshift_Jacobian
    logical      , intent(in) :: do_SurfPress_Jacobian

!  Gas control

    logical    , dimension ( maxgases), INTENT(IN)  :: do_gases
    logical    , dimension ( maxgases), INTENT(IN)  :: do_gas_wfs

!  Aerosol control

    logical, INTENT(IN)  :: do_aerosols

!  Cloud control

    logical, INTENT(IN)  :: do_clouds_V2

!  number of gases, aerosol/cloud scattering moments

    integer, INTENT(IN)  :: ngases
    integer, INTENT(IN)  :: aerosol_nscatmoms
    integer, INTENT(IN)  :: cloud_nscatmoms

!  Trace gas stuff
    !  Trace gas Layer column amounts
    !  Gas-column derivative with respect to Level VMRs. (needed for LAYERS-to_LEVELS Transformation)
    !  Gas-column derivative with respect to T-shift
    !  Gas-column derivative with respect to Surface Pressure
    !	LAYER gas columns are in mol/cm^2  or inverse [DU]

    real(fpk),    dimension ( maxlayers, maxgases, 2) , INTENT(IN) :: gascolumns
    real(fpk),    dimension ( maxlayers, maxgases, 2) , INTENT(IN) :: dGASCOLUMNS_dV
    real(fpk),    dimension ( maxlayers, maxgases, 2) , INTENT(IN) :: dGASCOLUMNS_dS
    real(fpk),    dimension ( maxgases)               , INTENT(IN) :: dGASCOLUMN_dP

!  Air density 
    !  Layer column amounts
    !  Air-column derivative with respect to T-shift
    !  Air-column derivative with respect to Surface Pressure
    !  Air columns are in mol/cm^2  or inverse [DU]

    real(fpk),    dimension ( maxlayers ) , INTENT(IN) :: aircolumns
    real(fpk),    dimension ( maxlayers ) , INTENT(IN) :: dAIRCOLUMNS_dS
    real(fpk),                              INTENT(IN) :: dAIRCOLUMN_dP

!  Tshift and Surface pressure

    real(fpk),    INTENT(IN) :: Tshift
    real(fpk),    INTENT(IN) :: Surfpress

!  Gas cross-sections, all levels, all gases;  Rayleigh stuff
    !	Units should be consistent, e.g:
    !	X-sections  are in cm^2/mol  or [DU]

    real(fpk),    dimension ( maxwav, 0:maxlayers, maxgases )      , INTENT(IN) :: gas_xsecs
    real(fpk),    dimension ( maxwav, 0:maxlayers, maxgases )      , INTENT(IN) :: dgasxsecs_dS
    real(fpk),    dimension ( maxwav, 0:maxlayers, maxgases )      , INTENT(IN) :: dgasxsecs_dP
    real(fpk),    dimension ( maxwav )   , INTENT(IN) :: Rayleigh_xsecs
    real(fpk),    dimension ( maxwav )   , INTENT(IN) :: Rayleigh_depols

!  Aerosol optical properties

    LOGICAL,      DIMENSION( maxlayers )              , INTENT(IN) :: AEROSOL_LAYERFLAGS 
    REAL(fpk),    DIMENSION( maxlayers, maxwav )      , INTENT(IN) :: AEROSOL_DELTAU
    REAL(fpk),    DIMENSION( maxwav )                 , INTENT(IN) :: AEROSOL_SSALBS 
    REAL(fpk),    DIMENSION( 6, 0:maxaermoms, maxwav ), INTENT(IN) :: AEROSOL_SCATMOMS 

!  Cloud optical properties

    LOGICAL,      DIMENSION( maxlayers )              , INTENT(IN) :: CLOUD_LAYERFLAGS 
    REAL(fpk),    DIMENSION( maxlayers, maxwav )      , INTENT(IN) :: CLOUD_DELTAU
    REAL(fpk),    DIMENSION( maxwav )                 , INTENT(IN) :: CLOUD_SSALBS 
    REAL(fpk),    DIMENSION( 6, 0:maxcldmoms, maxwav ), INTENT(IN) :: CLOUD_SCATMOMS 

!  VLIDORT OUTPUTS
!  ===============

!  Regular

    INTEGER, intent(OUT)      :: NGREEK_MOMENTS_INPUT

    REAL(fpk),    intent(out) :: DELTAU_VERT_INPUT     ( MAXLAYERS )
    REAL(fpk),    intent(out) :: OMEGA_TOTAL_INPUT     ( MAXLAYERS )
    REAL(fpk),    intent(out) :: GREEKMAT_TOTAL_INPUT  ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  Linearized

    REAL(fpk),    intent(out) :: L_DELTAU_VERT_INPUT    ( MAX_ATMOSWFS, MAXLAYERS )
    REAL(fpk),    intent(out) :: L_OMEGA_TOTAL_INPUT    ( MAX_ATMOSWFS, MAXLAYERS )
    REAL(fpk),    intent(out) :: L_GREEKMAT_TOTAL_INPUT ( MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )

!  Additional output (for Layers-to-Levels Transformation)
!  =================

!  Layer gas optical depths

    REAL(fpk),    intent(out) :: GASOPDEPS       ( maxlayers, maxgases )

!  Chain-Rule Transformation from Layer-GASOPDEP to Level-VMR Jacobians

    REAL(fpk),    intent(out) :: dGASOPDEPS_dV       ( maxlayers, maxgases, 2 )

!  Local variables
!  ===============

    logical      :: AERFLAG_N, CLDFLAG_N
    integer      :: N, N1, L, Q, G
    integer      :: k, ck, gk, cmask(8), gmask(8), smask(8), ngreekmat_entries
    data gmask /  1, 2, 5, 6, 11, 12, 15, 16 /
    data cmask /  1, 5, 5, 2, 3, 6, 6, 4 /
    data smask /  1, -1, -1, 1, 1, 1, -1, 1 /

    REAL(fpk)    :: RAY_ATMOS, GAS_ATMOS, AER_ATMOS, CLD_ATMOS
    REAL(fpk)    :: RAY_SCDEP, GAS_OPDEP, MOL_OPDEP
    REAL(fpk)    :: AER_OPDEP, AER_SCDEP
    REAL(fpk)    :: CLD_OPDEP, CLD_SCDEP
    REAL(fpk)    :: TOT_OPDEP, TOT_SCDEP
    REAL(fpk)    :: RAY_WT, AER_WT, CLD_WT, OMEGA
    REAL(fpk)    :: GAS_ABS, AER_SSALB, CLD_SSALB
    REAL(fpk)    :: PROBLEM_RAY(6,0:2), DEPOL, SK, SIG1, SIG2, RATIO, BETA_NL, FACTOR

    REAL(fpk)    :: dGASOPDEP_dS, dGASABS_dS, dMOLOPDEP_dS, dRAYSCDEP_dS 
    REAL(fpk)    :: dDELTA_dS, dOMEGA_dS, dBETANL_dS, dRAYWT_dS
    REAL(fpk)    :: dGASOPDEP_dP, dGASABS_dP, dMOLOPDEP_dP, dRAYSCDEP_dP 
    REAL(fpk)    :: dDELTA_dP, dOMEGA_dP, dBETANL_dP, dRAYWT_dP

!  Initialize Regular VLIDORT IOPS

    NGREEK_MOMENTS_INPUT = 0
    DELTAU_VERT_INPUT    = 0.0d0
    OMEGA_TOTAL_INPUT    = 0.0d0
    GREEKMAT_TOTAL_INPUT = 0.0d0
    GASOPDEPS      = 0.0d0

!  Initialize Linearized VLIDORT IOPS

    L_DELTAU_VERT_INPUT    = 0.0d0
    L_OMEGA_TOTAL_INPUT    = 0.0d0
    L_GREEKMAT_TOTAL_INPUT = 0.0d0

!  Set some variables

    NGREEK_MOMENTS_INPUT = 2
    if ( do_aerosols ) then
       NGREEK_MOMENTS_INPUT = aerosol_nscatmoms
    endif
    if ( do_clouds_v2 ) then
       NGREEK_MOMENTS_INPUT = max(NGREEK_MOMENTS_INPUT,cloud_nscatmoms)
    endif

!  Set masking limits

    if ( nstokes .eq. 1 ) ngreekmat_entries = 1
    if ( nstokes .eq. 3 ) ngreekmat_entries = 5
    if ( nstokes .eq. 4 ) ngreekmat_entries = 8

!  get the optical properties
!  ==========================

!  Rayleigh scattering law

   depol       = Rayleigh_depols(w)
   Problem_Ray = 0.0d0

   PROBLEM_RAY(1,0) =                    1.D0
   PROBLEM_RAY(4,1) =      3.D0 * ( 1.D0 - 2.D0*DEPOL ) / (2.D0 + DEPOL )
!   PROBLEM_RAY(3,1) =      3.D0 * ( 1.D0 - 2.D0*DEPOL ) / (2.D0 + DEPOL )
   PROBLEM_RAY(1,2) =                  ( 1.D0 - DEPOL ) / (2.D0 + DEPOL )
   PROBLEM_RAY(5,2) =  - DSQRT(6.D0) * ( 1.D0 - DEPOL ) / (2.D0 + DEPOL )
   PROBLEM_RAY(2,2) =           6.D0 * ( 1.D0 - DEPOL ) / (2.D0 + DEPOL )

!  Totals

   RAY_ATMOS = 0.0d0
   GAS_ATMOS = 0.0d0
   AER_ATMOS = 0.0d0
   CLD_ATMOS = 0.0d0

!  Start layer loop

   DO N = 1, NLAYERS

      N1 = N - 1

!  Rayeigh scattering optical depth

      RAY_SCDEP = RAYLEIGH_XSECS(W) * AIRCOLUMNS(N)
      RAY_ATMOS = RAY_ATMOS + RAY_SCDEP

!  Trace gases optical depths
!    Use Cross-sections defined at Levels.

      GAS_OPDEP = 0.0D0
      DO G = 1, NGASES
         GAS_ABS = 0.0d0
         IF (DO_GASES(G)) THEN
            sig1  = gas_xsecs(w,n1,g) ; sig2  = gas_xsecs(w,n,g) 
            GAS_ABS  = GASCOLUMNS(N,G,1) * sig1 + GASCOLUMNS(N,G,2) * sig2
            dGASOPDEPS_dV(N,G,1) =  dGASCOLUMNS_dV(N,G,1) * SIG1
            dGASOPDEPS_dV(N,G,2) =  dGASCOLUMNS_dV(N,G,2) * SIG2
            dGASOPDEPS_dV(N,G,1) =  GASCOLUMNS(N,G,1) * SIG1
            dGASOPDEPS_dV(N,G,2) =  GASCOLUMNS(N,G,2) * SIG2
         endif
         GASOPDEPS(N,G) = GAS_ABS
         GAS_OPDEP = GAS_OPDEP + GAS_ABS
      ENDDO
      GAS_ATMOS = GAS_ATMOS + GAS_OPDEP

!  Total molecular optical depth for layer

      MOL_OPDEP =  RAY_SCDEP + GAS_OPDEP

!  Total optical depth for layer
!  Also sets the scattering optical depths

      TOT_OPDEP = MOL_OPDEP
      TOT_SCDEP = RAY_SCDEP

!  Layer flags

      AERFLAG_N     = DO_AEROSOLS  .AND. AEROSOL_LAYERFLAGS(N)
      CLDFLAG_N     = DO_CLOUDS_V2 .AND. CLOUD_LAYERFLAGS(N)

!  Add aerosols if present in this layer

      AER_OPDEP = 0.0d0
      AER_SCDEP = 0.0d0
      IF ( AERFLAG_N ) THEN
         AER_OPDEP = AEROSOL_DELTAU(N,W)
         AER_SSALB = AEROSOL_SSALBS(W)
         AER_SCDEP = AER_OPDEP * AER_SSALB
         TOT_OPDEP = AER_OPDEP + TOT_OPDEP
         TOT_SCDEP = AER_SCDEP + TOT_SCDEP
      ENDIF
      AER_ATMOS = AER_ATMOS + AER_OPDEP

!  Add clouds if present in this layer
!    Only for cloud method 2

      CLD_OPDEP = 0.0d0
      CLD_SCDEP = 0.0d0
      IF ( CLDFLAG_N ) THEN
         CLD_OPDEP = CLOUD_DELTAU(N,W)
         CLD_SSALB = CLOUD_SSALBS(W)
         CLD_SCDEP = CLD_OPDEP * CLD_SSALB
         TOT_OPDEP = CLD_OPDEP + TOT_OPDEP
         TOT_SCDEP = CLD_SCDEP + TOT_SCDEP
      ENDIF
      CLD_ATMOS = CLD_ATMOS + CLD_OPDEP

!      write(*,'(i3,5f12.6)')n,ray_scdep,gas_opdep,CLD_OPDEP,CLD_SCDEP,tot_opdep

!  Scatter weighting

      RAY_WT                  = RAY_SCDEP / TOT_SCDEP
      IF ( AERFLAG_N ) AER_WT = AER_SCDEP / TOT_SCDEP
      IF ( CLDFLAG_N ) CLD_WT = CLD_SCDEP / TOT_SCDEP

!  VLIDORT Bulk properties. Use limiting value for OMEGA

      DELTAU_VERT_INPUT(N) = TOT_OPDEP
      OMEGA                = TOT_SCDEP / TOT_OPDEP
      IF ( OMEGA .GT. OMEGA_LIM ) THEN
         OMEGA_TOTAL_INPUT(N) = OMEGA_LIM
      ELSE
         OMEGA_TOTAL_INPUT(N) = TOT_SCDEP / TOT_OPDEP
      ENDIF

!  SCATTERING LAW moments, Rayleigh first

      DO K = 1, ngreekmat_entries
         GK = gmask(k)
         ck = cmask(k)
         DO L = 0, 2
            GREEKMAT_TOTAL_INPUT(L,n,gk) = RAY_WT * PROBLEM_RAY(ck,L)
         ENDDO
      ENDDO

!  Add aerosols

      IF ( AERFLAG_N ) THEN
         DO K = 1, ngreekmat_entries
            GK = gmask(k)
            sk = dble(smask(k))
            ck = cmask(k)
            DO L = 0, NGREEK_MOMENTS_INPUT
               GREEKMAT_TOTAL_INPUT(L,n,gk) = GREEKMAT_TOTAL_INPUT(L,n,gk)  &
                                            + AER_WT * SK * AEROSOL_SCATMOMS(ck,L,W)
            ENDDO
         ENDDO
      ENDIF

!  Add clouds

      IF ( CLDFLAG_N ) THEN
         DO K = 1, ngreekmat_entries
            GK = gmask(k)
            sk = dble(smask(k))
            ck = cmask(k)
            DO L = 0, NGREEK_MOMENTS_INPUT
               GREEKMAT_TOTAL_INPUT(L,n,gk) = GREEKMAT_TOTAL_INPUT(L,n,gk)  &
                                            + CLD_WT * SK * CLOUD_SCATMOMS(ck,L,W)
!            if (n.eq.29)write(*,*)n,l, GREEKMAT_TOTAL_INPUT(L,n,gk)
            ENDDO
         ENDDO
      ENDIF

!  Phase function normalization

      GREEKMAT_TOTAL_INPUT(0,n,1) = 1.0d0

!  Linearized VLIDORT properties
!  =============================

!  Profile Jacobians for Trace gases
!  ---------------------------------

      if ( do_Profile_Jacobians ) then
         Q = 0
         DO G = 1, NGASES
            IF ( DO_GASES(G).AND.DO_GAS_WFS(G) ) THEN
               Q = Q + 1
               RATIO = GASOPDEPS(N,G)  / TOT_OPDEP
               L_DELTAU_VERT_INPUT(Q,N) =   RATIO
               L_OMEGA_TOTAL_INPUT(Q,N) = - RATIO
            ENDIF
         ENDDO
      endif

!  Profile Jacobians for Aerosols
!  ------------------------------

      if ( do_AerOpdepProfile_Jacobians ) then
         if ( AERFLAG_N ) THEN
            Q = Q + 1
            RATIO  = AER_OPDEP  / TOT_OPDEP
            FACTOR = ( AER_SSALB / OMEGA ) - 1.0d0
            L_DELTAU_VERT_INPUT(Q,N) = RATIO
            L_OMEGA_TOTAL_INPUT(Q,N) = RATIO * FACTOR
            DO K = 1, ngreekmat_entries
               GK = gmask(k)
               sk = dble(smask(k))
               ck = cmask(k)
               DO L = 0, NGREEK_MOMENTS_INPUT
                  BETA_NL = GREEKMAT_TOTAL_INPUT(L,n,gk)
                  if ( BETA_NL .ne. 0.0d0 ) then
                     FACTOR = ( SK * AEROSOL_SCATMOMS(ck,L,W) / BETA_NL ) - 1.0d0
                     L_GREEKMAT_TOTAL_INPUT(Q,L,n,gk) = AER_WT * FACTOR
!                  if ( l.eq.10)write(*,*)n,q, BETA_NL * L_GREEKMAT_TOTAL_INPUT(Q,L,N,gk), BETA_NL    
                  endif
               ENDDO
               L_GREEKMAT_TOTAL_INPUT(Q,0,n,1) = 0.0d0   !  comes from phase function normalization
            ENDDO
         ENDIF
      endif

!  Column Jacobians: T-shift and  surface-pressue Jacobians
!  --------------------------------------------------------

      Q = 0

      if ( do_Tshift_Jacobian ) then

!  First column Jacobian

         Q = 1

!  Derivative of Rayleigh scattering optical depth

         dRAYSCDEP_dS = RAYLEIGH_XSECS(W) * dAIRCOLUMNS_dS(N)

!  Derivative of Absoprtion optical depth

         dGASOPDEP_dS = 0.0d0
         DO G = 1, NGASES
            dGASABS_dS = 0.0d0
            IF (DO_GASES(G)) THEN
               sig1  = gas_xsecs(w,n1,g) ; sig2  = gas_xsecs(w,n,g) 
               dGASABS_dS  =   dGASCOLUMNS_dS(N,G,1) * sig1 + GASCOLUMNS(N,G,1) * dgasxsecs_dS(w,n1,g) &
                             + dGASCOLUMNS_dS(N,G,2) * sig2 + GASCOLUMNS(N,G,2) * dgasxsecs_dS(w,n,g) 
            endif
            dGASOPDEP_dS = dGASOPDEP_dS + dGASABS_dS
         ENDDO

!  Derivative of Extinction optical depth

         dMOLOPDEP_dS =  dRAYSCDEP_dS + dGASOPDEP_dS

!  Derivative of DELTA and OMEGA

         dDELTA_dS = dMOLOPDEP_dS
         dOMEGA_dS = ( dRAYSCDEP_dS - OMEGA * dMOLOPDEP_dS ) / TOT_OPDEP

!  Normalized derivatives == VLIDORT inputs

         L_DELTAU_VERT_INPUT(Q,N) = TSHIFT * dDELTA_dS / TOT_OPDEP
         if ( OMEGA_TOTAL_INPUT(N) .ge. OMEGA_LIM ) then
            L_OMEGA_TOTAL_INPUT(Q,N) = 0.0_fpk
         else
            L_OMEGA_TOTAL_INPUT(Q,N) = TSHIFT * dOMEGA_dS / OMEGA   
         endif
         
!  First debug

!         write(998,*)n, L_DELTAU_VERT_INPUT(Q,N) * TOT_OPDEP, DELTAU_VERT_INPUT(N)

!  Derivative Greek moments, only for non-Rayleigh layers

         if ( CLDFLAG_N .or. AERFLAG_N ) THEN
            dRAYWT_dS  = dRAYSCDEP_dS / TOT_SCDEP
            DO K = 1, ngreekmat_entries
               GK = gmask(k) ; ck = cmask(k)
               DO L = 0, 2
                  BETA_NL = GREEKMAT_TOTAL_INPUT(L,n,gk)
                  if ( BETA_NL .ne. 0.0d0 ) then
                     dBETANL_dS = dRAYWT_dS * ( PROBLEM_RAY(ck,L) - BETA_NL )
                     L_GREEKMAT_TOTAL_INPUT(Q,L,n,gk) = TSHIFT * dBETANL_dS / BETA_NL
                  endif
               ENDDO
               DO L = 3, NGREEK_MOMENTS_INPUT
                  BETA_NL = GREEKMAT_TOTAL_INPUT(L,n,gk)
                  if ( BETA_NL .ne. 0.0d0 ) then
!                     dBETANL_dS = - dRAYWT_dS * BETA_NL
!                     L_GREEKMAT_TOTAL_INPUT(Q,L,n,gk) = TSHIFT * dBETANL_dS / BETA_NL
                     L_GREEKMAT_TOTAL_INPUT(Q,L,n,gk) = - TSHIFT * dRAYWT_dS
                  endif
               ENDDO
            ENDDO
         endif

!  End T-shift Jacobian

      endif

!  surface Pressure Jacobian

      if ( do_SurfPress_Jacobian .and. n.eq.nlayers ) then

!  First or second column Jacobian

         Q = Q + 1

!        write(*,*)dAIRCOLUMN_dP*SURFPRESS, AIRCOLUMNS(N)


!  Derivative of Absorption, Rayleigh scattering and Total molecular optical depth
!  Dont forget the Cross-section derivatives

         dRAYSCDEP_dP = RAYLEIGH_XSECS(W) * dAIRCOLUMN_dP
         dGASOPDEP_dP = 0.0d0
         DO G = 1, NGASES
            dGASABS_dP = 0.0d0
            IF (DO_GASES(G)) THEN
               dGASABS_dP  = dGASCOLUMN_dP(G) * gas_xsecs(w,n,g) + GASCOLUMNS(N,G,2) * dgasxsecs_dP(w,n,g)
            endif
            dGASOPDEP_dP = dGASOPDEP_dP + dGASABS_dP
         ENDDO
         dMOLOPDEP_dP = dRAYSCDEP_dP + dGASOPDEP_dP

!  Derivative of DELTA and OMEGA

         dDELTA_dP = dMOLOPDEP_dP
         dOMEGA_dP = ( dRAYSCDEP_dP - OMEGA * dMOLOPDEP_dP ) / TOT_OPDEP

!  Normalized derivatives == VLIDORT inputs

         L_DELTAU_VERT_INPUT(Q,N) = SURFPRESS * dDELTA_dP / TOT_OPDEP
         if ( OMEGA_TOTAL_INPUT(N) .ge. OMEGA_LIM ) then
            L_OMEGA_TOTAL_INPUT(Q,N) = 0.0_fpk
         else
            L_OMEGA_TOTAL_INPUT(Q,N) = SURFPRESS * dOMEGA_dP / OMEGA   
         endif

!  Derivative Greek moments, only for non-Rayleigh layers

         if ( CLDFLAG_N .or. AERFLAG_N ) THEN
            dRAYWT_dP  = dRAYSCDEP_dP / TOT_SCDEP
            DO K = 1, ngreekmat_entries
               GK = gmask(k) ; ck = cmask(k)
               DO L = 0, 2
                  BETA_NL = GREEKMAT_TOTAL_INPUT(L,n,gk)
                  if ( BETA_NL .ne. 0.0d0 ) then
                     dBETANL_dP = dRAYWT_dP * ( PROBLEM_RAY(ck,L) - BETA_NL )
                     L_GREEKMAT_TOTAL_INPUT(Q,L,n,gk) = SURFPRESS * dBETANL_dP / BETA_NL
                  endif
               ENDDO
               DO L = 3, NGREEK_MOMENTS_INPUT
                  BETA_NL = GREEKMAT_TOTAL_INPUT(L,n,gk)
                  if ( BETA_NL .ne. 0.0d0 ) then
                     L_GREEKMAT_TOTAL_INPUT(Q,L,n,gk) = - SURFPRESS * dRAYWT_dP
                  endif
               ENDDO
            ENDDO
         endif

!  End Surface-Pressure Jacobian

      endif

!  End layer loop

   ENDDO

!   pause'Checking L_DELTAU etc, Aero pRof'

!  Debug Output: Total optical depths (fort.33), Trace gas optical depths (fort.34)

   if ( do_clouds_V2) then
       write(33,567)w,wavnum,GAS_ATMOS,RAY_ATMOS,AER_ATMOS, CLD_ATMOS
   else
       write(33,567)w,wavnum,GAS_ATMOS,RAY_ATMOS,AER_ATMOS, 0.0d0
   endif
   write(34,567)w,wavnum,sum(GASOPDEPS(1:NLAYERS,1)),sum(GASOPDEPS(1:NLAYERS,2)),sum(GASOPDEPS(1:NLAYERS,3)),&
                         sum(GASOPDEPS(1:NLAYERS,4))
567 format(i5,f12.5,1p8e17.7)


!  DEBUG output of optical properties

   if ( IOP_debug_flag ) then
      if ( w.eq.1 ) then
         open(578,file=adjustl(trim(IOP_debug_filename)),status='replace')
!          open(579,file='fort.579',status='replace')
!          do n = 1, nlayers
!             write(579,'(i5,1p3e20.10)')n,height_grid(n), DELTAU_VERT_INPUT(n),OMEGA_TOTAL_INPUT(n)
!          enddo
!          close(579)
       endif
       if ( w.eq.1 ) then
          write(578,*)w
          do n = 1, nlayers
             write(578,'(i5,1p2e20.10)')n, DELTAU_VERT_INPUT(n), OMEGA_TOTAL_INPUT(n)
             if ( aerosol_layerflags(n) ) then
                do l = 0, ngreek_moments_input
                   write(578,'(i5,1p11e11.4)')l, (GREEKMAT_TOTAL_INPUT(L,n,k),k=1,11)
                enddo
             else
                do l = 0, 2
                   write(578,'(i5,1p11e11.4)')l, (GREEKMAT_TOTAL_INPUT(L,n,k),k=1,11)
                enddo
             endif
          enddo
       endif
       if ( w.eq.nwav) close(578)
 !      if ( w.eq.nwav) pause 'checker 1'
    endif

!    if (w.eq.4)stop'tempo stop @@@@@@@@@@@@@@@@@'

!  Finish

   RETURN
END subroutine GEMSTOOL_NSW_linearized_iops


!  End module

End  module GEMSTOOL_NSW_linearized_iops_m


