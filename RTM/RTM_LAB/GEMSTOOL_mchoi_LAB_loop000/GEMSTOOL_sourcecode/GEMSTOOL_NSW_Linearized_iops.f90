module GEMSTOOL_NSW_linearized_iops_m

!  routines public

!  @@ Rob fix 7/23/14  , add H2O scaling derivative
!  @@ Y.Jung fix 2/1/15, add CH4 scaling derivative
!  @@ Rob 10/20/16. n_scatmoms now has wavelength dependence

private
public :: GEMSTOOL_NSW_linearized_iops

contains

subroutine GEMSTOOL_NSW_linearized_iops &
     ( maxlayers, maxmoments_input, maxstokes_sq, max_atmoswfs,                  & ! VLIDORT dimensions
       maxwav, maxaermoms, maxcldmoms, maxgases, maxaerwfs,                      & ! GEMSTOOL dimensions
       IOP_debug_filename, IOP_Debug_flag,                                       & ! Debug IOP check 
       do_Profile_Jacobians, do_AerOpdepProfile_Jacobians,                       & ! Profile Jacobian control flags
       do_AerBulk_Jacobians, do_Tshift_Jacobian,                                 & ! Column/Bulk Jacobian control flags
       do_Surfpress_Jacobian, do_H2OScaling_Jacobian, do_CH4Scaling_Jacobian,    & ! Column/Bulk Jacobian control flags
       w, wavnum, omega_lim, do_gases, do_gas_wfs, ngases,                       & ! GEMSTOOL Inputs (Control)
       nstokes, nlayers, nwav, n_aerosol_wfs, aerosol_pars,                      & ! GEMSTOOL Inputs (Control)
       do_aerosols, aerosol_nscatmoms, do_clouds_V2, cloud_nscatmoms,            & ! GEMSTOOL Inputs (Aer/clouds)
       gascolumns, dGASCOLUMNS_dS, dGASCOLUMN_dP, dGASCOLUMNS_dV,                & ! GEMSTOOL Inputs (GasColumns + derivs)
       dH2OCOLUMN_dF, dCH4COLUMN_dF,                                             & ! GEMSTOOL Inputs (GasColumns + derivs)
       aircolumns, dAIRCOLUMNS_dS, dAIRCOLUMN_dP,                                & ! GEMSTOOL Inputs (AirColumns + derivs)
       Tshift, SurfPress, H2OScaling, gh2o, CH4Scaling, gch4,                    & ! GEMSTOOL Inputs (T-shift,SurfP,H2O+CH4Scaling)
       gas_xsecs, dgasxsecs_dS, dgasxsecs_dP, Rayleigh_xsecs, Rayleigh_depols,   & ! GEMSTOOL Inputs (Cross-sections, depol)
       aerosol_layerflags, aerosol_deltau, aerosol_ssalbs, aerosol_scatmoms,     & ! GEMSTOOL Inputs (Aerosol)
       L_aerosol_deltau, L_aerosol_ssalbs, L_aerosol_scatmoms,                   & ! GEMSTOOL Inputs (Aerosol linearized)
       cloud_layerflags, cloud_deltau, cloud_ssalbs, cloud_scatmoms,             & ! GEMSTOOL Inputs (Clouds)
       DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT,                                     & ! VLIDORT Outputs (Regular iops)
       NGREEK_MOMENTS_INPUT, GREEKMAT_TOTAL_INPUT,                               & ! VLIDORT Outputs (Regular iops)
       L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,         & ! VLIDORT Outputs (linearized)
       GASOPDEPS, dGASOPDEPS_dV )                                                  ! Layer-to-Level transformation

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

   integer, INTENT(IN)  ::  maxwav, maxaermoms, maxcldmoms, maxgases, maxaerwfs

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
!  @@ Rob fix 7/23/14, add H2O scaling derivative flag
!  @@ Y.Jung fix 1/2/15, add CH4 scaling derivative flag

    logical      , intent(in) :: do_Profile_Jacobians
    logical      , intent(in) :: do_AerOpdepProfile_Jacobians
    logical      , intent(in) :: do_AerBulk_Jacobians

    logical      , intent(in) :: do_Tshift_Jacobian
    logical      , intent(in) :: do_SurfPress_Jacobian

    logical      , intent(in) :: do_H2OScaling_Jacobian
    logical      , intent(in) :: do_CH4Scaling_Jacobian

!  Gas control

    logical    , dimension ( maxgases), INTENT(IN)  :: do_gases
    logical    , dimension ( maxgases), INTENT(IN)  :: do_gas_wfs

!  Aerosol control

    logical, INTENT(IN)  :: do_aerosols

!  Cloud control

    logical, INTENT(IN)  :: do_clouds_V2

!  number of gases, aerosol/cloud scattering moments

    integer, INTENT(IN)  :: ngases
    integer, INTENT(IN)  :: aerosol_nscatmoms (maxwav)
    integer, INTENT(IN)  :: cloud_nscatmoms

!  Number of aerosol weighting functions and parameter values
!   -- 10/20/16, number of aerosol moments now has wavelength dependence

    real(fpk), INTENT(IN) :: aerosol_pars (maxaerwfs)
    integer  , INTENT(IN) :: n_aerosol_wfs

!  Trace gas stuff    !  @@ Rob fix 7/23/14, add H2O scaling derivative
!  Trace gas stuff    !  @@ Y.Jung fix 1/2/15, add CH4 scaling derivative
    !  Trace gas Layer column amounts
    !  Gas-column derivative with respect to Level VMRs. (needed for LAYERS-to_LEVELS Transformation)
    !  Gas-column derivative with respect to T-shift
    !  Gas-column derivative with respect to Surface Pressure
    !  H2O-column derivative with respect to Scaling
    !  CH4-column derivative with respect to Scaling
    !	LAYER gas columns are in mol/cm^2  or inverse [DU]

    real(fpk),    dimension ( maxlayers, maxgases, 2) , INTENT(IN) :: gascolumns
    real(fpk),    dimension ( maxlayers, maxgases, 2) , INTENT(IN) :: dGASCOLUMNS_dV
    real(fpk),    dimension ( maxlayers, maxgases, 2) , INTENT(IN) :: dGASCOLUMNS_dS
    real(fpk),    dimension ( maxgases)               , INTENT(IN) :: dGASCOLUMN_dP
    REAL(fpk)   , dimension ( maxlayers, 2 )          , INTENT(IN) :: dH2OCOLUMN_dF
    REAL(fpk)   , dimension ( maxlayers, 2 )          , INTENT(IN) :: dCH4COLUMN_dF

!  Air density 
    !  Layer column amounts
    !  Air-column derivative with respect to T-shift
    !  Air-column derivative with respect to Surface Pressure
    !  Air columns are in mol/cm^2  or inverse [DU]

    real(fpk),    dimension ( maxlayers ) , INTENT(IN) :: aircolumns
    real(fpk),    dimension ( maxlayers ) , INTENT(IN) :: dAIRCOLUMNS_dS
    real(fpk),                              INTENT(IN) :: dAIRCOLUMN_dP

!  Tshift and Surface pressure
 !  @@ Rob fix 7/23/14, add H2O scaling
 !  @@ Y.Jung fix 2/1/15, add CH4 scaling

    real(fpk),    INTENT(IN) :: Tshift
    real(fpk),    INTENT(IN) :: Surfpress
    real(fpk),    INTENT(IN) :: H2OScaling
    integer  ,    INTENT(IN) :: gh2o
    real(fpk),    INTENT(IN) :: CH4Scaling
    integer  ,    INTENT(IN) :: gch4

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

!  Aerosol Linearized optical properties

    REAL(fpk),    DIMENSION ( maxlayers, maxwav, maxaerwfs )      , INTENT(IN) :: L_AEROSOL_DELTAU
    REAL(fpk),    DIMENSION ( maxwav, maxaerwfs )                 , INTENT(IN) :: L_AEROSOL_SSALBS 
    REAL(fpk),    DIMENSION ( 6, 0:maxaermoms, maxwav, maxaerwfs ), INTENT(IN) :: L_AEROSOL_SCATMOMS 

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
    integer      :: N, N1, L, Q, Q1, QOFFSET1, G
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
    REAL(fpk)    :: L_AER_OPDEP, L_AER_SCDEP(maxaerwfs), L_AER_SSALB, T1, T2

    REAL(fpk)    :: dGASOPDEP_dS, dGASABS_dS, dMOLOPDEP_dS, dRAYSCDEP_dS 
    REAL(fpk)    :: dDELTA_dS, dOMEGA_dS, dBETANL_dS, dRAYWT_dS
    REAL(fpk)    :: dGASOPDEP_dP, dGASABS_dP, dMOLOPDEP_dP, dRAYSCDEP_dP 
    REAL(fpk)    :: dDELTA_dP, dOMEGA_dP, dBETANL_dP, dRAYWT_dP
    REAL(fpk)    :: dDELTA_dF, dOMEGA_dF, dMOLOPDEP_dF, dGASOPDEP_dF

!  Initialize Regular VLIDORT IOPS

    NGREEK_MOMENTS_INPUT = 0
    DELTAU_VERT_INPUT    = 0.0d0
    OMEGA_TOTAL_INPUT    = 0.0d0
    GREEKMAT_TOTAL_INPUT = 0.0d0
    GASOPDEPS            = 0.0d0

!  Initialize Linearized VLIDORT IOPS
!   ( Bookkeeping for Jacobians is done outside this routine )

    L_DELTAU_VERT_INPUT    = 0.0d0
    L_OMEGA_TOTAL_INPUT    = 0.0d0
    L_GREEKMAT_TOTAL_INPUT = 0.0d0

!  Set some variables

    NGREEK_MOMENTS_INPUT = 2
    if ( do_aerosols ) then
       NGREEK_MOMENTS_INPUT = aerosol_nscatmoms(w)   ! 10/20/16 Wavelength dependence
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
            if (sig1 .ne. sig1) sig1=0.0d0 !mchoi nan control
            if (sig2 .ne. sig2) sig2=0.0d0 !mchoi nan control
 
            GAS_ABS  = GASCOLUMNS(N,G,1) * sig1 + GASCOLUMNS(N,G,2) * sig2
            dGASOPDEPS_dV(N,G,1) =  dGASCOLUMNS_dV(N,G,1) * SIG1
            dGASOPDEPS_dV(N,G,2) =  dGASCOLUMNS_dV(N,G,2) * SIG2
!            if (G.eq.1) write(31,*) w,n,g, gascolumns(N,G,1), GASCOLUMNS(N,G,2), sig1, sig2, GAS_ABS
!            if (G.eq.2) write(32,*) w,n,g, gascolumns(N,G,1), GASCOLUMNS(N,G,2), sig1, sig2, GAS_ABS
!            if (G.eq.3) write(33,*) w,n,g, gascolumns(N,G,1), GASCOLUMNS(N,G,2), sig1, sig2, GAS_ABS
           ! dGASOPDEPS_dV(N,G,1) =  GASCOLUMNS(N,G,1) * SIG1
           ! dGASOPDEPS_dV(N,G,2) =  GASCOLUMNS(N,G,2) * SIG2
         endif
         GASOPDEPS(N,G) = GAS_ABS
         GAS_OPDEP = GAS_OPDEP + GAS_ABS
      ENDDO
!      write(34,*) W, N, GAS_OPDEP
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

!  Always zero the count, for every layer
 
        Q = 0

!  Profile Jacobians for Trace gases
!  ---------------------------------

      if ( do_Profile_Jacobians ) then
!         Q = 0
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

!  Column Jacobians: H2O-Scaling, CH4-Scaling, T-shift and  surface-pressue Jacobians
!  =====================================================================

!  @@ Rob fix 7/23/14, add H2O scaling derivative

!  Always zero the count, for every layer

      Q = 0

!  H2O-Scaling Jacobian
!  --------------------

      if ( do_H2OScaling_Jacobian ) then

!  First column Jacobian

         Q = Q + 1

!  Derivative of Absorption and Extinction optical depths

         dGASOPDEP_dF =  dH2OCOLUMN_dF(N,1) * gas_xsecs(w,n1,gh2o) + dH2OCOLUMN_dF(N,2) * gas_xsecs(w,n,gh2o)
         dMOLOPDEP_dF =  dGASOPDEP_dF

!  Derivative of DELTA and OMEGA

         dDELTA_dF = dMOLOPDEP_dF
         dOMEGA_dF =  - OMEGA * dMOLOPDEP_dF / TOT_OPDEP

!  Normalized derivatives == VLIDORT inputs

         L_DELTAU_VERT_INPUT(Q,N) = H2OScaling * dDELTA_dF / TOT_OPDEP
         if ( OMEGA_TOTAL_INPUT(N) .ge. OMEGA_LIM ) then
            L_OMEGA_TOTAL_INPUT(Q,N) = 0.0_fpk
         else
            L_OMEGA_TOTAL_INPUT(Q,N) = H2OScaling * dOMEGA_dF / OMEGA   
         endif

!         write(35,*)W,Q,N,L_DELTAU_VERT_INPUT(Q,N)*DELTAU_VERT_INPUT(N),DELTAU_VERT_INPUT(N)
!         write(37,*)W,Q,N,L_OMEGA_TOTAL_INPUT(Q,N)*OMEGA_TOTAL_INPUT(N),OMEGA_TOTAL_INPUT(N)

!  End H2O Scaling Jacobian

      endif


!  @@ Y.Jung fix 2/1/15, add CH4 scaling derivative

!  CH4-Scaling Jacobian
!  --------------------

      if ( do_CH4Scaling_Jacobian ) then

!  Next column Jacobian

         Q = Q + 1

!  Derivative of Absorption and Extinction optical depths

         dGASOPDEP_dF =  dCH4COLUMN_dF(N,1) * gas_xsecs(w,n1,gch4) + dCH4COLUMN_dF(N,2) * gas_xsecs(w,n,gch4)
         dMOLOPDEP_dF =  dGASOPDEP_dF

!  Derivative of DELTA and OMEGA

         dDELTA_dF = dMOLOPDEP_dF
         dOMEGA_dF =  - OMEGA * dMOLOPDEP_dF / TOT_OPDEP

!  Normalized derivatives == VLIDORT inputs

         L_DELTAU_VERT_INPUT(Q,N) = CH4Scaling * dDELTA_dF / TOT_OPDEP
         if ( OMEGA_TOTAL_INPUT(N) .ge. OMEGA_LIM ) then
            L_OMEGA_TOTAL_INPUT(Q,N) = 0.0_fpk
         else
            L_OMEGA_TOTAL_INPUT(Q,N) = CH4Scaling * dOMEGA_dF / OMEGA   
         endif

!         write(35,*)W,Q,N,L_DELTAU_VERT_INPUT(Q,N)*DELTAU_VERT_INPUT(N),DELTAU_VERT_INPUT(N)
!         write(37,*)W,Q,N,L_OMEGA_TOTAL_INPUT(Q,N)*OMEGA_TOTAL_INPUT(N),OMEGA_TOTAL_INPUT(N)

!  End CH4 Scaling Jacobian

      endif


!  T-shift Jacobian
!  ----------------

      if ( do_Tshift_Jacobian ) then

!  Next column Jacobian

         Q = Q + 1

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
         
! debug
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
!  -------------------------
!  Next column Jacobian

!  Rob Fix. 26 March 2015. Must increase weighting function Count Q in all layers
!                          even though only N= NLAYERS has non-zero entries for Surf-pressure Jacobians

!      if ( do_SurfPress_Jacobian .and. n.eq.nlayers ) then
      if ( do_SurfPress_Jacobian ) then

!  First or second column Jacobian

         Q = Q + 1 ; if (n.ne.nlayers ) go to 3556    ! Increase Q-count, but skip if n /= nlayers

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


!  Continuation point for avoiding the calculation for NLAYERS

  3556  continue

!  End Surface-Pressure Jacobian

      endif

!  Offset so far for Bulk Jacobians

      qoffset1 = Q

!  Column/Bulk Jacobians for Aerosols (New, June 2014, R. Spurr)
!  -------------------------------------------------------------

      if ( do_AerBulk_Jacobians ) then
         if ( AERFLAG_N ) THEN
            do q = 1, n_aerosol_wfs
              Q1 = qoffset1 + Q
              L_AER_OPDEP    = L_AEROSOL_DELTAU(N,W,Q)
              L_AER_SSALB    = L_AEROSOL_SSALBS(W,Q)
              L_AER_SCDEP(Q) = L_AER_OPDEP * AER_SSALB + L_AER_SSALB * AER_OPDEP
              L_DELTAU_VERT_INPUT(Q1,N) = AEROSOL_PARS(Q) * L_AER_OPDEP / TOT_OPDEP
              L_OMEGA_TOTAL_INPUT(Q1,N) = AEROSOL_PARS(Q) * &
                       ( L_AER_SCDEP(Q) - OMEGA * L_AER_OPDEP ) / TOT_SCDEP
               DO K = 1, ngreekmat_entries
                  GK = gmask(k) ; sk = dble(smask(k)) ; ck = cmask(k)
                  DO L = 0, NGREEK_MOMENTS_INPUT
                     BETA_NL = GREEKMAT_TOTAL_INPUT(L,n,gk)
                     if ( BETA_NL.ne. 0.0d0 ) then
                        RATIO = aerosol_PARS(q) / BETA_NL
                        T1 = L_AER_SCDEP(Q) * &
                          ( SK * AEROSOL_SCATMOMS(ck,L,W) - GREEKMAT_TOTAL_INPUT(L,n,gk) )
                        T2 =   AER_SCDEP    * SK * L_AEROSOL_SCATMOMS(ck,L,W,Q) 
                        L_GREEKMAT_TOTAL_INPUT(Q1,L,n,gk) = RATIO * ( T1 + T2 ) / TOT_SCDEP
                     endif                
                  ENDDO
               ENDDO
               L_GREEKMAT_TOTAL_INPUT(Q,0,n,1) = 0.0d0        !  Normalization
            ENDDO
         ENDIF
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
!   write(34,567)w,wavnum,sum(GASOPDEPS(1:NLAYERS,1)),sum(GASOPDEPS(1:NLAYERS,2)),sum(GASOPDEPS(1:NLAYERS,3)),&
!                         sum(GASOPDEPS(1:NLAYERS,4))
  write(34,567)w,wavnum,sum(GASOPDEPS(1:NLAYERS,1)),sum(GASOPDEPS(1:NLAYERS,2))
567 format(i5,f12.5,1p8e17.7)

!    write(*,567)w,wavnum,sum(GASOPDEPS(1:NLAYERS,1)),sum(GASOPDEPS(1:NLAYERS,2)),sum(GASOPDEPS(1:NLAYERS,3)),&
!    sum(GASOPDEPS(1:NLAYERS,4))
! 567 format(i5,f12.5,1p8e17.7)

!  DEBUG output of optical properties

!   if ( IOP_debug_flag ) then
   !   if ( w.eq.1 ) then
   !      open(578,file=adjustl(trim(IOP_debug_filename)),status='replace')
   !       open(579,file='fort.579',status='replace')
   !       do n = 1, nlayers
   !          write(579,'(i5,1p3e20.10)')n,DELTAU_VERT_INPUT(n),OMEGA_TOTAL_INPUT(n) !,height_grid(n), 
   !       enddo
   !       close(579)
   !    endif
   !    if ( w.eq.1 ) then
   !       write(578,*)w
   !       do n = 1, nlayers
   !          write(578,'(i5,1p2e20.10)')n, DELTAU_VERT_INPUT(n), OMEGA_TOTAL_INPUT(n)
   !          if ( aerosol_layerflags(n) ) then
   !             do l = 0, ngreek_moments_input
   !                write(578,'(i5,1p11e11.4)')l, (GREEKMAT_TOTAL_INPUT(L,n,k),k=1,11)
   !             enddo
   !          else
   !             do l = 0, 2
   !                write(578,'(i5,1p11e11.4)')l, (GREEKMAT_TOTAL_INPUT(L,n,k),k=1,11)
   !             enddo
   !          endif
   !       enddo
   !    endif
   !    if ( w.eq.nwav) close(578)
!!       if ( w.eq.nwav) pause 'checker 1'
!    endif
!    if (w.eq.4)stop'tempo stop @@@@@@@@@@@@@@@@@'

!  Finish

   RETURN
END subroutine GEMSTOOL_NSW_linearized_iops


!  End module

End  module GEMSTOOL_NSW_linearized_iops_m


