module GEMSTOOL_UVN_regular_iops_m

!  routines public

!  @@ Rob 10/25/16. aerosol_n_scatmoms now has wavelength dependence

private
public :: GEMSTOOL_UVN_regular_iops

contains

subroutine GEMSTOOL_UVN_regular_iops                                             &
     ( maxlayers, maxmoments_input, maxstokes_sq,                                & ! VLIDORT dimensions
       maxwav, maxaermoms, maxcldmoms, maxgases,                                 & ! GEMSTOOL dimensions
       IOP_debug_filename, IOP_Debug_flag,                                       & ! Debug IOP check  
       w, wav, omega_lim, do_gases, which_gases, ngases, nstokes, nlayers, nwav, & ! GEMSTOOL Inputs
       do_aerosols, aerosol_nscatmoms, do_clouds_V2, cloud_nscatmoms,        & ! GEMSTOOL Inputs
       level_temps, gas_xsecs, O3C1_xsecs, O3C2_xsecs, H2O_xsecs, O2_xsecs,  & ! GEMSTOOL Inputs
       layer_gascolumns, layer_aircolumns, Rayleigh_xsecs, Rayleigh_depols,  & ! GEMSTOOL Inputs
       aerosol_layerflags, aerosol_deltau, aerosol_ssalbs, aerosol_scatmoms, & ! GEMSTOOL Inputs
       cloud_layerflags, cloud_deltau, cloud_ssalbs, cloud_scatmoms,         & ! GEMSTOOL Inputs
       DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT,                                 & ! VLIDORT Outputs
!       NGREEK_MOMENTS_INPUT, GREEKMAT_TOTAL_INPUT )                            ! VLIDORT Outputs
       NGREEK_MOMENTS_INPUT, GREEKMAT_TOTAL_INPUT,                           &  ! VLIDORT Outputs
       TOT_OPDEP_L,AER_OPDEP_L,AER_SCDEP_L,RAY_SCDEP_L,GAS_OPDEP_L)           !test value  

   implicit none

!  precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  VLIDORT INPUTS
!  ==============

!  Dimensioning

    integer, INTENT(IN)  :: maxlayers, maxmoments_input, maxstokes_sq

!  GEMSTOOL INPUTS
!  ==============

!  GEMSTOOL Dimensioning

   integer, INTENT(IN)  ::  maxwav, maxaermoms, maxcldmoms, maxgases

!  Debug Iop check

   logical      , intent(in) :: IOP_debug_flag
   character*(*), intent(in) :: IOP_debug_filename
  
!  Control

    integer, INTENT(IN)  ::  nstokes, nlayers

!  Wavelength index and value

    integer, INTENT(IN)   :: w, nwav
    real(fpk), INTENT(IN) :: wav

!  Limiting single scatter albedo
    !	This should be pre-set to coincide with the internal limit in
    !	the LIDORT models. Suggested value = 0.999999

    real(fpk),    INTENT(IN) :: omega_lim

!  Gas control

    logical    , dimension ( maxgases), INTENT(IN)  :: do_gases
    character*4, dimension ( maxgases), INTENT(IN)  :: which_gases

!  Aerosol control

    logical, INTENT(IN)  :: do_aerosols

!  Cloud control

    logical, INTENT(IN)  :: do_clouds_V2

!  number of gases, aerosol/cloud scattering moments
!   -- 10/25/16, number of aerosol moments now has wavelength dependence

    integer, INTENT(IN)  :: ngases
    integer, INTENT(IN)  :: aerosol_nscatmoms(maxwav)
    integer, INTENT(IN)  :: cloud_nscatmoms

!  Trace gas stuff
    !	Trace gas Layer column amounts and cross-sections
    !	Units should be consistent, e.g:
    !	X-sections  are in cm^2/mol  or [DU]
    !	LAYER gas columns are in mol/cm^2  or inverse [DU]

    real(fpk),    dimension ( 0:maxlayers )           , INTENT(IN) :: level_temps
    real(fpk),    dimension ( maxlayers, maxgases, 2) , INTENT(IN) :: layer_gascolumns
    real(fpk),    dimension ( maxwav, maxgases )      , INTENT(IN) :: gas_xsecs
    real(fpk),    dimension ( maxwav)                 , INTENT(IN) :: O3C1_xsecs
    real(fpk),    dimension ( maxwav)                 , INTENT(IN) :: O3C2_xsecs
    real(fpk),    dimension ( maxwav, 0:maxlayers )   , INTENT(IN) :: H2O_xsecs
    real(fpk),    dimension ( maxwav, 0:maxlayers )   , INTENT(IN) :: O2_xsecs

!  Air density and Rayleigh stuff
    !	Units should be consistent, e.g:
    !	X-sections  are in cm^2/mol  or [DU]
    !	Air columns are in mol/cm^2  or inverse [DU]

    real(fpk),    dimension ( maxlayers ), INTENT(IN) :: layer_aircolumns
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
    


!  Local variables
!  ===============

    logical      :: AERFLAG_N, CLDFLAG_N
    integer      :: N, N1, L, G
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
    REAL(fpk)    :: GAS_ABS, AER_SSALB, CLD_SSALB, TABS(maxlayers,maxgases)
    REAL(fpk)    :: PROBLEM_RAY(6,0:2), DEPOL, SK, SIG1, SIG2, TEMP1, TEMP1SQ, TEMP2, TEMP2SQ
    REAL(fpk), parameter :: TZERO = 273.15_fpk


    REAL(fpk)    :: TOT_OPDEP_L     ( MAXLAYERS)
    REAL(fpk)    :: AER_OPDEP_L     ( MAXLAYERS)
    REAL(fpk)    :: AER_SCDEP_L     ( MAXLAYERS)
    REAL(fpk)    :: RAY_SCDEP_L     ( MAXLAYERS)
    REAL(fpk)    :: GAS_OPDEP_L     ( MAXLAYERS)
    TOT_OPDEP_L = 0.0d0
    AER_OPDEP_L = 0.0d0
    AER_SCDEP_L = 0.0d0
    RAY_SCDEP_L = 0.0d0
    GAS_OPDEP_L = 0.0d0


!  Initialize Regular VLIDORT IOPS

    NGREEK_MOMENTS_INPUT = 0
    DELTAU_VERT_INPUT    = 0.0d0
    OMEGA_TOTAL_INPUT    = 0.0d0
    GREEKMAT_TOTAL_INPUT = 0.0d0

!  Set some variables

    NGREEK_MOMENTS_INPUT = 2
    if ( do_aerosols ) then
       NGREEK_MOMENTS_INPUT = aerosol_nscatmoms(w)   ! 10/25/16 Wavelength dependence
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
   TABS      = 0.0d0

!  Start layer loop

   DO N = 1, NLAYERS

      N1 = N - 1

!  Rayeigh scattering optical depth

      RAY_SCDEP = RAYLEIGH_XSECS(W) * LAYER_AIRCOLUMNS(N)
      RAY_ATMOS = RAY_ATMOS + RAY_SCDEP

!  Trace gases optical depths
!    Use Cross-sections defined at Levels.

      GAS_OPDEP = 0.0D0
      DO G = 1, NGASES
         IF (DO_GASES(G) .and. which_gases(g) .eq. 'O3  ' ) THEN
            if ( n.eq.1) then
               temp1 = Level_temps(n1) - TZERO ; temp1sq = temp1 * temp1
               sig1  = gas_xsecs(w,g) + temp1 * O3C1_xsecs(w) + temp1sq * O3C2_xsecs(w)
               temp2 = Level_temps(n)  - TZERO ; temp2sq = temp2 * temp2
               sig2  = gas_xsecs(w,g)  + temp2 * O3C1_xsecs(w) + temp2sq * O3C2_xsecs(w)
            else
               sig1  = sig2
               temp2 = Level_temps(n)  - TZERO ; temp2sq = temp2 * temp2
               sig2  = gas_xsecs(w,g)  + temp2 * O3C1_xsecs(w) + temp2sq * O3C2_xsecs(w)
            endif
            GAS_ABS      = LAYER_GASCOLUMNS(N,G,1) * sig1 + LAYER_GASCOLUMNS(N,G,2) * sig2
         else if (DO_GASES(G) .and. which_gases(g) .eq. 'O2  ' ) THEN
            sig1  = o2_xsecs(w,n1) ; sig2  = o2_xsecs(w,n) 
            GAS_ABS      = LAYER_GASCOLUMNS(N,G,1) * sig1 + LAYER_GASCOLUMNS(N,G,2) * sig2
         else if (DO_GASES(G) .and. which_gases(g) .eq. 'H2O ' ) THEN
            sig1  = h2o_xsecs(w,n1) ; sig2  = h2o_xsecs(w,n) 
            GAS_ABS      = LAYER_GASCOLUMNS(N,G,1) * sig1 + LAYER_GASCOLUMNS(N,G,2) * sig2
         else
            sig1  = gas_xsecs(w,g)
            GAS_ABS  = ( LAYER_GASCOLUMNS(N,G,1) + LAYER_GASCOLUMNS(N,G,2) ) * sig1
         endif
         TABS(N,G)    = GAS_ABS
         GAS_OPDEP    = GAS_OPDEP + GAS_ABS
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

!  debug Tshift -1

!      if ( n.gt.40.and.n.lt.45) then
!         write(*,*)n, AERFLAG_N,DELTAU_VERT_INPUT(N),OMEGA_TOTAL_INPUT(N)
!         if ( n.eq.44)pause'first 3 perturbed'
!      endif

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

!  debug Tshift 2

!      if ( n.gt.40.and.n.lt.45) then
!         write(*,*)n, GREEKMAT_TOTAL_INPUT(1,n,1) ,GREEKMAT_TOTAL_INPUT(2,n,1) 
!         if ( n.eq.44)pause'first 3 perturbed'
!      endif

!  Add clouds

      IF ( CLDFLAG_N ) THEN
         DO K = 1, ngreekmat_entries
            GK = gmask(k)
            sk = dble(smask(k))
            ck = cmask(k)
            DO L = 0, NGREEK_MOMENTS_INPUT
               GREEKMAT_TOTAL_INPUT(L,n,gk) = GREEKMAT_TOTAL_INPUT(L,n,gk)  &
                                            + CLD_WT * SK * CLOUD_SCATMOMS(ck,L,W)
!            if (n.eq.58)write(125,*)n,l,k,CLOUD_SCATMOMS(ck,L,W)
            ENDDO
         ENDDO
!         pause
      ENDIF

!  Phase function normalization

      GREEKMAT_TOTAL_INPUT(0,n,1) = 1.0d0

! debug Surfp

!         if ( n.gt.70) then
!            write(*,*)n, DELTAU_VERT_INPUT(N), OMEGA_TOTAL_INPUT(N), GAS_OPDEP
!            if ( n.eq.nlayers)pause'first pert'
!         endif



!  by MyungjeChoi
   TOT_OPDEP_L(N)=TOT_OPDEP
   AER_OPDEP_L(N)=AER_OPDEP
   AER_SCDEP_L(N)=AER_SCDEP
   RAY_SCDEP_L(N)=RAY_SCDEP
   GAS_OPDEP_L(N)=GAS_OPDEP

!  End layer loop

   ENDDO


!! output of Total optical depth for transmissivity calculation
!    open(10,file='./GEMSTOOL_UVN_Results/TOTAL_OPDEP.dat', status='replace')
!        write(10,'(f8.5)') sum(DELTAU_VERT_INPUT)
!    close(10)

!! output of Integrated Optical Depth by Myungje Choi 20171103 
!    open(10,file='./GEMSTOOL_UVN_Results/Info_OPDEP.dat', status='replace')
!        write(10,999)w,wav,sum(TOT_OPDEP_L(1:NLAYERS)),sum(AER_OPDEP_L(1:NLAYERS)),&
!	                   sum(AER_SCDEP_L(1:NLAYERS)),&
!	                   sum(RAY_SCDEP_L(1:NLAYERS)),sum(GAS_OPDEP_L(1:NLAYERS))
!        write(*,999)w,wav,sum(TOT_OPDEP_L(1:NLAYERS)),sum(AER_OPDEP_L(1:NLAYERS)),&
!	                   sum(AER_SCDEP_L(1:NLAYERS)),&
!	                   sum(RAY_SCDEP_L(1:NLAYERS)),sum(GAS_OPDEP_L(1:NLAYERS))
!999 format(i5,f10.4,f8.5,f8.5,f8.5,f8.5,f8.5)
!    close(10)



!  Trace gases optical depths

!   if ( do_clouds_V2) then
!       write(33,567)w,wav,GAS_ATMOS,RAY_ATMOS,AER_ATMOS, CLD_ATMOS
!   else
!       write(33,567)w,wav,GAS_ATMOS,RAY_ATMOS,AER_ATMOS, 0.0d0
!   endif
!567 format(i5,f10.4,1p8e17.7)

!   write(34,567)w,wav,sum(TABS(1:NLAYERS,1)),sum(TABS(1:NLAYERS,2)),sum(TABS(1:NLAYERS,3)),&
!                      sum(TABS(1:NLAYERS,4)),sum(TABS(1:NLAYERS,5)),sum(TABS(1:NLAYERS,6)),sum(TABS(1:NLAYERS,7))

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
END subroutine GEMSTOOL_UVN_regular_iops

!  End module

End  module GEMSTOOL_UVN_regular_iops_m


