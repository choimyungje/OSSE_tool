module GEMSTOOL_RTCALC_PLUS_m

!  Rob Fix, 10/18/16. Introduced SIF
!    GEMSTOOL_RTCALC_PLUS_Actual - SIF inputs added, VLIDORT variables set

!  Rob Fix, 10/25/16. Introduced BRDF, No linearization.
!    GEMSTOOL_RTCALC_PLUS_Actual - BRDF inputs added, VLIDORT variables set

!  Rob Fix 10/24/16
!    GEMSTOOL_RTCALC_PLUS_Final  - New SIF weighting function output.
!    GEMSTOOL_RTCALC_PLUS_Final  - Given UVN/NWS separation  --> SETJACOBIANS

!  Rob Fix, 10/25/16. SLEAVE Structure replace SIF (more general, allows for water-leaving)

!  Rob Fix, 11/30/16. Allowing for linear parameterization of SIF.

!  Module files for VLIDORT

   USE VLIDORT_PARS
   USE VLIDORT_IO_DEFS
   USE VLIDORT_MASTERS

   USE VLIDORT_LIN_IO_DEFS
   USE VLIDORT_LPS_MASTERS
   USE VLIDORT_LCS_MASTERS

!  Module files for the BRDF and SLEAVE inputs. Added 10/25/16.
!    * Note that linearization for BRDF inputs not yet available
!    * SIF amplitude linearization is trivial

   USE VBRDF_SUP_MOD
   USE VSLEAVE_SUP_MOD

public  :: GEMSTOOL_RTCALC_PLUS_Actual

contains

!   (  do_firstorder_option, FO_do_regular_ps,                        & ! FO control

subroutine GEMSTOOL_RTCALC_PLUS_Actual &
    (  do_Profile_Jacobians, do_AerOpdepProfile_Jacobians,                  & ! Main control
       do_SIF, Polynomial, dPolynomial_dSL, SIF_ExactOnly,                  & ! SIF control
       VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,                           & ! VLIDORT regular inputs
       VLIDORT_LinFixIn, VLIDORT_LinModIn, VLIDORT_LinSup,                  & ! VLIDORT linearized inputs
       nlayers, NGREEK_MOMENTS_INPUT, height_grid, lambertian_albedo,       & ! Control/Surface Proxies
       DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,          & ! Atmos-Optical Proxies
       L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,    & ! Atmos-Optical Proxies (linearized)
       GEMSTOOL_SLEAVE_Results, GEMSTOOL_BRDF_Results,                      & ! Surface Supplement Results
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

!  SIF control
!    11/30/16. Add SIF Polynomial factor and derivatives
!    11/30/16. Add SIF Exact-only approximation option

   logical  , intent(in) :: do_SIF
   REAL(fpk), intent(in) :: Polynomial, dPolynomial_dSL(2)
   Logical  , intent(in) :: SIF_ExactOnly
!   real(fpk), intent(in) :: SIF755_Amplitude

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

!  Proxy inputs for VLIDORT.
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

!  SLEAVE input structure, added 10/25/16. FIrst Attempt 10/18/16 was only for SIF....
!   REAL(fpk), intent(in) :: GEMSTOOL_SIF          ( MAX_GEOMETRIES )
!   REAL(fpk), intent(in) :: dGEMSTOOLSIF_dSCALING ( MAX_GEOMETRIES )

   TYPE(VSLEAVE_Sup_Outputs), INTENT(IN) :: GEMSTOOL_SLEAVE_Results

!  BRDF input structure, added 10/25/16

   TYPE(VBRDF_Sup_Outputs), INTENT(IN) :: GEMSTOOL_BRDF_Results

!  Outputs
!  -------

!  VLIDORT output structure

   TYPE(VLIDORT_Outputs), intent(inout)              :: VLIDORT_Out

!  VLIDORT linearized output structure

   TYPE(VLIDORT_LinOutputs), intent(inout)           :: VLIDORT_LinOut

!  Local. Additional quantities for SIF, 10/18/16.

   logical, parameter :: skip_vlidort = .false.
   integer   :: jj, g, L, nstreams, nmoments, ngeoms, qs, nsleave_wfs
   real(fpk) :: Iso_Value, LinIso_Value

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
   VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO     = LAMBERTIAN_ALBEDO    ! May be zero if BRDF in effect

!  proxies

   ngeoms   = VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS
   nstreams = VLIDORT_FixIn%Cont%TS_NSTREAMS
   nmoments = 2 * nstreams

!  Surface-leaving inputs
!  ----------------------

!  10/18/16, 10/25/16. Here is where you set the VLIDORT inputs
!     Only operating with unpolarized stuff (JJ = 1, indicates first Stokes component
!         VLIDORT supplement arrays are pre-zeroed.

!  11/30/16. Upgrade to include linear parameterization of Fluorescence
!    -- if flagged apply polynomial to basic 755 result.
!    -- Use of the Exact-only approximation to SI is in force.

   JJ = 1
   if ( VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING ) THEN
      if ( VLIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC ) THEN
         do g = 1, ngeoms
            Iso_Value = Polynomial * GEMSTOOL_SLEAVE_Results%SL_SLTERM_ISOTROPIC(JJ,g)
            VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC(JJ,g)        = Iso_Value
            VLIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES(JJ,g,g,g)   = Iso_Value
            if ( .not. SIF_ExactOnly ) then
               VLIDORT_Sup%SLEAVE%TS_SLTERM_F_0(0,JJ,1:nstreams,g) = Iso_Value
               VLIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0(0,JJ,g,g)     = Iso_Value
            endif
         enddo
      else ! This is rather academic as fluorescence is isotropic......
         do g = 1, ngeoms
            VLIDORT_Sup%SLEAVE%TS_SLTERM_ISOTROPIC(JJ,g)        = &
                          Polynomial * GEMSTOOL_SLEAVE_Results%SL_SLTERM_ISOTROPIC(JJ,g)
            VLIDORT_Sup%SLEAVE%TS_SLTERM_USERANGLES(JJ,g,g,g)   = &
                          Polynomial * GEMSTOOL_SLEAVE_Results%SL_SLTERM_USERANGLES(JJ,g,g,g)
            if ( .not. SIF_ExactOnly ) then
               do L = 0, nmoments
                  VLIDORT_Sup%SLEAVE%TS_SLTERM_F_0(L,JJ,1:nstreams,g) = &
                          Polynomial * GEMSTOOL_SLEAVE_Results%SL_SLTERM_F_0(L,JJ,1:nstreams,g)
                  VLIDORT_Sup%SLEAVE%TS_USER_SLTERM_F_0(L,JJ,g,g)     = &
                          Polynomial * GEMSTOOL_SLEAVE_Results%SL_USER_SLTERM_F_0(0,JJ,g,g) 
               enddo
            endif
         enddo
      endif
   endif

!  BRDF inputs. 10/25/16. Copy the GEMSTOOL_BRDF_Results variables to VLIDORT
!   VBRDF inputs are all pre-initialized; emissivities not requird

   ! mchoi --point to modify
   if ( .not.VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE ) THEN 
      VLIDORT_Sup%BRDF%TS_BRDF_F_0        = GEMSTOOL_BRDF_Results%BS_BRDF_F_0
      VLIDORT_Sup%BRDF%TS_BRDF_F          = GEMSTOOL_BRDF_Results%BS_BRDF_F
      VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0   = GEMSTOOL_BRDF_Results%BS_USER_BRDF_F_0
      VLIDORT_Sup%BRDF%TS_USER_BRDF_F     = GEMSTOOL_BRDF_Results%BS_USER_BRDF_F
      VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC = GEMSTOOL_BRDF_Results%BS_DBOUNCE_BRDFUNC
   endif

!  Linearized inputs

   VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT    = L_DELTAU_VERT_INPUT
   VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT    = L_OMEGA_TOTAL_INPUT
   VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT = L_GREEKMAT_TOTAL_INPUT

!  SIF, 10/18/16. Here is where you set the Linearized VLIDORT inputs
!   Currently defaults to the ISF case, though you don't know it.....
!                 --> Only 1 weighting function (qs = 1)

!  11/30/16. Upgrade to include linear parameterization of Fluorescence
!    -- if flagged apply polynomial to basic 755 result.
!    -- Use of the Exact-only approximation to SI is in force.

!  TEMPORARY FIX, only works for SIF

   nsleave_wfs = VLIDORT_LinFixIn%Cont%TS_N_SLEAVE_WFS
   if ( do_SIF ) then
     if ( VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING.and.VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS )   THEN 
       if ( VLIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC ) THEN
         do qs = 1, nsleave_wfs
           do g = 1, ngeoms
             LinIso_Value = dPolynomial_dSL(qs) * GEMSTOOL_SLEAVE_Results%SL_SLTERM_ISOTROPIC(JJ,g)
             VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_ISOTROPIC(qs,1,g)        = LinIso_Value
             VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_USERANGLES(qs,1,g,g,g)   = LinIso_Value
             if ( .not. SIF_ExactOnly ) then
                VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_F_0(qs,0,1,1:nstreams,g) = LinIso_Value
                VLIDORT_LinSup%SLEAVE%TS_LSSL_USER_SLTERM_F_0(qs,0,1,g,g)     = LinIso_Value
             endif
           enddo
         enddo
       else ! This is rather academic as fluorescence is isotropic......
         do qs = 1, nsleave_wfs
           do g = 1, ngeoms
             VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_ISOTROPIC(qs,JJ,g)        = &
                          dPolynomial_dSL(qs) * GEMSTOOL_SLEAVE_Results%SL_SLTERM_ISOTROPIC(JJ,g)
             VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_USERANGLES(qs,JJ,g,g,g)   = &
                          dPolynomial_dSL(qs) * GEMSTOOL_SLEAVE_Results%SL_SLTERM_USERANGLES(JJ,g,g,g)
             if ( .not. SIF_ExactOnly ) then
               do L = 0, nmoments
                  VLIDORT_LinSup%SLEAVE%TS_LSSL_SLTERM_F_0(qs,L,JJ,1:nstreams,g) = &
                          dPolynomial_dSL(qs) * GEMSTOOL_SLEAVE_Results%SL_SLTERM_F_0(L,JJ,1:nstreams,g)
                  VLIDORT_LinSup%SLEAVE%TS_LSSL_USER_SLTERM_F_0(qs,L,JJ,g,g)     = &
                          dPolynomial_dSL(qs) * GEMSTOOL_SLEAVE_Results%SL_USER_SLTERM_F_0(0,JJ,g,g) 
               enddo
             endif
           enddo
         enddo
       endif
     endif
   endif

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

!  End module

end module GEMSTOOL_RTCALC_PLUS_m
