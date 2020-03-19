Module LRRS_Phasfunc_Sup_Accessory_m

!   -- Rob mod 5/12/17 for 2p5a, New routine for generating PHASFUNC inputs from moments
!    - Given set of moments, range of angles; Legendre calculation is saved for repeated usage

public

contains

subroutine Phasfunc_Sup_Accessory &
        ( LOCAL_MAXMOMS, DO_UPWELLING, DO_DNWELLING,         & ! Input
          LOCAL_NMOMS, N_USER_ANGLES, N_USER_RELAZMS,        & ! Input
          THETA_BOA, USER_ANGLES, USER_RELAZMS, LOCAL_MOMS,  & ! Input
          DO_LEGCALC, DO_SCATANG, COSSCAT_UP, COSSCAT_DN,    & ! Modified input/output
          LEGCALC_UP, LEGCALC_DN, PHASFUNC_UP, PHASFUNC_DN )   ! Modified I/O and True output

   USE LRRS_PARS_m, ONLY : MAX_GEOMETRIES, MAX_USER_STREAMS, MAX_USER_RELAZMS, &
                           fpk, zero, one, deg_to_rad

   implicit none

!  Input Arguments
!  ===============

!  Local dimension

   integer, intent(in) :: LOCAL_MAXMOMS

!  Flags

   logical, intent(in) :: DO_UPWELLING, DO_DNWELLING

!  Local numbers

   integer, intent(in) :: LOCAL_NMOMS, N_USER_ANGLES, N_USER_RELAZMS

!  Geometry input

   real(fpk), intent(in) :: THETA_BOA, USER_ANGLES(MAX_USER_STREAMS), USER_RELAZMS(MAX_USER_RELAZMS)

!  Actual moments

   real(fpk), intent(in) :: LOCAL_MOMS (0:LOCAL_MAXMOMS)

!  output
!  ======

!  Flag for calculating Scattering angles

   logical  , intent(inout) :: DO_SCATANG

!  Legendre calculations

   real(fpk), intent(inout) :: COSSCAT_UP (MAX_GEOMETRIES)
   real(fpk), intent(inout) :: COSSCAT_DN (MAX_GEOMETRIES)

!  Flag for calculating Legendre functions

   logical  , intent(inout) :: DO_LEGCALC

!  Legendre calculations

   real(fpk), intent(inout) :: LEGCALC_UP (MAX_GEOMETRIES,0:LOCAL_MAXMOMS)
   real(fpk), intent(inout) :: LEGCALC_DN (MAX_GEOMETRIES,0:LOCAL_MAXMOMS)

!  True output phase functions

   real(fpk), intent(out) :: PHASFUNC_UP(MAX_GEOMETRIES)
   real(fpk), intent(out) :: PHASFUNC_DN(MAX_GEOMETRIES)

!  Local
!  =====

   integer   :: UA, UM, V, L, n_geometries
   REAL(FPK) :: ALPHA_BOA, PHI_BOA, COSSCAT
   REAL(FPK) :: SALPHA_BOA, STHETA_BOA, CALPHA_BOA, CTHETA_BOA, CPHI_BOA
   REAL(FPK) :: DF1(2:LOCAL_MAXMOMS), DF2(2:LOCAL_MAXMOMS), HELP, PHAS_E

!  Initialize

   PHASFUNC_UP = zero ; PHASFUNC_DN = zero
   n_geometries = N_USER_ANGLES * N_USER_RELAZMS

!  Calculate Scattering Angles
!  ===========================

   if ( .not. DO_SCATANG ) then

!  initialize

      COSSCAT_UP = zero ; COSSCAT_DN = zero

!  scattering cosines

      DO UM = 1, N_USER_ANGLES
         ALPHA_BOA = USER_ANGLES(UM)
         DO UA = 1, N_USER_RELAZMS
            PHI_BOA   = USER_RELAZMS(UA)
            V = N_USER_RELAZMS * (UM-1) + UA
            stheta_boa = sin(theta_boa * deg_to_rad) ; ctheta_boa = sqrt(one-stheta_boa*stheta_boa)
            salpha_boa = sin(alpha_boa * deg_to_rad) ; calpha_boa = sqrt(one-salpha_boa*salpha_boa)
            cphi_boa   = cos(phi_boa * deg_to_rad)
            cosscat_up (v) = - calpha_boa * ctheta_boa + salpha_boa * stheta_boa * cphi_boa
            cosscat_dn (v) = + calpha_boa * ctheta_boa + salpha_boa * stheta_boa * cphi_boa
         ENDDO
      ENDDO

!  re-set flag    

      DO_SCATANG = .true.

   endif

!  Calculate Legendre functions
!  ============================

   if ( .not. DO_LEGCALC ) then

!  initialize

      LEGCALC_UP = zero ; LEGCALC_DN = zero

!  scattering cosines

      DO UM = 1, N_USER_ANGLES
         ALPHA_BOA = USER_ANGLES(UM)
         DO UA = 1, N_USER_RELAZMS
            PHI_BOA   = USER_RELAZMS(UA)
            V = N_USER_RELAZMS * (UM-1) + UA
            stheta_boa = sin(theta_boa * deg_to_rad) ; ctheta_boa = sqrt(one-stheta_boa*stheta_boa)
            salpha_boa = sin(alpha_boa * deg_to_rad) ; calpha_boa = sqrt(one-salpha_boa*salpha_boa)
            cphi_boa   = cos(phi_boa * deg_to_rad)
            cosscat_up (v) = - calpha_boa * ctheta_boa + salpha_boa * stheta_boa * cphi_boa
            cosscat_dn (v) = + calpha_boa * ctheta_boa + salpha_boa * stheta_boa * cphi_boa
         ENDDO
      ENDDO

!  Help arrays

      DO L = 2, LOCAL_NMOMS
        HELP = DBLE(L) ; DF1(L) = DBLE(2*L-1)/HELP ; DF2(L) = DBLE(L-1)/HELP
      ENDDO

!  Upwelling

      IF ( DO_UPWELLING ) THEN
         DO V = 1, N_GEOMETRIES
            COSSCAT = cosscat_up (v) 
            LEGCALC_UP(V,0) = ONE
            LEGCALC_UP(V,1) = COSSCAT
            DO L = 2, LOCAL_NMOMS
               LEGCALC_UP(V,L) = DF1(L) * LEGCALC_UP(V,L-1) * COSSCAT - DF2(L) * LEGCALC_UP(V,L-2)
            ENDDO
         ENDDO
      ENDIF 

!  Downwelling

      IF ( DO_DNWELLING ) THEN
         DO V = 1, N_GEOMETRIES
            COSSCAT = cosscat_dn (v) 
            LEGCALC_DN(V,0) = ONE
            LEGCALC_DN(V,1) = COSSCAT
            DO L = 2, LOCAL_NMOMS
               LEGCALC_DN(V,L) = DF1(L) * LEGCALC_DN(V,L-1) * COSSCAT - DF2(L) * LEGCALC_DN(V,L-2)
            ENDDO
         ENDDO
      ENDIF 

!  re-set flag    

      DO_LEGCALC = .true.

   endif

!  Expansions
!  ==========

   IF ( DO_UPWELLING ) THEN
      DO V = 1, N_GEOMETRIES
         PHAS_E = DOT_PRODUCT(LOCAL_MOMS(0:LOCAL_NMOMS),LEGCALC_UP(V,0:LOCAL_NMOMS))
         PHASFUNC_UP(V) = PHAS_E
      ENDDO
   ENDIF

   IF ( DO_DNWELLING ) THEN
      DO V = 1, N_GEOMETRIES
         PHAS_E = DOT_PRODUCT(LOCAL_MOMS(0:LOCAL_NMOMS),LEGCALC_DN(V,0:LOCAL_NMOMS))
         PHASFUNC_DN(V) = PHAS_E
      ENDDO
   ENDIF

!  Finish

   RETURN
end subroutine Phasfunc_Sup_Accessory

!

subroutine Rayfunc_Sup_Accessory &
        ( DO_UPWELLING, DO_DNWELLING,                   & ! Input
          N_USER_ANGLES, N_USER_RELAZMS, NPOINTS,       & ! Input
          THETA_BOA, USER_ANGLES, USER_RELAZMS, RAYMOM, & ! Input
          DO_SCATANG, COSSCAT_UP, COSSCAT_DN,    & ! Modified input/output
          RAYFUNC_UP, RAYFUNC_DN )                 ! Modified I/O and True output

   USE LRRS_PARS_m, ONLY : MAX_GEOMETRIES, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS, &
                           fpk, zero, one, deg_to_rad

   implicit none

!  Input Arguments
!  ===============

!  Flags

   logical, intent(in)   :: DO_UPWELLING, DO_DNWELLING

!  Local numbers

   integer, intent(in)   :: N_USER_ANGLES, N_USER_RELAZMS, NPOINTS

!  Geometry input

   real(fpk), intent(in) :: THETA_BOA, USER_ANGLES(MAX_USER_STREAMS), USER_RELAZMS(MAX_USER_RELAZMS)

!  rayleigh depolarization

   real(fpk), intent(in) :: RAYMOM(MAX_POINTS)

!  output
!  ======

!  Flag for calculating Scattering angles

   logical  , intent(inout) :: DO_SCATANG

!  Legendre calculations

   real(fpk), intent(inout) :: COSSCAT_UP (MAX_GEOMETRIES)
   real(fpk), intent(inout) :: COSSCAT_DN (MAX_GEOMETRIES)

!  True output phase functions

   real(fpk), intent(out) :: RAYFUNC_UP (MAX_GEOMETRIES,MAX_POINTS)
   real(fpk), intent(out) :: RAYFUNC_DN (MAX_GEOMETRIES,MAX_POINTS)

!  Local
!  =====

   integer   :: UA, UM, V, W, n_geometries
   REAL(FPK) :: ALPHA_BOA, PHI_BOA, COSSCAT, LEGCALC_2
   REAL(FPK) :: SALPHA_BOA, STHETA_BOA, CALPHA_BOA, CTHETA_BOA, CPHI_BOA

!  Initialize

   RAYFUNC_UP = zero ; RAYFUNC_DN = zero
   n_geometries = N_USER_ANGLES * N_USER_RELAZMS

!  Calculate Scattering Angles
!  ===========================

   if ( .not. DO_SCATANG ) then

!  initialize

      COSSCAT_UP = zero ; COSSCAT_DN = zero

!  scattering cosines

      DO UM = 1, N_USER_ANGLES
         ALPHA_BOA = USER_ANGLES(UM)
         DO UA = 1, N_USER_RELAZMS
            PHI_BOA   = USER_RELAZMS(UA)
            V = N_USER_RELAZMS * (UM-1) + UA
            stheta_boa = sin(theta_boa * deg_to_rad) ; ctheta_boa = sqrt(one-stheta_boa*stheta_boa)
            salpha_boa = sin(alpha_boa * deg_to_rad) ; calpha_boa = sqrt(one-salpha_boa*salpha_boa)
            cphi_boa   = cos(phi_boa * deg_to_rad)
            cosscat_up (v) = - calpha_boa * ctheta_boa + salpha_boa * stheta_boa * cphi_boa
            cosscat_dn (v) = + calpha_boa * ctheta_boa + salpha_boa * stheta_boa * cphi_boa
         ENDDO
      ENDDO

!  re-set flag    

      DO_SCATANG = .true.

   endif

!  Expansions
!  ==========

   IF ( DO_UPWELLING ) THEN
      DO V = 1, N_GEOMETRIES
         COSSCAT = cosscat_up (v) 
         LEGCALC_2 = 1.5_fpk * COSSCAT * COSSCAT - 0.5_fpk
         DO W = 1, NPOINTS
            RAYFUNC_UP(V,W) = ONE + LEGCALC_2 * RAYMOM(W)
         ENDDO
      ENDDO      
   ENDIF

   IF ( DO_DNWELLING ) THEN
      DO V = 1, N_GEOMETRIES
         COSSCAT = cosscat_dn (v) 
         LEGCALC_2 = 1.5_fpk * COSSCAT * COSSCAT - 0.5_fpk
         DO W = 1, NPOINTS
            RAYFUNC_DN(V,W) = ONE + LEGCALC_2 * RAYMOM(W)
         ENDDO
      ENDDO
   ENDIF

!  Finish

   RETURN
end subroutine Rayfunc_Sup_Accessory

!  End module

end module LRRS_Phasfunc_Sup_Accessory_m
