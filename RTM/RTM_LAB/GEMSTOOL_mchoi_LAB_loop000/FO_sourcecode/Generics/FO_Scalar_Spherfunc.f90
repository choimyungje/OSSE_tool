module FO_Scalar_spherfuncs_m

public

contains

SUBROUTINE FO_Scalar_spherfuncs ( STARTER, MAXMOMS, NMOMS, DF1, DF2, MU, SS_PLEG )

   implicit none

!  parameter argument

   integer, parameter :: fpk = selected_real_kind(15)

!  I/O

      LOGICAL  , intent(inout) :: STARTER
      INTEGER  , intent(in)    :: MAXMOMS, NMOMS
      REAL(fpk), intent(in)    :: MU
      REAL(fpk), intent(out)   :: SS_PLEG(0:MAXMOMS)
      REAL(fpk), intent(inout) :: DF1(MAXMOMS)
      REAL(fpk), intent(inout) :: DF2(MAXMOMS)

!  Local

      integer  :: L
     
!  Help arrays

      IF ( STARTER ) THEN
         DF1(1) = 0.0d0 ; DF2(1) = 0.0d0
         DO L = 2, NMOMS
            DF1(L) = DBLE(2*L-1) / DBLE(L)
            DF2(L) = DBLE(L-1)   / DBLE(L)
         ENDDO
         STARTER = .false.
      ENDIF

!  Legendre

      SS_PLEG(0) = 1.0d0
      SS_PLEG(1) = MU
      DO L = 2, NMOMS
         SS_PLEG(L) = DF1(L) * SS_PLEG(L-1) * MU - DF2(L) * SS_PLEG(L-2)
      ENDDO

!  Finish

      RETURN
END SUBROUTINE FO_Scalar_spherfuncs

!  Finish

end module FO_Scalar_spherfuncs_m

