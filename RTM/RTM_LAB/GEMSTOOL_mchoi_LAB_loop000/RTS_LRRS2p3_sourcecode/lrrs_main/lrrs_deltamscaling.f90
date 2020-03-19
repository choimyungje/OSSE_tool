! ###############################################################
! #                                                             #
! #                    THE LIDORT_RRS MODEL                     #
! #                                                             #
! #      (LInearized Discrete Ordinate Radiative Transfer)      #
! #       --         -        -        -         -              #
! #                 (Rotational Raman Scatter)                  #
! #                  -          -     -                         #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Author :      Robert J. D. Spurr                           #
! #                                                             #
! #  Address :     RT SOLUTIONS Inc.                            #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                Tel: (617) 492 1183                          #
! #                                                             #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  Version      :  2.3                                        #
! #  Release Date :  March 2011                                 #
! #                                                             #
! ###############################################################

!    #########################################################
!    #                                                       #
!    #   This Version of LIDORT_RRS comes with a GNU-style   #
!    #   license. Please read the license carefully.         #
!    #                                                       #
!    #########################################################

! ###############################################################
! #                                                             #
! #      External stand-alone routines do Delta-M scaling:      #
! #                                                             #
! #         LRRS_DELTAM_SCALING_2p2                             #
! #         LRRS_DELTAM_SCALING_PLUS_2p2                        #
! #                                                             #
! ###############################################################

      MODULE lrrs_deltamscaling

      USE LRRS_PARS

      PRIVATE
      PUBLIC :: LRRS_DELTAM_SCALING_2P2,&
                LRRS_DELTAM_SCALING_PLUS_2P2

      CONTAINS

      SUBROUTINE LRRS_DELTAM_SCALING_2P2 &
         ( MAXLAYERS, MAXINPMOMS, MAXMOMENTS, MAXPOINTS, &
           DO_DELTAM_SCALING, NPOINTS, NLAYERS, NSTREAMS, &
           NINPMOMS, DELTAU_UNSCALED, OMEGAMOMS_UNSCALED, &
           NMOMENTS, TRUNC_FACTORS, DELTAU_SCALED, OMEGAMOMS_SCALED )

!  This is a self-contained stand-alone module for delta-M scaling all
!  elastic optical property inputs for the LRRS model, Version 2.2.

!  INPUTS
!  ======
!
!  The dimensioning parameters are
!
!      MAXLAYERS  = maximum number of atmospheric layers
!      MAXPOINTS  = maximum number of wavelength points
!      MAXMOMENTS = maximum number of Legendre moments (scaled)
!                   (Twice the maximum number of Streams)
!      MAXINPMOMS = maximum number of Legendre moments (unscaled)
!
!  Flag
!
!      DO_DELTAM_SCALING

!  The control inputs for number of points, layers and moments

!      NPOINTS        = number of points in window
!      NLAYERS        = actual number of atmospheric layers
!      NSTREAMS       = actual number of half-space Discrete Ordinates
!      NINPMOMS       = Actual number of phase function moments

!  The Input optical properties to be modified are
!    We use the product of the single-scatering albedo and phase functio

!      DELTAU_UNSCALED(n,wo)      = total optical depth for extinction.
!      OMEGAMOMS_UNSCALED(n,l,wo) = SSA x Phase function moments
!
!    n  is the layer index,            n  = 1, .... NLAYERS
!    wo is the wavelength index,       wo = 1, ...  NPOINTS_TOTAL
!    l  is the Legendre moment index,  l  = 0, .... NMOMENTS_TOTAL

!  The Output optical properties  are

!      NMOMENTS                 = Equals 2 * NSTREAMS - 1, scaled
!      TRUNC_FACTORS(n,wo)      = Actual deltaM scaling truncation facto
!      DELTAU_SCALED(n,wo)      = total optical depth for extinction.
!      OMEGAMOMS_SCALED(n,l,wo) = SSA x Phase function moments

      IMPLICIT NONE

!  ARGUMENT DECLARATIONS
!  =====================

!  dimensioning

      INTEGER, INTENT(IN) ::          MAXLAYERS
      INTEGER, INTENT(IN) ::          MAXPOINTS
      INTEGER, INTENT(IN) ::          MAXMOMENTS
      INTEGER, INTENT(IN) ::          MAXINPMOMS

!  Flag

      LOGICAL, INTENT(IN) ::          DO_DELTAM_SCALING

!  Number of discrete ordinates in the halfspace

      INTEGER, INTENT(IN) ::          NSTREAMS

!  number of layers

      INTEGER, INTENT(IN) ::          NLAYERS

!  spectral information

      INTEGER, INTENT(IN) ::          NPOINTS

!  optical properties (input)

      INTEGER, INTENT(IN) ::          NINPMOMS
      REAL(FPK), INTENT(IN) :: DELTAU_UNSCALED ( MAXLAYERS, MAXPOINTS )
      REAL(FPK), INTENT(IN) :: OMEGAMOMS_UNSCALED &
           ( MAXLAYERS, 0:MAXINPMOMS, MAXPOINTS )

!  optical properties (modified and output)

      INTEGER, INTENT(INOUT) ::          NMOMENTS
      REAL(FPK), INTENT(INOUT) :: TRUNC_FACTORS   ( MAXLAYERS, MAXPOINTS )
      REAL(FPK), INTENT(INOUT) :: DELTAU_SCALED   ( MAXLAYERS, MAXPOINTS )
      REAL(FPK), INTENT(INOUT) :: OMEGAMOMS_SCALED &
           ( MAXLAYERS, 0:MAXMOMENTS, MAXPOINTS )


!  LOCAL VARIABLES
!  ===============

!  help variables

      INTEGER ::          NM1, N, L, W
      REAL(FPK) :: DNM1, FDEL, FAC1, FAC2, DNL1, FDNL1, OMEGA

!  No deltam SCaling
!  -----------------

!  Copy inputs to output and return

      IF ( .not. DO_DELTAM_SCALING ) THEN
        NMOMENTS = MIN(NINPMOMS,2*NSTREAMS-1)
        DO W = 1, NPOINTS
          DO N = 1, NLAYERS
            DO L = 0, NMOMENTS
              OMEGAMOMS_SCALED(N,L,W) = OMEGAMOMS_UNSCALED(N,L,W)
            ENDDO
            DELTAU_SCALED(N,W) = DELTAU_UNSCALED(N,W)
            TRUNC_FACTORS(N,W) = 0.0D0
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Wih DEltam Scaling
!  ==================

!  Do not require Number of input moments.

!  Truncation moment

      NMOMENTS = 2*NSTREAMS - 1
      NM1  = NMOMENTS+1
      DNM1 = DBLE(2*NM1+1)

!  Start loops over spectrum and layers
!  ------------------------------------

      DO W = 1, NPOINTS
        DO N = 1, NLAYERS

!  overall truncation factor

          OMEGA = OMEGAMOMS_UNSCALED(N,0,W)
          FDEL  = OMEGAMOMS_UNSCALED(N,NM1,W) / OMEGA / DNM1
          TRUNC_FACTORS(N,W) = FDEL
          FAC2  = FDEL * OMEGA
          FAC1  = 1.0d0 - FAC2

!  scale optical depth grid

          DELTAU_SCALED(N,W) = DELTAU_UNSCALED(N,W)* FAC1

!  Scale phase moment entries, and finish results

          DO L = 0, NMOMENTS + 1
            DNL1  = DBLE(2*L + 1 )
            FDNL1 = FAC2 * DNL1
            OMEGAMOMS_SCALED(N,L,W) = &
                  ( OMEGAMOMS_UNSCALED(N,L,W) - FDNL1 ) / FAC1
          ENDDO

!  end layer loop

        ENDDO

!  End spectral loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LRRS_DELTAM_SCALING_2P2

!

      SUBROUTINE LRRS_DELTAM_SCALING_PLUS_2P2 &
         ( MAXLAYERS, MAXINPMOMS, MAXMOMENTS, MAXPOINTS, MAXVARS, &
           DO_DELTAM_SCALING, NPOINTS, NLAYERS, NSTREAMS, &
           NINPMOMS, NFLAGS, NVARS, &
           DELTAU_UNSCALED, OMEGAMOMS_UNSCALED, &
           L_DELTAU_UNSCALED, L_OMEGAMOMS_UNSCALED, &
           NMOMENTS, TRUNC_FACTORS, L_TRUNC_FACTORS, &
           DELTAU_SCALED, OMEGAMOMS_SCALED, &
           L_DELTAU_SCALED, L_OMEGAMOMS_SCALED )

!  This is a self-contained stand-alone module for delta-M scaling all
!  elastic optical property inputs for the LRRS model, Version 2.2.
!    Includes also the lineared EOPS.

!  INPUTS
!  ======
!
!  The dimensioning parameters are
!
!      MAXLAYERS  = maximum number of atmospheric layers
!      MAXPOINTS  = maximum number of wavelength points
!      MAXMOMENTS = maximum number of Legendre moments (scaled)
!                   (Twice the maximum number of Streams)
!      MAXINPMOMS = maximum number of Legendre moments (unscaled)
!      MAXVARS    = Maximum number of weighting functions
!
!  Flag
!
!      DO_DELTAM_SCALING

!  The control inputs for number of points,  layers and moments

!      NPOINTS        = number of points in window
!      NLAYERS        = actual number of atmospheric layers
!      NSTREAMS       = actual number of half-space Discrete Ordinates
!      NINPMOMS       = Actual number of phase function moments
!
!      NVARS(n)       = Number of weighting functions in layer n
!      NFLAGS(n)      = Flag for presence of weighting functions in laye

!  The Input optical properties to be modified are
!    We use the product of the single-scattering albedo and phase functi

!      DELTAU_UNSCALED(n,wo)      = total optical depth for extinction.
!      OMEGAMOMS_UNSCALED(n,l,wo) = SSA x Phase function moments
!
!      L_DELTAU_UNSCALED(q,n,wo)      = Linearized total optical depth f
!      L_OMEGAMOMS_UNSCALED(q,n,l,wo) = Linearized SSA x Phase function

!    n  is the layer index,            n  = 1, .... NLAYERS
!    wo is the wavelength index,       wo = 1, ...  NPOINTS_TOTAL
!    l  is the Legendre moment index,  l  = 0, .... NMOMENTS_TOTAL
!    q  is the Jacobian index          q  = 1, .....NVARS(N)

!  The Output optical properties  are

!      NMOMENTS                 = Equals 2 * NSTREAMS - 1, scaled
!      TRUNC_FACTORS(n,wo)      = Actual deltaM scaling truncation facto
!      L_TRUNC_FACTORS(q,n,wo)  = Linearized truncation factors

!      DELTAU_SCALED(n,wo)      = total optical depth for extinction.
!      OMEGAMOMS_SCALED(n,l,wo) = SSA x Phase function moments

!      L_DELTAU_SCALED(q,n,wo)      = Linearized total optical depth for
!      L_OMEGAMOMS_SCALED(q,n,l,wo) = Linearized SSA x Phase function mo

      IMPLICIT NONE

!  ARGUMENT DECLARATIONS
!  =====================

!  dimensioning

      INTEGER, INTENT(IN) ::          MAXLAYERS
      INTEGER, INTENT(IN) ::          MAXPOINTS
      INTEGER, INTENT(IN) ::          MAXMOMENTS
      INTEGER, INTENT(IN) ::          MAXINPMOMS
      INTEGER, INTENT(IN) ::          MAXVARS

!  Flag

      LOGICAL, INTENT(IN) ::          DO_DELTAM_SCALING

!  Number of discrete ordinates in the halfspace

      INTEGER, INTENT(IN) ::          NSTREAMS

!  number of layers

      INTEGER, INTENT(IN) ::          NLAYERS

!  spectral information

      INTEGER, INTENT(IN) ::          NPOINTS

!  Flagging of linearization

      LOGICAL, INTENT(IN) ::          NFLAGS(MAXLAYERS)

!  Number of linearization parameters

      INTEGER, INTENT(IN) ::          NVARS(MAXLAYERS)

!  optical properties (input)

      INTEGER, INTENT(IN) ::          NINPMOMS
      REAL(FPK), INTENT(IN) :: DELTAU_UNSCALED ( MAXLAYERS, MAXPOINTS )
      REAL(FPK), INTENT(IN) :: OMEGAMOMS_UNSCALED &
           ( MAXLAYERS, 0:MAXINPMOMS, MAXPOINTS )


!  LInearized opical properties (input)

      REAL(FPK), INTENT(IN) :: L_DELTAU_UNSCALED &
           ( MAXVARS, MAXLAYERS, MAXPOINTS )
      REAL(FPK), INTENT(IN) :: L_OMEGAMOMS_UNSCALED &
           ( MAXVARS, MAXLAYERS, 0:MAXINPMOMS, MAXPOINTS )

!  optical properties (modified and output)

      INTEGER, INTENT(INOUT) ::          NMOMENTS
      REAL(FPK), INTENT(INOUT) :: TRUNC_FACTORS   ( MAXLAYERS, MAXPOINTS )
      REAL(FPK), INTENT(INOUT) :: DELTAU_SCALED   ( MAXLAYERS, MAXPOINTS )
      REAL(FPK), INTENT(INOUT) :: OMEGAMOMS_SCALED &
           ( MAXLAYERS, 0:MAXMOMENTS, MAXPOINTS )

!  Linearized optical properties (output)

      REAL(FPK), INTENT(OUT) :: L_TRUNC_FACTORS &
           ( MAXVARS, MAXLAYERS, MAXPOINTS )
      REAL(FPK), INTENT(OUT) :: L_DELTAU_SCALED &
           ( MAXVARS, MAXLAYERS, MAXPOINTS )
      REAL(FPK), INTENT(OUT) :: L_OMEGAMOMS_SCALED &
           ( MAXVARS, MAXLAYERS, 0:MAXMOMENTS, MAXPOINTS )

!  LOCAL VARIABLES
!  ===============

!  help variables

      INTEGER ::          NM1, N, L, W, Q
      REAL(FPK) :: DNM1, FDEL, OMEGA, TERM1, L_OMEG
      REAL(FPK) :: FAC1, FAC2, DNL1, FDNL1, L_FAC2, L_FAC1

!  No deltam SCaling
!  -----------------

!  Copy inputs to output and return

      IF ( .not. DO_DELTAM_SCALING ) THEN
        NMOMENTS = MIN(NINPMOMS,2*NSTREAMS-1)
        DO W = 1, NPOINTS
         DO N = 1, NLAYERS
          DO L = 0, NMOMENTS
           OMEGAMOMS_SCALED(N,L,W) = OMEGAMOMS_UNSCALED(N,L,W)
          ENDDO
          DELTAU_SCALED(N,W) = DELTAU_UNSCALED(N,W)
          TRUNC_FACTORS(N,W) = 0.0D0
          DO Q = 1, NVARS(N)
           DO L = 0, NMOMENTS
            L_OMEGAMOMS_SCALED(Q,N,L,W) = L_OMEGAMOMS_UNSCALED(Q,N,L,W)
           ENDDO
           L_DELTAU_SCALED(Q,N,W) = L_DELTAU_UNSCALED(Q,N,W)
           L_TRUNC_FACTORS(Q,N,W) = 0.0D0
          ENDDO
         ENDDO
        ENDDO
        RETURN
      ENDIF

!  Wih DEltam Scaling
!  ==================

!  Truncation moment

      NMOMENTS = 2*NSTREAMS - 1
      NM1  = NMOMENTS+1
      DNM1 = DBLE(2*NM1+1)

!  Start loops over spectrum and layers
!  ------------------------------------

      DO W = 1, NPOINTS
        DO N = 1, NLAYERS

!  overall truncation factor

          OMEGA = OMEGAMOMS_UNSCALED(N,0,W)
          FDEL  = OMEGAMOMS_UNSCALED(N,NM1,W) / OMEGA / DNM1
          TRUNC_FACTORS(N,W) = FDEL
          FAC2  = FDEL * OMEGA
          FAC1  = 1.0d0 - FAC2

!  scale optical depth grid

          DELTAU_SCALED(N,W)   = DELTAU_UNSCALED(N,W)* FAC1

!  Scale phase moment entries, and finish results

          DO L = 0, NMOMENTS + 1
            DNL1  = DBLE(2*L + 1 )
            FDNL1 = FAC2 * DNL1
            OMEGAMOMS_SCALED(N,L,W) = &
                  ( OMEGAMOMS_UNSCALED(N,L,W) - FDNL1 ) / FAC1
          ENDDO

!  Linearized truncation factor

          IF ( NFLAGS(N) ) THEN
           DO Q = 1, NVARS(N)
            L_OMEG = L_OMEGAMOMS_UNSCALED(Q,N,0,W)
            L_FAC2 = L_OMEGAMOMS_UNSCALED(Q,N,NM1,W) / DNM1
            L_FAC1 = -L_FAC2
            L_DELTAU_SCALED(Q,N,W) = L_DELTAU_UNSCALED(Q,N,W) * FAC1 &
                                     + DELTAU_UNSCALED(N,W)   * L_FAC1
            L_TRUNC_FACTORS(Q,N,W) = 0.0d0
            IF ( FDEL .NE.0.0d0 ) THEN
             L_TRUNC_FACTORS(Q,N,W)=FDEL*((L_FAC2/FAC2)-(L_OMEG/OMEGA))
            ENDIF
            DO L = 0, NMOMENTS + 1
              DNL1  = DBLE(2*L + 1 )
              TERM1 = L_FAC2 * ( OMEGAMOMS_SCALED(N,L,W) - DNL1)
              TERM1 = TERM1 + L_OMEGAMOMS_UNSCALED(Q,N,L,W)
              L_OMEGAMOMS_SCALED(Q,N,L,W) = TERM1 / FAC1
           ENDDO
           ENDDO
          ENDIF

!  end layer loop

        ENDDO

!  End spectral loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE LRRS_DELTAM_SCALING_PLUS_2P2


      END MODULE lrrs_deltamscaling

