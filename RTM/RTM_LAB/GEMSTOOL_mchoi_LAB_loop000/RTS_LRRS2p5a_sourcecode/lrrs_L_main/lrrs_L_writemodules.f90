
! ###############################################################
! #                                                             #
! #                    THE LIDORT_RRS MODEL                     #
! #                                                             #
! #      (LInearized Discrete Ordinate Radiative Transfer)      #
! #       --         -        -        -         -              #
! #                 (Rotational Raman Scattering)               #
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
! #  Email   :     rtsolutions@verizon.net                      #
! #  Website :     www.rtslidort.com                            #
! #                                                             #
! #  Version  #   :  2.5                                        #
! #  Release Date :  March 2017                                 #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  --- History of the model ------------                      #
! #                                                             #
! #  Version 1.0 : 2005, Fortran 77                             #
! #  Version 1.1 : 2007, F77                                    #
! #                                                             #
! #  Version 2.1 : 2009, F77                                    #
! #       * Linearization for Atmospheric-property Jacobians    #
! #       * Single scatter corrections added                    #
! #                                                             #
! #  Version 2.3 : March 2011, Fortran 90                       #
! #       * Simplified Raman-setup procedure                    #
! #       * F90 Version with Type-structure I/O                 #
! #       * Test package developed for installation             #
! #                                                             #
! #  Version 2.5 : March 2017, F90                              #
! #       * Formal BRDF/SLEAVE supplements developed            #
! #       * New test-bed software for testing supplements       #
! #       * Thread-safe Code for OpenMP applications            #
! #       * Complete revision of Taylor-series modules          #
! #       * New User Guide and Review paper                     #
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
! #   SUBROUTINES :                                             #
! #                                                             #
! #             LRRS_WRITE_LIN_INPUT                            #
! #             LRRS_WRITE_LIN_SUP_BRDF_INPUT                   #
! #             LRRS_WRITE_LIN_SUP_SLEAVE_INPUT                 #
! #                                                             #
! ###############################################################

!  This is LRRS Version 2.5.
!      - Newly written for this version, by M. Christi

      MODULE lrrs_l_writemodules_m

!      USE LRRS_PARS_m, Only : LDU

      PRIVATE
      PUBLIC :: LRRS_WRITE_LIN_INPUT, &
                LRRS_WRITE_LIN_SUP_BRDF_INPUT, &
                LRRS_WRITE_LIN_SUP_SLEAVE_INPUT

      CONTAINS

      SUBROUTINE LRRS_WRITE_LIN_INPUT ( &
        NLAYERS,NPOINTS,N_USER_STREAMS,N_USER_RELAZMS,&
        DO_COLUMN_WFS,DO_PROFILE_WFS,DO_SURFACE_WFS,DO_SLEAVE_WFS,&
        DO_NORMALIZED_WFS,DO_AIRPROFILE_WFS,DO_TEMPPROFILE_WFS,DO_TEMPSHIFT_WF,&
        LAYER_VARY_FLAG,LAYER_VARY_NUMBER,&
        N_TOTALCOLUMN_WFS,N_SURFACE_WFS,N_SLEAVE_WFS,&
        LAYER_AIRCOLUMNS_DT,TEMPERATURES_UNSHIFTED,&
        L_DELTAU_INPUT_UNSCALED,L_OMEGAMOMS_ELASTIC_UNSCALED, &
        L_OMEGAPHASFUNC_ELASTIC_UP, L_OMEGAPHASFUNC_ELASTIC_DN )

!   -- Rob mod 5/12/17 for 2p5a, removed NMOMENTS_INPUT, Added L_OMEGAPHASFUNC inputs

      USE LRRS_PARS_m

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   NLAYERS
!      INTEGER, INTENT(IN) ::   NMOMENTS_INPUT
      INTEGER, INTENT(IN) ::   NPOINTS

      INTEGER  , INTENT(IN) ::  N_USER_STREAMS
      INTEGER  , INTENT(IN) ::  N_USER_RELAZMS

!  -----------------
!  Linearized Inputs
!  -----------------

      LOGICAL, INTENT(IN) ::     DO_COLUMN_WFS
      LOGICAL, INTENT(IN) ::     DO_PROFILE_WFS
      LOGICAL, INTENT(IN) ::     DO_SURFACE_WFS
      LOGICAL, INTENT(IN) ::     DO_SLEAVE_WFS

      LOGICAL, INTENT(IN) ::     DO_NORMALIZED_WFS
      LOGICAL, INTENT(IN) ::     DO_AIRPROFILE_WFS
      LOGICAL, INTENT(IN) ::     DO_TEMPPROFILE_WFS
      LOGICAL, INTENT(IN) ::     DO_TEMPSHIFT_WF

      INTEGER, INTENT(IN) ::     N_TOTALCOLUMN_WFS
      INTEGER, INTENT(IN) ::     N_SURFACE_WFS
      INTEGER, INTENT(IN) ::     N_SLEAVE_WFS

      LOGICAL, INTENT(IN) ::     LAYER_VARY_FLAG ( MAX_LAYERS )
      INTEGER, INTENT(IN) ::     LAYER_VARY_NUMBER ( MAX_LAYERS )

      REAL(fpk), INTENT(IN) ::   LAYER_AIRCOLUMNS_DT          ( MAX_LAYERS )
      REAL(fpk), INTENT(IN) ::   TEMPERATURES_UNSHIFTED       ( MAX_LAYERS )

!  -- Rob mod 5/12/17 for 2p5a, changed moment dimension in L_OMEGAMOMS_ELASTIC_UNSCALED

      REAL(fpk), INTENT(IN) ::   L_DELTAU_INPUT_UNSCALED      ( MAX_ATMOSWFS, MAX_LAYERS, MAX_POINTS )
      REAL(fpk), INTENT(IN) ::   L_OMEGAMOMS_ELASTIC_UNSCALED ( MAX_ATMOSWFS, MAX_LAYERS, 0:MAX_MOMENTS, MAX_POINTS )

!  -- Rob mod 5/12/17 for 2p5a, added these inputs for Write-up

      REAL(FPK), INTENT(IN) ::   L_OMEGAPHASFUNC_ELASTIC_UP ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )
      REAL(FPK), INTENT(IN) ::   L_OMEGAPHASFUNC_ELASTIC_DN ( MAX_ATMOSWFS, MAX_LAYERS, MAX_GEOMETRIES, MAX_POINTS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: LAY,MOM,PT,WF,V,NG
      INTEGER :: N_WFS
!      CHARACTER (LEN=9) :: DWFC1S

!  Open output file

      OUTUNIT = 111
      OPEN (OUTUNIT,file = 'LRRS_WRITE_LIN_INPUT.dbg',status = 'replace')

!  Define local variables

      N_WFS = 0
      IF (DO_COLUMN_WFS) THEN
        N_WFS = N_TOTALCOLUMN_WFS
      ELSE IF (DO_PROFILE_WFS) THEN
        N_WFS = MAX_ATMOSWFS
      END IF

!  Write all input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '-----------------'
      WRITE(OUTUNIT,'(A)') 'Linearized Inputs'
      WRITE(OUTUNIT,'(A)') '-----------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_COLUMN_WFS  = ',DO_COLUMN_WFS
      WRITE(OUTUNIT,DWFL)  'DO_PROFILE_WFS = ',DO_PROFILE_WFS
      WRITE(OUTUNIT,DWFL)  'DO_SURFACE_WFS = ',DO_SURFACE_WFS
      WRITE(OUTUNIT,DWFL)  'DO_SLEAVE_WFS  = ',DO_SLEAVE_WFS

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_NORMALIZED_WFS  = ',DO_NORMALIZED_WFS
      WRITE(OUTUNIT,DWFL)  'DO_AIRPROFILE_WFS  = ',DO_AIRPROFILE_WFS
      WRITE(OUTUNIT,DWFL)  'DO_TEMPPROFILE_WFS = ',DO_TEMPPROFILE_WFS
      WRITE(OUTUNIT,DWFL)  'DO_TEMPSHIFT_WF    = ',DO_TEMPSHIFT_WF

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'N_TOTALCOLUMN_WFS  = ',N_TOTALCOLUMN_WFS
      WRITE(OUTUNIT,DWFI)  'N_SURFACE_WFS      = ',N_SURFACE_WFS
      WRITE(OUTUNIT,DWFI)  'N_SLEAVE_WFS       = ',N_SLEAVE_WFS

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFL1)  'LAY = ',LAY,&
          ' LAYER_VARY_FLAG(LAY)   = ',LAYER_VARY_FLAG(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFI1)  'LAY = ',LAY,&
          ' LAYER_VARY_NUMBER(LAY) = ',LAYER_VARY_NUMBER(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' LAYER_AIRCOLUMNS_DT(LAY)    = ',LAYER_AIRCOLUMNS_DT(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' TEMPERATURES_UNSHIFTED(LAY) = ',TEMPERATURES_UNSHIFTED(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO LAY=1,NLAYERS
          DO WF=1,N_WFS
            WRITE(OUTUNIT,DWFR3) &
              'PT = ',PT,' LAY = ',LAY,' WF = ',WF,&
              ' L_DELTAU_INPUT_UNSCALED(WF,LAY,PT) = ',&
                L_DELTAU_INPUT_UNSCALED(WF,LAY,PT)
          END DO
        END DO
      END DO

      !THE 2ND & 3RD DIMENSIONS OF L_OMEGAMOMS_ELASTIC_UNSCALED ARE
      !"IN THE WRONG ORDER"??
!   -- Rob mod 5/12/17 for 2p5a, removed NMOMENTS_INPUT, output all moments to MAX_MOMENTS instead

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO LAY=1,NLAYERS
          DO MOM=0,MAX_MOMENTS      !  Rob mod 5/12/17 for 2p5a
            DO WF=1,N_WFS
              WRITE(OUTUNIT,DWFR4) &
                'PT = ',PT,' LAY = ',LAY,' MOM = ',MOM,' WF = ',WF,&
                ' L_OMEGAMOMS_ELASTIC_UNSCALED(WF,LAY,MOM,PT) = ',&
                  L_OMEGAMOMS_ELASTIC_UNSCALED(WF,LAY,MOM,PT)
            END DO
          END DO
        END DO
      END DO

!  Write the L_OMEGAPHASFUNC inputs
!   -- Rob mod 5/12/17 for 2p5a

      NG = N_USER_STREAMS * N_USER_RELAZMS
      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO LAY=1,NLAYERS
          DO V=1,NG
            DO WF=1,N_WFS
              WRITE(OUTUNIT,DWFR4)  'PT = ',PT,' LAY = ',LAY,' GEOM = ',V,' WF = ',WF,&
              ' L_OMEGAPHASFUNC_ELASTIC_UP(WF,LAY,GEOM,PT) = ',&
                L_OMEGAPHASFUNC_ELASTIC_UP(WF,LAY,V,PT) 
            END DO
          END DO
        END DO
      END DO
      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO LAY=1,NLAYERS
          DO V=1,NG
            DO WF=1,N_WFS
              WRITE(OUTUNIT,DWFR4)  'PT = ',PT,' LAY = ',LAY,' GEOM = ',V,' WF = ',WF,&
              ' L_OMEGAPHASFUNC_ELASTIC_DN(WF,LAY,GEOM,PT) = ',&
                L_OMEGAPHASFUNC_ELASTIC_DN(WF,LAY,V,PT) 
            END DO
          END DO
        END DO
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LRRS_WRITE_LIN_INPUT

!

      SUBROUTINE LRRS_WRITE_LIN_SUP_BRDF_INPUT ( &
        NSTREAMS,N_USER_STREAMS,N_USER_RELAZMS,NPOINTS,N_SURFACE_WFS,&
        LS_EXACTDB_BRDFUNC,LS_BRDF_F_0,LS_BRDF_F,LS_USER_BRDF_F_0,LS_USER_BRDF_F )

      USE LRRS_PARS_m

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   NSTREAMS
      INTEGER, INTENT(IN) ::   N_USER_STREAMS
      INTEGER, INTENT(IN) ::   N_USER_RELAZMS
      INTEGER, INTENT(IN) ::   NPOINTS
      INTEGER, INTENT(IN) ::   N_SURFACE_WFS

      REAL(fpk), INTENT(IN) :: LS_EXACTDB_BRDFUNC &
          ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS )

      REAL(fpk), INTENT(IN) :: LS_BRDF_F_0 &
          ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk), INTENT(IN) :: LS_BRDF_F &
          ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_STREAMS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk), INTENT(IN) :: LS_USER_BRDF_F_0 &
          ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS )
      REAL(fpk), INTENT(IN) :: LS_USER_BRDF_F &
          ( MAX_SURFACEWFS, 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_STREAMS, MAX_POINTS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: MOM,PT,STRM,STRMI,STRMJ,SWF,USTRM,URA
      INTEGER :: NMOMENTS

!  Open output file

      OUTUNIT = 112
      OPEN (OUTUNIT,file = 'LRRS_WRITE_LIN_SUP_BRDF_INPUT.dbg',&
            status = 'replace')

!  Define local variables

      NMOMENTS = 2*NSTREAMS

!  Write all BRDF input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '---------------------------------'
      WRITE(OUTUNIT,'(A)') 'Linearized BRDF Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '---------------------------------'

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO URA=1,N_USER_RELAZMS
          DO USTRM=1,N_USER_STREAMS
            DO SWF=1,N_SURFACE_WFS
              WRITE(OUTUNIT,DWFR4) &
                'PT = ',PT,' URA = ',URA,' USTRM = ',USTRM,&
                ' SWF = ',SWF,&
                ' LS_EXACTDB_BRDFUNC(SWF,USTRM,URA,PT) = ',&
                  LS_EXACTDB_BRDFUNC(SWF,USTRM,URA,PT)
            END DO
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO STRM=1,NSTREAMS
          DO MOM=0,NMOMENTS
            DO SWF=1,N_SURFACE_WFS
              WRITE(OUTUNIT,DWFR4) &
                'PT = ',PT,' STRM = ',STRM,' MOM = ',MOM,&
                ' SWF = ',SWF,&
                ' LS_BRDF_F_0(SWF,MOM,STRM,PT) = ',&
                  LS_BRDF_F_0(SWF,MOM,STRM,PT)
            END DO
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO STRMJ=1,NSTREAMS
          DO STRMI=1,NSTREAMS
            DO MOM=0,NMOMENTS
              DO SWF=1,N_SURFACE_WFS
                WRITE(OUTUNIT,DWFR5) &
                  'PT = ',PT,' STRMJ = ',STRMJ,' STRMI = ',STRMI,&
                  ' MOM = ',MOM,' SWF = ',SWF,&
                  ' LS_BRDF_F(SWF,MOM,STRMI,STRMJ,PT) = ',&
                    LS_BRDF_F(SWF,MOM,STRMI,STRMJ,PT)
              END DO
            END DO
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO USTRM=1,N_USER_STREAMS
          DO MOM=0,NMOMENTS
            DO SWF=1,N_SURFACE_WFS
              WRITE(OUTUNIT,DWFR4) &
                'PT = ',PT,' USTRM = ',USTRM,' MOM = ',MOM,&
                ' SWF = ',SWF,&
                ' LS_USER_BRDF_F_0(SWF,MOM,USTRM,PT) = ',&
                  LS_USER_BRDF_F_0(SWF,MOM,USTRM,PT)
            END DO
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO STRM=1,NSTREAMS
          DO USTRM=1,N_USER_STREAMS
            DO MOM=0,NMOMENTS
              DO SWF=1,N_SURFACE_WFS
                WRITE(OUTUNIT,DWFR5) &
                  'PT = ',PT,' STRM = ',STRM,' USTRM = ',USTRM,&
                  ' MOM = ',MOM,' SWF = ',SWF,&
                  ' LS_USER_BRDF_F(SWF,MOM,USTRM,STRM,PT) = ',&
                    LS_USER_BRDF_F(SWF,MOM,USTRM,STRM,PT)
              END DO
            END DO
          END DO
        END DO
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LRRS_WRITE_LIN_SUP_BRDF_INPUT

!

      SUBROUTINE LRRS_WRITE_LIN_SUP_SLEAVE_INPUT ( &
        NSTREAMS,N_USER_STREAMS,N_USER_RELAZMS,NPOINTS,N_SLEAVE_WFS,&
        LSSL_SLTERM_ISOTROPIC,LSSL_SLTERM_USERANGLES,&
        LSSL_SLTERM_F_0,LSSL_USER_SLTERM_F_0 )

      USE LRRS_PARS_m

      IMPLICIT NONE

      INTEGER, INTENT(IN) ::   NSTREAMS
      INTEGER, INTENT(IN) ::   N_USER_STREAMS
      INTEGER, INTENT(IN) ::   N_USER_RELAZMS
      INTEGER, INTENT(IN) ::   NPOINTS
      INTEGER, INTENT(IN) ::   N_SLEAVE_WFS

      REAL(fpk), INTENT(IN) :: LSSL_SLTERM_ISOTROPIC &
          ( MAX_SLEAVEWFS, MAX_POINTS )
      REAL(fpk), INTENT(IN) :: LSSL_SLTERM_USERANGLES &
          ( MAX_SLEAVEWFS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_POINTS )
      REAL(fpk), INTENT(IN) :: LSSL_SLTERM_F_0 &
          ( MAX_SLEAVEWFS, 0:MAX_MOMENTS, MAX_STREAMS, MAX_POINTS )
      REAL(fpk), INTENT(IN) :: LSSL_USER_SLTERM_F_0 &
          ( MAX_SLEAVEWFS, 0:MAX_MOMENTS, MAX_USER_STREAMS, MAX_POINTS )

!  Local variables

      INTEGER :: OUTUNIT
      INTEGER :: MOM,PT,STRM,SWF,USTRM,URA
      INTEGER :: NMOMENTS

!  Open output file

      OUTUNIT = 104
      OPEN (OUTUNIT,file = 'LRRS_WRITE_LIN_SUP_SLEAVE_INPUT.dbg',&
            status = 'replace')

!  Define local variable

      NMOMENTS = 2*NSTREAMS

!  Write all surface-leaving input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '-----------------------------------'
      WRITE(OUTUNIT,'(A)') 'Linearized SLEAVE Supplement Inputs'
      WRITE(OUTUNIT,'(A)') '-----------------------------------'

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO SWF=1,N_SLEAVE_WFS
          WRITE(OUTUNIT,DWFR2) &
            'PT = ',PT,' SWF = ',SWF,&
            ' LSSL_SLTERM_ISOTROPIC(SWF,PT) = ',&
              LSSL_SLTERM_ISOTROPIC(SWF,PT)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO URA=1,N_USER_RELAZMS
          DO USTRM=1,N_USER_STREAMS
            DO SWF=1,N_SLEAVE_WFS
              WRITE(OUTUNIT,DWFR4) &
                'PT = ',PT,' URA = ',URA,' USTRM = ',USTRM,' SWF = ',SWF,&
                ' LSSL_SLTERM_USERANGLES(SWF,USTRM,URA,PT) = ',&
                  LSSL_SLTERM_USERANGLES(SWF,USTRM,URA,PT)
            END DO
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO STRM=1,NSTREAMS
          DO MOM=0,NMOMENTS
            DO SWF=1,N_SLEAVE_WFS
              WRITE(OUTUNIT,DWFR4) &
                'PT = ',PT,' STRM = ',STRM,' MOM = ',MOM,' SWF = ',SWF,&
                ' LSSL_SLTERM_F_0(SWF,MOM,STRM,PT) = ',&
                  LSSL_SLTERM_F_0(SWF,MOM,STRM,PT)
            END DO
          END DO
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO PT=1,NPOINTS
        DO USTRM=1,N_USER_STREAMS
          DO MOM=0,NMOMENTS
            DO SWF=1,N_SLEAVE_WFS
              WRITE(OUTUNIT,DWFR4) &
                'PT = ',PT,' USTRM = ',USTRM,' MOM = ',MOM,' SWF = ',SWF,&
                ' LSSL_USER_SLTERM_F_0(SWF,MOM,USTRM,PT) = ',&
                  LSSL_USER_SLTERM_F_0(SWF,MOM,USTRM,PT)
            END DO
          END DO
        END DO
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE LRRS_WRITE_LIN_SUP_SLEAVE_INPUT

      END MODULE lrrs_l_writemodules_m
