module RTSMie_IO_readwrite

!  Mie modules

  use RTSMie_parameters_m

!  Auxiliary stand-alone routines for Reading and writing

!    RTSMie_read_configfile.          Reads the Regular configuration file
!    RTSMie_read_configbimodal.       Reads the Bimodal configuration file

!    RTSMie_write_standard.           Standard output
!    RTSMie_write_extended.           Standard + Linearized  output

!    RTSMie_write_bimodal_standard.   Standard output, bimodal
!    RTSMie_write_bimodal_extended.   Standard + Linearized output, bimodal

public

contains

subroutine RTSMie_read_configfile ( filename,               & ! File name
     Do_Expcoeffs, Do_Fmatrix, Do_Monodisperse,             & ! Gp 1 Inputs (Flags)
     Do_LinearRef,  Do_psd_Linearize,                       & ! Gp 1 Inputs (Flags)
     psd_Index, psd_pars, FixR1R2, R1, R2,                  & ! Gp 2 Inputs (PSD)
     R1R2_cutoff, nblocks, nweights,                        & ! Gp 2 Inputs (PSD)
     xparticle_limit, MonoRadius, lambda, n_real, n_imag,   & ! Gp 3 Inputs (Optical)
     n_Fmatrix_angles, Fmatrix_angles )                       ! Gp 4 Inputs (Fmatrix control)

   use RTSMie_parameters_m, ONLY : dp, max_Mie_angles

   implicit none

!  Filename
!  --------

   character*(*), intent(in)  :: filename

!  RTSMIE MASTER : INPUT ARGUMENTS
!  ===============================

!  Boolean Input arguments
!  -----------------------

!  Flags for Expansion Coefficient and Fmatrix calculations

   logical, intent(out)  :: Do_Expcoeffs
   logical, intent(out)  :: Do_Fmatrix

!  Logical flag for Monodisperse calculation

   logical, intent(out)  :: Do_monodisperse

!  linearization of monodisperse quantitites

   logical, intent(out)  :: Do_LinearRef    ! w.r.t. Real/Imag parts ref index

!  Linaerization for PSD parameters

   logical, intent(out)  :: Do_psd_linearize        ! PSD parameters

!  PSD inputs (Group 2)
!  --------------------

!  PSD index and parameters

   integer      , intent(out)   :: psd_Index
   real(kind=dp), intent(out)   :: psd_pars (3)

!  Flag for making an internal Fix of R1 and R2
!    ( Not relevant for the Old-style PSDs)

   logical     , intent(out)    :: FixR1R2

!  R1 and R2 (intent(inout) from the new program)

   real(kind=dp), intent(out)   :: R1, R2

!    R1R2_cutoff particle size for setting R1 and R2 internally

   REAL(KIND=dp), intent(out)  :: R1R2_cutoff

!    PSD integration ranges is divided into so many blocks.
!    For each block, integrate using Gaussian quadrature, with so many weights.

   INTEGER,       intent(out)  :: nblocks
   INTEGER      , intent(out)  :: nweights

!  Optical and Monodisperse (Group 3)
!  ----------------------------------

!  Limiting particle size value. Set to 10000.0 default.
!   If you exceed this, program will tell you to increase dimensioning.

   REAL(KIND=dp), intent(out)  :: xparticle_limit

!  Monodisperse radius

   real(kind=dp), intent(out)   :: Monoradius

!: Wavelength, refractive index

   real(kind=dp), intent(out)   :: lambda, n_real, n_imag

!  F-matrix Angular control input (Group 4)
!  ----------------------------------------

!  Calculate F-matrix at user-defined angles (do_Fmatrix flag MUST BE set)
!       n_Fmatrix_angles = number of user-defined angles.
!       Fmatrix_angles   = user-defined angles, in DEGREES between [0, 180]

   INTEGER      , intent(out) :: n_Fmatrix_angles
   REAL(KIND=dp), intent(out) :: Fmatrix_angles(max_Mie_angles)

!  Local
!  -----

   integer       :: j
   real(kind=dp) :: div

!  Open file and read
!  ------------------

   open(1,file=TRIM(ADJUSTL(filename)),status = 'old')
   read(1,*)
   read(1,*)
   read(1,*) ! First group (Boolean flags)
   read(1,*)
   read(1,*) Do_expcoeffs
   read(1,*) Do_Fmatrix
   read(1,*) Do_monodisperse
   read(1,*) Do_LinearRef
   read(1,*) Do_psd_Linearize
   read(1,*)
   read(1,*) ! Second group (PSD-related variables)
   read(1,*)
   read(1,*) psd_Index
   read(1,*) psd_pars(1)
   read(1,*) psd_pars(2)
   read(1,*) psd_pars(3)
   read(1,*) FixR1R2
   read(1,*) R1
   read(1,*) R2
   read(1,*) R1R2_cutoff
   read(1,*) nblocks
   read(1,*) nweights
   read(1,*)
   read(1,*) ! Third group (Monodisperse and optical inputs)
   read(1,*)
   read(1,*) xparticle_limit
   read(1,*) Monoradius
   read(1,*) lambda
   read(1,*) n_real
   read(1,*) n_imag
   read(1,*)
   read(1,*) ! Fourth group (Fmatrix inputs)
   read(1,*)
   read(1,*) n_Fmatrix_angles
   close(1)

!  Fmatrix angles are regular here between 0 and 180.

   if ( do_Fmatrix ) then
      div = 180.0d0 / dble( n_Fmatrix_angles - 1 )
      do j = 1, n_Fmatrix_angles
         Fmatrix_angles(j) = dble(j-1) * div
      enddo
   else
      Fmatrix_angles = 0.0d0
   endif

!  return

   return
end subroutine RTSMie_read_configfile

subroutine RTSMie_read_configbimodal ( filename,            & ! File name
     Do_Expcoeffs, Do_Fmatrix, Do_Monodisperse,             & ! Gp 1 Inputs (Flags)
     Do_LinearRef,  Do_psd_Linearize,                       & ! Gp 1 Inputs (Flags)
     psd_Index, psd_pars, FixR1R2, R1, R2,                  & ! Gp 2 Inputs (PSD)
     R1R2_cutoff, nblocks, nweights,                        & ! Gp 2 Inputs (PSD)
     xparticle_limit, MonoRadius, lambda, n_real, n_imag,   & ! Gp 3 Inputs (Optical)
     n_Fmatrix_angles, Fmatrix_angles, fraction )             ! Gp 4/5 Inputs (Fmatrix and fraction)

   use RTSMie_parameters_m, ONLY : dp, max_Mie_angles

   implicit none

!  Filename
!  --------

   character*(*), intent(in)  :: filename

!  RTSMIE MASTER : INPUT ARGUMENTS
!  ===============================

!  Boolean Input arguments
!  -----------------------

!  Flags for Expansion Coefficient and Fmatrix calculations

   logical, intent(out)  :: Do_Expcoeffs
   logical, intent(out)  :: Do_Fmatrix

!  Logical flag for Monodisperse calculation

   logical, intent(out)  :: Do_monodisperse

!  linearization of monodisperse quantitites

   logical, intent(out)  :: Do_LinearRef    ! w.r.t. Real/Imag parts ref index

!  Linaerization for PSD parameters

   logical, intent(out)  :: Do_psd_linearize        ! PSD parameters

!  PSD inputs (Group 2)
!  --------------------

!  PSD index and parameters

   integer      , intent(out)   :: psd_Index(2)
   real(kind=dp), intent(out)   :: psd_pars (3,2)

!  Flag for making an internal Fix of R1 and R2
!    ( Not relevant for the Old-style PSDs)

   logical     , intent(out)    :: FixR1R2(2)

!  R1 and R2 (intent(inout) from the new program)

   real(kind=dp), intent(out)   :: R1(2), R2(2)

!    R1R2_cutoff particle size for setting R1 and R2 internally

   REAL(KIND=dp), intent(out)   :: R1R2_cutoff(2)

!    PSD integration ranges is divided into so many blocks.
!    For each block, integrate using Gaussian quadrature, with so many weights.

   INTEGER,       intent(out)   :: nblocks(2)
   INTEGER      , intent(out)   :: nweights(2)

!  Optical and Monodisperse (Group 3)
!  ----------------------------------

!  Limiting particle size value. Set to 10000.0 default.
!   If you exceed this, program will tell you to increase dimensioning.

   REAL(KIND=dp), intent(out)   :: xparticle_limit

!  Monodisperse radius

   real(kind=dp), intent(out)   :: Monoradius

!: Wavelength, refractive index

   real(kind=dp), intent(out)   :: lambda, n_real(2), n_imag(2)

!  F-matrix Angular control input (Group 4)
!  ----------------------------------------

!  Calculate F-matrix at user-defined angles (do_Fmatrix flag MUST BE set)
!       n_Fmatrix_angles = number of user-defined angles.
!       Fmatrix_angles   = user-defined angles, in DEGREES between [0, 180]

   INTEGER      , intent(out) :: n_Fmatrix_angles
   REAL(KIND=dp), intent(out) :: Fmatrix_angles(max_Mie_angles)

!  Bimodal fractional weight for PSD 1 (Group 5)
!  ---------------------------------------------

   REAL(KIND=dp), intent(out)  :: fraction

!  Local
!  -----

   integer       :: j
   real(kind=dp) :: div

!  Open file and read
!  ------------------

   open(1,file=TRIM(ADJUSTL(filename)),status = 'old')
   read(1,*)
   read(1,*)
   read(1,*) ! First group (Boolean flags)
   read(1,*)
   read(1,*) Do_expcoeffs
   read(1,*) Do_Fmatrix
   read(1,*) Do_monodisperse
   read(1,*) Do_LinearRef
   read(1,*) Do_psd_Linearize
   read(1,*)
   read(1,*) ! Second group (PSD-related variables)
   read(1,*)
   read(1,*) psd_Index(1),   psd_Index(2)
   read(1,*) psd_pars(1,1),  psd_pars(1,2)
   read(1,*) psd_pars(2,1),  psd_pars(2,2)
   read(1,*) psd_pars(3,1),  psd_pars(3,2)
   read(1,*) FixR1R2(1),     FixR1R2(2)
   read(1,*) R1(1),          R1(2)
   read(1,*) R2(1),          R2(2)
   read(1,*) R1R2_cutoff(1), R1R2_cutoff(2)
   read(1,*) nblocks (1),    nblocks (2)
   read(1,*) nweights(1),    nweights(2)
   read(1,*)
   read(1,*) ! Third group (Monodisperse and optical inputs)
   read(1,*)
   read(1,*) xparticle_limit
   read(1,*) Monoradius
   read(1,*) lambda
   read(1,*) n_real(1),n_real(2)
   read(1,*) n_imag(1),n_imag(2)
   read(1,*)
   read(1,*) ! Fourth group (Fmatrix inputs)
   read(1,*)
   read(1,*) n_Fmatrix_angles
   read(1,*)
   read(1,*) ! Fifth group (fractional input)
   read(1,*)
   read(1,*) fraction
   close(1)

!  Some settings fixed

   do_monodisperse = .false.
   Monoradius = 0.0d0

!  Fmatrix angles are regular here between 0 and 180.

   if ( do_Fmatrix ) then
      div = 180.0d0 / dble( n_Fmatrix_angles - 1 )
      do j = 1, n_Fmatrix_angles
         Fmatrix_angles(j) = dble(j-1) * div
      enddo
   else
      Fmatrix_angles = 0.0d0
   endif

!  return

   return
end subroutine RTSMie_read_configbimodal


subroutine RTSMie_write_standard                              &
       ( filename, Do_Expcoeffs, Do_Fmatrix, Do_monodisperse, & ! Control Inputs
         n_Fmatrix_angles, Fmatrix_angles,                    & ! Control Inputs
         RTSMie_bulk, RTSMie_asymm, RTSMie_ncoeffs,           & ! RTSMie Results
         RTSMie_expcoeffs, RTSMie_Fmatrix, RTSMie_dist )        ! RTSMie Results

   use RTSMie_parameters_m, ONLY : dp, max_Mie_angles

   implicit none

!  routine inputs
!  ==============

!  Filename

   character*(*), intent(in)  :: filename

!  Flags for Expansion Coefficient and Fmatrix calculations

   logical, intent(in)  :: Do_Expcoeffs
   logical, intent(in)  :: Do_Fmatrix

!  Logical flag for Monodisperse calculation

   logical, intent(in)  :: Do_monodisperse

!  F-matrix Angular control input
!  Calculate F-matrix at user-defined angles (do_Fmatrix flag MUST BE set)
!       n_Fmatrix_angles = number of user-defined angles.
!       Fmatrix_angles   = user-defined angles, in DEGREES between [0, 180]

   INTEGER, intent(in)       :: n_Fmatrix_angles
   REAL(KIND=dp), intent(in) :: Fmatrix_angles(max_Mie_angles)

!  RTSMIE MASTER : OUTPUT ARGUMENTS
!  ================================

!  Bulk distribution parameters
!  ----------------------------

!    1 = Extinction Coefficient
!    2 = Scattering Coefficient
!    3 = Single scattering albedo

   real(kind=dp), intent(in)   :: RTSMie_bulk (3)

!  Expansion coefficients and Asymmetry parameter
!  ----------------------------------------------

   integer      , intent(in)   :: RTSMie_ncoeffs
   real(kind=dp), intent(in)   :: RTSMie_expcoeffs (6,0:max_Mie_angles)
   real(kind=dp), intent(in)   :: RTSMie_asymm

!  F-matrix,  optional output
!  --------------------------

   real(kind=dp), intent(in)   :: RTSMie_Fmatrix(4,max_Mie_angles)

!  Distribution parameters
!  -----------------------

!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(kind=dp), intent(in)   :: RTSMie_dist (5)

!  Local

   integer       :: L, k, cmask(6),fmask(4)

!  Old indexing

!   data cmask / 1, 3, 4, 6, 2, 5 /
!   data fmask / 1, 3, 2, 4 /

!  New indexing

   data cmask / 1, 2, 3, 4, 5, 6 /
   data fmask / 1, 2, 3, 4 /

!  open file

   open(35,file=TRIM(ADJUSTL(filename)),status = 'unknown' )

!  write distribution information

   if ( Do_monodisperse ) then
      write(35,'(/a)')' *** MONODISPERSE OUTPUT ONLY '
   else
      write(35,'(/a/)')' PSD output ------------- '
      write(35,'(A,T25,1pe15.6)')'Number density     = ',RTSMie_dist (1)
      write(35,'(A,T25,1pe15.6)')'Cross-section      = ',RTSMie_dist (2)
      write(35,'(A,T25,1pe15.6)')'Volume             = ',RTSMie_dist (3)
      write(35,'(A,T25,1pe15.6)')'Effective Radius   = ',RTSMie_dist (4)
      write(35,'(A,T25,1pe15.6)')'Effective Variance = ',RTSMie_dist (5)
   endif

!  write Bulk property output
!    ASYMM will be Zero is Do_expcoeffs not set.

   write(35,'(/a/)')' Bulk property output ------------- '
   write(35,'(A,T40,1pe20.11)')'Extinction Coefficient = ',RTSMie_bulk (1)
   write(35,'(A,T40,1pe20.11)')'Scattering Coefficient = ',RTSMie_bulk (2)
   write(35,'(A,T40,1pe20.11)')'Singlescattering albedo= ',RTSMie_bulk (3)
   write(35,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',RTSMie_asymm

!  Expansion coefficients
  
   if ( Do_expcoeffs ) then
      write(35,'(/a)')'Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
      write(35,'(a,T30,I5/)')' - Number of coefficients =  ',RTSMie_ncoeffs
      do L = 0, RTSMie_ncoeffs
        write(35,'(i5,1p6e20.11)')L,(RTSMie_expcoeffs (cmask(k),L),k=1,6)
      enddo
   endif

!  F-matrix

   if ( Do_Fmatrix ) then
      write(35,'(/a)')' F-matrix (F11=F22/F33=F44/F12/F34) @ angles (degs)'
      write(35,'(a,I5/)')' - Number of angles =  ',n_Fmatrix_angles
      do L = 1, n_Fmatrix_angles
        write(35,'(f6.2,4F20.11)')Fmatrix_angles(L),(RTSMie_Fmatrix(fmask(k),L),k=1,4)
      enddo
   endif

!  Close file

   close(35)

!  Finish

   return
end subroutine RTSMie_write_standard

subroutine RTSMie_write_extended                              &
       ( filename, Do_Expcoeffs, Do_Fmatrix, Do_monodisperse, & ! Control Inputs
         Do_LinearRef, Do_psd_Linearize, Donorm,              & ! Control Inputs
         psd_index, n_real, n_imag, psd_pars,                 & ! Control Inputs
         n_Fmatrix_angles, Fmatrix_angles,                    & ! Control Inputs
         RTSMie_bulk, RTSMie_asymm, RTSMie_ncoeffs,           & ! RTSMie Results
         RTSMie_expcoeffs, RTSMie_Fmatrix, RTSMie_dist,       & ! RTSMie Results
         LPSD_RTSMie_bulk, LPSD_RTSMie_asymm,                 & ! RTSMie Linearized Results
         LPSD_RTSMie_expcoeffs, LPSD_RTSMie_Fmatrix,          & ! RTSMie Linearized Results
         LRFE_RTSMie_bulk, LRFE_RTSMie_asymm,                 & ! RTSMie Linearized Results
         LRFE_RTSMie_expcoeffs, LRFE_RTSMie_Fmatrix,          & ! RTSMie Linearized Results
         LPSD_RTSMie_dist )                                     ! RTSMie Linearized Results

   use RTSMie_parameters_m, ONLY : dp, max_Mie_angles

   implicit none

!  routine inputs
!  ==============

!  Filename

   character*(*), intent(in)  :: filename

!  Flags for Expansion Coefficient and Fmatrix calculations

   logical, intent(in)  :: Do_Expcoeffs
   logical, intent(in)  :: Do_Fmatrix

!  Logical flag for Monodisperse calculation

   logical, intent(in)  :: Do_monodisperse

!  linearization w.r.t Ref Index

   logical, intent(in)  :: Do_LinearRef    ! w.r.t. Real/Imag parts ref index

!  Linearization for PSD parameters

   logical, intent(in)  :: Do_psd_linearize        ! PSD parameters

!  Flag for writing Normalized output
!    (Mie output is not normalized per se)

   logical, intent(in)  :: Donorm

!  PSD index

   integer, intent(in)  :: psd_Index

!  Refractive indices and PSD parameters.
!    --- Required for writing normalized output

   REAL(KIND=dp), intent(in) :: n_real, n_imag, psd_pars(3)

!  F-matrix Angular control input
!  Calculate F-matrix at user-defined angles (do_Fmatrix flag MUST BE set)
!       n_Fmatrix_angles = number of user-defined angles.
!       Fmatrix_angles   = user-defined angles, in DEGREES between [0, 180]

   INTEGER, intent(in)       :: n_Fmatrix_angles
   REAL(KIND=dp), intent(in) :: Fmatrix_angles(max_Mie_angles)

!  RTSMIE MASTER : OUTPUT ARGUMENTS
!  ================================

!  Bulk distribution parameters
!  ----------------------------

!    1 = Extinction Coefficient
!    2 = Scattering Coefficient
!    3 = Single scattering albedo

   real(kind=8), intent(in)   :: RTSMie_bulk (3)

!  linearizations w.r.t. PSD parameters

   real(kind=8), intent(in)  :: LPSD_RTSMie_bulk (3,3)

!  linearizations w.r.t. RefIdx parameters

   real(kind=8), intent(in)  :: LRFE_RTSMie_bulk (3,2)

!  Expansion coefficients and Asymmetry parameter
!  ----------------------------------------------

   integer     , intent(in)   :: RTSMie_ncoeffs
   real(kind=8), intent(in)   :: RTSMie_expcoeffs (6,0:max_mie_angles)
   real(kind=8), intent(in)   :: RTSMie_asymm

!  linearizations w.r.t. PSD parameters

   real(kind=8), intent(in)  :: LPSD_RTSMie_expcoeffs (6,0:max_mie_angles,3)
   real(kind=8), intent(in)  :: LPSD_RTSMie_asymm(3)

!  linearizations w.r.t. RefIdx parameters

   real(kind=8), intent(in)  :: LRFE_RTSMie_expcoeffs (6,0:max_mie_angles,2)
   real(kind=8), intent(in)  :: LRFE_RTSMie_asymm(2)

!  F-matrix,  optional output
!  --------------------------

   real(kind=8), intent(in)  :: RTSMie_Fmatrix (4,max_mie_angles)

!  Linearizations of F-matrix

   real(kind=8), intent(in)  :: LPSD_RTSMie_Fmatrix (4,max_mie_angles,3)
   real(kind=8), intent(in)  :: LRFE_RTSMie_Fmatrix (4,max_mie_angles,2)

!  Distribution parameters
!  -----------------------

!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(kind=8), intent(in)  :: RTSMie_dist (5)
   real(kind=8), intent(in)  :: LPSD_RTSMie_dist (5,3)

!  Local

   integer       :: L, k, M, MS, cmask(6),fmask(4)
   character*1   :: c1
   real(kind=8)  :: pf

!  Old indexing

!   data cmask / 1, 3, 4, 6, 2, 5 /
!   data fmask / 1, 3, 2, 4 /

!  New indexing

   data cmask / 1, 2, 3, 4, 5, 6 /
   data fmask / 1, 2, 3, 4 /

!  open file

   open(35,file=TRIM(ADJUSTL(filename)),status = 'unknown' )

!  write distribution information

   if ( Do_monodisperse ) then
      write(35,'(/a)')' *** MONODISPERSE OUTPUT ONLY '
   else
      write(35,'(/a/)')' PSD output ------------- '
      write(35,'(A,T25,1pe15.6)')'Number density     = ',RTSMie_dist (1)
      write(35,'(A,T25,1pe15.6)')'Cross-section      = ',RTSMie_dist (2)
      write(35,'(A,T25,1pe15.6)')'Volume             = ',RTSMie_dist (3)
      write(35,'(A,T25,1pe15.6)')'Effective Radius   = ',RTSMie_dist (4)
      write(35,'(A,T25,1pe15.6)')'Effective Variance = ',RTSMie_dist (5)
   endif

!  write Bulk property output

   write(35,'(/a/)')' Bulk property output ------------- '
   write(35,'(A,T40,1pe20.11)')'Extinction Coefficient = ',RTSMie_bulk (1)
   write(35,'(A,T40,1pe20.11)')'Scattering Coefficient = ',RTSMie_bulk (2)
   write(35,'(A,T40,1pe20.11)')'Singlescattering albedo= ',RTSMie_bulk (3)
   write(35,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',RTSMie_asymm

!  Expansion coefficients
  
   if ( Do_expcoeffs ) then
      write(35,'(/a)')'Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
      write(35,'(a,T30,I5/)')' - Number of coefficients =  ',RTSMie_ncoeffs
      do L = 0, RTSMie_ncoeffs
        write(35,'(i5,1p6e20.11)')L,(RTSMie_expcoeffs (cmask(k),L),k=1,6)
      enddo
   endif

!  F-matrix

   if ( Do_Fmatrix ) then
      write(35,'(/a)')' F-matrix (F11=F22/F33=F44/F12/F34) @ angles (degs)'
      write(35,'(a,I5/)')' - Number of angles =  ',n_Fmatrix_angles
      do L = 1, n_Fmatrix_angles
        write(35,'(f6.2,4F20.11)')Fmatrix_angles(L),(RTSMie_Fmatrix (fmask(k),L),k=1,4)
      enddo
   endif

!  Close file

   close(35)

!  write PSD-linearized results
!    Filenames derived from Input name, with '_LPSD_*' added (where * = 1, 2, or 3)
 
   if ( .not. do_monodisperse ) then
      if ( Do_psd_linearize ) then
         MS = 0 ; if (psd_index.eq.4) MS = 2 ; if (psd_index.ne.4) MS = 3
         do M = 1, MS
            pf = 1.0d0 ; if (donorm) pf = psd_pars(m)
            write(C1,'(I1)')M
            open(36,file=TRIM(ADJUSTL(filename))//'_LPSD_'//C1,status = 'unknown' )

            write(36,'(/a/)')' LPSD Distribution output ------------- '
            write(36,'(A,T25,1pe15.6)')'Number density     = ',pf*LPSD_RTSMie_dist (1,m) 
            write(36,'(A,T25,1pe15.6)')'Cross-section      = ',pf*LPSD_RTSMie_dist (2,m)
            write(36,'(A,T25,1pe15.6)')'Volume             = ',pf*LPSD_RTSMie_dist (3,m)
            write(36,'(A,T25,1pe15.6)')'Effective Radius   = ',pf*LPSD_RTSMie_dist (4,m)
            write(36,'(A,T25,1pe15.6)')'Effective Variance = ',pf*LPSD_RTSMie_dist (5,m)

            write(36,'(/a/)')' LPSD Bulk property output ------------- '
            write(36,'(A,T40,1pe20.11)')'Extinction Coefficient = ',pf*LPSD_RTSMie_bulk (1,m)
            write(36,'(A,T40,1pe20.11)')'Scattering Coefficient = ',pf*LPSD_RTSMie_bulk (2,m)
            write(36,'(A,T40,1pe20.11)')'Singlescattering albedo= ',pf*LPSD_RTSMie_bulk (3,m)
            write(36,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',pf*LPSD_RTSMie_asymm(m)

            if ( Do_expcoeffs ) then
               write(36,'(/a)')'LPSD Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
               write(36,'(a,I5/)')' - Number of coefficients =  ',RTSMie_ncoeffs
               do L = 0, RTSMie_ncoeffs
                  write(36,'(i5,1p6e20.11)')L,(pf*LPSD_RTSMie_expcoeffs (cmask(k),L,M),k=1,6)
               enddo
            endif

            if ( Do_Fmatrix ) then
               write(36,'(/a)')' LPSD F-matrix (F11=F22/F33=F44/F12/F34) @ angles (degs)'
               write(36,'(a,I5/)')' - Number of angles =  ',n_Fmatrix_angles
               do L = 1, n_Fmatrix_angles
                  write(36,'(f6.2,4F20.11)')Fmatrix_angles(L),(pf*LPSD_RTSMie_Fmatrix (fmask(k),L,m),k=1,4)
               enddo
            endif

            close(36)

         enddo
      endif
   endif

!  write RFE-linearized results
!    Filenames derived from Input name, with '_LRFE_*' added (where * = 1, 2)

   if ( Do_LinearRef ) then
      do M = 1, 2
         pf = 1.0d0 ; if (donorm.and.m.eq.1) pf = n_real ; if (donorm.and.m.eq.2) pf = n_imag

         write(C1,'(I1)')M
         open(37,file=TRIM(ADJUSTL(filename))//'_LRFE_'//C1,status = 'unknown' )

         write(37,'(/a/)')' LRFE Bulk property output ------------- '
         write(37,'(A,T40,1pe20.11)')'Extinction Coefficient = ',pf*LRFE_RTSMie_bulk (1,m)
         write(37,'(A,T40,1pe20.11)')'Scattering Coefficient = ',pf*LRFE_RTSMie_bulk (2,m)
         write(37,'(A,T40,1pe20.11)')'Singlescattering albedo= ',pf*LRFE_RTSMie_bulk (3,m)
         write(37,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',pf*LRFE_RTSMie_asymm(m)

         if ( Do_expcoeffs ) then
            write(37,'(/a)')'LRFE Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
            write(37,'(a,I5/)')' - Number of coefficients =  ',RTSMie_ncoeffs
            do L = 0, RTSMie_ncoeffs
               write(37,'(i5,1p6e20.11)')L,(pf*LRFE_RTSMie_expcoeffs (cmask(k),L,M),k=1,6)
            enddo
         endif

         if ( Do_Fmatrix ) then
            write(37,'(/a)')' LRFE F-matrix (F11=F22/F33=F44/F12/F34) @ angles (degs)'
            write(37,'(a,I5/)')' - Number of angles =  ',n_Fmatrix_angles
            do L = 1, n_Fmatrix_angles
               write(37,'(f6.2,4F20.11)')Fmatrix_angles(L),(pf*LRFE_RTSMie_Fmatrix (fmask(k),L,m),k=1,4)
            enddo
         endif

         close(37)

      enddo
   endif

!  Finish

   return
end subroutine RTSMie_write_extended

subroutine RTSMie_write_bimodal_standard                 &
       ( filename, Do_Expcoeffs, Do_Fmatrix,             & ! Control Inputs
         n_Fmatrix_angles, Fmatrix_angles, fraction,     & ! Control Inputs
         RTSMie_bulk, RTSMie_asymm, RTSMie_ncoeffs,      & ! RTSMie Results
         RTSMie_expcoeffs, RTSMie_Fmatrix, RTSMie_dist )   ! RTSMie Results

   use RTSMie_parameters_m, ONLY : dp, max_Mie_angles

!  writes the Bi-modal output
!   --- NO MONODISPERSE here          !!!!!!!!!!!!!!

!  writes out both distributions

   implicit none

!  routine inputs
!  ==============

!  Filename

   character*(*), intent(in)  :: filename

!  Flags for Expansion Coefficient and Fmatrix calculations

   logical, intent(in)  :: Do_Expcoeffs
   logical, intent(in)  :: Do_Fmatrix

!  F-matrix Angular control input
!  Calculate F-matrix at user-defined angles (do_Fmatrix flag MUST BE set)
!       n_Fmatrix_angles = number of user-defined angles.
!       Fmatrix_angles   = user-defined angles, in DEGREES between [0, 180]

   INTEGER, intent(in)       :: n_Fmatrix_angles
   REAL(KIND=dp), intent(in) :: Fmatrix_angles(max_Mie_angles)

!  Fraction

   REAL(KIND=dp), intent(in) :: fraction

!  RTSMIE MASTER : OUTPUT ARGUMENTS
!  ================================

!  Bulk distribution parameters
!  ----------------------------

!    1 = Extinction Coefficient
!    2 = Scattering Coefficient
!    3 = Single scattering albedo

   real(kind=8), intent(in)   :: RTSMie_bulk (3)

!  Expansion coefficients and Asymmetry parameter
!  ----------------------------------------------

   integer      , intent(in)   :: RTSMie_ncoeffs
   real(kind=dp), intent(in)   :: RTSMie_expcoeffs (6,0:max_Mie_angles)
   real(kind=dp), intent(in)   :: RTSMie_asymm

!  F-matrix,  optional output
!  --------------------------

   real(kind=dp), intent(in)   :: RTSMie_Fmatrix(4,max_Mie_angles)

!  Distribution parameters
!  -----------------------

!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(kind=8), intent(in)   :: RTSMie_dist (5,2)

!  Local

   integer       :: L, k, cmask(6),fmask(4)

!  Old indexing

!   data cmask / 1, 3, 4, 6, 2, 5 /
!   data fmask / 1, 3, 2, 4 /

!  New indexing

   data cmask / 1, 2, 3, 4, 5, 6 /
   data fmask / 1, 2, 3, 4 /

!  open file

   open(35,file=TRIM(ADJUSTL(filename)),status = 'unknown' )

!  write distribution information

   write(35,'(/a/)')' PSD output,  first mode  ------------- '
   write(35,'(A,T25,1pe15.6)')'Number density     = ',RTSMie_dist (1,1)
   write(35,'(A,T25,1pe15.6)')'Cross-section      = ',RTSMie_dist (2,1)
   write(35,'(A,T25,1pe15.6)')'Volume             = ',RTSMie_dist (3,1)
   write(35,'(A,T25,1pe15.6)')'Effective Radius   = ',RTSMie_dist (4,1)
   write(35,'(A,T25,1pe15.6)')'Effective Variance = ',RTSMie_dist (5,1)

   write(35,'(/a/)')' PSD output, second mode  ------------- '
   write(35,'(A,T25,1pe15.6)')'Number density     = ',RTSMie_dist (1,2)
   write(35,'(A,T25,1pe15.6)')'Cross-section      = ',RTSMie_dist (2,2)
   write(35,'(A,T25,1pe15.6)')'Volume             = ',RTSMie_dist (3,2)
   write(35,'(A,T25,1pe15.6)')'Effective Radius   = ',RTSMie_dist (4,2)
   write(35,'(A,T25,1pe15.6)')'Effective Variance = ',RTSMie_dist (5,2)

   write(35,'(/a/)')' PSD output, Bimodal fraction  ------------- '
   write(35,'(A,T25,f10.5)')  'Bi-modal fraction  = ',fraction

!  write Bulk property output
!    ASYMM will be Zero is Do_expcoeffs not set.

   write(35,'(/a/)')' Total Bulk property output ------------- '
   write(35,'(A,T40,1pe20.11)')'Extinction Coefficient = ',RTSMie_bulk (1)
   write(35,'(A,T40,1pe20.11)')'Scattering Coefficient = ',RTSMie_bulk (2)
   write(35,'(A,T40,1pe20.11)')'Singlescattering albedo= ',RTSMie_bulk (3)
   write(35,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',RTSMie_asymm

!  Expansion coefficients
  
   if ( Do_expcoeffs ) then
      write(35,'(/a)')'Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
      write(35,'(a,T30,I5/)')' - Number of coefficients =  ',RTSMie_ncoeffs
      do L = 0, RTSMie_ncoeffs
        write(35,'(i5,1p6e20.11)')L,(RTSMie_expcoeffs (cmask(k),L),k=1,6)
      enddo
   endif

!  F-matrix

   if ( Do_Fmatrix ) then
      write(35,'(/a)')' F-matrix (F11=F22/F33=F44/F12/F34) @ angles (degs)'
      write(35,'(a,I5/)')' - Number of angles =  ',n_Fmatrix_angles
      do L = 1, n_Fmatrix_angles
        write(35,'(f6.2,4F20.11)')Fmatrix_angles(L),(RTSMie_Fmatrix (fmask(k),L),k=1,4)
      enddo
   endif

!  Close file

   close(35)

!  Finish

   return
end subroutine RTSMie_write_bimodal_standard

subroutine RTSMie_write_bimodal_extended                   &
       ( filename, Do_Expcoeffs, Do_Fmatrix,               & ! Control Inputs
         Do_LinearRef, Do_psd_Linearize, psd_index,        & ! Control Inputs
         n_Fmatrix_angles, Fmatrix_angles, fraction,       & ! Control Inputs
         RTSMie_bulk, RTSMie_asymm, RTSMie_ncoeffs,        & ! RTSMie Results
         RTSMie_expcoeffs, RTSMie_Fmatrix, RTSMie_dist,    & ! RTSMie Results
         LPSD_RTSMie_bulk, LPSD_RTSMie_asymm,              & ! RTSMier Linearized Results
         LPSD_RTSMie_expcoeffs, LPSD_RTSMie_Fmatrix,       & ! RTSMie Linearized Results
         LRFE_RTSMie_bulk, LRFE_RTSMie_asymm,              & ! RTSMie Linearized Results
         LRFE_RTSMie_expcoeffs, LRFE_RTSMie_Fmatrix,       & ! RTSMie Linearized Results
         LFRC_RTSMie_bulk, LFRC_RTSMie_asymm,              & ! Outputs (Bimodal frac)
         LFRC_RTSMie_expcoeffs, LFRC_RTSMie_Fmatrix,       & ! Outputs (Bimodal frac)
         LPSD_RTSMie_dist )                                  ! RTSMie Linearized Results

   implicit none

!  routine inputs
!  ==============

!  Filename

   character*(*), intent(in)  :: filename

!  Flags for Expansion Coefficient and Fmatrix calculations

   logical, intent(in)  :: Do_Expcoeffs
   logical, intent(in)  :: Do_Fmatrix

!  linearization w.r.t Ref Index

   logical, intent(in)  :: Do_LinearRef    ! w.r.t. Real/Imag parts ref index

!  Linearization for PSD parameters

   logical, intent(in)  :: Do_psd_linearize        ! PSD parameters

!  PSD index

   integer, intent(in)  :: psd_Index(2)

!  F-matrix Angular control input
!  Calculate F-matrix at user-defined angles (do_Fmatrix flag MUST BE set)
!       n_Fmatrix_angles = number of user-defined angles.
!       Fmatrix_angles   = user-defined angles, in DEGREES between [0, 180]

   INTEGER, intent(in)       :: n_Fmatrix_angles
   REAL(KIND=dp), intent(in) :: Fmatrix_angles(max_Mie_angles)

!  Fraction

   REAL(KIND=dp), intent(in) :: fraction

!  RTSMIE MASTER : OUTPUT ARGUMENTS
!  ================================

!  Bulk distribution parameters
!  ----------------------------

!    1 = Extinction Coefficient
!    2 = Scattering Coefficient
!    3 = Single scattering albedo

   real(kind=8), intent(in)   :: RTSMie_bulk (3)

!  linearizations w.r.t. PSD parameters

   real(kind=8), intent(in)  :: LPSD_RTSMie_bulk (3,3,2)

!  linearizations w.r.t. RefIdx parameters

   real(kind=8), intent(in)  :: LRFE_RTSMie_bulk (3,2,2)

!  Expansion coefficients and Asymmetry parameter
!  ----------------------------------------------

   integer     , intent(in)   :: RTSMie_ncoeffs
   real(kind=8), intent(in)   :: RTSMie_expcoeffs (6,0:max_mie_angles)
   real(kind=8), intent(in)   :: RTSMie_asymm

!  linearizations w.r.t. PSD parameters

   real(kind=8), intent(in)  :: LPSD_RTSMie_expcoeffs (6,0:max_mie_angles,3,2)
   real(kind=8), intent(in)  :: LPSD_RTSMie_asymm(3,2)

!  linearizations w.r.t. RefIdx parameters

   real(kind=8), intent(in)  :: LRFE_RTSMie_expcoeffs (6,0:max_mie_angles,2,2)
   real(kind=8), intent(in)  :: LRFE_RTSMie_asymm(2,2)

!  F-matrix,  optional output
!  --------------------------

   real(kind=8), intent(in)  :: RTSMie_Fmatrix (4,max_mie_angles)

!  Linearizations of F-matrix

   real(kind=8), intent(in)  :: LPSD_RTSMie_Fmatrix (4,max_mie_angles,3,2)
   real(kind=8), intent(in)  :: LRFE_RTSMie_Fmatrix (4,max_mie_angles,2,2)

!  linearizations w.r.t fraction
!  -----------------------------

   real(kind=8), intent(in)  :: LFRC_RTSMie_bulk (3)
   real(kind=8), intent(in)  :: LFRC_RTSMie_expcoeffs (6,0:max_mie_angles)
   real(kind=8), intent(in)  :: LFRC_RTSMie_asymm
   real(kind=8), intent(in)  :: LFRC_RTSMie_Fmatrix (4,max_mie_angles)

!  Distribution parameters
!  -----------------------

!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(kind=8), intent(in)  :: RTSMie_dist (5,2)
   real(kind=8), intent(in)  :: LPSD_RTSMie_dist (5,3,2)

!  Local

   integer       :: L, k, M, MS, J, cmask(6),fmask(4)
   character*1   :: c1, CJ
   character*22  :: CPSD

!  Old indexing

!   data cmask / 1, 3, 4, 6, 2, 5 /
!   data fmask / 1, 3, 2, 4 /

!  New indexing

   data cmask / 1, 2, 3, 4, 5, 6 /
   data fmask / 1, 2, 3, 4 /

!  open file

   open(35,file=TRIM(ADJUSTL(filename)),status = 'unknown' )

!  write distribution information

   write(35,'(/a/)')' PSD output,  first mode  ------------- '
   write(35,'(A,T25,1pe15.6)')'Number density     = ',RTSMie_dist (1,1)
   write(35,'(A,T25,1pe15.6)')'Cross-section      = ',RTSMie_dist (2,1)
   write(35,'(A,T25,1pe15.6)')'Volume             = ',RTSMie_dist (3,1)
   write(35,'(A,T25,1pe15.6)')'Effective Radius   = ',RTSMie_dist (4,1)
   write(35,'(A,T25,1pe15.6)')'Effective Variance = ',RTSMie_dist (5,1)

   write(35,'(/a/)')' PSD output, second mode  ------------- '
   write(35,'(A,T25,1pe15.6)')'Number density     = ',RTSMie_dist (1,2)
   write(35,'(A,T25,1pe15.6)')'Cross-section      = ',RTSMie_dist (2,2)
   write(35,'(A,T25,1pe15.6)')'Volume             = ',RTSMie_dist (3,2)
   write(35,'(A,T25,1pe15.6)')'Effective Radius   = ',RTSMie_dist (4,2)
   write(35,'(A,T25,1pe15.6)')'Effective Variance = ',RTSMie_dist (5,2)

   write(35,'(/a/)')' PSD output, Bimodal fraction  ------------- '
   write(35,'(A,T25,f10.5)')  'Bi-modal fraction  = ',fraction

!  write Bulk property output
!    ASYMM will be Zero is Do_expcoeffs not set.

   write(35,'(/a/)')' Total Bulk property output ------------- '
   write(35,'(A,T40,1pe20.11)')'Extinction Coefficient = ',RTSMie_bulk (1)
   write(35,'(A,T40,1pe20.11)')'Scattering Coefficient = ',RTSMie_bulk (2)
   write(35,'(A,T40,1pe20.11)')'Singlescattering albedo= ',RTSMie_bulk (3)
   write(35,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',RTSMie_asymm

!  Expansion coefficients
  
   if ( Do_expcoeffs ) then
      write(35,'(/a)')'Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
      write(35,'(a,T30,I5/)')' - Number of coefficients =  ',RTSMie_ncoeffs
      do L = 0, RTSMie_ncoeffs
        write(35,'(i5,1p6e20.11)')L,(RTSMie_expcoeffs (cmask(k),L),k=1,6)
      enddo
   endif

!  F-matrix

   if ( Do_Fmatrix ) then
      write(35,'(/a)')' F-matrix (F11=F22/F33=F44/F12/F34) @ angles (degs)'
      write(35,'(a,I5/)')' - Number of angles =  ',n_Fmatrix_angles
      do L = 1, n_Fmatrix_angles
        write(35,'(f6.2,4F20.11)')Fmatrix_angles(L),(RTSMie_Fmatrix (fmask(k),L),k=1,4)
      enddo
   endif

!  Close file

   close(35)

!  write PSD-linearized results
!    Filenames derived from Input name, with '_LPSD_*' added (where * = 1, 2, or 3)
 
   if ( Do_psd_linearize ) then
      do J = 1, 2
         write(CJ,'(I1)')J
         MS = 0 ; if (psd_index(j).eq.4) MS = 2 ; if (psd_index(j).ne.4) MS = 3
         do M = 1, MS 
            write(C1,'(I1)')M
            CPSD = 'PSD # '//CJ//', parameter # '//C1
            open(36,file=TRIM(ADJUSTL(filename))//'_LPSD_output_PSD#'//CJ//'_PAR#'//C1,status='unknown')
          
            write(36,'(/a/)')' LPSD Distribution output, PSD # '//CJ//' ------------- '
            write(36,'(A,T25,1pe15.6)')'Number density     = ',LPSD_RTSMie_dist (1,m,j) 
            write(36,'(A,T25,1pe15.6)')'Cross-section      = ',LPSD_RTSMie_dist (2,m,j)
            write(36,'(A,T25,1pe15.6)')'Volume             = ',LPSD_RTSMie_dist (3,m,j)
            write(36,'(A,T25,1pe15.6)')'Effective Radius   = ',LPSD_RTSMie_dist (4,m,j)
            write(36,'(A,T25,1pe15.6)')'Effective Variance = ',LPSD_RTSMie_dist (5,m,j)

            write(36,'(/a/)')' LPSD Bulk property output, PSD # '//CJ//' ------------- '
            write(36,'(A,T40,1pe20.11)')'Extinction Coefficient = ',LPSD_RTSMie_bulk (1,m,j)
            write(36,'(A,T40,1pe20.11)')'Scattering Coefficient = ',LPSD_RTSMie_bulk (2,m,j)
            write(36,'(A,T40,1pe20.11)')'Singlescattering albedo= ',LPSD_RTSMie_bulk (3,m,j)
            write(36,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',LPSD_RTSMie_asymm(m,j)

            if ( Do_expcoeffs ) then
               write(36,'(/a)')'LPSD Expansion coefficient output (a1,a2,a3,a4,b1,b2), PSD # '//CJ
               write(36,'(a,I5/)')' - Number of coefficients =  ',RTSMie_ncoeffs
               do L = 0, RTSMie_ncoeffs
                  write(36,'(i5,1p6e20.11)')L,(LPSD_RTSMie_expcoeffs (cmask(k),L,M,J),k=1,6)
               enddo
            endif

            if ( Do_Fmatrix ) then
               write(36,'(/a)')' LPSD F-matrix (F11=F22/F33=F44/F12/F34) @ angles (degs), PSD # '//CJ
               write(36,'(a,I5/)')' - Number of angles =  ',n_Fmatrix_angles
               do L = 1, n_Fmatrix_angles
                  write(36,'(f6.2,4F20.11)')Fmatrix_angles(L),(LPSD_RTSMie_Fmatrix (fmask(k),L,m,j),k=1,4)
               enddo
            endif

            close(36)
         enddo
      enddo
   endif

!  write RFE-linearized results
!    Filenames derived from Input name, with '_LRFE_*' added (where * = 1, 2)

   if ( Do_LinearRef ) then
      do J = 1, 2
         write(CJ,'(I1)')J
         MS = 2
         do M = 1, MS
            write(C1,'(I1)')M
            CPSD = 'PSD # '//CJ//', parameter # '//C1
            open(37,file=TRIM(ADJUSTL(filename))//'_LRFE_output_PSD#'//CJ//'_PAR#'//C1,status = 'unknown' )
          
            write(37,'(/a/)')' LRFE Bulk property output, RFE # '//CJ//' ------------- '
            write(37,'(A,T40,1pe20.11)')'Extinction Coefficient = ',LRFE_RTSMie_bulk (1,m,j)
            write(37,'(A,T40,1pe20.11)')'Scattering Coefficient = ',LRFE_RTSMie_bulk (2,m,j)
            write(37,'(A,T40,1pe20.11)')'Singlescattering albedo= ',LRFE_RTSMie_bulk (3,m,j)
            write(37,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',LRFE_RTSMie_asymm(m,j)

            if ( Do_expcoeffs ) then
               write(37,'(/a)')'LRFE Expansion coefficient output (a1,a2,a3,a4,b1,b2), RFE # '//CJ
               write(37,'(a,I5/)')' - Number of coefficients =  ',RTSMie_ncoeffs
               do L = 0, RTSMie_ncoeffs
                  write(37,'(i5,1p6e20.11)')L,(LRFE_RTSMie_expcoeffs (cmask(k),L,M,j),k=1,6)
               enddo
            endif

            if ( Do_Fmatrix ) then
               write(37,'(/a)')' LRFE F-matrix (F11=F22/F33=F44/F12/F34) @ angles (degs), RFE # '//CJ
               write(37,'(a,I5/)')' - Number of angles =  ',n_Fmatrix_angles
               do L = 1, n_Fmatrix_angles
                  write(37,'(f6.2,4F20.11)')Fmatrix_angles(L),(LRFE_RTSMie_Fmatrix (fmask(k),L,m,j),k=1,4)
               enddo
            endif

            close(37)

         enddo
      enddo
   endif

!  write FRC-linearized results (w.r.t. fraction)

   open(38,file=TRIM(ADJUSTL(filename))//'_LFRC',status = 'unknown' )

   write(38,'(/a/)')' LFRC Bulk property output ------------- '
   write(38,'(A,T40,1pe20.11)')'Extinction Coefficient = ',LFRC_RTSMie_bulk(1)
   write(38,'(A,T40,1pe20.11)')'Scattering Coefficient = ',LFRC_RTSMie_bulk(2)
   write(38,'(A,T40,1pe20.11)')'Singlescattering albedo= ',LFRC_RTSMie_bulk(3)
   write(38,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',LFRC_RTSMie_asymm

   if ( Do_expcoeffs ) then
      write(38,'(/a)')'LFRC Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
      write(38,'(a,I5/)')' - Number of coefficients =  ',RTSMie_ncoeffs
      do L = 0, RTSMie_ncoeffs
         write(38,'(i5,1p6e20.11)')L,(LFRC_RTSMie_expcoeffs (cmask(k),L),k=1,6)
      enddo
   endif

   if ( Do_Fmatrix ) then
      write(38,'(/a)')' LFRC F-matrix (F11=F22/F33=F44/F12/F34) @ angles (degs)'
      write(38,'(a,I5/)')' - Number of angles =  ',n_Fmatrix_angles
      do L = 1, n_Fmatrix_angles
         write(38,'(f6.2,4F20.11)')Fmatrix_angles(L),(LFRC_RTSMie_Fmatrix (fmask(k),L),k=1,4)
      enddo
   endif

   close(38)

!  Finish

   return
end subroutine RTSMie_write_bimodal_extended

!  Finish Module

end module RTSMie_IO_readwrite

