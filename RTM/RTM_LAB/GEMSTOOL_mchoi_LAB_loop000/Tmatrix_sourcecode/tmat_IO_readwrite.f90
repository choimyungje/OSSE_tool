module tmat_IO_readwrite

!  Auxiliary stand-alone routines for Reading and writing

!    tmat_read_configfile.          Reads the Regular configuration file
!    tmat_read_configbimodal.       Reads the Bimodal configuration file

!    tmat_write_standard.           Standard output
!    tmat_write_extended.           Standard + Linearized  output

!    tmat_write_bimodal_standard.   Standard output, bimodal
!    tmat_write_bimodal_extended.   Standard + Linearized output, bimodal

public

contains

subroutine tmat_read_configfile ( filename,        & ! File name
     Do_Expcoeffs, Do_Fmatrix,                     & ! Gp 1 Inputs (Flags)
     Do_Monodisperse, Do_EqSaSphere,               & ! Gp 1   Inputs (Flags)
     Do_LinearRef, Do_LinearEps, Do_psd_Linearize, & ! Gp 1 Inputs (Flags)
     Do_psd_OldStyle, psd_Index, psd_pars,         & ! Gp 1/2 Inputs (PSD)
     MonoRadius, R1, R2, FixR1R2,                  & ! Gp 2   Inputs (PSD)
     np, nkmax, npna, ndgs, eps, accuracy,         & ! Gp 3 Inputs (General)
     lambda, n_real, n_imag )                        ! Gp 4 Inputs (Optical)

   implicit none

!  Filename
!  --------

   character*(*), intent(in)  :: filename

!  TMATRIX MASTER : INPUT ARGUMENTS
!  ================================

!  Boolean Input arguments
!  -----------------------

!  Flags for Expansion Coefficient and Fmatrix calculations

   logical, intent(out)  :: Do_Expcoeffs
   logical, intent(out)  :: Do_Fmatrix

!  Logical flag for Monodisperse calculation

   logical, intent(out)  :: Do_monodisperse

!  Logical flag for using equal-surface-area-specification (ESAS)

   logical, intent(out)  :: Do_EqSaSphere

!  style flag, PSD

   logical, intent(out)  :: Do_psd_OldStyle

!  linearization of monodisperse quantitites

   logical, intent(out)  :: Do_LinearRef   ! w.r.t. Real/Imag parts ref index
   logical, intent(out)  :: Do_LinearEps   ! w.r.t. Shape factor

!  Linaerization for PSD parameters

   logical, intent(out)  :: Do_psd_linearize        ! PSD parameters

!  General Input arguments
!  -----------------------

!  integers

   integer     , intent(out)   :: np, nkmax, ndgs, npna

!  Accuracy and aspect ratio

   real(kind=8), intent(out)   :: accuracy, eps

!  Optical: Wavelength, refractive index
!  -------------------------------------

   real(kind=8), intent(out)   :: lambda, n_real, n_imag

!  PSD inputs
!  ----------

!  Flag for making an internal Fix of R1 and R2
!    ( Not relevant for the Old-style PSDs)

   logical     , intent(out)   :: FixR1R2

!  Monodisperse radius

   real(kind=8), intent(out)   :: Monoradius

!  R1 and R2 (intent(inout) from the new program)

   real(kind=8), intent(out)   :: R1, R2

!  PSD index and parameters

   integer     , intent(out)   :: psd_Index
   real(kind=8), intent(out)   :: psd_pars (3)

!  Do it

   open(1,file=TRIM(ADJUSTL(filename)),status = 'old')
   read(1,*)
   read(1,*)
   read(1,*) ! First group (Boolean flags)
   read(1,*)
   read(1,*) Do_expcoeffs
   read(1,*) Do_Fmatrix
   read(1,*) Do_monodisperse
   read(1,*) Do_EqSaSphere
   read(1,*) Do_LinearRef
   read(1,*) Do_LinearEps
   read(1,*) Do_psd_Linearize
   read(1,*) Do_psd_OldStyle
   read(1,*)
   read(1,*) ! Second group (PSD-related variables)
   read(1,*)
   read(1,*) psd_Index
   read(1,*) psd_pars(1)
   read(1,*) psd_pars(2)
   read(1,*) psd_pars(3)
   read(1,*) Monoradius
   read(1,*) FixR1R2
   read(1,*) R1
   read(1,*) R2
   read(1,*)
   read(1,*) ! Third group (General inputs)
   read(1,*)
   read(1,*) np
   read(1,*) nkmax
   read(1,*) npna
   read(1,*) ndgs
   read(1,*) eps
   read(1,*) accuracy
   read(1,*)
   read(1,*) ! Fourth group (Optical inputs)
   read(1,*)
   read(1,*) lambda
   read(1,*) n_real
   read(1,*) n_imag
   close(1)

!  return

   return
end subroutine tmat_read_configfile

subroutine tmat_read_configbimodal ( filename,     & ! File name
     Do_Expcoeffs, Do_Fmatrix,                     & ! Gp 1 Inputs (Flags)
     Do_Monodisperse, Do_EqSaSphere,               & ! Gp 1 Inputs (Flags)
     Do_LinearRef, Do_LinearEps, Do_psd_Linearize, & ! Gp 1 Inputs (Flags)
     Do_psd_OldStyle, psd_Index, psd_pars,         & ! Gp 1/2 Inputs (PSD)
     MonoRadius, R1, R2, FixR1R2, fraction,        & ! Gp 2   Inputs (PSD)
     np, nkmax, npna, ndgs, eps, accuracy,         & ! Gp 3 Inputs (General)
     lambda, n_real, n_imag )                        ! Gp 4 Inputs (Optical)

   implicit none

!  Filename
!  --------

   character*(*), intent(in)  :: filename

!  TMATRIX MASTER : INPUT ARGUMENTS
!  ================================

!  Boolean Input arguments
!  -----------------------

!  Flags for Expansion Coefficient and Fmatrix calculations

   logical, intent(out)  :: Do_Expcoeffs
   logical, intent(out)  :: Do_Fmatrix

!  Logical flag for Monodisperse calculation

   logical, intent(out)  :: Do_monodisperse

!  Logical flag for using equal-surface-area-specification (ESAS)

   logical, intent(out)  :: Do_EqSaSphere

!  style flag, PSD

   logical, intent(out)  :: Do_psd_OldStyle

!  linearization of monodisperse quantitites

   logical, intent(out)  :: Do_LinearRef            ! Real part refractive index
   logical, intent(out)  :: Do_LinearEps            ! Imag part refractive index

!  Linaerization for PSD parameters

   logical, intent(out)  :: Do_psd_linearize        ! PSD parameters

!  General Input arguments
!  -----------------------

!  integers

   integer     , intent(out)   :: np, nkmax(2), ndgs(2), npna

!  Accuracy and aspect ratio

   real(kind=8), intent(out)   :: accuracy, eps(2)

!  Optical: Wavelength, refractive index, 2 distributions
!  -------------------------------------

   real(kind=8), intent(out)   :: lambda, n_real(2), n_imag(2)

!  PSD inputs
!  ----------

!  Flag for making an internal Fix of R1 and R2
!    ( Not relevant for the Old-style PSDs)

   logical     , intent(out)   :: FixR1R2(2)

!  Monodisperse radius

   real(kind=8), intent(out)   :: Monoradius

!  R1 and R2 (intent(inout) from the new program)  2 distributions

   real(kind=8), intent(out)   :: R1(2), R2(2)

!  PSD index and parameters. 2 distributions of same type.

   integer     , intent(out)   :: psd_Index(2)
   real(kind=8), intent(out)   :: psd_pars (3,2)

!  Fractional difference .

   real(kind=8), intent(out)   :: fraction

!  Do it

   open(1,file=TRIM(ADJUSTL(filename)),status = 'old')
   read(1,*)
   read(1,*)
   read(1,*) ! First group (Boolean flags)
   read(1,*)
   read(1,*) Do_expcoeffs
   read(1,*) Do_Fmatrix
   read(1,*) Do_monodisperse
   read(1,*) Do_EqSaSphere
   read(1,*) Do_LinearRef
   read(1,*) Do_LinearEps
   read(1,*) Do_psd_Linearize
   read(1,*) Do_psd_OldStyle
   read(1,*)
   read(1,*) ! Second group (PSD-related variables)
   read(1,*)
   read(1,*) psd_Index(1) , psd_index(2)
   read(1,*) psd_pars(1,1), psd_pars(1,2)
   read(1,*) psd_pars(2,1), psd_pars(2,2)
   read(1,*) psd_pars(3,1), psd_pars(3,2)
   read(1,*) Monoradius
   read(1,*) FixR1R2(2)
   read(1,*) R1(1),R1(2)
   read(1,*) R2(1),R2(2)
   read(1,*)
   read(1,*) ! Third group (General inputs)
   read(1,*)
   read(1,*) np
   read(1,*) nkmax(1),nkmax(2)
   read(1,*) npna
   read(1,*) ndgs(1),ndgs(2)
   read(1,*) eps(1),eps(2)
   read(1,*) accuracy
   read(1,*)
   read(1,*) ! Fourth group (Optical inputs)
   read(1,*)
   read(1,*) lambda
   read(1,*) n_real(1),n_real(2)
   read(1,*) n_imag(1),n_imag(2)
   read(1,*)
   read(1,*) ! Fifth group (fractional input)
   read(1,*)
   read(1,*) fraction
   close(1)

!  Some settings fixed

   do_monodisperse = .false.
   Do_psd_OldStyle = .false.
   FixR1R2    = .false.
   Monoradius = 0.0d0

!  return

   return
end subroutine tmat_read_configbimodal


subroutine tmat_write_standard                     &
       ( filename, Do_Expcoeffs, Do_Fmatrix,       & ! Control Inputs
         Do_monodisperse, maxnpa, npl1, npna,      & ! Control Inputs
         tmat_bulk, tmat_asymm, tmat_ncoeffs,      & ! Tmatrix Results
         tmat_expcoeffs, tmat_Fmatrix, Tmat_dist )   ! Tmatrix Results

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

!  Dimensioning

   integer, intent(in)  :: MAXNPA, NPL1

!  Number of scattering angles (Fmatrix only)

   integer, intent(in)  :: npna

!  TMATRIX MASTER : OUTPUT ARGUMENTS
!  =================================

!  Bulk distribution parameters
!  ----------------------------

!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(kind=8), intent(in)   :: Tmat_bulk (3)

!  Expansion coefficients and Asymmetry parameter
!  ----------------------------------------------

   integer     , intent(in)   :: tmat_ncoeffs
   real(kind=8), intent(in)   :: Tmat_expcoeffs (NPL1,6)
   real(kind=8), intent(in)   :: Tmat_asymm

!  F-matrix,  optional output
!  --------------------------

   real(kind=8), intent(in)   :: Tmat_Fmatrix (MAXNPA,6)

!  Distribution parameters
!  -----------------------

!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(kind=8), intent(in)   :: Tmat_dist (5)

!  Local

   integer       :: L, k
   real(kind=8)  :: tb, db

!  open file

   open(35,file=TRIM(ADJUSTL(filename)),status = 'unknown' )

!  write distribution information

   if ( Do_monodisperse ) then
      write(35,'(/a)')' *** MONODISPERSE OUTPUT ONLY '
   else
      write(35,'(/a/)')' PSD output ------------- '
      write(35,'(A,T25,1pe15.6)')'Number density     = ',Tmat_dist (1)
      write(35,'(A,T25,1pe15.6)')'Cross-section      = ',Tmat_dist (2)
      write(35,'(A,T25,1pe15.6)')'Volume             = ',Tmat_dist (3)
      write(35,'(A,T25,1pe15.6)')'Effective Radius   = ',Tmat_dist (4)
      write(35,'(A,T25,1pe15.6)')'Effective Variance = ',Tmat_dist (5)
   endif

!  write Bulk property output
!    ASYMM will be Zero is Do_expcoeffs not set.

   write(35,'(/a/)')' Bulk property output ------------- '
   write(35,'(A,T40,1pe20.11)')'Extinction coefficient = ',Tmat_bulk (1)
   write(35,'(A,T40,1pe20.11)')'Scattering coefficient = ',Tmat_bulk (2)
   write(35,'(A,T40,1pe20.11)')'Singlescattering albedo= ',Tmat_bulk (3)
   write(35,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',Tmat_asymm

!  Expansion coefficients
  
   if ( Do_expcoeffs ) then
      write(35,'(/a)')'Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
      write(35,'(a,T30,I5/)')' - Number of coefficients =  ',tmat_ncoeffs
      do L = 1, tmat_ncoeffs
        write(35,'(i5,1p6e20.11)')L-1,(Tmat_expcoeffs (L,k),k=1,6)
      enddo
   endif

!  F-matrix for npna angles (regularly spaced)

   if ( Do_Fmatrix ) then
      write(35,'(/a)')' F-matrix (F11/F22/F33/F44/F12/F34) @ angles (degs)'
      write(35,'(a,I5/)')' - Number of angles =  ',npna
      db = 180.0/dble(npna-1); tb = - db
      do L = 1, npna
        tb = tb + db
        write(35,'(f6.2,6F20.11)')tb,(tmat_Fmatrix (L,k),k=1,6)
      enddo
   endif

!  Close file

   close(35)

!  Finish

   return
end subroutine tmat_write_standard

subroutine tmat_write_extended                         &
       ( filename, Do_Expcoeffs, Do_Fmatrix,           & ! Control Inputs
         Do_Monodisperse, Do_LinearRef, Do_LinearEps,  & ! Gp 1 Inputs (Flags)
         Do_psd_OldStyle, Do_psd_Linearize, psd_index, & ! Gp 1/2 Inputs (PSD)
         maxnpa, npl1, npna,                           & ! Control Inputs
         tmat_bulk, tmat_asymm, tmat_ncoeffs,          & ! Tmatrix Results
         tmat_expcoeffs, tmat_Fmatrix, Tmat_dist,      & ! Tmatrix Results
         LPSD_tmat_bulk, LPSD_tmat_asymm,              & ! Tmatrix Linearized Results
         LPSD_tmat_expcoeffs, LPSD_tmat_Fmatrix,       & ! Tmatrix Linearized Results
         LRFE_tmat_bulk, LRFE_tmat_asymm,              & ! Tmatrix Linearized Results
         LRFE_tmat_expcoeffs, LRFE_tmat_Fmatrix,       & ! Tmatrix Linearized Results
         LPSD_Tmat_dist )                                ! Tmatrix Linearized Results

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

!  linearization of monodisperse quantitites

   logical, intent(in)  :: Do_LinearRef            ! Real part refractive index
   logical, intent(in)  :: Do_LinearEps            ! Imag part refractive index

!  style flag, PSD

   logical, intent(in)  :: Do_psd_OldStyle

!  Linearization for PSD parameters

   logical, intent(in)  :: Do_psd_linearize        ! PSD parameters

!  PSD index

   integer, intent(in)  :: psd_Index

!  Dimensioning

   integer, intent(in)  :: MAXNPA, NPL1

!  Number of scattering angles (Fmatrix only)

   integer, intent(in)  :: npna

!  TMATRIX MASTER : OUTPUT ARGUMENTS
!  =================================

!  Bulk distribution parameters
!  ----------------------------

!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(kind=8), intent(in)   :: Tmat_bulk (3)

!  linearizations w.r.t. PSD parameters

   real(kind=8), intent(in)  :: LPSD_Tmat_bulk (3,3)

!  linearizations w.r.t. RefIdx/Eps parameters

   real(kind=8), intent(in)  :: LRFE_Tmat_bulk (3,3)

!  Expansion coefficients and Asymmetry parameter
!  ----------------------------------------------

   integer     , intent(in)   :: tmat_ncoeffs
   real(kind=8), intent(in)   :: Tmat_expcoeffs (NPL1,6)
   real(kind=8), intent(in)   :: Tmat_asymm

!  linearizations w.r.t. PSD parameters

   real(kind=8), intent(in)  :: LPSD_Tmat_expcoeffs (NPL1,6,3)
   real(kind=8), intent(in)  :: LPSD_Tmat_asymm(3)

!  linearizations w.r.t. RefIdx/Eps parameters

   real(kind=8), intent(in)  :: LRFE_Tmat_expcoeffs (NPL1,6,3)
   real(kind=8), intent(in)  :: LRFE_Tmat_asymm(3)

!  F-matrix,  optional output
!  --------------------------

   real(kind=8), intent(in)   :: Tmat_Fmatrix (MAXNPA,6)

!  Linearizations of F-matrix

   real(kind=8), intent(in)  :: LPSD_Tmat_Fmatrix (MAXNPA,6,3)
   real(kind=8), intent(in)  :: LRFE_Tmat_Fmatrix (MAXNPA,6,3)

!  Distribution parameters
!  -----------------------

!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(kind=8), intent(in)  :: Tmat_dist (5)
   real(kind=8), intent(in)  :: LPSD_Tmat_dist (5,3)

!  Local

   integer       :: L, k, M, MS
   real(kind=8)  :: tb, db
   character*1   :: c1

!  open file

   open(35,file=TRIM(ADJUSTL(filename)),status = 'unknown' )

!  write distribution information

   if ( Do_monodisperse ) then
      write(35,'(/a)')' *** MONODISPERSE OUTPUT ONLY '
   else
      write(35,'(/a/)')' PSD output ------------- '
      write(35,'(A,T25,1pe15.6)')'Number density     = ',Tmat_dist (1)
      write(35,'(A,T25,1pe15.6)')'Cross-section      = ',Tmat_dist (2)
      write(35,'(A,T25,1pe15.6)')'Volume             = ',Tmat_dist (3)
      write(35,'(A,T25,1pe15.6)')'Effective Radius   = ',Tmat_dist (4)
      write(35,'(A,T25,1pe15.6)')'Effective Variance = ',Tmat_dist (5)
   endif

!  write Bulk property output
!    ASYMM will be Zero is Do_expcoeffs not set.

   write(35,'(/a/)')' Bulk property output ------------- '
   write(35,'(A,T40,1pe20.11)')'Extinction coefficient = ',Tmat_bulk (1)
   write(35,'(A,T40,1pe20.11)')'Scattering coefficient = ',Tmat_bulk (2)
   write(35,'(A,T40,1pe20.11)')'Singlescattering albedo= ',Tmat_bulk (3)
   write(35,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',Tmat_asymm

!  Expansion coefficients
  
   if ( Do_expcoeffs ) then
      write(35,'(/a)')'Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
      write(35,'(a,T30,I5/)')' - Number of coefficients =  ',tmat_ncoeffs
      do L = 1, tmat_ncoeffs
        write(35,'(i5,1p6e20.11)')L-1,(Tmat_expcoeffs (L,k),k=1,6)
      enddo
   endif

!  F-matrix for npna angles (regularly spaced)

   if ( Do_Fmatrix ) then
      write(35,'(/a)')' F-matrix (F11/F22/F33/F44/F12/F34) @ angles (degs)'
      write(35,'(a,I5/)')' - Number of angles =  ',npna
      db = 180.0/dble(npna-1); tb = - db
      do L = 1, npna
        tb = tb + db
        write(35,'(f6.2,6F20.11)')tb,(tmat_Fmatrix (L,k),k=1,6)
      enddo
   endif

!  Close file

   close(35)

!  write PSD-linearized results
!    Filenames derived from Input name, with '_LPSD_*' added (where * = 1, 2, or 3)
 
   if ( .not. do_monodisperse ) then
      if ( .not. Do_psd_OldStyle .and. Do_psd_linearize ) then
         MS = 0 ; if (psd_index.eq.4) MS = 2 ; if (psd_index.ne.4) MS = 3
         do M = 1, MS
            write(C1,'(I1)')M
            open(36,file=TRIM(ADJUSTL(filename))//'_LPSD_'//C1,status = 'unknown' )

            write(36,'(/a/)')' LPSD Distribution output ------------- '
            write(36,'(A,T25,1pe15.6)')'Number density     = ',LPSD_Tmat_dist (1,m) 
            write(36,'(A,T25,1pe15.6)')'Cross-section      = ',LPSD_Tmat_dist (2,m)
            write(36,'(A,T25,1pe15.6)')'Volume             = ',LPSD_Tmat_dist (3,m)
            write(36,'(A,T25,1pe15.6)')'Effective Radius   = ',LPSD_Tmat_dist (4,m)
            write(36,'(A,T25,1pe15.6)')'Effective Variance = ',LPSD_Tmat_dist (5,m)

            write(36,'(/a/)')' LPSD Bulk property output ------------- '
            write(36,'(A,T40,1pe20.11)')'Extinction coefficient = ',LPSD_Tmat_bulk (1,m)
            write(36,'(A,T40,1pe20.11)')'Scattering coefficient = ',LPSD_Tmat_bulk (2,m)
            write(36,'(A,T40,1pe20.11)')'Singlescattering albedo= ',LPSD_Tmat_bulk (3,m)
            write(36,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',LPSD_Tmat_asymm(m)

            if ( Do_expcoeffs ) then
               write(36,'(/a)')'LPSD Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
               write(36,'(a,I5/)')' - Number of coefficients =  ',tmat_ncoeffs
               do L = 1, tmat_ncoeffs
                  write(36,'(i5,1p6e20.11)')L-1,(LPSD_Tmat_expcoeffs (L,k,M),k=1,6)
               enddo
            endif

            if ( Do_Fmatrix ) then
               write(36,'(/a)')' LPSD F-matrix (F11/F22/F33/F44/F12/F34) @ angles (degs)'
               write(36,'(a,I5/)')' - Number of angles =  ',npna
               db = 180.0/dble(npna-1); tb = - db
               do L = 1, npna
                  tb = tb + db
                  write(36,'(f6.2,6F20.11)')tb,(LPSD_tmat_Fmatrix (L,k,m),k=1,6)
               enddo
            endif
            close(36)
         enddo
      endif
   endif

!  write RFE-linearized results
!    Filenames derived from Input name, with '_LRFE_*' added (where * = 1, 2, or 3)

   if ( Do_LinearRef.or.Do_LinearEps ) then
      MS = 0 ; if (Do_LinearRef) MS = MS + 2 ; if (Do_LinearEps) MS = MS + 1
      do M = 1, MS
         write(C1,'(I1)')M
         open(37,file=TRIM(ADJUSTL(filename))//'_LRFE_'//C1,status = 'unknown' )

         write(37,'(/a/)')' LRFE Bulk property output ------------- '
         write(37,'(A,T40,1pe20.11)')'Extinction coefficient = ',LRFE_Tmat_bulk (1,m)
         write(37,'(A,T40,1pe20.11)')'Scattering coefficient = ',LRFE_Tmat_bulk (2,m)
         write(37,'(A,T40,1pe20.11)')'Singlescattering albedo= ',LRFE_Tmat_bulk (3,m)
         write(37,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',LRFE_Tmat_asymm(m)

         if ( Do_expcoeffs ) then
            write(37,'(/a)')'LRFE Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
            write(37,'(a,I5/)')' - Number of coefficients =  ',tmat_ncoeffs
            do L = 1, tmat_ncoeffs
               write(37,'(i5,1p6e20.11)')L-1,(LRFE_Tmat_expcoeffs (L,k,M),k=1,6)
            enddo
         endif

         if ( Do_Fmatrix ) then
            write(37,'(/a)')' LRFE F-matrix (F11/F22/F33/F44/F12/F34) @ angles (degs)'
            write(37,'(a,I5/)')' - Number of angles =  ',npna
            db = 180.0/dble(npna-1); tb = - db
            do L = 1, npna
               tb = tb + db
               write(37,'(f6.2,6F20.11)')tb,(LRFE_tmat_Fmatrix (L,k,m),k=1,6)
            enddo
         endif

         close(37)
      enddo
   endif

!  Finish

   return
end subroutine tmat_write_extended

subroutine tmat_write_bimodal_standard             &
       ( filename, Do_Expcoeffs, Do_Fmatrix,       & ! Control Inputs
         maxnpa, npl1, npna, fraction,             & ! Control Inputs
         tmat_bulk, tmat_asymm, tmat_ncoeffs,      & ! Tmatrix Results
         tmat_expcoeffs, tmat_Fmatrix, Tmat_dist )   ! Tmatrix Results

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

!  Dimensioning

   integer, intent(in)  :: MAXNPA, NPL1

!  Number of scattering angles (Fmatrix only)

   integer, intent(in)  :: npna

!  TMATRIX MASTER : OUTPUT ARGUMENTS
!  =================================

!  Bulk distribution parameters
!  ----------------------------

!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(kind=8), intent(in)   :: Tmat_bulk (3)

!  Expansion coefficients and Asymmetry parameter
!  ----------------------------------------------

   integer     , intent(in)   :: tmat_ncoeffs
   real(kind=8), intent(in)   :: Tmat_expcoeffs (NPL1,6)
   real(kind=8), intent(in)   :: Tmat_asymm

!  F-matrix,  optional output
!  --------------------------

   real(kind=8), intent(in)   :: Tmat_Fmatrix (MAXNPA,6)

!  Distribution parameters
!  -----------------------

!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(kind=8), intent(in)   :: Tmat_dist (5,2)
   real(kind=8), intent(in)   :: fraction

!  Local

   integer       :: L, k
   real(kind=8)  :: tb, db

!  open file

   open(35,file=TRIM(ADJUSTL(filename)),status = 'unknown' )

!  write distribution information

   write(35,'(/a/)')' PSD output,  first mode  ------------- '
   write(35,'(A,T25,1pe15.6)')'Number density     = ',Tmat_dist (1,1)
   write(35,'(A,T25,1pe15.6)')'Cross-section      = ',Tmat_dist (2,1)
   write(35,'(A,T25,1pe15.6)')'Volume             = ',Tmat_dist (3,1)
   write(35,'(A,T25,1pe15.6)')'Effective Radius   = ',Tmat_dist (4,1)
   write(35,'(A,T25,1pe15.6)')'Effective Variance = ',Tmat_dist (5,1)

   write(35,'(/a/)')' PSD output, second mode  ------------- '
   write(35,'(A,T25,1pe15.6)')'Number density     = ',Tmat_dist (1,2)
   write(35,'(A,T25,1pe15.6)')'Cross-section      = ',Tmat_dist (2,2)
   write(35,'(A,T25,1pe15.6)')'Volume             = ',Tmat_dist (3,2)
   write(35,'(A,T25,1pe15.6)')'Effective Radius   = ',Tmat_dist (4,2)
   write(35,'(A,T25,1pe15.6)')'Effective Variance = ',Tmat_dist (5,2)

   write(35,'(/a/)')' PSD output, Bimodal fraction  ------------- '
   write(35,'(A,T25,f10.5)')  'Bi-modal fraction  = ',fraction

!  write Bulk property output
!    ASYMM will be Zero is Do_expcoeffs not set.

   write(35,'(/a/)')' Total Bulk property output ------------- '
   write(35,'(A,T40,1pe20.11)')'Extinction coefficient = ',Tmat_bulk (1)
   write(35,'(A,T40,1pe20.11)')'Scattering coefficient = ',Tmat_bulk (2)
   write(35,'(A,T40,1pe20.11)')'Singlescattering albedo= ',Tmat_bulk (3)
   write(35,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',Tmat_asymm

!  Expansion coefficients
  
   if ( Do_expcoeffs ) then
      write(35,'(/a)')'Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
      write(35,'(a,T30,I5/)')' - Number of coefficients =  ',tmat_ncoeffs
      do L = 1, tmat_ncoeffs
        write(35,'(i5,1p6e20.11)')L-1,(Tmat_expcoeffs (L,k),k=1,6)
      enddo
   endif

!  F-matrix for npna angles (regularly spaced)

   if ( Do_Fmatrix ) then
      write(35,'(/a)')' F-matrix (F11/F22/F33/F44/F12/F34) @ angles (degs)'
      write(35,'(a,I5/)')' - Number of angles =  ',npna
      db = 180.0/dble(npna-1); tb = - db
      do L = 1, npna
        tb = tb + db
        write(35,'(f6.2,6F20.11)')tb,(tmat_Fmatrix (L,k),k=1,6)
      enddo
   endif

!  Close file

   close(35)

!  Finish

   return
end subroutine tmat_write_bimodal_standard

subroutine tmat_write_bimodal_extended                 &
       ( filename, Do_Expcoeffs, Do_Fmatrix,           & ! Control Inputs
         Do_LinearRef, Do_LinearEps,                   & ! Gp 1 Inputs (Flags)
         Do_psd_Linearize, psd_index,                  & ! Gp 1/2 Inputs (PSD)
         maxnpa, npl1, npna, fraction,                 & ! Control Inputs
         tmat_bulk, tmat_asymm, tmat_ncoeffs,          & ! Tmatrix Results
         tmat_expcoeffs, tmat_Fmatrix, Tmat_dist,      & ! Tmatrix Results
         LPSD_tmat_bulk, LPSD_tmat_asymm,              & ! Tmatrix Linearized Results
         LPSD_tmat_expcoeffs, LPSD_tmat_Fmatrix,       & ! Tmatrix Linearized Results
         LRFE_tmat_bulk, LRFE_tmat_asymm,              & ! Tmatrix Linearized Results
         LRFE_tmat_expcoeffs, LRFE_tmat_Fmatrix,       & ! Tmatrix Linearized Results
         LFRC_tmat_bulk, LFRC_tmat_asymm,              & ! Outputs (Bimodal frac)
         LFRC_tmat_expcoeffs, LFRC_tmat_Fmatrix,       & ! Outputs (Bimodal frac)
         LPSD_Tmat_dist )                                ! Tmatrix Linearized Results

   implicit none

!  routine inputs
!  ==============

!  Filename

   character*(*), intent(in)  :: filename

!  Flags for Expansion Coefficient and Fmatrix calculations

   logical, intent(in)  :: Do_Expcoeffs
   logical, intent(in)  :: Do_Fmatrix

!  linearization of monodisperse quantitites

   logical, intent(in)  :: Do_LinearRef            ! Real part refractive index
   logical, intent(in)  :: Do_LinearEps            ! Imag part refractive index

!  Linearization for PSD parameters

   logical, intent(in)  :: Do_psd_linearize        ! PSD parameters

!  PSD index

   integer, intent(in)  :: psd_Index(2)

!  Dimensioning

   integer, intent(in)  :: MAXNPA, NPL1

!  Number of scattering angles (Fmatrix only)

   integer, intent(in)  :: npna

!  Fraction

   real(kind=8), intent(in)  :: fraction

!  TMATRIX MASTER : OUTPUT ARGUMENTS
!  =================================

!  Bulk distribution parameters
!  ----------------------------

!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(kind=8), intent(in)   :: Tmat_bulk (3)

!  linearizations w.r.t. PSD parameters

   real(kind=8), intent(in)  :: LPSD_Tmat_bulk (3,3,2)

!  linearizations w.r.t. RefIdx/Eps parameters

   real(kind=8), intent(in)  :: LRFE_Tmat_bulk (3,3,2)

!  Expansion coefficients and Asymmetry parameter
!  ----------------------------------------------

   integer     , intent(in)   :: tmat_ncoeffs
   real(kind=8), intent(in)   :: Tmat_expcoeffs (NPL1,6)
   real(kind=8), intent(in)   :: Tmat_asymm

!  linearizations w.r.t. PSD parameters

   real(kind=8), intent(in)  :: LPSD_Tmat_expcoeffs (NPL1,6,3,2)
   real(kind=8), intent(in)  :: LPSD_Tmat_asymm(3,2)

!  linearizations w.r.t. RefIdx/Eps parameters

   real(kind=8), intent(in)  :: LRFE_Tmat_expcoeffs (NPL1,6,3,2)
   real(kind=8), intent(in)  :: LRFE_Tmat_asymm(3,2)

!  F-matrix,  optional output
!  --------------------------

   real(kind=8), intent(in)  :: Tmat_Fmatrix (MAXNPA,6)

!  Linearizations of F-matrix

   real(kind=8), intent(in)  :: LPSD_Tmat_Fmatrix (MAXNPA,6,3,2)
   real(kind=8), intent(in)  :: LRFE_Tmat_Fmatrix (MAXNPA,6,3,2)

!  Fraction Jacobian
!  -----------------

   real(kind=8), intent(in)  :: LFRC_Tmat_bulk (3)
   real(kind=8), intent(in)  :: LFRC_Tmat_expcoeffs (NPL1,6)
   real(kind=8), intent(in)  :: LFRC_Tmat_asymm
   real(kind=8), intent(in)  :: LFRC_Tmat_Fmatrix (MAXNPA,6)

!  Distribution parameters
!  -----------------------

!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(kind=8), intent(in)  :: Tmat_dist (5,2)
   real(kind=8), intent(in)  :: LPSD_Tmat_dist (5,3,2)

!  Local

   integer       :: L, k, M, MS, J
   real(kind=8)  :: tb, db
   character*1   :: c1, CJ
   character*22  :: CPSD

!  open file

   open(35,file=TRIM(ADJUSTL(filename)),status = 'unknown' )

!  write distribution information

   write(35,'(/a/)')' PSD output,  first mode  ------------- '
   write(35,'(A,T25,1pe15.6)')'Number density     = ',Tmat_dist (1,1)
   write(35,'(A,T25,1pe15.6)')'Cross-section      = ',Tmat_dist (2,1)
   write(35,'(A,T25,1pe15.6)')'Volume             = ',Tmat_dist (3,1)
   write(35,'(A,T25,1pe15.6)')'Effective Radius   = ',Tmat_dist (4,1)
   write(35,'(A,T25,1pe15.6)')'Effective Variance = ',Tmat_dist (5,1)

   write(35,'(/a/)')' PSD output, second mode  ------------- '
   write(35,'(A,T25,1pe15.6)')'Number density     = ',Tmat_dist (1,2)
   write(35,'(A,T25,1pe15.6)')'Cross-section      = ',Tmat_dist (2,2)
   write(35,'(A,T25,1pe15.6)')'Volume             = ',Tmat_dist (3,2)
   write(35,'(A,T25,1pe15.6)')'Effective Radius   = ',Tmat_dist (4,2)
   write(35,'(A,T25,1pe15.6)')'Effective Variance = ',Tmat_dist (5,2)

   write(35,'(/a/)')' PSD output, Bimodal fraction  ------------- '
   write(35,'(A,T25,f10.5)')  'Bi-modal fraction  = ',fraction

!  write Bulk property output
!    ASYMM will be Zero is Do_expcoeffs not set.

   write(35,'(/a/)')' Total Bulk property output ------------- '
   write(35,'(A,T40,1pe20.11)')'Extinction coefficient = ',Tmat_bulk (1)
   write(35,'(A,T40,1pe20.11)')'Scattering coefficient = ',Tmat_bulk (2)
   write(35,'(A,T40,1pe20.11)')'Singlescattering albedo= ',Tmat_bulk (3)
   write(35,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',Tmat_asymm

!  Expansion coefficients
  
   if ( Do_expcoeffs ) then
      write(35,'(/a)')'Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
      write(35,'(a,T30,I5/)')' - Number of coefficients =  ',tmat_ncoeffs
      do L = 1, tmat_ncoeffs
        write(35,'(i5,1p6e20.11)')L-1,(Tmat_expcoeffs (L,k),k=1,6)
      enddo
   endif

!  F-matrix for npna angles (regularly spaced)

   if ( Do_Fmatrix ) then
      write(35,'(/a)')' F-matrix (F11/F22/F33/F44/F12/F34) @ angles (degs)'
      write(35,'(a,I5/)')' - Number of angles =  ',npna
      db = 180.0/dble(npna-1); tb = - db
      do L = 1, npna
        tb = tb + db
        write(35,'(f6.2,6F20.11)')tb,(tmat_Fmatrix (L,k),k=1,6)
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
            write(36,'(/a/)')' LPSD Distribution output, '//CPSD
            write(36,'(A,T25,1pe15.6)')'Number density     = ',LPSD_Tmat_dist (1,m,j) 
            write(36,'(A,T25,1pe15.6)')'Cross-section      = ',LPSD_Tmat_dist (2,m,j)
            write(36,'(A,T25,1pe15.6)')'Volume             = ',LPSD_Tmat_dist (3,m,j)
            write(36,'(A,T25,1pe15.6)')'Effective Radius   = ',LPSD_Tmat_dist (4,m,j)
            write(36,'(A,T25,1pe15.6)')'Effective Variance = ',LPSD_Tmat_dist (5,m,j)
            write(36,'(/a/)')' LPSD Bulk property output, '//CPSD
            write(36,'(A,T40,1pe20.11)')'Extinction coefficient = ',LPSD_Tmat_bulk (1,m,j)
            write(36,'(A,T40,1pe20.11)')'Scattering coefficient = ',LPSD_Tmat_bulk (2,m,j)
            write(36,'(A,T40,1pe20.11)')'Singlescattering albedo= ',LPSD_Tmat_bulk (3,m,j)
            write(36,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',LPSD_Tmat_asymm(m,j)
            if ( Do_expcoeffs ) then
               write(36,'(/a)')'LPSD Expansion coefficient output (a1,a2,a3,a4,b1,b2), '//CPSD
               write(36,'(a,I5/)')' - Number of coefficients =  ',tmat_ncoeffs
               do L = 1, tmat_ncoeffs
                  write(36,'(i5,1p6e20.11)')L-1,(LPSD_Tmat_expcoeffs (L,k,M,J),k=1,6)
               enddo
            endif
            if ( Do_Fmatrix ) then
               write(36,'(/a)')' LPSD F-matrix (F11/F22/F33/F44/F12/F34) @ angles (degs), '//CPSD
               write(36,'(a,I5/)')' - Number of angles =  ',npna
               db = 180.0/dble(npna-1); tb = - db
               do L = 1, npna
                  tb = tb + db
                  write(36,'(f6.2,6F20.11)')tb,(LPSD_tmat_Fmatrix (L,k,m,j),k=1,6)
               enddo
            endif
            close(36)
         enddo
      enddo
   endif

!  write RFE-linearized results
!    Filenames derived from Input name, with '_LRFE_*' added (where * = 1, 2, or 3)

   if ( Do_LinearRef.or.Do_LinearEps ) then
      do J = 1, 2
         write(CJ,'(I1)')J
         MS = 0 ; if (Do_LinearRef) MS = MS + 2 ; if (Do_LinearEps) MS = MS + 1
         do M = 1, MS
            write(C1,'(I1)')M
            CPSD = 'PSD # '//CJ//', parameter # '//C1
            open(37,file=TRIM(ADJUSTL(filename))//'_LRFE_output_PSD#'//CJ//'_PAR#'//C1,status = 'unknown' )
            write(37,'(/a/)')' LRFE Bulk property output, '//CPSD
            write(37,'(A,T40,1pe20.11)')'Extinction coefficient = ',LRFE_Tmat_bulk (1,m,j)
            write(37,'(A,T40,1pe20.11)')'Scattering coefficient = ',LRFE_Tmat_bulk (2,m,j)
            write(37,'(A,T40,1pe20.11)')'Singlescattering albedo= ',LRFE_Tmat_bulk (3,m,j)
            write(37,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',LRFE_Tmat_asymm(m,j)
            if ( Do_expcoeffs ) then
               write(37,'(/a)')'LRFE Expansion coefficient output (a1,a2,a3,a4,b1,b2), '//CPSD
               write(37,'(a,I5/)')' - Number of coefficients =  ',tmat_ncoeffs
               do L = 1, tmat_ncoeffs
                  write(37,'(i5,1p6e20.11)')L-1,(LRFE_Tmat_expcoeffs (L,k,M,j),k=1,6)
               enddo
            endif
            if ( Do_Fmatrix ) then
               write(37,'(/a)')' LRFE F-matrix (F11/F22/F33/F44/F12/F34) @ angles (degs), '//CPSD
               write(37,'(a,I5/)')' - Number of angles =  ',npna
               db = 180.0/dble(npna-1); tb = - db
               do L = 1, npna
                  tb = tb + db
                  write(37,'(f6.2,6F20.11)')tb,(LRFE_tmat_Fmatrix (L,k,m,j),k=1,6)
               enddo
            endif
            close(37)
         enddo
      enddo
   endif

!  write FRC-linearized results (w.r.t. fraction)

   open(38,file=TRIM(ADJUSTL(filename))//'_LFRC',status = 'unknown' )
   write(38,'(/a/)')' LFRC Bulk property output ------------- '
   write(38,'(A,T40,1pe20.11)')'Extinction coefficient = ',LFRC_Tmat_bulk(1)
   write(38,'(A,T40,1pe20.11)')'Scattering coefficient = ',LFRC_Tmat_bulk(2)
   write(38,'(A,T40,1pe20.11)')'Singlescattering albedo= ',LFRC_Tmat_bulk(3)
   write(38,'(A,T40,1pe20.11)')'Asymmetry parameter    = ',LFRC_Tmat_asymm
   if ( Do_expcoeffs ) then
      write(38,'(/a)')'LFRC Expansion coefficient output (a1,a2,a3,a4,b1,b2)'
      write(38,'(a,I5/)')' - Number of coefficients =  ',tmat_ncoeffs
      do L = 1, tmat_ncoeffs
         write(38,'(i5,1p6e20.11)')L-1,(LFRC_Tmat_expcoeffs (L,k),k=1,6)
      enddo
   endif
   if ( Do_Fmatrix ) then
      write(38,'(/a)')' LFRC F-matrix (F11/F22/F33/F44/F12/F34) @ angles (degs)'
      write(38,'(a,I5/)')' - Number of angles =  ',npna
      db = 180.0/dble(npna-1); tb = - db
      do L = 1, npna
         tb = tb + db
         write(38,'(f6.2,6F20.11)')tb,(LFRC_tmat_Fmatrix (L,k),k=1,6)
      enddo
   endif
   close(38)

!  Finish

   return
end subroutine tmat_write_bimodal_extended

!  Finish Module

end module tmat_IO_readwrite

