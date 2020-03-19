module Create_TriAx_for_GEMSTOOL_Mk2_m

   Use SMCConv_Interp_m , only : SMC_GAULEG, SMC_spline, SMC_splint
   Use SMCConv_Develop_m, only : SMCConv_Develop_Mk2
   Use SMCConv_Expand_m

   public  :: Create_TriAx_for_GEMSTOOL_Mk2
   private :: fpk

contains

subroutine Create_TriAx_for_GEMSTOOL_Mk2 &
    ( Max_Coeffs, Cutoff, RefWvl, Quadrature,     & ! Control
      TriAxData_Scatangles, TriAxData_MuellerMat, & ! Input data
      TriAxData_XsecsAsym,  TriAxData_Distpars,   & ! Input data
      TriAx_dist, TriAx_bulk, TriAx_asymm, TriAx_ncoeffs, TriAx_expcoeffs ) ! output for GEMSTOOL

  implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Control inputs
!  ==============

!  Number of coefficients

   integer, intent(in) :: Max_Coeffs

!  Cutoff

   real(fpk), intent(in) :: Cutoff

!  Logical reference flag

   logical, intent(in)   :: RefWvl

!  Single size flag and Radius
!   if set, then use just one size-parameter bin
!   logical  , intent(in) :: SingleSize
!   real(fpk), intent(in) :: SingleRadius

!  Quadrature of size parameters

   real(fpk), intent(in)  :: Quadrature(100)

!  TriAx data, extracted and input here
!  ====================================

!  Detailed arrays for use in a Type structure

   real(fpk), intent(in)  :: TriAxData_Scatangles(500)     ! Degrees
   real(fpk), intent(in)  :: TriAxData_MuellerMat(100,6,500)
   real(fpk), intent(in)  :: TriAxData_XsecsAsym(100,4)
   real(fpk), intent(in)  :: TriAxData_Distpars(100,4)

!  Mueller matrix quantities, indices 1-6 (middle dimension)

!  1:ln(P11)
!  2:P22/P11
!  3:P33/P11
!  4:P44/P11
!  5:P12/P11
!  6:P34/P11

!Triax_Efficiency(k,1) = extinction efficiency
!Triax_Efficiency(k,2) = absorption efficiency
!Triax_Efficiency(k,3) = scattering efficiency
!Triax_Efficiency(k,4) = single-scattering albedo

!Triax_XsecsAsym(k,1)  = extinction cross section
!Triax_XsecsAsym(k,2)  = absorption cross section
!Triax_XsecsAsym(k,3)  = absorption cross section
!Triax_XsecsAsym(k,4)  = asymmetry factor of phase function
	
!Triax_Distpars(k,1) = rproj  = the radius of the surface-equivalent spheres
!Triax_Distpars(k,2) = reff   = the radius of the volume-equivalent spheres
!Triax_Distpars(k,3) = parea  = the projected area
!Triax_Distpars(k,4) = volume = the volume of the particle

!  Output is the TriAx stuff as needed for GEMSTOOL
!  ================================================

!  Quadrature parameters

   real(fpk), intent(out) :: TriAx_dist(5)

!  Bulk optical parameters
!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(fpk), intent(out)  :: TriAx_bulk(3)

!  Asymmetry parameter and number of expansion coefficients

   real(fpk), intent(out)  :: TriAx_asymm

!  Number and value of expansion coefficients

   integer, intent(out)    :: TriAx_ncoeffs
   real(fpk), intent(out)  :: TriAx_expcoeffs(6,0:max_coeffs)

!  Local
!  =====

!  Number of angles for quadrature

   integer, parameter  :: Max_OutAngles = 5000
   integer, parameter  :: Max_InAngles  = 500  ! Data value

!  Quadrature

   Integer   :: N_OutAngles
   Real(fpk) :: OutAngles  ( Max_OutAngles )
   Real(fpk) :: OutCosines ( Max_OutAngles )
   Real(fpk) :: OutWeights ( Max_OutAngles )

!  Local Output Fmatrices (Interpolated)

   Real(fpk) :: OutScatMat ( Max_OutAngles, 6 )

!  Coefficients. Initialized locally inside "Develop" subroutine

   REAL(fpk) :: expcoeffs(0:max_OutAngles,6)

!  Parameters

   real(fpk), parameter :: zero = 0.0_fpk
   real(fpk), parameter :: one  = 1.0_fpk, mone  = - 1.0_fpk
   real(fpk), parameter :: pie2 = 2.0_fpk * acos(mone), dtr = pie2 / 360.0_fpk
   integer  , parameter :: nstokes = 4

!  Help variables

   INTEGER            :: m, q, qd, k, k1, k2, L, N_InAngles, ncoeffs
   REAL    (KIND=fpk) :: InCosines(Max_InAngles)
   REAL    (KIND=fpk) :: InAngles (Max_InAngles)
   REAL    (KIND=fpk) :: Phys_InScatMat(Max_InAngles,6)
   REAL    (KIND=fpk) :: Local_InScatMat(Max_InAngles,6)
   REAL    (KIND=fpk) :: Check_InScatMat(Max_InAngles,6)
   REAL    (KIND=fpk) :: maxd,x1,x2,w1,w2,x,y,y2(Max_InAngles,6),yp1,ypn,percheck(6)

!  Debug

   logical, parameter :: check_expansion = .false.
!   logical, parameter :: check_expansion = .true.

!  Initial section
!  ===============

!  initialize output

   TriAx_ncoeffs   = 0
   TriAx_expcoeffs = zero
   TriAx_bulk      = zero
   TriAx_asymm     = zero
   TriAx_dist      = zero

!  Quadrature - set number of angles

   N_OutAngles = Max_OutAngles

!  ORIGINAL SPECIFICATION ==>
!  select particle size bin and number of coefficients from input
!   if ( SingleSize ) then
!      xchoice = pie2 * SingleRadius / Wavelength
!      mindist = + 1.0d20
!      do n = 1, 100
!         dist = abs ( xchoice - TriAxData_Sizepars(n) )
!         mindist = min(mindist,dist)
!         if (mindist.eq. dist ) nc = n
!      enddo
!   endif
!   TriAx_ncoeffs = N_coeffs

!  set bulk and distpars

   TriAx_bulk(1) = Dot_Product(TriAxData_XsecsAsym(1:100,1),Quadrature(1:100))
   TriAx_bulk(2) = Dot_Product(TriAxData_XsecsAsym(1:100,3),Quadrature(1:100))
   TriAx_bulk(3) = TriAx_bulk(2) / TriAx_bulk(1)
   TriAx_asymm   = Dot_Product(TriAxData_XsecsAsym(1:100,4),Quadrature(1:100))
   do k1 = 1, 4
      TriAx_dist(k1)   = Dot_Product(TriAxData_Distpars(1:100,k1),Quadrature(1:100))
   enddo

!  If only reference wavelength, do not need to find coefficients

   if ( RefWvl ) return

!  Develop Quadrature

   Call SMC_GAULEG ( mone, one, OutCosines, OutWeights, N_OutAngles, Max_OutAngles )
   do k1 = 1, N_OutAngles / 2
      k2 = N_OutAngles + 1 - k1 
      x1 = OutCosines(k1) ; x2 = OutCosines(k2)
      w1 = OutWeights(k1) ; w2 = OutWeights(k2)
      OutCosines(k1) = x2 ; OutCosines(k2) = x1
      OutWeights(k1) = w2 ; OutWeights(k2) = w1
   enddo
   do k = 1, N_OutAngles
      OutAngles(k) = acos ( OutCosines(k) ) / dtr
!write(99,*)k,OutAngles(k) 
   enddo

!  Cosines of input angles.
!    -- Reverse directions for the Splining (Monotonically increasing)

   N_InAngles = 500
   do k = 1, N_InAngles
      k1 = N_InAngles + 1 - k
      InAngles(k1)            = TriAxData_Scatangles(k)
      InCosines(k1)           = cos ( TriAxData_Scatangles(k) * dtr )
   enddo

!  Find max value

   maxd = maxval(Quadrature(1:100))
   do q = 1, 100
      if ( maxd.eq.Quadrature(q)) qd = q
   enddo

!  start Loop over all Size parameters

   do q = 1, 100

      if ( q.gt.qd.and. Quadrature(q) .lt. 1.0d-6 ) go to 6778

!    -- Create local scatMat quantities from input data

      do k = 1, N_InAngles
         k1 = N_InAngles + 1 - k
         Local_Inscatmat(k1,1:6) = TriAxData_MuellerMat(q,1:6,k)
!         Local_Inscatmat(k1,1) = exp ( TriAxData_MuellerMat(q,1,k) )
!         Local_Inscatmat(k1,2:6) = TriAxData_MuellerMat(q,2:6,k)
         Phys_Inscatmat(k1,1)   = exp ( TriAxData_MuellerMat(q,1,k) )
         Phys_Inscatmat(k1,2:6) = TriAxData_MuellerMat(q,2:6,k) * Phys_Inscatmat(k1,1)
      enddo

!  Develop Splines
!    - Set the End-point gradient (YPN) = input gradient at forward-peak
!    - This make a HUGE difference to the accuracy

       do m = 1, 6
          yp1 = zero ; ypn = zero
          ypn = ( Local_InScatMat(N_InAngles,m) -  Local_InScatMat(N_InAngles-1,m) ) / &
                 ( InCosines(N_InAngles) -  InCosines(N_InAngles-1) )
          Call SMC_SPLINE(Max_InAngles,InCosines,Local_InScatMat(:,m),N_InAngles,yp1,ypn,y2(:,m))
       enddo

!  interpolate

      do k = 1, N_OutAngles
         k1 = N_OutAngles + 1 - k 
         x = OutCosines(k1)
         do m = 1, 6
            Call SMC_SPLINT(Max_InAngles,InCosines,Local_InScatMat(:,m),y2(:,m),N_InAngles,x,y)
            if ( m.eq.1 ) OutScatMat(k1,m) = exp(y)
!            if ( m.eq.1 ) OutScatMat(k1,m) = y
            if ( m.gt.1 ) OutScatMat(k1,m) = y * OutScatMat(k1,1)
         enddo
      enddo

!  Get the coefficients from the Develop routine.
!     Output = Expcoeffs,  which is local. Use CUTOFF this time......

      Call SMCConv_Develop_Mk2 &
         ( max_OutAngles, N_OutAngles, Cutoff, &
           OutCosines, OutWeights, OutScatMat, Expcoeffs, Ncoeffs )

!  save coefficients

      do L = 0, ncoeffs
         do m = 1, 6
            TriAx_expcoeffs(m,L) = TriAx_expcoeffs(m,L) + TriAxData_XsecsAsym(q,3) * Quadrature(q) * Expcoeffs(L,m)
         enddo
      enddo
      TriAx_Ncoeffs = max(TriAx_Ncoeffs,Ncoeffs)

!  progress

!      if (mod(q,10).eq.0)write(*,*)&
      if (mod(q,1).eq.0)write(*,*)&
           'Done Quadrature # , # coeffs, last coeff ',q,Quadrature(q),ncoeffs,Expcoeffs(ncoeffs,1)

!  debug 3/6/17.
!      write(4445,*)'Doing Quadrature # ',q,Quadrature(q),ncoeffs, expcoeffs(ncoeffs,1)
!      do L = 0, ncoeffs
!         write(4445,*)L,expcoeffs(L,1:6)
!      enddo

!  Check workings by calling Expand routine, just for one layer
!    - this should give results very close to the original Fmatrix input.  

      if ( check_expansion ) then
         Call SMCConv_Expand &
            ( max_InAngles, max_OutAngles, n_InAngles, Ncoeffs, nstokes, InCosines, expcoeffs, Check_InScatMat )
         do k = N_InAngles, 1, -1
            k1 = N_InAngles + 1 - k 
          write(771,55)k,k1,InAngles(k),InCosines(k), phys_Inscatmat(k,1:6)
          write(772,55)k,k1,InAngles(k),InCosines(k), Check_Inscatmat(k1,1:6)
          do m = 1, 6
             percheck(m) = zero
             if (Check_Inscatmat(k1,m).ne.zero) percheck(m) = ((phys_Inscatmat(k,m)/Check_Inscatmat(k1,m))-one)*100.0d0
          enddo
          write(773,56)k,k1,InAngles(k),InCosines(k), percheck(1:6)
       enddo
56     format(2i5,2f10.4,6f10.5)
55     format(2i5,2f10.4,1p6e15.7)
      endif

      if ( check_expansion ) then
         do k = 1, 500
            write(69,*)k, TriAxData_Scatangles(501-k),Log(phys_Inscatmat(k,1)), log(Check_Inscatmat(501-k,1)), &
             100.0d0 * ( ( phys_Inscatmat(k,1)/Check_Inscatmat(501-k,1) ) - 1.0d0 )
         enddo
      endif

!  continuation point for skipping

6778  continue

!  finish size parameter loop

   enddo

!  set the Gemstool output

   do L = 0, TriAx_ncoeffs
      do m = 1, 6
         TriAx_expcoeffs(m,L) = TriAx_expcoeffs(m,L) / TriAx_bulk(2)
      enddo
   enddo
   TriAx_expcoeffs(1,0) = 1.0d0

!  done

   return
end subroutine Create_TriAx_for_GEMSTOOL_Mk2

!  Finish module

end module Create_TriAx_for_GEMSTOOL_Mk2_m
