MODULE GEMSTOOL_numerical_m

!  Everything public, except as Indicated
!   Double precision. All routines stand-alone.

!  Rob Fix 10/26/16. Code added from the NASA depository.
!    -- Several more routines added, for convolution and RSRs

implicit none

!  Precision

integer, parameter :: fpk = SELECTED_REAL_KIND(15)

public  
private trapz, fpk
!public  asym_gauss_f2c, &
!        asym_gauss_f2ci0, &
!        spline, & !updated 11/20/2013
!        splint, & !updated 11/20/2013
!        wtd_mean, &
!        RSR_to_fine_grid, &
!        FlipRealVec, &
!        linear_trawler_ascending,&
!        lintp2,&
!        makechar3


CONTAINS


  SUBROUTINE LINEAR_TRAWLER_ASCENDING &
         ( MAXGRID, NGRID, GRID, X,   &
           F1, F2, K1, K2 )


      INTEGER          :: MAXGRID, NGRID
      real(kind=fpk)        :: GRID(MAXGRID), X
      INTEGER          :: K1, K2
      real(kind=fpk)        :: F1, F2
      LOGICAL          :: LOOP

      IF ( X .LT. GRID(2) ) THEN
        K1 = 1
      ELSE IF ( X .GE. GRID(NGRID-1) ) THEN
        K1 = NGRID - 1
      ELSE
        LOOP = .TRUE.
        K1 = 0
        DO WHILE (LOOP)
          K1 = K1 + 1
          IF ( X .GE. GRID(K1) .AND. X .LT. GRID(K1+1) ) LOOP=.FALSE.
        ENDDO
      ENDIF

      K2 = K1 + 1
      F1 = ( X - GRID(K1) ) / ( GRID(K2) - GRID(K1) )
      F2 = 1.0D0 - F1

      RETURN
  END SUBROUTINE LINEAR_TRAWLER_ASCENDING

!
SUBROUTINE asym_gauss_f2c (fwave, fspec, nf, nspec, hw1e, asym, cwave, cspec, nc)

! =========================================================================
!
! Convolves input spectrum with an asymmetric Gaussian slit function of
! specified HW1E (half-width at 1/e intensity)
!
!  g1=exp(-((lam_iter-lam_ch1(i)).^2)./(hwe(i).*(1+sign(lam_iter-lam_ch1(i)).*asym(i))).^2);
!
! The :asymmetric Gaussian g(x) is defined as
!                   _                                      _
!                  |               (x-x0)^2                 |
!      g(x) =  EXP | - ------------------------------------ |
!                  |_     (hw1e*(1+sign(x-x0)*asym))^2     _|
!
!  FWHM = 2.0 * sqrt(ln(2.0)) * hw1e = 1.66551 * hw1e
! =========================================================================

! convolve high-resolution spectra to low resolution spectra
! fwave: high-resolution wavelength grid
! fspec: high-resolution reference spectra (possibly more than 1)
! nf:    number of wavelengths at high-resolution
! nspec: number of spectra
! fwhm:  slit function in terms of FWHM (nm)
! cwave: low-resolution wavelength grid (could be the same as high-resolution grid)
! cspec: spectra at fwhm
! nc:    number of wavelengths for the low-resolution grid

  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER,                         INTENT (IN)         :: nc, nf, nspec
  REAL (KIND=FPK), DIMENSION (nf), INTENT (IN)         :: fwave
  REAL (KIND=FPK), DIMENSION (nf, nspec), INTENT (IN)  :: fspec
  REAL (KIND=FPK), DIMENSION (nc), INTENT (IN)         :: cwave,hw1e,asym

  REAL (KIND=FPK), DIMENSION (nc, nspec), INTENT (OUT) :: cspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                         :: i, j, midx, sidx, eidx, nhalf, qqq
  REAL (KIND=FPK)                 :: hw1esq, dfw, ssum, intf, cspeci
  REAL (KIND=FPK), DIMENSION (nf) :: slit, sign_asym, asym_factor2, fspeci

  dfw  = fwave(2) - fwave(1)
  nhalf  = INT((sum(hw1e) / MAX(1,size(hw1e))) / dfw * 4.0)

  ! Integration done at the fine wavelength grid (fwave) but slit function centred at the coarse wavelength grid (cwave)

  DO i = 1, nc
     ! Find the closest pixel
     midx = MINVAL(MAXLOC(fwave, MASK=(fwave <= cwave(i)))) + 1

     sidx = MAX(midx - nhalf, 1)
     eidx = MIN(nf, midx + nhalf)

     hw1esq = hw1e(i) ** 2

!mick fix 10/16/2013 - added loop & if block to avoid singularity
     !sign_asym(sidx:eidx) =  (fwave(sidx:eidx) - cwave(i)) / SQRT((fwave(sidx:eidx) - cwave(i))**2)
     DO j=sidx,eidx
       IF ((fwave(j) - cwave(i)) < 1.0E-16_FPK) THEN
         sign_asym(j) = 0.0_FPK
       ELSE
         sign_asym(j) = (fwave(j) - cwave(i)) / SQRT((fwave(j) - cwave(i))**2)
       END IF
       !write(*,*) 'j = ',j,' sign_asym(j) = ',sign_asym(j)
     ENDDO

!  Zero point (code by R. Spurr to avoid trouble)

     qqq = MINVAL(MAXLOC(fwave, MASK=(fwave .eq. cwave(i))))
     if ( qqq.gt.0) sign_asym(MINVAL(MAXLOC(fwave, MASK=(fwave .eq. cwave(i))))) = 1.0

     asym_factor2(sidx:eidx) = (1.0 + sign_asym(sidx:eidx) * asym(i)) ** 2

     slit(sidx:eidx) = EXP(-(cwave(i) - fwave(sidx:eidx))**2 / (hw1esq*asym_factor2(sidx:eidx)) )

     ! Trapezoidal integration
     call trapz(fwave(sidx:eidx),slit(sidx:eidx),eidx-sidx+1,ssum)

     DO j = 1, nspec
         fspeci(sidx:eidx) = fspec(sidx:eidx, j)

         ! Trapezoidal integration
         call trapz(fwave(sidx:eidx), fspeci(sidx:eidx)*slit(sidx:eidx), eidx-sidx+1, intf)

         ! Normalization by slit function area
         if (ssum.gt.0.0) then
             cspeci = intf / ssum
             cspec(i, j) = cspeci
         endif
     ENDDO
  ENDDO

END SUBROUTINE asym_gauss_f2c

! Same as above, except for correcting solar I0 effect using high resolution
! solar reference spectrum following the Beer's Law
! i0: high resolution solar reference
! scalex: scaling, normally number of molecules
SUBROUTINE asym_gauss_f2ci0(fwave, fspec, i0, nf, nspec, scalex, hw1e, asym, cwave, cspec, nc)

  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER,                                INTENT (IN)  :: nc, nf, nspec
  REAL (KIND=FPK),                        INTENT (IN)  :: scalex
  REAL (KIND=FPK), DIMENSION (nf), INTENT (IN)         :: fwave, i0
  REAL (KIND=FPK), DIMENSION (nf, nspec), INTENT (IN)  :: fspec
  REAL (KIND=FPK), DIMENSION (nc), INTENT (IN)         :: cwave, hw1e, asym
  REAL (KIND=FPK), DIMENSION (nc, nspec), INTENT (OUT) :: cspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                             :: i, ntemp
  REAL (KIND=FPK), DIMENSION(nc)        :: ci0
  REAL (KIND=FPK), DIMENSION(nc, nspec) :: cabspec
  REAL (KIND=FPK), DIMENSION(nf, nspec) :: abspec

  ntemp = 1
  CALL asym_gauss_f2c (fwave, i0, nf, ntemp, hw1e, asym, cwave, ci0, nc)

  ! Follow Beer's Law
  DO i = 1, nspec
     abspec(:, i) = i0 * EXP(-fspec(:, i) * scalex)
  ENDDO

  CALL asym_gauss_f2c (fwave, abspec, nf, nspec, hw1e, asym, cwave, cabspec, nc)

  DO i = 1, nspec
     cspec(:, i) = - LOG(cabspec(:, i)/ ci0) / scalex
  ENDDO

END SUBROUTINE asym_gauss_f2ci0


SUBROUTINE trapz (x,y,np,aire)

  ! Integration par la methode des trapezes
  !
  ! Dominique Lefebvre janvier 2007

  ! a = borne inferieure d'integration
  ! b = borne superieure d'integration
  ! n = nombre de pas

  ! aire = surface retournee
  IMPLICIT NONE

  INTEGER,                         INTENT (IN)  :: np
  REAL (KIND=FPK), DIMENSION (np), INTENT (IN)  :: x,y
  REAL (KIND=FPK),                 INTENT (OUT) :: aire

  real(KIND=FPK) a,b,h

  INTEGER i

  ! Boucle d'integration
  ! Initialisation des variables

  aire = 0.0

  ! Boucle de calcul
  ! h = (b-a)/n. Je sais que l'aire de chaque petit trapèze est Ai = (h/2)*(f(a+ih) + f(a+(i-1)h)).

  DO i = 1, np-1
    b = x(i+1)
    a = x(i)
    h = (b-a) 
    if (h.gt.0.0) aire = aire + (h / 2.0D0) * (y(i+1) + y(i))
  ENDDO
END SUBROUTINE trapz


SUBROUTINE spline(nmax,x,y,n,yp1,ypn,y2)
  !input
  INTEGER        :: nmax !max # of x & y points
  INTEGER        :: n !actual # of x & y points
  REAL(kind=fpk) :: x(nmax),y(nmax) !pts of orig func
  REAL(kind=fpk) :: yp1,ypn  !1st derivs of orig func at endpts
  !output
  REAL(kind=fpk) :: y2(nmax) !2nd derivs of the spline interp func
  !local
  INTEGER        :: i,k
  REAL(kind=fpk) :: p,qn,sig,un,u(nmax)

  if (yp1.gt..99e30) then
    y2(1)=0.0d0
    u(1)=0.0d0
  else
    y2(1)=-0.5d0
    u(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  endif

  do i=2,n-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2.
    y2(i)=(sig-1.)/p
    u(i)=(6.0d0*( (y(i+1)-y(i)) / (x(i+1)-x(i)) - &
                  (y(i)-y(i-1)) / (x(i)-x(i-1))   &
                ) / (x(i+1)-x(i-1)) - sig*u(i-1)  &
         )/p
  enddo

  if (ypn.gt..99d30) then
    qn=0.0d0
    un=0.d0
  else
    qn=0.5d0
    un=(3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  endif

  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
  do k=n-1,1,-1
    y2(k)=y2(k)*y2(k+1)+u(k)
  enddo

END SUBROUTINE spline


SUBROUTINE splint(nmax,xa,ya,y2a,n,x,y)
  !input
  INTEGER        :: nmax !max # of x & y points
  INTEGER        :: n !actual # of x & y points
  REAL(kind=fpk) :: xa(nmax),ya(nmax) !pts of orig func
  REAL(kind=fpk) :: y2a(nmax) !2nd derivs of the spline interp func
  REAL(kind=fpk) :: x !x value for which interp y value is desired
  !output
  REAL(kind=fpk) :: y !interp output value
  !local
  INTEGER        :: k,khi,klo
  REAL(kind=fpk) :: a,b,h

  klo=1
  khi=n

1 if (khi-klo.gt.1) then
    k=(khi+klo)/2
    if(xa(k).gt.x)then
      khi=k
    else
      klo=k
    endif
    goto 1
  endif
  h=xa(khi)-xa(klo)
!  if (h.eq.0.0d0) pause
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+ &
    ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0

 END SUBROUTINE splint


 FUNCTION wtd_mean(n,w,x)
   !input
   INTEGER, intent(in) :: n
   REAL(kind=fpk), dimension(:), intent(in) ::  w, x
   !output
   REAL(kind=fpk) :: wtd_mean

   wtd_mean = sum(w(1:n)*x(1:n))/sum(w(1:n))

 END FUNCTION wtd_mean


 SUBROUTINE RSR_to_fine_grid(rwave,rsr,fwave,frsr)
   !input
   REAL(kind=fpk), dimension(:) :: rwave,rsr,fwave
   !output
   REAL(kind=fpk), dimension(:) :: frsr
   !local
   INTEGER        :: i,j,nr,nf,istart
   REAL(kind=fpk) :: slope
   REAL(kind=fpk), dimension(:), allocatable :: rsr_der2

   LOGICAL :: debug=.false.
   !LOGICAL :: debug=.true.

   !define some dimensions
   nr = size(rwave)
   nf = size(fwave)

   if (debug) then
     write(*,*)
     write(*,*) 'in RSR_to_fine_grid:'
     write(*,*) 'nr = ',nr
     write(*,*) 'nf = ',nf
     write(*,*) 'rwave(1)  = ',rwave(1)
     write(*,*) 'rwave(nr) = ',rwave(nr)
     write(*,*) 'fwave(1)  = ',fwave(1)
     write(*,*) 'fwave(nf) = ',fwave(nf)
   end if

   allocate ( rsr_der2(nr) )

   !interpolate RSR
   istart = 1
   do j=1,nf !fine grid loop
     if ( (fwave(j) < rwave(1) ) .or. &
          (fwave(j) > rwave(nr)) ) then
       !fine grid point outside the boundaries of the RSR
       frsr(j) = 0.0_fpk
     else
       i = istart
       do !rsr grid loop
         if ( (rwave(i) <= fwave(j)) .and. &
              (fwave(j) <= rwave(i+1)) ) then

           if (debug) then
             !check - part 1: input
             write(*,*)
             write(*,*) '*************************************'
             write(*,*) 'input:'
             write(*,*) 'For i = ',i,' nr = ',nr
             write(*,*) 'For j = ',j,' nf = ',nf
             write(*,*) 'rwave(i)   fwave(j)   rwave(i+1)'
             write(*,*)  rwave(i) , fwave(j) , rwave(i+1)
             write(*,*) 'rsr(i)                rsr(i+1) '
             write(*,*)  rsr(i)   ,            rsr(i+1)
           end if

           !interpolate RSR to fine grid point - linearly
           !slope = (rsr(i+1) - rsr(i)) / (rwave(i+1) - rwave(i))
           !frsr(j) = slope*(fwave(j) - rwave(i)) + rsr(i)

           !interpolate RSR to fine grid point - using cubic spline
           call spline(nr,rwave,rsr,nr,0.0_fpk,0.0_fpk,rsr_der2)
           call splint(nr,rwave,rsr,rsr_der2,nr,fwave(j),frsr(j))

           if (debug) then
             !check - part 2: output
             write(*,*) 'output:'
             write(*,*) 'rsr(i)    frsr(j)     rsr(i+1) '
             write(*,*)  rsr(i)  , frsr(j)  ,  rsr(i+1)
             !if (mod(nf,5) == 0) read(*,*)
           end if

           exit !exit rsr grid loop
         end if
         i = i + 1
       end do
       !update starting point of rsr loop for the next pass
       istart = i
     end if
   end do

   deallocate ( rsr_der2 )

 END SUBROUTINE RSR_to_fine_grid


 SUBROUTINE FlipRealVec(n,vec)
   !input
   INTEGER, intent(in) :: n
   !input/output
   REAL(kind=fpk), dimension(n), intent(inout) :: vec
   !local
   INTEGER        :: i
   REAL(kind=fpk) :: temp

   do i=1,n/2
     temp       = vec(i)
     vec(i)     = vec(n+1-i)
     vec(n+1-i) = temp
   end do

 END SUBROUTINE FlipRealVec

!

  SUBROUTINE LINTP2 (N1,X1,Y1,N2,X2,Y2)

      INTEGER        ::     N1,N2
      REAL(kind=fpk)      :: X1(N1),Y1(N1),X2(N2),Y2(N2)
      INTEGER        :: INV, N2P, J, I2, I1, J2, I
      REAL(kind=fpk)     :: X, YP, XP, S

      IF(N1.LE.0.OR.N2.LE.0.OR.N1.GT.32767.OR.N2.GT.32767)GOTO 9
!
! --- Checking the condition for inverting arrays ---
!
      INV=0
      IF (N1.EQ.1.OR.N2.EQ.1) GOTO 500
      IF ((X1(N1)-X1(1))*(X2(N2)-X2(1))) 300,300,500
!
! --- Inversion of new grid ---
!

      N2P=0
  300 CONTINUE
      INV=1
      N2P=N2/2
      DO 301 J=1,N2P
      XP=X2(J)
      X2(J)=X2(N2-J+1)
      X2(N2-J+1)=XP
  301 CONTINUE
!
! --- Main block ---
!
  500 IF (N1.EQ.1) GOTO 7
      S=DSIGN(1.0D0,X1(N1)-X1(1))
      I2=1
      I1=2
    1 IF((X2(I2)-X1(1))*S) 2,2,3
    2 Y2(I2)=Y1(1)
      I2=I2+1
      IF(I2.GT.N2) GOTO 999
      GOTO 1
    3 IF((X2(I2)-X1(I1))*S) 4,4,5
    4 X=(X2(I2)-X1(I1-1))/(X1(I1)-X1(I1-1))
      Y2(I2)=Y1(I1)*X+Y1(I1-1)*(1.0D0-X)
      I2=I2+1
      IF(I2.GT.N2) GOTO 999
      GOTO 3
    5 I1=I1+1
      IF(I1.LE.N1) GOTO 3
      DO 6 J2=I2,N2
    6 Y2(J2)=Y1(N1)
      GOTO 999
    7 DO 8 I=1,N2
    8 Y2(I)=Y1(1)
  999 CONTINUE
!
! --- Checking the condition for back inversion ---
!
      IF (INV.NE.1) GOTO 1000
!
! --- Back inversion ---
!
      DO 302 J=1,N2P
      XP=X2(J)
      YP=Y2(J)
      X2(J)=X2(N2-J+1)
      Y2(J)=Y2(N2-J+1)
      X2(N2-J+1)=XP
      Y2(N2-J+1)=YP
  302 CONTINUE
!
! --- Exit block ---
!
 1000 RETURN
    9 PRINT 100,N1,N2
  100 FORMAT(10X,'  Error in subroutine LINTP2 !',2I15)
      STOP '  Error in subroutine LINTP2 !'

  END SUBROUTINE LINTP2

  subroutine makechar3(w,c3)
  integer, intent(in)          :: w
  character(LEN=3),intent(out) :: c3
  c3 = '000'
  if ( w.lt.10 ) then
     write(c3(3:3),'(I1)')W
  else if ( w.gt.99 ) then
     write(c3(1:3),'(I3)')W
  else
     write(c3(2:3),'(I2)')W
  endif
  end subroutine makechar3

!  End module

END MODULE GEMSTOOL_numerical_m
    
