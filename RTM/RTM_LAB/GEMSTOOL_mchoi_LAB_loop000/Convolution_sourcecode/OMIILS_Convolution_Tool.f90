Module OMIILS_Convolution_Tool_m

!  Convolution Tool for the OMIILS case
!    -- Programmed 5/4/16

!  Notes 21 september 2016
!    -- The Extract subroutine uses a fake-Gaussian. 
!    -- Need to replace by ILS-database extraction routine, HDF4 etc etc etc.....

!  Contains the following routines (PUBLIC):

!        Extract_OMIILS       Gets 505 OMIILS Slit function Values and wavelengths from data, for each coarse wavelength
!        Spline_OMIILS        Splines  OMIILS Slit function Values to fine-grid wavelengths,  for each coarse wavelength
!        Convolve_OMIILS_1    Applies convolution to a single spectrum
!        Convolve_OMIILS      Applies convolution to multiple spectra
!        Convolve_OMIILS_T    Applies convolution to multiple spectra, with both arrays transposed
!        Convolve_OMIILS_TR   Applies convolution to multiple spectra, with 1 of 2 arrays transposed

   implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  Everything public, except spline interpolation arrays and the parameter FPK

public  
private :: spline, splint, fpk

contains

Subroutine Extract_OMIILS &
         ( E_ndat_c, E_ndat_ILS, ndat_c, ndat_ILS, lambdas_c, &
           OMI_PixelIndex, OMI_Channel, OMI_Window,           &
           ILS_wavs, ILS_Values, fail_ILS, message_ILS )

!  Dimensions

   integer, intent(in) :: E_ndat_c, E_ndat_ILS

!  Numbers

   integer, intent(in)    :: ndat_c
   integer, intent(inout) :: ndat_ILS

!  Wavelengths

   real(fpk), intent(in)  :: lambdas_c  (E_ndat_c)

!  OMI control

   integer, intent(in)    :: OMI_PixelIndex, OMI_Channel, OMI_Window

!  ILS wavelengths and Values

   real(fpk), intent(out) :: ILS_Values  (E_ndat_ILS,E_ndat_c)
   real(fpk), intent(out) :: ILS_Wavs    (E_ndat_ILS,E_ndat_c)

!  Exception handling

   logical      , intent(out) :: fail_ILS
   character*(*), intent(out) :: message_ILS

!  Local

   integer :: nmid, nhalf, w, k, k1, k2
   real(fpk), parameter :: zero = 0.0_fpk
   real(fpk) :: res, sigma, diff, slit, norm

!  Initialize intent(out) quantities

   fail_ILS = .false.
   message_ils = ' '
   ILS_values = zero
   ILS_wavs   = zero

!  Exception handling

   if ( OMI_Channel.gt.2 ) then
      message_ILS = 'Extract OMIILS: OMI Channel not equal to 1 or 2' ; fail_ILS = .true. ; return
   endif
   if ( OMI_Window.gt.3 ) then
      message_ILS = 'Extract OMIILS: OMI window not equal to 1, 2 or 3' ; fail_ILS = .true. ; return
   endif
   if ( OMI_PixelIndex.eq.0 ) then
      message_ILS = 'Extract OMIILS: OMI PixelIndex cannot be 0' ; fail_ILS = .true. ; return
   endif

!  fake it with gaussians

   ndat_ILS = 505 ; nhalf = 252 ; nmid = 253
   res = 0.01 ; sigma = 0.625d0 / sqrt(log(2.0d0)) ; norm = sigma*sqrt(4.0d0 * atan(1.0d0))

   do w = 1, ndat_c

      ILS_values(nmid,w) = 1.0d0
      ILS_wavs(nmid,w)   = lambdas_c(w)

!  Fill up with gaussians

      do k = 1, nhalf
        k1 = nhalf - k + 1 ; k2 = nmid + k
        ILS_wavs(k1,w)   = lambdas_c(w) - res * dble(k)
        ILS_wavs(k2,w)   = lambdas_c(w) + res * dble(k)
        diff = res*dble(k)/sigma ; slit = exp(-diff*diff) 
        ILS_values(k2,w) = slit
        ILS_values(k1,w) = slit
      enddo

!  Normalize (Trapezium norm commented out here)
!      norm = res * ( 0.5d0 * ( ILS_values(1,w) + ILS_values(ndat_ILS,w) ) + sum(ILS_values(2:ndat_ILS-1,w)) )

      do k = 1, ndat_ILS  
         ILS_values(k,w) = ILS_values(k,w) / norm
      enddo

   enddo

!  Finish

   return
End Subroutine Extract_OMIILS

!

Subroutine Spline_OMIILS &
         ( E_ndat, E_ndat_c, E_ndat_ILS, ndat, lambdas, ndat_c,      &
           ndat_ILS, ILS_wavs, ILS_Values, ILS_offsets, ILS_npoints, &
           fail_ILS, message_ILS )

!  Dimensions

   integer, intent(in) :: E_ndat, E_ndat_c, E_ndat_ILS

!  Numbers

   integer, intent(in)    :: ndat, ndat_c, ndat_ILS

!  Wavelengths

   real(fpk), intent(in)  :: lambdas   (E_ndat)

!  ILS wavelengths and Values
!    -- Values will be regridded to the fine wavelengths

   real(fpk), intent(inout) :: ILS_Values  (E_ndat_ILS,E_ndat_c)
   real(fpk), intent(in)    :: ILS_Wavs    (E_ndat_ILS,E_ndat_c)

!  ILS Offsets and Npoints

   integer  , intent(out)   :: ILS_offsets (E_ndat_c)
   integer  , intent(out)   :: ILS_npoints (E_ndat_c)

!  Exception handling

   logical      , intent(out) :: fail_ILS
   character*(*), intent(out) :: message_ILS

!  Local

   real(fpk), parameter :: zero = 0.0_fpk
   real(fpk) :: wav, xfirst, xlast, yt, xt(E_ndat_ILS)
   real(fpk) :: x(ndat_ILS), y(ndat_ILS), y2d(ndat_ILS)
   logical   :: trawl
   integer   :: i, ic, istart, w, j
   
!  Initialize intent(out) quantities
!     --> ILS_Values will be re-assigned

   fail_ILS = .false.
   message_ils = ' '
   ILS_offsets = 0
   ILS_npoints = 0

!  Spline regridding for each coarse wavelength

   istart = 0
   do w = 1, ndat_c

!write(*,*) 'w = ',w,' istart = ',istart

!  Spline part

      xfirst =  ILS_wavs(1,w)
      xlast  =  ILS_wavs(ndat_ILS,w)
      x(1:ndat_ILS) = ILS_wavs(1:ndat_ILS,w)
      y(1:ndat_ILS) = ILS_values(1:ndat_ILS,w)
      call spline(x,y,ndat_ILS,zero,zero,y2d)

!  Find number of points to be regridded
!  mick note 10/11/2016 -
!     i  = index of fine (wvl) pts
!     ic = index of regridded ILS pts
      ic = 0 ; i = istart ; trawl = .true.
      do while(trawl .and. i.lt.ndat)
         i = i + 1 ; wav = lambdas(i)
         !see if fine pt is between 1st and last native ILS pts
         if ( wav.ge.xfirst ) then
            if ( wav.le.xlast ) then
               ic = ic + 1 ; xt(ic) = wav
            else
               trawl = .false.
            endif
         endif
      enddo

!   --> If no points, then exit with failure

      if ( ic .eq. 0 ) then
         fail_ILS = .true.
         message_ILS = 'No points in window - outside of range'
         return
      endif

!  set outputs

      ILS_offsets(w) = i - ic - 1
!mick fix 10/11/2016 - added if check
      if (ILS_offsets(w) < 0) ILS_offsets(w) = 0      
      ILS_npoints(w) = ic

!  Splint regridding

      do j = 1, ic
         call splint(x,y,y2d,ndat_ILS, xt(j), yt )
         ILS_values(j,w) = yt
      enddo

!  Zero remaining entries in reassigned array and upgrade starting point

      ILS_values(ic+1:E_ndat_ILS,w) = zero
      istart = ILS_offsets(w)

!  Done main loop 

   enddo

!  Done

   return
End Subroutine Spline_OMIILS

!

SUBROUTINE Convolve_OMIILS_1 ( E_ndat, E_ndat_c, E_ndat_ILS, ndat_c, &
                               ILS_Values, ILS_offsets, ILS_npoints, &
                               fwav, fspec, cspec )

!  Convolution with Single spectra : Fspec(wfine) ---> Cspec(wcoarse)

!  Dimensions

   integer, intent(in)   :: E_ndat, E_ndat_c, E_ndat_ILS

!  Numbers

   integer, intent(in)   :: ndat_c

!  ILS input

   integer  , intent(in) :: ILS_npoints (E_ndat_c)
   integer  , intent(in) :: ILS_offsets (E_ndat_c)
   real(fpk), intent(in) :: ILS_Values  (E_ndat_ILS,E_ndat_c)

!  Spectra

   real(fpk), intent(in)  :: fspec (E_ndat), fwav(E_ndat)
   real(fpk), intent(out) :: cspec (E_ndat_c)

!  Local

   integer   :: w, i1, i2, j, j1, j2
   real(fpk) :: conv, reshalf

!  Simple sum

   cspec = 0.0d0
   do w = 1, ndat_c
      conv = 0.0_fpk
      do j = 1, ILS_npoints(w) - 1
         j1 = j ; j2 = j1 + 1 ; i1 = ILS_offsets(w) + j1 ; i2 = i1 + 1
         reshalf = 0.5_fpk * ( fwav(i2) - fwav(i1) )
         conv = conv + reshalf * ( ILS_Values(j1,w) * fspec(i1) + ILS_Values(j2,w) * fspec(i2) )
      enddo
      cspec(w) = conv
   enddo

!  Finish

   return
end subroutine Convolve_OMIILS_1

SUBROUTINE Convolve_OMIILS ( E_ndat, E_ndat_c, E_ndat_ILS, ndat_c, nspec, &
                             ILS_Values, ILS_offsets, ILS_npoints,        &
                             fwav, fspec, cspec )

!  Convolution with Regular spectra : Fspec(wfine,S) ---> Cspec(wcoarse,S)

!  Dimensions

   integer, intent(in)   :: E_ndat, E_ndat_c, E_ndat_ILS

!  Numbers

   integer, intent(in)   :: ndat_c, nspec

!  ILS input

   integer  , intent(in) :: ILS_npoints (E_ndat_c)
   integer  , intent(in) :: ILS_offsets (E_ndat_c)
   real(fpk), intent(in) :: ILS_Values  (E_ndat_ILS,E_ndat_c)

!  Spectra

   real(fpk), intent(in)  :: fspec (E_ndat,nspec), fwav(E_ndat)
   real(fpk), intent(out) :: cspec (E_ndat_c,nspec)

!  Local

   integer   :: w, i1, i2, j, j1, j2,  k
   real(fpk) :: conv, reshalf

!  Simple sum

   cspec = 0.0d0
   do w = 1, ndat_c
      do k = 1, nspec
         conv = 0.0_fpk
         do j = 1, ILS_npoints(w) - 1
            j1 = j ; j2 = j1 + 1 ; i1 = ILS_offsets(w) + j1 ; i2 = i1 + 1
            reshalf = 0.5_fpk * ( fwav(i2) - fwav(i1) )
            conv = conv + reshalf * ( ILS_Values(j1,w) * fspec(i1,k) + ILS_Values(j2,w) * fspec(i2,k) )
         enddo
         cspec(w,k) = conv
      enddo
   enddo

!  Finish

   return
end subroutine Convolve_OMIILS

SUBROUTINE Convolve_OMIILS_T ( E_ndat, E_ndat_c, E_ndat_ILS, ndat_c, nspec, &
                               ILS_Values, ILS_offsets, ILS_npoints,        &
                               fwav, fspec, cspec )

!  Convolution with Transpose spectra : Fspec(S,wfine) ---> Cspec(S,wcoarse)

!  Dimensions

   integer, intent(in)   :: E_ndat, E_ndat_c, E_ndat_ILS

!  Numbers

   integer, intent(in)   :: ndat_c, nspec

!  ILS input

   integer  , intent(in) :: ILS_npoints (E_ndat_c)
   integer  , intent(in) :: ILS_offsets (E_ndat_c)
   real(fpk), intent(in) :: ILS_Values  (E_ndat_ILS,E_ndat_c)

!  Spectra

   real(fpk), intent(in)  :: fspec (nspec,E_ndat), fwav(E_ndat)
   real(fpk), intent(out) :: cspec (nspec,E_ndat_c)

!  Local

   integer   :: w, i1, i2, j, j1, j2,  k
   real(fpk) :: conv, reshalf

!  Simple sum

   cspec = 0.0d0
   do w = 1, ndat_c
      do k = 1, nspec
         conv = 0.0d0
         do j = 1, ILS_npoints(w) - 1
            j1 = j ; j2 = j1 + 1 ; i1 = ILS_offsets(w) + j1 ; i2 = i1 + 1
            reshalf = 0.5_fpk * ( fwav(i2) - fwav(i1) )
            conv = conv + reshalf * ( ILS_Values(j1,w) * fspec(k,i1) + ILS_Values(j2,w) * fspec(k,i2) )
         enddo
         cspec(k,w) = conv
      enddo
   enddo

!  Finish

   return
end subroutine Convolve_OMIILS_T

SUBROUTINE Convolve_OMIILS_TR ( E_ndat, E_ndat_c, E_ndat_ILS, ndat_c, nspec, &
                                ILS_Values, ILS_offsets, ILS_npoints,        &
                                fwav, fspec, cspec )

!  Convolution with Reverse-Transpose spectra : Fspec(S,wfine) ---> Cspec(wcoarse,S)

!  Dimensions

   integer, intent(in)   :: E_ndat, E_ndat_c, E_ndat_ILS

!  Numbers

   integer, intent(in)   :: ndat_c, nspec

!  ILS input

   integer  , intent(in) :: ILS_npoints (E_ndat_c)
   integer  , intent(in) :: ILS_offsets (E_ndat_c)
   real(fpk), intent(in) :: ILS_Values  (E_ndat_ILS,E_ndat_c)

!  Spectra

   real(fpk), intent(in)  :: fspec (nspec,E_ndat), fwav(E_ndat)
   real(fpk), intent(out) :: cspec (E_ndat_c,nspec)

!  Local

   integer   :: w, i1, i2, j, j1, j2,  k
   real(fpk) :: conv, reshalf

!  Simple sum

   cspec = 0.0d0
   do w = 1, ndat_c
      do k = 1, nspec
         conv = 0.0d0
         do j = 1, ILS_npoints(w) - 1
            j1 = j ; j2 = j1 + 1 ; i1 = ILS_offsets(w) + j1 ; i2 = i1 + 1
            reshalf = 0.5_fpk * ( fwav(i2) - fwav(i1) )
            conv = conv + reshalf * ( ILS_Values(j1,w) * fspec(k,i1) + ILS_Values(j2,w) * fspec(k,i2) )
         enddo
         cspec(w,k) = conv
      enddo
   enddo

!  Finish

   return
end subroutine Convolve_OMIILS_TR

!******************************************************************************
  SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER            ::  n
      REAL(kind=fpk)       :: yp1,ypn,x(n),y(n),y2(n)
      INTEGER            :: i,k
      REAL(kind=fpk)       :: p,qn,sig,un,u(n)


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
                     (y(i)-y(i-1)) / (x(i)-x(i-1))    &
                   ) / (x(i+1)-x(i-1)) - sig*u(i-1)  )/p
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

      return
  END SUBROUTINE SPLINE

!******************************************************************************

  SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER            ::  n
      REAL(kind=fpk)       :: x, y, xa(n),ya(n),y2a(n)
      INTEGER            ::  k,khi,klo
      REAL(kind=fpk)       :: a,b,h

      klo=1
      khi=n

 1    if (khi-klo.gt.1) then
          k=(khi+klo)/2
          if(xa(k).gt.x)then
            khi=k
          else
            klo=k
          endif
       goto 1
       endif
       h=xa(khi)-xa(klo)
       if (h.eq.0.0d0) stop 'GRONK GRONK GRONK SPLINT'
       a=(xa(khi)-x)/h
       b=(x-xa(klo))/h
       y=a*ya(klo)+b*ya(khi)+ &
        ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0

      return
  END SUBROUTINE SPLINT

!  End module

end Module OMIILS_Convolution_Tool_m

