program gausscheck
   implicit none
   integer, parameter :: fpk = selected_real_kind(15)
   integer    :: k, w
   REAL(FPK) :: FL_DataGAUSSIANS(3,2)
   real(fpk) :: FsSum, FsSum2, lin, grad, ampli, lamda, sigma, var, arg, Gauss, FL_Wavelength

      data FL_DataGAUSSIANS(1,1) / 1.445d0 /
      data FL_DataGAUSSIANS(2,1) / 736.8d0 /
      data FL_DataGAUSSIANS(3,1) / 21.2d0  /
      data FL_DataGAUSSIANS(1,2) / 0.868d0 /
      data FL_DataGAUSSIANS(2,2) / 685.2d0 /
      data FL_DataGAUSSIANS(3,2) / 9.55d0  /


!  Straightline

         FL_Wavelength = 755.0d0
         FsSum2 = 0.0_fpk
         do k = 1, 2
            ampli = FL_DataGaussians(1,k)
            lamda = FL_DataGaussians(2,k)
            sigma = FL_DataGaussians(3,k)
            var = 0.5d0/sigma/sigma
            arg = ( FL_Wavelength - lamda ) * ( FL_Wavelength - lamda ) * var
            Gauss = 0.0_fpk
            if ( arg.lt.88.0d0 ) gauss = ampli * exp ( - arg )
            FsSum2 = FsSum2 + Gauss
         enddo
         FL_Wavelength = 775.0d0
         FsSum = 0.0_fpk
         do k = 1, 2
            ampli = FL_DataGaussians(1,k)
            lamda = FL_DataGaussians(2,k)
            sigma = FL_DataGaussians(3,k)
            var = 0.5d0/sigma/sigma
            arg = ( FL_Wavelength - lamda ) * ( FL_Wavelength - lamda ) * var
            Gauss = 0.0_fpk
            if ( arg.lt.88.0d0 ) gauss = ampli * exp ( - arg )
            FsSum = FsSum + Gauss
         enddo

         grad = ( FsSum - FsSum2 ) / (1.0d0 - ( 775.0_fpk / 755.0_fpk ) )
write(*,*)grad, 1.0d0/grad
!  gaussian

      do w = 640, 840

         FL_Wavelength = real(w,fpk)
         FsSum = 0.0_fpk
         do k = 1, 2
            ampli = FL_DataGaussians(1,k)
            lamda = FL_DataGaussians(2,k)
            sigma = FL_DataGaussians(3,k)
            var = 0.5d0/sigma/sigma
            arg = ( FL_Wavelength - lamda ) * ( FL_Wavelength - lamda ) * var
            Gauss = 0.0_fpk
            if ( arg.lt.88.0d0 ) gauss = ampli * exp ( - arg )
            FsSum = FsSum + Gauss
         enddo
         Lin = FsSum2 + grad * ( 1.0d0 - FL_Wavelength/755.0d0 )
         write(777,'(f10.5,3f12.6)')FL_Wavelength, FsSum, Lin, 100.0d0*(1.0-(Lin/FsSum))

     enddo

   stop
end program gausscheck
