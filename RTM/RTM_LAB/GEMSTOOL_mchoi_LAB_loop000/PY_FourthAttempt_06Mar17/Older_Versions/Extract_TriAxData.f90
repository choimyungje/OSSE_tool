module Extract_TriAxData_m

public
private :: fpk, getblank, getcommas

contains

subroutine Extract_TriAxData ( Datafile, &
      TriAxData_Scatangles, TriAxData_Sizepars,   TriAxData_MuellerMat, & ! Output data
      TriAxData_XsecsAsym,  TriAxData_Efficiency, TriAxData_Distpars   )  ! output data

 implicit none

!  Precision

   integer, parameter :: fpk = SELECTED_REAL_KIND(15)

!  TriAx data, extracted and output here
!  =====================================

!  data file

   character*(*), intent(in) :: Datafile

!  Detailed arrays for use in a Type structure

   real(fpk), intent(out) :: TriAxData_Scatangles(500)     ! Degrees
   real(fpk), intent(out) :: TriAxData_Sizepars  (100)     

   real(fpk), intent(out)  :: TriAxData_MuellerMat(100,6,500)
   real(fpk), intent(out)  :: TriAxData_XsecsAsym(100,4)
   real(fpk), intent(out)  :: TriAxData_Efficiency(100,4)
   real(fpk), intent(out)  :: TriAxData_Distpars(100,4)

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

!  Local
!  =====

!  holding arrays

   real(fpk) :: Array(100,6,505)
   real(fpk) :: Scatangles(500)     ! Degrees
   real(fpk) :: Sizepars  (100)     

!  Local

   character*1   :: blank, c1
   character*100 :: regline, startline
   integer       :: k, kang, kpar, kpos, ms, mf, ncom, bpos, blen, lc
   logical       :: trawl

!  Data statements
!  ===============

Data Scatangles &
/0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,&
0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,&
0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.40,0.41,0.42,0.43,0.44,0.45,0.46,0.47,&
0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60,0.61,0.62,0.63,&
0.64,0.65,0.66,0.67,0.68,0.69,0.70,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,&
0.80,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,&
0.96,0.97,0.98,0.99,1.00,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,1.10,1.11,&
1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,&
1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,&
1.44,1.45,1.46,1.47,1.48,1.49,1.50,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,&
1.60,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.70,1.71,1.72,1.73,1.74,1.75,&
1.76,1.77,1.78,1.79,1.80,1.81,1.82,1.83,1.84,1.85,1.86,1.87,1.88,1.89,1.90,1.91,&
1.92,1.93,1.94,1.95,1.96,1.97,1.98,1.99,2.00,2.05,2.10,2.15,2.20,2.25,2.30,2.35,&
2.40,2.45,2.50,2.55,2.60,2.65,2.70,2.75,2.80,2.85,2.90,2.95,3.00,3.05,3.10,3.15,&
3.20,3.25,3.30,3.35,3.40,3.45,3.50,3.55,3.60,3.65,3.70,3.75,3.80,3.85,3.90,3.95,&
4.00,4.05,4.10,4.15,4.20,4.25,4.30,4.35,4.40,4.45,4.50,4.55,4.60,4.65,4.70,4.75,&
4.80,4.85,4.90,4.95,5.00,5.10,5.20,5.30,5.40,5.50,5.60,5.70,5.80,5.90,6.00,6.10,&
6.20,6.30,6.40,6.50,6.60,6.70,6.80,6.90,7.00,7.10,7.20,7.30,7.40,7.50,7.60,7.70,&
7.80,7.90,8.00,8.10,8.20,8.30,8.40,8.50,8.60,8.70,8.80,8.90,9.00,9.10,9.20,9.30,&
9.40,9.50,9.60,9.70,9.80,9.90,10.00,10.50,11.00,11.50,12.00,12.50,13.00,13.50,14.00,&
14.50,15.00,16.00,17.00,18.00,19.00,20.00,21.00,22.00,23.00,24.00,25.00,26.00,27.00,&
28.00,29.00,30.00,31.00,32.00,33.00,34.00,35.00,36.00,37.00,38.00,39.00,40.00,41.00,&
42.00,43.00,44.00,45.00,46.00,47.00,48.00,49.00,50.00,51.00,52.00,53.00,54.00,55.00,&
56.00,57.00,58.00,59.00,60.00,61.00,62.00,63.00,64.00,65.00,66.00,67.00,68.00,69.00,&
70.00,71.00,72.00,73.00,74.00,75.00,76.00,77.00,78.00,79.00,80.00,81.00,82.00,83.00,&
84.00,85.00,86.00,87.00,88.00,89.00,90.00,91.00,92.00,93.00,94.00,95.00,96.00,97.00,&
98.00,99.00,100.00,101.00,102.00,103.00,104.00,105.00,106.00,107.00,108.00,109.00,110.00,&
111.00,112.00,113.00,114.00,115.00,116.00,117.00,118.00,119.00,120.00,121.00,122.00,123.00,&
124.00,125.00,126.00,127.00,128.00,129.00,130.00,131.00,132.00,133.00,134.00,135.00,136.00,&
137.00,138.00,139.00,140.00,141.00,142.00,143.00,144.00,145.00,146.00,147.00,148.00,149.00,&
150.00,151.00,152.00,153.00,154.00,155.00,156.00,157.00,158.00,159.00,160.00,161.00,162.00,&
163.00,164.00,165.00,166.00,167.00,168.00,169.00,170.00,171.00,172.00,173.00,174.00,175.00,&
175.50,175.75,176.00,176.25,176.50,176.75,177.00,177.25,177.50,177.75,178.00,178.25,178.50,&
178.75,179.00,179.25,179.50,179.75,180.00 /

Data Sizepars &
/0.025000,0.027824,0.030968,0.034466,0.038360,&
0.042694,0.047517,0.052886,0.058861,0.065510,&
0.072912,0.081148,0.090316,0.100519,0.111876,&
0.124515,0.138582,0.154239,0.171663,0.191057,&
0.212642,0.236666,0.263402,0.293160,0.326280,&
0.363142,0.404168,0.449829,0.500649,0.557209,&
0.620160,0.690222,0.768200,0.854988,0.951581,&
1.059085,1.178735,1.311903,1.460115,1.625072,&
1.808665,2.012999,2.240418,2.493530,2.775237,&
3.088769,3.437724,3.826100,4.258355,4.739443,&
5.274882,5.870812,6.534068,7.272255,8.093839,&
9.008241,10.025948,11.158631,12.419279,13.822348,&
15.383930,17.121931,19.056282,21.209169,23.605276,&
26.272085,29.240178,32.543590,36.220205,40.312188,&
44.866461,49.935255,55.576696,61.855480,68.843610,&
76.621225,85.277517,94.911755,105.634421,117.568482,&
130.850795,145.633679,162.086660,180.398419,200.778952,&
223.461977,248.707620,276.805392,308.077514,342.882607,&
381.619809,424.733351,472.717651,526.122983,585.561788,&
651.715698,725.343355,807.289105,898.492684,1000.000000/

!  Zero holding array

   Array = 0.0d0

!  Open file

   open(1,file=Trim(Datafile),status='old')

!  Read header and initial set of blanks

   do k = 1, 157
     read(1,*)
   enddo

   lc = 157

!  Loop over particles
!  -------------------

   do kpar = 1, 100

!  read initial set of blanks

     do k = 1, 21
       read(1,*)
     enddo
     lc = lc + 21

!  Position loop

     do 55 kpos = 1, 6
       read(1,'(a)')startline ; ms = 1 ; lc = lc + 1
       call getcommas(startline,ncom) ; mf = ms + ncom - 2
       backspace(1)
       read(1,*,err=909)blank,array(kpar,kpos,ms:mf) ; ms = mf + 1
       trawl = .true.
       do while (trawl)
         read(1,'(a)')regline ; lc = lc + 1
         call getblank(regline,bpos,blen)
write(77,*)kpar,kpos,lc,Trim(regline),bpos,blen
    if ( lc.eq.27306 ) stop'now'
         
         if ( kpar.eq.100.and.kpos.eq.6.and.blen.eq.1) goto 677
         if ( bpos.eq.5 ) then
           trawl = .false.
write(*,*)'hello3',kpar,kpos,bpos
if ( lc.eq.27305 ) stop'gggggggggggggggggggggg, bpos 5'
         else if ( bpos.gt.15 ) then
           call getcommas(regline,ncom) ; mf = ms + ncom - 5
           backspace(1)
if ( lc.eq.27305 ) pause'gggggggggggggggggggggg, bpos >15'
           read(1,*,err=910)array(kpar,kpos,ms:mf),c1,c1,c1,c1
           trawl = .false.
write(*,*)'hello1',kpar,kpos,bpos
         else if ( bpos.eq.3 ) then
           backspace(1) ; trawl = .false. ; lc = lc - 1
write(*,*)'hello2',kpar, kpos,bpos
         else if ( bpos.eq.0 ) then
           call getcommas(regline,ncom) ; mf = ms + ncom - 1
if ( lc.eq.27305 ) stop'gggggggggggggggggggggg, bpos 0'
           backspace(1)
           read(1,*,err=911)array(kpar,kpos,ms:mf) ; ms = mf + 1
         endif

       enddo

! if ( kpos.eq.1 ) then
!     do k = 1, 505
!        write(566,*)k,array(kpar,kpos,k)
!     enddo
! endif

55   continue

!  finish

   enddo

!  Normal finish

677 close(1)

!  Copy to the TriAx Arrays

   do k = 1, 100
      TriAxData_DistPars  (k,1:4) = Array(k,3,501:504)
      TriAxData_XsecsAsym (k,1:4) = Array(k,2,501:504)
      TriAxData_Efficiency(k,1:4) = Array(k,1,501:504)
      do kang = 1, 500
         TriAxData_MuellerMat(k,1:6,kang) = Array(k,1:6,kang)
      enddo
   enddo
   TriAxData_Scatangles = ScatAngles
   TriAxData_Sizepars = Sizepars  
   return

!  debug error finishes

909 write(*,*) 'Error',lc ; return
910 write(*,*) 'Error',lc ; return
911 write(*,*) 'Error',lc

!  done

return
end subroutine Extract_TriAxData

subroutine getcommas(line,ncom)
character*(*), intent(in)  :: line
integer      , intent(out) :: ncom
ncom = 0
do n = 1, len(line)
  if ( line(n:n).eq.',' )ncom=ncom+1
enddo
return
end subroutine getcommas

subroutine getblank(line,bpos,blen)
character*(*), intent(in)  :: line
integer      , intent(out) :: bpos,blen
bpos = 0 ; blen=len(trim(line))
do n = 1, blen
  if ( line(n:n).eq.'_' )bpos = n
  if ( bpos.ne.0 ) go to 56
enddo
56 continue
return
end subroutine getblank

!  End module

end module Extract_TriAxData_m
