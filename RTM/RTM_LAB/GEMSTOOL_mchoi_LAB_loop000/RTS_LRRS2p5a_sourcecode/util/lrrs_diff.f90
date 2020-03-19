      PROGRAM LRRS_DIFF

!THIS PROGRAM PERFORMS A SPECIALIZED DIFF OF LIDORT-RRS RESULT FILES FOR THE
!FILES SPECIFIED

!PROGRAMMER: MATT CHRISTI
!DATE LAST MODIFIED: 12/20/2015

!DATA DICTIONARY****************************************************************
!
! VARIABLE  = DESCRIPTION
!
!*******************************************************************************

!INTRINSIC SUBPROGRAMS USED BY LRRS_DIFF****************************************
!     TRIM,ADJUSTL,LEN
!*******************************************************************************

!EXTERNAL SUBPROGRAMS USED BY LRRS_DIFF*****************************************
!     GETARG
!*******************************************************************************

      !USE F90_UNIX_ENV !FOR NAG F90/95

      IMPLICIT NONE
!INPUT VARIABLES
      CHARACTER (LEN=1) :: &
        Dummy
      CHARACTER (LEN=120), DIMENSION(2) :: &
        InputFile
!INPUT/OUTPUT VARIABLES
      CHARACTER (LEN=330), DIMENSION(:,:), ALLOCATABLE :: &
        Line
!OUTPUT VARIABLES
      INTEGER, DIMENSION(:), ALLOCATABLE :: &
        LineNumDiffs
!INTERNAL VARIABLES
      INTEGER :: &
        ExpoNum1,ExpoNum2,&
        I,IARG,IOS,ISTAT,J,K,NumLineDiffs,NumResDiffs,&
        NumSavedLines,OutputUnit,StrLen,&
        MaxExtent,Num1,Num2,ThresholdDiffCheck3b,ThresholdDiffCheck3c
      INTEGER, DIMENSION(2) :: &
        InputUnit,FileSize
      CHARACTER (LEN=80) :: &
        OutputFile

      INTEGER, DIMENSION(:), ALLOCATABLE :: &
        SavedLineNum

      INTEGER :: &
        MaxNumLineDiffs=6,& !not used at present
        MaxNumResDiffs=3,&  !should be >= 2
        MagThreshold=-9
      LOGICAL :: &
        TEST=.FALSE.
        !TEST=.TRUE.

!START PROGRAM
      !WRITE(*,*)
      !WRITE(*,*) '*****************************************************'
      !WRITE(*,*) '*****************************************************'
      !WRITE(*,*) 'STARTING LRRS_DIFF'
      !WRITE(*,*)

!mick fix 11/11/2016 - added checks for validity of max inputs
!                    - added thresholds for two added diff checks
!CHECK VALIDITY OF MAX INPUTS
      IF (MaxNumResDiffs < 2) THEN
        WRITE(*,*)
        WRITE(*,*) 'ERROR: SETTING OF MaxNumResDiffs MUST BE >= 2'
        STOP
      END IF

!DEFINE THRESHOLDS FOR TWO ADDED DIFF CHECKS
      ThresholdDiffCheck3b = 5*10**(MaxNumResDiffs-2)
      ThresholdDiffCheck3c = 10**MaxNumResDiffs

!READ TWO INPUT FILES TO DIFF

!OPEN INPUT FILES IF THEY EXIST AND DETERMINE THEIR FILE SIZES
      DO K=1,2
        IARG = K
        CALL GETARG(IARG,InputFile(IARG))
        IF (TEST) THEN
          IF (IARG == 1) THEN
            !InputFile(IARG) = &
            !  'saved_results/gfortran/results_brdf_full_tester.all'
            InputFile(IARG) = &
              'saved_results/gfortran/results_thermal_lps_tester.all'
          ELSE IF (IARG == 2) THEN
            !InputFile(IARG) = &
            !  'results_brdf_full_tester.all'
            InputFile(IARG) = &
              'results_thermal_lps_tester.all'
          END IF
        END IF

        IF (TEST) THEN
          WRITE(*,'(A,I1,3A)') 'K = ',K,&
            ' TRIM(ADJUSTL(InputFile(K))) = |', &
              TRIM(ADJUSTL(InputFile(K))), '|'
          READ(*,*)
        END IF
        InputUnit(K) = 10+K
        OPEN(UNIT=InputUnit(K),&
             FILE=TRIM(ADJUSTL(InputFile(K))), &
             STATUS='OLD',&
             ACTION='READ',&
             IOSTAT=IOS)

        IF (TEST) WRITE(*,*) 'IOS = ',IOS
        IF (IOS /= 0) THEN
          WRITE(*,*) 'ERROR: FILE ' // TRIM(ADJUSTL(InputFile(K))) // ' DOES NOT EXIST'
          STOP
        END IF

        !DETERMINE FILE SIZE
        FileSize(K) = 0
        DO
          !READ FROM CURRENT LINE
          READ(InputUnit(K),'(A)',IOSTAT=ISTAT) Dummy
          IF (ISTAT == -1) EXIT
          FileSize(K) = FileSize(K) + 1
        END DO
        REWIND (InputUnit(K))
      END DO

!CHECK FILE SIZES FOR COMPATIBILITY.  IF PROBLEM, THEN ABORT.
      IF (FileSize(1) /= FileSize(2)) THEN
        WRITE(*,*) 'FILES SIZES NOT EQUAL FOR FILES:'
        WRITE(*,*) TRIM(ADJUSTL(InputFile(1)))
        WRITE(*,*) TRIM(ADJUSTL(InputFile(2)))
        STOP
      END IF

!OPEN OUTPUT FILE
      OutputFile = 'diff_' // Trim(InputFile(2))  ! bug fix 9/13/17 Rob, add Trim 
      OutputUnit = 13
      OPEN(UNIT=OutputUnit,&
           FILE=TRIM(ADJUSTL(OutputFile)),&
           STATUS='REPLACE',&
           ACTION='WRITE')

!ALLOCATE SOME ARRAYS
      ALLOCATE ( Line(2,FileSize(1)),&
                 LineNumDiffs(0:MaxNumLineDiffs),&
                 SavedLineNum(FileSize(1)) )

!DO A MORE DETAILED DIFF OF LINES FROM THE INPUT FILES THAN THE
!STANDARD UNIX DIFF UTILITY AND THEN SEND THE RESULTS TO THE OUTPUT FILE
      NumSavedLines = 0
      DO I=1,FileSize(1)
        !READ FROM CURRENT LINE IN EACH INPUT FILE
        READ(InputUnit(1),'(A)',IOSTAT=IOS) Line(1,I)
        IF (IOS /= 0) THEN
          WRITE(*,*) 'ERROR: IOS = ',IOS
          READ(*,*)
        END IF

        READ(InputUnit(2),'(A)',IOSTAT=IOS) Line(2,I)
        IF (IOS /= 0) THEN
          WRITE(*,*) 'ERROR: IOS = ',IOS
          READ(*,*)
        END IF

        !DETERMINE IF CURRENT LINE HAS SIGNIFICANT DIFFS (I.E. IT CONSISTS OF ANY
        !RESULTS WHOSE NUMBER OF DIFFS ARE MORE THAN MaxNumResDiffs AND WHOSE
        !MAGNITUDES ARE MORE THAN MagThreshold
        StrLen = LEN(TRIM(Line(1,I)))

        !#1: LINE PRE-CHECK
        NumLineDiffs = 0
        DO J=1,StrLen
          IF (Line(1,I)(J:J) /= Line(2,I)(J:J)) NumLineDiffs = NumLineDiffs + 1
        END DO

        IF (NumLineDiffs > 0) THEN
          DO J=1,StrLen
!mick fix 8/9/2014 - adjust valid result discriminator
            !IF (Line(1,I)(J:J) == 'E') THEN
            IF ( (Line(1,I)(J:J) == 'E') .AND. &
                 ( (Line(1,I)(J+1:J+1) == '+') .OR. &
                   (Line(1,I)(J+1:J+1) == '-') ) &
               ) THEN

              READ(Line(1,I)(J+1:J+3),'(I3)') ExpoNum1
              READ(Line(2,I)(J+1:J+3),'(I3)') ExpoNum2

              !#2: CHECK MAGNITUDE OF CURRENT RESULT
              IF ((ExpoNum1 > MagThreshold) .AND. (ExpoNum2 > MagThreshold)) THEN
                !#3a: CHECK NUMBER OF INDIVIDUAL RESULT DIFFS
                NumResDiffs = 0
                DO K=J-8,J+3
                  IF (Line(1,I)(K:K) /= Line(2,I)(K:K)) NumResDiffs = NumResDiffs + 1
                END DO

!mick fix 11/11/2016 - added checks #3b & #3c for two unusual conditions 
                !CHECK FOR TWO UNUSUAL CONDITIONS
                !WRITE(*,*)
                MaxExtent = MaxNumResDiffs + 1
                READ(Line(1,I)(J-MaxExtent:J-1),'(I4)') Num1
                READ(Line(2,I)(J-MaxExtent:J-1),'(I4)') Num2
                !WRITE(*,*) 'Num1 = ',Num1,' Num2 = ',Num2
                IF (NumResDiffs > MaxNumResDiffs) THEN
                  !#3b: TYPE #1 ERROR (FALSE POSITIVE) - ANOMALOUS CLOSE CONDITION
                  IF (ABS(Num1 - Num2) < ThresholdDiffCheck3b) &
                    NumResDiffs = NumResDiffs - 1
                ELSE IF (NumResDiffs <= MaxNumResDiffs) THEN
                  !#3c: TYPE #2 ERROR (FALSE NEGATIVE) - ANOMALOUS FAR AWAY CONDITION
                  IF (ABS(Num1 - Num2) > ThresholdDiffCheck3c) &
                    NumResDiffs = MaxExtent
                END IF

                !SAVE DIFF INFO FOR OUTPUT
                IF (NumResDiffs > MaxNumResDiffs) THEN
                  NumSavedLines = NumSavedLines + 1
                  SavedLineNum(NumSavedLines) = I
                  EXIT !EXIT SEARCH ON CURRENT PAIR OF LINES
                END IF
              END IF
            END IF
          END DO
        END IF

      END DO

!SEND DIFF INFO TO OUTPUT FILE
      WRITE(OutputUnit,'(A,I4)') 'TOTAL NUMBER OF LINES IN RESULT FILE: ',FileSize(1)

      WRITE(OutputUnit,*)
      WRITE(OutputUnit,'(A)') 'NUMBER OF LINES CONTAINING RESULTS OF SIGNIFICANT MAGNITUDE'
      WRITE(OutputUnit,'(A,I3)') '  WITH SIGNIFICANT DIFFS: ',NumSavedLines
      WRITE(OutputUnit,'(A,I3,A)') '  (MAGNITUDE DEFINED TO BE OF ORDER > 10^',MagThreshold,')'
      WRITE(OutputUnit,'(A,I1,A)') '  (NUMBER OF DIFFS DEFINED TO BE    >    ',MaxNumResDiffs,')'

      IF (NumSavedLines > 0) THEN
        WRITE(OutputUnit,*)
        DO I=1,NumSavedLines
          WRITE(OutputUnit,'(A,I4,A)') 'LINE ',SavedLineNum(I),':'
          WRITE(OutputUnit,'(2A)') '< ',TRIM(ADJUSTL( Line(1,SavedLineNum(I)) ))
          WRITE(OutputUnit,'(2A)') '> ',TRIM(ADJUSTL( Line(2,SavedLineNum(I)) ))
        END DO
      END IF

!CLOSE FILES
      CLOSE(InputUnit(1))
      CLOSE(InputUnit(2))
      CLOSE(OutputUnit)

!DEALLOCATE SOME ARRAYS
      DEALLOCATE ( Line,LineNumDiffs,SavedLineNum )

      IF (TEST) THEN
        WRITE(*,*) 'FINISHED WITH CURRENT CALL TO LRRS_DIFF'
        READ(*,*)
      END IF

!END PROGRAM
      !WRITE(*,*)
      !WRITE(*,*) 'LRRS_DIFF DONE'
      !WRITE(*,*) '*****************************************************'
      !WRITE(*,*) '*****************************************************'

      END PROGRAM LRRS_DIFF

!*******************************************************************************
!*******************************************************************************
