#!/bin/bash

# Sample commands:
# lrrs_run_subset.bash gfortran EB
# lrrs_run_subset.bash ifort CR
# lrrs_run_subset.bash gfortran   (test 10 or 11 --> i.e. EB or CR designation only required for tests 1-9)


# Select desired LIDORT-RRS test(s) to run
# Note to user: just set to 1 to activate desired test(s),
#               then execute the script

# 2/9/17. Added Test 11 (BRDF FULL RAMAN)

test[1]=0  #lrrs_o3mono_M_tester # not_working
test[2]=0  #lrrs_LPmono_M_tester # not_working
test[3]=0  #lrrs_LCSmono_M_tester # not_working

test[4]=1  #lrrs_o3bin_M_tester
test[5]=0  #lrrs_LPbin_M_tester
test[6]=0  #lrrs_LCSbin_M_tester # not_working

test[7]=0  #lrrs_brdf_full_LIDVal_tester
test[8]=0  #lrrs_sleave_full_LIDVal_tester
test[9]=0  #lrrs_brdf_full_Raman_tester

test[10]=0 #lrrs_brdf_self_tester
test[11]=0 #lrrs_sleave_self_tester


#### Run the selected test ####

# Define test names:

testname[1]='lrrs_o3mono_M_tester'
testname[2]='lrrs_LPmono_M_tester'
testname[3]='lrrs_LCSmono_M_tester'

testname[4]='lrrs_o3bin_M_tester'
testname[5]='lrrs_LPbin_M_tester'
testname[6]='lrrs_LCSbin_M_tester'

testname[7]='lrrs_brdf_full_LIDVal_tester'
testname[8]='lrrs_sleave_full_LIDVal_tester'
testname[9]='lrrs_brdf_full_Raman_tester'

testname[10]='lrrs_brdf_self_tester'
testname[11]='lrrs_sleave_self_tester'


# Give feedback with regard to LRRS compute mode chosen

echo '# = '$#
arg2="ZZ"
if [ $# -ge 2 ] ; then
  if [ $2 = "EB"  ] || [ $2 = "CR"  ] ; then
    arg2=$2
  fi
fi
#echo 'arg2 = |'$arg2'|'

if [ $arg2 = "EB"  ] ; then
  echo
  echo '*** Testing LRRS in Energy-Balance Mode ***'
elif [ $arg2 = "CR"  ] ; then
  echo
  echo '*** Testing LRRS in Cabannes-Raman Mode ***'
else
  # Check to see if doing one or more main tests.  If so, exit with error msg.
  do_MAIN_tests=0
  for ((i=1 ; i<=9 ; i++)) ; do
    if [ "${test[i]}" = "1" ] ; then
      do_MAIN_tests=1
    fi
  done
  if [ "${do_MAIN_tests}" = "1" ] ; then
    echo
    echo 'Do not recognize 3rd command line argument for chosen mono, bin, or full test(s)'
    echo '  (should be either EB or CR).'
    echo
    exit 1
  fi
fi

ulimit -c 10   # Maximum size of core file created if a problem occurs (in KB)

# Run chosen LIDORT-RRS test(s)

if [ -e makefile ] ; then
  echo
  #check to see if active makefile is current
  if [ $(diff -q lrrs_test/makefile_REENG makefile | wc -l) != "0" ] ; then
    echo 'makefile has changed - copying new version up to script directory ...'
    cp lrrs_test/makefile_REENG makefile
  else
    echo 'makefile in script directory is up to date ...'
  fi
else
  echo 'copying current test makefile up to script directory ...'
  cp lrrs_test/makefile_REENG makefile
fi

for ((i=1 ; i<=11 ; i++)) ; do
  #echo "test[$i] = "${test[i]}
  if [ "${test[i]}" = "1" ] ; then

    # Replace current "lrrs_pars.f90" file with the one required
    # for the current test if necessary
    if [ $i -ge 1 ] && [ $i -le 3 ] ; then
      # Mono tests
      #echo
      #echo 'mock making clean & replacing original pars file'
      #echo 'replacing original pars file'
      if [ $(diff -q lrrs_def/lrrs_pars.f90_NPTS_234 lrrs_def/lrrs_pars.f90 | wc -l) != "0" ] ; then
        make clean
        cp lrrs_def/lrrs_pars.f90_NPTS_234 lrrs_def/lrrs_pars.f90
      fi
    elif [ $i -ge 4 ] && [ $i -le 6 ] ; then
      # Bin tests
      #echo
      #echo 'mock making clean & replacing original pars file'
      #echo 'replacing original pars file'
      if [ $(diff -q lrrs_def/lrrs_pars.f90_NPTS_100 lrrs_def/lrrs_pars.f90 | wc -l) != "0" ] ; then
        make clean
        cp lrrs_def/lrrs_pars.f90_NPTS_100 lrrs_def/lrrs_pars.f90
      fi
    elif [ $i -eq 7 ] || [ $i -eq 10 ] ; then
      # BRDF full & self tests
      #echo
      #echo 'mock making clean & replacing original pars file'
      #echo 'replacing original pars file'
      if [ $(diff -q lrrs_def/lrrs_pars.f90_NPTS_100 lrrs_def/lrrs_pars.f90 | wc -l) != "0" ] ; then
        make clean
        cp lrrs_def/lrrs_pars.f90_NPTS_100 lrrs_def/lrrs_pars.f90
      fi
    elif [ $i -eq 8 ] || [ $i -eq 11 ] ; then
      # SLEAVE full & self tests
      #echo
      #echo 'mock making clean & replacing original pars file'
      #echo 'replacing original pars file'
      if [ $(diff -q lrrs_def/lrrs_pars.f90_sleave_test lrrs_def/lrrs_pars.f90 | wc -l) != "0" ] ; then
        make clean
        cp lrrs_def/lrrs_pars.f90_sleave_test lrrs_def/lrrs_pars.f90
      fi
    elif [ $i -eq 9 ] ; then
      # BRDF full Raman test, need 234 for Mono calculation
      #echo
      #echo 'mock making clean & replacing original pars file'
      #echo 'replacing original pars file'
      if [ $(diff -q lrrs_def/lrrs_pars.f90_NPTS_234 lrrs_def/lrrs_pars.f90 | wc -l) != "0" ] ; then
        make clean
        cp lrrs_def/lrrs_pars.f90_NPTS_234 lrrs_def/lrrs_pars.f90
      fi
    fi

    # Provide config file for the current test.
    #    Some renaming of supplement tests 2/9/17

    if [ $i -le 6 ] ; then
      # Mono or bin test
      cp lrrs_test/config/${testname[i]}.cfg_$arg2         lrrs_test/${testname[i]}.cfg
      cp lrrs_test/config/lrrs_BRDF_ReadInput_LRRSVal.cfg  lrrs_test/.
    elif [ $i -eq 7 ] ; then
      # BRDF full test, LIDVal
      cp lrrs_test/config/brdffulltest_LIDVal_ReadInput.cfg       lrrs_test/.
      cp lrrs_test/config/brdffulltest_BRDF_LIDVal_ReadInput.cfg  lrrs_test/.
    elif [ $i -eq 8 ] ; then
      # SLEAVE full test, LIDVal
      cp lrrs_test/config/sleavefulltest_LIDVal_ReadInput.cfg         lrrs_test/.
      cp lrrs_test/config/sleavefulltest_SLEAVE_LIDVal_ReadInput.cfg  lrrs_test/.
    elif [ $i -eq 9 ] ; then
      # BRDF full test, Full Raman
      cp lrrs_test/config/brdffulltest_Raman_ReadInput.cfg       lrrs_test/.
      cp lrrs_test/config/brdffulltest_BRDF_Raman_ReadInput.cfg  lrrs_test/.
    elif [ $i -eq 10 ] ; then
      # BRDF self test, LIDVal
      cp lrrs_test/config/brdfselftest_LIDVal.cfg  lrrs_test/.
    elif [ $i -eq 11 ] ; then
      # SLEAVE self test, LIDVal
      cp lrrs_test/config/sleaveselftest_*_LIDVal.cfg  lrrs_test/.
    fi

    echo
    echo "Making ${testname[i]} ..."
    echo
    make ${testname[i]}.exe FC=$1 $3 $4
    if [ ! -e ${testname[i]}.exe ] ; then
      echo
      echo 'Error in making file: '${testname[i]}.exe
      echo 'Check compilation or making of the file'
      echo
      exit 1
    fi
    echo
    echo "Running ${testname[i]} ..."
    echo
    ./${testname[i]}.exe

    # Remove current config file(s) & executable
    #rm lrrs_test/*.cfg ./${testname[i]}.exe

    # Move result files to results dir
    #mv *.resg *.plot *.all results

    echo
  fi
done

#make clean
#rm makefile

echo
echo 'done'
echo

exit 0
