#!/bin/bash

# This script runs an LRRS-reengineered package tests which have been set up to
# run using OpenMP (see http://openmp.org for further information about OpenMP)

# Sample commands:
# lrrs_run_REENG.bash gfortran EB
# lrrs_run_REENG.bash gfortran CR

# Note: gfortran is the only OpenMP 3.0 compatible compiler set up in the makefile
#       at this time


# Select desired LRRS OpenMP test(s) to run

# Note to user: (1) set to 1 to activate the desired test(s)
#               (2) insure enough memory is set aside for each OpenMP thread using
#                   the environment variable "OMP_STACKSIZE" 
#               (3) execute the script

# OMP tests:

test[1]=0 #lrrs_o3bin_M_REENG_tester_OMP  --> use stacksize 500M, see below
test[2]=1 #lrrs_LPbin_M_REENG_tester_OMP  --> use stacksize 500M, see below    STILL TO BE DETERMINED

# Do some computer memory set up

ulimit -s unlimited         # Maximum stack size of the main thread
#ulimit -c unlimited         # Maximum size of core file created if a problem occurs
ulimit -c 5000

# OMP_STACKSIZE -->         # Maximum stack size of OpenMP-spawned threads

# gfortran requirements for:
#export OMP_STACKSIZE=500M    #o3bin_REENG - 2 thread
export OMP_STACKSIZE=500M    #LPbin_REENG - 2 thread

#### Run the selected test(s) ####

# Define test names:

testname[1]='lrrs_o3bin_M_REENG_tester'
testname[2]='lrrs_LPbin_M_REENG_tester'

# make clean before EVERY test

# make clean

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
fi

# Run chosen LRRS tests using OpenMP

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

for ((i=1 ; i<=2 ; i++)) ; do
   #echo "test[$i] = "${test[i]}
   if [ "${test[i]}" = "1" ] ; then

      # Replace current "lrrs_pars.f90" file with the one required
      # for the current bin_REENG test if necessary
      if [ $(diff -q lrrs_def/lrrs_pars.f90_NPTS_100 lrrs_def/lrrs_pars.f90 | wc -l) != "0" ] ; then
        make clean
        cp lrrs_def/lrrs_pars.f90_NPTS_100 lrrs_def/lrrs_pars.f90
      fi

      # Provide config file for the current test      
      configfile=${testname[i]}.cfg
      cp lrrs_test/config/${configfile}_$arg2 lrrs_test/${configfile}
      cp lrrs_test/config/lrrs_BRDF_ReadInput_LRRSVal.cfg lrrs_test/.

      # Do selected test
      exec=${testname[i]}.exe
      echo
      echo "making ${testname[i]} ..."
      echo
      make $exec FC=$1 OPENMP=t $3
      echo
      echo "running ${testname[i]} ..."
      echo
      ./$exec
      echo

      # Remove current config file(s) from test dir
      #rm lrrs_test/*.cfg

      # Move result files to results dir
      #mv *.resg* *.plot* *.all results
   fi
done

#make clean

#rm makefile

echo
echo 'done'
echo

exit 0
