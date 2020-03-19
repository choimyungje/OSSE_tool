#!/bin/bash

# Sample command:
# lrrs_check2.bash dirname
#   where "dirname" is the directory name in subdirectory "saved_results"
#   one wishes to compare their package test results against

# Run LIDORT-RRS checks (using tailored lrrs diff utility)

cd lrrs_test

file_count=0
echo
echo 'checking results of testers...'
for filename in results_*.* ; do
  #echo
  #echo '  filename     = '$filename
  if [ -e $filename ]; then
    #echo '  filename exists'
    ref_filename=saved_results/$1/$filename
    #echo '  ref_filename = '$ref_filename
    if [ -e $ref_filename ]; then
      #echo '  ref_filename exists'
      ../lrrs_diff $ref_filename $filename
      echo $((++file_count)) files processed
    fi
  fi
done

cd ..

echo
echo 'done'
echo

exit 0

