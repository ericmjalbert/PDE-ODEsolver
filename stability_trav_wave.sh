#!/bin/bash

list=$(echo parameterOutput/out* |\
  sed 's/[^0-9 ]//g' |\
  fmt -1 |\
  sort -n |\
  sed 's#^#parameterOutput/out##g' |\
  sed 's#$#.dat#')

rm wave_stddev.dat 2> /dev/null

for i in $list
do
  cp $i ./stability.dat
  python stability_stddev.py >> wave_stddev.dat
done

mv wave_stddev.dat *Output/
rm stability.dat
