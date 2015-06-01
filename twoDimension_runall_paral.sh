#!/bin/bash
TEMP=$1
TEMP2=`echo $TEMP | cut -d '.' -f 1` # Takes off the '.txt'
DIRNAME=${TEMP2}ParalOutput  

gfortran -O3 -fopenmp -fdefault-real-8 PDE-ODEsolver.f90 CGdiag.f90 init_random_seed.f90
echo $1 | ./a.out

./seperateOutput2D.py

# The 2> /dev/null ignores error messages 
# (since it complains if the directory already exists)
mkdir $DIRNAME 2> /dev/null 
mv peakInfo.dat $DIRNAME/
mv total.dat $DIRNAME/
mv statReport.dat $DIRNAME/
cp $1 $DIRNAME/
mv 2D_output.dat ./$DIRNAME/out.dat
mv COprod.dat ./$DIRNAME/

# Counts files, lines, and words from input (the pipe "|")
find 2D_out* | wc > temp    

#saves files lines and words from wc
read fileName fileCount wordCount < temp  

COUNTER=1
while [ $COUNTER -lt $(($fileCount+1)) ]    # while ( count < filescount)
do
  n=$((1000+$COUNTER));  n=${n#1}    #adds leading zero by removing the 1 from 1XX 
  
  #Opens gnuplot and << EOF write the commands to gnuplot 
  gnuplot << EOF 2> /dev/null
  set terminal post enhanced color eps
  set xr[0:1]
  set yr[0:1.1]
  set xlabel "x"
  set output "biot$n.eps"
  plot "./2D_out$COUNTER.dat" using 1:2 with lines, "./2D_out$COUNTER.dat" using 1:3 with lines
EOF
  mv 2D_out$COUNTER.dat $DIRNAME/
  mv biot$n.eps $DIRNAME/
  let COUNTER=COUNTER+1   #increment counter
done

#Combines all .eps files into a .pdf
gs -sDEVICE=pdfwrite -dDEVICEWIDTHPOINTS=400 -dDEVICEHEIGHTPOINTS=300 -dFIXEDMEDIAswitch -dNOPAUSE -dBATCH -dSAFER -sOutputFile=$DIRNAME/biomassGraphs.pdf $DIRNAME/biot* > /dev/null
rm temp
rm $DIRNAME/biot*
