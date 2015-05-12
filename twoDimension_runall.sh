#!/bin/bash
TEMP=$1
TEMP2=`echo $TEMP | cut -d '.' -f 1` # Takes off the '.txt'
DIRNAME=${TEMP2}Output  

gfortran -O3 -fdefault-real-8 PDE-ODEsolver.f90 CGdiag.f90
echo $1 | ./a.out

# The 2> /dev/null ignores error messages 
# (since it complains if the directory already exists)
mkdir $DIRNAME 2> /dev/null 
mv statReport.dat $DIRNAME/
cp $1 $DIRNAME/

