#!/bin/bash
TEMP=$1
TEMP2=`echo $TEMP | cut -d '.' -f 1` # Takes off the '.txt'
DIRNAME=${TEMP2}Output  

# Check if parallel program has been specified
if [ "$#" == 1 ] ; then #Default parallel library
  PARA=omp
elif [ "$2" == "omp" ] || [ "$2" == "acc" ] ; then
  PARA=$2
else
  printf "[!] Invalid Arguments; expected format is:\n[!]  bash ------.sh parameter_file_name.txt [{omp, acc}] \n"
  exit -1
fi

gfortran -O3 -fdefault-real-8 PDE-ODEsolver.f90 ${PARA}cg.f90 ${PARA}amuxd.f90 ${PARA}distdot.f90
echo $1 | ./a.out


# The 2> /dev/null ignores error messages 
# (since it complains if the directory already exists)
mkdir $DIRNAME 2> /dev/null 
mv statReport.dat $DIRNAME/
cp $1 $DIRNAME/

