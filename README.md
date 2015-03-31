# PDE-ODEsolver
A numerical solver for coupled PDE-ODE systems modelling biofilm growth

## TO RUN
```
bash twoDimension_runall.sh FILENAME.txt [{omp, acc}]
```
  
Notes:
  1. `FILENAME.txt` is the parameter file you use. You can copy the original `parameter.txt` and use that. The order of variables needs to be preserved.
  2. `[{omp,acc}]` is an optional argument, it specifies which set of parallel-compatible library you want to use. `omp` is for OpenMP, `acc` is for OpenACC. The default is `omp` if nothing is specified.

## OUTPUTTED
Running the bash script will print the parameters to the screen and then it will print the following after every *nOuts* timesteps:

    time     avgIters     maxIter     avgNit     maxNit     avgM     avgC

Once fortran is done, the script will move all the output files into another folder called 'FILNAMEOutput'. This folder will contain the following:

    |-------------------|-------------------------------------------|
    | 2D_outXX.dat      | The solution at outputted time XX.        |
    |                   | Format: x, M(x), C(x)                     |
    |-------------------|-------------------------------------------|
    | biomassGraphs.pdf | All of the 2D_outXX.dat's in a pdf.       |
    |-------------------|-------------------------------------------|
    | out.dat           | All of 2D_outXX.dat's in a single file    |
    |                   | Format: x, M(x), C(x)                     |
    |-------------------|-------------------------------------------|
    | FILENAME.txt      | The parameter file used in simulation     |
    |-------------------|-------------------------------------------|
    | peakInfo.dat      | Data on the Max M value (x,y) and         |
    |                   | interface location (furthest x where      |
    |                   | M > 0.1)                                  |
    |                   | Frm.: t, x at MaxM(x), Max M(x), interface|
    |-------------------|-------------------------------------------|
    | statReport.dat    | A file reporting computation time and     |
    |                   | avg/max iters for lin solv and between    |
    |                   | solution iterations                       |
    |-------------------|-------------------------------------------|
    | total.dat         | The average M and C at each time....      |
    |                   | Formate: t, avg M(x), avg C(x)            | 
    |-------------------|-------------------------------------------|

## NOTES
Some parameter values make the problem too stiff to solve and computations begin to take a **VERY** long time. 

Also there is not any measures for when the iterations between two solutions *M* and *C* fails to converge. So if the simulator gets stuck and doesn't move, its probably because the parameter set chosen results in a failed convergence (never documented this, so please do if it happens)
