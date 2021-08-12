#!/bin/tcsh

echo READ DATA 

set NSUB = 216

set iSC = 1

module load impi

while ( $iSC <= $NSUB )

mpiexec -np 1 ./exeU2 54 $iSC 1000 10 1000000

@ iSC++

end

