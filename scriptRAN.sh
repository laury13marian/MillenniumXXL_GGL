#!/bin/bashrc

# USE THE CURRENT WORKING DIRECTORY
#$ -cwd

# NAME OF JOB
#$ -N RAN_new_bins_zmax_100_N1e+6

# MAIL ADDRESS FOR COMMUNICATIONS
#S -M lmarian@mpa-garching.mpg.de

# SEND MAIL AT BEGINNING AND END OF JOB
#S -m be

# REDIRECT OUTPUT TO THE FOLLOWING STREAM
#$ -o /u/lmarian/LauraTreeCode/XXL_CALCULATIONS/LOGS/RANDOM

# JOIN ERROR STREAM TO THE OUTPUT STREAM
#$ -j yes

# RUN AN MPI JOB
#$ -pe impi_hydra 32

# SET THE TIME LIMIT [HRS:MINS:SECS]
#$ -l h_rt=24:00:00

# REQUEST RESOURCES
#$ -R y

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# COMPUTE THE FIDUCIAL POWER SPECTRA

module load impi

mpiexec -np ${NSLOTS} ./exeRAN 1000000 100.
