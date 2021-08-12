#!/bin/bashrc

# USE THE CURRENT WORKING DIRECTORY
#$ -cwd

# NAME OF JOB
#$ -N JK_GG_RG
# -N GG_RG_P_NRAN_1e+6_zmax_100_part1

# MAIL ADDRESS FOR COMMUNICATIONS
#S -M lmarian@mpa-garching.mpg.de

# SEND MAIL AT BEGINNING AND END OF JOB
#S -m be

# REDIRECT OUTPUT TO THE FOLLOWING STREAM
#$ -o /u/lmarian/LauraTreeCode/XXL_CALCULATIONS/LOGS/JK/GG_RG/

# JOIN ERROR STREAM TO THE OUTPUT STREAM
#$ -j yes

# RUN AN MPI JOB
#$ -pe impi_hydra 32

# SET THE TIME LIMIT [HRS:MINS:SECS]
#$ -l h_rt=24:00:00

# REQUEST RESOURCES
#$ -R y

# RUN AN ARRAY JOB
#$ -t 11-20

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# COMPUTE THE FIDUCIAL POWER SPECTRA

module load impi

module load intel

#mpiexec -np ${NSLOTS} ./exe_GG_RG_p11 54 ${SGE_TASK_ID} 1000000 100

mpiexec -np ${NSLOTS} ./exeJK_GG 101 ${SGE_TASK_ID}

