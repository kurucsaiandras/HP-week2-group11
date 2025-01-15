#!/bin/bash
# 02614 - High-Performance Computing, January 2022
#
#
# Author: Bernd Dammann <bd@cc.dtu.dk>
#
#BSUB -J mm_batch
#BSUB -o mm_batch_%J.out
#BSUB -q hpcintro
#BSUB -n 24
#BSUB -R "rusage[mem=2048]"
#BSUB -R "select[model == XeonE5_2650v4]"
#BSUB -W 25
# uncomment the following line, if you want to assure that your job has
# a whole CPU for itself (shared L3 cache)
#BSUB -R "span[hosts=1] affinity[socket(1)]"

#
EXECUTABLE=poisson_j

GRIDSIZE=1024
MAXITER=10000
TOL=0.01
T=0
OUTPUT=

CC=${1-"gcc"}

# start the collect command with the above settings
./$EXECUTABLE $GRIDSIZE $MAXITER $TOL $T $OUTPUT | grep -v CPU >> test.txt