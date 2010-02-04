#!/bin/bash
#PBS -l nodes=2:ppn=8
#PBS -q default
#PBS -N RNAflow
#PBS -r n

START=$(date +%s)

qsub --env OMP_NUM_THREADS=3 -A ParallelSampling -n 512 --mode smp -q default -t 30 flow_ga.x

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Total time of job: $DIFF seconds"
