#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --time=06:00:00
#SBATCH --job-name=170070050Project
#SBATCH --error=job.%J.err_node_40
#SBATCH --output=job.%J.out_node_40
#SBATCH --partition=hm

g++ -o FWAlgoOpenMP FloydWarshall.cpp -fopenmp
for N in 10 20 40 80 160 320 640 1280 2560 5120
do
  for iter in 1 2 3 4 5 6 7 8 9 10 11
  do
    for numThreads in 1 2 4 8 10 12 16
    do
	  sleep 0.1
      ./FWAlgoOpenMP $N $numThreads
	done
  done
done
