#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --time=06:00:00
#SBATCH --job-name=170070050Project
#SBATCH --error=job.%J.err_node_40
#SBATCH --output=job.%J.out_node_40
#SBATCH --partition=hm

module load compiler/intel/2018.2.199
cd ~/ME766/Path-Planner-ME766/
g++ -o DijkstraOpenMP Dijkstra.cpp -fopenmp

for numThreads in 1 2 4 8 16
do
  sleep 0.1
  export OMP_NUM_THREADS=$numThreads
  rm USA_NY_out.txt
  rm USA_NY_parallel_out.txt
  ./DijkstraOpenMP
done
