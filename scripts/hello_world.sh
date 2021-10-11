#!/bin/bash

#SBATCH -t 14-00:00:00 # assign time limit

#SBATCH -J trial # assign job name

#SBATCH -N 1 # node number

#SBATCH -o out%tnseq # assign standard output file name

#SBATCH -e err%tnseq # assign standard error file name

#SBATCH -n 1 # assign CPU number (the owner cluster should have 128 threads for 64 cores)

#SBATCH -A mpag-np  # group name (the owner node name)

#SBATCH -p mpag-shared-np  # microbe, plant and animal genetics (the partition name). Note: mpag-np (to have the entire node) or mpag-shared-np (to request a portion of the node)

#SBATCH --mail-type=begin,end # email you when begin, end, fail, requeue

#SBATCH --mail-user=aubrey.hawks@utah.edu # your email address

#SBATCH --mem=5000 # request memory in Gb  (the owner node have maximum 500 GB memory, assign it if you know the exact memory you want to use)

# The following lines are commented out as they work with the notch-peak freecycle
## SBATCH -A karasov-group1 # my group name

## SBATCH --partition notchpeak-freecycle # request notchpeak-freecycle for script testing

sleep 60

echo "Hello World"
