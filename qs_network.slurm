#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 1:00:00
#SBATCH -J f19
#SBATCH --mail-type=ALL
#SBATCH --output="slurm-%j.out"
#SBATCH --mem=100G

echo '-------------------------'
cd ${SLURM_SUBMIT_DIR}
echo ${SLURM_SUBMIT_DIR}
echo Running on host ${hostname}
echo Time is ${date}
echo SLURM_NODES are $(echo ${SLURM_NODELIST})
echo '-------------------------'
echo -e '\n\n'

export PROCS=${SLURM_CPUS_ON_NODE}

python qs_network.py