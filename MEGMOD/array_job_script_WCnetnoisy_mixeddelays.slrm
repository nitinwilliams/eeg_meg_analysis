#!/bin/bash -l
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-200
#SBATCH -o slurm-%A_%a.out

module load matlab

srun matlab_multithread -nodisplay -nosplash -r "WCnet_mixeddelays_noisy_posteriorpredictivechecks($SLURM_ARRAY_TASK_ID) ; exit(0)"
