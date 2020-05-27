#!/usr/bin/env bash

#name job, output to slurm file, use partition all, run for 60 minutes and use 4GB of ram
#SBATCH -J 'Wahoo'
#SBATCH -o out/STnmf_output_%j.out
#SBATCH -p all
#SBATCH -t 5
#SBATCH --mem-per-cpu=16G  
#SBATCH --array=1000-2501
#SBATCH --exclude=redshirt-n[12-49]
#SBATCH --mail-type=FAIL                         # Event(s) that triggers email notification (BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=<camdenm@princeton.edu>    # Destination email address

# run the latest matlab
module load matlab/R2018a

