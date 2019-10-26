#!/usr/bin/env bash

#name job, output to slurm file, use partition all, run for 60 minutes and use 4GB of ram
#SBATCH -J 'dynomite'
#SBATCH -o out/dynomitescript_output_%j.out
#SBATCH -p all
#SBATCH -t 15
#SBATCH --mem-per-cpu=16G  
#SBATCH --mail-type=END                     # Event(s) that triggers email notification (BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=<temp@princeton.edu>    # Destination email address


# run the latest matlab
module load matlab/R2018a

#move to the repo
cd "/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Spock/"

#run a paranoid version of matlab that crashes gracefully and records why just incase.
# in addition, xvfb-run also creates a virtual desktop so that you can plot easily
