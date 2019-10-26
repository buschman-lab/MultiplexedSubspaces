#!/usr/bin/env bash

#name job, output to slurm file, use partition all, run for 60 minutes and use 4GB of ram
#SBATCH -J 'Swiffer'
#SBATCH -o out/Swiffer_output_%j.out
#SBATCH -p all
#SBATCH -t 720
#SBATCH --mem-per-cpu=16G  
#SBATCH --array=32
#SBATCH --mail-type=FAIL                         # Event(s) that triggers email notification (BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=<camdenm@princeton.edu>    # Destination email address

# run the latest matlab
module load matlab/R2018a

#move to the repo
cd "/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Spock/"

#run a paranoid version of matlab that crashes gracefully and records why just incase.
# in addition, xvfb-run also creates a virtual desktop so that you can plot easily
xvfb-run -d matlab -nosplash -nodisplay -nodesktop -r Spock_Preprocessing_Pipeline('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\RestingStateHemo\Mouse431_10_16_2019\431_10_16_2019_1\431_10_16_2019_1_MMStack_Pos0.ome.tif','Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\RestingStateHemo\Mouse431_10_16_2019\431_10_16_2019_1\prepro_log.mat');