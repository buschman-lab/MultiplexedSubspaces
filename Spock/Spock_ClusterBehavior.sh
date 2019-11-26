#!/usr/bin/env bash

#name job, output to slurm file, use partition all, run for 60 minutes and use 4GB of ram
#SBATCH -J 'pheno'
#SBATCH -o out/cluster_behavior_%j.out
#SBATCH -p all
#SBATCH -t 5
#SBATCH --mem-per-cpu=10G  
#SBATCH --array=1
#SBATCH --mail-type=FAIL                         # Event(s) that triggers email notification (BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=<camdenm@princeton.edu>    # Destination email address

# run the latest matlab
module load matlab/R2018a

#move to the repo
cd "/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/Spock/"

#run a paranoid version of matlab that crashes gracefully and records why just incase.
# in addition, xvfb-run also creates a virtual desktop so that you can plot easily
xvfb-run -d matlab -nosplash -nodisplay -nodesktop -r "try; Spock_ClusterBehavior($SLURM_ARRAY_TASK_ID,30); catch me; fprintf('Error: %s / %s\n',me.identifier,me.message); for k=1:length(me.stack), me.stack(k), end; end; exit"
