#!/usr/bin/env bash
#SBATCH -J 'xvalid'
#SBATCH -o out/DeconvMethods%j.out
#SBATCH -p all
#SBATCH -t 239
#SBATCH --exclude=redshirt-n[12-49]
#SBATCH --mem-per-cpu=64G
#SBATCH --array=1-7
#SBATCH --mail-type=END
#SBATCH --mail-user=<temp@princeton.edu>

cd "/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis/Deconvolution"
module load matlab/R2019b
matlab -batch "CompareDeconvolutionMethods(1,$SLURM_ARRAY_TASK_ID)"
