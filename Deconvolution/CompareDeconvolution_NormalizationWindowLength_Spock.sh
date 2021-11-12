#!/usr/bin/env bash
#SBATCH -J 'xvalid'
#SBATCH -o out/DeconvMethods%j.out
#SBATCH -p all
#SBATCH -t 239
#SBATCH --exclude=redshirt-n[12-49]
#SBATCH --mem-per-cpu=16G
#SBATCH --array=1-48
#SBATCH --mail-type=END
#SBATCH --mail-user=<temp@princeton.edu>

cd "/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis/Deconvolution"
module load matlab/R2019b
matlab -batch "CompareDeconvolution_NormalizationWindowLength($SLURM_ARRAY_TASK_ID,'std')"

# Camden this should be array 1-48 if you ran 6 recs and 8 windows
