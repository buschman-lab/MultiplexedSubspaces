#!/usr/bin/env bash
#SBATCH -J 'fit_GLM'
#SBATCH -o out/Fit_GLM%j.out
#SBATCH -p all
#SBATCH -t 239
#SBATCH --exclude=redshirt-n[12-49]
#SBATCH --mem-per-cpu=64G
#SBATCH --array=1
#SBATCH --mail-type=END
#SBATCH --mail-user=<temp@princeton.edu>

cd "/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis/Deconvolution"
module load matlab/R2019b
matlab -batch "GenerateDeconvolutionGLM('std',1)"
