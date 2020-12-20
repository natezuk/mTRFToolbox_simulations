#!/bin/bash
#SBATCH -t 2-00:00:00
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -n 1

module load matlab
matlab -nodesktop -nosplash -r "Test_TRFestimate_SNRnTr_delta"
