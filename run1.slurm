#!/bin/bash
#SBATCH --job-name="FinalNegBin1"
#SBATCH --time=144:00:00
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:0
#SBATCH -a 1
#3888
source ~/.bashrc
export TMPDIR=$WORKDIR
module load SC 
module load gcc/9.2.0
export STAN_NUM_THREADS=-1
time ./FinalNegBinMultiPhi sample num_samples=50000 num_warmup=20000 thin=1 adapt delta=0.8 random seed=1  data file=MultiRegionData output file=/scratch1/gri311/FinalModels/NegBin/output1.csv


