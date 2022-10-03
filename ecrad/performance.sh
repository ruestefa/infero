#!/usr/bin/env bash
#SBATCH --job-name=perf_infero 
#SBATCH --output=perf.log   
#SBATCH --error=perf.log   

set -eu

datadir=$1

model=$datadir/patch_shrNhbrWeightslvl0_all_heights_BN_mse_datav2IG_SSD_MT_featNormv2_fit_seq_MergedLwupLwdnSwupSwdn/serialized_model_steps_182400_238640_80560_188480

# run inference
srun ./performance.exe $model


