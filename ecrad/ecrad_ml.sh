#!/usr/bin/env bash

set -eu

data_dir="${1:-"/scratch/snx3000/gbertoli/icon/regress-radiation/saved_models"}"

default_model="patch_shrNhbrWeightslvl0_all_heights_BN_mse_datav2IG_SSD_MT_featNormv2_fit_seq_MergedLwupLwdnSwupSwdn/serialized_model_steps_182400_238640_80560_188480"
model_dir="${data_dir}/${default_model}"


eng_type=tf_c

# run inference
cmd=(./multiple_fields_model.exe "${model_dir}" "${eng_type}")
echo "${cmd[@]}"
srun "${cmd[@]}"
