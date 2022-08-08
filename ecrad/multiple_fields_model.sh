#!/usr/bin/env bash

set -eu

datadir=/scratch/snx3000/juckerj/WG_1_tasks/ecrad/regress-radiation/saved_models
exedir=/scratch/snx3000/juckerj/WG_1_tasks/ecrad/infero/ecrad


eng_type=tf_c

model_path=$datadir/patch_shrNhbrWeightslvl1_all_heights_BN_mse_energyconsv1W1e-4_lw_datav2IG_SSD_MT_featNormv2/serialized_model_step_427000

# run inference
$exedir/multiple_fields_model.exe \
  $model_path \
  $eng_type \


