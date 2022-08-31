#!/usr/bin/env bash

set -eu

datadir=$1
eng_type=tf_c

lw_model=$datadir/patch_shrNhbrWeightslvl1_all_heights_BN_mse_energyconsv1W1e-4_lw_datav2IG_SSD_MT_featNormv2/serialized_model_step_427000
sw_model=$datadir/patch_shrNhbrWeightslvl1_all_heights_BN_mse_energyconsv1W1e-4_sw_datav2IG_SSD_MT_featNormv2/serialized_model_step_491000

# run inference
./multiple_fields_model.exe \
  $lw_model \
  $sw_model \
  $eng_type \


