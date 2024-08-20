#!/bin/bash

# Batch size was 50 on the ACAN run
batch_size=50
mem=5
output="acan_replicate"
region_file="$WORKSPACE_BUCKET/saraj/vntr_reference_panel/ACAN_region_1k_each_side.bed"
vntr_id="290964"

time python genotype_aou.py \
  --dataset lrwgs \
  --sample-count 4 \
  --batch-size 2 \
  --output-name $output \
  --region-file $region_file \
  --mem $mem \
  --vntr-id $vntr_id
