#!/bin/bash

PROTEIN_MPNN_PATH="${HOME}/turbo/opt/ProteinMPNN"

folder_with_pdbs="data/"
output_dir="intermediate_data/ProteinMPNN/884_T0.1_5000"
mkdir -p ${output_dir}

path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
python ${PROTEIN_MPNN_PATH}/helper_scripts/parse_multiple_chains.py \
       --input_path=$folder_with_pdbs \
       --output_path=$path_for_parsed_chains

python ${PROTEIN_MPNN_PATH}/protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --out_folder $output_dir \
        --num_seq_per_target 5000 \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1

# chain_id_jsonl is NOT loaded
# ----------------------------------------
# fixed_positions_jsonl is NOT loaded
# ----------------------------------------
# pssm_jsonl is NOT loaded
# ----------------------------------------
# omit_AA_jsonl is NOT loaded
# ----------------------------------------
# bias_AA_jsonl is NOT loaded
# ----------------------------------------
# tied_positions_jsonl is NOT loaded
# ----------------------------------------
# bias by residue dictionary is not loaded, or not provided
# ----------------------------------------
# discarded {'bad_chars': 0, 'too_long': 0, 'bad_seq_length': 0}
# ----------------------------------------
# Number of edges: 48
# Training noise level: 0.2A
# Generating sequences for: Cauris_AF-A0A2H1A4M0-F1-model_v4
# 5000 sequences of length 111 generated in 4677.9136 seconds
