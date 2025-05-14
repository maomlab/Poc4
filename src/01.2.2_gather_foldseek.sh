

FOLDSEEK_M8_PATH="data/foldseek/Cauris_884_AF2_mmseqs_results_tt3_maeaoRTmWfE9qU_0Cs202D9oywnNW6T1XA"
FOLDSEEK_FASTA_PATH="intermediate_data/foldseek/Cauris_884_AF2_mmseqs_results_tt3_maeaoRTmWfE9qU_0Cs202D9oywnNW6T1XA"

mkdir -p ${FOLDSEEK_FASTA_PATH}

for m8_fname in $(find "${FOLDSEEK_M8_PATH}" -name "*.m8");
do
  fasta_fname="${FOLDSEEK_FASTA_PATH}/$(basename ${m8_fname%.*}).fasta"
  awk -f src/m8_to_fasta.awk "${m8_fname}" > "${fasta_fname}"
done
