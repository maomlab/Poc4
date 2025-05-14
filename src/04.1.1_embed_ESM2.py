

import os
import pathlib

from src import embed_ESM2

import pyarrow.parquet
import pyarrow as pa

model_name = 'esm2_t33_650M_UR50D'



# ProteinMPNN Cauris_884_AF2 T1.0 5000
fasta_file = pathlib.Path(
    "intermediate_data/ProteinMPNN/884_T0.1_5000/seqs/Cauris_AF-A0A2H1A4M0-F1-model_v4.fa")
embedding_path = pathlib.Path(
    "intermediate_data/ESM2_embeddings/ProteinMPNN_Cauris_884_AF2_T0.1_5000.parquet")
embedding_path.parent.mkdir(parents=True, exist_ok=True)
dataset = embed_ESM2.extract_embeddings(model_name, fasta_file, verbose=True)
pa.parquet.write_table(
    table = pa.Table.from_pandas(dataset),
    where = embedding_path)



# Frame2Seq Cauris_884_AF2 T1.0 5000
fasta_file = pathlib.Path(
    "intermediate_data/frame2seq_design/Cauris_884_AF2_Frame2seq_T1.0_5000_seqs_clean.fasta")
embedding_path = pathlib.Path(
    "intermediate_data/ESM2_embeddings/frame2seq_Cauris_884_AF2_T1.0_5000.parquet")
embedding_path.parent.mkdir(parents=True, exist_ok=True)
dataset = embed_ESM2.extract_embeddings(model_name, fasta_file, verbose=True)
pa.parquet.write_table(
    table = pa.Table.from_pandas(dataset),
    where = embedding_path)

# Frame2Seq Cauris_884_irc25_AF3 5000
for temp in ["10", "1.0", "0.1", "0.01"]:
    print(f"Embedding frame2seq for temp {temp} ...")
    fasta_file = pathlib.Path(
        f"intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T{temp}_5000_seqs_clean.fasta")
    embedding_path = pathlib.Path(
        f"intermediate_data/ESM2_embeddings/frame2seq_Cauris_884_irc25_AF3_T{temp}_5000.parquet")
    embedding_path.parent.mkdir(parents=True, exist_ok=True)
    dataset = embed_ESM2.extract_embeddings(model_name, fasta_file, verbose=True)
    pa.parquet.write_table(
        table = pa.Table.from_pandas(dataset),
        where = embedding_path)
        


# Frame2Seq Cauris_884_AF2 fix interface T1.0 5000
fasta_file = pathlib.Path(
    "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_5000_seqs_clean.fasta")
embedding_path = pathlib.Path(
    "intermediate_data/ESM2_embeddings/frame2seq_Cauris_884_irc25_AF3_fixinterface_T1.0_5000.parquet")
embedding_path.parent.mkdir(parents=True, exist_ok=True)
dataset = embed_ESM2.extract_embeddings(model_name, fasta_file, verbose=True)
pa.parquet.write_table(
    table = pa.Table.from_pandas(dataset),
    where = embedding_path)

# Frame2Seq Cauris_884_AF2 fix conserved T1.0 50000
fasta_file = pathlib.Path(
    "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixconserved_T1.0_50000_seqs_clean.fasta")
embedding_path = pathlib.Path(
    "intermediate_data/ESM2_embeddings/frame2seq_Cauris_884_irc25_AF3_fixconserved_T1.0_50000.parquet")
embedding_path.parent.mkdir(parents=True, exist_ok=True)
dataset = embed_ESM2.extract_embeddings(model_name, fasta_file, verbose=True)
pa.parquet.write_table(
    table = pa.Table.from_pandas(dataset),
    where = embedding_path)




# FoldSeek Cauris_884_AF2
for fasta_file in pathlib.Path(
    "intermediate_data/foldseek/Cauris_884_AF2_mmseqs_results_tt3_maeaoRTmWfE9qU_0Cs202D9oywnNW6T1XA/").glob("*fasta"):
    embedding_path = pathlib.Path(
        f"intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/{fasta_file.stem}.parquet")
    embedding_path.parent.mkdir(parents=True, exist_ok=True)
    dataset = embed_ESM2.extract_embeddings(model_name, fasta_file, verbose=True)
    pa.parquet.write_table(
        table = pa.Table.from_pandas(dataset),
        where = embedding_path)

# PFAM PF10448
fasta_file = pathlib.Path(
    "data/pfam/protein-matching-PF10448.fasta")
embedding_path = pathlib.Path(
    "intermediate_data/EMS2_embeddings/pfam_PF10448.parquet")
embedding_path.parent.mkdir(parents=True, exist_ok=True)
dataset = embed_ESM2.extract_embeddings(model_name, fasta_file, verbose=True)
pa.parquet.write_table(
    table = pa.Table.from_pandas(dataset),
    where = embedding_path)


