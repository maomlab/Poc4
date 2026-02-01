

import os
import pathlib
from src import embed_ESM2
import pyarrow.parquet
import pyarrow as pa
model_name = 'esm2_t33_650M_UR50D'


fasta_file = pathlib.Path(
    "product/frame2seq_design/designs_fixinterface_T1.0_af2_v2_pldd2_gt96_score_lt0.6_20251103/designs.fasta")
embedding_path = pathlib.Path(
    "intermediate_data/ESM2_embeddings/designs_fixinterface_T1.0_af2_v2_pldd2_gt96_score_lt0.6_20251103.parquet")
embedding_path.parent.mkdir(parents=True, exist_ok=True)
dataset = embed_ESM2.extract_embeddings(model_name, fasta_file, verbose=True)
pa.parquet.write_table(
    table = pa.Table.from_pandas(dataset),
    where = embedding_path)



fasta_file = pathlib.Path(
    "intermediate_data/AlphaFold3/afs_dockq.fasta")
embedding_path = pathlib.Path(
    "intermediate_data/AlphaFold3/afs_dockq_ems_embedding.parquet")
embedding_path.parent.mkdir(parents=True, exist_ok=True)
dataset = embed_ESM2.extract_embeddings(model_name, fasta_file, verbose=True)
pa.parquet.write_table(
    table = pa.Table.from_pandas(dataset),
    where = embedding_path)


fasta_file = pathlib.Path(
    "intermediate_data/AlphaFold3/afs_dockq_test.fasta")
embedding_path = pathlib.Path(
    "intermediate_data/AlphaFold3/afs_dockq_test_ems_embedding.parquet")
embedding_path.parent.mkdir(parents=True, exist_ok=True)
dataset = embed_ESM2.extract_embeddings(model_name, fasta_file, verbose=True)
pa.parquet.write_table(
    table = pa.Table.from_pandas(dataset),
    where = embedding_path)
