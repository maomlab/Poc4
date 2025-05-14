
import shutil
from frame2seq import Frame2seqRunner
runner = Frame2seqRunner()

# Temperature 10
runner.design(
    pdb_file = "data/AF3/fold_884_irc25/fold_884_irc25_model_0.pdb",
    chain_id = "A",
    temperature = 10,
    num_samples = 5000,
    omit_AA = ["C"],
    fixed_positions=None,
    save_indiv_seqs=False,
    save_indiv_neg_pll=False,
    verbose=True)

shutil.move(
    src = "frame2seq_outputs/seqs/seqs.fasta",
    dst = "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T10_5000_seqs.fasta")


# Temperature 1
runner.design(
    pdb_file = "data/AF3/fold_884_irc25/fold_884_irc25_model_0.pdb",
    chain_id = "A",
    temperature = 1.0,
    num_samples = 5000,
    omit_AA = ["C"],
    fixed_positions=None,
    save_indiv_seqs=False,
    save_indiv_neg_pll=False,
    verbose=True)

shutil.move(
    src = "frame2seq_outputs/seqs/seqs.fasta",
    dst = "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T1.0_5000_seqs.fasta")



# Temperature 0.1
runner.design(
    pdb_file = "data/AF3/fold_884_irc25/fold_884_irc25_model_0.pdb",
    chain_id = "A",
    temperature = 0.1,
    num_samples = 5000,
    omit_AA = ["C"],
    fixed_positions=None,
    save_indiv_seqs=False,
    save_indiv_neg_pll=False,
    verbose=True)

shutil.move(
    src = "frame2seq_outputs/seqs/seqs.fasta",
    dst = "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T0.1_5000_seqs.fasta")


# Temperature 0.01
runner.design(
    pdb_file = "data/AF3/fold_884_irc25/fold_884_irc25_model_0.pdb",
    chain_id = "A",
    temperature = 0.01,
    num_samples = 5000,
    omit_AA = ["C"],
    fixed_positions=None,
    save_indiv_seqs=False,
    save_indiv_neg_pll=False,
    verbose=True)

shutil.move(
    src = "frame2seq_outputs/seqs/seqs.fasta",
    dst = "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T0.01_5000_seqs.fasta")


#####################################################3

runner.design(
    pdb_file = "data/AF3/fold_884_irc25/fold_884_irc25_model_0.pdb",
    chain_id = "A",
    temperature = 1.0,
    num_samples = 50000,
    omit_AA = ["C"],
    fixed_positions=[18, 24, 26, 42, 44, 52, 57, 67, 76, 78, 79, 83, 85, 86, 89],
    save_indiv_seqs=False,
    save_indiv_neg_pll=False,
    verbose=True)

shutil.move(
    src = "frame2seq_outputs/seqs/seqs.fasta",
    dst = "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixconserved_T1.0_5000_seqs.fasta")
