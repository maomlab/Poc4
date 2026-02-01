
library(seqinr)

frame2seq_add_index <- function(
    input_path,
    output_path) {
    cat("Reading sequences from '", input_path, "'\n", sep = "")
    cat("Writing sequences to '", output_path, "'\n", sep = "")    
    seqs <- seqinr::read.fasta(file = input_path, seqtype = "AA")
    seqs |> seqinr::write.fasta(
        file.out = output_path,
        names = paste0("seq_", 1:length(seqs)))
}

frame2seq_metadata <- function(path) {
    cat("Reading sequences from '", path, "'\n", sep = "")
    metadata <- seqinr::read.fasta(file = path, seqtype = "AA") |>
        purrr::map_df(
            .f = function(seq) {
                data.frame(
                    Annot = seqinr::getAnnot(seq) |>
                        stringr::str_replace("^>", "")) |>
                    tidyr::separate_rows(
                        Annot,
                    sep = " ") |>
                    tidyr::separate_wider_delim(
                        Annot,
                        delim = "=",
                        names = c("key", "value")) |>
                    tidyr::pivot_wider(
                        names_from = "key")
            })

    metadata <- metadata |>
        dplyr::mutate(
            recovery_percent = recovery |>
                stringr::str_replace("%", "") |>
                as.numeric(),
            score = score |> as.numeric()) |>
        dplyr::select(recovery_percent, score) |>
        dplyr::mutate(
            seq_id = paste0("seq_", 1:dplyr::n()),
            .before = 1)
    metadata
}

####################
# Structure: AlphaFold2 Cauris Poc4
# Temperature: 1.0
# Number of designs: 5000
metadata <- frame2seq_metadata(
    path = "intermediate_data/frame2seq_design/Cauris_884_AF2_Frame2seq_T1.0_5000_seqs.fasta")
metadata |>
    readr::write_tsv(
        "intermediate_data/frame2seq_design/Cauris_884_AF2_Frame2seq_T1.0_5000_metadata.tsv")

####################
# Structure: AlphaFold3 Cauris Poc4 Irc25 dimer
#            data/AF3/fold_884_irc25/fold_884_irc25_model_0.pdb
# Temperature: 10-0.1
# Number of designs: 5000
for (temp in c("10", "1.0", "0.1", "0.01")) {
   frame2seq_add_index(
        input_path = paste0(
            "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T", temp, "_5000_seqs.fasta"),
        output_path = paste0(
            "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T", temp, "_5000_seqs_clean.fasta"))

    metadata <- frame2seq_metadata(
        path = paste0(
            "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T", temp, "_5000_seqs.fasta"))
    metadata |>
        readr::write_tsv(
            paste0("intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T", temp, "_5000_metadata.tsv"))
}




#########################
# Structure: AlphaFold3 Cauris Poc4 Irc25 dimer
#            data/AF3/fold_884_irc25/fold_884_irc25_model_0.pdb
# Fixed conserved residues
#   * [18, 24, 26, 42, 44, 52, 57, 67, 76, 78, 79, 83, 85, 86, 89]
# Temperature: 1.0
# Number of designs: 50000

frame2seq_add_index(
     input_path = paste0(
         "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixconserved_T1.0_50000_seqs.fasta"),
     output_path = paste0(
         "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixconserved_T1.0_50000_seqs_clean.fasta"))

metadata <- frame2seq_metadata(
    path = paste0(
        "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixconserved_T1.0_50000_seqs.fasta"))
metadata |>
    readr::write_tsv(
        paste0("intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixconserved_T1.0_50000_metadata.tsv"))


#########################
# Structure: AlphaFold3 Cauris Poc4 Irc25 dimer
#            data/AF3/fold_884_irc25/fold_884_irc25_model_0.pdb
# Fixed interface residues
#    [1, 3, 5, 12, 14, 16, 18, 24, 25, 26, 28, 30, 32, 33, 34, 35, 36, 37, 38, 39, 42, 44, 55, 56, 57, 58, 83, 85, 87]
# Temperature: 1.0
# Number of designs: 5000
# First round of predictions, which we used for experimental testing
frame2seq_add_index(
     input_path = paste0(
         "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_5000_seqs.fasta"),
     output_path = paste0(
         "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_5000_seqs_clean.fasta"))
metadata <- frame2seq_metadata(
    path = paste0(
        "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_5000_seqs.fasta"))
metadata |>
    readr::write_tsv(
        paste0("intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_5000_metadata.tsv"))


#########################
# Structure: AlphaFold3 Cauris Poc4 Irc25 dimer
#            data/AF3/fold_884_irc25/fold_884_irc25_model_0.pdb
# Fixed interface residues
#    [1, 3, 5, 12, 14, 16, 18, 24, 25, 26, 28, 30, 32, 33, 34, 35, 36, 37, 38, 39, 42, 44, 55, 56, 57, 58, 83, 85, 87]
# Temperature: 1.0
# Number of designs: 500000
# Second round of predictions ("v2"), which we used for sequence-structure analysis
frame2seq_add_index(
     input_path = paste0(
         "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_500000_v2_seqs.fasta"),
     output_path = paste0(
         "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_500000_v2_seqs_clean.fasta"))
metadata <- frame2seq_metadata(
    path = paste0(
        "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_500000_v2_seqs.fasta"))
metadata |>
    readr::write_tsv(
        paste0("intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_500000_v2_metadata.tsv"))

#######################################################
# Filter round 2 designs for the top designs => "good"
# Frame2Seq score < 0.6
# lowest 10000 by percent recovery
metadata <- readr::read_tsv(
    paste0("intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_500000_v2_metadata.tsv"),
    show_col_types = FALSE)
sequences <- seqinr::read.fasta(
    file = "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_500000_v2_seqs_clean.fasta",
    seqtype = "AA")

metadata_good_score <- metadata |>
    dplyr::filter(score < 0.6) |>
    dplyr::arrange(recovery_percent) |>
    head(10000)
sequences_good_score <- sequences[names(sequences) %in% metadata_good_score$seq_id]

metadata_good_score |>
    readr::write_tsv(
        paste0("intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_500000_v2_metadata_good_score.tsv"))
sequences_good_score |>
    seqinr::write.fasta(
        file.out = "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_500000_v2_seqs_good_score.fasta",
        names = names(sequences_good_score))
