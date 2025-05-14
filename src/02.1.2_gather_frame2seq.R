
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
    

metadata <- frame2seq_metadata(
    path = "intermediate_data/frame2seq_design/Cauris_884_AF2_Frame2seq_T1.0_5000_seqs.fasta")
            
metadata |>
    readr::write_tsv(
        "intermediate_data/frame2seq_design/Cauris_884_AF2_Frame2seq_T1.0_5000_metadata.tsv")

for (temp in c("10", "1.0", "0.1", "0.01")) {
   frame2seq_add_index(
        input_path = paste0(
            "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T", temp, "_5000_seqs.fasta"),
        output_path = paste0(
            "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T", temp, "_5000_seqs_clean.fasta"))
}


for (temp in c("10", "1.0", "0.1", "0.01")) {
    metadata <- frame2seq_metadata(
        path = paste0(
            "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T", temp, "_5000_seqs.fasta"))
    metadata |>
        readr::write_tsv(
            paste0("intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T", temp, "_5000_metadata.tsv"))
}


########

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

