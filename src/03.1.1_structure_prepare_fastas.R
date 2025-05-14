

library(seqinr)


####
# Predict 884 variants alone
prepare_fastas <- function(
    fasta_fname,
    sequence_labels_fname,
    dataset_tag,
    output_path,
    verbose = FALSE) {

    fasta_files <- data.frame(
        fname = fasta_fname)

    sequence_labels <- dplyr::bind_rows(
        readr::read_tsv(
            sequence_labels_fname,
            show_col_types = FALSE) |>
        dplyr::mutate(dataset = dataset_tag))
    

    sequences <- fasta_files |>
        dplyr::rowwise() |>
        dplyr::do({
            data <- .
            sequences <- seqinr::read.fasta(
                file = data$fname[1],
                seqtype = "AA")
    
            tibble::tibble(
                sequence = sequences |> purrr::map_chr(~paste0(.x, collapse = ""))) |>
                dplyr::mutate(
                    sequence = sequence |> stringr::str_replace("[*]", "")) |>
                dplyr::mutate(
                    fname = data$fname[1],
                    .before = 1)
        }) 
    
    sequences <- dplyr::bind_cols(
        sequences,
        sequence_labels)
    
    if (!dir.exists(output_path)) {
        cat("Creating '", output_path, "' ...\n", sep = "")
        dir.create(output_path, recursive = TRUE)
    }
    
    sequences |>
        dplyr::rowwise() |>
        dplyr::do({
            data <- .
    
            fasta_path <- paste0(
                output_path,
                data$dataset[1], "_",
                data$seq_id[1])
            if (!dir.exists(fasta_path)) {
                cat("Creating fasta path folder '", fasta_path, "'\n", sep = "")
                dir.create(fasta_path, recursive = TRUE)
            }
    
            cat(
                "> ", data$dataset[1], "_", data$seq_id[1], "\n",
                data$sequence[1], "\n",
                file = paste0(fasta_path, "/input.fasta"),
                sep = "")
            data.frame()
        })

}

# Predict 884 variants alone
prepare_fastas(
    fasta_fname = "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T0.1_5000_seqs_clean.fasta",
    sequence_labels_fname = "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T0.1_5000_metadata.tsv",
    dataset_tag = "Frame2Seq_T0.1",
    output_path = "intermediate_data/parafold/input/",
    verbose = TRUE)


# Predict 884 variants alone
prepare_fastas(
    fasta_fname = "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_5000_seqs_clean.fasta",
    sequence_labels_fname = "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_5000_metadata.tsv",
    dataset_tag = "Frame2Seq_fixinterface_T1.0",
    output_path = "intermediate_data/parafold/Frame2Seq_fixinterface_T1.0/input",
    verbose = TRUE)

