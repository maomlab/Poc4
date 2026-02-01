
metadata <- seqinr::read.fasta(
    file = "intermediate_data/ProteinMPNN/884_T0.1_5000/seqs/Cauris_AF-A0A2H1A4M0-F1-model_v4.fa",
    seqtype = "AA")

# the first one is just the input sequence, doesn't fit h the format

ref_metadata <- data.frame(
    Annot = metadata[[1]] |>
        seqinr::getAnnot() |>
        stringr::str_replace("^>", "")) |>
    tidyr::separate_rows(
        Annot,
        sep = ", ") |>
    dplyr::filter(dplyr::row_number() > 1) |> # first one is just the seq_id
    tidyr::separate_wider_delim(
        Annot,
        delim = "=",
        names = c("key", "value")) |>
    tidyr::pivot_wider(
        names_from = "key") |>
    dplyr::mutate(
        score = score |> as.numeric(),
        global_score = global_score |> as.numeric(),
        seed = seed |> as.numeric())
    
design_metadata <- metadata[2:length(metadata)] |>
    purrr::map_df(
        .f = function(seq) {
            
        data.frame(
            Annot = seqinr::getAnnot(seq) |>
                stringr::str_replace("^>", "")) |>
            tidyr::separate_rows(
                Annot,
                sep = ", ") |>
            tidyr::separate_wider_delim(
                Annot,
                delim = "=",
                names = c("key", "value")) |>
            dplyr::mutate(
                value = value |> as.numeric()) |>
            tidyr::pivot_wider(
                names_from = "key")
    })

# add in a row at the top to make it have the right length
metadata <- dplyr::bind_rows(
    ref_metadata |>
    dplyr::transmute(
        T = NA,
        sample = 0,
        score = score,
        global_score = global_score,
        seq_recover = 1),
    design_metadata)
            


metadata |>
    readr::write_tsv(
        "intermediate_data/ProteinMPNN/884_T0.1_5000/seqs/Cauris_AF-A0A2H1A4M0-F1-model_v4_metadata.tsv")

        
    
