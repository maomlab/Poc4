
load_af2_scores <- function(
    metadata_path,
    fasta_path,
    scores_path,
    job_tag) {
    
    metadata <- readr::read_tsv(
        metadata_path,
        show_col_types = FALSE)
    metadata <- metadata |>
        dplyr::mutate(
            structure_id = paste0(job_tag, "_", seq_id))

    seqs <- seqinr::read.fasta(file = fasta_path, seqtype = "AA") |>
        purrr::map_df(
            .f = function(seq) {
                data.frame(
                    seq_id = seqinr::getAnnot(seq) |> stringr::str_replace("^>", ""),
                    fasta = paste0(seq, collapse = ""))
            })

    af2_scores <- readr::read_tsv(
        scores_path,
        show_col_types = FALSE)
    
    af2_scores <- af2_scores |>
        dplyr::left_join(
            metadata,
            by = "structure_id") |>
        dplyr::mutate(
            seq_id = structure_id |> stringr::str_extract("seq_[0-9]+$")) |>
        dplyr::left_join(
            seqs,
            by = "seq_id")

    af2_scores <- af2_scores |>
        dplyr::mutate(
            job_tag = job_tag,
            .before = 1)
    af2_scores
}

af2_scores <- dplyr::bind_rows(
    load_af2_scores(
        metadata_path = "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T0.1_5000_metadata.tsv",
        fasta_path = "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T0.1_5000_seqs_clean.fasta",
        scores_path = "intermediate_data/parafold/af2_scores.tsv",
        job_tag = "Frame2Seq_T0.1"),
    load_af2_scores(
        metadata_path = "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_5000_metadata.tsv",
        fasta_path = "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_5000_seqs_clean.fasta",
        scores_path = "intermediate_data/parafold/Frame2Seq_fixinterface_T1.0/af2_scores.tsv",
        job_tag = "Frame2Seq_fixinterface_T1.0"))

af2_scores |>
    readr::write_tsv(
        "product/frame2seq_design/af2_scores_20250506.tsv")



plot <- ggplot2::ggplot(
    data = af2_scores) +
    ggplot2::theme_bw() +
    ggplot2::geom_histogram(
        mapping = ggplot2::aes(
            x = af2_plddt),
        bins = 50) +
    ggplot2::scale_x_continuous(
        "AlphaFold2 pLDDT") +
    ggplot2::facet_grid(
        cols = dplyr::vars(rank),
        rows = dplyr::vars(job_tag)) +
    ggplot2::theme(
        legend.position = "bottom") +
    ggplot2::ggtitle(
        "Cauris 884 with IRC25 AF3",
        subtitle = "By Rank and Design Protocols")

ggplot2::ggsave(
    filename = "product/frame2seq_design/plddt_histogram_rank_protocol_20241027.pdf",
    width = 7,
    height = 4,
    useDingbats = FALSE)




plot <- ggplot2::ggplot(
    data = af2_scores) +
    ggplot2::theme_bw() +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = score,
            y = af2_plddt,
            color = recovery_percent),
        alpha = 0.8,
        shape = 16) +
    ggplot2::facet_wrap(
        facets = dplyr::vars(job_tag)) +
    ggplot2::scale_x_continuous(
        "Frame2Seq Score") +
    ggplot2::scale_y_continuous(
        "AlphaFold2 pLDDT") +
    viridis::scale_color_viridis(
        "Frame2Seq SeqRecovery") +
    ggplot2::theme(
        legend.position = "bottom") +
    ggplot2::ggtitle(
        "Cauris 884 with IRC25 AF3")

ggplot2::ggsave(
    filename = "product/frame2seq_design/plddt_vs_score_20241027.pdf",
    width = 7,
    height = 4,
    useDingbats = FALSE)



plot <- ggplot2::ggplot(
    data = af2_scores) +
    ggplot2::theme_bw() +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = recovery_percent,
            y = af2_plddt,
            color = score),
        alpha = 0.8,
        shape = 16) +
    ggplot2::facet_wrap(
        facets = dplyr::vars(job_tag)) +
    ggplot2::scale_x_continuous(
        "Frame2Seq Recovery Percent") +
    ggplot2::scale_y_continuous(
        "AlphaFold2 pLDDT") +
    viridis::scale_color_viridis(
        "Frame2Seq Score") +
    ggplot2::theme(
        legend.position = "bottom") +
    ggplot2::ggtitle(
        "Cauris 884 with IRC25 AF3")

ggplot2::ggsave(
    filename = "product/frame2seq_design/plddt_vs_recovery_20241027.pdf",
    width = 7,
    height = 4,
    useDingbats = FALSE)


plot <- ggplot2::ggplot(
    data = af2_scores |>
        dplyr::filter(af2_plddt > 85)) +
    ggplot2::theme_bw() +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = score,
            y = recovery_percent,
            color = af2_plddt),
        alpha = 0.8,
        shape = 16) +
    ggplot2::facet_wrap(
        facets = dplyr::vars(job_tag)) +
    ggplot2::scale_x_continuous(
        "Frame2Seq Score") +
    ggplot2::scale_y_continuous(
        "Frame2Seq Recovery Percent") +
    viridis::scale_color_viridis(
        "AlphaFold2 pLDDT") +
    ggplot2::theme(
        legend.position = "bottom") +
    ggplot2::ggtitle(
        "Cauris 884 with IRC25 AF3")

ggplot2::ggsave(
    filename = "product/frame2seq_design/score_vs_recovery_20241027.pdf",
    width = 7,
    height = 4,
    useDingbats = FALSE)


######################
collate_designs <- function(
    designs,
    output_path) {

    if (!dir.exists(output_path)) {
        cat("creating ", output_path, " ...\n", sep = "")
        dir.create(output_path)
    }
    
    wt_884_fasta <- readr::read_file(
        file = "data/B9J08_000884.fasta")
    
    
    
    con <- file(paste0(output_path, "/designs.fasta"), "w")
    writeLines(wt_884_fasta, con = con)
    designs$fasta |>
        paste0(collapse = "\n") |>
        writeLines(con = con)
    close(con)
    
    
    designs |>
        dplyr::select(-fasta) |>
        readr::write_tsv(
            paste0(output_path, "/summary.tsv"))
    
    file.copy(
        from = "intermediate_data/AF3/fold_884_irc25/fold_884_irc25_model_0.pdb",
        to = paste0(output_path, "/Cauris_884_irc25_AF3.pdb"))
    
    file.copy(
        from = "intermediate_data/AF3/fold_884_irc25/fold_884_irc25_full_data_0.json",
        to = paste0(output_path, "/fold_884_irc25_full_data_0.json"))
    
    file.copy(
        from = "intermediate_data/AF3/fold_884_irc25/fold_884_irc25_summary_confidences_0.json",
        to = paste0(output_path, "/fold_884_irc25_summary_confidences_0.json"))
    
    
    designs |>
        dplyr::mutate(
            pdb_path = af2_results_fname |>
                stringr::str_replace("result_model", "unrelaxed_model") |>
                stringr::str_replace(".pkl", ".pdb")) |>
        dplyr::rowwise() |>
        dplyr::do({
            pdb_path <- .$pdb_path[1]
            structure_id <- .$structure_id[1]
            file.copy(
                from = pdb_path,
                to = paste0(output_path, "/", structure_id, ".pdb"))
            data.frame()
        })
}



# top_designs <- af2_scores |>
#     dplyr::filter(
#         job_tag == "Frame2Seq_T0.1",
#         af2_plddt > 96,
#         score < 0.6) |>
#     dplyr::arrange(desc(af2_plddt)) |>
#     dplyr::distinct(seq_id, .keep_all = TRUE) |>
#     dplyr::arrange(recovery_percent) |>
#     head(30)
# output_path <- "product/frame2seq_design/designs_T0.1_af2_pldd2_gt96_score_lt0.6_top30_20241022"
# 
# collate_designs(
#     top_designs = top_designs,
#     output_path = output_path)


top_designs <- af2_scores |>
    dplyr::filter(
        job_tag == "Frame2Seq_fixinterface_T1.0",
        af2_plddt > 96,
        score < 0.6) |>
    dplyr::arrange(desc(af2_plddt)) |>
    dplyr::distinct(seq_id, .keep_all = TRUE) |>
    dplyr::arrange(recovery_percent) |>
    head(30)
output_path <- "product/frame2seq_design/designs_fixinterface_T1.0_af2_pldd2_gt96_score_lt0.6_top30_20241022"

collate_designs(
    designs = top_designs,
    output_path = output_path)




