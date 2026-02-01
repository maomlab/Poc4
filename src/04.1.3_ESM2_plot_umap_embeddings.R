
# this needs to exactly match the that was used ot generate the embedding to make sure
# the ids line up
data <- dplyr::bind_rows(
    # PFAM PF10448",
    readr::read_tsv(
        "data/pfam/protein-matching-PF10448.tsv",
        show_col_types = FALSE) |>
    dplyr::transmute(
        dataset = "pfam",
        seq_id = Accession,
        prob = NA,
        score = NA,
        taxon = `Tax ID`,
        organism = `Tax Name`),

    # Foldseek 884
    readr::read_tsv(
        "intermediate_data/foldseek/Cauris_884_AF2_mmseqs_results_tt3_maeaoRTmWfE9qU_0Cs202D9oywnNW6T1XA/foldseek_metadata.tsv",
        show_col_types = FALSE) |>
    dplyr::transmute(
        dataset,
        seq_id = target,
        prob,
        score,
        taxon,
        organism)) |>
    dplyr::mutate(
        seq_label = dplyr::case_when(
            seq_id == "A0A0L0P8N2" ~ "CuIRC25",
            seq_id == "AF-A0A510P858-F1-model_v4 Hypothetical_protein" ~ "Cu884",
            dataset == "foldseek_afdb50" & seq_id == "AF-A0A202G644-F1-model_v4 Uncharacterized protein" ~ "ClPOC4",
            dataset == "foldseek_afdb50" & seq_id == "AF-A0A4P9ZGT9-F1-model_v4 Uncharacterized protein" ~ "MbPOC4",
            seq_id == "A0A1D8PFG6" ~ "CaPOC3",
            dataset == "foldseek_cath50" & seq_id == "af_A0A1D8PTK4_4_109_3.30.230.100" ~ "CaPOC4",
            dataset == "foldseek_afdb_proteome" & seq_id == "AF-Q12245-F1-model_v4 Proteasome chaperone 4" ~ "ScPOC4",
            TRUE ~ NA)) |>
    dplyr::bind_cols(
        arrow::read_parquet(
            "intermediate_data/pfam,foldseeka=1,b=1.8,n_neighbors=100,metric=euclidean_APDP_20240507/umap_embedding.parquet"))

output_path <- "product/UMAP_embedding/pfam,foldseeka=1,b=1.8,n_neighbors=100,metric=euclidean_APDP_20240507/"

if (!dir.exists(output_path)) {
    cat("creating dir ", output_path, "\n", sep = "")
    dir.create(output_path)
}

plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::coord_equal() +
    ggplot2::geom_point(
        data = data |> dplyr::filter(dataset |> stringr::str_detect("foldseek")),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            size = score^3,
            color = dataset),
        alpha = 0.9,
        shape = 16) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(dataset |> stringr::str_detect("pfam")),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            color = dataset),
        alpha = 0.9,
        shape = 16) +
    ggplot2::geom_point(
        data = data |>
            dplyr::filter(!is.na(seq_label)),            
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            color = dataset),
        size = 1.5,
        alpha = 1,
        shape = 16) +    
    ggrepel::geom_text_repel(
        data = data |>
            dplyr::filter(!is.na(seq_label)),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            label = seq_label),
        force = 10,
        min.segment.length = 0) +
    ggplot2::scale_x_continuous("UMAP1") +
    ggplot2::scale_y_continuous("UMAP2") +
    ggplot2::scale_size_continuous(
        "FoldSeek P-value",
        breaks = c(220, 410, 908)^3,
        labels = c("1e-6", "1e-10", "1e-22")) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::ggtitle(
        label = "UMAP Projection of EMS2 embeddings for Natural POC3 and POC4 Sequences",
        subtitle = "PFAM PF10448 and Cu884 FoldSeek Hits")


ggplot2::ggsave(
    filename = paste0(output_path, "/umap_embedding_20240507.pdf"),
    plot = plot,
    width = 10,
    height = 10,
    useDingbats = FALSE)

########################################


data <- dplyr::bind_rows(
    readr::read_tsv(
        "data/pfam/protein-matching-PF10448.tsv",
        show_col_types = FALSE) |>
    dplyr::transmute(
        dataset = "pfam",
        seq_id = Accession,
        prob = NA,
        score = NA,
        taxon = `Tax ID`,
        organism = `Tax Name`),
    readr::read_tsv(
        "intermediate_data/foldseek/Cauris_884_AF2_mmseqs_results_tt3_maeaoRTmWfE9qU_0Cs202D9oywnNW6T1XA/foldseek_metadata.tsv",
        show_col_types = FALSE) |>
    dplyr::transmute(
        dataset,
        seq_id = target,
        prob,
        score,
        taxon,
        organism),
    readr::read_tsv(
        "intermediate_data/Cauris_884_AF2_Frame2seq_T1.0_5000_metadata.tsv",
        show_col_types = FALSE) |>
    dplyr::transmute(
        dataset = "Frame2Seq",
        seq_id,
        prob = recovery_percent,
        score,
        taxon = NA_integer_,
        organism = NA_character_)) |>
    dplyr::mutate(
        seq_label = dplyr::case_when(
            seq_id == "A0A0L0P8N2" ~ "CuIRC25",
            seq_id == "AF-A0A510P858-F1-model_v4 Hypothetical_protein" ~ "Cu884",
            dataset == "foldseek_afdb50" & seq_id == "AF-A0A202G644-F1-model_v4 Uncharacterized protein" ~ "ClPOC4",
            dataset == "foldseek_afdb50" & seq_id == "AF-A0A4P9ZGT9-F1-model_v4 Uncharacterized protein" ~ "MbPOC4",
            seq_id == "A0A1D8PFG6" ~ "CaPOC3",
            dataset == "foldseek_cath50" & seq_id == "af_A0A1D8PTK4_4_109_3.30.230.100" ~ "CaPOC4",
            dataset == "foldseek_afdb_proteome" & seq_id == "AF-Q12245-F1-model_v4 Proteasome chaperone 4" ~ "ScPOC4",
            TRUE ~ NA)) |>
    dplyr::bind_cols(
        arrow::read_parquet(
            "intermediate_data/pfam,foldseek,frame2seq_monomera=1,b=1.8,n_neighbors=100,metric=euclidean_APDP_20240506/umap_embedding.parquet"))

output_path <- "product/UMAP_embedding/pfam,foldseek,frame2seqa=1,b=1.8,n_neighbors=100,metric=cosine_APDP_20240506/"

if (!dir.exists(output_path)) {
    cat("creating dir ", output_path, "\n", sep = "")
    dir.create(output_path)
}

plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::coord_equal() +
    ggplot2::geom_point(
        data = data |> dplyr::filter(dataset |> stringr::str_detect("foldseek")),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            size = score^3,
            color = dataset),
        alpha = 0.9,
        shape = 16) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(dataset %in% c("pfam", "Frame2Seq")),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            color = dataset),
        alpha = 0.9,
        shape = 16) +
    ggplot2::geom_point(
        data = data |>
            dplyr::filter(!is.na(seq_label)),            
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            color = dataset),
        size = 1.5,
        alpha = 1,
        shape = 16) +    
    ggrepel::geom_text_repel(
        data = data |>
            dplyr::filter(!is.na(seq_label)),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            label = seq_label),
        force = 10,
        min.segment.length = 0) +
    ggplot2::scale_x_continuous("UMAP1") +
    ggplot2::scale_y_continuous("UMAP2") +
    ggplot2::scale_size_continuous(
        "FoldSeek P-value",
        breaks = c(220, 410, 908)^3,
        labels = c("1e-6", "1e-10", "1e-22")) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::ggtitle(
        label = "UMAP Projection of EMS2 embeddings for Natural and Designed POC3 and POC4 Sequences",
        subtitle = "PFAM PF10448, Cu884 FoldSeek Hits, Cu884 Monomer Frame2Seq Design")
    


ggplot2::ggsave(
    filename = paste0(output_path, "/umap_embedding_natural_designed_20240506.pdf"),
    plot = plot,
    width = 10,
    height = 10,
    useDingbats = FALSE)


############################################################

data <- dplyr::bind_rows(
    readr::read_tsv(
        "data/pfam/protein-matching-PF10448.tsv",
        show_col_types = FALSE) |>
    dplyr::transmute(
        dataset = "pfam",
        seq_id = Accession,
        prob = NA,
        score = NA,
        taxon = `Tax ID`,
        organism = `Tax Name`),
    readr::read_tsv(
        "intermediate_data/foldseek/Cauris_884_AF2_mmseqs_results_tt3_maeaoRTmWfE9qU_0Cs202D9oywnNW6T1XA/foldseek_metadata.tsv",
        show_col_types = FALSE) |>
    dplyr::transmute(
        dataset,
        seq_id = target,
        prob,
        score,
        taxon,
        organism),
    dplyr::bind_rows(
        readr::read_tsv(
            "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T0.01_5000_metadata.tsv",
            show_col_types = FALSE) |>
            dplyr::mutate(dataset = "Frame2Seq_T0.01"),
        readr::read_tsv(
            "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T0.1_5000_metadata.tsv",
            show_col_types = FALSE) |>
            dplyr::mutate(dataset = "Frame2Seq_T0.1"),
        readr::read_tsv(
            "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T1.0_5000_metadata.tsv",
            show_col_types = FALSE) |>
        dplyr::mutate(dataset = "Frame2Seq_T1.0"),
        readr::read_tsv(
            "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_T10_5000_metadata.tsv",
            show_col_types = FALSE) |>
        dplyr::mutate(dataset = "Frame2Seq_T10")) |>
    dplyr::transmute(
        dataset,
        seq_id,
        prob = recovery_percent,
        score,
        taxon = NA_integer_,
        organism = NA_character_)) |>
    dplyr::mutate(
        seq_label = dplyr::case_when(
            seq_id == "A0A0L0P8N2" ~ "CuIRC25",
            seq_id == "AF-A0A510P858-F1-model_v4 Hypothetical_protein" ~ "Cu884",
            dataset == "foldseek_afdb50" & seq_id == "AF-A0A202G644-F1-model_v4 Uncharacterized protein" ~ "ClPOC4",
            dataset == "foldseek_afdb50" & seq_id == "AF-A0A4P9ZGT9-F1-model_v4 Uncharacterized protein" ~ "MbPOC4",
            seq_id == "A0A1D8PFG6" ~ "CaPOC3",
            dataset == "foldseek_cath50" & seq_id == "af_A0A1D8PTK4_4_109_3.30.230.100" ~ "CaPOC4",
            dataset == "foldseek_afdb_proteome" & seq_id == "AF-Q12245-F1-model_v4 Proteasome chaperone 4" ~ "ScPOC4",
            TRUE ~ NA)) |>
    dplyr::bind_cols(
        arrow::read_parquet(
            "intermediate_data/pfam,foldseek,frame2seq_monomera=1,b=1.8,n_neighbors=100,metric=cosine_20240521/umap_embedding.parquet"))

data <- data |>
    dplyr::mutate(
        dist_884 = sqrt((UMAP_1 - -24.82538)^2 + (UMAP_2 -  6.462263)^2))



output_path <- "product/UMAP_embedding/pfam,foldseek,frame2seq_monomera=1,b=1.8,n_neighbors=100,metric=cosine_20240521/"

if (!dir.exists(output_path)) {
    cat("creating dir ", output_path, "\n", sep = "")
    dir.create(output_path)
}




plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::coord_equal() +
    ggplot2::geom_point(
        data = data |> dplyr::filter(dataset |> stringr::str_detect("foldseek")),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            size = score^3,
            color = dataset),
        alpha = 0.9,
        shape = 16) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(dataset ==  "pfam"),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            color = dataset),
        alpha = 0.9,
        shape = 16) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(
            dataset |> stringr::str_detect("Frame2Seq")),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            color = dataset),
        size = 0.8,
        alpha = 0.9,
        shape = 17) +
    ggplot2::geom_point(
        data = data |>
            dplyr::filter(!is.na(seq_label)),            
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            color = dataset),
        size = 1.5,
        alpha = 1,
        shape = 16) +    
    ggrepel::geom_text_repel(
        data = data |>
            dplyr::filter(!is.na(seq_label)),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            label = seq_label),
        force = 10,
        min.segment.length = 0) +
    ggplot2::scale_x_continuous("UMAP1") +
    ggplot2::scale_y_continuous("UMAP2") +
    ggplot2::scale_size_continuous(
        "FoldSeek P-value",
        breaks = c(220, 410, 908)^3,
        labels = c("1e-6", "1e-10", "1e-22")) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::ggtitle(
        label = "UMAP Projection of EMS2 embeddings for Natural and Designed POC3 and POC4 Sequences",
        subtitle = "PFAM PF10448, Cu884 FoldSeek Hits, Cu884 Monomer Frame2Seq Design")
    


ggplot2::ggsave(
    filename = paste0(output_path, "/umap_embedding_natural_designed_20240521.pdf"),
    plot = plot,
    width = 12,
    height = 12,
    useDingbats = FALSE)


# zoom in
plot <- plot +
    ggplot2::coord_equal(
        xlim = c(-28, -24),
        ylim = c(4, 9.5))

ggplot2::ggsave(
    filename = paste0(output_path, "/umap_embedding_natural_designed_zoom_20240521.pdf"),
    plot = plot,
    width = 12,
    height = 12,
    useDingbats = FALSE)



######################################3

# fix interface


data <- dplyr::bind_rows(
    readr::read_tsv(
        "data/pfam/protein-matching-PF10448.tsv",
        show_col_types = FALSE) |>
    dplyr::transmute(
        dataset = "pfam",
        seq_id = Accession,
        prob = NA,
        score = NA,
        taxon = `Tax ID`,
        organism = `Tax Name`),
    readr::read_tsv(
        "intermediate_data/foldseek/Cauris_884_AF2_mmseqs_results_tt3_maeaoRTmWfE9qU_0Cs202D9oywnNW6T1XA/foldseek_metadata.tsv",
        show_col_types = FALSE) |>
    dplyr::transmute(
        dataset,
        seq_id = target,
        prob,
        score,
        taxon,
        organism),
    dplyr::bind_rows(
        readr::read_tsv(
            "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_5000_metadata.tsv",
            show_col_types = FALSE) |>
        dplyr::mutate(dataset = "Frame2Seq_T1.0")) |>
    dplyr::transmute(
        dataset,
        seq_id,
        prob = recovery_percent,
        score,
        taxon = NA_integer_,
        organism = NA_character_)) |>
    dplyr::mutate(
        seq_label = dplyr::case_when(
            seq_id == "A0A0L0P8N2" ~ "CuIRC25",
            seq_id == "AF-A0A510P858-F1-model_v4 Hypothetical_protein" ~ "Cu884",
            dataset == "foldseek_afdb50" & seq_id == "AF-A0A202G644-F1-model_v4 Uncharacterized protein" ~ "ClPOC4",
            dataset == "foldseek_afdb50" & seq_id == "AF-A0A4P9ZGT9-F1-model_v4 Uncharacterized protein" ~ "MbPOC4",
            seq_id == "A0A1D8PFG6" ~ "CaPOC3",
            dataset == "foldseek_cath50" & seq_id == "af_A0A1D8PTK4_4_109_3.30.230.100" ~ "CaPOC4",
            dataset == "foldseek_afdb_proteome" & seq_id == "AF-Q12245-F1-model_v4 Proteasome chaperone 4" ~ "ScPOC4",
            dataset == "Frame2Seq_T1.0" & seq_id == "seq_10" ~ "Design10",
            TRUE ~ NA)) |>
    dplyr::bind_cols(
        arrow::read_parquet(
            "intermediate_data/pfam,foldseek,frame2seq_dimera=1,b=1.8,n_neighbors=100,metric=cosine_20240522/umap_embedding.parquet"))

#data <- data |>
#    dplyr::mutate(
#        dist_884 = sqrt((UMAP_1 - -24.82538)^2 + (UMAP_2 -  6.462263)^2))




########################

output_path <- "product/UMAP_embedding/pfam,foldseek,frame2seq_dimera=1,b=1.8,n_neighbors=100,metric=cosine_20240522/"

if (!dir.exists(output_path)) {
    cat("creating dir ", output_path, "\n", sep = "")
    dir.create(output_path)
}




plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::coord_equal() +
    ggplot2::geom_point(
        data = data |> dplyr::filter(dataset |> stringr::str_detect("foldseek")),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            size = score^3,
            color = dataset),
        alpha = 0.9,
        shape = 16) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(dataset ==  "pfam"),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            color = dataset),
        alpha = 0.9,
        shape = 16) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(
            dataset |> stringr::str_detect("Frame2Seq")),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            color = dataset),
        size = 0.8,
        alpha = 0.9,
        shape = 17) +
    ggplot2::geom_point(
        data = data |>
            dplyr::filter(!is.na(seq_label)),            
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            color = dataset),
        size = 1.5,
        alpha = 1,
        shape = 16) +    
    ggrepel::geom_text_repel(
        data = data |>
            dplyr::filter(!is.na(seq_label)),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            label = seq_label),
        force = 10,
        min.segment.length = 0) +
    ggplot2::scale_x_continuous("UMAP1") +
    ggplot2::scale_y_continuous("UMAP2") +
    ggplot2::scale_size_continuous(
        "FoldSeek P-value",
        breaks = c(220, 410, 908)^3,
        labels = c("1e-6", "1e-10", "1e-22")) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::ggtitle(
        label = "UMAP Projection of EMS2 embeddings for Natural and Designed POC3 and POC4 Sequences",
        subtitle = "PFAM PF10448, Cu884 FoldSeek Hits, Cu884 Monomer Frame2Seq Design")
    


ggplot2::ggsave(
    filename = paste0(output_path, "/umap_embedding_natural_designed_20240522.pdf"),
    plot = plot,
    width = 12,
    height = 12,
    useDingbats = FALSE)


# zoom in
plot <- plot +
    ggplot2::coord_equal(
        xlim = c(-28, -24),
        ylim = c(4, 9.5))

ggplot2::ggsave(
    filename = paste0(output_path, "/umap_embedding_natural_designed_zoom_20240521.pdf"),
    plot = plot,
    width = 12,
    height = 12,
    useDingbats = FALSE)



##############################
######################################3

# fix conserved

data_pfam <- readr::read_tsv(
    "data/pfam/protein-matching-PF10448.tsv",
    show_col_types = FALSE) |>
    dplyr::transmute(
        dataset = "pfam",
        seq_id = Accession,
        prob = NA,
        score = NA,
        taxon = `Tax ID`,
        organism = `Tax Name`)

data_foldseek <- readr::read_tsv(
    "intermediate_data/foldseek/Cauris_884_AF2_mmseqs_results_tt3_maeaoRTmWfE9qU_0Cs202D9oywnNW6T1XA/foldseek_metadata.tsv",
    show_col_types = FALSE) |>
    dplyr::transmute(
        dataset,
        seq_id = target,
        prob,
        score,
        taxon,
        organism)


data_frame2seek <- dplyr::bind_rows(
    readr::read_tsv(
        "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_5000_metadata.tsv",
            show_col_types = FALSE) |>
        dplyr::mutate(dataset = "Frame2Seq_fI_T1.0"),
    readr::read_tsv(
        "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixconserved_T1.0_50000_metadata.tsv",
        show_col_types = FALSE) |>
        dplyr::mutate(dataset = "Frame2Seq_fC_T1.0")) |>
    dplyr::transmute(
        dataset,
        seq_id,
        prob = recovery_percent,
        score,
        taxon = NA_integer_,
        organism = NA_character_)

data <- dplyr::bind_rows(
    data_pfam,
    data_foldseek,
    data_frame2seek)

# add some seq_label
data <- data |>
    dplyr::mutate(
        seq_label = dplyr::case_when(
            seq_id == "A0A0L0P8N2" ~ "CuIRC25",
            seq_id == "AF-A0A510P858-F1-model_v4 Hypothetical_protein" ~ "Cu884",
            dataset == "foldseek_afdb50" & seq_id == "AF-A0A202G644-F1-model_v4 Uncharacterized protein" ~ "ClPOC4",
            dataset == "foldseek_afdb50" & seq_id == "AF-A0A4P9ZGT9-F1-model_v4 Uncharacterized protein" ~ "MbPOC4",
            seq_id == "A0A1D8PFG6" ~ "CaPOC3",
            dataset == "foldseek_cath50" & seq_id == "af_A0A1D8PTK4_4_109_3.30.230.100" ~ "CaPOC4",
            dataset == "foldseek_afdb_proteome" & seq_id == "AF-Q12245-F1-model_v4 Proteasome chaperone 4" ~ "ScPOC4",
            dataset == "Frame2Seq_T1.0" & seq_id == "seq_10" ~ "Design10",
            TRUE ~ NA))

# add umap coords
data <- data |>
    dplyr::bind_cols(
        arrow::read_parquet(
            "intermediate_data/pfam,foldseek,frame2seq_dimer,fixinterface,fixconserved,a=1,b=1.8,n_neighbors=100,metric=cosine_20240522/umap_embedding.parquet"))

#data <- data |>
#    dplyr::mutate(
#        dist_884 = sqrt((UMAP_1 - -24.82538)^2 + (UMAP_2 -  6.462263)^2))




output_path <- "product/UMAP_embedding/pfam,foldseek,frame2seq_dimer,fixinterface,fixconserved,a=1,b=1.8,n_neighbors=100,metric=cosine_20240522"

if (!dir.exists(output_path)) {
    cat("creating dir ", output_path, "\n", sep = "")
    dir.create(output_path)
}




plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::coord_equal() +
    ggplot2::geom_point(
        data = data |> dplyr::filter(dataset |> stringr::str_detect("foldseek")),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
#            size = score^3,
            color = dataset),
        alpha = 0.9,
        shape = 16) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(dataset ==  "pfam"),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            color = dataset),
        alpha = 0.9,
        shape = 16) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(
            dataset |> stringr::str_detect("Frame2Seq")),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            color = dataset,
            size = prob),
        size = 0.8,
        alpha = 0.9,
        shape = 17) +
    ggplot2::geom_point(
        data = data |>
            dplyr::filter(!is.na(seq_label)),            
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            color = dataset),
        size = 1.5,
        alpha = 1,
        shape = 16) +    
    ggrepel::geom_text_repel(
        data = data |>
            dplyr::filter(!is.na(seq_label)),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            label = seq_label),
        force = 10,
        min.segment.length = 0) +
    ggplot2::scale_x_continuous("UMAP1") +
    ggplot2::scale_y_continuous("UMAP2") +
    ggplot2::scale_size_continuous(
        "Frame2Seq Prob") +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::ggtitle(
        label = "UMAP Projection of EMS2 embeddings for Natural and Designed POC3 and POC4 Sequences",
        subtitle = "PFAM PF10448, Cu884 FoldSeek Hits, Cu884 Monomer Frame2Seq Design")
    


ggplot2::ggsave(
    filename = paste0(output_path, "/umap_embedding_natural_designed_20240522.pdf"),
    plot = plot,
    width = 12,
    height = 12,
    useDingbats = FALSE)


# zoom in
plot <- plot +
    ggplot2::coord_equal(
        xlim = c(22, 24.5),
        ylim = c(-4, -2.5))

ggplot2::ggsave(
    filename = paste0(output_path, "/umap_embedding_natural_designed_zoom_20240521.pdf"),
    plot = plot,
    width = 12,
    height = 12,
    useDingbats = FALSE)



##########################################3


# fix conserved

data_pfam <- readr::read_tsv(
    "data/pfam/protein-matching-PF10448.tsv",
    show_col_types = FALSE) |>
    dplyr::transmute(
        dataset = "pfam",
        seq_id = Accession,
        prob = NA,
        score = NA,
        taxon = `Tax ID`,
        organism = `Tax Name`)

data_foldseek <- readr::read_tsv(
    "intermediate_data/foldseek/Cauris_884_AF2_mmseqs_results_tt3_maeaoRTmWfE9qU_0Cs202D9oywnNW6T1XA/foldseek_metadata.tsv",
    show_col_types = FALSE) |>
    dplyr::transmute(
        dataset,
        seq_id = target,
        prob,
        score,
        taxon,
        organism)


data_frame2seek <- dplyr::bind_rows(
    readr::read_tsv(
        "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_5000_metadata.tsv",
            show_col_types = FALSE) |>
        dplyr::mutate(dataset = "Frame2Seq_fI_T1.0"),
    readr::read_tsv(
        "intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixconserved_T1.0_50000_metadata.tsv",
        show_col_types = FALSE) |>
        dplyr::mutate(dataset = "Frame2Seq_fC_T1.0")) |>
    dplyr::transmute(
        dataset,
        seq_id,
        prob = recovery_percent,
        score,
        taxon = NA_integer_,
        organism = NA_character_)

data_ProteinMPNN <- readr::read_tsv(
        "intermediate_data/ProteinMPNN/884_T0.1_5000/seqs/Cauris_AF-A0A2H1A4M0-F1-model_v4_metadata.tsv",
        show_col_types = FALSE) |>
    dplyr::transmute(
        dataset = "ProteinMPNN_T0.1",
        seq_id = sample |> as.character(),
        prob = seq_recovery,
        score,
        taxon = NA_integer_,
        organism = NA_character_)
    


data <- dplyr::bind_rows(
    data_pfam,
    data_foldseek,
    data_frame2seek,
    data_ProteinMPNN)

designs <- readr::read_tsv(
    paste0(
        "product/frame2seq_design/",
        "designs_fixinterface_T1.0_af2_pldd2_gt96_score_lt0.6_top30_20241022/",
        "summary.tsv"),
    show_col_types = FALSE) |>
    dplyr::filter(!is.na(`Selected for testing`))


# add some seq_label
data <- data |>
    dplyr::mutate(
        seq_label = dplyr::case_when(
            seq_id == "A0A0L0P8N2" ~ "CuIRC25",
            seq_id == "AF-A0A510P858-F1-model_v4 Hypothetical_protein" ~ "Cu884",
            dataset == "foldseek_afdb50" & seq_id == "AF-A0A202G644-F1-model_v4 Uncharacterized protein" ~ "ClPOC4",
            dataset == "foldseek_afdb50" & seq_id == "AF-A0A4P9ZGT9-F1-model_v4 Uncharacterized protein" ~ "MbPOC4",
            seq_id == "A0A1D8PFG6" ~ "CaPOC3",
            dataset == "foldseek_cath50" & seq_id == "af_A0A1D8PTK4_4_109_3.30.230.100" ~ "CaPOC4",
            dataset == "foldseek_afdb_proteome" & seq_id == "AF-Q12245-F1-model_v4 Proteasome chaperone 4" ~ "ScPOC4",
            dataset == "Frame2Seq_fI_T1.0" & (seq_id %in% designs$seq_id) ~ "Design",
            TRUE ~ NA))

# add umap coords
data <- data |>
    dplyr::bind_cols(
        arrow::read_parquet(
            "intermediate_data/pfam,foldseek,frame2seq_dimer,fixinterface,fixconserved,ProteinMPNN,a=1,b=1.8,n_neighbors=100,metric=cosine_20240719/umap_embedding.parquet"))

#data <- data |>
#    dplyr::mutate(
#        dist_884 = sqrt((UMAP_1 - -24.82538)^2 + (UMAP_2 -  6.462263)^2))




output_path <- "product/UMAP_embedding/pfam,foldseek,frame2seq_dimer,fixinterface,fixconserved,ProteinMPNN,a=1,b=1.8,n_neighbors=100,metric=cosine_202401104"

if (!dir.exists(output_path)) {
    cat("creating dir ", output_path, "\n", sep = "")
    dir.create(output_path)
}




plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::coord_equal() +
    ggplot2::geom_point(
        data = data |> dplyr::filter(dataset |> stringr::str_detect("foldseek")),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
#            size = score^3,
            color = dataset),
        alpha = 0.9,
        shape = 16) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(dataset ==  "pfam"),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            color = dataset),
        alpha = 0.9,
        shape = 16) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(
            dataset |> stringr::str_detect("Frame2Seq")),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            color = dataset,
            size = prob),
        size = 0.8,
        alpha = 0.9,
        shape = 17) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(
            dataset |> stringr::str_detect("ProteinMPNN")),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            color = dataset,
            size = prob),
        size = 0.8,
        alpha = 0.9,
        shape = 17) +
    ggplot2::geom_point(
        data = data |>
            dplyr::filter(!is.na(seq_label)),            
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            color = dataset),
        size = 1.5,
        alpha = 1,
        shape = 16) +    
    ggrepel::geom_text_repel(
        data = data |>
            dplyr::filter(!is.na(seq_label)),
        mapping = ggplot2::aes(
            x = UMAP_1,
            y = UMAP_2,
            label = seq_label),
        force = 10,
        min.segment.length = 0) +
    ggplot2::scale_x_continuous("UMAP1") +
    ggplot2::scale_y_continuous("UMAP2") +
    ggplot2::scale_size_continuous(
        "Frame2Seq Prob") +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::ggtitle(
        label = "UMAP Projection of EMS2 embeddings for Natural and Designed POC3 and POC4 Sequences",
        subtitle = "PFAM PF10448, Cu884 FoldSeek Hits, Cu884 Monomer Frame2Seq Design, ProteinMPNN")
    


ggplot2::ggsave(
    filename = paste0(output_path, "/umap_embedding_natural_designed_20241105.pdf"),
    plot = plot,
    width = 12,
    height = 12,
    useDingbats = FALSE)


# zoom in
plot <- plot +
    ggplot2::coord_equal(
        xlim = c(22, 24.5),
        ylim = c(-4, -2.5))

ggplot2::ggsave(
    filename = paste0(output_path, "/umap_embedding_natural_designed_zoom_20241105.pdf"),
    plot = plot,
    width = 12,
    height = 12,
    useDingbats = FALSE)



