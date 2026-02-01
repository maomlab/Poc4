# to get the path to the umap_embed program
source("parameters.R")



date_code <- "20240507"
tag <- "ESM2_embeddings"
expand.grid(
    a = c(1),
    b = c(1.8),
    n_neighbors = c(100),
    metric = c("euclidean")) |>
    dplyr::rowwise() |>
    dplyr::do({
        params <- .
        alt_tag <- paste0(
            "pfam,foldseek",
            "a=", params$a, ",",
            "b=", params$b, ",",
            "n_neighbors=", params$n_neighbors, ",",
            "metric=", params$metric, "_",
            "APDP_",
            date_code)
        feature_columns <- paste0(
            "intermediate_data/", tag, "/feature_columns.tsv ")
        cat("Computing embedding for: ", alt_tag, "\n", sep = "")
        command <- paste0(
            parameters$embed_umap_program, " ",
            "--dataset ",
            "intermediate_data/ESM2_embeddings/pfam_PF10448.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_afdb50.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_afdb-proteome.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_afdb-swissprot.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_cath50.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_pdb100.parquet ",
#            "intermediate_data/ESM2_embeddings/frame2seq_Cauris_884_AF2_T1.0_5000.parquet ",
            "--pca_n_components 100 ",
            "--feature_columns ", feature_columns, " ",
            "--tag ", alt_tag, " ",
            "--umap_a ", params$a, " ",
            "--umap_b ", params$b, " ",
            "--umap_n_neighbors ", params$n_neighbors, " ",
            "--umap_metric ", params$metric, " ",
            "--umap_n_epochs 2000 ",
            "--verbose",
        sep = "")
        cat(command, sep = "")
        system(command)
        data.frame()
    })



date_code <- "20240506"
tag <- "ESM2_embeddings"
expand.grid(
    a = c(1),
    b = c(1.8),
    n_neighbors = c(100),
    metric = c("cosine")) |>
    dplyr::rowwise() |>
    dplyr::do({
        params <- .
        alt_tag <- paste0(
            "pfam,foldseek,frame2seq_monomer",
            "a=", params$a, ",",
            "b=", params$b, ",",
            "n_neighbors=", params$n_neighbors, ",",
            "metric=", params$metric, "_",
            "APDP_",
            date_code)
        feature_columns <- paste0(
            "intermediate_data/", tag, "/feature_columns.tsv ")
        cat("Computing embedding for: ", alt_tag, "\n", sep = "")
        command <- paste0(
            parameters$embed_umap_program, " ",
            "--dataset ",
            "intermediate_data/ESM2_embeddings/pfam_PF10448.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_afdb50.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_afdb-proteome.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_afdb-swissprot.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_cath50.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_pdb100.parquet ",
            "intermediate_data/ESM2_embeddings/frame2seq_Cauris_884_AF2_T1.0_5000.parquet ",
            "--pca_n_components 100 ",
            "--feature_columns ", feature_columns, " ",
            "--tag ", alt_tag, " ",
            "--umap_a ", params$a, " ",
            "--umap_b ", params$b, " ",
            "--umap_n_neighbors ", params$n_neighbors, " ",
            "--umap_metric ", params$metric, " ",
            "--umap_n_epochs 2000 ",
            "--verbose",
        sep = "")
        cat(command, sep = "")
        system(command)
        data.frame()
    })


date_code <- "20240521"
tag <- "ESM2_embeddings"
expand.grid(
    a = c(1),
    b = c(1.8),
    n_neighbors = c(100),
    metric = c("cosine")) |>
    dplyr::rowwise() |>
    dplyr::do({
        params <- .
        alt_tag <- paste0(
            "pfam,foldseek,frame2seq_monomer",
            "a=", params$a, ",",
            "b=", params$b, ",",
            "n_neighbors=", params$n_neighbors, ",",
            "metric=", params$metric, "_",
            date_code)
        feature_columns <- paste0(
            "intermediate_data/", tag, "/feature_columns.tsv ")
        cat("Computing embedding for: ", alt_tag, "\n", sep = "")
        command <- paste0(
            parameters$embed_umap_program, " ",
            "--dataset ",
            "intermediate_data/ESM2_embeddings/pfam_PF10448.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_afdb50.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_afdb-proteome.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_afdb-swissprot.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_cath50.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_pdb100.parquet ",
#            "intermediate_data/ESM2_embeddings/frame2seq_Cauris_884_AF2_T1.0_5000.parquet ",
            "intermediate_data/ESM2_embeddings/frame2seq_Cauris_884_irc25_AF3_T0.01_5000.parquet ",
            "intermediate_data/ESM2_embeddings/frame2seq_Cauris_884_irc25_AF3_T0.1_5000.parquet ",
            "intermediate_data/ESM2_embeddings/frame2seq_Cauris_884_irc25_AF3_T1.0_5000.parquet ",
            "intermediate_data/ESM2_embeddings/frame2seq_Cauris_884_irc25_AF3_T10_5000.parquet ",
            "--pca_n_components 100 ",
            "--feature_columns ", feature_columns, " ",
            "--tag ", alt_tag, " ",
            "--umap_a ", params$a, " ",
            "--umap_b ", params$b, " ",
            "--umap_n_neighbors ", params$n_neighbors, " ",
            "--umap_metric ", params$metric, " ",
            "--umap_n_epochs 2000 ",
            "--verbose",
        sep = "")
        cat(command, sep = "")
        system(command)
        data.frame()
    })



##################3
date_code <- "20240522"
tag <- "ESM2_embeddings"
expand.grid(
    a = c(1),
    b = c(1.8),
    n_neighbors = c(100),
    metric = c("cosine")) |>
    dplyr::rowwise() |>
    dplyr::do({
        params <- .
        alt_tag <- paste0(
            "pfam,foldseek,frame2seq_dimer",
            "a=", params$a, ",",
            "b=", params$b, ",",
            "n_neighbors=", params$n_neighbors, ",",
            "metric=", params$metric, "_",
            date_code)
        feature_columns <- paste0(
            "intermediate_data/", tag, "/feature_columns.tsv ")
        cat("Computing embedding for: ", alt_tag, "\n", sep = "")
        command <- paste0(
            parameters$embed_umap_program, " ",
            "--dataset ",
            "intermediate_data/ESM2_embeddings/pfam_PF10448.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_afdb50.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_afdb-proteome.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_afdb-swissprot.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_cath50.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_pdb100.parquet ",
            "intermediate_data/ESM2_embeddings/frame2seq_Cauris_884_irc25_AF3_fixinterface_T1.0_5000.parquet ",
            "--pca_n_components 100 ",
            "--feature_columns ", feature_columns, " ",
            "--tag ", alt_tag, " ",
            "--umap_a ", params$a, " ",
            "--umap_b ", params$b, " ",
            "--umap_n_neighbors ", params$n_neighbors, " ",
            "--umap_metric ", params$metric, " ",
            "--umap_n_epochs 2000 ",
            "--verbose",
        sep = "")
        cat(command, sep = "")
        system(command)
        data.frame()
    })


##################3
date_code <- "20240522"
tag <- "ESM2_embeddings"
expand.grid(
    a = c(1),
    b = c(1.8),
    n_neighbors = c(100),
    metric = c("cosine")) |>
    dplyr::rowwise() |>
    dplyr::do({
        params <- .
        alt_tag <- paste0(
            "pfam,foldseek,frame2seq_dimer,fixinterface,fixconserved,",
            "a=", params$a, ",",
            "b=", params$b, ",",
            "n_neighbors=", params$n_neighbors, ",",
            "metric=", params$metric, "_",
            date_code)
        feature_columns <- paste0(
            "intermediate_data/", tag, "/feature_columns.tsv ")
        cat("Computing embedding for: ", alt_tag, "\n", sep = "")
        command <- paste0(
            parameters$embed_umap_program, " ",
            "--dataset ",
            "intermediate_data/ESM2_embeddings/pfam_PF10448.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_afdb50.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_afdb-proteome.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_afdb-swissprot.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_cath50.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_pdb100.parquet ",
            "intermediate_data/ESM2_embeddings/frame2seq_Cauris_884_irc25_AF3_fixinterface_T1.0_5000.parquet ",
            "intermediate_data/ESM2_embeddings/frame2seq_Cauris_884_irc25_AF3_fixconserved_T1.0_50000.parquet ",
            "--pca_n_components 100 ",
            "--feature_columns ", feature_columns, " ",
            "--tag ", alt_tag, " ",
            "--umap_a ", params$a, " ",
            "--umap_b ", params$b, " ",
            "--umap_n_neighbors ", params$n_neighbors, " ",
            "--umap_metric ", params$metric, " ",
            "--umap_n_epochs 2000 ",
            "--verbose",
        sep = "")
        cat(command, sep = "")
        system(command)
        data.frame()
    })


##########################################
date_code <- "20240719"
tag <- "ESM2_embeddings"
expand.grid(
    a = c(1),
    b = c(1.8),
    n_neighbors = c(100),
    metric = c("cosine")) |>
    dplyr::rowwise() |>
    dplyr::do({
        params <- .
        alt_tag <- paste0(
            "pfam,foldseek,frame2seq_dimer,fixinterface,fixconserved,",
            "a=", params$a, ",",
            "b=", params$b, ",",
            "n_neighbors=", params$n_neighbors, ",",
            "metric=", params$metric, "_",
            date_code)
        feature_columns <- paste0(
            "intermediate_data/", tag, "/feature_columns.tsv ")
        cat("Computing embedding for: ", alt_tag, "\n", sep = "")
        command <- paste0(
            parameters$embed_umap_program, " ",
            "--dataset ",
            "intermediate_data/ESM2_embeddings/pfam_PF10448.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_afdb50.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_afdb-proteome.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_afdb-swissprot.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_cath50.parquet ",
            "intermediate_data/ESM2_embeddings/foldseek/Cauris_2044_AF2/alis_pdb100.parquet ",
            "intermediate_data/ESM2_embeddings/frame2seq_Cauris_884_irc25_AF3_fixinterface_T1.0_5000.parquet ",
            "intermediate_data/ESM2_embeddings/frame2seq_Cauris_884_irc25_AF3_fixconserved_T1.0_50000.parquet ",
            "intermediate_data/ESM2_embeddings/ProteinMPNN_Cauris_884_AF2_T0.1_5000.parquet ",
            "--pca_n_components 100 ",
            "--feature_columns ", feature_columns, " ",
            "--tag ", alt_tag, " ",
            "--umap_a ", params$a, " ",
            "--umap_b ", params$b, " ",
            "--umap_n_neighbors ", params$n_neighbors, " ",
            "--umap_metric ", params$metric, " ",
            "--umap_n_epochs 2000 ",
            "--verbose",
        sep = "")
        cat(command, sep = "")
        system(command)
        data.frame()
    })
