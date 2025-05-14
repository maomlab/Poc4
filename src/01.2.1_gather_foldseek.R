
#########################
# parse foldseek output #
#########################

foldseek_M8_path <- "data/foldseek/Cauris_884_AF2_mmseqs_results_tt3_maeaoRTmWfE9qU_0Cs202D9oywnNW6T1XA"
foldseek_column_names <- c(
    "query",
    "target",
    "qca",
    "tca",
    "alntmscore",
    "qtmscore",
    "ttmscore",
    "u",
    "t",
    "lddt",
    "lddtfull",
    "prob",
    "score",
    paste0("X", 13:18),
    "taxon",
    "organism")

foldseek_metadata <- dplyr::bind_rows(
    readr::read_tsv(
        paste0(foldseek_m8_path, "/alis_afdb50.m8"),
        col_names = foldseek_column_names,
        show_col_types = FALSE) |>
    dplyr::mutate(dataset = "foldseek_afdb50", .before = 1),
    readr::read_tsv(
        paste0(foldseek_m8_path, "/alis_afdb-proteome.m8"),
        col_names = foldseek_column_names,
        show_col_types = FALSE) |>
    dplyr::mutate(dataset = "foldseek_afdb_proteome", .before = 1),
    readr::read_tsv(
        paste0(foldseek_m8_path, "/alis_afdb-swissprot.m8"),
        col_names = foldseek_column_names,
        show_col_types = FALSE) |>
    dplyr::mutate(dataset = "foldseek_afdb_swissprot", .before = 1),
    readr::read_tsv(
        paste0(foldseek_m8_path, "/alis_cath50.m8"),
        col_names = foldseek_column_names,
        show_col_types = FALSE) |>
    dplyr::mutate(dataset = "foldseek_cath50", .before = 1),
    readr::read_tsv(
        paste0(foldseek_m8_path, "/alis_pdb100.m8"),
        col_names = foldseek_column_names,
        show_col_types = FALSE) |>
    dplyr::mutate(dataset = "foldseek_pdb100", .before = 1))


#   dataset                     n
#   <chr>                   <int>
# 1 foldseek_afdb50          1000
# 2 foldseek_afdb_proteome     81
# 3 foldseek_afdb_swissprot   103
# 4 foldseek_cath50            90
# 5 foldseek_pdb100            56

foldseek_metadata |>
    readr::write_tsv("intermediate_data/foldseek/Cauris_884_AF2_mmseqs_results_tt3_maeaoRTmWfE9qU_0Cs202D9oywnNW6T1XA/foldseek_metadata.tsv")
