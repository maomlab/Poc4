#############
# Set paths #
#############

python_path <- "/nfs/turbo/umms-maom/maom/opt/miniforge3/envs/plip/bin/python"

intermediate_path <- "intermediate_data/dockq_analysis/structure_model"
if (!dir.exists(intermediate_path)) {
    cat("Creating directory '", intermediate_path, "'\n", sep = "")
    dir.create(intermediate_path, recursive = TRUE)
}


product_path <- "product/dockq_analysis/structure_model"
if (!dir.exists(product_path)) {
    cat("Creating directory '", product_path, "'\n", sep = "")
    dir.create(product_path, recursive = TRUE)
}

date_code <- "20260118"


#############
# Load Data #
#############

# load the training data features to make sure they match the feature types
afs_dockq_features_ml <- readr::read_csv(
    file = paste0(
        "intermediate_data/dockq_analysis/structure_model/afs_dockq_features_ml.csv"),
    show_col_types = FALSE)

load("intermediate_data/AlphaFold3/afs_dockq_test.Rdata")

afs_dockq_test_features_input <- afs_dockq_test |>
    dplyr::transmute(
        design_id = structure_id,
        complex_cif = cif_fname,
        target_chain = "D",
        design_chain = "K",
        DockQ2 = pDockQ2_Poc4_Pup2) |>
    dplyr::arrange(dplyr::desc(DockQ2))

afs_dockq_test_features_input |>
    readr::write_csv(
        file = paste0(
            intermediate_path,
            "/afs_dockq_test_features_input.csv"))

cmd <- paste0(
    python_path, " src/interaction_features_extract.py ",
    "--input '", intermediate_path, "/afs_dockq_test_features_input.csv' ",
    "--output '", intermediate_path, "/afs_dockq_test_features.csv' ",
    "--n_jobs 1")
cat(cmd, "\n", sep = "")
system(cmd)

afs_dockq_test_features <- readr::read_csv(
    file = paste0(intermediate_path, "/afs_dockq_test_features.csv"),
    show_col_types = FALSE)

afs_dockq_test_features_ml <- afs_dockq_test_features |>
    dplyr::transmute(
        design_id, cif_path, DockQ2,
        feature_id = paste(
            target_chain, target_resnr, target_restype, interaction_type,
            sep = "_")) |>
    dplyr::count(
        design_id, cif_path, DockQ2,
        feature_id) |>
    tidyr::pivot_wider(
        id_cols = c(design_id, cif_path, DockQ2),
        names_from = "feature_id",
        values_from = "n",
        values_fill = 0)

# filter the features to match the training data, filling with zeros where they are missing
feature_columns <- names(afs_dockq_features_ml)
missing_features <- setdiff(
    feature_columns,
    names(afs_dockq_test_features_ml))
missing_feature_values <- matrix(
    data = 0,
    nrow = nrow(afs_dockq_test_features_ml),
    ncol = length(missing_features),
    dimnames = list(NULL, missing_features))
afs_dockq_test_features_ml <- afs_dockq_test_features_ml |>
    dplyr::bind_cols(missing_feature_values) |>
    dplyr::select(tidyselect::all_of(feature_columns))

afs_dockq_test_features_ml |>
    readr::write_csv(
        file = paste0(
            intermediate_path, "/afs_dockq_test_features_ml.csv"))

afs_dockq_combined_features_ml <- dplyr::bind_rows(
    afs_dockq_features_ml,
    afs_dockq_test_features_ml)

afs_dockq_combined_features_ml |>
    readr::write_csv(
        file = paste0(
            intermediate_path, "/afs_dockq_combined_features_ml.csv"))


# Write out the features table to facilite inspecting structures and features
afs_dockq_test_features_ml_inspect <- afs_dockq_test_features_ml |>
    dplyr::group_by(design_id) |>
    dplyr::arrange(dplyr::desc(DockQ2)) |>
    dplyr::mutate(
        rank = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::arrange(rank, design_id)

afs_dockq_test_features_ml_inspect |>
    readr::write_tsv(
        file = paste0(
            product_path,
            "/afs_dockq_feastures_ml_inspect_", date_code, ".tsv"))
