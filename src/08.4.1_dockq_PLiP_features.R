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

date_code <- "20260103"


#############
# Load Data #
#############

load("intermediate_data/growth_data_summary.Rdata")
load("intermediate_data/AlphaFold3/afs_dockq.Rdata")

afs_dockq_features_input <- afs_dockq |>
    dplyr::transmute(
        design_id = structure_id,
        complex_cif = cif_fname,
        target_chain = "D",
        design_chain = "K",
        DockQ2 = pDockQ2_Poc4_Pup2) |>
    dplyr::arrange(dplyr::desc(DockQ2))

afs_dockq_features_input |>
    readr::write_csv(
        file = paste0(
            intermediate_path,
            "/afs_dockq_features_input.csv"))

cmd <- paste0(
    python_path, " src/interaction_features_extract.py ",
    "--input '", intermediate_path, "/afs_dockq_features_input.csv' ",
    "--output '", intermediate_path, "/afs_dockq_features.csv' ",
    "--n_jobs 1")
cat(cmd, "\n", sep = "")
system(cmd)

afs_dockq_features <- readr::read_csv(
    file = paste0(intermediate_path, "/afs_dockq_features.csv"),
    show_col_types = FALSE)

afs_dockq_features_ml <- afs_dockq_features |>
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

afs_dockq_features_ml |>
    readr::write_csv(
        file = paste0(
            intermediate_path, "/afs_dockq_features_ml.csv"))

# Write out the features table to facilite inspecting structures and features
afs_dockq_features_ml_inspect <- afs_dockq_features_ml |>
    dplyr::group_by(design_id) |>
    dplyr::arrange(dplyr::desc(DockQ2)) |>
    dplyr::mutate(
        rank = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::arrange(rank, design_id) |>
    dplyr::left_join(
        growth_data_summary |>
        dplyr::select(
            design_id = Strain,
            recovery_score = score,
            recovery_class = class),
        by = "design_id") |>
    dplyr::filter(
        !is.na(recovery_class)) |>
    dplyr::arrange(rank, design_id)

afs_tested_features |>
    readr::write_tsv(
        file = paste0(
            product_path,
            "/afs_dockq_feastures_ml_inspect_", date_code, ".tsv"))
