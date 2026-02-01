#############
# Set paths #
#############

intermediate_path <- "intermediate_data/dockq_analysis/hinge_model"
if (!dir.exists(intermediate_path)) {
    cat("Creating directory '", intermediate_path, "'\n", sep = "")
    dir.create(intermediate_path, recursive = TRUE)
}


product_path <- "product/dockq_analysis/hinge_model"
if (!dir.exists(product_path)) {
    cat("Creating directory '", product_path, "'\n", sep = "")
    dir.create(product_path, recursive = TRUE)
}

date_code <- "20260113"

#### Load Data ###

load("intermediate_data/growth_data_summary.Rdata")

growth_data_summary <- growth_data_summary |>
    dplyr::select(
        structure_id = Strain,
        recovery_score = score,
        recovery_class = class) |>
    dplyr::mutate(
        recovery_class = recovery_class |>
            stringr::str_replace("Parial", "Partial"))


load("intermediate_data/AlphaFold3/afs_dockq.Rdata")
load("intermediate_data/AlphaFold3/afs_dockq_test.Rdata")

hinge_geometries <- dplyr::bind_rows(
    readr::read_tsv(
        file = paste0(intermediate_path, "/design_hinge_features_20260121.tsv"),
        show_col_types = FALSE),
    readr::read_tsv(
        file = paste0(intermediate_path, "/design_hinge_features_test_20260121.tsv"),
        show_col_types = FALSE)) |>
    dplyr::rename(hinge_distance = distance)


plip_train_predicted_cv <- readr::read_tsv(
    file = paste0("product/dockq_analysis/structure_model/train_predicted_cv_20260118.tsv"),
    show_col_types = FALSE)

plip_test_predicted <- readr::read_tsv(
    file = paste0("product/dockq_analysis/structure_model/test_predicted_20260118.tsv"),
    show_col_types = FALSE)

hinge_geometries <- hinge_geometries |>
    dplyr::inner_join(
        dplyr::bind_rows(
            afs_dockq,
            afs_dockq_test),
        by = "cif_fname") |>
    dplyr::left_join(
        dplyr::bind_rows(
            plip_train_predicted_cv,
            plip_test_predicted) |>
            dplyr::select(
                cif_fname,
                plip_predict = predict),
        by = "cif_fname")

plot <- ggplot2::ggplot(data = hinge_geometries) +
    ggplot2::theme_bw() +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = sqrt(hinge_distance),
            y = pDockQ2_Poc4_Pup2,
            color = plip_predict),
        alpha = 0.8,
        shape = 16) +
    ggplot2::ggtitle(
        "DockQ2 vs Hinge Distance") +
    viridis::scale_color_viridis(
        option = "C",
        begin = 0.1,
        end = 0.9,
        guide = ggplot2::guide_colourbar(barwidth = 10, barheight = 1)) +
    ggplot2::scale_x_continuous(
        breaks = sqrt(c(20, 25, 30, 35, 40)),
        labels = c("20A", "25A", "30A", "35A", "40A")) +
    ggplot2::labs(
        x = "Hinge Distance",
        y = "pDockQ2 Poc4 vs. Pup2",
        color = "PLiP Model Prediction") +
    ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(
    filename = paste0(product_path, "/DockQ2_vs_hinge_distance_", date_code, ".pdf"),
    width = 5,
    height = 5,
    useDingbats = FALSE)
