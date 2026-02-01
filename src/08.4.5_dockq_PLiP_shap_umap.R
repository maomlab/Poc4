library(h2o)


#############
# Set paths #
#############

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

date_code <- "20260130"

h2o::h2o.init()

#### Load Data ####

load("intermediate_data/growth_data_summary.Rdata")

growth_data_summary <- growth_data_summary |>
    dplyr::select(
        design_id = Strain,
        recovery_score = score,
        recovery_class = class) |>
    dplyr::mutate(
        recovery_class = recovery_class |>
            stringr::str_replace("Parial", "Partial"))


#train_data_path <- "intermediate_data/dockq_analysis/structure_model/afs_dockq_features_ml.csv"
#load("intermediate_data/AlphaFold3/afs_dockq.Rdata")
#
#test_data_path <- "intermediate_data/dockq_analysis/structure_model/afs_dockq_test_features_ml.csv"
#load("intermediate_data/AlphaFold3/afs_dockq_test.Rdata")
#
#train_data <- h2o::h2o.importFile(
#    path = train_data_path,
#    destination_frame = "train_data")
#test_data <- h2o::h2o.importFile(
#    path = test_data_path,
#    destination_frame = "test_data")



combined_data_path <- "intermediate_data/dockq_analysis/structure_model/afs_dockq_combined_features_ml.csv"

combined_data <- h2o::h2o.importFile(
    path = combined_data_path,
    destination_frame = "combined_data")

combined_df <- readr::read_csv(
    file = combined_data_path,
    show_col_types = FALSE)


gbm <- h2o.upload_model(
    path = paste0(
        "intermediate_data/dockq_analysis/structure_model/",
        "combined_model/GBM_2_AutoML_2_20260129_150314"))


###############################################
## UMAP embedding of shap normalized features #
###############################################

contributions <- gbm |>
    h2o::h2o.predict_contributions(combined_data) |>
    as.data.frame()

tidyr::expand_grid(
    a = c(0.5),
    b = c(12)) |>
    dplyr::rowwise() |>
    dplyr::do({
        x <- .

        cat("Plotting umap with a=", x$a, " b=", x$b, "\n", sep = "")
        umap <- uwot::umap(
            X = contributions,
            a = x$a,
            b = x$b)

        umap <- umap |> as.data.frame()
        names(umap) <- c("UMAP1", "UMAP2")

        umap <- umap |>
            dplyr::bind_cols(
                combined_df)

        umap <- umap |>
            dplyr::left_join(
                growth_data_summary,
                by = "design_id")


        umap <- umap |>
            dplyr::group_by(design_id) |>
            dplyr::arrange(dplyr::desc(DockQ2)) |>
            dplyr::mutate(
                rank = dplyr::row_number()) |>
            dplyr::ungroup()

        # sort the points by what should be on top of other points
        umap <- umap |>
            dplyr::arrange(DockQ2)
        #umap <- umap |>
        #    dplyr::sample_n(size = dplyr::n(), replace = FALSE)

        plot <- ggplot2::ggplot(
            data = umap) +
            ggplot2::theme_bw() +
            ggplot2::geom_point(
                data = umap |> dplyr::filter(rank > 1),
                mapping = ggplot2::aes(
                    x = UMAP1,
                    y = UMAP2),
                color = "lightgray",
                alpha = 0.8,
                shape = 16) +
            ggplot2::geom_point(
                data = umap |> dplyr::filter(rank == 1),
                mapping = ggplot2::aes(
                    x = UMAP1,
                    y = UMAP2,
                    color = DockQ2),
                alpha = 0.8,
                shape = 16) +
           ggplot2::geom_point(
               data = umap |>
                   dplyr::filter(
                       rank == 1,
                       recovery_class == "Not Recovered"),
               mapping = ggplot2::aes(
                   x = UMAP1,
                   y = UMAP2),
               color = "yellow",
               alpha = 0.8,
               shape = 4) +
           ggplot2::geom_point(
               data = umap |> dplyr::filter(
                   rank == 1,
                   recovery_class == "Partial"),
               mapping = ggplot2::aes(
                   x = UMAP1,
                   y = UMAP2),
               color = "orange",
               alpha = 0.8,
               shape = 4) +
           ggplot2::geom_point(
               data = umap |> dplyr::filter(
                   rank == 1,
                   recovery_class == "Recovered"),
               mapping = ggplot2::aes(
                   x = UMAP1,
                   y = UMAP2),
               color = "red",
               alpha = 0.8,
               shape = 4) +
           ggrepel::geom_label_repel(
               data = umap |> dplyr::filter(
                   rank == 1,
                   !is.na(recovery_class)),
               mapping = ggplot2::aes(
                   x = UMAP1,
                   y = UMAP2,
                   label = design_id),
               size = .6) +
            ggplot2::coord_equal() +
            viridis::scale_color_viridis(
                "DockQ2",
                begin = 0.2,
                end = 0.8,
                option = "D") +
            ggplot2::theme(
                legend.position = "bottom") +
            ggplot2::ggtitle(
                "UMAP of SHAP normalized PLiP Features",
                subtitle = "Poc4 vs Pup2")

        ggplot2::ggsave(
            filename = paste0(
                product_path, "/shap_umap_",
                "a=", x$a, "_b=", x$b, "_", date_code, ".pdf"),
            plot = plot,
            width = 7,
            height = 7,
            useDingbats = FALSE)

        data.frame()
    })
