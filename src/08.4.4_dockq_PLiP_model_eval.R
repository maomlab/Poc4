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

date_code <- "20260118"

#############
# Load Data #
#############

train_predicted_cv <- readr::read_tsv(
    file = paste0(product_path, "/train_predicted_cv_20260118.tsv"),
    show_col_types = FALSE)

test_predicted <- readr::read_tsv(
    file = paste0(product_path, "/test_predicted_20260118.tsv"),
    show_col_types = FALSE)



#################################
# Compare predicted vs observed #
#################################

all_predicted <- dplyr::bind_rows(
    train_predicted_cv,
    test_predicted) |>
    dplyr::mutate(
        dataset_id = dplyr::case_when(
            dataset_id == "AFS" ~ "Cross Validated",
            dataset_id == "AFS_test" ~ "Prospective"))


lm_fit <- lm(
    formula = pDockQ2_Poc4_Pup2 ~ predict,
    data = all_predicted)
r_squared_label <- summary(lm_fit)$r.squared |> signif(2)


plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_smooth(
        data = all_predicted,
        mapping = ggplot2::aes(
            x = predict,
            y = pDockQ2_Poc4_Pup2),
        formula = y ~ x,
        method = "lm") +
    ggplot2::geom_point(
        data = all_predicted,
        mapping = ggplot2::aes(
            x = predict,
            y = pDockQ2_Poc4_Pup2,
            color = dataset_id),
        alpha = 0.8,
        shape = 16) +
    viridis::scale_color_viridis(
        option = "inferno",
        begin = 0.2,
        end = 0.8,
        discrete = TRUE) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::ggtitle(
        "Structure ML DockQ2 Prediction",
        subtitle = paste0(
            "based on Poc4 PuP2 PLiP Interaction Features ",
            "(R2 = ", r_squared_label, ")")) +
    ggplot2::labs(
        x = "h2o AutoML Ensemble Predictor Score",
        y = "pDockQ2 between Poc4 vs. Pup2 (max)",
        color = "Dataset")

ggplot2::ggsave(
    filename =
        paste0(product_path, "/afs_dockq_prediction_accuracy_", date_code, ".pdf"),
    width = 5,
    height = 4,
    useDingbats = FALSE)
