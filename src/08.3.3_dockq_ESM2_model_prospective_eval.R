#############
# Set paths #
#############

intermediate_path <- "intermediate_data/dockq_analysis/sequence_model"
if (!dir.exists(intermediate_path)) {
    cat("Creating directory '", intermediate_path, "'\n", sep = "")
    dir.create(intermediate_path, recursive = TRUE)
}


product_path <- "product/dockq_analysis/sequence_model"
if (!dir.exists(product_path)) {
    cat("Creating directory '", product_path, "'\n", sep = "")
    dir.create(product_path, recursive = TRUE)
}

date_code <- "20260127"

#### Load Data ###

train_predicted_cv <- readr::read_tsv(
    file = paste0(product_path, "/train_predicted_cv_20260117.tsv"),
    show_col_types = FALSE)

load("intermediate_data/AlphaFold3/afs_dockq_test.Rdata")
test_predicted <- readr::read_tsv(
    file = paste0(
        product_path,
        "/designs_fixinterface_T1.0_af2_v2_pldd2_gt96_score_lt0.6_20251103_predicted_20260113.tsv"),
    show_col_types = FALSE)

test_predicted <- afs_dockq_test |>
    dplyr::left_join(
        test_predicted |>
        dplyr::rename(full_structure_id = structure_id),
        by = c("Poc4_sequence" = "fasta"))


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
        "ML DockQ2 Prediction",
        subtitle = paste0("based on Poc4 ESM Embedding (R2 = ", r_squared_label, ")")) +
    ggplot2::labs(
        x = "h2o AutoML Ensemble Predictor Score",
        y = "pDockQ2 between Poc4 vs. Pup2 (max)",
        color = "Dataset")

ggplot2::ggsave(
    filename = "product/dockq_analysis/sequence_model/afs_esm_ml_prediction_and_test_20260117.pdf",
    width = 5,
    height = 4,
    useDingbats = FALSE)


#################
median_scores <- all_predicted |>
    dplyr::group_by(dataset_id) |>
    dplyr::summarize(
        pDockQ2_Poc4_Pup2 = median(pDockQ2_Poc4_Pup2),
        .groups = "drop")

plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_vline(
        data = median_scores,
        mapping = ggplot2::aes(
            xintercept = pDockQ2_Poc4_Pup2)) +
    ggplot2::geom_density(
        data = all_predicted,
        mapping = ggplot2::aes(
            x = pDockQ2_Poc4_Pup2,
            group = dataset_id,
            fill = dataset_id)) +
    viridis::scale_fill_viridis(
        option = "inferno",
        begin = 0.2,
        end = 0.8,
        alpha = 0.6,
        discrete = TRUE) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::ggtitle(
        "ML DockQ2 Prediction",
        subtitle = "based on Poc4 ESM Embedding") +
    ggplot2::labs(
        x = "h2o AutoML Ensemble Predictor Score",
        y = "pDockQ2 between Poc4 vs. Pup2 (max)",
        fill = "Dataset")

ggplot2::ggsave(
    filename = paste0(
        product_path,
        "/afs_esm_ml_prediction_and_test_density_", date_code, ".pdf"),
    width = 5,
    height = 4,
    useDingbats = FALSE)


roc_obj <- pROC::roc(
  response  = all_predicted$dataset_id,
  predictor = all_predicted$pDockQ2_Poc4_Pup2,
  levels    = c("Prospective", "Cross Validated"),
  direction = ">")

pROC::auc(roc_obj)
pROC::ci.auc(roc_obj, method = "bootstrap", boot.n = 2000)
