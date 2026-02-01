library(patchwork)

#############
# Set paths #
#############

intermediate_path <- "intermediate_data/dockq_analysis"
if (!dir.exists(intermediate_path)) {
    cat("Creating directory '", intermediate_path, "'\n", sep = "")
    dir.create(intermediate_path, recursive = TRUE)
}


product_path <- "product/dockq_analysis"
if (!dir.exists(product_path)) {
    cat("Creating directory '", product_path, "'\n", sep = "")
    dir.create(product_path, recursive = TRUE)
}

date_code <- "20260124"

##################
#### Load Data ###
##################

load("intermediate_data/growth_data_summary.Rdata")

growth_data_summary <- growth_data_summary |>
    dplyr::select(
        structure_id = Strain,
        recovery_score = score,
        recovery_class = class) |>
    dplyr::mutate(
        recovery_class = recovery_class |>
            stringr::str_replace("Parial", "Partial"))

# DockQ2 Scores
load("intermediate_data/AlphaFold3/afs_dockq.Rdata")
load("intermediate_data/AlphaFold3/afs_dockq_test.Rdata")

dockq_scores <- data <- dplyr::bind_rows(
    afs_dockq,
    afs_dockq_test) |>
    dplyr::mutate(
        cif_fname = cif_fname |>
            stringr::str_replace("/nfs/turbo/umms-tromeara", "/home/maom/turbo_tromeara")) |>
    dplyr::select(
        dataset_id,
        structure_id,
        Poc4_sequence,
        cif_fname,
        pDockQ2_Poc4_Pup2)


# ESM2 Sequence model predictions
ems2_model_train_predicted_cv <- readr::read_tsv(
    file = "product/dockq_analysis/sequence_model/train_predicted_cv_20260117.tsv",
    show_col_types = FALSE) |>
    dplyr::select(
        Poc4_sequence,
        esm2_predict = predict)

esm2_model_test_predicted <- readr::read_tsv(
    file = paste0(
        "product/dockq_analysis/sequence_model",
        "/designs_fixinterface_T1.0_af2_v2_pldd2_gt96_score_lt0.6_20251103_predicted_20260113.tsv"),
    show_col_types = FALSE) |>
    dplyr::select(
        Poc4_sequence = fasta,
        esm2_predict = predict)

esm2_model_predictions <- dplyr::bind_rows(
    ems2_model_train_predicted_cv,
    esm2_model_test_predicted)

# PLiP Structure model predictions
plip_model_train_predicted_cv <- readr::read_tsv(
    file = "product/dockq_analysis/structure_model/train_predicted_cv_20260118.tsv",
    show_col_types = FALSE) |>
    dplyr::select(
        Poc4_sequence,
        plip_predict = predict)

plip_model_test_predicted <- readr::read_tsv(
    file = "product/dockq_analysis/structure_model/test_predicted_20260118.tsv",
    show_col_types = FALSE) |>
    dplyr::select(
        Poc4_sequence,
        plip_predict = predict)

plip_model_predictions <- dplyr::bind_rows(
    plip_model_train_predicted_cv,
    plip_model_test_predicted)

# Hinge geometry model predictions
hinge_model_train_predicted_cv <- readr::read_tsv(
    file = "intermediate_data/dockq_analysis/hinge_model/design_hinge_features_20260121.tsv",
    show_col_types = FALSE) |>
    dplyr::rename(hinge_distance = distance)

hinge_model_test_predicted <- readr::read_tsv(
    file = "intermediate_data/dockq_analysis/hinge_model/design_hinge_features_test_20260121.tsv",
    show_col_types = FALSE) |>
    dplyr::rename(hinge_distance = distance)

hinge_model_predictions <- dplyr::bind_rows(
    hinge_model_train_predicted_cv,
    hinge_model_test_predicted)

data <- dockq_scores |>
    dplyr::left_join(
        esm2_model_predictions,
        by = "Poc4_sequence") |>
    dplyr::left_join(
        plip_model_predictions,
        by = "Poc4_sequence") |>
    dplyr::left_join(
        hinge_model_predictions,
        by = "cif_fname") |>
    dplyr::left_join(
        growth_data_summary,
        by = "structure_id")

base_plot <- ggplot2::ggplot(
    data = data,
    mapping = ggplot2::aes(
        y = pDockQ2_Poc4_Pup2,
        color = dataset_id)) +
    ggplot2::theme_bw() +
    viridis::scale_color_viridis(
        option = "inferno",
        begin = 0.2,
        end = 0.8,
        discrete = TRUE) +
    ggplot2::scale_y_continuous("pDockQ2 between Poc4 vs. Pup2")

plot1 <- base_plot +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = esm2_predict),
        alpha = 0.8,
        shape = 16) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(!is.na(recovery_score), recovery_class == "Not Recovered"),
        mapping = ggplot2::aes(
            x = esm2_predict),
        color = "orange",
        shape = 1) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(!is.na(recovery_score), recovery_class == "Partial"),
        mapping = ggplot2::aes(
            x = esm2_predict),
        color = "yellow",
        shape = 1) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(!is.na(recovery_score), recovery_class == "Recovered"),
        mapping = ggplot2::aes(
            x = esm2_predict),
        color = "green",
        shape = 1) +
    ggrepel::geom_label_repel(
        data = data |> dplyr::filter(!is.na(recovery_score)),
        mapping = ggplot2::aes(
            x = esm2_predict,
            label = structure_id),
        size = 3,
        force = 10,
        min.segment.length = 0) +
    ggplot2::scale_x_continuous("ESM2 ML Predicted pDockQ2")

plot2 <- base_plot +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = plip_predict),
        alpha = 0.8,
        shape = 16) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(!is.na(recovery_score), recovery_class == "Not Recovered"),
        mapping = ggplot2::aes(
            x = plip_predict),
        color = "orange",
        shape = 1) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(!is.na(recovery_score), recovery_class == "Partial"),
        mapping = ggplot2::aes(
            x = plip_predict),
        color = "yellow",
        shape = 1) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(!is.na(recovery_score), recovery_class == "Recovered"),
        mapping = ggplot2::aes(
            x = plip_predict),
        color = "green",
        shape = 1) +
    ggrepel::geom_label_repel(
        data = data |> dplyr::filter(!is.na(recovery_score)),
        mapping = ggplot2::aes(
            x = esm2_predict,
            label = structure_id),
        size = 3,
        force = 10,
        min.segment.length = 0) +
    ggplot2::scale_x_continuous("PLiP ML Predicted pDockQ2")

plot3 <- base_plot +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = sqrt(hinge_distance)),
        alpha = 0.8,
        shape = 16) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(!is.na(recovery_score), recovery_class == "Not Recovered"),
        mapping = ggplot2::aes(
            x = sqrt(hinge_distance)),
        color = "orange",
        shape = 1) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(!is.na(recovery_score), recovery_class == "Partial"),
        mapping = ggplot2::aes(
            x = sqrt(hinge_distance)),
        color = "yellow",
        shape = 1) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(!is.na(recovery_score), recovery_class == "Recovered"),
        mapping = ggplot2::aes(
            x = sqrt(hinge_distance)),
        color = "green",
        shape = 1) +
    ggrepel::geom_label_repel(
        data = data |> dplyr::filter(!is.na(recovery_score)),
        mapping = ggplot2::aes(
            x = esm2_predict,
            label = sqrt(hinge_distance)),
        size = 3,
        force = 10,
        min.segment.length = 0) +
    ggplot2::scale_x_continuous(
        "Poc4 Pup2 Hinge Distance",
        breaks = sqrt(c(20, 25, 30, 35)),
        limits = sqrt(c(15, 35)),
        labels = c("20A", "25A", "30A", "35A"),
        trans = "reverse")


plot <- patchwork::wrap_plots(
    plot1, plot2, plot3,
    nrow = 1,
    axes = "collect_y",
    guides = "collect")

ggplot2::ggsave(
    filename = paste0(product_path, "/DockQ2_by_predictors_", date_code, ".pdf"),
    width = 9,
    height = 3.5,
    useDingbats = FALSE)

data |>
    readr::write_tsv(
        file = paste0(product_path, "/DockQ2_by_predictors_", date_code, ".tsv"))



##############################

base_plot <- ggplot2::ggplot(
    data = data,
    mapping = ggplot2::aes(
        y = pDockQ2_Poc4_Pup2)) +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous("pDockQ2 between Poc4 vs. Pup2")

plot3 <- base_plot +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = sqrt(hinge_distance)),
        size = 1.5,
        alpha = 0.4,
        shape = 16) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(!is.na(recovery_score), recovery_class == "Not Recovered"),
        mapping = ggplot2::aes(
            x = sqrt(hinge_distance)),
        color = "#E05F4B",
        shape = 16,
        alpha = 0.8,
        size = 3) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(!is.na(recovery_score), recovery_class == "Partial"),
        mapping = ggplot2::aes(
            x = sqrt(hinge_distance)),
        color = "#E0AC50",
        shape = 16,
        alpha = 0.8,
        size = 3) +
    ggplot2::geom_point(
        data = data |> dplyr::filter(!is.na(recovery_score), recovery_class == "Recovered"),
        mapping = ggplot2::aes(
            x = sqrt(hinge_distance)),
        color = "#DEE157",
        shape = 16,
        alpha = 0.8,
        size = 3) +
    ggrepel::geom_label_repel(
        data = data |> dplyr::filter(!is.na(recovery_score)),
        mapping = ggplot2::aes(
            x = sqrt(hinge_distance),
            label = structure_id),
        size = 3,
        force = 10,
        min.segment.length = 0) +
    ggplot2::scale_x_continuous(
        "Poc4 Pup2 Hinge Distance",
        breaks = sqrt(c(20, 25, 30, 35)),
        limits = sqrt(c(15, 35)),
        labels = c("20A", "25A", "30A", "35A"),
        trans = "reverse")

ggplot2::ggsave(
    filename = paste0(product_path, "/DockQ2_by_hinge_", date_code, ".pdf"),
    plot = plot3,
    width = 5,
    height = 5,
    useDingbats = FALSE)

data |>
    readr::write_tsv(
        file = paste0(product_path, "/DockQ2_by_hinge_", date_code, ".tsv"))
