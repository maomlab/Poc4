load("intermediate_data/growth_data_summary.Rdata")

load("intermediate_data/AlphaFold3/AFS_20251003.Rdata")
load("intermediate_data/AlphaFold3/AFS_20251210.Rdata")
load("intermediate_data/AlphaFold3/AFS_20251224.Rdata")



output_path <- "product/dockq_analysis"
if (!dir.exists(output_path)) {
    cat("Creating output directory '", output_path, "'\n", sep = "")
    dir.create(output_path)
}

date_code <- "20251230"


# Compare Dataset
score_means <- dplyr::bind_rows(
        AFS_20251003,
        AFS_20251210,
        AFS_20251224) |>
    dplyr::mutate(
        dataset_id = paste0(dataset_id, "_mean")) |>
    dplyr::group_by(
        dataset_id, structure_id) |>
    dplyr::summarize(
        score_mean = mean(pDockQ2_Poc4_Pup2),
        .groups = "drop") |>
    tidyr::pivot_wider(
        id_cols = c("structure_id"),
        names_from = "dataset_id",
        values_from = "score_mean")

score_sds <- dplyr::bind_rows(
        AFS_20251003,
        AFS_20251210,
        AFS_20251224) |>
    dplyr::mutate(
        dataset_id = paste0(dataset_id, "_sd")) |>
    dplyr::group_by(
        dataset_id, structure_id) |>
    dplyr::summarize(
        score_sd = sd(pDockQ2_Poc4_Pup2),
        .groups = "drop") |>
    tidyr::pivot_wider(
        id_cols = c("structure_id"),
        names_from = "dataset_id",
        values_from = "score_sd")

plot_data <- dplyr::inner_join(
    score_means,
    score_sds,
    by = c("structure_id"))

plot_data <- plot_data |>
    dplyr::left_join(
        growth_data_summary |>
        dplyr::select(
            structure_id = Strain,
            recovery_score = score,
            recovery_class = class),
        by = "structure_id")


# AFS_20251210 vs AFS_20251003
plot <- ggplot2::ggplot(data = plot_data) +
    ggplot2::theme_bw() +
    ggplot2::geom_smooth(
        mapping = ggplot2::aes(
            x = AFS_20251003_mean,
            y = AFS_20251210_mean),
        formula = y ~ x,
        method = "lm") +
    ggplot2::geom_crossbar(
        mapping = ggplot2::aes(
            x = AFS_20251003_mean,
            y = AFS_20251210_mean,
            xmin = AFS_20251003_mean - AFS_20251003_sd,
            xmax = AFS_20251003_mean + AFS_20251003_sd),
        color = "darkgray") +
    ggplot2::geom_errorbar(
        mapping = ggplot2::aes(
            x = AFS_20251003_mean,
            y = AFS_20251210_mean,
            ymin = AFS_20251210_mean - AFS_20251210_sd,
            ymax = AFS_20251210_mean + AFS_20251210_sd),
        color = "darkgray") +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = AFS_20251003_mean,
            y = AFS_20251210_mean,
            color = recovery_class),
        alpha = 0.8,
        shape = 16) +
    ggplot2::coord_equal(
        xlim = c(0, .35),
        ylim = c(0, .35)) +
    ggplot2::ggtitle("pDockQ2 Pup2 Poc4 FrameSeq Designs")

ggplot2::ggsave(
    filename = paste0(output_path, "/AFS_20251210_vs_AFS_20251003_", date_code, ".pdf"),
    width = 5,
    height = 5,
    useDingbats = FALSE)

#################################

load("intermediate_data/AlphaFold3/afs_dockq.Rdata")
load("intermediate_data/AlphaFold3/afs_dockq_test.Rdata")

afs_modeling_data <- dplyr::bind_rows(
    afs_dockq,
    afs_dockq_test) |>
    dplyr::group_by(structure_id) |>
    dplyr::arrange(dplyr::desc(pDockQ2_Poc4_Pup2)) |>
    dplyr::slice_head(n = 1) |>
    dplyr::left_join(
        growth_data_summary |>
        dplyr::select(
            structure_id = Strain,
            recovery_score = score,
            recovery_class = class),
        by = "structure_id")

plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_histogram(
        data = afs_modeling_data,
        mapping = ggplot2::aes(
            x = pDockQ2_Poc4_Pup2),
        bins = 50) +
    ggplot2::geom_rug(
        data = afs_modeling_data |>
            dplyr::filter(!is.na(recovery_class)),
        mapping = ggplot2::aes(
            x = pDockQ2_Poc4_Pup2,
            color = recovery_class)) +
    ggplot2::ggtitle(
        "Distribution of Max DockQ2 Score Per Design",
        paste0(
            nrow(afs_modeling_data),
            " Designed Poc4 vs. Pup2 in AFS computed structures")) +
    ggplot2::labs(
        x = "DockQ2 Score")

ggplot2::ggsave(
    filename = paste0(output_path, "/AFS_scores_histogram_", date_code, ".pdf"),
    width = 6,
    height = 4,
    useDingbats = FALSE)
