

frame2seq <- readr::read_tsv(
    paste0("intermediate_data/frame2seq_design/Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_5000.tsv"),
    show_col_types = FALSE)


frame2seq_metadata <- data.frame(
    path = c(
        "Cauris_884_AF2_Frame2seq_T1.0_5000_metadata.tsv",
        "Cauris_884_irc25_AF3_Frame2seq_fixconserved_T1.0_50000_metadata.tsv",
        "Cauris_884_irc25_AF3_Frame2seq_fixinterface_T1.0_5000_metadata.tsv",
        "Cauris_884_irc25_AF3_Frame2seq_T0.01_5000_metadata.tsv",
        "Cauris_884_irc25_AF3_Frame2seq_T0.1_5000_metadata.tsv",
        "Cauris_884_irc25_AF3_Frame2seq_T1.0_5000_metadata.tsv",
        "Cauris_884_irc25_AF3_Frame2seq_T10_5000_metadata.tsv")) |>
    tidyr::separate_wider_regex(
        cols = path,
        patterns = c(
            "Cauris",                    "_",
            structure = ".*"     ,       "_",
            "Frame2seq",                 "_",
            constraint = ".*",
            "T", temperature = "[^_]+",  "_",
            n_samples = "[0-9]+",        "_",
            "metadata.tsv"),
        cols_remove = FALSE) |>
    dplyr::mutate(
        constraint = constraint |> stringr::str_replace("_$", ""),
        temperature = temperature |> as.numeric(),
        n_samples = n_samples |> as.numeric()) |>
    dplyr::rowwise() |>
    dplyr::do({
        data <- .
        data.frame(
            data,
            readr::read_tsv(
                paste0("intermediate_data/frame2seq_design/", data$path[1]),
                show_col_types = FALSE))
    }) |>
    dplyr::mutate(
        set_label = paste0(
            structure, " ", constraint, " T", temperature))


plot <- ggplot2::ggplot(
    data = frame2seq_metadata |>
        dplyr::filter(temperature < 10) |>
        dplyr::sample_n(size = dplyr::n(), replace = FALSE)) +
    ggplot2::theme_bw() +
    ggplot2::geom_jitter(
        mapping = ggplot2::aes(
            x = score,
            y = recovery_percent,
            color = set_label),
        width = 0.05,
        height = 0.8,
        shape = 16,
        size = 0.4,
        alpha = 0.4) +
    ggplot2::geom_smooth(
        mapping = ggplot2::aes(
            x = score,
            y = recovery_percent),
        formula = y ~ x,
        method = "lm",
        color = "blue") +
    ggplot2::scale_x_continuous("Frame2seq Score") +
    ggplot2::scale_y_continuous("Recovery Percent") +
    ggplot2::scale_color_discrete("Prediction Set") +
    ggplot2::ggtitle(
        label = "Frame2seq Sequence Recovery by Score")

ggplot2::ggsave(
    filename = "product/frame2seq_design/score_vs_recovery_20240902.pdf",
    plot = plot,
    width = 8,
    height = 8,
    useDingbats = FALSE)
