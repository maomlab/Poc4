growth_data <- readr::read_csv(
    file = "data/growth_data_20251216.csv",
    show_col_types = FALSE) |>
    dplyr::filter(!is.na(Strain)) |>
    dplyr::select(-tidyselect::starts_with("...")) |>
    tidyr::pivot_longer(
        cols = c(-Strain),
        names_to = "condition_measurement_replica",
        values_to = "AUC") |>
    dplyr::filter(!is.na(AUC)) |>
    tidyr::separate_wider_delim(
        cols = "condition_measurement_replica",
        delim = "_",
        names = c("temperature_C", "measurement", "replica")) |>
    dplyr::mutate(
        Strain = Strain |>
            stringr::str_replace("^Seq([0-9]+)$", "seq_\\1"),
        temperature_C = temperature_C |>  as.factor(),
        replica = replica |> as.factor()) |>
    dplyr::select(
        Strain,
        temperature_C,
        replica,
        AUC)

growth_data_summary <- growth_data |>
    dplyr::group_by(Strain, temperature_C) |>
    dplyr::summarize(
        AUC_mean = mean(AUC),
        .groups = "drop") |>
    tidyr::pivot_wider(
        id_cols = Strain,
        names_from = temperature_C,
        values_from = AUC_mean) |>
    dplyr::mutate(
        AUC_diff = `42C` - `30C`)

growth_data_summary <- growth_data_summary |>
    dplyr::mutate(
        PC = growth_data_summary |>
            dplyr::filter(Strain == "AR0382") |>
            dplyr::pull("AUC_diff"),
        NC = growth_data_summary |>
            dplyr::filter(Strain |> stringr::str_detect("Del_poc4")) |>
            dplyr::summarize(AUC_diff_mean = mean(AUC_diff)) |>
            dplyr::pull("AUC_diff_mean"),
        score = (AUC_diff - NC) / (PC - NC),
        class = dplyr::case_when(
            score >= 0.5 ~ "Recovered",
            score >= 0.2 ~ "Parial",
            TRUE ~ "Not Recovered") |>
            as.factor())



save(growth_data, file = "intermediate_data/growth_data.Rdata")
save(growth_data_summary, file = "intermediate_data/growth_data_summary.Rdata")
