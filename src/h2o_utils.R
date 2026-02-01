library(h2o)


h2o_model_report <- function(model, output_path, date_code) {

    attributes(model)$model$model_summary |>
        as.data.frame() |>
        readr::write_tsv(
            file = paste0(
                output_path, "/model_summary_", date_code, ".tsv"))

    attributes(model)$model$training_metrics |>
        show() |>
        capture.output() |>
        cat(
            sep = "\n",
            file = paste0(
                output_path, "/training_metrics_", date_code, ".txt"))

    attributes(model)$model$cross_validation_metrics |>
        show() |>
        capture.output() |>
        cat(
            sep = "\n",
            file = paste0(
                output_path, "/cross_validation_metrics_", date_code, ".txt"))

    attributes(model)$model$metalearner_model |>
        show() |>
        capture.output() |>
        cat(
            sep = "\n",
            file = paste0(
                output_path, "/metalearner_model_", date_code, ".txt"))
}




h2o_shap_summary_plot <- function(
    model,
    newdata,
    output_path,
    date_code,
    ...) {

    plot_data <- h2o.predict_contributions(
        object = model,
        newdata = newdata,
        ...) |>
        as.data.frame()

    plot_data |>
        readr::write_tsv(
            file = paste0(output_path, "/shap_", date_code, ".tsv"))

    shap_plot <- h2o.shap_summary_plot(
        model, newdata,
        sample_size = Inf,
        ...)

    ggplot2::ggsave(
        filename = paste0(output_path, "/shap_", date_code, ".pdf"),
        plot = shap_plot,
        width = 6,
        height = 6)
}
