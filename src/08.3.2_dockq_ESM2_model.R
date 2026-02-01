source("src/h2o_utils.R")


library(h2o)

#############
# Load Data #
#############

#Train data
load("intermediate_data/AlphaFold3/afs_dockq.Rdata")
afs_dockq_esm_embeddings <- arrow::read_parquet(
    file = "intermediate_data/AlphaFold3/afs_dockq_ems_embedding.parquet")

data_train <- afs_dockq |>
    dplyr::select(
        structure_id,
        pDockQ2_Poc4_Pup2) |>
    dplyr::bind_cols(
        afs_dockq_esm_embeddings |> dplyr::select(-id))


# Test data
designs_round2_summary <- readr::read_tsv(
    file = "product/frame2seq_design/designs_fixinterface_T1.0_af2_v2_pldd2_gt96_score_lt0.6_20251103/summary.tsv",
    show_col_types = FALSE) |>
    dplyr::distinct(fasta, .keep_all = TRUE)

frame2seq_esm_embeddings <- arrow::read_parquet(
    file = "intermediate_data/ESM2_embeddings/designs_fixinterface_T1.0_af2_v2_pldd2_gt96_score_lt0.6_20251103.parquet")


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

date_code <- "20260117"


# Start the H2O cluster (locally)
h2o::h2o.init()


#############
# Fit Model #
#############

train <- h2o::as.h2o(data_train, destination_frame = "data_train")
y <- "pDockQ2_Poc4_Pup2"
ids <- c("structure_id")
x <- setdiff(names(train), c(y, ids))

# Run AutoML for 20 base models
aml <- h2o::h2o.automl(
    x = x,
    y = y,
    training_frame = train,
    max_models = 20,
    seed = 1,
    keep_cross_validation_predictions = TRUE)


# View the AutoML Leaderboard
lb <- aml@leaderboard
print(lb, n = nrow(lb))


leader <- aml@leader

leader |> h2o::h2o.download_model(
    path = intermediate_path)


h2o_model_report(
    model = leader,
    output_path = product_path,
    date_code = date_code)


# cross_validation_predictions appears to only work for single models
# and resturns a list of H2OFrame objects of length nfold, where the
# predictions for each fold that the point is in are zeroed out.
glm <- gbm <- h2o.get_best_model(aml, algorithm = "GBM")

train_predicted_cv <- glm |>
    h2o::h2o.cross_validation_predictions() |>
    purrr::map_dfc(as.data.frame) |>
    dplyr::transmute(
        predict = rowSums(dplyr::across(tidyselect::everything())))

train_predicted_cv <- dplyr::bind_cols(
    afs_dockq,
    train_predicted_cv)

readr::write_tsv(
    train_predicted_cv,
    file = paste0(product_path, "/train_predicted_cv_", date_code, ".tsv"))




test <- h2o::as.h2o(
    frame2seq_esm_embeddings,
    destination_frame = "data_test")

test_predicted <- h2o::h2o.predict(
    object = leader,
    test) |>
    as.data.frame()

test_predicted <- dplyr::bind_cols(
    designs_round2_summary,
    test_predicted)

readr::write_tsv(
    test_predicted,
    file = paste0(product_path, "/designs_fixinterface_T1.0_af2_v2_pldd2_gt96_score_lt0.6_20251103_predicted_", date_code, ".tsv"))


plot <- ggplot2::ggplot(data = train_predicted) +
    ggplot2::theme_bw() +
    ggplot2::geom_abline(
        intercept = 0,
        slope = 1,
        color = "darkgray") +
    ggplot2::geom_point(
        data = train_predicted,
        mapping = ggplot2::aes(
            x = predict,
            y = pDockQ2_Poc4_Pup2),
        alpha = 0.8,
        shape = 16) +
    ggplot2::geom_rug(
        data = test_predicted,
        mapping = ggplot2::aes(
            x = predict)) +
    ggplot2::coord_equal() +
    ggplot2::ggtitle(
        "ML DockQ2 Prediction",
        subtitle = "based on Poc4 ESM Embedding") +
    ggplot2::labs(
        x = "h2o AutoML Ensemble Predictor Score",
        y = "pDockQ2 between Poc4 vs. Pup2 (max)")

ggplot2::ggsave(
    filename = "product/dockq_analysis/afs_esm_ml_prediction_20260112.pdf",
    width = 5,
    height = 5,
    useDingbats = FALSE)


test_predicted_new <- test_predicted |>
    dplyr::anti_join(
        train_predicted,
        by = c("fasta" = "Poc4_sequence")) |>
    dplyr::arrange(dplyr::desc(predict)) |>
    head(150)


####### Write out AF3 jsons for the top ESM2 predicted sequences ####
template <- readLines("src/AFS_refold_input_template.json")

test_predicted_new |>
    dplyr::mutate(
        batch_id = paste0("ESM2_DockQ2_predicted_set", rep(1:5, times = 30))) |>
    dplyr::rowwise() |>
    dplyr::do({
        x <- .

        batch_path <- paste0(intermediate_path, "/", x$batch_id)
        if (dir.exists(batch_path)) {
            warning("Batch path '", batch_path, "' exists.\n", sep = "")
        } else {
            dir.create(batch_path, recursive = TRUE)
        }

        input <- template |>
            purrr::map(~stringr::str_replace(.x, "DESIGN_NAME", x$seq_id)) |>
            purrr::map(~stringr::str_replace(.x, "DESIGN_SEQUENCE", x$fasta)) |>
            readr::write_lines(
                file =  paste0(batch_path, "/", x$seq_id, ".json"))
        data.frame()
    })
