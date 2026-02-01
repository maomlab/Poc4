library(h2o)
source("src/h2o_utils.R")

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

train_data_path <- "intermediate_data/dockq_analysis/structure_model/afs_dockq_features_ml.csv"
load("intermediate_data/AlphaFold3/afs_dockq.Rdata")

test_data_path <- "intermediate_data/dockq_analysis/structure_model/afs_dockq_test_features_ml.csv"
load("intermediate_data/AlphaFold3/afs_dockq_test.Rdata")

#############
# Fit Model #
#############

# Start the H2O cluster (locally)
h2o::h2o.init()


train_data <- h2o::h2o.importFile(
    path = train_data_path,
    destination_frame = "train_data")
test_data <- h2o::h2o.importFile(
    path = test_data_path,
    destination_frame = "test_data")

# Identify predictors and response
y <- "DockQ2"
ids <- c("design_id", "cif_path")
x <- setdiff(names(train_data), c(y, ids))

#####################
#  Train data Model #
#####################
# Run AutoML for 20 base models
aml <- h2o.automl(
    x = x,
    y = y,
    training_frame = train_data,
    max_models = 20,
    seed = 1,
    keep_cross_validation_predictions = TRUE)

# View the AutoML Leaderboard and pick top single model
lb <- aml@leaderboard
print(lb, n = nrow(lb))

gbm <- h2o.get_best_model(aml, algorithm = "GBM")
gbm_path <- h2o.download_model(
    model = gbm,
    path = intermediate_path)

gbm <- h2o.upload_model(
    path = paste0(intermediate_path, "/GBM_5_AutoML_1_20260118_82422"))

# save model reports
h2o_model_report(
    model = gbm,
    output_path = product_path,
    date_code = date_code)

h2o_shap_summary_plot(
    model = gbm,
    newdata = train_data,
    output_path = product_path,
    date_code = date_code)

# Save cross validated predictions
train_predicted_cv <- gbm |>
    h2o::h2o.cross_validation_predictions() |>
    purrr::map_dfc(as.data.frame) |>
    dplyr::transmute(
        predict = rowSums(dplyr::across(tidyselect::everything())))

train_predicted_cv <- afs_dockq |>
    dplyr::left_join(
        dplyr::bind_cols(
            train_data |>
                as.data.frame() |>
                dplyr::select(cif_fname = cif_path),
            train_predicted_cv),
        by = c("cif_fname"))

readr::write_tsv(
    train_predicted_cv,
    file = paste0(product_path, "/train_predicted_cv_", date_code, ".tsv"))


# Predict on held out test data
test_predicted <- h2o::h2o.predict(
    object = gbm,
    test_data) |>
    as.data.frame()

test_predicted <- dplyr::bind_cols(
    afs_dockq_test,
    test_predicted)

readr::write_tsv(
    test_predicted,
    file = paste0(product_path, "/test_predicted_", date_code, ".tsv"))


########################
#  Combined data Model #
########################

combined_data_path <- "intermediate_data/dockq_analysis/structure_model/afs_dockq_combined_features_ml.csv"

combined_data <- h2o::h2o.importFile(
    path = combined_data_path,
    destination_frame = "combined_data")

# Identify predictors and response
y <- "DockQ2"
ids <- c("design_id", "cif_path")
x <- setdiff(names(combined_data), c(y, ids))



# Run AutoML for 20 base models
aml <- h2o.automl(
    x = x,
    y = y,
    training_frame = combined_data,
    max_models = 20,
    seed = 1,
    keep_cross_validation_predictions = TRUE,
    keep_cross_validation_fold_assignment = TRUE)

# View the AutoML Leaderboard and pick top single model
lb <- aml@leaderboard
print(lb, n = nrow(lb))

gbm <- h2o.get_best_model(aml, algorithm = "GBM")

dir.create(
    path = paste0(intermediate_path, "/combined_model"),
    showWarning = FALSE,
    recursive = TRUE)

gbm_path <- h2o.download_model(
    model = gbm,
    path = paste0(intermediate_path, "/combined_model"))

gbm_path <- "intermediate_data/dockq_analysis/structure_model/combined_model/GBM_2_AutoML_2_20260129_150314"
gbm <- h2o.upload_model(
    path = "intermediate_data/dockq_analysis/structure_model/combined_model/GBM_2_AutoML_2_20260129_150314")

# save model reports

dir.create(
    path = paste0(product_path, "/combined_model"),
    showWarning = FALSE,
    recursive = TRUE)

h2o_model_report(
    model = gbm,
    output_path = paste0(product_path, "/combined_model"),
    date_code = date_code)

h2o_shap_summary_plot(
    model = gbm,
    newdata = combined_data,
    output_path = paste0(product_path, "/combined_model"),
    date_code = date_code)

# Save cross validated predictions
combined_predicted_cv <- gbm |>
    h2o::h2o.cross_validation_predictions() |>
    purrr::map_dfc(as.data.frame) |>
    dplyr::transmute(
        predict = rowSums(dplyr::across(tidyselect::everything())))

combined_predicted_cv <- dplyr::bind_rows(
    afs_dockq,
    afs_dockq_test) |>
    dplyr::left_join(
        dplyr::bind_cols(
            train_data |>
                as.data.frame() |>
                dplyr::select(cif_fname = cif_path),
            combined_predicted_cv),
        by = c("cif_fname"))

readr::write_tsv(
    combined_predicted_cv,
    file = paste0(product_path, "combined_model/train_predicted_cv_", date_code, ".tsv"))
