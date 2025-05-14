# Submit AlphaFold2 prediction jobs after having pre-generated MSAs

# 1) find all the structures that have pre-generated MSAs
# 2) find all the structures with that have generated AF2 predictions
# 3) save the remaining ones in intermediate_data/parafold/structures.todo
# 4) submit to slurm compute jobs that iterate through the structures.todo

prepare_AF2_jobs <- function(
    input_path,
    output_path,
    logs_path,
    slurm_account,
    slurm_partition,
    expected_files = c("full", "minimal"),
    array_range = NULL,
    verbose = FALSE) {


    if (!dir.exists(logs_path)) {
        cat("Creating '", logs_path, "' ...\n", sep = "")
        dir.create(logs_path, recursive = TRUE)
    }
    
    expected_files <- match.arg(expected_files)
    
    if (expected_files == "full") {
        n_models <- 5
        expected_files <- tibble::tibble(
            fname = c(
                paste0("ranked_", 0:(n_models - 1), ".pdb"),
                "ranking_debug.json",
                paste0("result_model_", 1:n_models, "_pred_0.pkl"),
                "timings.json",
                paste0("unrelaxed_model_", 1:n_models, "_pred_0.pdb")))
    } else if (expected_files == "minimal") {
        n_models <- 4
        expected_files <- tibble::tibble(
            fname = c(
                paste0("result_model_", 1:n_models, "_pred_0.pkl"),
                paste0("unrelaxed_model_", 1:n_models, "_pred_0.pdb")))
    }

    all_structure_meta <- tibble::tibble(
        path = Sys.glob(
            paste0(output_path, "/*/input/features.pkl")) |>
            stringr::str_replace("/features.pkl$", ""),
        structure_id = path |>
            stringr::str_extract(
                pattern = paste0("^", output_path, "/(.*)/input"),
                group = 1))
    
    output_meta <- tibble::tibble(
        structure_id = list.dirs(
            output_path,
            full.names = FALSE,
            recursive = FALSE)) |>
        dplyr::mutate(
            structure_path = paste0(output_path, "/", structure_id, "/input"))
    
    verbose <- FALSE
    prepared_structures_meta <- output_meta |>
        dplyr::rowwise() |>
        dplyr::do({
            structure_id <- .$structure_id[1]
            structure_path <- .$structure_path[1]
            structure_files <- tibble::tibble(
                fname = list.files(path = structure_path, recursive = TRUE))
            missing_files <- expected_files |>
                dplyr::anti_join(structure_files, by = "fname")
            if (nrow(missing_files) > 0) {
                if (verbose) {
                    cat("WARNING: output for structure_id '", structure_id, "' exists but is incomplete\n", sep = "")
                    print(missing_files)
                }
                result <- data.frame()
            } else {
                result <- data.frame(structure_id = structure_id)
            }
            result
        })

    if (nrow(prepared_structures_meta) > 0) {
        structures_todo_meta <- all_structure_meta |>
            dplyr::select(structure_id) |>
            dplyr::anti_join(
                prepared_structures_meta,
                by = "structure_id")
    } else {
        structures_todo_meta <- all_structure_meta |>
            dplyr::select(structure_id)
    }
    
    
    structures_todo_meta |>
        readr::write_tsv(
            file = paste0(output_path, "/structures.todo"),
            col_names = FALSE)
    
    n_structures_remaining <- readr::read_tsv(
        paste0(output_path, "/structures.todo"),
        col_names = "structure_id",
        show_col_types = FALSE) |>
        nrow()
    cat("N structures remaining: ", n_structures_remaining, "\n", sep = "")

    if (is.null(array_range)) {
        array_range <- paste0("1-", n_structures_remaining)
    }
    
    submit_sbatch <- paste0("sbatch ",
        "--account=", slurm_account, " ",
        "--partition=", slurm_partition, " ",
        "--array=", array_range, " ",
        "--ntasks-per-node=1 ",
        "--cpus-per-task=2 ",
        "--time=3:00:00 ",
        "--mem-per-cpu=40GB ",        
        "--gres=gpu:1 ",
        "--output='", logs_path, "/%x_%j.out' ",
        "--error='", logs_path, "/%x_%j.err' ",
        "--export=",
          "INPUT_PATH='", input_path, "',",
          "OUTPUT_PATH='", output_path, "' ",
        "src/submit_alphafold_monomer.slurm")
    cat(submit_sbatch, "\n", sep = "")
    submit_sbatch
}

submit_sbatch <- prepare_AF2_jobs(
    input_path = "intermediate_data/parafold/Frame2Seq_fixinterface_T1.0/input",
    output_path = "intermediate_data/parafold/Frame2Seq_fixinterface_T1.0/output",
    logs_path = "intermediate_data/parafold/Frame2Seq_fixinterface_T1.0/logs_AF2",
    slurm_account = "tromeara99",
    slurm_partition = "gpu",
    array_range = "1-3",
    verbose = TRUE)

system(submit_sbatch)


# sbatch --account=tromeara99 --partition=gpu --array=1-5000 --ntasks-per-node=1 --cpus-per-task=2 --time=3:00:00 --mem-per-cpu=3GB --gres=gpu:1 src/subit_alphafold_monomer.slurm
