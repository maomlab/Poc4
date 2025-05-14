
prepare_MSA_jobs <- function(
    input_path,
    output_path,
    logs_path,
    slurm_account,
    slurm_partition,
    array_range = NULL,
    verbose = FALSE) {

    if (!dir.exists(output_path)) {
        cat("Creating '", output_path, "' ...\n", sep = "")
        dir.create(output_path, recursive = TRUE)
    }
    if (!dir.exists(logs_path)) {
        cat("Creating '", logs_path, "' ...\n", sep = "")
        dir.create(logs_path, recursive = TRUE)
    }

    expected_msa_files <- tibble::tibble(fname = c(
        "input/features.pkl",
        "input/msas/bfd_uniref_hits.a3m",
        "input/msas/mgnify_hits.sto",
        "input/msas/pdb_hits.hhr",
        "input/msas/uniref90_hits.sto"))
    
    all_structure_meta <- tibble::tibble(
        structure_id = list.files(
            input_path,
            pattern = ".*",
            full.names = FALSE,
            recursive = FALSE))
    
    output_meta <- tibble::tibble(
        structure_id = list.dirs(
            output_path,
            full.names = FALSE,
            recursive = FALSE)) |>
        dplyr::mutate(
            msa_path = paste0(output_path, "/", structure_id))
    
    if (nrow(output_meta) == 0) {
        structures_todo_meta <- all_structure_meta
    } else {
        prepared_structures_meta <- output_meta |>
            dplyr::rowwise() |>
            dplyr::do({
                structure_id <- .$structure_id[1]
                msa_path <- .$msa_path[1]
                msa_files <- tibble::tibble(
                    fname = list.files(path = msa_path, recursive = TRUE))
                missing_msa_files <- expected_msa_files |>
                    dplyr::anti_join(msa_files, by = "fname")
                if (nrow(missing_msa_files) > 0) {
                    cat("WARNING: output for structure_id '", structure_id, "' exists but is incomplete\n", sep = "")
                    print(missing_msa_files)
                    result <- data.frame()
                } else {
                    result <- data.frame(structure_id = structure_id)
                }
                result
            })
    
        structures_todo_meta <- all_structure_meta |>
            dplyr::anti_join(
                prepared_structures_meta,
                by = "structure_id")
    }
    
    structures_todo_meta |>
        readr::write_tsv(
            file = paste0(output_path, "/msa.todo"),
            col_names = FALSE)
    
    n_structures_remaining <- readr::read_tsv(
        paste0(output_path, "/msa.todo"),
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
        "--cpus-per-task=8 ",
        "--time=10:00:00 ",
        "--mem-per-cpu=6GB ",
        "--output='", logs_path, "/%x_%j.out' ",
        "--error='", logs_path, "/%x_%j.err' ",
        "--export=",
          "INPUT_PATH='", input_path, "',",
          "OUTPUT_PATH='", output_path, "' ",
        "src/submit_alphafold_monomer_msa.slurm")

    cat(submit_sbatch, "\n", sep = "")
    submit_sbatch
}

# submit_sbatch <- prepare_MSA_jobs(
#     input_path = "intermediate_data/parafold/input",
#     output_path = "intermediate_data/parafold/output",
#     slurm_account = "tromeara99",
#     slurm_partition = "standard",
#     verbose = TRUE)
# 
# system(submit_sbatch)



submit_sbatch <- prepare_MSA_jobs(
    input_path = "intermediate_data/parafold/Frame2Seq_fixinterface_T1.0/input",
    output_path = "intermediate_data/parafold/Frame2Seq_fixinterface_T1.0/output",
    logs_path = "intermediate_data/parafold/Frame2Seq_fixinterface_T1.0/logs",
    slurm_account = "tromeara99",
    slurm_partition = "standard",
    array_range = "1-2000",
    verbose = TRUE)

system(submit_sbatch)


# sbatch --account=tromeara99 --partition=standard --array=1-5000 --ntasks-per-node=1 --cpus-per-task=8 --time=10:00:00 --mem-per-cpu=6GB src/sub_alphafold_multimer_msa.slurm
