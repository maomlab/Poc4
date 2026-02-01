compute_dockq_af3 <- function(
    json_fname,
    cif_fname,
    pae_cutoff = 10,
    dist_cutoff = 10,
    overwrite = FALSE,
    future = FALSE,
    verbose = FALSE) {

    if (!file.exists(json_fname)) {
        cat("WARNING: json_fname doesn't exist '", json_fname, "'\n")
    }
    if (!file.exists(cif_fname)) {
        cat("WARNING: cif_fname doesn't exist '", cif_fname, "'\n")
    }

    summary_fname <- cif_fname |>
        stringr::str_replace("[.]cif", paste0("_", pae_cutoff, "_", dist_cutoff, ".txt"))

    if (!overwrite && file.exists(summary_fname)) {
        if (verbose) {
            cat(
                "output file '", summary_fname, "' already exists, ",
                "and overwrite=FALSE, skipping computing it...\n",
                sep = "")
        }
    } else {
        cmd <- paste0(
            "python src/IPSAE/ipsae.py ",
            json_fname, " ",
            cif_fname, " ",
            pae_cutoff, " ",
            dist_cutoff)

        if (verbose) {
            cat(cmd, "\n")
        }

        system(cmd, intern=TRUE)
    }

    summary <- readr::read_table(
        file = summary_fname,
        skip = 1,
        show_col_types = FALSE)
    summary
}


compute_all_dockq_hpc <- function(
    results_path,
    pae_cutoff = 10,
    dist_cutoff = 10,
    overwrite = FALSE,
    verbose = FALSE) {

    data.frame(
        results_path = results_path) |>
    dplyr::rowwise() |>
    dplyr::do({
        x <- .
        results_path <- x$results_path[1]
        structure_id <- x$results_path |> basename()

        scores <- data.frame(
            cif_fname = Sys.glob(
                path = paste0(results_path, "/*/*.cif"))) |>
            dplyr::mutate(
                sample_index = cif_fname |> stringr::str_extract("sample-([0-9]+)", group = 1),
                json_fname = cif_fname |> stringr::str_replace("model.cif", "confidences.json")) |>
            dplyr::rowwise() |>
            dplyr::do({
                x <- .
                compute_dockq_af3(
                        json_fname = x$json_fname,
                        cif_fname = x$cif_fname,
                        pae_cutoff = pae_cutoff,
                        dist_cutoff = dist_cutoff,
                        overwrite = overwrite,
                        verbose = verbose) |>
                        dplyr::mutate(
                            structure_id = structure_id,
                            sample_index = x$sample_index,
                            .before = 1)
            })
        scores
    })
}

compute_all_dockq_afs <- function(
    results_path,
    pae_cutoff = 10,
    dist_cutoff = 10,
    overwrite = FALSE,
    verbose = FALSE) {

    data.frame(results_path = results_path) |>
    dplyr::rowwise() |>
    dplyr::do({
        x <- .
        results_path <- x$results_path[1]
        structure_id <- x$results_path |> basename()

        scores <- data.frame(
            cif_fname = Sys.glob(
                path = paste0(results_path, "/*.cif"))) |>
            dplyr::mutate(
                sample_index = cif_fname |> stringr::str_extract("model_([0-9]+)", group = 1),
                json_fname = cif_fname |> stringr::str_replace("model_([0-9]+).cif", "full_data_\\1.json")) |>
            dplyr::rowwise() |>
            dplyr::do({
                x <- .
                compute_dockq_af3(
                    json_fname = x$json_fname,
                    cif_fname = x$cif_fname,
                    pae_cutoff = pae_cutoff,
                    dist_cutoff = dist_cutoff,
                    verbose = TRUE) |>
                    dplyr::mutate(
                        structure_id = structure_id,
                        sample_index = x$sample_index,
                        .before = 1)
            })

        job_request_path <- Sys.glob(paste0(results_path, "/*job_request.json"))
        if (length(job_request_path) != 1 || !file.exists(job_request_path)) {
            cat("WARNING: job_request path is not valid: '", job_request_path, "'\n", sep = "")
        }


        job_request <- jsonlite::read_json(
            path = Sys.glob(paste0(results_path, "/*job_request.json")))

        sequences <- job_request[[1]]$sequences |>
            purrr::map_dfr(function(chain) {
                data.frame(
                    sequence = chain$proteinChain$sequence)
            }) |>
            dplyr::mutate(
                chain_id = LETTERS[1L : dplyr::n()],
                .before = 1)

        scores <- scores |>
            dplyr::left_join(
                sequences |>
                dplyr::rename(
                    Chn1 = chain_id,
                    Seq1 = sequence),
                by = "Chn1") |>
            dplyr::left_join(
                sequences |>
                dplyr::rename(
                    Chn2 = chain_id,
                    Seq2 = sequence),
                by = "Chn2")

        scores
    }) |>
    dplyr::ungroup()
}

submit_all_dockq_afs <- function(
    results_path,
    pae_cutoff = 10,
    dist_cutoff = 10,
    overwrite = FALSE,
    verbose = FALSE) {

    data.frame(
        results_path = results_path) |>
    dplyr::rowwise() |>
    dplyr::do({
        x <- .
        results_path <- x$results_path[1]
        structure_id <- x$results_path |> basename()


        scores <- data.frame(
            cif_fname = Sys.glob(
                path = paste0(results_path, "/*.cif"))) |>
            dplyr::mutate(
                sample_index = cif_fname |> stringr::str_extract("model_([0-9]+)", group = 1),
                json_fname = cif_fname |> stringr::str_replace("model_([0-9]+).cif", "full_data_\\1.json")) |>
            dplyr::rowwise() |>
            dplyr::do({
                x <- .
                result <- future::future({
                    compute_dockq_af3(
                        json_fname = x$json_fname,
                        cif_fname = x$cif_fname,
                        pae_cutoff = pae_cutoff,
                        dist_cutoff = dist_cutoff,
                        overwrite = overwrite,
                        verbose = verbose) |>
                        dplyr::mutate(
                            structure_id = structure_id,
                            sample_index = x$sample_index,
                            .before = 1)
                })
                tibble::tibble(job = list(result))
            })
        scores
    })
}
