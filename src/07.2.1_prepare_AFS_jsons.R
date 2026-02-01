library(jsonlite)
library(seqinr)


template <- readLines("src/AFS_refold_input_template.json")

designs <- readr::read_tsv(
    "product/frame2seq_design/designs_fixinterface_T1.0_af2_v2_pldd2_gt96_score_lt0.6_20251103/summary.tsv",
    show_col_types = FALSE)

output_path <- "product/frame2seq_design/designs_fixinterface_T1.0_af2_v2_pldd2_gt96_score_lt0.6_20251103/AFS_inputs"

if (!dir.exists(output_path)) {
    cat("Creating path '", output_path, "'\n", sep = "")
    dir.create(output_path, recursive = TRUE)
}

designs <- tibble::tibble(
    structure_id = names(designs),
    sequence = designs)

designs |>
    dplyr::rowwise() |>
    dplyr::do({
        x <- .

        input <- template |>
            purrr::map(~stringr::str_replace(.x, "DESIGN_NAME", x$seq_id)) |>
            purrr::map(~stringr::str_replace(.x, "DESIGN_SEQUENCE", x$fasta)) |>
            readr::write_lines(
                file = paste0(output_path, "/", x$seq_id, ".json"))
        data.frame()
    })
