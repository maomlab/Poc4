source("src/ipsea_util.R")

future::plan(
    future.batchtools::batchtools_slurm(
        template = system.file(
            "hpc",
            "batchtools.greatlakes.tmpl",
            package = "CalCEN"),
        resources = list(
            ncpus = 1,
            memory = 30000)))


#Pup2
# >B9J08_005129 COORDS:PEKT02000010_C_auris_B8441:509141-509899W, translated using codon table 12 (252 amino acids)
# MFLTRSEYDRGVSTFSPEGRLFQVEYSLEAIKLGSTAIGVSTSEGVILGVEKRVTSPLLE
# SSSIEKIVEIDTHIGCAMSGLTADARSMIDHARVASLSHNLYYDEAMQVESLTQSVCDLA
# LRFGEGASGEKRLMSRPFGVALLIAGIDENGPQLYHAEPSGTFYRYDAKAIGSGSEGAQT
# ELQNEYHKSLTLKEAELLTLKILKQVMEEKLDCKNAQLASITKEGGFQVYDDEKTDKLIK
# ELNETVTDDTQI*

#POC4
# >B9J08_000884 COORDS:PEKT02000002_C_auris_B8441:821389-821054C, translated using codon table 12 (111 amino acids)
# MDSAVLPLAAGEYTVSVTKGDAPKKPNVVFVTENSSTDIGSYIYAVPGKKDVFTAVLQGS
# ADSRVQDVATRIGTLLVKKSGCPSYVCVSGAVSEMDEMELVREVLTKSRSS*


################
# AFS_20251003 #
################

submit_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_20251003/seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE)

AFS_20251003 <- compute_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_20251003/fold_seq*"),
    pae_cutoff = 10,
    dist_cutoff = 10,
    verbose = TRUE) |>
    dplyr::filter(
        Chn1 == "E", #Pup2
        Chn2 == "K", #Poc4
        Type == "max") |>
    dplyr::transmute(
        dataset_id = "AFS_20251003",
        structure_id = structure_id |> stringr::str_replace("fold_", ""),
        sample_index = sample_index |> as.numeric(),
        pDockQ2_Poc4_Pup2 = pDockQ2)

save(AFS_20251003, file = "intermediate_data/AlphaFold3/AFS_20251003.Rdata")

################
# AFS_20251210 #
################

AFS_20251210 <- compute_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_20251210/fold_seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE) |>
    dplyr::filter(
        Chn1 == "D", #Pup2
        Chn2 == "K", #Poc4
        Type == "max") |>
    dplyr::transmute(
        dataset_id = "AFS_20251210",
        structure_id = structure_id |> stringr::str_replace("fold_", ""),
        sample_index = sample_index |> as.numeric(),
        pDockQ2_Poc4_Pup2 = pDockQ2)

save(AFS_20251210, file = "intermediate_data/AlphaFold3/AFS_20251210.Rdata")


################
# AFS_20251224 #
################

AFS_20251224 <- compute_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_20251224/seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE) |>
    dplyr::filter(
        Chn1 == "D", #Pup2
        Chn2 == "K", #Poc4
        Type == "max") |>
    dplyr::transmute(
        dataset_id = "AFS_20251224",
        structure_id = structure_id |> stringr::str_replace("fold_", ""),
        sample_index = sample_index |> as.numeric(),
        pDockQ2_Poc4_Pup2 = pDockQ2)

save(AFS_20251224, file = "intermediate_data/AlphaFold3/AFS_20251224.Rdata")



AFS_20251229 <- submit_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_20251229/seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE) |>
    dplyr::filter(
        Chn1 == "D", #Pup2
        Chn2 == "K", #Poc4
        Type == "max") |>
    dplyr::transmute(
        dataset_id = "AFS_2025129",
        structure_id = structure_id |> stringr::str_replace("fold_", ""),
        sample_index = sample_index |> as.numeric(),
        pDockQ2_Poc4_Pup2 = pDockQ2)


AFS_20251229 <- compute_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_20251229/seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE) |>
    dplyr::filter(
        Chn1 == "D", #Pup2
        Chn2 == "K", #Poc4
        Type == "max") |>
    dplyr::transmute(
        dataset_id = "AFS_2025129",
        structure_id = structure_id |> stringr::str_replace("fold_", ""),
        sample_index = sample_index |> as.numeric(),
        pDockQ2_Poc4_Pup2 = pDockQ2)

save(AFS_20251229, file = "intermediate_data/AlphaFold3/AFS_20251229.Rdata")


submit_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_set1_20251229/seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE)

submit_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_set2_20251229/seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE)

submit_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_set3_20251229/seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE)

submit_all_dockq_afs(
     results_path = Sys.glob(
         path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_set4_20251229/seq*"),
     pae_cutoff = 15,
     dist_cutoff = 15,
     verbose = TRUE)

submit_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_set5_20251229/seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE)

submit_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_set6_20251229/seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE)

submit_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_set7_20251231/seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE)


submit_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_set8_20251231/seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE)


submit_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_set9_20251231/seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE)


submit_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_set10_20251231/seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE)


submit_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_set_extra1_20251231/seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE)


submit_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_set_extra2_20251231/seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE)


submit_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_set11-12_20251230/seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE)

submit_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_set13_20251230/seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE)


submit_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_results_set16_20251231/seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE)


afs_dockq <- compute_all_dockq_afs(
    results_path = c(
        Sys.glob(
            path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_*/seq*"),
        Sys.glob(
            path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/af3_afs_*/fold_seq*")),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE) |>
    dplyr::filter(
        ifelse(
            Model |> stringr::str_detect("af3_afs_results_20251003"),
            Chn1 == "E",
            Chn1 == "D"), #Pup2
        Chn2 == "K", #Poc4
        Type == "max") |>
    dplyr::mutate(
        #Fix path names due to renaming folders after computing DockQ2 scores
        Model = Model |>
            stringr::str_replace(
                "af3_afs_results_20251229",
                "af3_afs_results_set1_20251229")) |>
    dplyr::transmute(
        dataset_id = "AFS",
        cif_fname = paste0(Model, ".cif"),
        structure_id = structure_id |> stringr::str_extract("seq_[0-9]+"),
        sample_index = sample_index |> as.numeric(),
        pDockQ2_Poc4_Pup2 = pDockQ2,
        Pup2_sequence = Seq1,
        Poc4_sequence = Seq2) |>
    dplyr::arrange(dplyr::desc(pDockQ2_Poc4_Pup2)) |>
    dplyr::distinct(Poc4_sequence, .keep_all = TRUE)


#check cif files exist
afs_dockq |>
    dplyr::rowwise() |>
    dplyr::do({
        x <- .
        if (!file.exists(x$cif_fname)) {
            cat("WARNING: cif file '", x$cif_fname, "' does not extst.\n", sep = "")
        }
        data.frame()
    })

save(afs_dockq, file = "intermediate_data/AlphaFold3/afs_dockq.Rdata")

afs_dockq |>
    readr::write_tsv(
        file = "intermediate_data/AlphaFold3/afs_dockq.tsv")


afs_dockq |>
    dplyr::mutate(
        fasta = paste0("> ", structure_id, "-", sample_index, "\n", Poc4_sequence, "\n")) |>
    dplyr::pull(fasta) |>
    cat(sep = "", file = "intermediate_data/AlphaFold3/afs_dockq.fasta")


###############################

## Test set
submit_all_dockq_afs(
    results_path = Sys.glob(
        path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/afs_outputs_january/*/seq*"),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE)



afs_dockq_test <- compute_all_dockq_afs(
    results_path = c(
        Sys.glob(
            path = "/home/maom/turbo_tromeara/CauProteosome/intermediate_data/AF3/afs_outputs_january/*/seq*")),
    pae_cutoff = 15,
    dist_cutoff = 15,
    verbose = TRUE) |>
    dplyr::filter(
        Chn1 == "D", #Pup2
        Chn2 == "K", #Poc4
        Type == "max") |>
    dplyr::transmute(
        dataset_id = "AFS_test",
        cif_fname = paste0(Model, ".cif"),
        structure_id = structure_id |> stringr::str_extract("seq_[0-9]+"),
        sample_index = sample_index |> as.numeric(),
        pDockQ2_Poc4_Pup2 = pDockQ2,
        Pup2_sequence = Seq1,
        Poc4_sequence = Seq2) |>
    dplyr::arrange(dplyr::desc(pDockQ2_Poc4_Pup2)) |>
    dplyr::distinct(Poc4_sequence, .keep_all = TRUE)


#check cif files exist
afs_dockq_test |>
    dplyr::rowwise() |>
    dplyr::do({
        x <- .
        if (!file.exists(x$cif_fname)) {
            cat("WARNING: cif file '", x$cif_fname, "' does not extst.\n", sep = "")
        }
        data.frame()
    })

save(afs_dockq_test, file = "intermediate_data/AlphaFold3/afs_dockq_test.Rdata")

afs_dockq_test |>
    readr::write_tsv(
        file = "intermediate_data/AlphaFold3/afs_dockq_test.tsv")


afs_dockq_test |>
    dplyr::mutate(
        fasta = paste0("> ", structure_id, "-", sample_index, "\n", Poc4_sequence, "\n")) |>
    dplyr::pull(fasta) |>
    cat(sep = "", file = "intermediate_data/AlphaFold3/afs_dockq_test.fasta")
