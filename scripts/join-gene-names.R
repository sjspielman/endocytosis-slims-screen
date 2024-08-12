#!/usr/bin/env Rscript

dat <- readr::read_csv("temp.csv")
genes <- readr::read_csv("../reference/human_uniprot_genenames.csv")

dat |>
    dplyr::inner_join(
        genes,
        by = c("uniprot_id" = "sequence_id")
    ) |>
    dplyr::select(gene_name, uniprot_id, motif_sequence, start_position) |>
    dplyr::arrange(gene_name) |>
    readr::write_csv("../results/motif-results.csv")