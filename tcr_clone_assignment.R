#!/usr/bin/R
# =============================================================================
# TCR Clone Assignment from Cell Ranger VDJ Output
# Handles dual alpha/beta chains per cell barcode
# Author: Hermes Agent (for Thiago Y Oliveira)
# Date: 2026-07-02
# =============================================================================

library(tidyverse)
library(igraph)

# =============================================================================
# 1. LOAD AND FILTER CONTIGS
# =============================================================================

load_contigs <- function(all_contig_csv,
                         filtered_contig_csv = NULL,
                         min_umi = 5,
                         min_read = 20) {
  message(">>> Loading all_contig_annotations.csv...")

  contigs <- read.csv(all_contig_csv, stringsAsFactors = FALSE)

  # Basic filtering: cell, high confidence, productive, TRA/TRB only
  contigs <- contigs %>%
    filter(is_cell == TRUE) %>%
    filter(high_confidence == TRUE) %>%
    filter(productive == TRUE) %>%
    filter(chain %in% c("TRA", "TRB"))

  # Quality filter: minimum UMI and read counts
  # Note: adjust column names — Cell Ranger uses 'umis' and 'reads'
  # Some versions use 'umi_count' and 'read_count'
  if ("umis" %in% colnames(contigs)) {
    contigs <- contigs %>% filter(umis >= min_umi)
  } else if ("umi_count" %in% colnames(contigs)) {
    contigs <- contigs %>% rename(umis = umi_count) %>% filter(umis >= min_umi)
  }

  if ("reads" %in% colnames(contigs)) {
    contigs <- contigs %>% filter(reads >= min_read)
  } else if ("read_count" %in% colnames(contigs)) {
    contigs <- contigs %>% rename(reads = read_count) %>% filter(reads >= min_read)
  }

  message(sprintf("  %d productive contigs from %d unique cell barcodes",
                  nrow(contigs), length(unique(contigs$barcode))))

  return(contigs)
}


# =============================================================================
# 2. RANK CHAINS PER CELL BARCODE
# =============================================================================

rank_chains <- function(contigs) {
  message(">>> Ranking chains by UMI per barcode...")

  # For each barcode + chain type, rank by UMI (descending)
  ranked <- contigs %>%
    group_by(barcode, chain) %>%
    mutate(rank = row_number(desc(umis))) %>%
    ungroup()

  # Separate primary (rank 1) and secondary (rank 2+)
  primary <- ranked %>%
    filter(rank == 1) %>%
    select(barcode, chain, v_gene, d_gene, j_gene, c_gene,
           cdr3, cdr3_nt, umis, reads) %>%
    rename_with(~ paste0(.x, "_prim"), -c(barcode, chain)) %>%
    rename(chain_prim = chain)

  secondary <- ranked %>%
    filter(rank == 2) %>%
    select(barcode, chain, v_gene, d_gene, j_gene, c_gene,
           cdr3, cdr3_nt, umis, reads) %>%
    rename_with(~ paste0(.x, "_sec"), -c(barcode, chain)) %>%
    rename(chain_sec = chain)

  # Count chains per barcode
  n_chains <- contigs %>%
    group_by(barcode, chain) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = chain, values_from = n,
                names_prefix = "n_", values_fill = 0)

  # Flag dual-chain cells
  n_chains <- n_chains %>%
    mutate(
      dual_alpha = n_TRA >= 2,
      dual_beta  = n_TRB >= 2,
      chain_status = case_when(
        dual_alpha & dual_beta  ~ "dual_alpha_beta",
        dual_alpha              ~ "dual_alpha",
        dual_beta               ~ "dual_beta",
        TRUE                    ~ "single"
      )
    )

  message(sprintf("  Chain status breakdown:"))
  print(table(n_chains$chain_status))

  return(list(
    primary   = primary,
    secondary = secondary,
    n_chains  = n_chains
  ))
}


# =============================================================================
# 3. FILTER DUAL TCR BY QUALITY THRESHOLD
# =============================================================================

filter_dual_tcr <- function(primary, secondary, min_ratio = 0.1, min_abs_umi = 5) {
  message(sprintf(">>> Filtering dual-TCR cells (min_ratio=%.0f%%, min_abs_umi=%d)...",
                  min_ratio * 100, min_abs_umi))

  # Join primary and secondary by barcode + chain
  merged <- secondary %>%
    inner_join(
      primary %>% select(barcode, chain_prim, umis_prim),
      by = "barcode"
    ) %>%
    # Only compare same chain type (TRA prim with TRA sec, TRB prim with TRB sec)
    filter(chain_sec == chain_prim) %>%
    mutate(
      umi_ratio = umis_sec / umis_prim,
      trusted_dual = umi_ratio >= min_ratio & umis_sec >= min_abs_umi
    )

  n_trusted <- sum(merged$trusted_dual)
  n_total   <- nrow(merged)

  message(sprintf("  %d / %d secondary chains pass quality filter (%.1f%%)",
                  n_trusted, n_total, n_trusted / n_total * 100))

  return(merged)
}


# =============================================================================
# 4. BUILD CLONE TABLE
# =============================================================================

build_clones <- function(primary, secondary_filtered, clone_definition = "TRB") {
  message(sprintf(">>> Building clone table (definition: %s)...", clone_definition))

  # Wide format: one row per barcode, columns split by chain
  prim_wide <- primary %>%
    pivot_wider(
      id_cols    = barcode,
      names_from = chain_prim,
      values_from = c(v_gene_prim, d_gene_prim, j_gene_prim, c_gene_prim,
                      cdr3_prim, cdr3_nt_prim, umis_prim, reads_prim),
      values_fill = NA
    )

  # Build clone ID based on definition
  if (clone_definition == "TRB") {
    # Clone = TRB CDR3 amino acid + V + J gene (standard)
    prim_wide <- prim_wide %>%
      mutate(clone_key = paste(
        cdr3_prim_TRB %||% "NA",
        v_gene_prim_TRB %||% "NA",
        j_gene_prim_TRB %||% "NA",
        sep = "|"
      ))
  } else if (clone_definition == "TRA_TRB") {
    # Clone = TRA + TRB CDR3 combo (more stringent)
    prim_wide <- prim_wide %>%
      mutate(clone_key = paste(
        cdr3_prim_TRA %||% "NA", v_gene_prim_TRA %||% "NA",
        cdr3_prim_TRB %||% "NA", v_gene_prim_TRB %||% "NA",
        sep = "|"
      ))
  } else if (clone_definition == "TRB_cdr3_only") {
    # Clone = TRB CDR3 aa only (most permissive)
    prim_wide <- prim_wide %>%
      mutate(clone_key = cdr3_prim_TRB %||% "NA")
  }

  # Assign clone IDs (ordered by frequency)
  clone_freq <- prim_wide %>%
    count(clone_key, sort = TRUE) %>%
    mutate(clone_id = paste0("clon", row_number()))

  prim_wide <- prim_wide %>%
    left_join(clone_freq, by = "clone_key") %>%
    select(-clone_key)

  # Add dual-chain annotation
  prim_wide <- prim_wide %>%
    mutate(
      has_dual_alpha = barcode %in% secondary_filtered$barcode[secondary_filtered$chain_sec == "TRA" & secondary_filtered$trusted_dual],
      has_dual_beta  = barcode %in% secondary_filtered$barcode[secondary_filtered$chain_sec == "TRB" & secondary_filtered$trusted_dual]
    )

  message(sprintf("  %d clonotypes from %d cells (top clone: %d cells)",
                  nrow(clone_freq), nrow(prim_wide), clone_freq$n[1]))

  return(list(
    cells       = prim_wide,
    clonotypes  = clone_freq
  ))
}

# Null-coalescing helper
`%||%` <- function(a, b) if (is.null(a) || all(is.na(a))) b else a


# =============================================================================
# 5. QC METRICS AND SUMMARY
# =============================================================================

qc_summary <- function(cells_df, n_chains, clonotypes) {
  message(">>> QC Summary")
  message("")

  n_total_cells <- nrow(cells_df)

  # Chain pairing stats
  message("=== Chain pairing ===")
  n_both  <- sum(!is.na(cells_df$cdr3_prim_TRA) & !is.na(cells_df$cdr3_prim_TRB))
  n_alpha_only <- sum(!is.na(cells_df$cdr3_prim_TRA) & is.na(cells_df$cdr3_prim_TRB))
  n_beta_only  <- sum(is.na(cells_df$cdr3_prim_TRA) & !is.na(cells_df$cdr3_prim_TRB))

  message(sprintf("  Cells with alpha+beta:   %d (%.1f%%)", n_both, n_both/n_total_cells*100))
  message(sprintf("  Cells with alpha only:   %d (%.1f%%)", n_alpha_only, n_alpha_only/n_total_cells*100))
  message(sprintf("  Cells with beta only:    %d (%.1f%%)", n_beta_only, n_beta_only/n_total_cells*100))
  message("")

  # Dual-chain stats
  message("=== Dual-chain prevalence ===")
  dual_stats <- table(n_chains$chain_status)
  for (status in names(dual_stats)) {
    message(sprintf("  %-16s  %d (%.1f%%)", status,
                    dual_stats[status], dual_stats[status]/sum(dual_stats)*100))
  }
  message("")

  # Expected: ~30% dual-alpha, ~5-10% dual-beta
  pct_dual_alpha <- mean(n_chains$dual_alpha) * 100
  message(sprintf("  Dual-alpha rate: %.1f%% (expected ~25-35%%)", pct_dual_alpha))
  if (pct_dual_alpha > 50) {
    message("  WARNING: Dual-alpha rate >50%% — possible doublet contamination!")
  }
  message("")

  # Clonotype stats
  message("=== Clonotype diversity ===")
  n_clones <- nrow(clonotypes)
  message(sprintf("  Total clonotypes:     %d", n_clones))
  message(sprintf("  Singletons (n=1):     %d (%.1f%%)",
                  sum(clonotypes$n == 1), sum(clonotypes$n == 1)/n_clones*100))
  message(sprintf("  Top 10 clone share:   %.1f%% of cells",
                  sum(head(clonotypes$n, 10)) / sum(clonotypes$n) * 100))

  # Shannon diversity
  props <- clonotypes$n / sum(clonotypes$n)
  shannon <- -sum(props * log(props))
  pielou   <- shannon / log(n_clones)
  message(sprintf("  Shannon diversity:    %.2f", shannon))
  message(sprintf("  Pielou evenness:      %.2f", pielou))
  message("")

  return(invisible(list(
    n_cells       = n_total_cells,
    n_clones      = n_clones,
    shannon       = shannon,
    pielou        = pielou,
    dual_alpha_pct = pct_dual_alpha,
    dual_beta_pct  = mean(n_chains$dual_beta) * 100
  )))
}


# =============================================================================
# 6. PLOTTING FUNCTIONS
# =============================================================================

plot_clone_size_dist <- function(clonotypes, top_n = 50) {
  ggplot(head(clonotypes, top_n), aes(x = reorder(clone_id, n), y = n)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    coord_flip() +
    scale_y_log10() +
    labs(
      title    = sprintf("Top %d TCR Clonotypes", top_n),
      subtitle = "Clone size (log scale)",
      x        = "Clonotype",
      y        = "Number of cells"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

plot_chain_status <- function(n_chains) {
  status_df <- as.data.frame(table(n_chains$chain_status))
  colnames(status_df) <- c("status", "count")

  ggplot(status_df, aes(x = "", y = count, fill = status)) +
    geom_col(width = 1) +
    coord_polar("y") +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "TCR Chain Status per Cell Barcode",
      fill  = "Status"
    ) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

plot_dual_umi_ratio <- function(secondary_filtered) {
  ggplot(secondary_filtered, aes(x = umi_ratio, fill = trusted_dual)) +
    geom_histogram(bins = 50, alpha = 0.8, color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0.1, linetype = "dashed", color = "red") +
    scale_x_log10(labels = scales::percent) +
    scale_fill_manual(values = c(`TRUE` = "#2ca25f", `FALSE` = "#d73027"),
                      labels = c(`TRUE` = "Trusted dual", `FALSE` = "Filtered out")) +
    labs(
      title    = "UMI Ratio: Secondary / Primary Chain",
      subtitle = "Red dashed line = 10% threshold",
      x        = "UMI ratio (log scale)",
      y        = "Count",
      fill     = ""
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}


# =============================================================================
# 7. scRepertoire COMPATIBILITY (optional)
# =============================================================================

to_screpertoire <- function(contigs) {
  message(">>> Formatting for scRepertoire compatibility...")

  # scRepertoire expects: barcode, chain, v_gene, d_gene, j_gene, c_gene,
  #                       cdr3, cdr3_nt, reads, umis, productive, is_cell, etc.
  df <- contigs %>%
    select(
      barcode, chain,
      v_gene, d_gene, j_gene, c_gene,
      cdr3, cdr3_nt,
      reads, umis,
      productive, is_cell, high_confidence
    ) %>%
    mutate(sample = "sample1")  # adjust if multiple samples

  message(sprintf("  Output: %d contigs, ready for scRepertoire::combineTCR()"))
  return(df)
}


# =============================================================================
# 7b. BATCH LOADING (multiple samples — barcode de-duplication)
# =============================================================================

load_contigs_batch <- function(sample_list,
                               min_umi = 5,
                               min_read = 20) {
  # sample_list: named list of paths to all_contig_annotations.csv
  #   e.g. list(CTRL = "ctrl/vdj_t/outs/all_contig_annotations.csv",
  #             PATIENT = "patient/vdj_t/outs/all_contig_annotations.csv")
  #
  # Returns a single data.frame with barcodes prefixed by sample ID
  # to avoid collisions when multiple samples share the same 10x whitelist.

  message(">>> Loading and combining ", length(sample_list), " samples...")

  all_contigs <- lapply(names(sample_list), function(sample_id) {
    message("  - ", sample_id, ": ", sample_list[[sample_id]])

    df <- read.csv(sample_list[[sample_id]], stringsAsFactors = FALSE)

    # Basic filtering
    df <- df %>%
      filter(is_cell == TRUE) %>%
      filter(high_confidence == TRUE) %>%
      filter(productive == TRUE) %>%
      filter(chain %in% c("TRA", "TRB"))

    # Normalize UMI/read column names
    if ("umis" %in% colnames(df)) {
      df <- df %>% filter(umis >= min_umi)
    } else if ("umi_count" %in% colnames(df)) {
      df <- df %>% rename(umis = umi_count) %>% filter(umis >= min_umi)
    }

    if ("reads" %in% colnames(df)) {
      df <- df %>% filter(reads >= min_read)
    } else if ("read_count" %in% colnames(df)) {
      df <- df %>% rename(reads = read_count) %>% filter(reads >= min_read)
    }

    # Prefix barcode with sample ID to prevent collisions across samples
    # 10x barcodes come from a shared ~3.4M whitelist — same sequence in
    # different runs = different physical cells
    df$barcode <- paste0(sample_id, "_", df$barcode)
    df$sample  <- sample_id

    df
  }) |> dplyr::bind_rows()

  message(sprintf("  Combined: %d productive contigs from %d unique cell barcodes",
                  nrow(all_contigs), length(unique(all_contigs$barcode))))

  return(all_contigs)
}


# =============================================================================
# 7c. CORE FUNCTION — assign clones to a data.frame (no files, no plots)
# =============================================================================

assign_clones <- function(contigs,
                          clone_def   = "TRB",
                          min_ratio   = 0.1,
                          min_abs_umi = 5) {
  # Core function: takes filtered contigs (from load_contigs or load_contigs_batch),
  # assigns clonotypes, returns a single data.frame.
  # No files written, no plots generated. Pure transformation.
  #
  # Input:  contigs data.frame with columns: barcode, chain, v_gene, j_gene,
  #         cdr3, cdr3_nt, umis, productive, high_confidence, is_cell
  # Output: data.frame with one row per barcode + clone_id column

  ranked        <- rank_chains(contigs)
  dual_filtered <- filter_dual_tcr(ranked$primary, ranked$secondary,
                                   min_ratio = min_ratio, min_abs_umi = min_abs_umi)
  clones        <- build_clones(ranked$primary, dual_filtered,
                                clone_definition = clone_def)

  return(clones$cells)
}


# =============================================================================
# 8. MAIN PIPELINE (single sample)
# =============================================================================

run_tcr_pipeline <- function(
    all_contig_csv,
    output_dir     = "tcr_results",
    clone_def      = "TRB",
    min_umi        = 5,
    min_read       = 20,
    min_ratio      = 0.1,
    min_abs_umi    = 5,
    make_plots     = TRUE
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  message("=" |> strrep(70))
  message("  TCR Clone Assignment Pipeline")
  message("  Clone definition: ", clone_def)
  message("=" |> strrep(70))
  message("")

  # Step 1: Load
  contigs <- load_contigs(all_contig_csv,
                          min_umi = min_umi,
                          min_read = min_read)

  # Step 2: Rank chains
  ranked <- rank_chains(contigs)

  # Step 3: Filter dual-TCR
  dual_filtered <- filter_dual_tcr(ranked$primary,
                                   ranked$secondary,
                                   min_ratio   = min_ratio,
                                   min_abs_umi = min_abs_umi)

  # Step 4: Build clones
  clones <- build_clones(ranked$primary, dual_filtered,
                         clone_definition = clone_def)

  # Step 5: QC
  qc <- qc_summary(clones$cells, ranked$n_chains, clones$clonotypes)

  # Step 6: Plots
  if (make_plots) {
    message(">>> Generating plots...")

    p1 <- plot_clone_size_dist(clones$clonotypes)
    p2 <- plot_chain_status(ranked$n_chains)
    p3 <- plot_dual_umi_ratio(dual_filtered)

    ggsave(file.path(output_dir, "clone_size_distribution.pdf"), p1,
           width = 8, height = 10)
    ggsave(file.path(output_dir, "chain_status_pie.pdf"), p2,
           width = 6, height = 6)
    ggsave(file.path(output_dir, "dual_umi_ratio.pdf"), p3,
           width = 8, height = 5)
  }

  # Step 7: Save outputs
  message(">>> Saving outputs...")
  write.csv(clones$cells,
            file.path(output_dir, "tcr_cells_with_clones.csv"), row.names = FALSE)
  write.csv(clones$clonotypes,
            file.path(output_dir, "clonotype_table.csv"), row.names = FALSE)
  write.csv(dual_filtered,
            file.path(output_dir, "dual_chain_filtered.csv"), row.names = FALSE)
  write.csv(to_screpertoire(contigs),
            file.path(output_dir, "contigs_screpertoire_format.csv"), row.names = FALSE)

  message("")
  message("Done! Outputs saved to: ", output_dir, "/")
  message("  - tcr_cells_with_clones.csv   (main table: barcode x clonotype)")
  message("  - clonotype_table.csv          (clone frequencies)")
  message("  - dual_chain_filtered.csv      (trusted dual-TCR annotations)")
  message("  - contigs_screpertoire_format.csv (for scRepertoire)")

  return(list(
    contigs         = contigs,
    ranked          = ranked,
    dual_filtered   = dual_filtered,
    cells           = clones$cells,
    clonotypes      = clones$clonotypes,
    qc              = qc
  ))
}

# =============================================================================
# 8b. MAIN PIPELINE (batch — multiple samples combined)
# =============================================================================

run_tcr_pipeline_batch <- function(
    sample_list,
    output_dir     = "tcr_results_batch",
    clone_def      = "TRB",
    min_umi        = 5,
    min_read       = 20,
    min_ratio      = 0.1,
    min_abs_umi    = 5,
    make_plots     = TRUE
) {
  # sample_list: named list where names = sample IDs, values = paths to CSVs
  #   e.g. list(CTRL = "ctrl/outs/all_contig_annotations.csv",
  #             PATIENT = "patient/outs/all_contig_annotations.csv")
  #
  # Combines all samples into a single analysis with unique barcodes
  # (prefixed by sample ID), enabling cross-sample clonotype comparison.

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  message("=" |> strrep(70))
  message("  TCR Clone Assignment Pipeline — BATCH MODE")
  message("  Samples: ", paste(names(sample_list), collapse = ", "))
  message("  Clone definition: ", clone_def)
  message("=" |> strrep(70))
  message("")

  # Step 1: Load and combine samples with unique barcodes
  contigs <- load_contigs_batch(sample_list,
                                min_umi  = min_umi,
                                min_read = min_read)

  # Step 2: Rank chains (operates on combined data)
  ranked <- rank_chains(contigs)

  # Step 3: Filter dual-TCR
  dual_filtered <- filter_dual_tcr(ranked$primary,
                                   ranked$secondary,
                                   min_ratio   = min_ratio,
                                   min_abs_umi = min_abs_umi)

  # Step 4: Build clones (combined)
  clones <- build_clones(ranked$primary, dual_filtered,
                         clone_definition = clone_def)

  # Step 5: QC
  qc <- qc_summary(clones$cells, ranked$n_chains, clones$clonotypes)

  # Step 5b: Per-sample breakdown
  message("")
  message(">>> Per-sample clonotype summary:")
  per_sample <- clones$cells %>%
    group_by(sample) %>%
    summarise(
      n_cells    = n(),
      n_clones   = n_distinct(clone_id),
      top_clone_n = max(table(clone_id)),
      top_clone_pct = round(max(table(clone_id)) / n() * 100, 1),
      .groups = "drop"
    ) %>%
    arrange(desc(n_cells))

  print(per_sample)

  # Step 6: Cross-sample overlap analysis
  message("")
  message(">>> Cross-sample clonotype overlap:")
  overlap <- clones$cells %>%
    distinct(sample, clone_id) %>%
    group_by(clone_id) %>%
    summarise(
      n_samples = n_distinct(sample),
      samples   = paste(sort(unique(sample)), collapse = ", "),
      .groups   = "drop"
    ) %>%
    filter(n_samples >= 2) %>%
    arrange(desc(n_samples))

  message(sprintf("  %d clonotypes shared across 2+ samples",
                  nrow(overlap)))
  if (nrow(overlap) > 0) {
    message("  Top shared clonotypes:")
    print(head(overlap, 10))
  }

  # Step 7: Plots
  if (make_plots) {
    message("")
    message(">>> Generating plots...")

    # Standard plots
    p1 <- plot_clone_size_dist(clones$clonotypes)
    p2 <- plot_chain_status(ranked$n_chains)
    p3 <- plot_dual_umi_ratio(dual_filtered)

    ggsave(file.path(output_dir, "clone_size_distribution.pdf"), p1,
           width = 8, height = 10)
    ggsave(file.path(output_dir, "chain_status_pie.pdf"), p2,
           width = 6, height = 6)
    ggsave(file.path(output_dir, "dual_umi_ratio.pdf"), p3,
           width = 8, height = 5)

    # Batch-specific: per-sample clone counts
    p4 <- ggplot(per_sample, aes(x = reorder(sample, -n_cells), y = n_cells)) +
      geom_col(fill = "steelblue", alpha = 0.8) +
      labs(
        title    = "Cells per Sample",
        x        = "Sample",
        y        = "Number of cells"
      ) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))

    # Batch-specific: clonotype diversity per sample
    p5 <- ggplot(per_sample, aes(x = reorder(sample, -n_clones), y = n_clones)) +
      geom_col(fill = "darkorange", alpha = 0.8) +
      labs(
        title    = "Unique Clonotypes per Sample",
        x        = "Sample",
        y        = "Number of clonotypes"
      ) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))

    # Batch-specific: top clone expansion per sample
    p6 <- ggplot(per_sample, aes(x = reorder(sample, -top_clone_pct), y = top_clone_pct)) +
      geom_col(fill = "#d73027", alpha = 0.8) +
      labs(
        title    = "Largest Clone per Sample (% of cells)",
        x        = "Sample",
        y        = "Top clone size (%)"
      ) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))

    ggsave(file.path(output_dir, "batch_cells_per_sample.pdf"), p4,
           width = 6, height = 4)
    ggsave(file.path(output_dir, "batch_clonotypes_per_sample.pdf"), p5,
           width = 6, height = 4)
    ggsave(file.path(output_dir, "batch_top_clone_expansion.pdf"), p6,
           width = 6, height = 4)
  }

  # Step 8: Save outputs
  message("")
  message(">>> Saving outputs...")
  write.csv(clones$cells,
            file.path(output_dir, "tcr_cells_with_clones.csv"), row.names = FALSE)
  write.csv(clones$clonotypes,
            file.path(output_dir, "clonotype_table.csv"), row.names = FALSE)
  write.csv(dual_filtered,
            file.path(output_dir, "dual_chain_filtered.csv"), row.names = FALSE)
  write.csv(to_screpertoire(contigs),
            file.path(output_dir, "contigs_screpertoire_format.csv"), row.names = FALSE)
  write.csv(per_sample,
            file.path(output_dir, "per_sample_summary.csv"), row.names = FALSE)
  write.csv(overlap,
            file.path(output_dir, "cross_sample_overlap.csv"), row.names = FALSE)

  message("")
  message("Done! Outputs saved to: ", output_dir, "/")
  message("  - tcr_cells_with_clones.csv       (main table: barcode x clonotype)")
  message("  - clonotype_table.csv              (clone frequencies)")
  message("  - dual_chain_filtered.csv          (trusted dual-TCR annotations)")
  message("  - contigs_screpertoire_format.csv  (for scRepertoire)")
  message("  - per_sample_summary.csv           (per-sample stats)")
  message("  - cross_sample_overlap.csv         (shared clonotypes)")

  return(list(
    contigs       = contigs,
    ranked        = ranked,
    dual_filtered = dual_filtered,
    cells         = clones$cells,
    clonotypes    = clones$clonotypes,
    per_sample    = per_sample,
    overlap       = overlap,
    qc            = qc
  ))
}


# =============================================================================
# 9. USAGE
# =============================================================================

# Single sample:
# results <- run_tcr_pipeline(
#   all_contig_csv = "vdj_t/outs/all_contig_annotations.csv",
#   output_dir     = "tcr_results",
#   clone_def      = "TRB",     # "TRB", "TRA_TRB", or "TRB_cdr3_only"
#   min_umi        = 5,
#   min_read       = 20,
#   min_ratio      = 0.1,       # secondary chain must be ≥10% of primary UMI
#   min_abs_umi    = 5,
#   make_plots     = TRUE
# )
#
# Multiple samples — BATCH MODE (cross-sample comparison):
# samples <- list(
#   CTRL    = "ctrl/vdj_t/outs/all_contig_annotations.csv",
#   PATIENT = "patient/vdj_t/outs/all_contig_annotations.csv"
# )
# results <- run_tcr_pipeline_batch(
#   sample_list = samples,
#   output_dir  = "tcr_results_batch",
#   clone_def   = "TRB"
# )
#
# # Per-sample independent processing (no cross-sample comparison):
# # all_results <- lapply(names(samples), function(s) {
# #   run_tcr_pipeline(samples[[s]], output_dir = paste0("tcr_results/", s))
# # })
# # names(all_results) <- names(samples)
#
# # For scRepertoire integration:
# library(scRepertoire)
# combined <- combineTCR(
#   list(results$contigs %>% mutate(sample = "sample1")),
#   samples = NULL,
#   ID      = NULL
# )
# combined <- clonalHomeostasis(combined)
# combined <- clonalProportion(combined)
# combined <- clonalOverlap(combined, method = "morisita")
