source('../haplotype_analysis.R')

bc_to_sample = read_excel("sample_list.xlsx", col_types = 'numeric') %>%
  arrange(
    `Ratio NF54`, `replicate`
  ) %>%
  mutate(
    Sample = factor(sprintf("NF54_%03.1f__%1.0f", `Ratio NF54`*100, `replicate`), levels = sprintf("NF54_%03.1f__%1.0f", `Ratio NF54`*100, `replicate`)),
    Barcode = as.character(barcode),
  )

dhps = init_phare(
  folder = 'export/dhps',
  pattern = '*.tsv',
  samples = bc_to_sample, # a DATA FRAME! with one column containing `Source` (== barcode) and another titeled `Sample`
  fasta_ref = '../templates/dhps.fasta',
  gb_report = '../templates/dhps.jsonl',
  gene_name = 'DHPS',
  threshold = 0.05,
  filter_silent = TRUE,
  min_reads = 50
) %>%
  make_snp_df() %>%
  filter_dupl_snps() %>%
  phare_translate() %>%
  phare_frequency() %>%
  apply_threshold() %>%
  make_plot_df() %>%
  phare_plot() %>%
  save_plot(sup='')