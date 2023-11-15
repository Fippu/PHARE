source('workflow/scripts/haplotype_analysis.R')




# PARSE COMMAND LINE ARGUMENTS
library(optparse)

option_list <- list(
  make_option("--folder", type="character", help="Specify the folder"),
  make_option("--pattern", type="character", help="Specify the pattern"),
  # make_option("--samples", type="character", help="Specify the path to the samples data frame"),
  make_option("--fasta_ref", type="character", help="Specify the path to the FASTA reference file"),
  make_option("--gb_report", type="character", help="Specify the path to the GenBank report file"),
  make_option("--gene_name", type="character", help="Specify the gene name"),
  make_option("--threshold", type="numeric", help="Specify the threshold value"),
  make_option("--filter_silent", type="numeric", help="Filter silent option [default]"),
  make_option("--min_reads", type="numeric", help="Specify the minimum number of reads")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser)
# print(args)



bc_to_sample = read_excel("config/sample_list.xlsx", col_types = 'text') %>%
  # arrange(
  #   `Ratio NF54`, `replicate`
  # ) %>%
  mutate(
    Sample = factor(as.character(sample)),
    # Sample = factor(sprintf("NF54_%03.1f__%1.0f", `Ratio NF54`*100, `replicate`), levels = sprintf("NF54_%03.1f__%1.0f", `Ratio NF54`*100, `replicate`)),
    Barcode = as.character(barcode),
  )


print(bc_to_sample$Sample)
# ANALYZE THE DATA AND GENERATE A PLOT
gene = init_phare(
  folder = args$folder,
  pattern = args$pattern,
  samples = bc_to_sample, # a DATA FRAME! with one column containing `Source` (== barcode) and another titeled `Sample`
  fasta_ref = args$fasta_ref,
  gb_report = args$gb_report,
  gene_name = args$gene_name,
  threshold = args$threshold,
  filter_silent = args$filter_silent,
  min_reads = args$min_reads
) %>%
  make_snp_df() %>%
  filter_dupl_snps() %>%
  phare_translate() %>%
  phare_frequency() %>%
  apply_threshold() %>%
  make_plot_df() %>%
  phare_plot() %>%
  save_plot(dir='results', sup='')
