library(readxl)
library(seqinr)
library(reshape2)
library(rjson)
library(tidyverse)
library(ggpubr)
library(data.table)



# _______________________________SETUP___________________________________



# examples:
# setwd("C:/Users/uname/Documents/SNPs")


# dhps = init_phare(
#   folder = 'car/dhps',
#   pattern = '*.tsv',
#   samples = bc_to_sample_car, # a DATA FRAME! with one column containing `Source` (== barcode) and another titeled `Sample`
#   fasta_ref = 'dhps.fasta',
#   gb_report = 'dhps.jsonl',
#   gene_name = 'DHPS',
#   threshold = 0.09,
#   filter_silent = TRUE,
#   min_reads = 50
# ) %>%
#   make_snp_df() %>%
#   filter_dupl_snps() %>%
#   phare_translate() %>%
#   phare_frequency() %>%
#   apply_threshold() %>%
#   make_plot_df() %>%
#   phare_plot() %>%
#   save_plot(sup='minqual10')
# 
# dhps$fig
# 
# 
# 
# dhfr = init_phare(
#   folder = 'car/dhfr',
#   pattern = '*.tsv',
#   samples = bc_to_sample_car, # a DATA FRAME! with one column containing `Source` (== barcode) and another titeled `Sample`
#   fasta_ref = 'dhfr.fasta',
#   gb_report = 'dhfr.jsonl',
#   gene_name = 'DHFR',
#   threshold = 0.09,
#   filter_silent = TRUE,
#   min_reads = 100
# ) %>%
#   make_snp_df() %>%
#   filter_dupl_snps() %>%
#   phare_translate() %>%
#   phare_frequency() %>%
#   apply_threshold() %>%
#   make_plot_df() %>%
#   phare_plot() %>%
#   save_plot(sup='minqual10')
# dhfr$fig
# 
# 
# 
# k13 = init_phare(
#   folder = 'car/k13',
#   pattern = '*.tsv',
#   samples = bc_to_sample_car, # a DATA FRAME! with one column containing `Source` (== barcode) and another titeled `Sample`
#   fasta_ref = 'k13.fasta',
#   gb_report = 'k13.jsonl',
#   gene_name = 'K13',
#   threshold = 0.09,
#   filter_silent = TRUE,
#   min_reads = 100
# ) %>%
#   make_snp_df() %>%
#   filter_dupl_snps() %>%
#   phare_translate() %>%
#   phare_frequency() %>%
#   apply_threshold() %>%
#   make_plot_df() %>%
#   phare_plot() %>%
#   save_plot(sup='minqual10')







# _______________________________FUNCTIONS________________________________________________________________________________________________________


# A FUNCTION TO DO THE INITIAL SETUP. ALL VARIABLES ARE STORED IN A LIST FOR EASY ACCESS (and to simulate a kind of class like structure)
init_phare = function(
    folder, # a folder with input reads
    pattern = '*.tsv',
    samples, # a DATA FRAME! with one column containing `Barcode` (NEEDS A LEADING 0 IF A BARCODE BELOW 10) and another titeled `Sample`
    fasta_ref,
    gb_report,
    gene_name,
    threshold,
    filter_silent,
    min_reads = 100
){
  
  # CHECKS
  if (! length(list.files(path=folder, pattern=pattern) > 0)){ print(sprintf('The folder was not found or does not contnain any files matching the pattern %s', pattern)); return }
  
  phare = list()
  
  # ASSIGN STUFF
  phare$ref = read.fasta(fasta_ref, seqtype = 'DNA')[[1]] %>%  toupper()
  
  phare$samples = samples
  phare$gb_rep = rjson::fromJSON(file=gb_report)
  phare$gene = gene_name
  phare$threshold = threshold
  phare$filter_silent = filter_silent
  phare$min_reads = min_reads
  
  # INITIAL OPERATIONS
  phare$df = read_folder(folder, pattern)
  phare = phare %>% make_hapl_df()
  
  phare
}



### A FUNCTION TO READ ALL FILES IN A PATH
# it assumes three numbers the file name: first, the barcode/Source, then the subset and lastly the replicate number for in silico replicates
read_folder = function(path, pattern = '*.tsv') {
  #create a list of the files from your target directory
  file_list <- list.files(path=path, pattern=pattern)
  print(file_list)
  
  #initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable
  raw <- data.frame()
  
  
  for (i in 1:length(file_list)){
    source = str_extract_all(file_list[i], "(?<=[barcode|_|\\.])\\d+")[[1]]
    temp_data <- read.csv(paste(path, file_list[i], sep='/'), sep='\t', header=TRUE, check.names=FALSE, colClasses='character') %>% #each file will be read in, specify which columns you need read in to avoid any errors
      mutate(
        Source = ifelse(length(source > 0), source[1], 'unknown'),
        sub = ifelse(length(source > 1), source[2], 'unknown'),
        in_silico_rep = ifelse(length(source > 2), source[3], 'unknown'),
      )
    raw <- dplyr::bind_rows(raw, temp_data) #for each iteration, bind the new data to the building dataset
  }
  tibble(raw %>% dplyr::rename(`haplotype` = hapl)) #return raw
}




### A FUNCTION TO JOIN THE SAMPLES TO THE DATA FRAME AND CALCULATE THE TOTAL NUMBER OF READS
# hapl: the df as returned by read_folder
# samples: a dataframe linking a barcode/Source to a sample ID
# min_reads: The minimal number of reads required in each sample
make_hapl_df = function(phare){
  hapl = phare$df
  snp_list = colnames(hapl)[1:grep("hapl", colnames(hapl))-1]
  
  hapl = hapl %>%
    mutate(
      Source = as.factor(Source),
      count = as.numeric(count)
    ) 
  totals = hapl %>% group_by(Source, sub, in_silico_rep) %>% summarize(total_reads = sum(count))
  
  phare$df = hapl %>%
    left_join(
      totals,
      by=c('Source', 'sub', 'in_silico_rep')
    ) %>%
    group_by(Source) %>% 
    filter(
      total_reads > phare$min_reads
    ) %>%
    left_join(
      phare$samples,
      by=c(Source='Barcode')
    ) %>%
    mutate(
      Sample_name = factor(Sample, levels=levels(as.factor(phare$samples$Sample))),
      Sample = as.factor(paste(Sample, sub, in_silico_rep, sep='-'))
    )

  phare$df$Sample_name = droplevels(phare$df$Sample_name)
  # maybe different because it's grouped?
  # phare$df = phare$df %>% mutate(
  #   Sample_name = droplevels(Sample_name)
  # )
  
  phare
}












# CREATE SNP DATAFRAME AND A REFERENCE SEQUENCE
# gb_report: A reportfile downloaded from GenBank in .json format, which corresponds to the reference sequence. It is used to calculate exon and intron positions.
# hapl: The dataframe with all the SNPs => used to create the output df
make_snp_df = function(phare) {
  # conveniene
  gb_report = phare$gb_rep
  hapl = phare$df
  
  # extract exons from the genbank file
  ref_begin = as.numeric(gb_report[['transcripts']][[1]][['genomicRange']][['range']][[1]][['begin']])
  exons = gb_report[['transcripts']][[1]][['exons']][['range']] %>% 
    map_dfr(~ .x) %>% 
    mutate(
      begin = as.numeric(begin) - ref_begin + 1,
      end = as.numeric(end) - ref_begin + 1,
      length = end - begin
    )
  
  # filter snps for snps in exons (should already be done in python), calculate the total intron space and make an exon reference sequence
  oldsnp_list = as.numeric(colnames(hapl)[1:grep("haplotype", colnames(hapl))[[1]]-1])
  
  snps = tibble(`snp` = c(), `expos` = c(), `snp_ex` = c(), `aapos` = c(), `tripi` = c())
  ref_ex = c()
  for (row in 1:nrow(exons)){
    ref_ex = c(ref_ex, phare$ref[exons[[row, 'begin']]:exons[[row, 'end']]])
    intron_sum =  ifelse(row > 1, exons[row-1, 'intron_sum']+exons[row,'begin']-exons[row-1,'end']-1, 0)[[1]]
    exons[row,'intron_sum'] = intron_sum
    snps_in_exon = oldsnp_list[oldsnp_list > exons[[row,'begin']] & oldsnp_list < exons[[row,'end']]]
    
    snps = rbind(snps, tibble(
      `snp` = snps_in_exon,
      `exon` = row,
      `expos` = (snps_in_exon - intron_sum),
      `aapos` = ceiling((snps_in_exon - intron_sum) / 3),
      `tripi` = ((snps_in_exon - intron_sum)-1) %% 3 + 1
    ))
  }
  
  # SET NEW STUFF
  phare$ref_ex = ref_ex
  phare$ref_aa = translate(ref_ex)
  phare$snps = snps
  
  phare
}
# tmp = copy(dhps) %>% make_snp_df()




# FILTER DUPLICATE SNPs
# if we have a duplicate SNP (same AA position but different nucleotide)
filter_dupl_snps = function(phare){
  snps = phare$snps # SOMEHOW NOT PASSED AS REFERENCE??
  
  if (nrow(snps) > 1){
    snp_dupl = c()
    for (row in 2:nrow(snps)){
      if (snps[row-1,]$aapos == snps[row,]$aapos){
        snp_dupl = c(snp_dupl, row)
      }
    }
    if (length(snp_dupl) > 0) { phare$snps = snps[-snp_dupl,]; print(snp_dupl) }
  }
  phare
}
# tmp2 = copy(tmp) %>% filter_dupl_snps()


# TRANSLATE NUCLEOTIDE SNPS TO AMINO ACIDS
phare_translate = function(phare){
  hapl = data.table(phare$df)
  # TRANSLATE NUCLEOTIDE SNPS TO AMINO ACIDS
  for (snp_i in 1:nrow(phare$snps)){
    snp = phare$snps[snp_i,]
    ref_triplet = phare$ref_ex[(1:3) + snp$expos - snp$tripi]
    reftripstr = paste(ref_triplet, collapse='')

    print(sprintf("processing nucleotide position %s", snp$snp))
    
    hapl[get(as.character(snp$snp)) == reftripstr, paste(snp$aapos, 'aa', sep='_') := translate(ref_triplet)]
    hapl[is.na(get(paste(snp$aapos, 'aa', sep='_'))), paste(snp$aapos, 'aa', sep='_') := tolower(translate(strsplit(get(as.character(snp$snp)), split='')[[1]])), by=seq_len(length(Source))]
    
    hapl[ ,paste(snp$aapos, 'reftrip', sep = '_') := reftripstr]
  }
  phare$df = tibble(hapl)
  phare
}


# CALCULATE THE FREQUENCY OF THE AA VARIANTS
phare_frequency = function(phare){
  hapl = phare$df
  # concatenate amino acid snps to a haplotype
  hapl$aa_hapl = hapl[paste(phare$snps$aapos, 'aa', sep='_')] %>% apply(1, paste, collapse='')
  
  
  # group by aa haplotypes in each sample and calculate their frequency (as part of all reads)
  hapl = hapl %>% 
    group_by(Sample, aa_hapl) %>%
    filter(haplotype == first(haplotype)) %>%
    select(- count) %>%
    left_join(
      hapl %>% group_by(Sample, aa_hapl) %>% summarize(count = sum(count))
    ) %>%
    mutate(
      frequency = count / total_reads,
    )
  
  phare$df = hapl
  phare
}
# tmp4 = copy(tmp3) %>% phare_frequency()





# APPLY A TRESHOLD AND THEN RECALCULATE HAPLOTYPES AND THEIR RELATIVE FREQUENCIES FOR ALL THE HAPLOTYPES WHICH STILL EXIST (ARE NOT CONSIDERED NOISE)
apply_threshold = function(phare) {
  hapl = phare$df
  snps = phare$snps
  
  # FOR EACH SNP CREATE A LIST OF VARIANTS
  snp_silent = c()
  for (snp_i in 1:nrow(snps)){
    snp = snps[snp_i,]
    wt_aa = phare$ref_aa[[snp$aapos]]
    variants = unique(c(wt_aa, hapl[hapl$frequency > phare$threshold,][[paste(snp$aapos, 'aa', sep='_')]])) # wildtype plus all AA variants occurring in a snp
    snps[snp_i, 'variants'] = paste(variants, collapse='')
    snps[snp_i, 'wt'] = wt_aa
    if (length(unique(toupper(variants))) == 1){ if (unique(toupper(variants)) == wt_aa){
      snp_silent = c(snp_silent, snp_i)
    }}
  }
  # print(snps)
  # FILTER OUT SILENT MUTATIONS
  if (phare$filter_silent){
    if (length(snp_silent) > 0) { snps = snps[-snp_silent,] }
  }
  
  # print(snp_silent)
  
  # NOW THAT WE HAVE DEFINED ALL SNPs we can do the haplotypes again
  hapl$aa_hapl_plot = hapl[paste(snps$aapos, 'aa', sep='_')] %>% apply(1, paste, collapse='')
  
  hapl = hapl %>%
    group_by(Sample, aa_hapl_plot) %>%
    filter(haplotype == first(haplotype)) %>%
    select(- count) %>%
    left_join(
      hapl %>% group_by(Sample, aa_hapl_plot) %>% summarize(count = sum(count))
    ) %>%
    mutate(
      frequency = count / total_reads,
    ) %>%
    ungroup()
  
  


  # calculate haplotype frequency (as part of reads not considered noise)
  hapl = hapl %>%
    left_join(
      hapl  %>% filter(frequency > phare$threshold) %>% group_by(Sample) %>% summarize(total_over_threshold = sum(count)),
      # by = c('Sample', 'sub', 'in_silico_rep')
    ) %>%
    mutate(
      freq_plot = ifelse(frequency > phare$threshold, count/total_over_threshold, frequency - phare$threshold),
      aa_hapl_plot = ifelse(aa_hapl_plot %in% unique(filter(hapl, frequency > phare$threshold)$aa_hapl_plot), aa_hapl_plot, ''),
      aa_hapl_plot = factor(aa_hapl_plot, levels=sort(unique(aa_hapl_plot)))
    )

  phare$snps = snps
  phare$df = hapl
  phare
}
# tmp5 = copy(tmp4) %>% apply_threshold()










make_plot_df = function(phare){
  hapl = phare$df %>% arrange(Sample_name, freq_plot)
  snps = phare$snps
  
  
  # add plot-specific data frame columns
  for (row in 1:nrow(hapl)){
    hapl[row, 'xmin'] = ifelse(row == 1,
       0,
       ifelse (hapl[row-1, 'freq_plot'] > 0,
               hapl[row-1, 'xmin']+hapl[row-1, 'freq_plot'],
               hapl[row-1, 'xmin']  # entries with a freq_plot < 0 are not to be displayed
       )
    )
  }
  
  # CREATE A LIST OF MUTATIONS
  mutations = tibble(aapos=c(), variants=c(), mut=c(), wt=c(), var_i=c())
  for (row in 1:nrow(snps)){
    snp = snps[row,]
    for(i in 1:nchar(snp$variants)){
      if(substr(snp$variants,i,i) == snp$wt){
        mutations = rbind(mutations, tibble(aapos=snp$aapos, variants=snp$variants, mut='', wt=snp$wt, var_i=i))
      }else{
        mutations = rbind(mutations, tibble(aapos=snp$aapos, variants=snp$variants, mut=substr(snp$variants,i,i), wt='', var_i=i))
      }
    }
  }
  mutations = mutations %>% mutate(aapos = factor(aapos))
  
  
  melted = hapl %>%
    filter(freq_plot > 0) %>%
    arrange(Sample, -freq_plot) %>%
    select(c('aa_hapl_plot', 'Sample', 'Sample_name', 'freq_plot', 'xmin', all_of(paste(snps$aapos, 'aa', sep='_')))) %>%
    reshape2::melt(id.vars=c('aa_hapl_plot', 'Sample', 'Sample_name', 'freq_plot', 'xmin')) %>%
    mutate(aa_hapl_plot = factor(aa_hapl_plot, levels = sort(as.character(unique(hapl$aa_hapl_plot))))) %>%
    left_join(mutate(snps, variable = factor(paste(aapos, 'aa', sep='_'))), by='variable') %>% 
    rowwise %>% mutate(
      variant_i = regexpr(value, variants)[[1]]
    )
  
  phare$mutations = mutations
  phare$melted = melted
  phare$plotdf = hapl
  phare
}
# tmp6 = copy(tmp5) %>% make_plot_df()








































phare_plot = function(phare, supl=''){
  hapl = phare$df
  melted = phare$melted
  mutations = phare$mutations
  snps = phare$snps
  gene = phare$gene
  threshold = phare$threshold
  
  p1 = (
    ggplot(
      melted,
      aes(
        # y= as.numeric(Sample)*1.2,
      )
    )
    + theme_classic()
    + theme(
      # axis.text.y = element_text(angle = 90, vjust = 1, hjust=0),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      legend.text = element_text(family='mono'),
      legend.title = element_blank(),
      # legend.position = 'none'
    )
    + scale_colour_brewer(
      type='seq',
      palette='PuBuGn',
      aesthetics = c('colour', 'fill')
    )
    + scale_x_continuous(position='top', breaks=1:length(unique(melted$variable)), labels=as.character(snps$aapos), limits=c(NA, length(unique(melted$variable)) + .5))
    + scale_y_continuous(breaks=sort(unique(as.numeric(hapl$Sample_name)))*1.2-.55,
                         labels=(summarize(group_by(hapl[c('Sample_name', 'total_reads')], Sample_name), tot=max(total_reads)) %>% mutate(lab = paste(as.character(Sample_name), tot, sep='\nn=')))$lab,
                         position = 'left',)
    + coord_cartesian(ylim = c(0, max(as.numeric(hapl$Sample_name))*1.2+.5-.55), clip = 'off')
    + geom_rect(
      aes(
        ymin = xmin + as.numeric(Sample_name)*.2,
        ymax = xmin+freq_plot + as.numeric(Sample_name)*.2,
        xmin = as.numeric(variable) -.2 +
          .4/nchar(variants) * (variant_i-1),
        xmax = as.numeric(variable) -.2 +
          .4/nchar(variants) * variant_i,
        fill = aa_hapl_plot
        # fill=variant_i == 1
      ),
      color=NA,
      # fill=value
    )
    
    
    
    
    ## ANNOTATIONS
    
    # DIVIDER RECTS BETWEEN SAMPLES
    + geom_rect(
      aes(xmin = 1 - .2,
          xmax = length(unique(melted$variable)) + .5,
          ymin = yax,
          ymax = yax + .2,
      ),
      fill = 'grey',
      alpha = .15,
      data = tibble(yax = unique(as.numeric(melted$Sample_name))*1.2)
    )
    
    # V-LINES
    + geom_vline(xintercept=c(1:length(unique(melted$variable))) -.2 + (.4/nchar(subset(mutations, wt != '')$variants)))
    
    # HAPLOTYPE LABEL
    + geom_text(
      aes(
        y=max(as.numeric(hapl$Sample_name)) * 1.2 + .2,
        x=as.numeric(aapos) - .2 + .4/nchar(variants) * (var_i - .5),
        label = ifelse(mut == '', wt, tolower(mut))
      ),
      data = mutations,
      family='mono'
    )
  )
  p1
  p2= (
    ggplot(
      hapl,
      aes(
        # y=as.numeric(Sample)*1.2
      ),
    )
    + theme_classic()
    + theme(
      legend.title = element_text(),
      # legend.position = 'top',
      axis.title.y =  element_blank(),
      axis.text.x.top = element_text(hjust=-0.3, lineheight = 0.8, family='mono'),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor.x = element_line(linewidth=0.1, linetype=5)
    )
    + scale_colour_brewer(
      type='seq',
      palette='PuBuGn',
      aesthetics = c('colour', 'fill')
    )
    + guides(fill='none', color='none')
    # + scale_color_brewer(palette='Paired', drop=FALSE)
    # + scale_fill_brewer(palette='Paired', drop=FALSE)
    + geom_vline(xintercept=c(1:max(as.integer(hapl$aa_hapl_plot))))
    
    # + scale_y_continuous(breaks=1:length(unique(hapl$Sample))*1.2, labels=levels(hapl$Sample))
    + scale_x_continuous(position='top', expand=c(0.00,0), limits=c(1,max(as.integer(hapl$aa_hapl_plot))+1), breaks=c(1:max(as.numeric(hapl$aa_hapl_plot))), minor_breaks=seq(1.1,1.9, by=(1.9-1.1)/5), labels=paste(c('Noise', levels(hapl$aa_hapl_plot)[2:length(levels(hapl$aa_hapl_plot))]), '\n', sep=''), name='')
    + coord_cartesian(ylim=c(.55,max(as.numeric(hapl$Sample_name))*1.2+.5), clip = 'off')
    
    # HAPLOTYPE BARS
    + geom_rect(aes(
      ymax=as.integer(Sample_name)*1.2+.2,
      ymin=as.integer(Sample_name)*1.2-.2,
      xmax=as.integer(aa_hapl_plot) + 0.1 + freq_plot*0.6,
      xmin=as.integer(aa_hapl_plot) + 0.1,
      fill=aa_hapl_plot,
    ),
    data=subset(hapl, freq_plot > 0),
    color=NA)
    
    # NOISE POINTS
    + geom_point(aes(
      color = aa_hapl_plot,
      y = as.numeric(Sample_name)*1.2,
      x = frequency * 1/threshold * 0.8 + 0.1 + 1
    ), data=subset(hapl, freq_plot < 0 & aa_hapl_plot != ''))
    + geom_point(aes(
      y = as.numeric(Sample_name)*1.2,
      x = frequency * 1/threshold * 0.8 + 0.1 + 1
    ), color='grey', data=subset(hapl, freq_plot < 0 & aa_hapl_plot == ''))
    
    # PERCENT TEXT
    + geom_text(
      aes(
        x = as.integer(aa_hapl_plot) + .1 + freq_plot*0.6,
        y = as.integer(Sample_name)*1.2,
        label = paste(round(freq_plot*100, digits=0), '%', sep='')
      ),
      hjust=0,
      size=3.2,
      data=subset(hapl, frequency > threshold))
    
    # NOISE PERCENT
    + geom_text(
      aes(x=posx, y=max(as.integer(hapl$Sample_name))*1.2+.5+.8, label=paste(lab, '%', sep='')),
      data=data.frame(posx=seq(1.1,1.9, by=(1.9-1.1)/2), lab=as.character(seq(0, threshold*100, by=threshold*100/2))),
      size=2
    )
    
    # DIVIDER LINES BETWEEN SAMPLES
    + geom_rect(
      aes(xmin = 1,
          xmax = length(unique(hapl$aa_hapl_plot)) + 1,
          ymin = yax + .6,
          ymax = yax + .8,
      ),
      fill = 'grey',
      alpha = .15,
      data = tibble(yax = unique(as.numeric(hapl$Sample_name))*1.2)
    )
  )
  p2
  fig = ggarrange(p1, p2, labels=c('A', 'B'), ncol=2, align='h', common.legend=TRUE, legend='top', widths=c(2,3))
  fig = annotate_figure(fig,
                  top = text_grob(paste(gene, supl, sep=''))) + bgcolor('White')
  
  phare$fig = fig
  phare
}



save_plot = function(phare, width=320, height=180, dir='', sup='phare'){
  sup = ifelse(sup == '', '', paste('_', sup))
  ggsave(sprintf('%s/%s_%s%s.png', dir, phare$gene, phare$threshold, sup), phare$fig, width=width, height=height, units='mm')
  phare
}


## UTILITY
get_small_df = function(phare){
  snpdf = (phare %>% make_snp_df())$snps
  phare$df %>% select(-paste(snpdf$aapos, 'reftrip', sep='_'), -as.character(snpdf$snp), -sub, -in_silico_rep, -haplotype)
}
get_small_df_no_snp = function(phare){
  snpdf = (phare %>% make_snp_df())$snps
  phare$df %>% select(-paste(snpdf$aapos, 'reftrip', sep='_'), -paste(snpdf$aapos, 'aa', sep='_'), -as.character(snpdf$snp), -sub, -in_silico_rep, -haplotype)
}

get_noise = function(phare, expected){
  tmp = phare %>% get_small_df() %>%
    left_join(
      expected
    ) %>% rowwise() %>%
    mutate(
      rightpos=ifelse(freq_plot > 0, ifelse(aa_hapl_plot %in% aa_var_list, 1, 0),0),
      falsepos=ifelse(freq_plot > 0, ifelse(! aa_hapl_plot %in% aa_var_list, 1, 0),0),
      rightneg=ifelse(freq_plot < 0, ifelse(! aa_hapl_plot %in% aa_var_list, 1, 0),0),
      falseneg=ifelse(freq_plot < 0, ifelse(aa_hapl_plot %in% aa_var_list, 1, 0),0),
    ) %>%
    # print(n=100)
    group_by(Sample_name) %>% summarize(
      totalreads=sum(count),
      noise=sum(ifelse(rightneg==1, count, ifelse(falsepos==1, count,0))),
      correctreads=sum(ifelse(rightpos==1, count, ifelse(falseneg==1, count,0))),
      rightpos=sum(rightpos),
      falsepos=sum(falsepos),
      rightneg=sum(rightneg),
      falseneg=sum(falseneg),
      # mistakes=sum(mistakes),
    )
}