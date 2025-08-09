#!/usr/bin/env Rscript
library(optparse)
library(tidyverse)  # dplyr, tidyr, stringr, forcats, readr

set.seed(42)

##################### Options ########################
opts=list(
  make_option("--metadata", default=NULL, type="character",
              help="CSV file with columns: sample_id,filtered_file (required)"),
  make_option("--counts_col", default="scaled", type="character",
              help="Column name of the counts (default: %default)")
)
params=parse_args(OptionParser(option_list=opts))

if (is.null(params$metadata)) {
  stop("Error: --metadata is required. Provide a CSV with sample_id and filtered_file columns.")
}

countcolumn = params$counts_col

message("Options: \n")
print(params)

#########################################################


knownjunctions = c("annotated")
neojunctions = c("novel_acceptor", "novel_combo", "novel_donor", 
                 "novel_exon_skip", "unannotated", "ambig_gene")
dasperannots = c(knownjunctions, neojunctions)

results  = data.frame()

# Load metadata CSV (sample_id, filtered_file)
metadata <- read_csv(params$metadata)

for (i in seq_len(nrow(metadata))) {
  
  sample_id <- metadata$sample_id[i]
  file <- metadata$filtered_file[i]
  message("Calculating neojunction for sample: ", sample_id, " (file: ", file, ")")

  # Load dasper filtered file
  dasper = read.delim(file) %>% as_tibble()
  
  # Check if the count column exists
  if (!countcolumn %in% colnames(dasper)) {
    stop(paste0("Column '", countcolumn, "' not found in ", file))
  }
  
  # Add unique intron identifier
  # Create 'knownneo' column that is either 'knownjunction' or 'neojunction'
  dasper = dasper %>%
    mutate(intron = paste0(dasper$seqnames, ":", dasper$start, ":",
                           dasper$end, ":", dasper$strand, ":",
                           seq_along(dasper$seqnames)),
           type = factor(type),
           knownneo = fct_collapse(dasper$type,
                                   knownjunction = !!knownjunctions,
                                   neojunction = !!neojunctions))
  print(table(dasper$type, useNA="always"))
  print(table(dasper$type, dasper$knownneo, useNA="always"))
  
  
  # Neojunction calculation
  #         totalCts known known% neo neo% c5 c5% c3 c3% crypt crypt% novelannot novelannot% 
  # sample1 100      80    80%    20  20%  10 10% 5  5%  3     3%      2    2% (10+5+3+2 = neo)
  
  # Get the 'knownjunction' and 'neojunction' sums and percentage for each sample 
  knownneo = dasper %>% group_by(knownneo) %>%
    summarize(sum = sum(!!as.name(countcolumn))) %>%
    mutate(perc = sum/sum(sum)) %>%
    pivot_wider(names_from = knownneo, values_from = c(sum, perc))
  
  # Get the Type (Dasper annotations) sums and percentage for each sample 
  #   use .drop = FALSE to keep all types, even if they are not present in the data
  type = dasper %>% group_by(type, .drop = FALSE) %>%
    summarize(sum = sum(!!as.name(countcolumn))) %>%
    mutate(perc = sum/sum(sum)) %>%
    pivot_wider(names_from = type, values_from = c(sum, perc))
  
  sums = cbind(knownneo, type)
  
  # Data sanity checks
  #   Show that the Knonwjunction value is equivalent to the Annotated value 
  if (sums$sum_knownjunction != sums$sum_annotated) {
    message("Sum of Knownjunction does not equal Sum of Annotated")
  }
  #   Show that the Neojunction value is equivalent to the sum of all neojunctions
  sumofneo = sums %>%
    select(!contains("perc") & !contains("junction") & !contains("sum_annotated")) %>%
    rowSums()
  if (abs(sums$sum_neojunction - sumofneo) > 1e-5) {
    message("STOP")
  }
  
  sums$sample <- sample_id
  results = rbind(results, sums)
}

results = results %>% select(sample, everything()) %>%
  arrange(sample)

resultsfilename = paste0("neojunction_results.tsv")
write_tsv(results, resultsfilename)


sessionInfo()


