#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(optparse)

set.seed(42)

##################### Options ########################
opts=list(
  make_option("--dasperannotfile", default="", type="character",
              help="Dasper annotation file from dasper.R (default: none)"),
  make_option("--librarysize", default="", type="integer", help="STAR library size"),
  make_option("--outdir", default="dasperannots", type="character",
              help="Output directory for dasper annotation files (default: dasperannots)"),
  make_option("--min", default=50, type="integer", help="Min intron length (default: 50)"),
  make_option("--max", default=500000, type="integer", help="Max intron length (default: 500,000)"),
  make_option("--sf", default=1e7, type="integer", help="scale factor (default: 1e7)")
)
params=parse_args(OptionParser(option_list=opts))
message("Options: \n")
print(params)

basefile = basename(params$dasperannotfile)
libsize = params$librarysize
message("Library size for ", basefile, " is ", libsize)

dir.create(params$outdir, showWarnings = FALSE)
########################################################

# read Dasper annotation file
d = read_delim(params$dasperannotfile, delim="\t")

# scale Raw junction counts to the STAR library size and
# then filter junctions based on min and max intron size
d2 = d %>%
  mutate(scaled = junctioncounts/libsize * params$sf) %>% 
  filter(width >= params$min & width <= params$max) %>%
  select(seqnames, start, end, strand, type, junctioncounts, scaled, gene_id_start)
d2

write_delim(d2, file.path(params$outdir, paste0(basefile, ".filtered.tsv")), delim="\t")

sessionInfo()


