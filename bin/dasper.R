#!/usr/bin/env Rscript
library(dasper)
library(SummarizedExperiment)
library(readr)
library(optparse)

# This script takes a STAR SJ.out.tab file and adds annotation to each
# junction. The output is as follows for the 20 columns:
# col 1:  chromosome (equivalent to col 1 from SJ.out.tab)
# col 2:  first base of the intron (1-based) (equivalent to col 2 from SJ.out.tab)
# col 3:  last base of the intron (1-based) (equivalent to col 3 from SJ.out.tab)
# cols 4 to 18: Dasper various annotations for Gene, Transcript and Exon
# col 19: Dasper junction annotation ('type')
# col 20:  number of uniquely mapping reads crossing the junction (equivalent to col 7 from SJ.out.tab)
 

set.seed(42)

##################### Options ########################
opts=list(
  make_option("--prefix", default="", type="character", help="Prefix used for output file (default: none)"),
  make_option("--sj", default="", type="character", help="STAR SJ.out.tab file (default: none)"),
  make_option("--gtf", default="", type="character", help="GTF file (default: none)"),
  make_option("--outdir", default="dasperannots", type="character", help="Output directory for dasper annotation files (default: dasperannots)")
)
params=parse_args(OptionParser(option_list=opts))
message("Options: \n")
print(params)


# Validate required parameters
if (params$prefix == "") {
  stop("Please provide output prefix using --prefix")
}

if (params$sj == "") {
  stop("Please provide STAR SJ.out.tab file using --sj")
}
if (!file.exists(params$sj)) {
  stop(paste0("STAR SJ.out.tab file '", params$sj, "' does not exist"))
}

if (params$gtf == "") {
  stop("Please provide GTF file using --gtf")
}
if (!file.exists(params$gtf)) {
  stop(paste0("GTF file '", params$gtf, "' does not exist"))
}

dir.create(params$outdir, showWarnings = FALSE)
#########################################################

# Read STAR SJ file to take a look at it
sj = read.delim(params$sj, header = F, sep="\t")
dim(sj)
str(sj)
head(sj, n=10)

# Load SJ junction file
d = junction_load(params$sj)
d
d@rowRanges
rowRanges(d)
d@assays
dassay = assays(d)[['raw']]
assays(d)[['raw']][1:10,]  # dasper is only keeping Column 7 from SJ.out.tab (num uniquely mapped reads to junction)


# Annotate SJ junctions
d2 = junction_annot(d, ref = params$gtf)
d2
rowRanges(d2)
d2assay = assays(d2)[['raw']]
identical(dassay, d2assay)


# Convert rowRanges to a data.frame
df = as.data.frame(d2@rowRanges)
table(df$type)
lapply(df, typeof)
#str(df)
head(df)
dim(df)
#write.table(df, "df.tsv", sep="\t")  # doesn't work since lists are expanded
head(df$gene_id_start)


# Collapse list columns of the data.frame to character
# 
# https://www.biostars.org/p/271154/
# https://stackoverflow.com/questions/22377713/concatenating-the-list-elements-in-r
df2 = as.data.frame(lapply(df, function (x) {
  if (is.list(x)) {
    sapply(x, paste, collapse = ",")    
  }
  else { return(x) }
}))
lapply(df2, typeof)
#df2 %>% filter(gene_id_start == "ENSG00000279457.2")


# Output annotated junctions 
df2$junctioncounts = dassay  # add number of mapped reads for each junction
write_delim(df2, file.path(params$outdir, paste0(params$prefix, ".dasper.tsv")), delim="\t")

sessionInfo()



