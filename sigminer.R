suppressPackageStartupMessages({
library(sigminer)
library(dplyr)
library(readr)
library(reshape2)
library(kableExtra)
library(ggplot2)
library(gridExtra)
library(NMF)
library(BSgenome.Hsapiens.UCSC.hg38)
})

args <- commandArgs(trailingOnly = TRUE)


nmfruns <- as.integer(args[1])
nsigs <- as.integer(args[2])
threshold <- as.double(args[3])
method <- args[4]

ref_genome_lib <- "BSgenome.Hsapiens.UCSC.hg38"

# Get the current working directory
wd <- getwd()
# Specify the new directory path relative to the current working directory
data_dir <- file.path(wd, "RawData")
# Set working directory to data folder
setwd(data_dir)
# Get list of VCF files in current working directory
vcf_files <- list.files(pattern = "\\.vcf$", full.names = TRUE)
# Read all VCF files into a list
vcf_data <- lapply(vcf_files, read_delim, delim = "\t", comment = "##", col_names = TRUE, show_col_types = FALSE)
# Form MAF object needed for Sig Tally function
maf <- read_vcf(vcf_files, genome_build = "hg38", keep_only_pass = FALSE)

mt_tally <- sig_tally(maf, mode= "SBS",ref_genome = ref_genome_lib, genome_build = "hg38",useSyn = TRUE)

# vytiahnut idealny pocet signatur !!! a potom poslat do sig extractu
# mt est ma vyplut idealny pocet signatur co budeme posuvat do sig extract nsigs
mt_est <- sig_estimate(mt_tally$nmf_matrix,
  range = 2:10,
  nrun = 20, # increase this value if you wana a more stable estimation
  use_random = FALSE, # if TRUE, add results from randomized input
  seed = 123456,
  cores = 4,
  verbose = FALSE,
  save_plots = TRUE
)
# De novo
mt_sig <- sig_extract(mt_tally$nmf_matrix,
  n_sig = nsigs,
  nrun = nmfruns,
  seed = 123456,
  method = method,
  cores = 4)

mat_obj <- mt_tally$nmf_matrix %>% t()
# Fit
fit_res <- sig_fit(
  mat_obj,
  mt_sig,
  sig_db = "SBS_hg38",
  sig_index = "ALL",
  db_type =  "human-exome",
  show_index = TRUE,
  return_class = "data.table",
  method = "NNLS",
  auto_reduce = FALSE,
  return_error = TRUE,
  rel_threshold = threshold,
  mode = "SBS")

# Specify the new directory path relative to the current working directory
output_dir <- file.path(wd, "SigMinerOutput")
setwd(output_dir)

# Loop through each matrix in the list and save them individually
for (i in 1:4) {
  # Get the name of the matrix
  matrix_name <- names(mt_sig)[i]
  # Write the matrix to text file
  write.table(mt_sig[[i]], file = paste0(matrix_name, ".txt"), sep = "\t", quote = FALSE)
}

output_tally(mt_tally, output_dir, mut_type = "SBS")

sig_output <- file.path(output_dir, "Signatures")

setwd(sig_output)
output_sig(mt_sig, sig_output, mut_type = "SBS")

fit_output <- file.path(output_dir, "Refitting")

setwd(fit_output)
output_fit(fit_res, fit_output, mut_type = "SBS", sig_db = "SBS_hg38")

quit(status = 0)
