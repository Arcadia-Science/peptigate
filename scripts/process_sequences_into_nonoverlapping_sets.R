library(readr)
library(dplyr)
source("scripts/parse_sequence_information.R")

# read in data sets ------------------------------------------------------

# KC: the 'sequence' column is actually a sequence IDs (or 'contig names')
fai_col_names <- c("sequence", "length", "offset", "linebases", "linewidth")
clusters <- read_tsv(snakemake@input[['clusters']], col_names = c("rep", "cluster_member"))

# KC: these are the coding and noncoding datasets from the RNAChallenge
coding_validation <- read_tsv(unlist(snakemake@input[['validation_fai']])[1], col_names = fai_col_names)
noncoding_validation <- read_tsv(unlist(snakemake@input[['validation_fai']])[2], col_names = fai_col_names)
all_fai <- read_tsv(snakemake@input[['all_fai']], col_names = fai_col_names) %>%
  parse_sequence_information_from_seqkit_fai()

metadata <- read_tsv(snakemake@input[['metadata']]) %>%
  select(organism, genome = genome_abbreviation, set_name)

# homology reduction ------------------------------------------------------

# filter to homology reduced sequences -- keep sequences that had =<80% homology to other sequences
all_fai_filtered <- all_fai %>%
  mutate(sequence_id = gsub(" .*", "", sequence)) %>%
  filter(sequence_id %in% clusters$rep) %>%
  select(-sequence_id)

# data set separation -----------------------------------------------------

traintest <- all_fai_filtered %>%
  filter(!sequence %in% c(coding_validation$sequence, noncoding_validation$sequence)) %>% # remove sequences that are in the validation sets
  parse_sequence_information_from_seqkit_fai(dataset_type = FALSE) %>%
  left_join(metadata, by = c("genome"))

train <- traintest %>%
  filter(set_name == "train")
coding_train <- train %>% filter(sequence_type == "cdna")
noncoding_train <- train %>% filter(sequence_type == "ncrna")

test <- traintest %>%
  filter(set_name == "test")
coding_test <- test %>% filter(sequence_type == "cdna")
noncoding_test <- test %>% filter(sequence_type == "ncrna")

coding_validation_filtered <- all_fai_filtered %>%
  filter(sequence %in% coding_validation$sequence) %>% # filter to coding sequences
  filter(!sequence %in% traintest$sequence_id)         # remove sequences that overlap with the training/testing data

noncoding_validation_filtered <- all_fai_filtered %>%
  filter(sequence %in% noncoding_validation$sequence) %>% # filter to noncoding sequences
  filter(!sequence %in% traintest$sequence_id)            # remove sequences that overlap with the training/testing data

# augment the noncoding traintest set to balance sizes  ------------------

# To avoid biases in classification, it's better to have the same number of sequences in the different classes.
# This could also be done with changing the weights on the activation function of the model building in tensorflow,
# but that would require altering the source code of the RNAsamba tool.
# This approach should produce equivalent results.

num_rows_target <- nrow(coding_train)

noncoding_train <- noncoding_train %>%
  slice_sample(n = num_rows_target, replace = TRUE)

# write out contigs names for each set ------------------------------------

snakemake_output <- unlist(snakemake@output)
write.table(coding_train$sequence, snakemake_output[1], quote = F, col.names = F, row.names = F)
write.table(coding_test$sequence, snakemake_output[2], quote = F, col.names = F, row.names = F)
write.table(coding_validation_filtered$sequence, snakemake_output[3], quote = F, col.names = F, row.names = F)
write.table(noncoding_train$sequence, snakemake_output[4], quote = F, col.names = F, row.names = F)
write.table(noncoding_test$sequence, snakemake_output[5], quote = F, col.names = F, row.names = F)
write.table(noncoding_validation_filtered$sequence, snakemake_output[6], quote = F, col.names = F, row.names = F)
