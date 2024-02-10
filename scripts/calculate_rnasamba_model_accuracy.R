library(tidyverse)
library(caret)

files <- unlist(snakemake@input)

# read in and format model results for each dataset type
model_predictions <- files %>%
  set_names() %>%
  map_dfr(read_tsv, .id = "tmp") %>%
  mutate(tmp = gsub(".tsv", "", basename(tmp))) %>%
  separate(tmp, into = c("coding_type", "dataset_type"), sep = "_", remove = T) %>%
  mutate(coding_type = as.factor(coding_type),
         classification = as.factor(classification))  

# determine model performance
model_confusion_matrix <- confusionMatrix(data = model_predictions$classification, 
                                          reference = model_predictions$coding_type)

# convert the model performance metrics into a data frame
model_confusion_matrix_df <- bind_rows(data.frame(value = model_confusion_matrix$overall), 
                                       data.frame(value = model_confusion_matrix$byClass)) %>%
  rownames_to_column("metric") %>%
  mutate(dataset_type = model_predictions$dataset_type[1]) # set dataset type as column even tho it will be in the file name too

# convert the confusion matrix into a data frame
model_confusion_matrix_freq_df <- model_confusion_matrix %>% as.table() %>% data.frame()

# KC: this is 'confusionmatrix_{dataset_type}.tsv'
write_tsv(model_confusion_matrix_df, snakemake@output[['metrics']])

# KC: this is 'accuracy_metrics_{dataset_type}.tsv'
write_tsv(model_confusion_matrix_freq_df, snakemake@output[['freq']])
