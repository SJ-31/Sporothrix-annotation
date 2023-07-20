library(glue)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
go_csv <- args[2]
stat_file <- args[1]
go_mapping <- read.csv(args[2]) %>%
  filter(GO != "None")
stat_header <- strsplit(readLines(stat_file)[2], split = "\t") %>% unlist()
match <- regexec(".*/(.*)_.*", stat_file)
sample_name <- unlist(regmatches(stat_file, match), use.names = FALSE)[2]

population <- read.table(stat_file, sep = "\t") %>%
  `colnames<-`(stat_header) %>%
  filter(GeneId %in% go_mapping$Gene_Name)
writeLines(population$GeneId, glue("{sample_name}_population.txt"))

sample <- population %>%
  filter(variants_impact_HIGH > 0)
writeLines(sample$GeneId, glue("{sample_name}_interest.txt"))

header <- "GoStat IDs Format Version 1.0"
formatted_go <- go_mapping$GO %>%
  lapply(., gsub,
    pattern = "_.*;",
    replacement = ","
  ) %>%
  lapply(., sub, pattern = "_.*", replacement = "") %>%
  unlist(use.names = FALSE)
formatted_mappings <- data.frame(
  ids = go_mapping$Gene_Name,
  GO = formatted_go
)

output <- "ontologizer_mappings.ids"
writeLines(header, output)
write_delim(formatted_mappings, output, append = TRUE,
  col_names = FALSE, delim = "\t"
)
