library(phangorn)
library(stringr)
library(glue)
library(TreeDist)
args <- commandArgs(trailingOnly = TRUE)

ref <- read.tree(args[1])
query <- NJ(dist.ml(read.dna(args[2], format = "fasta")))
name <- args[2] %>%
  str_replace_all("\\./", "") %>%
  str_replace("_.*", "")
query$tip.label <- query$tip.label %>% str_replace("\\|.*", "")
score <- TreeDistance(query, ref)
write.tree(query, glue("{name}_NJ.nwk"))
cat(glue("{name}\t{score}"), file = glue("{name}_score.txt"))
# TreeDistance is based on Generalized RobinsonFoulds distance
