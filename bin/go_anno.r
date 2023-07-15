# Add annotations from a mapping file to a SnpEff-annotated vcf file
# The annotations will be mapped by having the same gene names

args <- commandArgs(trailingOnly = TRUE)
read <- args[1]
reference <- read.csv(args[2])
output <- paste0(gsub(".*/|.vcf", "", read), "_GO.vcf")


return_match <- function(gene_id, reference) {
  match <- reference[reference$Gene_Name == gene_id, ]
  go <- match$GO
  if (go == "None") {
    go <- "GO:None"
  }
  uniprot <- paste0("UniProt:", match$UniProtKB_ID)
  return(
    paste0(uniprot, "|", go)
  )
}

name_from_anno <- function(annotation) {
  name <- strsplit(annotation, split = "|", fixed = TRUE) |> unlist()
  return(name[4])
}

add_anno <- function(row_num, vcf, ref) {
  current <- vcf[row_num, ]
  current$V8 <- paste0(
    current$V8, "|",
    return_match(name_from_anno(current$V8), ref)
  )
  entries <- current |> as.character() |> paste(collapse = "\t")
  return(entries)
}

all <- readLines(read)
header <- all[grep("#", all)]
vars <- read.table(read, sep = "\t")
annotated <- lapply(1:length(vars$V2), add_anno,
  vcf = vars,
  ref = reference
) |> unlist(use.names = FALSE)
writeLines(header, output)
cat(annotated, file = output, sep = "\n", append = TRUE)
