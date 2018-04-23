library(dplyr)

all_lines <- readLines("./data/sequences_e.coli.txt")

gene_id <- cumsum(grepl("^>", all_lines))

all_genes <- split(all_lines, gene_id)

gene_name <- unlist(lapply(all_genes, first))

genes_list <- lapply(all_genes, function(i) unlist(strsplit(i[-1], "")))

nucleotide_count <- lapply(genes_list, function(i) 
  factor(i, levels = c("A", "C", "G", "T")) %>% 
    table %>% 
    data.frame %>% 
    pull(Freq)) %>% 
  unlist(., use.names = FALSE) %>% 
  matrix(., ncol = 4, byrow = TRUE) %>% 
  data.frame()

colnames(nucleotide_count) <- c("A", "C", "G", "T")

all_codons <- expand.grid(c("A", "C", "G", "T"), c("A", "C", "G", "T"), c("A", "C", "G", "T")) %>% 
  apply(1, paste0, collapse = "")

all_codons_per_gene <- lapply(genes_list, function(i) 
  factor(apply(matrix(i, ncol = 3, byrow = TRUE), 1, paste0, collapse = ""),
         levels = all_codons) 
)


all_codons_counts <- lapply(all_codons_per_gene, function(i)
  table(i) %>% 
    data.frame %>% 
    pull(Freq)
) %>% 
  unlist(., use.names = FALSE) %>% 
  matrix(., ncol = 64, byrow = TRUE) %>% 
  data.frame()

colnames(all_codons_counts) <- all_codons


part_codons_counts <- lapply(all_codons_per_gene, function(i) {
  no_start_end <- i[-c(1, length(i))]
  
  
  
  which((floor(1L:(length(no_start_end))/(length(no_start_end) %/% 3 + 1)) ==
  ceiling(1L:(length(no_start_end))/(length(no_start_end) %/% 3 + 1)) - 1))
  
  
  browser()
})



#nucleotide_count <- data.frame(gene = gene_name, nucleotide_count)

