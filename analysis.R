library(dplyr)
library(WriteXLS)

all_lines <- readLines("./data/sequences_e.coli.txt")

gene_id <- cumsum(grepl("^>", all_lines))

all_genes <- split(all_lines, gene_id)

gene_name <- unlist(lapply(all_genes, first))

genes_list <- lapply(all_genes, function(i) unlist(strsplit(i[-1], "")))

# do not count codon with atypical nucleotides

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

all_codons_per_gene <- lapply(genes_list, function(i) {
  # length(i) %/% 3 * 3 - only full codons
  factor(apply(matrix(i[1L:(length(i) %/% 3 * 3)], ncol = 3, byrow = TRUE), 1, paste0, collapse = ""),
         levels = all_codons) 
})

# codon position composition
# strsplit(as.character(i[1L:start1]), "") %>% 
#   do.call(rbind, .) %>% 
#   apply(2, function(ith_col) data.frame(table(factor(ith_col, c("A", "C", "G", "T")))))

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
  
  len <- length(i)
  
  group_mod <- len %% 3
  
  # fuzzy splits - border codons in the case of unequal group belong
  # to two groups
  
  group_size <- len %/% 3 + group_mod
  
  start3 <- (len - group_size + group_mod)
  
  group_codons_list <- list(i[1L:group_size],
                              i[(group_size + 1 - group_mod):start3],
                              i[(start3 + 1 - group_mod):len]) %>% 
    lapply(function(ith_group) 
      table(ith_group) %>% 
        data.frame %>% 
        pull(Freq)) 
  
  group_codons_counts <- unlist(group_codons_list)
  names(group_codons_counts) <- unlist(lapply(paste0("r", 1L:3, "_"), function(ith_region_name) 
    paste0(ith_region_name, all_codons)))
  
  group_codons_fraction <- unlist(lapply(group_codons_list, function(ith_region) 
    ith_region/sum(ith_region)))
  names(group_codons_fraction) <- unlist(lapply(paste0("f", 1L:3, "_"), function(ith_region_name) 
    paste0(ith_region_name, all_codons)))
  
  c(group_codons_counts, group_codons_fraction)
}) %>% do.call(rbind, .)


count_codon_df <- data.frame(gene = gene_name, nucleotide_count, all_codons_counts, part_codons_counts)

WriteXLS(count_codon_df, ExcelFileName = "count_codon_table.xlsx")
