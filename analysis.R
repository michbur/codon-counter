library(dplyr)
library(reshape2)
library(WriteXLS)

all_codons <- expand.grid(c("A", "C", "G", "T"), c("A", "C", "G", "T"), c("A", "C", "G", "T")) %>% 
  apply(1, paste0, collapse = "")

codon_types <- rbind(data.frame(type = "normal", codon = all_codons),
                     data.frame(type = "dangerous_codons", codon = c("TGG", "TAC", "TAT", "TCA", "TTA", 
                                                                     "TGC", "TGT", "GAA", "GAG", "AAA", 
                                                                     "AAG", "CAA", "CAG", "TCG", "TTG", 
                                                                     "AGA", "CGA", "GGA")),
                     data.frame(type = "with_alt", codon = c("TCA", "TCG", "TTA", "TTG", "AGA", "CGA", 
                                                             "GGA")),
                     data.frame(type = "double_danger", codon = c("TGG", "TAC", "TAT", "TCA", "TTA")),
                     data.frame(type = "no_alt", codon = c("TGG", "TGT", "TGC", "TAC", "TAT", "CAA", 
                                                           "CAG", "GAA", "GAG", "AAA", "AAG"))) %>% 
  dcast(codon ~ type, value.var = "codon") %>% 
  group_by(codon) %>% 
  mutate_all(function(i) !is.na(i)) %>% 
  select(-normal) %>% 
  ungroup %>% 
  mutate(codon = as.character(codon))


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

all_codons_per_gene_raw <- lapply(genes_list, function(i) {
  # length(i) %/% 3 * 3 - only full codons
  factor(apply(matrix(i[1L:(length(i) %/% 3 * 3)], ncol = 3, byrow = TRUE), 1, paste0, collapse = ""),
         levels = all_codons) 
})

all_codons_per_gene <- all_codons_per_gene_raw[lengths(all_codons_per_gene_raw) %in% 100L:500]

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

# approach 1 ----------------------------

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


#count_codon_df <- data.frame(gene = gene_name, nucleotide_count, all_codons_counts, part_codons_counts)

# approach 2 ---------------------------------------


all_lengths <- lengths(all_codons_per_gene)

n_seq <- sapply(1L:max(all_lengths), function(i)
  sum(i <= all_lengths)) %>% 
  data.frame(pos = 1L:max(all_lengths), n_seq = .)

# codon type plot
lapply(colnames(codon_types)[-1], function(ith_codon_type) 
  lapply(all_codons_per_gene, function(i)
    data.frame(pos = 1L:length(i), codon = i, stringsAsFactors = FALSE)) %>% 
    bind_rows %>% 
    inner_join(codon_types, by = c("codon" = "codon")) %>% 
    group_by_(ith_codon_type, "pos") %>% 
    summarise(n_codon = length(codon)) %>% 
    ungroup %>% 
    inner_join(n_seq, by = c("pos" = "pos")) %>% 
    mutate(norm_n_codon = n_codon/n_seq)) %>% 
  bind_rows %>% 
  melt(id.vars = c("pos", "n_codon", "n_seq", "norm_n_codon")) %>% 
  na.omit %>% 
  ggplot(aes(x = pos, y = norm_n_codon, color = variable)) +
  #geom_point() +
  geom_smooth(method = "loess") +
  facet_wrap(~ variable + value, scales = "free_y")


library(ggplot2)
#filter(codon_pos_df, codon == "CAA") %>% 
ggplot(codon_pos_df, aes(x = pos, y = norm_n_codon, color = dangerous_codons)) +
  geom_smooth(method = "loess") 





# saving results ----------------------------------

WriteXLS(count_codon_df, ExcelFileName = "count_codon_table.xlsx")


