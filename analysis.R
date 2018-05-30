library(dplyr)
library(reshape2)
library(WriteXLS)
library(magrittr)

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

all_res <- pblapply(list.files("./data/", full.names = TRUE), function(file_name) {
  
  just_name <- strsplit(file_name, "/")[[1]] %>% last 
  all_lines <- readLines(file_name)
  
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
  
  #all_codons_per_gene <- all_codons_per_gene_raw[lengths(all_codons_per_gene_raw) %in% 1L:50]
  all_codons_per_gene <- all_codons_per_gene_raw
  
  if(length(all_codons_per_gene) > 0) {
    
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
    
    group_codons_counts <- lapply(names(all_codons_per_gene), function(ith_name) {
      i <- all_codons_per_gene[[ith_name]]
      
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
      
      group_counts_list <- lapply(1L:3, function(single_seq_part_id) 
        lapply(c("dangerous_codons", "with_alt", "double_danger", "no_alt"), 
               function(single_codon_type) {
                 data.frame(in_group = codon_types[[single_codon_type]], 
                            counts = group_codons_list[[single_seq_part_id]]) %>% 
                   group_by(in_group) %>% 
                   summarise(count = sum(counts)) %>%
                   mutate(type = single_codon_type,
                          freq = count/sum(count),
                          region = single_seq_part_id,
                          name = ith_name)
               }) %>% bind_rows()
      ) %>% bind_rows()

      browser()
    }) %>% bind_rows() %>% 
      mutate(file_name = just_name)
  } else {
    NULL
  }
}) %>% 
  bind_rows() %>% 
  select(file_name, name, region, type, in_group, count, freq)

WriteXLS(all_res, ExcelFileName = "./results/group_counts.xls")

all_res <- data.table::fread("./results/group_counts.csv", data.table = FALSE)

p <- group_by(all_res, file_name, region, type, in_group) %>% 
  summarise(freq = mean(freq)) %>% 
  filter(in_group) %>% 
  ggplot(aes(x = file_name, y = freq, fill = factor(region), 
             label = formatC(freq, 4))) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(position = position_dodge(width = 0.9), hjust = "right") +
  facet_wrap( ~ type, scales = "free_x", nrow = 1) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom") +
  

png("./results/plot.png", width = 2550, 
    height = 100 + length(unique(all_res[["file_name"]])) * 50)
p
dev.off()
#count_codon_df <- data.frame(gene = gene_name, nucleotide_count, all_codons_counts, part_codons_counts)

filter(all_res, 
       file_name == "GCF_000825665.1_Borrelia_crocidurae_str._03-02_cds_from_genomic.fna",
       type == "dangerous_codons",
       in_group,
       region != 1) %>% 
  select(region, freq, name) %>% 
  mutate(region = paste0("r", region)) %>% 
  dcast(name ~ region, value.var = "freq") %$%
  wilcox.test(r2, r3, paired = TRUE)

filter(all_res, 
       file_name == "GCF_000825665.1_Borrelia_crocidurae_str._03-02_cds_from_genomic.fna",
       type == "dangerous_codons",
       in_group,
       region != 2) %>% 
  select(region, freq, name) %>% 
  mutate(region = paste0("r", region)) %>% 
  ggplot(aes(x = freq, fill = region)) +
  geom_density(alpha = 0.3)
