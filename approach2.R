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
