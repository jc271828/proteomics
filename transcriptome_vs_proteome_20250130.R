
# load libraries ----------------------------------------------------------

library(tidyverse)
library(magrittr)
library(readxl)
library(edgeR)
library(arrangements)
library(reshape2)
library(writexl)
library(eulerr)
library(SuperExactTest)


# load files --------------------------------------------------------------

wbps19 <- read_tsv("wbps19_20240320_cleaned.txt", col_names = T)
ws286 <- read_csv("ws286_simple.csv", col_names = T)
rna <- read_csv("webster2018_all_detected_genes_edgeR_output_20231108.csv")
protein <- read_xlsx("10455_SupplementalData_050224.xlsx", sheet = "Table 5 Fed_v_Starv")

# FC in the core-generated spreadsheet is mean(fed normalized intensity)/mean(starved normalized intensity)
mean(as.numeric(protein[1, 6:9]))/mean(as.numeric(protein[1, 10:13]))
protein[1, 20]

# pval in the core-generated spreadsheet is two-sided unpaired variance-not-pooled t-test using log2(fed) and log2(starved) instead of fed and starved
t.test(log2(as.numeric(protein[1, 6:9])), log2(as.numeric(protein[1, 10:13])), alternative = "two.sided", paired = F, var.equal = F)$p.value
protein[1, 21]

# supposed to be FDR-corrected, not sure why it doesn't match. Re-adjust
p.adjust(protein$Fed_v_Starv_pval, method = "fdr")[1]
protein[1, 22]

# we want this
protein$log2FC_jingxian <- NA
for (i in 1:nrow(protein)) {
  protein$log2FC_jingxian[i] <- log2(mean(as.numeric(protein[i, 10:13]))/mean(as.numeric(protein[i, 6:9])))
}

protein$padj_jingxian <- p.adjust(protein$Fed_v_Starv_pval, method = "fdr")

core_pval_cutoff <- protein %>% mutate(row_number = 1:nrow(protein)) %>% filter(Fed_v_Starv_pval < 0.05 & abs(Fed_v_Starv_FC) >= 1.5)
core_stats <- protein %>% mutate(row_number = 1:nrow(protein)) %>% filter(Stats %in% c("Up", "Down"))

temp <- protein[setdiff(core_pval_cutoff$row_number, core_stats$row_number), ] # core's Up and Stats annotation is incomplete

core_padj_cutoff <- protein %>% mutate(row_number = 1:nrow(protein)) %>% filter(Fed_v_Starv_adjpval < 0.05 & abs(Fed_v_Starv_FC) >= 1.5)
jingxian_fdr_cutoff <- protein %>% mutate(row_number = 1:nrow(protein)) %>% filter(padj_jingxian < 0.05 & abs(Fed_v_Starv_FC) >= 1.5)


# match protein with WormBase ID ----------------------------------------------

# fed: fresh S-complete (<= 1 week old) with 1x fresh HB101 (<= 1 month old), 1 worm/ul, collected 18h post bleach, N2
# starved: virgin S-basal, 1 worm/ul, collected 24h post bleach, N2
conditions <- c("fed", "starved")

# clean proteomics dataframe
raw <- protein %>% select(c("Genes", "Accession", "113710 Fed Normalized", "113712 Fed Normalized", "113714 Fed Normalized", "113716 Fed Normalized", "113711 Starv Normalized", "113713 Starv Normalized", "113715 Starv Normalized", "113717 Starv Normalized"))
colnames(raw) <- c("identifier", "uniprot", "fed_rep1", "fed_rep2", "fed_rep3", "fed_rep4", "starved_rep1", "starved_rep2", "starved_rep3", "starved_rep4")
raw$row_number <- 1:nrow(raw)

raw_one_uniprot <- raw[str_detect(raw$uniprot, pattern = ";", negate = T), ]
raw_multiple_uniprot <- raw[str_detect(raw$uniprot, pattern = ";"), ]
df <- data.frame(matrix(nrow = 1, ncol = ncol(raw)))
colnames(df) <- colnames(raw)

for (i in 1:nrow(raw_multiple_uniprot)) {
  temp <- data.frame(identifier = raw_multiple_uniprot$identifier[i],
                     uniprot = unlist(str_split(raw_multiple_uniprot$uniprot[i], ";")),
                     fed_rep1 = raw_multiple_uniprot$fed_rep1[i],
                     fed_rep2 = raw_multiple_uniprot$fed_rep2[i],
                     fed_rep3 = raw_multiple_uniprot$fed_rep3[i],
                     fed_rep4 = raw_multiple_uniprot$fed_rep4[i],
                     starved_rep1 = raw_multiple_uniprot$starved_rep1[i],
                     starved_rep2 = raw_multiple_uniprot$starved_rep2[i],
                     starved_rep3 = raw_multiple_uniprot$starved_rep3[i],
                     starved_rep4 = raw_multiple_uniprot$starved_rep4[i],
                     row_number = raw_multiple_uniprot$row_number[i])
  df %<>% rbind(temp)
}
df <- df [-1, ]

raw <- rbind(raw_one_uniprot, df)

# match proteomics id with gene_id
raw$gene_id <- NA
raw$identifier %<>% str_remove("^CELE_")

for (i in 1:nrow(raw)) {
  match_by_uniprot <- wbps19 %>% filter(uniprot == raw$uniprot[i]) %>% .$gene_id %>% unique()
  match_by_identifier <- wbps19 %>% filter(identifier == raw$identifier[i]) %>% .$gene_id %>% unique()
  match_by_locus <- wbps19 %>% filter(locus == raw$identifier[i]) %>% .$gene_id %>% unique()
  
  if (length(match_by_uniprot) == 1) {
    raw$gene_id[i] <- match_by_uniprot
  } else if (length(match_by_identifier) == 1) {
    raw$gene_id[i] <- match_by_identifier
  } else if (length(match_by_locus) == 1) {
    raw$gene_id[i] <- match_by_locus
  } else if (length(match_by_uniprot) > 1) {
    raw$gene_id[i] <- wbps19 %>% filter(uniprot == raw$uniprot[i]) %>% .$deal_with_dup_uniprot %>% unique()
  }
}

# spot-check on duplicate gene_id
raw[which(str_detect(raw$gene_id, ",")), ]
raw[which(raw$identifier == "eft-4"), "gene_id"] <- "WBGene00001169" # WormBase

raw[which(str_detect(raw$gene_id, ",")), ]
raw[which(raw$uniprot == "G5EG90"), "gene_id"] <- "WBGene00016648" # UniProt

raw[which(str_detect(raw$gene_id, ",")), ]
raw[which(raw$uniprot == "G5EGJ4"), "gene_id"] <- "WBGene00016646" # UniProt

raw[which(str_detect(raw$gene_id, ",")), ]
raw[which(raw$uniprot == "G5ECN6"), "gene_id"] <- "WBGene00008359" # UniProt

raw[which(str_detect(raw$gene_id, ",")), ]
raw[which(raw$uniprot == "G5ECP8"), "gene_id"] <- "WBGene00008356" # UniProt

raw[which(str_detect(raw$gene_id, ",")), ]
raw[which(raw$uniprot == "G5EF55"), "gene_id"] <- "WBGene00009637" # UniProt

raw[which(str_detect(raw$gene_id, ",")), ]
raw[which(raw$identifier == "his-68;his-35"), "gene_id"] <- "WBGene00001942,WBGene00001909" # WormBase

# spot-check on empty gene_id
raw[which(is.na(raw$gene_id)), ] %>% print(n = 50)

# remove non-worm proteins
raw <- raw[str_detect(raw$uniprot, "_", negate = T), ]

raw[which(is.na(raw$gene_id)), ]
raw[which(raw$identifier == "retr-1"), "gene_id"] <- "WBGene00018416" # WormBase

raw[which(is.na(raw$gene_id)), ]
raw[which(raw$identifier == "F55A11.4"), "gene_id"] <- "WBGene00010077" # WormBase

raw[which(is.na(raw$gene_id)), ]
raw[which(raw$identifier == "tag-53"), "gene_id"] <- "WBGene00006432" # WormBase

raw[which(is.na(raw$gene_id)), ]
raw[which(raw$identifier == "Y82E9BR.16"), "gene_id"] <- "WBGene00022348" # WormBase

raw[which(is.na(raw$gene_id)), ]
raw[which(raw$identifier == "Y39E4B.10"), "gene_id"] <- "WBGene00012719" # WormBase

raw[which(is.na(raw$gene_id)), ]
raw[which(raw$uniprot == "C0HLB3"), "gene_id"] <- "WBGene00006911" # UniProt

# check uniqueness of gene_id
temp <- raw[which(duplicated(raw$gene_id)), ] %>% .$gene_id
raw %>% filter(gene_id %in% temp)

# rows with duplicated gene_id also have duplicated row_num, can delete duplicated entries
raw <- raw[!(duplicated(raw$gene_id)), ]

raw_merged <- merge(raw, ws286, by = "gene_id", sort = FALSE, all.x = T, all.y = F)
raw_merged[which(is.na(raw_merged$gene_name)), ]
raw_merged[which(raw_merged$gene_id == "WBGene00306133"), "gene_name"] <- "azyx-1" # WormBase
raw_merged[which(raw_merged$gene_id == "WBGene00018416"), "gene_name"] <- "cerv-1" # WormBase
raw_merged[which(raw_merged$gene_id == "WBGene00001942,WBGene00001909"), "gene_name"] <- "his-68,his-35" # WormBase

raw_merged[which(raw_merged$gene_id == "WBGene00306133"), "gene_biotype"] <- "protein_coding"
raw_merged[which(raw_merged$gene_id == "WBGene00018416"), "gene_biotype"] <- "protein_coding"
raw_merged[which(raw_merged$gene_id == "WBGene00001942,WBGene00001909"), "gene_biotype"] <- "protein_coding"

# all WormBase ID correspond to protein_coding genes, which is what they are supposed to do 
raw_merged$gene_biotype %>% unique()
raw_merged_wide <-
  raw_merged %>%
  mutate(id = paste(row_number, gene_id, sep = "_")) %>%
  select(c("id", "fed_rep1", "fed_rep2", "fed_rep3", "fed_rep4", "starved_rep1", "starved_rep2", "starved_rep3", "starved_rep4"))

# back to analyses
intensity_long <-
  raw_merged %>%
  mutate(id = paste(row_number, gene_id, sep = "_")) %>%
  select(c("id", "fed_rep1", "fed_rep2", "fed_rep3", "fed_rep4", "starved_rep1", "starved_rep2", "starved_rep3", "starved_rep4")) %>%
  pivot_longer(cols = -c("id"), values_to = "normalized_intensity", names_to = "sample") %>%
  as.data.frame()

# distribution of intensity
ggplot(data = intensity_long,
       aes(x = log10(normalized_intensity), color = sample))+
  geom_density()+
  labs(x = "log10(normalized_intensity)")+
  ggtitle("distribution of log10(normalized_intensity)")+
  theme_classic()+
  theme(aspect.ratio = 1)

# pearson correlation matrix
intensity_for_pearson <- raw_merged_wide
count_cormatrix <- round(cor(log2(intensity_for_pearson[, 2:9] + 1),
                             use = "all.obs",
                             method = "pearson"),
                         digits = 3)
melted_cormatrix <- melt(count_cormatrix)

# plot unordered Pearson correlation heatmap
ggplot(data = melted_cormatrix,
       aes(x = Var1, y = Var2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#900C3F",
                       high = "#6B33FF",
                       mid = "white",
                       midpoint = 0.96,
                       limit = c(0.91, 1))+
  theme_classic(base_size = 10)+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 3)+
  ggtitle("Pearson correlation matrix")+
  theme(aspect.ratio = 1)

# pca using covariance matrix
intensity_for_pca <- 
  intensity_for_pearson %>%
  select(-c("id"))

# Amy's way of transforming the data
intensity_for_pca$mean <- rowMeans(intensity_for_pca)
pca_using_cov <- prcomp(t(log2(intensity_for_pca[, 1:8]/intensity_for_pca$mean + 1)),
                        center = TRUE,
                        scale. = FALSE)
summary(pca_using_cov)

# getting the score matrix
pca.cov.df <- as.data.frame(pca_using_cov$x)
pca.cov.df <- cbind(pca.cov.df, rownames(pca.cov.df))
colnames(pca.cov.df)[9] <- "condition"
pca.cov.df$condition <- c(rep("fed", 4), rep("starved", 4)) %>% as.factor()
pca.cov.df$rep <- c(1:4, 1:4) %>% as.factor()

# PCA with confidence region
ggplot(pca.cov.df, mapping = aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3, mapping = aes(shape = rep)) +
  stat_ellipse(level = 0.95) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  labs(title = "rep-level PCA using cov matrix (95% confidence ellipses)",
       x = "PC1 (73.86 % of variance)",
       y = "PC2 (11.58% of variance)")


# # can we regenerate proteomics' stats like we do with RNASeq ----------------

# actually don't, because edgeR is designed for count-based data

# # put data in a DGEList object
# df <- DGEList(counts = raw_merged_wide[, c(2:9)],
#               genes = raw_merged_wide[, 1],
#               group = c(rep("fed", 4), rep("starved", 4)))
# df <- df[order(rowSums(df$counts), decreasing = TRUE), ]
# rownames(df$counts) <- df$genes$genes
# 
# # differential expression analysis
# df <- estimateCommonDisp(df)
# df <- estimateTagwiseDisp(df)
# # edgeR does log2(mean(intensity)) -- arithmetic means
# df$AveLogCPM %>% range()
# 
# # pairwise comparison using exact test
# # 1st group in the pair is baseline
# # so if the pair is c("A","B") then the comparison is B - A
# # so genes with positive logFC are up-regulated in B compared to A
# df$samples$group %<>% fct_relevel(c("fed", "starved"))
# levels(df$samples$group)
# 
# # DE analysis for starved vs fed
# DE.tag <- exactTest(df, df$tagwise.dispersion, pair = c("fed", "starved"))
# DE.tag$table$logFC %>% range() # edgeR does log2(mean(FC)) -- arithmetic means
# DE.tag$table$logCPM %>% range() # values in the logCPM column are exactly the same as AveLogCPM in df -- log2(mean(intensity))
# 
# # default sorting method of topTags is by PValue
# DE.tag_sort <- topTags(DE.tag, n = nrow(DE.tag$table))$table
# colnames(DE.tag_sort)[1:3] <- c("id" ,"log2FC", "log2_mean_norm_intens")
# 
# # add gene_name into result
# DE.tag_merge <- merge(select(mutate(raw_merged, id = paste(row_number, gene_id, sep = "_")), c("id", "gene_name")), DE.tag_sort, by = "id")
# DE.tag_merge$id %>% unique() %>% length()
# DE.tag_merge$gene_name %>% unique() %>% length()
# DE.tag_merge[which(str_detect(DE.tag_merge$gene_name, "his-")), ]
# 
# # visualize DEGs
# DE.tag_top_0.05 <- rownames(DE.tag_sort)[DE.tag_sort$FDR < 0.05] # 3799
# DE.tag_top_0.1 <- rownames(DE.tag_sort)[DE.tag_sort$FDR < 0.1] # 4330
# 
# how many DE proteins per the core's standards?
# if p-value < 0.05 and abs(FC) >= 1.5, then depending on the sign of FC, Stats is Up or Down (up means fed is greater than starved and vice versa).
# core_pval_cutoff %>% nrow() # 2627 but some entries correspond to multiple WormBase ID
# core_DE_protein <- raw_merged %>% mutate(id = paste(row_number, gene_id, sep = "_")) %>% filter(row_number %in% core_pval_cutoff$row_number) %>% .$id

# # overlap between EdgeR and core-determined DE proteins
# set.seed(1)
# plot(euler(list(EdgeR_FDR0.05_DEP = DE.tag_top_0.05,
#                 core_determined_DEP = core_DE_protein),
#            shape = "ellipse"),
#      quantities = TRUE) # do not do UP and DOWN DEG separately here, bc they will show overlap -- some genes are UP in certain tissues while DOWN in others
# length(intersect(DE.tag_top_0.05, core_DE_protein))/length(union(DE.tag_top_0.05, core_DE_protein)) # jaccard 0.632507
# 
# set.seed(1)
# plot(euler(list(EdgeR_FDR0.1_DEP = DE.tag_top_0.1,
#                 core_determined_DEP = core_DE_protein),
#            shape = "ellipse"),
#      quantities = TRUE) # do not do UP and DOWN DEG separately here, bc they will show overlap -- some genes are UP in certain tissues while DOWN in others
# length(intersect(DE.tag_top_0.1, core_DE_protein))/length(union(DE.tag_top_0.1, core_DE_protein)) # jaccard 0.5771598
# 
# # FDR 0.05 gives a better jaccard, proceed with it
# 
# # Smear plot
# plotSmear(DE.tag, de.tags = DE.tag_top_0.05, main = "starved vs fed at FDR 0.05 cutoff")
# 
# # MA plot
# ggplot()+
#   geom_point(data = filter(DE.tag_merge, FDR >= 0.05),
#              aes(x = log2_mean_norm_intens, y = log2FC),
#              alpha = 0.1)+
#   geom_point(data = filter(DE.tag_merge, FDR < 0.05),
#              aes(x = log2_mean_norm_intens, y = log2FC, color = "FDR < 0.05"),
#              alpha = 0.25)+
#   theme_classic()+
#   theme(aspect.ratio = 1)+
#   labs(x = "log2(normalized intensity)", y = "log2FC(starved/fed)")
# 
# # volcano plot
# ggplot()+
#   geom_point(data = filter(DE.tag_merge, FDR >= 0.05),
#              aes(x = log2FC, y = -log10(FDR)),
#              alpha = 0.2, size = 1.5)+
#   geom_point(data = filter(DE.tag_merge, FDR < 0.05),
#              aes(x = log2FC, y = -log10(FDR), color = "FDR < 0.05"),
#              alpha = 0.4, size = 1.75)+
#   theme_classic()+
#   labs(x = "log2(FC(starved/fed))", y = "-log10(FDR)")


# compare proteomics vs transcriptomics -----------------------------------

protein_for_comparison <- protein %>% mutate(row_number = 1:nrow(protein)) %>% select(c("row_number", "log2FC_jingxian", "padj_jingxian", "Peptides", "Detected_Imputed"))
colnames(protein_for_comparison)[c(2:5)] <- c("log2FC_protein", "FDR_protein", "unique_peptide_counts", "detected1_vs_imputed0_4fed_then_4starved_samples")
protein_for_comparison <- merge(raw[, c(1:2, 11:12)], protein_for_comparison, by = "row_number", all.x = T, all.y = F) # exclude non-worm protein

# deal with duplicated gene_id
temp <- protein_for_comparison[which(str_detect(protein_for_comparison$gene_id, ",")), ]
temp %<>% rbind(., .)
temp
temp$gene_id <- c("WBGene00001942", "WBGene00001909")

protein_for_comparison <- rbind(protein_for_comparison[which(str_detect(protein_for_comparison$gene_id, ",", negate = T)), ], temp)

# compare
DEP_FDR0.05 <- protein_for_comparison %>% filter(FDR_protein < 0.05 & abs(log2FC_protein) >= log2(1.5)) %>% .$gene_id %>% unique()
UP_DEP <- protein_for_comparison %>% filter(FDR_protein < 0.05 & log2FC_protein >= log2(1.5)) %>% .$gene_id %>% unique()
DOWN_DEP <- protein_for_comparison %>% filter(FDR_protein < 0.05 & log2FC_protein <= log2(1/1.5)) %>% .$gene_id %>% unique()

rna_for_comparison <- rna %>% select(c("gene_id", "sequence", "symbol", "logFC", "logCPM", "FDR"))
colnames(rna_for_comparison)[4:6] <- c("log2FC_transcript", "log2CPM", "FDR_transcript")
DEG_FDR0.05 <- rna_for_comparison %>% filter(FDR_transcript < 0.05) %>% .$gene_id %>% unique() %>% na.omit()
UP_DEG <- rna_for_comparison %>% filter(FDR_transcript < 0.05 & log2FC_transcript > 0) %>% .$gene_id %>% unique() %>% na.omit()
DOWN_DEG <- rna_for_comparison %>% filter(FDR_transcript < 0.05 & log2FC_transcript < 0) %>% .$gene_id %>% unique() %>% na.omit()

# overlap of detected transcripts vs proteins
bg_protein <- unique(protein_for_comparison$gene_id)
bg_gene <- unique(rna_for_comparison$gene_id) %>% na.omit()

set.seed(1)
plot(euler(list(detected_proteins = bg_protein,
                detected_transcripts = bg_gene),
           shape = "ellipse"),
     quantities = TRUE)
length(intersect(bg_protein, bg_gene))/length(union(bg_protein, bg_gene)) # jaccard 0.6752944

# hypergeometric test p-value < 2.2e-16
MSET(list(bg_protein,
          bg_gene),
     nrow(filter(ws286, gene_biotype == "protein_coding")),
     lower.tail = FALSE)

# overlap of DE transcripts vs proteins
set.seed(1)
plot(euler(list(DE_proteins = DEP_FDR0.05,
                DE_transcripts = DEG_FDR0.05),
           shape = "ellipse"),
     quantities = TRUE)
length(intersect(DEP_FDR0.05, DEG_FDR0.05))/length(union(DEP_FDR0.05, DEG_FDR0.05)) # jaccard 0.1688701

# hypergeometric test p-value = 1
MSET(list(DEP_FDR0.05,
          DEG_FDR0.05),
     length(union(bg_protein, bg_gene)),
     lower.tail = FALSE)

# overlap of DE transcripts vs proteins considering directions
set.seed(1)
plot(euler(list(UP_protein = UP_DEP,
                DOWN_protein = DOWN_DEP,
                UP_gene = UP_DEG,
                DOWN_gene = DOWN_DEG),
           shape = "ellipse"),
     quantities = TRUE)

# UP protein & UP gene hypergeometric test p-value = 2.75097e-30
MSET(list(UP_DEP,
          UP_DEG),
     length(union(bg_protein, bg_gene)),
     lower.tail = FALSE)

# UP protein & DOWN gene hypergeometric test p-value = 1
MSET(list(UP_DEP,
          DOWN_DEG),
     length(union(bg_protein, bg_gene)),
     lower.tail = FALSE)

# DOWN protein & UP gene hypergeometric test p-value = 1
MSET(list(DOWN_DEP,
          UP_DEG),
     length(union(bg_protein, bg_gene)),
     lower.tail = FALSE)

# DOWN protein & DOWN gene hypergeometric test p-value = 1.857956e-05
MSET(list(DOWN_DEP,
          DOWN_DEG),
     length(union(bg_protein, bg_gene)),
     lower.tail = FALSE)

# what about genes' whose transcript shows changes but protein level lags?
set.seed(1)
plot(euler(list(DE_gene_but_not_protein = setdiff(DEG_FDR0.05, DEP_FDR0.05),
                detected_protein = bg_protein),
           shape = "ellipse"),
     quantities = TRUE)

# what's the cause of a gene not changing in protein level but in transcript level? protein detected but not changing or not detected -- hypergeometric test p-value < 2.2e-16
MSET(list(setdiff(DEG_FDR0.05, DEP_FDR0.05),
          bg_protein),
     nrow(filter(ws286, gene_biotype == "protein_coding")),
     lower.tail = FALSE)

# if a gene's abundance changes at protein level but not transcript level (though detected at transcript level), does that imply regulatory mechanism?
set.seed(1)
plot(euler(list(DE_protein_but_not_gene = setdiff(DEP_FDR0.05, DEG_FDR0.05),
                detected_gene = bg_gene),
           shape = "ellipse"),
     quantities = TRUE)

# hypergeometric test p-value 0.9905243 --
# when a gene changes at protein but not transcript level,
# it's usually not because it's detected but transcript level doesn't change,
# but rather, it's minimal and not detected at transcript level 
MSET(list(setdiff(DEP_FDR0.05, DEG_FDR0.05),
          bg_gene),
     nrow(filter(ws286, gene_biotype == "protein_coding")),
     lower.tail = FALSE)

# generate combined result
comparison <- merge(protein_for_comparison, rna_for_comparison, by = "gene_id", all = T)
comparison %<>% select(c("row_number", "unique_peptide_counts", "detected1_vs_imputed0_4fed_then_4starved_samples", "gene_id", "identifier", "sequence", "symbol", "uniprot", "log2FC_protein", "FDR_protein", "log2FC_transcript", "FDR_transcript"))

for(i in 1:nrow(comparison)) {
  if(is.na(comparison$symbol[i])) {
    if (!(is.na(comparison$identifier[i]))) {
      comparison$symbol[i] <- comparison$identifier[i]
    }
  }
}

for(i in 1:nrow(comparison)) {
  if(is.na(comparison$symbol[i])) {
    if (!(is.na(comparison$sequence[i]))) {
      comparison$symbol[i] <- comparison$sequence[i]
    }
  }
}

comparison %<>% select(-c("identifier", "sequence"))

comparison %<>% arrange(gene_id)

# write_xlsx(comparison, "proteomics_vs_transcriptomics_20250130.xlsx")
