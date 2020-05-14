library(Seurat) #install.packages('Seurat')
library(dplyr) #install.packages('dplyr')
library(tidyr) #install.packages('tidyr')

adult_circulating <- readRDS("temp/dropseq.mtcut_50.norm.Rds")

adult_circulating_cluster <- data.frame(cluster=Seurat::Idents(adult_circulating)) %>%
  dplyr::add_rownames("UMI")

adult_circulating_exp <- data.frame(GetAssayData(object = adult_circulating, slot = 'data')) %>%
  dplyr::add_rownames(var="gene") %>%
  tidyr::gather(key="UMI", value="expression", -gene) %>%
  dplyr::select(UMI,gene, expression) %>%
  dplyr::mutate(UMI = gsub("\\.","-",UMI)) %>%
  dplyr::left_join(., adult_circulating_cluster, by='UMI') 

adult_markers <- readRDS("temp/dropseq.mtcut_50.norm.markers.Rds") %>%
  dplyr::mutate(log10_padj=log10(p_val_adj)) %>%
  dplyr::select(cluster, gene, Fraction_target = pct.1, Fraction_others = pct.2, Mean_log_FC = avg_logFC, 
                P_value=p_val, P_adjust = p_val_adj, log10_P_adjust=log10_padj)
