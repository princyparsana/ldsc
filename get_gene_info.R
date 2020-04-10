library(recount)
library(biomaRt)

human <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl");

gene_info = getBM(
    attributes= c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand"), 
    values = TRUE, mart=human)

gene_info$start_position_bystrand <- ifelse(gene_info$strand==1, 
                            gene_info$start_position, gene_info$end_position)

saveRDS(gene_info, file = "data/gene_metadata.Rds")
