library(tidyverse)
inputargs <- commandArgs(TRUE)
annot_fn <- inputargs[1]
window_size <- as.numeric(inputargs[2])+1
save_fn_prefix <- inputargs[3]

#annot_fn <- "/work-zfs/abattle4/prashanthi/recount_networks/files_for_princy/degree.GeneSet"
gene_info <- readRDS("data/gene_metadata.Rds")

# Read annotation file. The first column corresponds to gene ensembl ID and second column is the network annotation

dat <- read_csv(annot_fn, col_names = F)

colnames(dat) <- c("ensembl_gene_id", "network_annotation")

gene_of_interest <- gene_info %>% inner_join(dat, by = c("ensembl_gene_id" = "ensembl_gene_id")) %>% 
				filter(chromosome_name!="X") %>% arrange(chromosome_name, start_position) %>%
				mutate(START = pmax(0,start_position-window_size), END = end_position+window_size)


class(gene_of_interest$chromosome_name) <- "numeric"

for(eachchr in unique(gene_of_interest$chromosome_name)){
	cat("starting chr", eachchr,"\n")
	bimfn <- paste0("data/1000G_EUR_Phase3_plink/1000G.EUR.QC.",
		eachchr, ".bim")
	bimfile <- read_tsv(bimfn, col_names = F)
	thischr <- gene_of_interest %>% filter(chromosome_name == eachchr)
		tmp <- left_join(bimfile, thischr, 
			by = c("X1"="chromosome_name")) %>% select(X1,X2,X3,X4,X5,network_annotation, START, END) %>%
		filter(X4 >= START, X4 <= END) %>% 
		group_by(X2) %>% mutate(ANNO = max(network_annotation)) %>% 
		select(-c("network_annotation", "START", "END")) %>% ungroup() %>% distinct()
		out_thischr <- left_join(bimfile, tmp, 
			by = c("X2" = "X2", "X1" = "X1", "X4" = "X4")) %>% 
		select("X1", "X2", "X4", "ANNO") %>% 
		rename(CHR = X1, SNP =X2, BP = X4, AN1 = ANNO) %>%
		arrange(BP) %>% mutate(AN1 = ifelse(is.na(AN1), 0, AN1)) %>% select(AN1)
		write_tsv(out_thischr, path = paste0(data/save_fn_prefix, eachchr, ".annot.gz"))
		cat("Done chr", eachchr,"\n")
}
