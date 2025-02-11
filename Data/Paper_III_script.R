rm(list=ls())


# Functions -------------------------------------------------------------------------
eDNAindex<- function(df, Sample_column, OTU_column, Counts_column, Biological.replicate ){ # 
	require(dplyr)
	require(rlang)
	Sample_column <- rlang::enquo(Sample_column)
	OTU_column    <- rlang::enquo(OTU_column)
	Counts_column <- rlang::enquo(Counts_column)
	Biological.replicate <- rlang::enquo(Biological.replicate)
	
	if (quo_name(Biological.replicate) %in% colnames(df)){
		print ("Averaging ratios between Biological replicates")
		
		df %>%
			group_by(!!Sample_column, !!OTU_column, !!Biological.replicate) %>%
			summarise ( sumreads = sum(!!Counts_column)) %>%  # This sums technical replicates
			group_by(!!Sample_column,!!Biological.replicate) %>% 
			mutate (Tot = sum(sumreads),
							Row.prop = sumreads / Tot)  %>%                      # This creates the proporttion on each biological replicate    
			group_by(!!Sample_column) %>% 
			mutate (nreps = n_distinct(!!Biological.replicate)) %>% 
			group_by(!!Sample_column, !!OTU_column) %>% 
			mutate(Row.prop=if_else(is.nan(Row.prop),0,Row.prop)) %>% 
			summarise (mean.prop = sum (Row.prop) / max(nreps))   %>% 
			mutate(mean.prop=if_else(is.nan(mean.prop),0,mean.prop)) %>% 
			group_by (!!OTU_column) %>%
			mutate (Colmax = max (mean.prop),
							Normalized.reads = mean.prop / Colmax) %>% 
			dplyr::select( -Colmax, -mean.prop) -> output 
		return(output)
	}
}


collapse <- function(tax, int, otu) {
	l1 <- tax_to_vector(tax) %>% duplicated()
	l2 <- tax_to_vector(tax)[l1]
	l2 <- unique(l2)
	
	for (i in 1:length(l2)) {
		l4 <- l2[i]
		l3 <- tax_to_vector(tax) %in% l4
		
		otu[which(tax_to_vector(tax) %in% l4)[1], ] <-
			otu[l3, ] %>% colSums()
		otu <- otu[-which(tax_to_vector(tax) %in% l4)[-1], ]
		
		int[which(tax_to_vector(tax) %in% l4)[1], ] <-
			int[l3, ] %>% sapply(., function(x)
				paste(x, collapse = " | "))
		int <- int[-which(tax_to_vector(tax) %in% l4)[-1], ]
		
		tax <- tax[-which(tax_to_vector(tax) %in% l4)[-1], ]
	}
	
	return(cbind(tax, int, otu))
}

collapse_otu2 <- function (otu, metadata = mmn,fun="Sum") 
{
	cat("Is collapse.name column created in the metadata?")
	bb <- as.data.frame(otu)
	colnames(bb) <- metadata[, "collapse.name"]
	otu_coll <- data.frame(matrix(0, nrow(bb), length(unique(colnames(bb)))))
	colnames(otu_coll) <- unique(colnames(bb))
	if (fun=="Sum") {
		for (i in unique(colnames(bb))) {
			otu_coll[, colnames(otu_coll) == i] <- rowSums(bb[colnames(bb) == i])
		} 
	}
	if (fun=="Mean") {
		for (i in unique(colnames(bb))) {
			otu_coll[, colnames(otu_coll) == i] <- rowMeans(bb[colnames(bb) == i])
		} 
	}
	return(as.data.frame(otu_coll))
}

# Libraries -------------------------------------------------------------------------


#Install libraries
library(ggu.base.fun)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(vegan)
library(MoMAColors)
library(indicspecies)

overview <- data.frame(
	Year = c(2018, 2019, 2019, 2020, 2020, 2021, 2021, 2021, 2022, 2022),
	Time = c("October", "May", "October", "June", "October", "March", "June", "October", "March", "September"),
	Cruise_1 = c("2018625", "2019631", "2019629", "2020611", "2020620", "2021607", "2021614", "2021624", "2022607", "2022624"),
	Cruise_2 = c("", "", "", "", "", "", "", "2021512", "", ""),
	Balsfjord = c(4, 4, 4, 3, 4, 10, 5, 4, 7, 0),
	Frakfjord = c(0, 0, 0, 3, 0, 4, 3, 4, 4, 4),
	Olderfjord = c(0, 0, 0, 1, 0, 3, 1, 3, 3, 3),
	Bergsfjord = c(0, 0, 0, 6, 0, 8, 6, 0, 6, 6)
)

overview

# setwd("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III")
# 
# # mm <- read.csv("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Documents & Protocolls/ALL_FISHDIV_samples_metadata_curated.csv")
# mm <- read.csv("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/Metadata.csv")
# 
# exp_idx <- read.csv("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/exposure_index.csv")
# exp_idx <- read.csv("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/exposure_index_100km.csv")
# exp_idx <- exp_idx %>% group_by(St_name) %>% 
# 	summarize(exp_idx=sum(Distance)) %>% 
# 	as.data.frame()
# 
# setwd("./5_collapsed_blasted_annotated_collapsed")
# list.files(pattern = "blasted_annotated_collapsed.csv")
# 
# Plate_A <- read.csv("Plate_A_collapsed_blasted_annotated_collapsed.csv")
# tx1 <- Plate_A[colnames(Plate_A)%in%c("kingdom","phylum","class","order","family","genus","species")] %>% mutate(plate="Plate_A")
# li1 <- paste0(tx1[,1],tx1[,2],tx1[,3],tx1[,4],tx1[,5],tx1[,6],tx1[,7])
# 
# Plate_B <- read.csv("Plate_B_collapsed_blasted_annotated_collapsed.csv")
# tax_merge <- Plate_B[colnames(Plate_B)%in%c("kingdom","phylum","class","order","family","genus","species")]%>% mutate(plate="Plate_B")
# li2 <- paste0(tax_merge[,1],tax_merge[,2],tax_merge[,3],tax_merge[,4],tax_merge[,5],tax_merge[,6],tax_merge[,7])
# 
# tax_merge <-	tx1 %>% filter(li1%notin%li2) %>% rbind(.,tax_merge) %>% arrange(class,order,family,species)
# 
# 
# Plate_CD <- read.csv("Plate_CD_collapsed_blasted_annotated_collapsed.csv")
# tx1 <- Plate_CD[colnames(Plate_CD)%in%c("kingdom","phylum","class","order","family","genus","species")] %>% mutate(plate="Plate_CD")
# li1 <- paste0(tx1[,1],tx1[,2],tx1[,3],tx1[,4],tx1[,5],tx1[,6],tx1[,7])
# 
# li2 <- paste0(tax_merge[,1],tax_merge[,2],tax_merge[,3],tax_merge[,4],tax_merge[,5],tax_merge[,6],tax_merge[,7])
# tax_merge <-	tx1 %>% filter(li1%notin%li2) %>% rbind(.,tax_merge) %>% arrange(class,order,family,species)
# 
# 
# 	
# Plate_NANB <- read.csv("Plate_NANB_collapsed_blasted_annotated_collapsed.csv")
# tx1 <- Plate_NANB[colnames(Plate_NANB)%in%c("kingdom","phylum","class","order","family","genus","species")] %>% mutate(plate="Plate_NANB")
# li1 <- paste0(tx1[,1],tx1[,2],tx1[,3],tx1[,4],tx1[,5],tx1[,6],tx1[,7])
# 
# li2 <- paste0(tax_merge[,1],tax_merge[,2],tax_merge[,3],tax_merge[,4],tax_merge[,5],tax_merge[,6],tax_merge[,7])
# tax_merge <-	tx1 %>% filter(li1%notin%li2) %>% rbind(.,tax_merge) %>% arrange(class,order,family,species)
# 
# 
# 
# Plate_NE <- read.csv("Plate_NE_collapsed_blasted_annotated_collapsed.csv")
# tx1 <- Plate_NE[colnames(Plate_NE)%in%c("kingdom","phylum","class","order","family","genus","species")] %>% mutate(plate="Plate_NE")
# li1 <- paste0(tx1[,1],tx1[,2],tx1[,3],tx1[,4],tx1[,5],tx1[,6],tx1[,7])
# 
# li2 <- paste0(tax_merge[,1],tax_merge[,2],tax_merge[,3],tax_merge[,4],tax_merge[,5],tax_merge[,6],tax_merge[,7])
# tax_merge <-	tx1 %>% filter(li1%notin%li2) %>% rbind(.,tax_merge) %>% arrange(class,order,family,species)
# 
# 
# 
# Plate_NF <- read.csv("Plate_NF_collapsed_blasted_annotated_collapsed.csv")
# tx1 <- Plate_NF[colnames(Plate_NF)%in%c("kingdom","phylum","class","order","family","genus","species")] %>% mutate(plate="Plate_NF")
# li1 <- paste0(tx1[,1],tx1[,2],tx1[,3],tx1[,4],tx1[,5],tx1[,6],tx1[,7])
# 
# li2 <- paste0(tax_merge[,1],tax_merge[,2],tax_merge[,3],tax_merge[,4],tax_merge[,5],tax_merge[,6],tax_merge[,7])
# tax_merge <-	tx1 %>% filter(li1%notin%li2) %>% rbind(.,tax_merge) %>% arrange(class,order,family,species)
# 
# 
# 
# # Plate_NGNH <- read.csv("Plate_NGNH_collapsed_blasted_annotated_collapsed.csv")
# # tx1 <- Plate_NGNH[colnames(Plate_NGNH)%in%c("kingdom","phylum","class","order","family","genus","species")] %>% mutate(plate="Plate_NGNH")
# # li1 <- paste0(tx1[,1],tx1[,2],tx1[,3],tx1[,4],tx1[,5],tx1[,6],tx1[,7])
# # 
# # li2 <- paste0(tax_merge[,1],tax_merge[,2],tax_merge[,3],tax_merge[,4],tax_merge[,5],tax_merge[,6],tax_merge[,7])
# # tax_merge <-	tx1 %>% filter(li1%notin%li2) %>% rbind(.,tax_merge) %>% arrange(class,order,family,species)
# 
# 
# 
# Plate_NI <- read.csv("Plate_NI_collapsed_blasted_annotated_collapsed.csv")
# tx1 <- Plate_NI[colnames(Plate_NI)%in%c("kingdom","phylum","class","order","family","genus","species")] %>% mutate(plate="Plate_NI")
# li1 <- paste0(tx1[,1],tx1[,2],tx1[,3],tx1[,4],tx1[,5],tx1[,6],tx1[,7])
# 
# li2 <- paste0(tax_merge[,1],tax_merge[,2],tax_merge[,3],tax_merge[,4],tax_merge[,5],tax_merge[,6],tax_merge[,7])
# tax_merge <-	tx1 %>% filter(li1%notin%li2) %>% rbind(.,tax_merge) %>% arrange(class,order,family,species)
# 
# 
# 
# Plate_NJ <- read.csv("Plate_NJ_collapsed_blasted_annotated_collapsed.csv")
# tx1 <- Plate_NJ[colnames(Plate_NJ)%in%c("kingdom","phylum","class","order","family","genus","species")] %>% mutate(plate="Plate_NJ")
# li1 <- paste0(tx1[,1],tx1[,2],tx1[,3],tx1[,4],tx1[,5],tx1[,6],tx1[,7])
# 
# li2 <- paste0(tax_merge[,1],tax_merge[,2],tax_merge[,3],tax_merge[,4],tax_merge[,5],tax_merge[,6],tax_merge[,7])
# tax_merge <-	tx1 %>% filter(li1%notin%li2) %>% rbind(.,tax_merge) %>% arrange(class,order,family,species)
# 
# 
# 
# Plate_NK <- read.csv("Plate_NK_collapsed_blasted_annotated_collapsed.csv")
# tx1 <- Plate_NK[colnames(Plate_NK)%in%c("kingdom","phylum","class","order","family","genus","species")] %>% mutate(plate="Plate_NK")
# li1 <- paste0(tx1[,1],tx1[,2],tx1[,3],tx1[,4],tx1[,5],tx1[,6],tx1[,7])
# 
# li2 <- paste0(tax_merge[,1],tax_merge[,2],tax_merge[,3],tax_merge[,4],tax_merge[,5],tax_merge[,6],tax_merge[,7])
# tax_merge <-	tx1 %>% filter(li1%notin%li2) %>% rbind(.,tax_merge) %>% arrange(class,order,family,species)
# 
# 
# 
# Plate_NL <- read.csv("Plate_NL_collapsed_blasted_annotated_collapsed.csv")
# tx1 <- Plate_NL[colnames(Plate_NL)%in%c("kingdom","phylum","class","order","family","genus","species")] %>% mutate(plate="Plate_NL")
# li1 <- paste0(tx1[,1],tx1[,2],tx1[,3],tx1[,4],tx1[,5],tx1[,6],tx1[,7])
# 
# li2 <- paste0(tax_merge[,1],tax_merge[,2],tax_merge[,3],tax_merge[,4],tax_merge[,5],tax_merge[,6],tax_merge[,7])
# tax_merge <-	tx1 %>% filter(li1%notin%li2) %>% rbind(.,tax_merge) %>% arrange(class,order,family,species)
# 
# 
# 
# Plate_NG01 <- read.csv("Plate_NG01_collapsed_blasted_annotated_collapsed.csv")
# tx1 <- Plate_NG01[colnames(Plate_NG01)%in%c("kingdom","phylum","class","order","family","genus","species")] %>% mutate(plate="Plate_NG01")
# li1 <- paste0(tx1[,1],tx1[,2],tx1[,3],tx1[,4],tx1[,5],tx1[,6],tx1[,7])
# 
# li2 <- paste0(tax_merge[,1],tax_merge[,2],tax_merge[,3],tax_merge[,4],tax_merge[,5],tax_merge[,6],tax_merge[,7])
# tax_merge <-	tx1 %>% filter(li1%notin%li2) %>% rbind(.,tax_merge) %>% arrange(class,order,family,species)
# 
# 
# 
# Plate_NG02 <- read.csv("Plate_NG02_collapsed_blasted_annotated_collapsed.csv")
# tx1 <- Plate_NG02[colnames(Plate_NG02)%in%c("kingdom","phylum","class","order","family","genus","species")] %>% mutate(plate="Plate_NG02")
# li1 <- paste0(tx1[,1],tx1[,2],tx1[,3],tx1[,4],tx1[,5],tx1[,6],tx1[,7])
# 
# li2 <- paste0(tax_merge[,1],tax_merge[,2],tax_merge[,3],tax_merge[,4],tax_merge[,5],tax_merge[,6],tax_merge[,7])
# tax_merge <-	tx1 %>% filter(li1%notin%li2) %>% rbind(.,tax_merge) %>% arrange(class,order,family,species)
# 
# 
# 
# Plate_NH01 <- read.csv("Plate_NH01_collapsed_blasted_annotated_collapsed.csv")
# tx1 <- Plate_NH01[colnames(Plate_NH01)%in%c("kingdom","phylum","class","order","family","genus","species")] %>% mutate(plate="Plate_NH01")
# li1 <- paste0(tx1[,1],tx1[,2],tx1[,3],tx1[,4],tx1[,5],tx1[,6],tx1[,7])
# 
# li2 <- paste0(tax_merge[,1],tax_merge[,2],tax_merge[,3],tax_merge[,4],tax_merge[,5],tax_merge[,6],tax_merge[,7])
# tax_merge <-	tx1 %>% filter(li1%notin%li2) %>% rbind(.,tax_merge) %>% arrange(class,order,family,species)
# 
# 
# 
# Plate_NH02 <- read.csv("Plate_NH02_collapsed_blasted_annotated_collapsed.csv")
# tx1 <- Plate_NH02[colnames(Plate_NH02)%in%c("kingdom","phylum","class","order","family","genus","species")] %>% mutate(plate="Plate_NH02")
# li1 <- paste0(tx1[,1],tx1[,2],tx1[,3],tx1[,4],tx1[,5],tx1[,6],tx1[,7])
# 
# li2 <- paste0(tax_merge[,1],tax_merge[,2],tax_merge[,3],tax_merge[,4],tax_merge[,5],tax_merge[,6],tax_merge[,7])
# tax_merge <-	tx1 %>% filter(li1%notin%li2) %>% rbind(.,tax_merge) %>% arrange(class,order,family,species)
# 
# 
# # Plate_A$species[Plate_A$species%in%"Somateria fischeri"] <- "Somateria spp."
# # Plate_A$species[Plate_A$id%in%"GL01_000340199"] <- "Ammodytes marinus|Hyperoplus lanceolatus|Ammodytes tobianus"
# # Plate_A$species[Plate_A$id%in%"GL01_000000013"] <- "Anarhichas spp."
# # Plate_A$family[Plate_A$id%in%"GL01_000000051"] <- "Gadidae"
# # write.csv(Plate_A,"Plate_A_collapsed_blasted_annotated_collapsed.csv",row.names = F)
# # 
# # Plate_NANB$genus[Plate_NANB$id%in%"GG01_001663732"] <- "Poecilia"
# # Plate_NANB$species[Plate_NANB$id%in%"GG01_001663732"] <- "Poecilia sphenops"
# # Plate_NANB$species[Plate_NANB$id%in%"GG01_000324077"] <- "Poecilia sphenops"
# # Plate_NANB$species[Plate_NANB$id%in%"GG01_002657274"] <- "Gadus morhua"
# # Plate_NANB$order[Plate_NANB$id%in%"GG01_002637929"] <- "Scorpaeniformes"
# # Plate_NANB$family[Plate_NANB$id%in%"GG01_002637929"] <- "Scorpaenidae"
# # Plate_NANB$genus[Plate_NANB$id%in%"GG01_002637929"] <- "Sebastes"
# # Plate_NANB$species[Plate_NANB$id%in%"GG01_002637929"] <- "Sebastes spp."
# # write.csv(Plate_NANB,"Plate_NANB_collapsed_blasted_annotated_collapsed.csv",row.names = F)
# # 
# # Plate_NK$class[Plate_NK$id%in%"GG01_002480775"] <- "Echinoidea"
# # write.csv(Plate_NK,"Plate_NK_collapsed_blasted_annotated_collapsed.csv",row.names = F)
# # 
# # Plate_NF$species[Plate_NF$id%in%"GG01_008064849"] <- "Clupea harengus"
# # Plate_NF$[Plate_NF$id%in%"GG01_001338213"] <- "Boreogadus saida"
# # Plate_NF$genus[Plate_NF$id%in%"GG01_001338213"] <- "Boreogadus"
# # Plate_NF$species[Plate_NF$id%in%"GG01_001338213"] <- "Boreogadus saida"
# # Plate_NF$order[Plate_NF$id%in%"GG01_001338213"] <- "Perciformes"
# # Plate_NF$species[Plate_NF$id%in%"GG01_001014766"] <- "Gadus morhua"
# # Plate_NF$order[Plate_NF$id%in%"GG01_001338213"] <- "Gadiformes"
# # Plate_NF$order[Plate_NF$id%in%"GG01_005498532"] <- "Perciformes"
# # Plate_NF$order[Plate_NF$id%in%"GG01_003221632"] <- "Scorpaeniformes"
# # Plate_NF$family[Plate_NF$id%in%"GG01_003221632"] <- "Scorpaenidae"
# # Plate_NF$genus[Plate_NF$id%in%"GG01_003221632"] <- "Sebastes"
# # Plate_NF$species[Plate_NF$id%in%"GG01_003221632"] <- "Sebastes spp."
# # write.csv(Plate_NF,"Plate_NF_collapsed_blasted_annotated_collapsed.csv",row.names = F)
# # 
# # Plate_NI$order[Plate_NI$id%in%"GG01_003811191"] <- "Scorpaeniformes"
# # Plate_NI$family[Plate_NI$id%in%"GG01_003811191"] <- "Scorpaenidae"
# # Plate_NI$genus[Plate_NI$id%in%"GG01_003811191"] <- "Sebastes"
# # Plate_NI$species[Plate_NI$id%in%"GG01_003811191"] <- "Sebastes spp."
# # write.csv(Plate_NI,"Plate_NI_collapsed_blasted_annotated_collapsed.csv",row.names = F)
# # 
# # Plate_NE$order[Plate_NE$id%in%"GG01_000859848"] <- "Scorpaeniformes"
# # Plate_NE$family[Plate_NE$id%in%"GG01_000859848"] <- "Scorpaenidae"
# # Plate_NE$genus[Plate_NE$id%in%"GG01_000859848"] <- "Sebastes"
# # Plate_NE$species[Plate_NE$id%in%"GG01_000859848"] <- "Sebastes spp."
# # write.csv(Plate_NE,"Plate_NE_collapsed_blasted_annotated_collapsed.csv",row.names = F)
# # 
# # Plate_NJ$order[Plate_NJ$id%in%"GG01_006625862"] <- "Perciformes|Scorpaeniformes"
# # Plate_NJ$genus[Plate_NJ$id%in%"GG01_000178580"] <- "Clupea"
# # Plate_NJ$species[Plate_NJ$id%in%"GG01_000178580"] <- "Clupea harengus"
# # Plate_NJ$species[Plate_NJ$id%in%"GG01_003597153"] <- "Anarhichas spp."
# # Plate_NJ$species[Plate_NJ$id%in%"GG01_004693283"] <- "Sebastes spp."
# # Plate_NJ$order[Plate_NJ$id%in%"GG01_004693283"] <- "Scorpaeniformes"
# # Plate_NJ$family[Plate_NJ$id%in%"GG01_004693283"] <- "Scorpaenidae"
# # write.csv(Plate_NJ,"Plate_NJ_collapsed_blasted_annotated_collapsed.csv",row.names = F)
# # 
# # Plate_NL$order[Plate_NL$id%in%"GG01_000917475"] <- "Pleuronectidae"
# # Plate_NL$genus[Plate_NL$id%in%"GG01_000917475"] <- "Pleuronectes"
# # Plate_NL$species[Plate_NL$id%in%"GG01_000917475"] <- "Pleuronectes platessa"
# # Plate_NL$order[Plate_NL$id%in%"GG01_000917475"] <- "Pleuronectiformes"
# # Plate_NL$species[Plate_NL$id%in%"GG01_002829028"] <- "Anarhichas spp."
# # Plate_NL$order[Plate_NL$id%in%"GG01_000550564"] <- "Scorpaeniformes"
# # Plate_NL$family[Plate_NL$id%in%"GG01_000550564"] <- "Scorpaenidae"
# # Plate_NL$species[Plate_NL$id%in%"GG01_000550564"] <- "Sebastes spp."
# # write.csv(Plate_NL,"Plate_NL_collapsed_blasted_annotated_collapsed.csv",row.names = F)
# # 
# # Plate_NK$genus[Plate_NK$id%in%"GG01_002960057"] <- "Clupea"
# # Plate_NK$species[Plate_NK$id%in%"GG01_002960057"] <- "Clupea harengus"
# # 
# # Plate_NK$order[Plate_NK$id%in%"GG01_002376339"] <- "Scorpaeniformes"
# # Plate_NK$family[Plate_NK$id%in%"GG01_002376339"] <- "Scorpaenidae"
# # Plate_NK$genus[Plate_NK$id%in%"GG01_002376339"] <- "Sebastes"
# # Plate_NK$species[Plate_NK$id%in%"GG01_002376339"] <- "Sebastes spp."
# # write.csv(Plate_NK,"Plate_NK_collapsed_blasted_annotated_collapsed.csv",row.names = F)
# 
# # Plate_NG01$species[Plate_NG01$id%in%"GG01_000372424"] <- "Glyptocephalus cynoglossus"
# # Plate_NG01$species[Plate_NG01$id%in%"GG01_000572092"] <- "Liparis fabricii|Liparis inquilinus|Liparis liparis"
# # Plate_NG01$species[Plate_NG01$id%in%"GG01_000382676"] <- "Sebastes spp."
# # Plate_NG01[Plate_NG01$id%in%"GG01_001813697",1:7] <- c("Animalia","Chordata", "Teleostei", "Scorpaeniformes","Scorpaenidae","Sebastes","Sebastes spp.")
# # write.csv(Plate_NG01,"Plate_NG01_collapsed_blasted_annotated_collapsed.csv",row.names = F)
# 
# # Plate_NG02$species[Plate_NG02$id%in%"GG02_000659545"] <- "Gadus morhua"
# # Plate_NG02[Plate_NG02$id%in%"GG02_000596027",1:7] <- c("Animalia","Chordata", "Teleostei", "Scorpaeniformes","Scorpaenidae","Sebastes","Sebastes spp.")
# # write.csv(Plate_NG02,"Plate_NG02_collapsed_blasted_annotated_collapsed.csv",row.names = F)
# 
# # Plate_NH01$species[Plate_NH01$id%in%"GG03_000007358"] <- "Myoxocephalus scorpius"
# # Plate_NH01[Plate_NH01$id%in%"GG03_000570544",1:7] <- c("Animalia","Chordata", "Teleostei", "Scorpaeniformes","Scorpaenidae","Sebastes","Sebastes spp.")
# # write.csv(Plate_NH01,"Plate_NH01_collapsed_blasted_annotated_collapsed.csv",row.names = F)
# 
# 
# # Plate_NH02$species[Plate_NH02$id%in%"GG04_000228285"] <- "Liparis fabricii|Liparis inquilinus|Liparis liparis"
# # Plate_NH02[Plate_NH02$id%in%"GG04_001201577",1:7] <- c("Animalia","Chordata", "Teleostei", "Scorpaeniformes","Scorpaenidae","Sebastes","Sebastes spp.")
# # write.csv(Plate_NH02,"Plate_NH02_collapsed_blasted_annotated_collapsed.csv",row.names = F)
# 
# 
# 
# Plate_A_tax <- Plate_A %>% 
# 	select(c("kingdom","phylum","class","order","family","genus","species"))
# Plate_A_int <- Plate_A %>% 
# 	select(c("id", "annotated_pident", "best_identity"))
# Plate_A_otu <- Plate_A %>% 
# 	select(-c("kingdom","phylum","class","order","family","genus","species","id", "annotated_pident", "best_identity"))
# Plate_A <- collapse(Plate_A_tax,Plate_A_int,Plate_A_otu)
# 
# Plate_NF_tax <- Plate_NF %>% 
# 	select(c("kingdom","phylum","class","order","family","genus","species"))
# Plate_NF_int <- Plate_NF %>% 
# 	select(c("id", "annotated_pident", "best_identity"))
# Plate_NF_otu <- Plate_NF %>% 
# 	select(-c("kingdom","phylum","class","order","family","genus","species","id", "annotated_pident", "best_identity"))
# Plate_NF <- collapse(Plate_NF_tax,Plate_NF_int,Plate_NF_otu)
# 
# Plate_NI_tax <- Plate_NI %>% 
# 	select(c("kingdom","phylum","class","order","family","genus","species"))
# Plate_NI_int <- Plate_NI %>% 
# 	select(c("id", "annotated_pident", "best_identity"))
# Plate_NI_otu <- Plate_NI %>% 
# 	select(-c("kingdom","phylum","class","order","family","genus","species","id", "annotated_pident", "best_identity"))
# Plate_NI <- collapse(Plate_NI_tax,Plate_NI_int,Plate_NI_otu)
# 
# Plate_NJ_tax <- Plate_NJ %>% 
# 	select(c("kingdom","phylum","class","order","family","genus","species"))
# Plate_NJ_int <- Plate_NJ %>% 
# 	select(c("id", "annotated_pident", "best_identity"))
# Plate_NJ_otu <- Plate_NJ %>% 
# 	select(-c("kingdom","phylum","class","order","family","genus","species","id", "annotated_pident", "best_identity"))
# Plate_NJ <- collapse(Plate_NJ_tax,Plate_NJ_int,Plate_NJ_otu)
# 
# Plate_NK_tax <- Plate_NK %>% 
# 	select(c("kingdom","phylum","class","order","family","genus","species"))
# Plate_NK_int <- Plate_NK %>% 
# 	select(c("id", "annotated_pident", "best_identity"))
# Plate_NK_otu <- Plate_NK %>% 
# 	select(-c("kingdom","phylum","class","order","family","genus","species","id", "annotated_pident", "best_identity"))
# Plate_NK <- collapse(Plate_NK_tax,Plate_NK_int,Plate_NK_otu)
# 
# Plate_NL_tax <- Plate_NL %>% 
# 	select(c("kingdom","phylum","class","order","family","genus","species"))
# Plate_NL_int <- Plate_NL %>% 
# 	select(c("id", "annotated_pident", "best_identity"))
# Plate_NL_otu <- Plate_NL %>% 
# 	select(-c("kingdom","phylum","class","order","family","genus","species","id", "annotated_pident", "best_identity"))
# Plate_NL <- collapse(Plate_NL_tax,Plate_NL_int,Plate_NL_otu)
# 
# Plate_NG01_tax <- Plate_NG01 %>% 
# 	select(c("kingdom","phylum","class","order","family","genus","species"))
# Plate_NG01_int <- Plate_NG01 %>% 
# 	select(c("id", "annotated_pident", "best_identity"))
# Plate_NG01_otu <- Plate_NG01 %>% 
# 	select(-c("kingdom","phylum","class","order","family","genus","species","id", "annotated_pident", "best_identity"))
# Plate_NG01 <- collapse(Plate_NG01_tax,Plate_NG01_int,Plate_NG01_otu)
# 
# Plate_NG02_tax <- Plate_NG02 %>% 
# 	select(c("kingdom","phylum","class","order","family","genus","species"))
# Plate_NG02_int <- Plate_NG02 %>% 
# 	select(c("id", "annotated_pident", "best_identity"))
# Plate_NG02_otu <- Plate_NG02 %>% 
# 	select(-c("kingdom","phylum","class","order","family","genus","species","id", "annotated_pident", "best_identity"))
# Plate_NG02 <- collapse(Plate_NG02_tax,Plate_NG02_int,Plate_NG02_otu)
# 
# Plate_NH01_tax <- Plate_NH01 %>% 
# 	select(c("kingdom","phylum","class","order","family","genus","species"))
# Plate_NH01_int <- Plate_NH01 %>% 
# 	select(c("id", "annotated_pident", "best_identity"))
# Plate_NH01_otu <- Plate_NH01 %>% 
# 	select(-c("kingdom","phylum","class","order","family","genus","species","id", "annotated_pident", "best_identity"))
# Plate_NH01 <- collapse(Plate_NH01_tax,Plate_NH01_int,Plate_NH01_otu)
# 
# Plate_NH02_tax <- Plate_NH02 %>% 
# 	select(c("kingdom","phylum","class","order","family","genus","species"))
# Plate_NH02_int <- Plate_NH02 %>% 
# 	select(c("id", "annotated_pident", "best_identity"))
# Plate_NH02_otu <- Plate_NH02 %>% 
# 	select(-c("kingdom","phylum","class","order","family","genus","species","id", "annotated_pident", "best_identity"))
# Plate_NH02 <- collapse(Plate_NH02_tax,Plate_NH02_int,Plate_NH02_otu)
# 
# 
# 
# merging <- c("kingdom"="kingdom",
# 						 "phylum"="phylum",
# 						 "class"="class",
# 						 "order"="order",
# 						 "family"="family",
# 						 "genus"="genus",
# 						 "species"="species")
# 
# x <- left_join(tax_merge,Plate_A,by=merging) %>% 
# 	left_join(.,Plate_B,by=merging) %>%
# 	left_join(.,Plate_CD,by=merging) %>% 
# 	left_join(.,Plate_NANB,by=merging) %>% 
# 	left_join(.,Plate_NE,by=merging) %>% 
# 	left_join(.,Plate_NF,by=merging) %>% 
# 	left_join(.,Plate_NG01,by=merging) %>% 
# 	left_join(.,Plate_NG02,by=merging) %>% 
# 	left_join(.,Plate_NH01,by=merging) %>% 
# 	left_join(.,Plate_NH02,by=merging) %>% 
# 	left_join(.,Plate_NI,by=merging) %>% 
# 	left_join(.,Plate_NJ,by=merging) %>% 
# 	left_join(.,Plate_NK,by=merging) %>% 
# 	left_join(.,Plate_NL,by=merging)
# 
# 
# 
# xn <- colnames(x)
# 
# xn <- xn %>% gsub("18KB70656","2018625_",.)
# xn <- xn %>% gsub("X2018625_0","X2018625_",.)
# xn <- xn %>% gsub("X19KB55241","X2019629_",.)
# xn <- xn %>% gsub("X2019629_0","X2019629_",.)
# xn <- xn %>% gsub("X2020611_0","X2020611_",.)
# xn <- xn %>% gsub("X2021607_00","X2021607_",.)
# xn <- xn %>% gsub("X2021607_0","X2021607_",.)
# xn <- gsub("(X)([1-2]\\.[LBEB]{2})", "X2018625_\\2", xn)
# xn <- gsub("(X)([1-2]EB[0-9]{3}521)","X2021607_\\2",xn)
# xn <- xn %>% gsub("\\.[xy]{1}","",.)
# xn[xn=="PCRblank_D1"] <- c("NTC_A_D1","NTC_B_D1")
# xn[xn=="PCRblank_B1"] <- c("NTC_A_B1","NTC_B_B1")
# 
# 
# xn <- make.unique.2(xn,sep = "-")
# 
# o <- xn[grepl("-",xn)][!grepl("id|annotated_pident|best_identity|tot_reads",xn[grepl("-",xn)])]
# o <- o[!grepl("\\-1",o)]
# 
# o <-
# 	o %>% gsub("\\-[0-9]","",.) %>% 
# 	as.data.frame(.) %>% setNames("Filter_name") %>% 
# 	left_join(.,mm,by="Filter_name") %>% 
# 	select(colnames(mm))
# 	
# mm <- rbind(mm,o)
# mm$Filter_name <- make.unique.2(mm$Filter_name,"-")
# 
# mm <- mm[mm$Filter_name%in%xn,]
# 
# 
# # xn[xn%in%mm$Filter_name] %>% length()
# 
# colnames(x) <- xn
# 
# # mm$Filter_name[mm$Filter_name%in%xn]
# # mm[mm$Filter_name%notin%xn,]
# 
# x <-
# 	x %>% 
# 	filter(class%in%c("Teleostei","Chondrichthyes")) %>% 
# 		filter(is.na(species) | species != "Poecilia sphenops") %>% 
# 		filter(is.na(species) | species != "Trigonostigma espei") %>%
# 		filter(is.na(species) | species != "Epiplatys dageti")
# 
# 
# tax <- x %>% 
# 	select(c("kingdom","phylum","class","order","family","genus","species"))
# otu <- x %>% 
# 	select(-colnames(.)[multi_grep(c("kingdom","phylum","class","order","family","genus","species","id", "annotated_pident", "best_identity","plate","tot_reads"),colnames(.))])
# 
# otu[is.na(otu)] <- 0
# 
# mm <- mm %>% mutate(Station_name=if_else(grepl("\\-2",Filter_name),paste0(Station_name,"_2"),Station_name))
# 
# # xx <- mm %>% filter(grepl("\\-1",Filter_name))
# # xx2 <- mm %>% filter(grepl("\\-2",Filter_name))
# # 
# # 
# # a_l <- which(substr(xx$Filter_name,0,nchar(xx$Filter_name)-2)%in%gsub("X2019629_0","X2019629_",gsub("X19KB55241","X2019629_",colnames(Plate_A))))
# # xx[a_l,]
# # xx$Plate[a_l] <- gsub("Plate_NB \\| ","",xx$Plate[a_l])
# # xx$Plate[a_l] <- gsub("Plate_NE \\| ","",xx$Plate[a_l])
# # xx$Plate[a_l] <- gsub("Plate_NC \\| ","",xx$Plate[a_l])
# # 
# # cd_l <- which(substr(xx$Filter_name,0,nchar(xx$Filter_name)-2)%in%colnames(Plate_CD))
# # xx[cd_l,]
# # xx$Plate[cd_l] <- gsub("Plate_NA \\| ","",xx$Plate[cd_l])
# # xx$Plate[cd_l] <- gsub("Plate_NF \\| ","",xx$Plate[cd_l])
# # 
# # nanb_l <- which(substr(xx$Filter_name,0,nchar(xx$Filter_name)-2)%in%colnames(Plate_NANB))
# # xx[nanb_l,]
# # xx$Plate[nanb_l] <- gsub("Plate_NA \\| ","",xx$Plate[nanb_l])
# # 
# # a2_l <- which(substr(xx$Filter_name,0,nchar(xx$Filter_name)-2)%in%gsub("X2018625_0","X2018625_",gsub("18KB70656","2018625_",colnames(Plate_A))))
# # xx[a2_l,]
# # xx$Plate[a2_l] <- gsub("Plate_NF \\| ","",xx$Plate[a2_l])
# # xx$Plate[a2_l] <- gsub("Plate_NE \\| ","",xx$Plate[a2_l])
# # xx$Plate[a2_l] <- gsub("Plate_ND \\| ","",xx$Plate[a2_l])
# # 
# # 
# # ne_l <- which(substr(xx$Filter_name,0,nchar(xx$Filter_name)-2)%in%colnames(Plate_NE))
# # xx[ne_l,]
# # xx$Plate[ne_l] <- gsub("Plate_NE \\| ","",xx$Plate[ne_l])
# # 
# # nf_l <- which(substr(xx$Filter_name,0,nchar(xx$Filter_name)-2)%in%colnames(Plate_NF))
# # xx[nf_l,]
# # xx$Plate[nf_l] <- gsub("Plate_NF \\| ","",xx$Plate[nf_l])
# # 
# # a_l <- which(substr(xx2$Filter_name,0,nchar(xx2$Filter_name)-2)%in%gsub("X2019629_0","X2019629_",gsub("X19KB55241","X2019629_",colnames(Plate_A))))
# # xx2[a_l,]
# # xx2$Plate[a_l] <- gsub(" \\| Plate_A","",xx2$Plate[a_l])
# # xx2$Plate[a_l] <- gsub("Plate_NC \\| ","",xx2$Plate[a_l])
# # 
# # cd_l <- which(substr(xx2$Filter_name,0,nchar(xx2$Filter_name)-2)%in%colnames(Plate_CD))
# # xx2[cd_l,]
# # xx2$Plate[cd_l] <- gsub(" \\| Plate_D","",xx2$Plate[cd_l])
# # xx2$Plate[cd_l] <- gsub(" \\| Plate_C","",xx2$Plate[cd_l])
# # 
# # nanb_l <- which(substr(xx2$Filter_name,0,nchar(xx2$Filter_name)-2)%in%colnames(Plate_NANB))
# # xx2[nanb_l,]
# # xx2$Plate[nanb_l] <- gsub(" \\| Plate_D","",xx2$Plate[nanb_l])
# # xx2$Plate[nanb_l] <- gsub(" \\| Plate_C","",xx2$Plate[nanb_l])
# # 
# # a2_l <- which(substr(xx2$Filter_name,0,nchar(xx2$Filter_name)-2)%in%gsub("X2018625_0","X2018625_",gsub("18KB70656","2018625_",colnames(Plate_A))))
# # xx2[a2_l,]
# # xx2$Plate[a2_l] <- gsub(" \\| Plate_A","",xx2$Plate[a2_l])
# # xx2$Plate[a2_l] <- gsub("Plate_ND \\| ","",xx2$Plate[a2_l])
# # 
# # 
# # ne_l <- which(substr(xx2$Filter_name,0,nchar(xx2$Filter_name)-2)%in%colnames(Plate_NE))
# # xx2[ne_l,]
# # xx2$Plate[ne_l] <- gsub(" \\| Plate_A","",xx2$Plate[ne_l])
# # 
# # nf_l <- which(substr(xx2$Filter_name,0,nchar(xx2$Filter_name)-2)%in%colnames(Plate_NF))
# # xx2[nf_l,]
# # xx2$Plate[nf_l] <- gsub(" \\| Plate_A","",xx2$Plate[nf_l])
# # xx2$Plate[nf_l] <- gsub(" \\| Plate_C","",xx2$Plate[nf_l])
# # 
# # xx %>% distinct(Plate)
# # xx2 %>% distinct(Plate)
# # 
# # mm[grepl("\\-1",mm$Filter_name),] <- xx
# # mm[grepl("\\-2",mm$Filter_name),] <- xx2
# # 
# # 
# # 
# # xx <- mm %>% filter(grepl(" \\| ",Plate))
# # xx$Plate <- gsub("Plate_NC \\| ","",xx$Plate)
# # xx$Plate <- gsub("Plate_ND \\| ","",xx$Plate)
# # xx$Plate <- gsub(" \\| Plate_E","",xx$Plate)
# # xx$Plate[xx$Filter_name%in%c("X2018625_6","X2018625_30")] <- gsub("Plate_NE \\| ","",xx$Plate[xx$Filter_name%in%c("X2018625_6","X2018625_30")])
# # 
# # mm[grepl(" \\| ",mm$Plate),] <- xx
# # xx <- mm %>% filter()
# # 
# # write.csv(mm,"metadata_curated.csv",row.names = F)
# 
# mm <- read.csv("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/Metadata_curated.csv")
# mm <-
# 	mm %>% 
# 	mutate(Season=substr(start = 4,stop = 5,Date)) %>% 
# 	mutate(Season=if_else(Season=="03","Spring",Season)) %>% 
# 	mutate(Season=if_else(Season=="05","Summer",Season)) %>% 
# 	mutate(Season=if_else(Season=="06","Summer",Season)) %>% 
# 	mutate(Season=if_else(Season=="09","Autumn",Season)) %>% 
# 	mutate(Season=if_else(Season=="10","Autumn",Season)) %>% 
# 	mutate(Season = factor(Season, levels = c("Spring", "Summer", "Autumn"))) %>% 
# 	mutate(Location=substr(start = 0,stop = 2,Station_name)) %>% 
# 	mutate(Location=if_else(grepl("\\-2",Filter_name),paste0(Location,"_2"),Location)) %>% 
# 	mutate(Location = factor(Location, levels = c("BA", "FR", "OL","BE","BA_2","FR_2"))) %>%
# 	mutate(Month=substr(start = 4,stop = 5,Date)) %>% 
# 	mutate(Year=substr(start = 7,stop = 10,Date)) %>% 
# 	mutate(Date_2=paste0(Year,"_",Month))


# Import data -----------------------------------------------------------------------

otu <- read.csv("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/6_OTU_TAX_META/otu.csv",check.names = FALSE)
tax <- read.csv("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/6_OTU_TAX_META/tax.csv")
mm <- read.csv("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/6_OTU_TAX_META/Metadata_curated.csv")


# Prepare the data ------------------------------------------------------------------

otu_coll <- collapse_otu(otu %>% select(mm$Filter_name),
															mm %>% mutate(collapse.name=paste0(Cruise,"_",Station_name,"_",Sample_type)))
mm_coll <- collapse_meta(mm %>% mutate(collapse.name=paste0(Cruise,"_",Station_name,"_",Sample_type)))

#Filter the duplicated stations (BA_2 and FR_2)
mm_coll <- mm_coll %>% 
	filter(!(Location=="BA"&Cruise_no==2019629)) %>% 
	filter(!(Location=="BA"&Cruise_no==2018625)) %>% 
	filter(!(Location=="FR"&Cruise_no==2021607))

otu_coll <- otu_coll %>% select(mm_coll$collapse.name)

mm_coll <- mm_coll %>% 
	mutate(Location=if_else(grepl("_2",Location),substr(Location,0,2),Location)) 

	
# Alpha diversity analysis ----------------------------------------------------------
alpha <- pa_conversion(otu_coll) %>% colSums() %>% as.data.frame() %>% 
	rownames_to_column("Sample") %>% setNames(c("Sample","Richness"))

alpha <- mm_coll %>% mutate(Rich=merr2(collapse.name,alpha$Sample,alpha$Richness)) %>% 
	mutate(Rich = as.numeric(Rich))

# Filter the blanks
alpha <- alpha %>% filter(sample_type=="e_sample")

alpha_summary <-
	alpha %>%
	group_by(Location,Season, Year) %>%
	summarise(richness_mean = mean(Rich, na.rm = TRUE),
						q_25=quantile(Rich,0.25),
						q_75=quantile(Rich,0.75), .groups = 'drop')


# # Figure 2 - Alpha metrics ----------------------------------------------------------
# 
# cc <- c("BA" = "deepskyblue2",
# 				"FR" = "orange2",
# 				"OL" = "tomato2",
# 				"BE" = "green3")
# 
# p1 <-
# 	ggplot() +
# 	geom_jitter(data = alpha %>% filter(Season=="Autumn"),
# 							aes(x = as.character(Year), y = Rich, color=Location), position=position_dodge(width=0.4),alpha=0.5) +
# 	theme_classic() +
# 	geom_errorbar(data = alpha_summary %>% 	filter(Season=="Autumn"),
# 								aes(x = as.character(Year), ymin = q_25, ymax = q_75,color=Location),width = 0.3,
# 								position = position_dodge(width = 0.4))+
# 	geom_point(data = alpha_summary %>% 	filter(Season=="Autumn"),
# 						 aes(x = as.character(Year), y = richness_mean, color=Location),
# 						 position = position_dodge(width = 0.4), size = 6) + # Nudge points to the right
# 	ylim(3, 43)+
# 	ylab("Species richness")+
# 	xlab("")+
# 	scale_color_manual(values=cc)+
# 	facet_wrap(~Season)+
# 	theme(legend.position = "none",
# 				axis.title.y = element_text(size = 19),
# 				plot.margin = margin(0.1, 0, 0.5, 0.5, "cm"),
# 				axis.title.x = element_text(size = 19),
# 				legend.title = element_text(size=17),
# 				legend.key.height = unit(1.0, 'cm'),
# 				legend.text = element_text(size=16),
# 				axis.text.x = element_text(size = 17),
# 				strip.text = element_text(size = 17, color="white"),
# 				strip.background = element_rect(fill = "deepskyblue4", color = "black", size = 0.5, linetype = "dotted"),
# 				axis.text.y=element_text(size=15))
# p2 <-
# 	ggplot() +
# 	geom_jitter(data = alpha %>% filter(Season=="Spring"),
# 							aes(x = as.character(Year), y = Rich, color=Location), position=position_dodge(width=0.4),alpha=0.5) +
# 	theme_classic() +
# 	geom_errorbar(data = alpha_summary %>% 	filter(Season=="Spring"),
# 								aes(x = as.character(Year), ymin = q_25, ymax = q_75,color=Location),width = 0.3,
# 								position = position_dodge(width = 0.4))+
# 	geom_point(data = alpha_summary %>% 	filter(Season=="Spring"),
# 						 aes(x = as.character(Year), y = richness_mean, color=Location),
# 						 position = position_dodge(width = 0.4), size = 6) + # Nudge points to the right
# 	ylim(3, 43)+
# 	ylab("Species richness")+
# 	xlab("")+
# 	scale_color_manual(values=cc)+
# 	facet_wrap(~Season)+
# 	theme(legend.position = "none",
# 				axis.title.y = element_blank(),
# 				plot.margin = margin(0.1, 0.0, 0.5, 0.5, "cm"),
# 				axis.title.x = element_text(size = 19),
# 				legend.title = element_text(size=17),
# 				legend.key.height = unit(1.0, 'cm'),
# 				legend.text = element_text(size=16),
# 				axis.text.x = element_text(size = 17),
# 				axis.text.y = element_blank(),
# 				axis.line.y = element_blank(),
# 				strip.background = element_rect(fill = "orange2", color = "black", size = 0.5, linetype = "dotted"),
# 				strip.text = element_text(size = 17),
# 				axis.ticks.y = element_blank())
# 
# p3 <-
# 	ggplot() +
# 	geom_jitter(data = alpha %>% filter(Season=="Summer"),
# 							aes(x = as.character(Year), y = Rich, color=Location), position=position_dodge(width=0.4),alpha=0.5) +
# 	theme_classic() +
# 	geom_errorbar(data = alpha_summary %>% filter(Season=="Summer"),
# 								aes(x = as.character(Year), ymin = q_25, ymax = q_75,color=Location),width = 0.3,
# 								position = position_dodge(width = 0.4))+
# 	geom_point(data = alpha_summary %>% filter(Season=="Summer"),
# 						 aes(x = as.character(Year), y = richness_mean, color=Location),
# 						 position = position_dodge(width = 0.4), size = 6) + # Nudge points to the right
# 	ylim(3, 43)+
# 	ylab("Species richness")+
# 	xlab("")+
# 	scale_color_manual(values=cc)+
# 	facet_wrap(~Season)+
# 	theme(legend.position = "none",
# 				axis.title.y = element_blank(),
# 				plot.margin = margin(0.1, -0.5, 0.5, 0.5, "cm"),
# 				axis.title.x = element_text(size = 19),
# 				legend.title = element_text(size=17),
# 				legend.key.height = unit(1.0, 'cm'),
# 				legend.text = element_text(size=16),
# 				axis.text.x = element_text(size = 17),
# 				axis.text.y = element_blank(),
# 				axis.line.y = element_blank(),
# 				strip.text = element_text(size = 17),
# 				strip.background = element_rect(fill = "lightblue", color = "black", size = 0.5, linetype = "dotted"),
# 				axis.ticks.y = element_blank())
# 
# pp1 <- cowplot::plot_grid(p1,p2,p3,align = "h",rel_widths = c(5,2,3),ncol = 3)
# 
# cc1 <- c("Balfsfjord" = "deepskyblue2",
# 				"Frakkfjord" = "orange2",
# 				"Olderfjord" = "tomato2",
# 				"Bergsfjord" = "green3")
# 
# p1_leg <-
# 	ggplot() +
# 	theme_classic() +
# 	geom_errorbar(data = alpha_summary %>%
# 									mutate(Location=if_else(Location=="BA","Balfsfjord",Location)) %>%
# 									mutate(Location=if_else(Location=="FR","Frakkfjord",Location)) %>%
# 									mutate(Location=if_else(Location=="OL","Olderfjord",Location)) %>%
# 									mutate(Location=if_else(Location=="BE","Bergsfjord",Location)), aes(x = as.character(Year), ymin = q_25, ymax = q_75,color=Location))+
# 	geom_point(data = alpha_summary %>%
# 						 	mutate(Location=if_else(Location=="BA","Balfsfjord",Location)) %>%
# 						 	mutate(Location=if_else(Location=="FR","Frakkfjord",Location)) %>%
# 						 	mutate(Location=if_else(Location=="OL","Olderfjord",Location)) %>%
# 						 	mutate(Location=if_else(Location=="BE","Bergsfjord",Location)),
# 						 aes(x = as.character(Year), y = richness_mean, color=Location),size = 6) + # Nudge points to the right
# 	scale_color_manual("Location",values=cc1)+
# 	theme(plot.margin = margin(1, 1, 1, 0, "cm"),
# 				legend.key.width = unit(1.2,"cm"),
# 				legend.title = element_text(size=18),
# 				legend.key.height = unit(0.8, 'cm'),
# 				legend.text = element_text(size=16))
# 
# p1_leg_leg <- cowplot::get_legend(p1_leg+
# 																		theme(legend.justification = c(0.5,0.5)))
# 
# # pp2 <- cowplot::plot_grid(pp1,p1_leg_leg,rel_widths = c(8,1.4),ncol = 2)
# # ggsave(paste0("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/Plots",
# # 							"/Figure_2_Alpha diversity.jpg"),pp2,width=35,height=25,units="cm")



# Multivariate analysis -------------------------------------------------------------

#Filter for empty columns
otu_coll_nmds <- otu_coll %>% select(colnames(.)[colSums(.)>0])
mm_coll_nmds <- mm_coll[mm_coll$collapse.name%in%colnames(otu_coll),]

# 
# dist.rel_0<-vegdist(t(sqrt(sqrt(otu_coll_nmds))), method = "bray")
# nmds.rel_0 <- metaMDS(dist.rel_0, k = dime, autotransform = T)
# mm_coll_nmds_0 <- mm_coll_nmds

#Filter for only environmental samples
mm_coll_nmds <- mm_coll_nmds %>% filter((grepl("e_sample",sample_type)))
otu_coll_nmds <- otu_coll_nmds %>% select(colnames(.)[colnames(.)%in%mm_coll_nmds$collapse.name])

#Filter for species that have no occurances
tax_nmds <- tax %>% filter(rowSums(otu_coll_nmds)>0)
otu_coll_nmds <- otu_coll_nmds %>% filter(rowSums(.)>0)

# Add grouping for Locaition - Season - Year
mm_coll_nmds <- mm_coll_nmds %>% mutate(group=paste0(Location,"_",Season,"_",Year))

# otu_coll_nmds <- otu_coll_nmds[rowSums(otu_coll_nmds)>1000,]


# If you want to use biological replicates as an individual sample
totu <- trans(otu_coll_nmds) %>%  # this line transposes the otu table with X Y Z columns
	mutate(Rep=1) %>%
	setNames(c("Species","Sample","Reads","Rep")) #this line renames the columns to "Species","Sample", and "Reads"

# If you want to pool/collapse the biological replicates
totu_coll <- trans(otu_coll_nmds) %>% #this line transposes the otu table with X Y Z columns
	mutate(Rep=substr(Y,nchar(Y),nchar(Y))) %>% #this line removes the last character that indicates the replicates (e.g., "_1","_2","_3")
	mutate(Y=substr(Y,0,nchar(Y)-2)) %>% #this line removes the last character that indicates the replicates (e.g., "_1","_2","_3")
	setNames(c("Species","Sample","Reads","Rep")) #this line renames the columns to "Species","Sample", and "Reads"


# Run the eDNA index function -------------------------------------------------------

# This line runs converts the reads into eDNA index
totu_indexed_nmds <- eDNAindex(df = totu,Sample_column = Sample,OTU_column = Species,Counts_column = Reads,Biological.replicate = Rep)


# Re-transpose the otu table to rows=species and column=station
otu_indexed_nmds <- inv.trans(input.X = totu_indexed_nmds$Species,
												 input.Y = totu_indexed_nmds$Sample,
												 input.Z = totu_indexed_nmds$Normalized.reads)

dime=2

# dist.sq<-vegdist(t(sqrt(otu_coll_nmds)), method = "bray")
# nmds.sq <- metaMDS(dist.sq, k = dime, autotransform = T)

dist.sq.sq<-vegdist(t(sqrt(sqrt(otu_coll_nmds))), method = "bray")
nmds.sq.sq <- metaMDS(dist.sq.sq, k = dime, autotransform = T,trymax = 99,try = 2)

# dist.rel.sq<-vegdist(t(rel.ab(sqrt(otu_coll_nmds))), method = "bray")
# nmds.rel.sq <- metaMDS(dist.rel.sq, k = dime, autotransform = T)
# 
dist.rel<-vegdist(t(rel.ab(otu_coll_nmds)), method = "bray")
nmds.rel <- metaMDS(dist.rel, k = dime, autotransform = T,trymax = 999)
# 
dist.idx<-vegdist(t(otu_indexed_nmds), method = "bray")
nmds.idx <- metaMDS(dist.idx, k = dime, autotransform = T,trymax = 99)
# 
# dist.lg<-vegdist(t(log(otu_coll_nmds+1)), method = "bray")
# nmds.lg <- metaMDS(dist.lg, k = dime, autotransform = T)


par(mar=c(1,1,1,1))


# # Figurexxx - e_sample vs blanks -----------------------------------------------------
# 
# jpeg(paste0("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/Plots/",
# 						"Difference between blank and samples_bc.jpg"),
# 		 width = 2000, height = 1900)
# 
# group <- mm_coll_nmds %>% distinct(sample_type) %>% mutate(group=if_else(sample_type=="e_sample","Sample","Blank"))
# label <- group %>% pull(group)
# 
# co <- c("tomato2","deepskyblue2")
# co1 <- co[factor(group$group)]
# group <- cbind(group,co1)
# 
# plot(nmds.sq.sq, display="si", type="n", ylab=NA,xlab=NA,xaxt="n",yaxt="n")
# 
# for (i in unique(group$group)) {
# 	g1 <- nmds.sq.sq$points %>% as.data.frame() %>% filter(mm_coll_nmds$sample_type%in%group$sample_type[group$group==i])
# 	points(g1, display="si", col=group$co1[group$group==i], cex=1,pch=16)
# 	ordispider(g1, groups = rep(group$group[group$group==i] %>% unique(),nrow(g1)), col=group$co1[group$group==i],  lwd = 1, lty=4, label = T,cex=1)
# }
# dev.off()

# 
# # Figurexxx - Location (duplicates exposed) ------------------------------------------
# group <- mm_coll_nmds %>% distinct(Location) %>% mutate(group=Location)
# label <- group %>% pull(group)
# 
# # co <- c("tomato2","orange","springgreen2","#008d98","deepskyblue2")
# co <- c("deepskyblue2","tomato2","olivedrab4","orange","deepskyblue3","tomato")
# co1 <- co[factor(group$group)]
# group <- cbind(group,co1)
# 
# # jpeg(paste0("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/Plots/",
# # 						"Year difference all fjords.jpg"),
# # 		 width = 2000, height = 1900)
# 
# plot(nmds.sq.sq, display="si", type="n", ylab=NA,xlab=NA,xaxt="n",yaxt="n")
# 
# for (i in unique(group$group)) {
# 	g1 <- nmds.sq.sq$points %>% as.data.frame() %>% filter(mm_coll_nmds$Location%in%group$Location[group$group==i])
# 	points(g1, display="si", col=group$co1[group$group==i], cex=1,pch=16)
# 	ordispider(g1, groups = rep(group$group[group$group==i] %>% unique(),nrow(g1)), col=group$co1[group$group==i],  lwd = 1, lty=4, label = T,cex=1)
# }
# # par(mar=c(1,1,80,130))
# # boxplot(betadisper(dist.sq.sq,group = merr2(mm_coll_nmds$Location,group$Location,group$group),type = "centroid"),col=group$co1)
# legend("topright",inset=c(0.,0.02), legend = group %>% filter(!duplicated(group)) %>% pull(group), title="Location",
# 			 col = group$co1, pch = c(16), bty = "n", cex = 4)
# # dev.off()
# 
# # Figurexxx - Location (duplicates merged) ------------------------------------------
# group <- mm_coll_nmds %>% distinct(Location) %>% mutate(group=substr(Location,0,2))
# label <- group %>% pull(group)
# 
# # co <- c("tomato2","orange","springgreen2","#008d98","deepskyblue2")
# co <- c("deepskyblue2","orange","tomato2","olivedrab4")
# co1 <- co[factor(group$group)]
# group <- cbind(group,co1)
# 
# # jpeg(paste0("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/Plots/",
# # 						"Year difference all fjords.jpg"),
# # 		 width = 2000, height = 1900)
# 
# plot(nmds.sq.sq, display="si", type="n", ylab=NA,xlab=NA,xaxt="n",yaxt="n")
# 
# for (i in unique(group$group)) {
# 	g1 <- nmds.sq.sq$points %>% as.data.frame() %>% filter(mm_coll_nmds$Location%in%group$Location[group$group==i])
# 	points(g1, display="si", col=group$co1[group$group==i], cex=1,pch=16)
# 	ordispider(g1, groups = rep(group$group[group$group==i] %>% unique(),nrow(g1)), col=group$co1[group$group==i],  lwd = 1, lty=4, label = T,cex=1)
# }
# # par(mar=c(1,1,80,130))
# # boxplot(betadisper(dist.sq.sq,group = merr2(mm_coll_nmds$Location,group$Location,group$group),type = "centroid"),col=group$co1)
# legend("topright",inset=c(0.,0.02), legend = group %>% filter(!duplicated(group)) %>% pull(group), title="Location",
# 			 col = group$co1, pch = c(16), bty = "n", cex = 4)
# # dev.off()
# 
# 
# # Figurexxx - Location effect ------------------------------------------------------
# 
# group <- mm_coll_nmds %>% distinct(Location) %>% mutate(group=substr(Location,0,2)) %>% 
# 	mutate(group=if_else(group=="BA","BA","FR-OL-BE"))
# 
# co <- c("deepskyblue2","tomato2")
# co1 <- co[factor(group$group)]
# group <- cbind(group,co1)
# 
# # jpeg(paste0("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/Plots/",
# # 						"Year difference all fjords.jpg"),
# # 		 width = 2000, height = 1900)
# 
# plot(nmds.sq.sq, display="si", type="n", ylab=NA,xlab=NA,xaxt="n",yaxt="n")
# 
# for (i in unique(group$group)) {
# 	g1 <- nmds.sq.sq$points %>% as.data.frame() %>% filter(mm_coll_nmds$Location%in%group$Location[group$group==i])
# 	points(g1, display="si", col=group$co1[group$group==i], cex=1,pch=16)
# 	ordispider(g1, groups = rep(group$group[group$group==i] %>% unique(),nrow(g1)), col=group$co1[group$group==i],  lwd = 1, lty=4, label = T,cex=1)
# }
# # par(mar=c(1,1,80,130))
# # boxplot(betadisper(dist.sq.sq,group = merr2(mm_coll_nmds$Location,group$Location,group$group),type = "centroid"),col=group$co1)
# legend("topright",inset=c(0.,0.02), legend = group %>% filter(!duplicated(group)) %>% pull(group), title="Location",
# 			 col = group$co1, pch = c(16), bty = "n", cex = 4)
# # dev.off()
# 
# 
# 
# # Figurexxx Season - Location (BA) -------------------------------------------------
# group <- mm_coll_nmds %>% distinct(Year,Location,Season) %>% mutate(group=paste0(substr(Location,0,2),Season,Year)) %>% arrange(group) %>% filter(Location=="BA") %>% 
# 	mutate(label=Year)
# 
# co <- c("orange",
# 				"deepskyblue",
# 				"deepskyblue4")
# # co <- c("deepskyblue","deepskyblue2","deepskyblue3","deepskyblue4","dodgerblue4",
# # 				"orange","orange2","orange4",
# # 				"tomato","tomato2","tomato3",
# # 				"olivedrab1","olivedrab3","olivedrab4")
# co1 <- co[factor(group$Season)]
# group <- cbind(group,co1)
# 
# # jpeg(paste0("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/Plots/",
# # 						"Year difference all fjords.jpg"),
# # 		 width = 2000, height = 1900)
# 
# plot(nmds.sq.sq, display="si", type="n", ylab=NA,xlab=NA,xaxt="n",yaxt="n")
# 
# for (i in unique(group$group)) {
# 	g1 <- nmds.sq.sq$points %>% as.data.frame() %>% filter(mm_coll_nmds$Location%in%group$Location[group$group==i]&mm_coll_nmds$Season%in%group$Season[group$group==i]&mm_coll_nmds$Year%in%group$Year[group$group==i])
# 	points(g1, display="si", col=group$co1[group$group==i], cex=1,pch=16)
# 	ordispider(g1, groups = rep(group$label[group$group==i] %>% unique(),nrow(g1)), col=group$co1[group$group==i],  lwd = 1, lty=4, label = T,cex=1)
# }
# # par(mar=c(1,1,80,130))
# # boxplot(betadisper(dist.sq.sq,group = merr2(mm_coll_nmds$Location,group$Location,group$group),type = "centroid"),col=group$co1)
# legend("topright",inset=c(0.,0.02), legend = group %>% filter(!duplicated(Season)) %>% pull(Season), title="Location",
# 			 col = group %>% filter(!duplicated(Season)) %>% pull(co1), pch = c(16), bty = "n", cex = 4)
# # dev.off()
# 
# 
# # Figurexxx Season - Location (FR,OL,BE) -------------------------------------------
# 
# group <- mm_coll_nmds %>% distinct(Year,Location,Season) %>% mutate(group=paste0(substr(Location,0,2),Season,Year)) %>% arrange(group) %>% filter(Location!="BA"&Location!="BA_2") %>% 
# 	mutate(label=Year) %>% 
# 	mutate(Location2=substr(Location,0,2))
# 
# co <- c("orange",
# 				"deepskyblue",
# 				"deepskyblue4")
# co1 <- co[factor(group$Season)]
# sy <- c(15,17,18)
# sy1 <- sy[factor(substr(group$Location2,0,2))]
# group <- cbind(group,co1,sy1)
# 
# 
# 
# # jpeg(paste0("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/Plots/",
# # 						"Year difference all fjords.jpg"),
# # 		 width = 2000, height = 1900)
# 
# plot(nmds.sq.sq, display="si", type="n", ylab=NA,xlab=NA,xaxt="n",yaxt="n")
# 
# for (i in unique(group$group)) {
# 	g1 <- nmds.sq.sq$points %>% as.data.frame() %>% filter(mm_coll_nmds$Location%in%group$Location[group$group==i]&mm_coll_nmds$Season%in%group$Season[group$group==i]&mm_coll_nmds$Year%in%group$Year[group$group==i])
# 	points(g1, display="si", col=group$co1[group$group==i], cex=1,pch=group$sy[group$group==i])
# 	ordispider(g1, groups = rep(group$label[group$group==i] %>% unique(),nrow(g1)), col=group$co1[group$group==i],  lwd = 1, lty=4, label = T,cex=1)
# }
# # par(mar=c(1,1,80,130))
# # boxplot(betadisper(dist.sq.sq,group = merr2(mm_coll_nmds$Location,group$Location,group$group),type = "centroid"),col=group$co1)
# legend("topright",inset=c(0.,0.02), legend = group %>% filter(!duplicated(Season)) %>% pull(Season), title="Season",
# 			 col = group %>% filter(!duplicated(Season)) %>% pull(co1), pch = c(16), bty = "n", cex = 4)
# legend("bottomright",inset=c(0.,0.02), legend = group %>% filter(!duplicated(Location2)) %>% pull(Location2), title="Location",
# 			 col="black",pch = group %>% filter(!duplicated(Location2)) %>% pull(sy1), bty = "n", cex = 4)
# # dev.off()


# Figure 3 - nMDS Season - Location - Year -------------------------------------------

group <- mm_coll_nmds %>% distinct(Year,Location,Season) %>% mutate(group=paste0(substr(Location,0,2),Season,Year)) %>% arrange(group) %>% 
	mutate(label=Year) %>% 
	mutate(Location2=substr(Location,0,2)) %>% 
	mutate(Location2=if_else(Location2=="BA","Balfsfjord",Location2)) %>%
	mutate(Location2=if_else(Location2=="FR","Frakkfjord",Location2)) %>%
	mutate(Location2=if_else(Location2=="OL","Olderfjord",Location2)) %>%
	mutate(Location2=if_else(Location2=="BE","Bergsfjord",Location2))

co <- c("deepskyblue4",
				"orange",
				"deepskyblue")
co1 <- co[factor(group$Season)]
sy <- c(16,15,17,18)
sy1 <- sy[factor(substr(group$Location2,0,2))]
group <- cbind(group,co1,sy1)


par(mar=c(0.4,0.4,0.4,0.4))

# jpeg(paste0("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/Plots/",
# 						"Figure_3_nMDS_idx.jpg"),
# 		 width = 2000, height = 1900)

# pdf(paste0("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/Plots/",
# 						"Figure_3_nMDS.pdf"),
# 		width = 28, height = 26.5)

plot(nmds.sq.sq, display="si", type="n", ylab=NA,xlab=NA,xaxt="n",yaxt="n")

for (i in unique(group$group)) {
	g1 <- nmds.sq.sq$points %>% as.data.frame() %>% filter(mm_coll_nmds$Location%in%group$Location[group$group==i]&mm_coll_nmds$Season%in%group$Season[group$group==i]&mm_coll_nmds$Year%in%group$Year[group$group==i])
	points(g1, display="si", col=group$co1[group$group==i], cex=2,pch=group$sy[group$group==i])
}
for (i in unique(group$group)) {
	g1 <- nmds.sq.sq$points %>% as.data.frame() %>% filter(mm_coll_nmds$Location%in%group$Location[group$group==i]&mm_coll_nmds$Season%in%group$Season[group$group==i]&mm_coll_nmds$Year%in%group$Year[group$group==i])
	ordispider(g1, groups = rep(group$label[group$group==i] %>% unique(),nrow(g1)), col=group$co1[group$group==i],  lwd = 1.5, lty=1, label = F,cex=2)
	# points(x=g1 %>% pull(MDS1) %>% mean(),y=g1 %>% pull(MDS2) %>% mean(),
	# 			 col=group$co1[group$group==i], cex=3,pch=group$sy[group$group==i])
}
for (i in unique(group$group)) {
	g1 <- nmds.sq.sq$points %>% as.data.frame() %>% filter(mm_coll_nmds$Location%in%group$Location[group$group==i]&mm_coll_nmds$Season%in%group$Season[group$group==i]&mm_coll_nmds$Year%in%group$Year[group$group==i])
	rect(g1 %>% pull(MDS1) %>% mean() - 0.04, g1 %>% pull(MDS2) %>% mean() - 0.020,  
			 g1 %>% pull(MDS1) %>% mean() + 0.04, g1 %>% pull(MDS2) %>% mean() + 0.020,
			 col=group$co1[group$group==i])
	text(x=g1 %>% pull(MDS1) %>% mean(),y=g1 %>% pull(MDS2) %>% mean(),labels=group$Year[group$group==i],
				 col="white", cex=2,pch=group$sy[group$group==i])
}
# par(mar=c(1,1,80,130))
# boxplot(betadisper(dist.sq.sq,group = merr2(mm_coll_nmds$Location,group$Location,group$group),type = "centroid"),col=group$co1)
legend("topleft",inset=c(-0.04,-0.02), legend = paste0("Stress = ",round(nmds.sq.sq$stress,3)), bty = "n", cex = 4)
legend("topright",inset=c(0.,0.02), legend = group %>% filter(!duplicated(Season)) %>% pull(Season), title="Season",
			 col = group %>% filter(!duplicated(Season)) %>% pull(co1), pch = c(16), bty = "n", cex = 4)
legend("bottomright",inset=c(0.,0.02), legend = group %>% filter(!duplicated(Location2)) %>% pull(Location2), title="Location",
			 col="black",pch = group %>% filter(!duplicated(Location2)) %>% pull(sy1), bty = "n", cex = 4)
# dev.off()


par(mar=c(4,4,1,1))
accum_curve <- specaccum(t(otu_coll_nmds %>% select(mm_coll_nmds$collapse.name[mm_coll_nmds$Location=="BE"])))

# Plot the species accumulation curve
plot(accum_curve, xlab = "Number of stations visited", ylab = "Number of species observed",
		 main = "Species Accumulation Curve")




# Initialize an empty vector to store the cumulative number of species observed
accumulated_species <- numeric(ncol(otu_coll_nmds))

# Iterate through each station
for (i in 1:ncol(otu_coll_nmds)) {
	# Calculate the cumulative number of species observed at this station
	accumulated_species[i] <- length(unique(otu_coll_nmds[, i][otu_coll_nmds[, i] > 0]))
}

# Plot the species accumulation curve
plot(1:ncol(otu_coll_nmds), accumulated_species, type = "l", 
		 xlab = "Number of stations visited", ylab = "Number of species observed",
		 main = "Species Accumulation Curve")



# # Figurexxx Season - Location (All) -------------------------------------------
# mm_coll_nmds <- mm_coll_nmds %>% mutate(nmds_group=if_else(Location!="BA"&Location!="BA_2"&Season=="Autumn","1",NA)) %>%
# 	mutate(nmds_group=if_else(Season=="Spring","2",nmds_group)) %>%
# 	mutate(nmds_group=if_else((Location=="BA"|Location=="BA_2")&(Season=="Autumn"|Season=="Summer"),"3",nmds_group)) %>%
# 	mutate(nmds_group=if_else((Location!="BA"&Location!="BA_2")&(Season=="Summer")&Year==2020,"4",nmds_group)) %>%
# 	mutate(nmds_group=if_else((Location!="BA"&Location!="BA_2")&(Season=="Summer")&Year==2021,"5",nmds_group))
# 
# 
# group <- mm_coll_nmds %>% distinct(nmds_group,Year,Location,Season) %>% mutate(group=paste0(substr(Location,0,2),Season,Year)) %>% arrange(nmds_group) %>% 
# 	mutate(label=Year) %>% 
# 	mutate(Location2=substr(Location,0,2))
# 
# co <- c("orange",
# 				"deepskyblue",
# 				"tomato",
# 				"olivedrab",
# 				"deepskyblue4")
# co1 <- co[factor(group$nmds_group)]
# sy <- c(16,15,17,18)
# sy1 <- sy[factor(substr(group$Location2,0,2))]
# group <- cbind(group,co1,sy1)
# 
# # jpeg(paste0("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/Plots/",
# # 						"Main_Figure2.jpg"),
# # 		 width = 2000, height = 1900)
# 
# plot(nmds.sq.sq, display="si", type="n", ylab=NA,xlab=NA,xaxt="n",yaxt="n")
# 
# for (i in unique(group$group)) {
# 	g1 <- nmds.sq.sq$points %>% as.data.frame() %>% filter(mm_coll_nmds$Location%in%group$Location[group$group==i]&mm_coll_nmds$Season%in%group$Season[group$group==i]&mm_coll_nmds$Year%in%group$Year[group$group==i])
# 	points(g1, display="si", col=group$co1[group$group==i], cex=2,pch=group$sy[group$group==i])
# }
# for (i in unique(group$group)) {
# 	g1 <- nmds.sq.sq$points %>% as.data.frame() %>% filter(mm_coll_nmds$Location%in%group$Location[group$group==i]&mm_coll_nmds$Season%in%group$Season[group$group==i]&mm_coll_nmds$Year%in%group$Year[group$group==i])
# 	ordispider(g1, groups = rep(group$label[group$group==i] %>% unique(),nrow(g1)), col=group$co1[group$group==i],  lwd = 1.5, lty=1, label = F,cex=2)
# 	# points(x=g1 %>% pull(MDS1) %>% mean(),y=g1 %>% pull(MDS2) %>% mean(),
# 	# 			 col=group$co1[group$group==i], cex=3,pch=group$sy[group$group==i])
# }
# for (i in unique(group$group)) {
# 	g1 <- nmds.sq.sq$points %>% as.data.frame() %>% filter(mm_coll_nmds$Location%in%group$Location[group$group==i]&mm_coll_nmds$Season%in%group$Season[group$group==i]&mm_coll_nmds$Year%in%group$Year[group$group==i])
# 	rect(g1 %>% pull(MDS1) %>% mean() - 0.04, g1 %>% pull(MDS2) %>% mean() - 0.020,  
# 			 g1 %>% pull(MDS1) %>% mean() + 0.04, g1 %>% pull(MDS2) %>% mean() + 0.020,
# 			 col=group$co1[group$group==i])
# 	text(x=g1 %>% pull(MDS1) %>% mean(),y=g1 %>% pull(MDS2) %>% mean(),labels=group$Year[group$group==i],
# 				 col="white", cex=2,pch=group$sy[group$group==i])
# }
# # par(mar=c(1,1,80,130))
# # boxplot(betadisper(dist.sq.sq,group = merr2(mm_coll_nmds$Location,group$Location,group$group),type = "centroid"),col=group$co1)
# legend("topright",inset=c(0.,0.02), legend = group %>% filter(!duplicated(Season)) %>% pull(Season), title="Season",
# 			 col = group %>% filter(!duplicated(Season)) %>% pull(co1), pch = c(16), bty = "n", cex = 4)
# legend("bottomright",inset=c(0.,0.02), legend = group %>% filter(!duplicated(Location2)) %>% pull(Location2), title="Location",
# 			 col="black",pch = group %>% filter(!duplicated(Location2)) %>% pull(sy1), bty = "n", cex = 4)
# dev.off()



# Fugre 4 - Bubble plot -----------------------------------------------------------------------

# Are the otu column name in the same order an metadata?
identical(mm_coll_nmds$collapse.name,colnames(otu_coll_nmds))

tax_nmds_bub <- tax_nmds

# Collapse otu and mm into Location - Season - Year
otu_coll_nmds_bub <- otu_coll_nmds %>% `colnames<-`(paste0(mm_coll_nmds$Location,"_",mm_coll_nmds$Season,"_",mm_coll_nmds$Year))

otu_coll_nmds_bub <- collapse_otu2(((otu_coll_nmds_bub)),
																	 mm_coll_nmds %>% mutate(collapse.name=paste0(Location,"_",Season,"_",Year)),fun = "Mean")
otu_coll_nmds_bub <- ((otu_coll_nmds_bub))

mm_coll_nmds_bub <- collapse_meta(mm_coll_nmds %>% mutate(collapse.name=paste0(Location,"_",Season,"_",Year)))


tax_nmds_bub$species[is.na(tax_nmds_bub$species)] <- "-"
tax_nmds_bub$species <- make.unique.2(tax_nmds_bub$species)

totu_bub <- trans(otu_coll_nmds_bub) %>% setNames(c("Species","Sample","Reads"))
totu_bub$Species <- merr2(totu_bub$Species,rownames(tax_nmds_bub),tax_nmds_bub$species)

totu_bub$Location <- merr2(totu_bub$Sample,mm_coll_nmds_bub$group,as.character(mm_coll_nmds_bub$Location))
totu_bub$Season <- merr2(totu_bub$Sample,mm_coll_nmds_bub$group,as.character(mm_coll_nmds_bub$Season))
totu_bub$Year <- merr2(totu_bub$Sample,mm_coll_nmds_bub$group,mm_coll_nmds_bub$Year)

x_axis <- totu_bub %>% distinct(Season,Year) %>% arrange(desc(Season),desc(Year))
lo <- levels(factor(totu$Location, levels = c("BA", "FR", "OL", "BE", "BA_2", "FR_2")))
sp_l <- c("Gadus morhua",
					"Mallotus villosus",
					"Salmo salar",
					"Clupea harengus",
					"Melanogrammus aeglefinus|Merlangius merlangus",
					"Pleuronectes platessa",
					"Hippoglossoides platessoides",
					"Pollachius virens",
					"Anarhichas spp.",
					"Leptoclinus maculatus")


co <- c("deepskyblue",
				"orange",
				"tomato",
				"deepskyblue3",
				"darkorange2",
				"olivedrab",
				"cornsilk3",
				"thistle2",
				"paleturquoise1",
				"deepskyblue4")


jpeg(paste0("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/Plots/",
						"Figure_4_Bubble_plot.jpg"),
		 width = 1500, height = 1000)
pdf(paste0("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/Plots/",
						"Figure_4_Bubble_plot.pdf"),
		 width = 15*1.4, height = 10*1.4)
{par(mar = c(4, 10, 4, 1))
	plot(c(1, 4.8), c(1, nrow(x_axis)), xlab = "", ylab = "", type = "n", las = 1, axes = FALSE, bty = "o")
	axis(side = 2, at = c(1:nrow(x_axis)), labels = paste0(x_axis$Season, " ", x_axis$Year), las = 2, tick = FALSE,cex.axis = 1.6)  
	axis(side = 1, at = c(1:4) + 0.4, labels = c("Balsfjord", "Frakkfjord", "Olderfjord", "Bergsfjord"), tick = FALSE, cex.axis = 2.0)
	abline(h=5.5);abline(h=3.5)
	abline(v=0.92);abline(v=1.92);abline(v=2.92);abline(v=3.92)

	for (i in 1:4) {
		for (k in 1:(length(sp_l))) {
			loc <- lo[i]
			
			bub <- totu_bub %>% 
				filter(Species==sp_l[k]) %>% 
				filter(Location==loc)
			
			bub <- x_axis %>% left_join(.,bub,by=c("Season"="Season","Year"="Year"))
			
			##Circle values
			for (j in 1:nrow(bub)){
				par(new=T)
				points(i+((k-0.5)*0.09),j, pch=21, cex=0.5*sqrt(sqrt(as.numeric(bub$Reads[j]))),bg=co[k],xpd=NA)
			}
		}
	}

# Apply italic formatting to each element of sp_l
italic_species <- lapply(sp_l, function(x) substitute(italic(y), list(y = x)))

par(new=T)  # Create a new plot
par(mar = c(0, 0, 0, 0))  # Set margins to zero for the legend plot

# Plot legend with italicized species names
legend("topright", inset = c(0, -0.05), legend = italic_species,
			 title = "Species",
			 col = co, pch = 19, cex = 1.8,
			 bg = "white", box.lwd = 1, box.lty = 1, box.col = "black")
}

{par(new=T)  # Create a new plot
par(mar = c(0, 0, 0, 0))  # Set margins to zero for the legend plot

rect(2.2, 7.8,  
		 3.3, 10.85,
		 col="white",border="black")

points(2.4,10.3, pch=21, cex=sqrt(sqrt(1))*0.5, bg="white",xpd=NA)
points(2.4,10.1, pch=21, cex=sqrt(sqrt(100))*0.5, bg="white",xpd=NA)
points(2.4,9.8, pch=21, cex=sqrt(sqrt(10000))*0.5, bg="white",xpd=NA)
points(2.4,9.3, pch=21, cex=sqrt(sqrt(100000))*0.5, bg="white",xpd=NA)
points(2.4,8.5, pch=21, cex=sqrt(sqrt(1000000))*0.5, bg="white",xpd=NA)

text(x=2.6,y=c(10.3,10.1,9.8,9.3,8.5),cex=1.5,
		 labels=c("1","100","10000","100000","1000000"),xpd=NA,pos = 4,adj=c(1,0.5),font=1)
text(x=2.2,y=c(10.6),cex=1.8,
		 labels=c("Reads (mean per station)"),xpd=NA,pos = 4,adj=c(1,0.5),font=1)
}

dev.off()


# # Figurexxx - Plate - Location ---------------------------------------------------------
# 
# group <- mm_coll_nmds %>% distinct(Plate,Location) %>% mutate(group=paste0(substr(Location,0,2),Plate)) %>% arrange(group) %>% 
# 	mutate(label=substr(Plate,7,nchar(.)))
# co <- c("deepskyblue",
# 				"orange",
# 				"tomato",
# 				"olivedrab")
# co1 <- co[factor(substr(group$Location,0,2))]
# group <- cbind(group,co1)
# group
# # jpeg(paste0("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/Plots/",
# # 						"Year difference all fjords.jpg"),
# # 		 width = 2000, height = 1900)
# 
# plot(nmds.sq.sq, display="si", type="n", ylab=NA,xlab=NA,xaxt="n",yaxt="n")
# 
# for (i in unique(group$group)) {
# 	g1 <- nmds.sq.sq$points %>% as.data.frame() %>% filter(mm_coll_nmds$Location%in%group$Location[group$group==i]&mm_coll_nmds$Plate%in%group$Plate[group$group==i])
# 	points(g1, display="si", col=group$co1[group$group==i], cex=1,pch=16)
# 	ordispider(g1, groups = rep(group$label[group$group==i] %>% unique(),nrow(g1)), col=group$co1[group$group==i],  lwd = 1, lty=4, label = T,cex=0.7)
# }
# # par(mar=c(1,1,80,130))
# # boxplot(betadisper(dist.sq.sq,group = merr2(mm_coll_nmds$Location,group$Location,group$group),type = "centroid"),col=group$co1)
# legend("topright",inset=c(0.,0.02), legend = group %>% filter(!duplicated(group)) %>% pull(group), title="Location",
# 			 col = group %>% filter(!duplicated(group)) %>% pull(co1), pch = c(16), bty = "n", cex = 4)
# # dev.off()
# 
# sp.vec <- envfit(nmds.sq.sq, t(sqrt(sqrt(otu_coll_nmds)))) # species scores
# plot(sp.vec,p.max = 0.001,col="mediumspringgreen",cex=0.5)




# Figure 5 - Heatmap ----------------------------------------------------------------
mm_coll_heat <- mm_coll_nmds %>% 
	mutate(nmds_group=if_else(Location!="BA"&Location!="BA_2"&Season=="Autumn","1",NA)) %>%
	mutate(nmds_group=if_else(Season=="Spring","2",nmds_group)) %>%
	mutate(nmds_group=if_else((Location=="BA"|Location=="BA_2")&(Season=="Autumn"|Season=="Summer"),"3",nmds_group)) %>%
	mutate(nmds_group=if_else((Location!="BA"&Location!="BA_2")&(Season=="Summer")&Year==2020,"4",nmds_group)) %>%
	mutate(nmds_group=if_else((Location!="BA"&Location!="BA_2")&(Season=="Summer")&Year==2021,"5",nmds_group))

mm_coll_heat <- mm_coll_heat %>% arrange(nmds_group,Location,Season,Year)
otu_coll_heat <- otu_coll_nmds %>% 
	select(mm_coll_heat %>% pull(collapse.name)) %>% 
	filter(!is.na(tax_nmds$species))
tax_heat <- tax_nmds[!is.na(tax_nmds$species),]

jpeg(paste0("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/Plots/",
						"Figure_5_Heatmap.jpg"),
		 width = 1550, height = 1200)
pdf(paste0("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/Plots/",
						"Figure_5_Heatmap.pdf"),
		 width = 15.5*1.5, height = 12*1.5)
heatmap(sqrt(sqrt(otu_coll_heat %>% filter(rowSums(.)>1000))) %>% as.matrix(), cexRow=1.7,cexCol=1.4,Rowv=NA, Colv=NA, col = hcl.colors(13,palette = "Blues",rev=T), scale="none",labCol = paste0(mm_coll_heat$Season,"_",mm_coll_heat$Location),
				margins=c(8,30),labRow=as.expression(lapply((tax_heat$species), function(a) bquote(italic(.(a))))))

legend(x="left",title=expression(sqrt(sqrt("eDNA reads"))), legend=round(seq(0,sqrt(sqrt(max(otu_coll_heat))),length.out = 13),0),
			 fill=hcl.colors(13,palette = "Blues",rev=T),cex = 1.4)
dev.off()


#

d <- otu_coll_nmds_loseye[!is.na(tax_nmds$species)&tax_nmds$species=="Salmo salar",]


indsp <- multipatt(t(sqrt(sqrt(otu_coll_nmds))),mm_coll_nmds$nmds_group,control = how(nperm=999),func="IndVal.g")
indsp$sign %>% mutate(id=row_number()) %>% filter(p.value<0.01)
summary(indsp) 

heat_mm <- mm_coll_nmds %>% arrange(nmds_group,Location,Season,Year)
heat_otu <- otu_coll_nmds %>% select(heat_mm$collapse.name)

jpeg(paste0("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/Plots/",
						"Main_Heatmap.jpg"),
		 width = 1500, height = 1150)
heatmap(sqrt(sqrt(heat_otu)) %>% as.matrix(), cexRow=1.2,cexCol=1.1,Rowv=NA, Colv=NA, col = hcl.colors(13,palette = "Blues",rev=T), scale="none",labCol = paste0(heat_mm$Location,"_",heat_mm$Season,"_",heat_mm$Year),
				margins=c(10,15),labRow=as.expression(lapply(rownames(heat_otu), function(a) bquote(italic(.(a))))))
legend(x="left",title=expression(sqrt(sqrt("eDNA reads"))), legend=round(seq(0,sqrt(sqrt(max(heat_otu))),length.out = 13),0),
			 fill=hcl.colors(13,palette = "Blues",rev=T))
dev.off()



# Figure 6 - Species accumulation curves --------------------------------------------

otu_coll_accum <- otu_coll
mm_coll_accum <- mm_coll

# mm_coll_accum <- mm_coll_accum %>% filter(Season=="Spring"&Year=="2021")
# otu_coll_accum <- otu_coll_accum %>% select(mm_coll_accum$collapse.name)

iter=100

list_ba <- which(mm_coll_accum$Location=="BA")
list_fr <- which(mm_coll_accum$Location=="FR")
list_ol <- which(mm_coll_accum$Location=="OL")
list_be <- which(mm_coll_accum$Location=="BE")

# list <- list(list_ba,list_fr,list_ol,list_be)
list <- list(list_be,list_fr,list_ol)
li_le <- length(list)
x_columns <- sapply(list, length) %>% sum()
x <- matrix(NA,ncol = x_columns,nrow = iter)

# co <- rep(c("tomato2","orange2","deepskyblue2","green3"),sapply(list,length))
c1 <- rep(c("Olderfjord","Frakkfjord","Bergsfjord"),sapply(list,length))
co <- rep(c("tomato2","orange2","green3"),sapply(list,length))

for (j in 1:iter) {
	
	for (i in 1:x_columns) {
		
		if (li_le>=1) {
			if (i < sapply(list, length)[1] %>% sum()+1) {
				selected_samples <- sample(list[[1]],i)
			} }
		
		if (li_le>=2) {
			if (i > sapply(list, length)[1] %>% sum() & i < sapply(list, length)[1:2] %>% sum()+1) {
				selected_samples <- sample(list[[1]],length(list[[1]]))
				selected_samples <- c(selected_samples,sample(list[[2]],i-length(list[[1]])))
			} }
		
		if (li_le>=3) {
			if (i > sapply(list, length)[1:2] %>% sum() & i < sapply(list, length)[1:3] %>% sum()+1) {
				selected_samples <- sample(list[[1]],length(list[[1]]))
				selected_samples <- c(selected_samples,sample(list[[2]],length(list[[2]])))
				selected_samples <- c(selected_samples,sample(list[[3]],i-sapply(list, length)[1:2] %>% sum()))
			} }
		
		if (li_le>=4) {
			if (i > sapply(list, length)[1:3] %>% sum() & i < sapply(list, length)[1:4] %>% sum()+1) {
				selected_samples <- sample(list[[1]],length(list[[1]]))
				selected_samples <- c(selected_samples,sample(list[[2]],length(list[[2]])))
				selected_samples <- c(selected_samples,sample(list[[3]],length(list[[3]])))
				selected_samples <- c(selected_samples,sample(list[[4]],i-sapply(list, length)[1:3] %>% sum()))
			} }
		
		x[j,i] <- otu_coll_accum[,selected_samples] %>% as.data.frame() %>% rowSums() %>% as.data.frame() %>% 
			pa_conversion() %>% colSums() %>% setNames("n_sp")
	}
	
}

acc_cur <- matrix(NA,ncol = ncol(x),nrow=3)
acc_cur[1,] <- colMeans(x)
acc_cur[2,] <- apply(x, 2, quantile, probs = 0.25)
acc_cur[3,] <- apply(x, 2, quantile, probs = 0.75)
acc_cur <- t(acc_cur) %>% as.data.frame() %>% rownames_to_column("Sample") %>% 
	cbind(.,c1 %>% as.data.frame()) %>% 
	setNames(c("Sample","Mean","q25","q75","Location")) %>% 
	mutate(Sample=as.numeric(Sample))

ggplot()+
	geom_errorbar(data=acc_cur,aes(x=Sample,ymin=q25,ymax=q75,color=Location),width=0.8)+
	geom_point(data=acc_cur,aes(x=Sample,y=Mean,color=Location),size=2)+
	theme_bw()



#
#
# 
# 
# 
# 

boxplot(betadisper(dist.sq,group = mm_coll_nmds$Year,type = "centroid"),ylim=c(0,1))
boxplot(betadisper(dist.sq,group = paste0(mm_coll_nmds$Year,mm_coll_nmds$Season,type = "centroid")),ylim=c(0,1))
boxplot(betadisper(dist.sq,group = paste0(mm_coll_nmds$Year,mm_coll_nmds$Location,type = "centroid")),ylim=c(0,1))
boxplot(betadisper(dist.sq,group = paste0(mm_coll_nmds$Year,mm_coll_nmds$Season,mm_coll_nmds$Location,type = "centroid")),ylim=c(0,1))

x <- betadisper(dist.rel,group = mm_coll_nmds$sample_type,type = "centroid")
x1 <- betadisper(dist.rel,group = mm_coll_nmds$Year,type = "centroid")
x2 <- betadisper(dist.rel,group = paste0(mm_coll_nmds$Year,mm_coll_nmds$Season,type = "centroid"))
x3 <- betadisper(dist.rel,group = paste0(mm_coll_nmds$Year,mm_coll_nmds$Location,type = "centroid"))
x4 <- betadisper(dist.rel,group = paste0(mm_coll_nmds$Year,mm_coll_nmds$Season,mm_coll_nmds$Location,type = "centroid"))
x5 <- betadisper(dist.rel,group = paste0(mm_coll_nmds$Year,mm_coll_nmds$Season,mm_coll_nmds$Station_name,type = "centroid"))
x6 <- betadisper(dist.rel,group = paste0(mm_coll_nmds$Year,mm_coll_nmds$Plate,type = "centroid"))



permdisp <- rbind(x$distances %>% as.data.frame() %>% mutate(group="1.All"),
									x1$distances%>% as.data.frame() %>% mutate(group="2.Year"),
									x2$distances%>% as.data.frame()%>% mutate(group="3.Year-Season"),
									x3$distances%>% as.data.frame()%>% mutate(group="4.Year-Location"),
									x4$distances%>% as.data.frame()%>% mutate(group="5.Year-Season-Location"),
									# x6$distances%>% as.data.frame()%>% mutate(group="7.Year-Plate"),
									x5$distances%>% as.data.frame()%>% mutate(group="6.Year-Season-Location-Station")) %>% setNames(c("distance","Grouping"))


# p <- 
ggplot() +
	geom_density(data=permdisp,aes(x=distance,fill=Grouping),colour="black")+
	theme_bw()+
	xlim(0,1)+
	scale_fill_manual(values = c("black","deepskyblue2","orange","olivedrab4","tomato2","mediumblue","violet"))
ggsave(paste0("/Users/a36142/Library/CloudStorage/OneDrive-Havforskningsinstituttet/IMR/Data/eDNA/Paper_III/Plots/",
							"Centroid distances.jpg"), 
			 plot=p, device="jpg", width=10, height=10, units="in", dpi=300)