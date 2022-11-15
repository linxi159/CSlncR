################################################################################################
############## Running scripts for exploring cell-specific lncRNA regulation ##############
################################################################################################

## Load required R packages, please firstly install the following R packages before running scripts
library(pracma)
library(WGCNA)
library(igraph)
library(miRspongeR)
library(biclique)
library(corrplot)
library(dendextend)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(clusterProfiler)
library(msigdbr)
library(vroom)
library(doParallel)

## Load utility functions
source("CSlncR.R")

##############################################################################################################################
############################################## <<1>> Loading data  ###########################################################
##############################################################################################################################
## Load prepared datasets
ASD <- read.csv("./data/ASD_gene_lncRNAs_mRNAs.csv")
ASD <- as.matrix(ASD)

## Load the preprocessing single-cell sequencing data including filtering genes with expression <= 10 cells, and log e(x+1).
lncRNA_scRNA_raw  <- read.csv("./data/GSE71315/GSE71315_scell_ncounts_genes_thresh_lncRNA.csv")
rownames(lncRNA_scRNA_raw ) <- lncRNA_scRNA_raw[,1]
lncRNA_scRNA_raw  <- lncRNA_scRNA_raw[,-1]
lncRNA_scRNA_raw   <- t(lncRNA_scRNA_raw)
lncRNA_scRNA_raw   <- as.data.frame(lncRNA_scRNA_raw)

mRNA_scRNA_raw <- read.csv("./data/GSE71315/GSE71315_scell_ncounts_genes_thresh_mRNA.csv")
rownames(mRNA_scRNA_raw) <- mRNA_scRNA_raw[,1]
mRNA_scRNA_raw <- mRNA_scRNA_raw[,-1]
mRNA_scRNA_raw  <- t(mRNA_scRNA_raw)
mRNA_scRNA_raw <- as.data.frame(mRNA_scRNA_raw)

## compute the average expression values of duplicate genes and remove genes with constant expression values in all cells
# Transformation using log2(x+1)
lncRNA_scRNA_norm <- log2(lncRNA_scRNA_raw+1)
mRNA_scRNA_norm <- log2(mRNA_scRNA_raw+1) 

# Compute the average expression values of duplicate genes
lncRNA_scRNA_norm_average <- Averg_Duplicate(lncRNA_scRNA_norm)
mRNA_scRNA_norm_average <- Averg_Duplicate(mRNA_scRNA_norm)

# Remove genes with constant expression values in all cells
lncRNA_scRNA_norm_sd <- unlist(lapply(seq(dim(lncRNA_scRNA_norm_average)[2]), function(i) sd(lncRNA_scRNA_norm_average[, i])))
lncRNA_scRNA_norm_filter <- lncRNA_scRNA_norm_average[, which(lncRNA_scRNA_norm_sd > 0)]
mRNA_scRNA_norm_sd <- unlist(lapply(seq(dim(mRNA_scRNA_norm_average)[2]), function(i) sd(mRNA_scRNA_norm_average[, i])))
mRNA_scRNA_norm_filter <- mRNA_scRNA_norm_average[, which(mRNA_scRNA_norm_sd > 0)]

# loading prediction lcnRNA-target to filter to human brain lncRNA-mRNA 
lncRNA_mRNA_pre_1 <- unique(read.csv("./data/lncRNA_mRNA-pre.csv"))
lncRNA_mRNA_pre_lncRNA <- as.vector(as.matrix(lncRNA_mRNA_pre_1[1]))
lncRNA_name <- colnames(lncRNA_scRNA_norm_filter)
lncRNA_mRNA_pre_2 <- lncRNA_mRNA_pre_1[which(lncRNA_mRNA_pre_lncRNA %in% lncRNA_name),]
lncRNA_mRNA_pre_mRNA <- as.vector(as.matrix(lncRNA_mRNA_pre_2[2]))
mRNA_name <- colnames(mRNA_scRNA_norm_filter)
lncRNA_mRNA_pre_3 <- lncRNA_mRNA_pre_2[which(lncRNA_mRNA_pre_mRNA %in% mRNA_name),]
lncRNA_mRNA_pre <- lncRNA_mRNA_pre_3 
# loading experimental lcnRNA-target
lncRNA_mRNA_exp <- read.csv("./data/lncRNA_mRNA-exp.csv")
lncRNA_mRNA_exp <- as.matrix(lncRNA_mRNA_exp)
lncRNA_mRNA_exp <- unique(lncRNA_mRNA_exp)

# save    
# save(lncRNA_scRNA_norm, mRNA_scRNA_norm, lncRNA_scRNA_norm_filter,mRNA_scRNA_norm_filter,lncRNA_mRNA_pre,lncRNA_mRNA_exp,file = "Exp_247_lncRNAs_10208_mRNAs_276_single_cells_GSE71315.RData")
load("Exp_247_lncRNAs_10208_mRNAs_276_single_cells_GSE71315.RData")


##############################################################################################################################
####################################### <<2>> Exploring cell-specific lncRNA regulation ######################################
##############################################################################################################################
## Discovering cell-specific lncRNA-mRNA regulatory network   
CSlncR_network_null <- CSlncR_net(lncRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter)
load("CSlncR_network_null_GSE71315.RData")
CSlncR_network_null = CSlncR_network_null[[1]]
prior_graph <- make_graph(c(t(lncRNA_mRNA_pre[, 1:2])), directed = TRUE)
CSlncR_network_null_graph <- lapply(seq(CSlncR_network_null), function(i) make_graph(c(t(CSlncR_network_null[[i]][, 1:2])), directed = TRUE))
CSlncR_network <- lapply(seq(CSlncR_network_null), function(i) as_data_frame(CSlncR_network_null_graph[[i]] %s% prior_graph))
# without prior knowledge
CSlncR_network_ <- lapply(seq(CSlncR_network_null), function(i) as_data_frame(CSlncR_network_null_graph[[i]]))

# save    
# save(CSlncR_network_null_graph,file = "CSlncR_network_null_graph_GSE71315.RData")
# save(CSlncR_network,file = "CSlncR_network_GSE71315.RData")
# save(CSlncR_network_,file = "CSlncR_network_GSE71315_.RData")
# load("CSlncR_network_null_graph_GSE71315.RData")       
load("CSlncR_network_GSE71315.RData")  # performing experiments to obtain CSlncR_network with prior knowledge
load("CSlncR_network_GSE71315_.RData") # performing experiments to obtain CSlncR_network_ without prior knowledge


##############################################################################################################################
#################### <<3>> single-cell lncRNA regulation network for each cell type (development stages) #####################
##############################################################################################################################
# GW16  GW21     A1      A2       S44      S46
# GW16  GW21     GW20.5  GW20.5   GW19.5   GW23.5
# 1:26  27:50    51:115  116:173  174:199  200:276
#  26    24            123           26       77
# [0]GW16 to GW23.5
# Merging single-cell networks from 5 different developmental stages, separately, (CSlncR_network)
CSlncR_network_src <- CSlncR_network
CSlncR_network_GW16 <- CSlncR_network_src[c(1:26)]
CSlncR_network_GW19.5 <- CSlncR_network_src[c(174:199)]
CSlncR_network_GW20.5 <- CSlncR_network_src[c(51:173)]
CSlncR_network_GW21 <- CSlncR_network_src[c(27:50)]
CSlncR_network_GW23.5 <- CSlncR_network_src[c(200:276)]

# Adjusted GW16 to GW23.5
tmp_1 <- append(CSlncR_network_GW16,CSlncR_network_GW19.5)
tmp_1 <- append(tmp_1,CSlncR_network_GW20.5)
tmp_1 <- append(tmp_1,CSlncR_network_GW21)
CSlncR_network_GW16_to_GW23.5 <-  append(tmp_1,CSlncR_network_GW23.5)

# Integrate expression values in the order GW16 to GW23.5
lncRNA_scRNA_norm_filter_GW16_to_GW23.5 <- rbind(lncRNA_scRNA_norm_filter[1:26,],lncRNA_scRNA_norm_filter[174:199,],lncRNA_scRNA_norm_filter[51:173,],lncRNA_scRNA_norm_filter[27:50,],lncRNA_scRNA_norm_filter[200:276,])
mRNA_scRNA_norm_filter_GW16_to_GW23.5 <- rbind(mRNA_scRNA_norm_filter[1:26,],mRNA_scRNA_norm_filter[174:199,],mRNA_scRNA_norm_filter[51:173,],mRNA_scRNA_norm_filter[27:50,],mRNA_scRNA_norm_filter[200:276,])
lncRNA_mRNA_scRNA_norm_filter_GW16_to_GW23.5 <- cbind(lncRNA_scRNA_norm_filter_GW16_to_GW23.5,mRNA_scRNA_norm_filter_GW16_to_GW23.5)

#rm(list= "CSlncR_network_src","CSlncR_network_null_graph","CSlncR_network_null")

# [1-2-3-4-5:GW16,GW21,GW20.5,GW19.5,GW23.5]: selection option [1-5]
#[1] GW16
CSlncR_network_GW16 <- CSlncR_network[c(1:26)]
lncRNA_scRNA_norm_filter <- lncRNA_scRNA_norm_filter[1:26,]
mRNA_scRNA_norm_filter <- mRNA_scRNA_norm_filter[1:26,]
CSlncR_network  <- CSlncR_network_GW16
cell_number <- 26
# without prior knowledge
CSlncR_network_GW16_ <- CSlncR_network_[c(1:26)]
CSlncR_network_  <- CSlncR_network_GW16_
lncRNA_mRNA_scRNA_norm_filter <- cbind(lncRNA_scRNA_norm_filter,mRNA_scRNA_norm_filter)

# #[2] GW21
# CSlncR_network_GW21 <- CSlncR_network[c(27:50)]
# lncRNA_scRNA_norm_filter <- lncRNA_scRNA_norm_filter[27:50,]
# mRNA_scRNA_norm_filter <- mRNA_scRNA_norm_filter[27:50,]
# CSlncR_network  <- CSlncR_network_GW21
# cell_number <- 24
# # without prior knowledge
# CSlncR_network_GW21_ <- CSlncR_network_[c(27:50)]
# CSlncR_network_  <- CSlncR_network_GW21_
# lncRNA_mRNA_scRNA_norm_filter <- cbind(lncRNA_scRNA_norm_filter,mRNA_scRNA_norm_filter)

# #[3] GW20.5
# CSlncR_network_GW20.5 <- CSlncR_network[c(51:173)]
# lncRNA_scRNA_norm_filter <- lncRNA_scRNA_norm_filter[51:173,]
# mRNA_scRNA_norm_filter <- mRNA_scRNA_norm_filter[51:173,]
# CSlncR_network  <- CSlncR_network_GW20.5
# cell_number <- 123
# # without prior knowledge
# CSlncR_network_GW20.5_ <- CSlncR_network_[c(51:173)]
# CSlncR_network_  <- CSlncR_network_GW20.5_
# lncRNA_mRNA_scRNA_norm_filter <- cbind(lncRNA_scRNA_norm_filter,mRNA_scRNA_norm_filter)

# #[4] GW19.5
# CSlncR_network_GW19.5 <- CSlncR_network[c(174:199)]
# lncRNA_scRNA_norm_filter <- lncRNA_scRNA_norm_filter[174:199,]
# mRNA_scRNA_norm_filter <- mRNA_scRNA_norm_filter[174:199,]
# CSlncR_network  <-  CSlncR_network_GW19.5
# cell_number <- 26
# # without prior knowledge
# CSlncR_network_GW19.5_ <- CSlncR_network_[c(174:199)]
# CSlncR_network_  <- CSlncR_network_GW19.5_
# lncRNA_mRNA_scRNA_norm_filter <- cbind(lncRNA_scRNA_norm_filter,mRNA_scRNA_norm_filter)

# #[5] GW23.5
# CSlncR_network_GW23.5 <- CSlncR_network[c(200:276)]
# lncRNA_scRNA_norm_filter <- lncRNA_scRNA_norm_filter[200:276,]
# mRNA_scRNA_norm_filter <- mRNA_scRNA_norm_filter[200:276,]
# CSlncR_network  <- CSlncR_network_GW23.5
# cell_number <- 77
# # without prior knowledge
# CSlncR_network_GW23.5_ <- CSlncR_network_[c(200:276)]
# CSlncR_network_  <- CSlncR_network_GW23.5_
# lncRNA_mRNA_scRNA_norm_filter <- cbind(lncRNA_scRNA_norm_filter,mRNA_scRNA_norm_filter)

##############################################################################################################################
########################################## <<4>> Generating comparison data  #################################################
##############################################################################################################################
############ CSlncR_network_random ###########
# random method as following # Extraction without replacement
# The implementation is as following.

############ CSlncR_network_LncRNA2Target ###########
# LncRNA2Target
CSlncR_network_LncRNA2Target_1 <- unique(read.csv("./data/lncRNA_mRNA-LncRNA2TargetPrediction.csv"))
LncRNA2Target_lncRNA <- as.vector(as.matrix(CSlncR_network_LncRNA2Target_1[1]))
lncRNA_name <- colnames(lncRNA_scRNA_norm_filter)
CSlncR_network_LncRNA2Target_2 <- CSlncR_network_LncRNA2Target_1[which(LncRNA2Target_lncRNA %in% lncRNA_name),]
LncRNA2Target_mRNA <- as.vector(as.matrix(CSlncR_network_LncRNA2Target_2[2]))
mRNA_name <- colnames(mRNA_scRNA_norm_filter)
CSlncR_network_LncRNA2Target_3 <- CSlncR_network_LncRNA2Target_2[which(LncRNA2Target_mRNA %in% mRNA_name),]
CSlncR_network_LncRNA2Target_ <- CSlncR_network_LncRNA2Target_3 
CSlncR_network_LncRNA2Target <-list()
for (i in 1:cell_number){
  CSlncR_network_LncRNA2Target[[i]] <- CSlncR_network_LncRNA2Target_
}
# save(CSlncR_network_LncRNA2Target,file = "CSlncR_network_LncRNA2Target.RData")
# load("CSlncR_network_LncRNA2Target.RData")

############ CSlncR_network_NPInterPrediction  ###########
# NPInterPrediction
lncRNA_gene_NPInterPrediction <- read.csv("./data/lncRNA_mRNA-NPInterPrediction.csv")
CSlncR_network_NPInterPrediction_1 <-unique(lncRNA_gene_NPInterPrediction)
NPInterPrediction_lncRNA <- as.vector(as.matrix(CSlncR_network_NPInterPrediction_1[1]))
lncRNA_name <- colnames(lncRNA_scRNA_norm_filter)
CSlncR_network_NPInterPrediction_2 <- CSlncR_network_NPInterPrediction_1[which(NPInterPrediction_lncRNA %in% lncRNA_name),]
NPInterPrediction_mRNA <- as.vector(as.matrix(CSlncR_network_NPInterPrediction_2[2]))
mRNA_name <- colnames(mRNA_scRNA_norm_filter)
CSlncR_network_NPInterPrediction_3 <- CSlncR_network_NPInterPrediction_2[which(NPInterPrediction_mRNA %in% mRNA_name),]
CSlncR_network_NPInterPrediction_ <- CSlncR_network_NPInterPrediction_3 
CSlncR_network_NPInterPrediction <-list()
for (i in 1:cell_number){
  CSlncR_network_NPInterPrediction[[i]] <- CSlncR_network_NPInterPrediction_1
}
# save(CSlncR_network_NPInterPrediction,file = "CSlncR_network_NPInterPrediction.RData")
# load("CSlncR_network_NPInterPrediction.RData")


###########################################################################################################################
################################## <<5>> 15 downstream computation  ->  5 Result analysis #################################
###########################################################################################################################
# downstream [1] validated
## Experimentally validated cell-specific lncRNA-mRNA interactions, the ground-truth (lncRTarget variable) is from the literature of NPInter database.
# remove missing value
    lncRNA_mRNA_exp <- lncRNA_mRNA_exp[!is.na(lncRNA_mRNA_exp[1:length(lncRNA_mRNA_exp[,1]),1]),]  
    lncRTarget_graph <- make_graph(c(t(lncRNA_mRNA_exp[, 1:2])), directed = TRUE)
    CSlncR_network_graph <- lapply(seq(CSlncR_network), function(i) make_graph(c(t(CSlncR_network[[i]][, 1:2])), directed = TRUE))
    CSlncR_network_validated <- lapply(seq(CSlncR_network), function(i) as_data_frame(CSlncR_network_graph[[i]] %s% lncRTarget_graph))
    # save(CSlncR_network_validated,file = "CSlncR_network_validated.RData")
    # load("CSlncR_network_validated.RData") 
    
    # GW16_to_GW23.5
    CSlncR_network_GW16_to_GW23.5_graph <- lapply(seq(CSlncR_network_GW16_to_GW23.5), function(i) make_graph(c(t(CSlncR_network_GW16_to_GW23.5[[i]][, 1:2])), directed = TRUE))
    
    # without prior knowledge 
    CSlncR_network_graph_ <- lapply(seq(CSlncR_network_), function(i) make_graph(c(t(CSlncR_network_[[i]][, 1:2])), directed = TRUE))
    CSlncR_network_validated_ <- lapply(seq(CSlncR_network_), function(i) as_data_frame(CSlncR_network_graph_[[i]] %s% lncRTarget_graph))
    # save(CSlncR_network_validated_,file = "CSlncR_network_validated_.RData")
    # load("CSlncR_network_validated_.RData")  

    # random
    net_num <- c()
    for (i in 1:cell_number){
      net_num <- append(net_num,length(rownames(CSlncR_network[[i]])))
    }
    CSlncR_network_validated_random <-random.net.validate(net_num,lncRTarget_graph,lncRNA_scRNA_norm_filter,mRNA_scRNA_norm_filter)
    # save(CSlncR_network_validated_random,file = "CSlncR_network_validated_random.RData")
    # load("CSlncR_network_validated_random.RData")

    # LncRNA2Target 
    CSlncR_network_graph_LncRNA2Target <- lapply(seq(CSlncR_network_LncRNA2Target), function(i) make_graph(c(t(CSlncR_network_LncRNA2Target[[i]][, 1:2])), directed = TRUE))
    CSlncR_network_validated_LncRNA2Target <- lapply(seq(CSlncR_network_LncRNA2Target), function(i) as_data_frame(CSlncR_network_graph_LncRNA2Target[[i]] %s% lncRTarget_graph))
    # save(CSlncR_network_validated_LncRNA2Target,file = "CSlncR_network_validated_LncRNA2Target.RData")
    # load("CSlncR_network_validated_LncRNA2Target.RData") 

    # NPInterPrediction 
    CSlncR_network_graph_NPInterPrediction  <- lapply(seq(CSlncR_network_NPInterPrediction ), function(i) make_graph(c(t(CSlncR_network_NPInterPrediction [[i]][, 1:2])), directed = TRUE))
    CSlncR_network_validated_NPInterPrediction  <- lapply(seq(CSlncR_network_NPInterPrediction ), function(i) as_data_frame(CSlncR_network_graph_NPInterPrediction [[i]] %s% lncRTarget_graph))
    # save(CSlncR_network_validated_NPInterPrediction,file = "CSlncR_network_validated_NPInterPrediction.RData")
    # load("CSlncR_network_validated_NPInterPrediction.RData") 

# downstream [2] ASD
## ASD-related cell-specific lncRNA-mRNA interactions. the list of ASD-related lncRNAs and mRNAs (ASD variable) is from SFARI tools, respectively    
    CSlncR_network_ASD <- lapply(seq(CSlncR_network), function(i) CSlncR_network[[i]][intersect(which(CSlncR_network[[i]][, 1] %in% as.matrix(ASD)), which(CSlncR_network[[i]][, 2] %in% as.matrix(ASD))), ])

# downstream [3] interactions  ####[1] Result analysis: The lncRNA regulation in each cell is unique####
    ## Overlap of cell-specific lncRNA-mRNA interactions across cells  
    Overlap_network <- Overlap_net_interaction(CSlncR_network, Intersect_num = round(length(CSlncR_network)*0.9)) #default: 0.9
    Overlap_network_union <- Overlap_net_interaction(CSlncR_network, Intersect_num = 1)
    Overlap_network_two <- Overlap_net_interaction(CSlncR_network, Intersect_num = 2) 
    
    Overlap_network_union_graph <- make_graph(c(t(Overlap_network_union[, 1:2])), directed = TRUE)
    Overlap_network_two_graph <- make_graph(c(t(Overlap_network_two[, 1:2])), directed = TRUE)
    Overlap_network_rewired <- as_data_frame(Overlap_network_union_graph %m% Overlap_network_two_graph)
    
    Overlap_network_graph <- make_graph(c(t(Overlap_network[, 1:2])), directed = TRUE)
    Overlap_network_conserved_validated <- as_data_frame(Overlap_network_graph %s% lncRTarget_graph)
    Overlap_network_rewired_graph <- make_graph(c(t(Overlap_network_rewired[, 1:2])), directed = TRUE)
    Overlap_network_rewired_validated <- as_data_frame(Overlap_network_rewired_graph %s% lncRTarget_graph)
    
    # write.csv(Overlap_network_rewired,"./tmp.csv")

# downstream [4] ASD-related regulatory network
    ## ASD-related conserved and rewired lncRNA-mRNA regulatory network    
    Overlap_network_ASD <- Overlap_network[intersect(which(Overlap_network[, 1] %in% as.matrix(ASD)), which(Overlap_network[, 2] %in% as.matrix(ASD))), ]
    Overlap_network_rewired_ASD <- Overlap_network_rewired[intersect(which(Overlap_network_rewired[, 1] %in% as.matrix(ASD)), which(Overlap_network_rewired[, 2] %in% as.matrix(ASD))), ]
    
# downstream [5] hub  ####[1] Result analysis: The lncRNA regulation in each cell is unique####
    ## Identifying cell-specific hub lncRNAs    
    CSlncR_network_outdegree <- lapply(seq(CSlncR_network), function(i) degree(CSlncR_network_graph[[i]], mode="out"))
    hub_lncRNAs <- lapply(seq(CSlncR_network), function(i) names(sort(CSlncR_network_outdegree[[i]][which(CSlncR_network_outdegree[[i]]!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSlncR_network_outdegree[[i]]!=0)))])
    
    # GW16_to_GW23.5
    CSlncR_network_GW16_to_GW23.5_outdegree <- lapply(seq(CSlncR_network_GW16_to_GW23.5), function(i) degree(CSlncR_network_GW16_to_GW23.5_graph[[i]], mode="out"))
    hub_lncRNAs_GW16_to_GW23.5 <- lapply(seq(CSlncR_network_GW16_to_GW23.5), function(i) names(sort(CSlncR_network_GW16_to_GW23.5_outdegree[[i]][which(CSlncR_network_GW16_to_GW23.5_outdegree[[i]]!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSlncR_network_GW16_to_GW23.5_outdegree[[i]]!=0)))])
    
    ## Overlap of cell-specific hub lncRNAs across cells    
    Overlap_hub_lncRNAs_conserved <- Overlap_hub(hub_lncRNAs, Intersect_num = round(length(hub_lncRNAs)*0.9)) #default: 0.9
    Overlap_hub_lncRNAs_union <- Overlap_hub(hub_lncRNAs, Intersect_num = 1)
    Overlap_hub_lncRNAs_two <- Overlap_hub(hub_lncRNAs, Intersect_num = 2)
    Overlap_hub_lncRNAs_rewired <- setdiff(Overlap_hub_lncRNAs_union, Overlap_hub_lncRNAs_two)
    
    # write.csv(Overlap_hub_lncRNAs_rewired,"./tmp.csv")
    
# downstream [6] ASD-related hub       
    ## ASD-related cell-specific hub lncRNAs    
    hub_lncRNAs_ASD <- lapply(seq(hub_lncRNAs), function(i) hub_lncRNAs[[i]][which(hub_lncRNAs[[i]] %in% as.matrix(ASD))])
    
    ## ASD-related conserved and rewired hub lncRNAs    
    Overlap_hub_lncRNAs_conserved_ASD <- Overlap_hub_lncRNAs_conserved[which(Overlap_hub_lncRNAs_conserved %in% as.matrix(ASD))]    
    Overlap_hub_lncRNAs_rewired_ASD <- Overlap_hub_lncRNAs_rewired[which(Overlap_hub_lncRNAs_rewired %in% as.matrix(ASD))]
    
# downstream [7]  similarity matrix    
    ## Calculating similarity matrix of cell-specific lncRNA-mRNA regulatory network across cells    
    CSlncR_network_Sim <- Sim.network(CSlncR_network, CSlncR_network, directed = TRUE)
    #CSlncR_network_GW16_to_GW23.5_Sim <- Sim.network(CSlncR_network_GW16_to_GW23.5, CSlncR_network_GW16_to_GW23.5, directed = TRUE)
    
    ## Calculating similarity matrix of cell-specific hub lncRNAs across cells    
    CSlncR_hub_Sim <- Sim.hub(hub_lncRNAs, hub_lncRNAs)
    #CSlncR_hub_GW16_to_GW23.5_Sim <- Sim.hub(hub_lncRNAs_GW16_to_GW23.5, hub_lncRNAs_GW16_to_GW23.5)
    
    ## Calculating similarity matrix of expression across cells
    CSlncR_expression_Sim <- cor(t(lncRNA_mRNA_scRNA_norm_filter))
    #CSlncR_expression_Sim_GW16_to_GW23.5 <- cor(t(lncRNA_mRNA_scRNA_norm_filter_GW16_to_GW23.5))
    
    #save(CSlncR_network_Sim,CSlncR_hub_Sim,CSlncR_expression_Sim,file = "sim.RData")
    #load("sim.RData")
    #save(CSlncR_network_GW16_to_GW23.5_Sim,CSlncR_hub_GW16_to_GW23.5_Sim,CSlncR_expression_Sim_GW16_to_GW23.5,file = "GW16_to_GW23.5_Sim.RData")
    #load("GW16_to_GW23.5_Sim.RData")
    
##############################################################################################################################
######## downstream [8-9-10] cell analysis ###### [5] Result analysis: CSlncR helps to understand cell-cell crosstalk  #######
# downstream [8]  cell-cell crosstalk      
## Identifying cell-cell crosstalk network in terms of network similarity matrix    
    CSlncR_network_adjacency_matrix <- ifelse(CSlncR_network_Sim > median(CSlncR_network_Sim[lower.tri(CSlncR_network_Sim)]), 1, 0)
    diag(CSlncR_network_adjacency_matrix) <- 0
    colnames(CSlncR_network_adjacency_matrix) <- rownames(CSlncR_network_adjacency_matrix) <- rownames(lncRNA_scRNA_norm_filter) 
    CSlncR_network_adjacency_matrix_graph <- graph_from_adjacency_matrix(CSlncR_network_adjacency_matrix, mode = "undirected")

## Identifying cell-cell crosstalk network in terms of hub lncRNA similarity matrix    
    CSlncR_hub_adjacency_matrix <- ifelse(CSlncR_hub_Sim > median(CSlncR_hub_Sim[lower.tri(CSlncR_hub_Sim)]), 1, 0)
    diag(CSlncR_hub_adjacency_matrix) <- 0
    colnames(CSlncR_hub_adjacency_matrix) <- rownames(CSlncR_hub_adjacency_matrix) <- rownames(lncRNA_scRNA_norm_filter)
    CSlncR_hub_adjacency_matrix_graph <- graph_from_adjacency_matrix(CSlncR_hub_adjacency_matrix, mode = "undirected")

    adjacency_ <-  CSlncR_network_adjacency_matrix 
    # adjacency_ <-  CSlncR_hub_adjacency_matrix  
    crosstalk_ <- data.frame(cell1=0,cell2=0)
    cnt <- 0
    for(i in seq(dim(adjacency_)[1])){
      for (j in seq(i)) {
        if (adjacency_[i,j] == 1){
          cnt <- cnt + 1
          crosstalk_[cnt,1] <- rownames(adjacency_)[i]
          crosstalk_[cnt,2] <- colnames(adjacency_)[j]
        }
        
      }
    }
    # write.csv(crosstalk_,"./tmp.csv")
    
# downstream [9]  hub cells   
## Identifying hub cells in terms of network similarity    
    CSlncR_network_cell_degree <- degree(CSlncR_network_adjacency_matrix_graph)
    names(CSlncR_network_cell_degree) <- rownames(lncRNA_scRNA_norm_filter)
    CSlncR_network_hub_cells <- names(sort(CSlncR_network_cell_degree[which(CSlncR_network_cell_degree!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSlncR_network_cell_degree!=0)))]

    # write.csv(CSlncR_network_cell_degree[CSlncR_network_hub_cells],"./tmp.csv")

## Identifying hub cells in terms of hub lncRNA similarity    
    CSlncR_hub_cell_degree <- degree(CSlncR_hub_adjacency_matrix_graph)
    names(CSlncR_hub_cell_degree) <- rownames(lncRNA_scRNA_norm_filter)
    CSlncR_hub_hub_cells <- names(sort(CSlncR_hub_cell_degree[which(CSlncR_hub_cell_degree!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSlncR_hub_cell_degree!=0)))]

    # write.csv(CSlncR_hub_cell_degree[CSlncR_hub_hub_cells],"./tmp.csv")
    
# downstream [10]  cell-cell crosstalk modules  
## Identifying cell-cell crosstalk modules in terms of network similarity matrix   
    CSlncR_network_cell_module <- netModule(CSlncR_network_adjacency_matrix_graph %>% as_data_frame)

## Identifying cell-cell crosstalk modules in terms of hub lncRNA similarity matrix    
    CSlncR_hub_cell_module <- netModule(CSlncR_hub_adjacency_matrix_graph %>% as_data_frame)
   
    # save id
    cell_module <- CSlncR_network_cell_module[[1]]
    # cell_module <- CSlncR_hub_cell_module[[1]]
    tmp_module = vector()
    for (ii in seq(length(cell_module))) {
      #tmp_module <- append(tmp_module,as.numeric(strsplit(cell_module[ii],"_")[[1]][2])) # GW16: "_";
      tmp_module <- append(tmp_module,as.character(strsplit(cell_module[ii],"[.]")[[1]][2])) # GW19.5,GW23.5: "[.]";
      #tmp_module <- append(tmp_module,as.character(substr(cell_module[ii],5,100))) # GW21
      #tmp_module <- append(tmp_module,as.character(cell_module[ii])) # GW20.5
    }
    tmp_module <- sort(tmp_module) # GW16,GW19.5,GW20.5,GW21
    #write.csv(tmp_module,"./tmp.csv") 
#######################################################################################################################

########################################################################################################################    
############## downstream [11] cell type-specific lncRNAs analysis 细胞类型特异的lncRNAs分析##############################
############ [2] Result analysis: The cell type-specific regulation across single‑cells of human brain region NCX  ######
## downstream [11] The cell type-specific lncRNAs analysis    
    # 关于细胞类型及marker的来源：原始数据的文章
    # endothelia "LINC00339", "TRIM52-AS1"
    lncR_endothelia  <- c("LINC-MILR1-3","SLC38A3","LINC00152","RP11-401P9.4","MIR4435-1HG",
                          "LINC00339","RP11-483C6.1","AP000459.4","AC127904.2","RP11-161M6.2",
                          "RP11-417F21.1", "TRIM52-AS1","CTD-2081C10.7","RP11-296I10.3","RP11-532M24.1")
    
    # radial glia "LINC00943","MAGI2-AS3", "RUSC1-AS1"
    lncR_radial_glia <- c("Z83001.1", "RP11-731J8.2","LINC00943", "RP3-418C23.2", "RP11-1002K11.1",
                          "MAGI2-AS3","RP11-421L21.3", "LINC-FZD8-3", "LINC-FZD8-1", "LINC00263", 
                          "EIF3J-AS1", "LOC646329","LINC-KREMEN1-1","RUSC1-AS1","DGKK")
    
    # dividing radial glia "THAP9-AS1"
    lncR_dividing_radial_glia <- c("UHRF1", "CTRR-175P5.4","RP11-138A9.1", "RP11-849F2.9", "RP11-143K11.1",
                          "AC004447.2","SNORA59B", "CTC-503J8.6", "RP11-138A9.2", "RP11-95D17.1", 
                          "THAP9-AS1", "SNHG1","CTD-2017D11.1","RP11-58B17.2","DYNLL1-AS1")
    #放射状胶质细胞特异的lncRNAs分析
    lncR_radial_glia_family <- c("LINC00943","MAGI2-AS3", "RUSC1-AS1", "THAP9-AS1")
    
    # intermediate progenitors "DGCR11"
    lncR_intermediate_progenitors <- c("LINC-TMEM200C-1","RP11-798G7.8","RP11-351J23.1-AS1","RP3-326L13.3","CTD-2245E15.3",
                                       "C1orf132","AC084018.1","RP11-73O6.3","RP11-594N15.3","RP11-436D23.1",
                                       "AC083884.8","DGCR11","RP11-456K23.1","RP6-24A23.3","RP1-20C7.6")
    
    # newborn neurons "INHBA-AS1","MYT1L-AS1","KIF9-AS1",
    lncR_newborn_neurons <- c("RP5-1024G6.8","LINC-PTCHD2-3","RP11-513M16.8","RP11-661O13.1","RP11-524C21.2",
                              "RP11-356K23.1","LINC01105","INHBA-AS1","MYT1L-AS1","KIF9-AS1",
                              "RP11-1006G14.4","RP11-296O14.3","RP11-452L6.5","CTD-3099C6.9","RP11-452H21.4")
    
    # maturing excitatory "MIR137HG","PWAR6","SIK3-IT1","NAV2-AS3","DAPK1-IT1",
    lncR_maturing_excitatory_neurons <- c("MIR137HG","LINC00599","PWAR6","SIK3-IT1","RP11-53O19.3",
                                  "RP11-402L6.1","RP11-18I14.10","RP11-486F17.1","NAV2-AS3","DAPK1-IT1",
                                  "RP11-397O4.1","RP11-64K12.10","LINC00643","RP3-462E2.5","LINC-TMEM182-5")
    
    # inhibitory interneurons  "DLX6-AS1","SOX2-OT","MEG3", 
    lncR_inhibitory_interneurons <- c("DLX6-AS1","RP11-588P7.1","SOX2-OT","GS1-18A18.1","MEG3",
                                      "LINC-DKFZP761K2322-2","GRIP2","AC087393.1","LINC00966","RP11-450H6.3",
                                      "RP13-514E23.1","RP11-379H18.1","RP11-69E11.4","AC012358.8","LINC-TBCC-1")
    
    
    lncR_cell_type_ <- c("LINC00339", "TRIM52-AS1","LINC00943","MAGI2-AS3", "RUSC1-AS1",
                        "THAP9-AS1","DGCR11","INHBA-AS1","MYT1L-AS1","KIF9-AS1",
                        "MIR137HG","PWAR6","SIK3-IT1","NAV2-AS3","DAPK1-IT1","DLX6-AS1","SOX2-OT","MEG3")
    
    lncR_cell_type <- c(lncR_endothelia,lncR_radial_glia,lncR_dividing_radial_glia,lncR_intermediate_progenitors,lncR_newborn_neurons,lncR_maturing_excitatory_neurons,lncR_inhibitory_interneurons)

    # Extracting lncR_cell_type related cell-specific lncRNA-mRNA regulatory networks
    CSlncR_network_lncRcelltype <- lapply(seq(CSlncR_network), function(i) CSlncR_network[[i]][which(CSlncR_network[[i]][, 1] %in% as.matrix(lncR_cell_type)), ])
    
    # Extracting conserved and rewired lncRNA-mRNA regulatory networks associated with the lncR_cell_type
    Overlap_network_lncRcelltype <- Overlap_network[which(Overlap_network[, 1] %in% as.matrix(lncR_cell_type)), ]
    Overlap_network_rewired_lncRcelltype <- Overlap_network_rewired[which(Overlap_network_rewired[, 1] %in% as.matrix(lncR_cell_type)), ]
    
    # write.csv(Overlap_network_lncRcelltype,"./tmp.csv")
    
    # Extracting ASD-related conserved and rewired lncRNA-mRNA regulatory networks associated with the lncR_cell_type
    CSlncR_network_lncRcelltype_ASD <- lapply(seq(CSlncR_network_lncRcelltype), function(i) CSlncR_network_lncRcelltype[[i]][intersect(which(CSlncR_network_lncRcelltype[[i]][, 1] %in% as.matrix(ASD)), which(CSlncR_network_lncRcelltype[[i]][, 2] %in% as.matrix(ASD))), ])
    
    # Experimentally validated lncRNA-mRNA interactions associated with the lncRcelltype
    CSlncR_network_lncRcelltype_graph <- lapply(seq(CSlncR_network_lncRcelltype), function(i) make_graph(c(t(CSlncR_network_lncRcelltype[[i]][, 1:2])), directed = TRUE))
    CSlncR_network_lncRcelltype_validated <- lapply(seq(CSlncR_network_lncRcelltype), function(i) as_data_frame(CSlncR_network_lncRcelltype_graph[[i]] %s% lncRTarget_graph))
    
    ## Difference of pvalue from the lncRcelltype regulation    
    CSlncR_pvalue <- function(tmp,cell_number) {
    res_pvalue <- matrix(NA, nrow = cell_number, ncol = cell_number)
    for (i in seq(cell_number)) {
      for (j in seq(cell_number)) {
      
        cell_1 <- tmp[[i]]
        cell_2 <- tmp[[j]]

        cell_1_ <- paste(cell_1[,1], cell_1[,2], sep = "&")
        cell_2_ <- paste(cell_2[,1], cell_2[,2], sep = "&")
        cell_1_2_union <- union(cell_1_,cell_2_)
        cell_1_2_intersect <- intersect(cell_1_,cell_2_)
        cell_1_remaining_num <-  length(cell_1_) -  length(cell_1_2_intersect)
        cell_2_remaining_num <-  length(cell_2_) -  length(cell_1_2_intersect)
        
        intersect_ <- rep(c(1),length(cell_1_2_intersect))
       
        cell_1_remaining_1_ <- rep(c(1),cell_1_remaining_num)
        cell_1_remaining_0_ <- rep(c(0),cell_1_remaining_num)
        
        cell_2_remaining_1_ <- rep(c(1),cell_2_remaining_num)
        cell_2_remaining_0_ <- rep(c(0),cell_2_remaining_num)
        
        cell_1_regulation <- c(intersect_,cell_1_remaining_1_,cell_2_remaining_0_)
        cell_2_regulation <- c(intersect_,cell_1_remaining_0_,cell_2_remaining_1_)
        
        res_cell1 <- cell_1_regulation
        res_cell2 <- cell_2_regulation

        res_pvalue[i, j] <- ks.test(res_cell1, res_cell2)$p.value
      }
    }
    
    return(res_pvalue) } 
    
    res_predict_pvalue <- CSlncR_pvalue(tmp= CSlncR_network_lncRcelltype,cell_number=cell_number) 
    #res_validated_pvalue <- CSlncR_pvalue(tmp= CSlncR_network_lncRcelltype_validated,cell_number=cell_number) 
    res_ASD_pvalue <- CSlncR_pvalue(tmp= CSlncR_network_lncRcelltype_ASD,cell_number=cell_number) 

    rownames(res_predict_pvalue) <- colnames(res_predict_pvalue) <- paste("Cell",c(1:cell_number),sep=" ")
    #rownames(res_validated_pvalue) <- colnames(res_validated_pvalue) <- paste("Cell",c(1:cell_number),sep=" ")
    rownames(res_ASD_pvalue) <- colnames(res_ASD_pvalue) <- paste("Cell",c(1:cell_number),sep=" ")
    test.value <- matrix(0, nrow = cell_number, ncol = cell_number)
    rownames(test.value) <- colnames(test.value) <- paste("Cell",c(1:cell_number),sep=" ")
    
    library(corrplot)
    par(mfrow=c(1,2)) # 形成1行、2列的图形矩阵
    corrplot(test.value, p.mat = res_predict_pvalue, method = "square", diag = FALSE, type = "upper", title="A  Difference in predicted targets of cell type-specific",
                cl.pos="n", sig.level = c(.001, .01, .05), pch.cex = 1, mar=c(0,0,1.5,0), insig = "label_sig", pch.col = "black")

    corrplot(test.value, p.mat = res_predict_pvalue, method = "square", diag = FALSE, type = "upper", title="A  Difference in predicted targets of cell type-specific",
             cl.pos="n", sig.level = c(.001, .01, .05), pch.cex = 1, mar=c(0,0,1.5,0), insig = "label_sig", pch.col = "black", tl.pos = 'n')
    
    corrplot(test.value, p.mat = res_ASD_pvalue, method = "square", diag = FALSE, type = "upper", title="B  Difference in ASD-related targets of cell type-specific",
             cl.pos="n", sig.level = c(.001, .01, .05), pch.cex = 1, mar=c(0,0,1.5,0), insig = "label_sig", pch.col = "black")

    corrplot(test.value, p.mat = res_ASD_pvalue, method = "square", diag = FALSE, type = "upper", title="B  Difference in ASD-related targets of cell type-specific",
             cl.pos="n", sig.level = c(.001, .01, .05), pch.cex = 1, mar=c(0,0,1.5,0), insig = "label_sig", pch.col = "black", tl.pos = 'n')
    
    pred_ <- data.frame(id=c("conserved","rewired"),value=c(dim(Overlap_network_lncRcelltype)[1],dim(Overlap_network_rewired_lncRcelltype)[1]))
    ggplot(pred_, aes(x = id, y = value)) +
      geom_bar(fill = ifelse(pred_$id=="conserved","blue","red"), stat = "identity", width = 0.6) +
      xlab("lncRNA regulation type") +
      ylab("#Predicted targerts of cell type-specific") +
      geom_text(aes(label=value,y=value+5),size=3) +
      labs(title='E') # GW16: A, GW19.5: B, GW20.5: C, GW21: D, GW23.5: E
    
    # save         11 5.5 inches-pdf, 1100 550 tiff # GW16, GW19.5,GW21
    # save         16 8.0 inches-pdf, 1600 800 tiff # GW20.5, GW23.5
    # each plot is  5.5 5.5           550 550
    
    ## Enrichment analysis of conserved and rewired the lncRcelltype regulation   
    CSlncR_network_lncRcelltype_conserved = list(Overlap_network_lncRcelltype)
    CSlncR_network_lncRcelltype_rewired = list(Overlap_network_rewired_lncRcelltype)
    
    lncRcelltype_conserved_gene_list=c()
    for (i in seq(nrow(CSlncR_network_lncRcelltype_conserved[[1]]))) {
      for (j in seq(ncol(CSlncR_network_lncRcelltype_conserved[[1]]))) {
        lncRcelltype_conserved_gene_list=append(x=lncRcelltype_conserved_gene_list, CSlncR_network_lncRcelltype_conserved[[1]][i,j])
      }
      
    }
    lncRcelltype_conserved_gene_list <- list(unique(lncRcelltype_conserved_gene_list))
    
    lncRcelltype_rewired_gene_list=c()
    for (i in seq(nrow(CSlncR_network_lncRcelltype_rewired[[1]]))) {
      for (j in seq(ncol(CSlncR_network_lncRcelltype_rewired[[1]]))) {
        lncRcelltype_rewired_gene_list=append(x=lncRcelltype_rewired_gene_list, CSlncR_network_lncRcelltype_rewired[[1]][i,j])
      }
      
    }
    lncRcelltype_rewired_gene_list <- list(unique(lncRcelltype_rewired_gene_list))
    
    # GO, KEGG and Reactome enrichment analysis
    lncRcelltype_conserved_FEA <- moduleFEA(lncRcelltype_conserved_gene_list, padjustvaluecutoff = 0.05)
    lncRcelltype_rewired_FEA <- moduleFEA(lncRcelltype_rewired_gene_list, padjustvaluecutoff = 0.05)
    
    # ASD enrichment analysis   
    lncRcelltype_conserved_ASD_EA <- module_ASD_EA(lncRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, ASD, lncRcelltype_conserved_gene_list)
    lncRcelltype_rewired_ASD_EA <- module_ASD_EA(lncRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, ASD, lncRcelltype_rewired_gene_list)
    
    # Hallmark enrichment analysis
    m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, human_gene_symbol)
    
    lncRcelltype_conserved_Hallmark <- lapply(seq(lncRcelltype_conserved_gene_list), function(i) enricher(lncRcelltype_conserved_gene_list[[i]], TERM2GENE=m_t2g, minGSSize=1) %>% as.data.frame)
    lncRcelltype_rewired_Hallmark <- lapply(seq(lncRcelltype_rewired_gene_list), function(i) enricher(lncRcelltype_rewired_gene_list[[i]], TERM2GENE=m_t2g, minGSSize=1) %>% as.data.frame)
    
    # Cell marker enrichment analsyis
    # cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
    #   tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
    #   dplyr::select(cellMarker, geneSymbol)    
    cell_markers  <- read.csv("./data/Cell_marker_Human.csv") %>% tidyr::unite("Marker", Tissue.type, Cancer.type, Cell.name, sep=", ") %>% dplyr::select(Marker, Symbol)   
    
    lncRcelltype_conserved_Cellmarker <- lapply(seq(lncRcelltype_conserved_gene_list), function(i) enricher(lncRcelltype_conserved_gene_list[[i]], TERM2GENE=cell_markers, minGSSize=1) %>% as.data.frame)
    lncRcelltype_rewired_Cellmarker <- lapply(seq(lncRcelltype_rewired_gene_list), function(i) enricher(lncRcelltype_rewired_gene_list[[i]], TERM2GENE=cell_markers, minGSSize=1) %>% as.data.frame)

    # write.csv(lncRcelltype_rewired_FEA[[3]][[1]]@result,"./tmp.csv")
    # write.csv(lncRcelltype_rewired_Hallmark[[1]],"./tmp.csv")
    # write.csv(lncRcelltype_rewired_Cellmarker[[1]],"./tmp.csv")
    

################################################################################################################################################      
################################################################################################################################################    
# downstream [12]  Hierarchical cluster analysis  #### [4] Result analysis: CSlncR provides a novel strategy for clustering single-cells ####
    par(mfrow=c(2,2)) # 形成2行、2列的图形矩阵
## Hierarchical cluster analysis of cell-specific lncRNA-mRNA regulatory network    
    rownames(CSlncR_network_GW16_to_GW23.5_Sim) <- colnames(CSlncR_network_GW16_to_GW23.5_Sim) <- rownames(lncRNA_mRNA_scRNA_norm_filter)
    hclust_res_network <- hclust(as.dist(1-CSlncR_network_GW16_to_GW23.5_Sim), "complete")
    
    dend_network <- as.dendrogram(hclust_res_network)
   
    # 5 developments
    dend_network %>% set("branches_k_color", value = c("red", "green", "blue","pink","skyblue"), k=5) %>% 
      set("labels_col", value = c("red", "green", "blue","pink","skyblue"), k=5) %>% 
      set("labels_cex", value = 0.001) %>%  #set("labels_cex", value = 0.7) %>% 
      plot(main="A  Hierarchical clustering analysis using interaction similarity",ylab = "Height")
    
    pre_network = cutree(dend_network,k=5)
    
## Hierarchical cluster analysis of cell-specific hub lncRNAs    
    rownames(CSlncR_hub_GW16_to_GW23.5_Sim) <- colnames(CSlncR_hub_GW16_to_GW23.5_Sim) <- rownames(lncRNA_mRNA_scRNA_norm_filter)
    hclust_res_hub <- hclust(as.dist(1-CSlncR_hub_GW16_to_GW23.5_Sim), "complete")
 
    dend_hub <- as.dendrogram(hclust_res_hub)
    
    # 5 developments
    dend_hub %>% set("branches_k_color", value = c("red", "green", "blue","pink","skyblue"), k=5) %>% 
      set("labels_col", value = c("red", "green", "blue","pink","skyblue"), k=5) %>% 
      set("labels_cex", value = 0.001) %>%  #set("labels_cex", value = 0.7) %>% 
      plot(main="B  Hierarchical clustering analysis using hub lncRNA similarity",ylab = "Height")
    
    pre_hub = cutree(hclust_res_hub,k=5)
    
## Hierarchical cluster analysis of expression similarity  
    rownames(CSlncR_expression_Sim_GW16_to_GW23.5) <- colnames(CSlncR_expression_Sim_GW16_to_GW23.5) <- rownames(lncRNA_mRNA_scRNA_norm_filter_GW16_to_GW23.5)
    hclust_res_expression <- hclust(as.dist(1-CSlncR_expression_Sim_GW16_to_GW23.5), "complete")

    dend_expression <- as.dendrogram(hclust_res_expression)
    
    # 5 developments
    dend_expression %>% set("branches_k_color", value = c("red", "green", "blue","pink","skyblue"), k=5) %>% 
      set("labels_col", value = c("red", "green", "blue","pink","skyblue"), k=5) %>% 
      set("labels_cex", value = 0.001) %>% #set("labels_cex", value = 0.7) %>% 
      plot(main="C  Hierarchical clustering analysis using expression similarity",ylab = "Height")
    
    pre_expression = cutree(hclust_res_expression,k=5)

    # label of GW16 to GW23.5: 1 2 3 4 5
    labels_ <- as.vector(c(rep(1,26),rep(2,26),rep(3,123),rep(4,24),rep(5,77)))
    pre_network <- as.vector(pre_network)
    pre_hub <- as.vector(pre_hub)
    pre_expression <- as.vector(pre_expression)
    
    library(aricode)
    ari_network <- ARI(pre_network,labels_)
    nid_network <- NID(pre_network,labels_)
    nvi_network <- NVI(pre_network,labels_)
    
    ari_hub <- ARI(pre_hub,labels_)
    nid_hub <- NID(pre_hub,labels_)
    nvi_hub <- NVI(pre_hub,labels_)
    
    ari_expression <- ARI(pre_expression,labels_)
    nid_expression <- NID(pre_expression,labels_)
    nvi_expression <- NVI(pre_expression,labels_)
  
## barplot similarity      
    Count1_A <- c(ari_network,nid_network,nvi_network)
    Count2_B <- c(ari_hub,nid_hub,nvi_hub)
    Count3_C <- c(ari_expression,nid_expression,nvi_expression)
    data_clustering  <- t(cbind(cbind(as.data.frame(Count1_A),Count2_B),Count3_C))

    #非堆积柱状图
    barplot(as.matrix(data_clustering),
            names.arg = c("ARI","NID","NVI"),
            main="D    Clustering results using interaction - hub lncRNA - expression similarity",
            xlab="Clustering Measures",ylab="Value",
            col=c("red","green","blue"),
            #legend=c("Interaction","Hub lncRNA","Expression"),args.legend = list(x = "topleft"),
            ylim=c(0,1),
            beside=TRUE)
    # Add the legend to the chart.
    legend("left", c("Interaction","Hub lncRNA","Expression"), cex=1, fill=c("red","green","blue"))
    box() #边框
    
    # save 13 7.5  inches-pdf, 1300 750 tiff
################################################################################################################################################      
################################################################################################################################################ 
    
# downstream [13] Similarity network plot   #### [1] Result analysis: The lncRNA regulation in each cell is unique ####   
## Similarity network plot in terms of cell-specific lncRNA-mRNA regulatory netowork    
    par(mfrow=c(1,2)) # 形成1行、2列的图形矩阵
    rownames(CSlncR_network_Sim) <- colnames(CSlncR_network_Sim) <- paste("Cell",c(1:cell_number),sep=" ")
    corrplot(CSlncR_network_Sim, method = "pie", type = "upper",title = "A   Cell-specific lncRNA-mRNA interactions", mar=c(0,0,2,0),diag = FALSE, cl.lim = c(0, 1), tl.cex = 1)
    corrplot(CSlncR_network_Sim, method = "pie", type = "upper",title = "A   Cell-specific lncRNA-mRNA interactions", mar=c(0,0,2,0),diag = FALSE, cl.lim = c(0, 1), tl.pos = 'n')

## Similarity network plot in terms of cell-specific hub lncRNAs
    rownames(CSlncR_hub_Sim) <- colnames(CSlncR_hub_Sim) <- paste("Cell",c(1:cell_number),sep=" ")
    corrplot(CSlncR_hub_Sim, method = "pie", type = "upper", title =    "B       Cell-specific hub lncRNAs", diag = FALSE, mar=c(0,0,2,0),cl.lim = c(0, 1), tl.cex = 1)
    corrplot(CSlncR_hub_Sim, method = "pie", type = "upper", title =    "B       Cell-specific hub lncRNAs", diag = FALSE, mar=c(0,0,2,0),cl.lim = c(0, 1), tl.pos = 'n')

    # save 9.5 4.5 inches-pdf, 950 450 tiff # GW16,GW19.5,GW21
    # save 16 8.0 inches-pdf, 1600 800 tiff # GW20.5, GW23.5

# downstream [14] ggplot plot patchwork   #### [1] Result analysis: The lncRNA regulation in each cell is unique ####   
## Stem plots
index_net <- data.frame(value = unlist(lapply(seq(CSlncR_network), function(i) nrow(CSlncR_network[[i]]))), id = seq(cell_number))

col1 <- rep("#FF9999", cell_number)
p1 <- ggplot(index_net, aes(x = id, y = value)) +
    geom_point(aes(color = col1), size = 5) +
    geom_bar(aes(fill = col1), stat = "identity", width = 0.2) +
    #theme_bw(base_family = "Times") +
    xlab("Single-cell ID") +
    ylab("#Predicted lncRNA-mRNA interactions") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position="none", 
	  panel.border = element_blank(),
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) +
    scale_x_continuous(breaks = seq(1, cell_number, 1)) 

index_validated_net <- data.frame(value = unlist(lapply(seq(CSlncR_network_validated), function(i) nrow(CSlncR_network_validated[[i]])/nrow(CSlncR_network[[i]])*100)), id = seq(cell_number))

col2 <- rep("plum4", cell_number)
p2 <- ggplot(index_validated_net, aes(x = id, y = value)) +
    geom_point(aes(color = col2), size = 5) +
    geom_bar(aes(fill = col2), stat = "identity", width = 0.2) +
    scale_fill_manual(values=c("plum4"), aesthetics = "fill") +
    scale_colour_manual(values=c("plum4"), aesthetics = "colour") +
    xlab("Single-cell ID") +
    ylab("%Validated lncRNA-mRNA interactions") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position="none", 
	  panel.border = element_blank(),
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) + 
    scale_x_continuous(breaks = seq(1, cell_number, 1)) 

index_ASD_net <- data.frame(value = unlist(lapply(seq(CSlncR_network_ASD), function(i) nrow(CSlncR_network_ASD[[i]])/nrow(CSlncR_network[[i]])*100)), id = seq(cell_number))

col3 <- rep("blue", cell_number)
p3 <- ggplot(index_ASD_net, aes(x = id, y = value)) +
    geom_point(aes(color = col3), size = 5) +
    geom_bar(aes(fill = col3), stat = "identity", width = 0.2) +
    scale_fill_manual(values=c("blue"), aesthetics = "fill") +
    scale_colour_manual(values=c("blue"), aesthetics = "colour") +
    xlab("Single-cell ID") +
    ylab("%ASD-related lncRNA-mRNA interactions") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position="none", 
	  panel.border = element_blank(),
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) +
    scale_x_continuous(breaks = seq(1, cell_number, 1))

index_ASD_hub <- data.frame(value = unlist(lapply(seq(hub_lncRNAs_ASD), function(i) length(hub_lncRNAs_ASD[[i]])/length(hub_lncRNAs[[i]])*100)), id = seq(cell_number))

col4 <- rep("green", cell_number)
p4 <- ggplot(index_ASD_hub, aes(x = id, y = value)) +
    geom_point(aes(color = col4), size = 5) +
    geom_bar(aes(fill = col4), stat = "identity", width = 0.2) +
    scale_fill_manual(values=c("green"), aesthetics = "fill") +
    scale_colour_manual(values=c("green"), aesthetics = "colour") +
    xlab("Single-cell ID") +
    ylab("%ASD-related hub lncRNAs") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position="none", 
	  panel.border = element_blank(),
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) +
    scale_x_continuous(breaks = seq(1, cell_number, 1))

library(patchwork)
(p1+p2)/(p3+p4) + plot_annotation(tag_levels = 'A')
# save 10 10 inches-pdf, 1000 1000 tiff  # GW16, GW19.5, GW21
# save 13 13 inches-pdf, 1300 1300 tiff  # GW20.5,GW23.5

# downstream [15] comparison of CSlncR-Random-TargetScan-NPInter #### [3] Result analysis: CSlncR is effective in predicting cell-specific lncRNA targets  ####   
## Stem plots
index_validated_net <- data.frame(value = unlist(lapply(seq(CSlncR_network_validated), function(i) nrow(CSlncR_network_validated[[i]])/nrow(CSlncR_network[[i]])*100)), id = seq(cell_number))
index_validated_net_ <- data.frame(value = unlist(lapply(seq(CSlncR_network_validated_), function(i) nrow(CSlncR_network_validated_[[i]])/nrow(CSlncR_network_[[i]])*100)), id = seq(cell_number))
index_validated_net_random <- index_validated_net
index_validated_net_random[1] <- CSlncR_network_validated_random*100
index_validated_net_LncRNA2Target <- data.frame(value = unlist(lapply(seq(CSlncR_network_validated_LncRNA2Target), function(i) nrow(CSlncR_network_validated_LncRNA2Target[[i]])/nrow(CSlncR_network_LncRNA2Target[[i]])*100)), id = seq(cell_number))
index_validated_net_NPInterPrediction <- data.frame(value = unlist(lapply(seq(CSlncR_network_validated_NPInterPrediction), function(i) nrow(CSlncR_network_validated_NPInterPrediction[[i]])/nrow(CSlncR_network_NPInterPrediction[[i]])*100)), id = seq(cell_number))

data1 <- rbind(index_validated_net,index_validated_net_)
data1$Methods <- c(rep("CSlncR",cell_number), rep("CSlncR without prior knowledge",cell_number))
data2 <- rbind(index_validated_net,index_validated_net_random)
data2$Methods <- c(rep("CSlncR",cell_number), rep("Random",cell_number))
data3 <- rbind(index_validated_net,index_validated_net_LncRNA2Target)
data3$Methods <- c(rep("CSlncR",cell_number), rep("LncRNA2Target",cell_number))
data4 <- rbind(index_validated_net,index_validated_net_NPInterPrediction)
data4$Methods <- c(rep("CSlncR",cell_number), rep("NPInterPrediction",cell_number))

p1 <- ggplot(data1,aes(x=id,y=value,group=Methods,color=Methods))+
  geom_line(size=0.5)+
  geom_point(aes(shape=Methods,color=Methods),size=2)+
  xlab("Single-cell ID") +
  ylab("%Validated numbers") +
  theme(legend.position="top",
        axis.text.x = element_text(face = "bold",size=10),
        axis.text.y = element_text(face = "bold",size=10),
        axis.title.x = element_text(face = "bold",size=12),
        axis.title.y = element_text(face = "bold",size=12)) 

p2 <- ggplot(data2,aes(x=id,y=value,group=Methods,color=Methods))+
  geom_line(size=0.5)+
  geom_point(aes(shape=Methods,color=Methods),size=2)+
  xlab("Single-cell ID") +
  ylab("%Validated numbers") +
  theme(legend.position="top",
        axis.text.x = element_text(face = "bold",size=10),
        axis.text.y = element_text(face = "bold",size=10),
        axis.title.x = element_text(face = "bold",size=12),
        axis.title.y = element_text(face = "bold",size=12)) 

p3 <- ggplot(data3,aes(x=id,y=value,group=Methods,color=Methods))+
  geom_line(size=0.5)+
  geom_point(aes(shape=Methods,color=Methods),size=2)+
  xlab("Single-cell ID") +
  ylab("%Validated numbers") +
  theme(legend.position="top",
      axis.text.x = element_text(face = "bold",size=10),
      axis.text.y = element_text(face = "bold",size=10),
      axis.title.x = element_text(face = "bold",size=12),
      axis.title.y = element_text(face = "bold",size=12)) 

p4 <- ggplot(data4,aes(x=id,y=value,group=Methods,color=Methods))+
  geom_line(size=0.5)+
  geom_point(aes(shape=Methods,color=Methods),size=2)+
  xlab("Single-cell ID") +
  ylab("%Validated numbers") +
  theme(legend.position="top",
        axis.text.x = element_text(face = "bold",size=10),
        axis.text.y = element_text(face = "bold",size=10),
        axis.title.x = element_text(face = "bold",size=12),
        axis.title.y = element_text(face = "bold",size=12)) 

library(patchwork)
(p1+p2)/(p3+p4) + plot_annotation(tag_levels = 'A')
# save 10 7.2 inches-pdf, 1000 720 tiff

#Supplement 1:average values and paired t-test
cn <- cell_number+1
cn2 <- cell_number*2
mean(data1[1][1:cell_number,]) # 0.4123549
mean(data1[1][cn:cn2,]) # 0.02520622
t.test(data1[1][1:cell_number,],data1[1][cn:cn2,],paired = TRUE) # p-value = 2.20E-16
mean(data2[1][1:cell_number,]) # 0.4123549
mean(data2[1][cn:cn2,]) # 0.01902152
t.test(data2[1][1:cell_number,],data2[1][cn:cn2,],paired = TRUE) # p-value = 2.2E-16
mean(data3[1][1:cell_number,]) # 0.4123549
mean(data3[1][cn:cn2,]) #  0.339597
t.test(data3[1][1:cell_number,],data3[1][cn:cn2,],paired = TRUE) # p-value = 5.36E-05
mean(data4[1][1:cell_number,]) # 0.4123549
mean(data4[1][cn:cn2,]) # 0.01225273
t.test(data4[1][1:cell_number,],data4[1][cn:cn2,],paired = TRUE) # p-value = 2.2E-16

#Supplement 2: deduplication, get the relevant proportions of conserved and rewired
Netlist = CSlncR_network
tmp = NULL
for (i in seq(Netlist)){
  Interin <- Netlist[[i]]
  Interin_paste <- paste(Interin[, 1], Interin[, 2], sep = "_")
  tmp <- c(tmp,Interin_paste)
}
num = unique(tmp) 
percentage_conserved = length(rownames(Overlap_network)) / length(num) 
percentage_rewired = length(rownames(Overlap_network_rewired)) /length(num) 

Netlist = hub_lncRNAs
tmp = NULL
for (i in seq(Netlist)){
  Interin <- Netlist[[i]]
  for (j in seq(Interin)) {
    tmp <- c(tmp,Interin[j])
  }
}
num = unique(tmp)
percentage_conserved = length(Overlap_hub_lncRNAs_conserved) / length(num) 
percentage_rewired = length(Overlap_hub_lncRNAs_rewired) /length(num) 

