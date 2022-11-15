############################################################################################
############## Utility functions for exploring cell-type-specific lncRNA regulation ########
############################################################################################

## The original version of the function is written in Matlab at https://github.com/wys8c764/CSN 
## (Dai H, Li L, Zeng T, Chen L. Cell-specific network constructed by single-cell RNA sequencing data. Nucleic Acids Res. 2019, doi: 10.1093/nar/gkz172.), 
## We reimplement the function in R for single cell RNA sequencing data.  
# gx and gy: Gene expression values of gene x (a vector) and gene y (a vector) in n cells 
# boxsize: Size of neighborhood (0.1 in default)
# Output: res is a vector, the normalized statistic of edge gx-gy in n cells
csn_edge <- function(gx, gy, boxsize = 0.1){

# Define the neighborhood of each plot
    n <- length(gx)
    upper <- zeros(1, n)
    lower <- zeros(1, n)
    a <- zeros(2, n)
    B <- list()

    for (i in seq_len(2)){
        g <- gx*(i==1)+gy*(i==2)
        s1 <- sort(g, index.return = TRUE)[[1]]
        s2 <- sort(g, index.return = TRUE)[[2]]
	n0 <- n - sum(sign(s1))
        h <- round(boxsize/2*sum(sign(s1))+eps(1))
	k <- 1
	while (k <= n){
	    s <- 0
	    while ( (n >= k+s+1) && (s1[k+s+1] == s1[k]) ) {
	    s <- s+1
	    }
	    if (s >= h){
	        upper[s2[k:(k+s)]] <- g[s2[k]]
                lower[s2[k:(k+s)]] <- g[s2[k]]
	    } else {	    	    
                upper[s2[k:(k+s)]] = g[s2[min(n,k+s+h)]]
                lower[s2[k:(k+s)]] = g[s2[max(n0*(n0>h)+1,k-h)]]
	    }	
	    k <- k+s+1
	}

	B[[i]] <- (do.call(cbind, lapply(seq_len(n), function(i) g <= upper[i]))) & 
	          (do.call(cbind, lapply(seq_len(n), function(i) g >= lower[i])))
	a[i,] <- colSums(B[[i]])
}

# Calculate the normalized statistic of edge gx-gy
    res <- (colSums(B[[1]] & B[[2]])*n-a[1,]*a[2,])/sqrt(a[1,]*a[2,]*(n-a[1,])*(n-a[2,])/(n-1)+eps(1))

    return(res)

 }

## Function for computing the average expression values of duplicate genes
# Exp_scRNA: Gene expression values of lncRNAs or mRNAs in single cells, rows are cells and columns are lncRNAs or mRNAs
# Output: temp is single cell expression data without duplicate genes
Averg_Duplicate <- function(Exp_scRNA){
    
    uniqueNameList <- unique(colnames(Exp_scRNA))
    noOfgenes <- length(uniqueNameList)
    temp <- matrix(0, nrow = nrow(Exp_scRNA), ncol = noOfgenes)
    colnames(temp) <- uniqueNameList
    rownames(temp) <- rownames(Exp_scRNA)
    for(c in 1:noOfgenes){
        GeneList <- which(colnames(Exp_scRNA) == colnames(temp)[c])
    for(r in 1:nrow(temp)) {
        temp[r, c] <- mean(as.numeric(Exp_scRNA[r, GeneList]))  
  }
}
    return(temp)
}


## Function for discovering cell-specific lncRNA-mRNA regulatory network
# lncR and mR: Gene expression values of lncRNAs and mRNAs in single cells, rows are cells and columns are lncRNAs or mRNAs
# boxsize: Size of neighborhood (0.1 in default)
# p.value.cutoff: Significance p-value for identifying cell-specific lncRNA-mRNA regulatory network
# Output: res_list is a list of cell-specific lncRNA-mRNA regulatory network
CSlncR_net <- function(lncR, mR, boxsize = 0.1, p.value.cutoff = 0.05) {
  
  lncRs_num <- ncol(lncR)
  mRs_num <- ncol(mR)
  cell_num <- nrow(lncR)
  res <- matrix(NA, nrow = lncRs_num*mRs_num, ncol = cell_num + 2)
  for (i in seq(lncRs_num)){
    for (j in seq(mRs_num)){
      res[(i-1)*mRs_num+j, 1] <- colnames(lncR)[i]
      res[(i-1)*mRs_num+j, 2] <- colnames(mR)[j]
      res[(i-1)*mRs_num+j, 3:(cell_num + 2)] <- csn_edge(lncR[, i], mR[, j], boxsize = boxsize)                                                              
    }
  }
  
  q <- -qnorm(p.value.cutoff)
  
  res_list <- lapply(seq(cell_num), function(i) res[which(as.numeric(res[, i+2]) > q), seq(2)])
  
  return(res_list)
}


## Function for identifying the overlap of cell-type-specific lncRNA-mRNA regulatory network
# Netlist: list object, a list of cell-type-specific lncRNA-mRNA regulatory network
# Intersect_num??The least number of different cell types intersected for overlap.
# The value of 1 means the union of cell-type-specific lncRNA-mRNA interactions from different cell types.
# Output: Overlap_res is the overlap of cell-type-specific lncRNA-mRNA regulatory network in all cell types
Overlap_net_interaction <- function(Netlist, Intersect_num) {
  
  if (length(Netlist) >= 2 & length(Netlist) >= Intersect_num) {
    
    tmp = NULL
    for (i in seq(Netlist)){
      Interin <- Netlist[[i]]
      Interin_paste <- paste(Interin[, 1], Interin[, 2], sep = "_")
      tmp <- c(tmp,Interin_paste)
    }

  Interin_table <- table(tmp)    
  if(Intersect_num == 1 | Intersect_num == 2 ){
    Interin_names <- names(Interin_table)[which(Interin_table ==  Intersect_num)]
  }else{
    Interin_names <- names(Interin_table)[which(Interin_table >=  Intersect_num)]
  }
  
  from <- c()
  to <- c()
  tmpp  <- data.frame(from,to)
  
  for (i in seq(length(Interin_names))){
    tmpp[i,"from"]= strsplit(Interin_names[i],split = "_")[[1]][1]
    tmpp[i,"to"]= strsplit(Interin_names[i],split = "_")[[1]][2]
  }
  
  Overlap_res <- tmpp
  return(Overlap_res)
  
  } else {
  stop("Please check your input!\n")}
}
  
 
## Function for identifying the overlap of cell-type-specific hub lncRNAs
# Netlist: list object, a list of cell-type-specific hub lncRNAs
# Intersect_num??The least number of different cell types intersected for overlap.
# The value of 1 means the union of cell-type-specific hub lncRNAs from different cell types.
# Output: Overlap_res is the overlap of cell-type-specific hub lncRNAs in all cell types
Overlap_hub <- function(hublist, Intersect_num) {
  
  if (length(hublist) >= 2 & length(hublist) >= Intersect_num) {
    
    tmp = NULL
    for (i in seq(hublist)){
      Interin <- hublist[[i]]
      for (j in seq(Interin)) {
        Interin_paste <- Interin[j]
        tmp <- c(tmp,Interin_paste)
      }
    }
    Interin_table <- table(tmp)
    if(Intersect_num == 1 | Intersect_num == 2 ){
      Interin_names <- names(Interin_table)[which(Interin_table ==  Intersect_num)]
    }else{
      Interin_names <- names(Interin_table)[which(Interin_table >=  Intersect_num)]
    }

    Overlap_res <- Interin_names
    return(Overlap_res)
    
  } else {
    stop("Please check your input!\n")
  }
  
}


## Function for calculating similarity matrix between two list of networks
# net1: List object, the first list of network
# net2: List object, the second list of network
# directed: Logical value, network directed (TRUE) or undirected (FALSE)
# Output: Sim is a similarity matrix between two list of networks
Sim.network <- function(net1, net2, directed = TRUE){

    if(class(net1)!="list" | class(net2)!="list") {
    stop("Please check your input network! The input network should be list object! \n")
    }

    m <- length(net1)
    n <- length(net2)
    Sim <- matrix(NA, m, n)
    for (i in seq(m)){
        for (j in seq(n)){
	    net1_graph_interin <- make_graph(c(t(net1[[i]][, 1:2])), directed = directed)
            net2_graph_interin <- make_graph(c(t(net2[[j]][, 1:2])), directed = directed)
	    overlap_interin <- nrow(as_data_frame(net1_graph_interin %s% net2_graph_interin))
	    Sim[i, j] <- overlap_interin/min(nrow(net1[[i]]), nrow(net2[[j]]))
	}
    }

    return(Sim)
}

## Function for calculating similarity matrix between two list of hubs
# hub1: List object, the first list of hub
# hub2: List object, the second list of hub
# Output: Sim is a similarity matrix between two list of hubs
Sim.hub <- function(hub1, hub2){

    if(class(hub1)!="list" | class(hub2)!="list") {
    stop("Please check your input hub! The input hub should be list object! \n")
    }

    m <- length(hub1)
    n <- length(hub2)
    Sim <- matrix(NA, m, n)
    for (i in seq(m)){
        for (j in seq(n)){	    
	    overlap_interin <- length(intersect(hub1[[i]], hub2[[j]]))
	    Sim[i, j] <- overlap_interin/min(length(hub1[[i]]), length(hub2[[j]]))
	}
    }

    return(Sim)
}

## ASD enrichment analysis using hypergeometric distribution test
# lncRExp and mRExp: Gene expression values of lncRNAs and mRNAs in single cells, rows are cells and columns are lncRNAs or mRNAs 
# ASDgenes: ASD-related genes (lncRNAs and mRNAs)
# Modulelist: List object, a list of lncRNA-mRNA bicliques or networks
# Output: A list of significance p-values enriched in ASD
module_ASD_EA <- function(lncRExp, mRExp, ASDgenes, Modulelist) {

    ExpData <- cbind(lncRExp, mRExp)      

    B <- ncol(ExpData)
    N <- length(intersect(colnames(ExpData), as.matrix(ASDgenes)))
    M <- unlist(lapply(seq_along(Modulelist), function(i) length(Modulelist[[i]])))
    x <- unlist(lapply(seq_along(Modulelist), function(i) length(intersect(Modulelist[[i]], as.matrix(ASDgenes)))))    
    p.value <- 1 - phyper(x - 1, N, B - N, M)
    
    names(p.value) <- names(Modulelist)
    return(p.value)
}


random.net.validate <- function(net_num, lncRTarget_graph, lncRExp, mRExp, perm = 100) {
  lncRmR.comb <- expand.grid(colnames(lncRExp), colnames(mRExp))
  colnames(lncRmR.comb) <- c("lncRNA", "mRNA") 
  net_mean <- net_num
  
  for (i in seq(length(net_num))) {
    interin <- 0
    for (j in seq(perm)){
      lncRmR.random <- lncRmR.comb[sample(seq(nrow(lncRmR.comb)), net_num[i]),]
      lncRmR.random.graph <-make_graph(c(t(lncRmR.random[, 1:2])), directed = TRUE)
      lncRmR.random.validated <- as_data_frame(lncRmR.random.graph %s% lncRTarget_graph)
      interin <- interin + nrow(lncRmR.random.validated)/nrow(lncRmR.random)
    }
    interin.mean <- interin/perm
    net_mean[i] <- interin.mean
  }
  
  return(net_mean)
}





