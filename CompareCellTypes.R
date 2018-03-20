library(ape)
library(phytools)
library(phangorn)
library(ggtree)

# Function to simulate ZINB counts
simulate_counts_from_ZINB <- function(nc,ng){
  ## set seed for mu, r , and d per gene ##
  mu_vec <- 10^(runif(ng,0,4))
  r_vec <- rgamma(n = ng,shape = 1,rate = 0.5)
  d_vec <- runif(ng,0,0.9)
  
  ## build simulation row by row using gene specific values of mu r and d ##
  counts <- t(sapply(1:ng,function(g){
    
    # Generate counts from negative binomial #
    gNB <- rnbinom(nc, size = r_vec[g], mu = mu_vec[g])
    
    # Generate technical dropouts by simulating from binomial; 1=dropout #
    gDO <- rbinom(nc, size = 1, prob = d_vec[g])
    
    # Generate zero-inflated counts #
    gcount <- gNB * (1-gDO)
    return(gcount)
  }))
  
  row.names(counts) <- paste('gene',1:ng,sep='')
  return(counts)
}

## Fitting ZINB ##

calculate_summary_values <- function(counts) {
  if (sum(counts < 0) >0) {stop("Expression matrix contains negative values! Please provide raw UMI counts!")}
  #if ( sum(counts >= 1) != sum(counts > 0) ) {stop("Error: Expression matrix is not integers! Please provide raw UMI counts.")}
  #        if (sum(!is.integer(counts)) >0) {stop("Expression matrix is not integers! Please provide a matrix (not data.frame) raw UMI counts!")}
  
  tjs <- rowSums(counts, na.rm=T) # Total molecules/gene
#  if (sum(tjs <= 0) > 0) {stop("Error: all genes must have at least one detected molecule.")}
  tis <- colSums(counts, na.rm=T) # Total molecules/cell
  if (sum(tis <= 0) > 0) {stop("Error: all cells must have at least one detected molecule.")}
  djs <- ncol(counts)-rowSums(counts > 0, na.rm=T) # Observed Dropouts per gene
  dis <- nrow(counts)-colSums(counts > 0, na.rm=T) # Observed Dropouts per cell
  nc <- length(counts[1,]) # Number of cells
  ng <- length(counts[,1]) # Number of genes
  total <- sum(tis, na.rm=T) # Total molecules sampled
  return(list(tis=tis, tjs=tjs, dis=dis, djs=djs, total=total,nc=nc,ng=ng));
}

fit_zero_inflated_negative_binomial <- function(obs, g_row, vals, e=0.00001) {
  if (sum(vals$tjs[g_row]) <= 0) {return(list(mu=NA, r=NA, d=NA))}
  l_is <- vals$tis*vals$nc/vals$total # cell-specific weights due to difference in library size
  d0 <- vals$djs[g_row]/vals$nc # dropout-rate, initially all zeros due to dropouts
  max_r <- 10^10 # r = size parameter for R's nbinom
  
  d_curr <- d0; 
  d_prev <- -100;
  while( abs(d_curr-d_prev) > e ) {
    
    # Fit the NB #
    mu_j <- vals$tjs[g_row]/(vals$nc-d_curr*vals$nc) # dropouts only affect zeros to just change the denominator of the mean
    mu_ijs <- mu_j*l_is; #account for different library sizes
    weights <- rep(1, times=length(mu_ijs)) # rather than remove some columns 
    
    weights[obs == 0] <- (1-d_curr*vals$nc/vals$djs[g_row]) # lets down-weight the contribution of 
    # zeros for fitting the NB
    # Note: error = n*variance 
    # weight-adjusted observed error correcting for cell-specific library sizes
    #obs_err <- (obs - mu_ijs)^2*weights
    obs_err <- sum((obs - mu_ijs)^2*weights)
    
    # fit the dispersion as:
    # Observed error = sum of variances of cell-specific NB distributions
    # variance of NB = mean+mean^2/dispersion
    # rearrange to solve for dispersion (r)
    r_j <- sum( mu_ijs^2*weights )/(obs_err - sum(mu_ijs*weights))
    if (r_j <= 0) {r_j <- max_r} # BUG
    
    # probability of a zero in each cell given the NB
    p_ijs <- (1 + mu_ijs/r_j)^(-r_j)
    d_exp <- mean(p_ijs)
    
    # Update dropout rate (d) estimate
    d_prev <- d_curr
    d_curr <- (vals$djs[g_row] - d_exp*vals$nc)/vals$nc
    if (d_curr <= 0) {d_curr <- d_prev}
  }
  return(list(mu=mu_j, r=r_j, d=d_prev,gene=g_row));
}

fit_ZINB_to_Matrix <- function(counts) {
  # passes count matrix through fitting fuction
  vals <- calculate_summary_values(counts);
  estim <- lapply(1:vals$ng, function(g) {
    fit_zero_inflated_negative_binomial(counts[g,], g, vals)
  })
  mu <- sapply(estim,function(x){unlist(x)[1]})
  r <- sapply(estim,function(x){unlist(x)[2]})
  d <- sapply(estim,function(x){unlist(x)[3]})
  
  # return data frame with the mu, r, d for each gene 
  return(data.frame(gene=row.names(counts),mu=mu,r=r,d=d))
}


## Fit ZINB to each gene in a cell type
fit_clusters <- function(counts, cluster_assignments) {
  clusters <- sort(unique(cluster_assignments))
  fit.list <- lapply(1:length(clusters), function(x){
    cluster.df <- counts[,which(cluster_assignments == clusters[x])]
    return(fit_ZINB_to_Matrix(cluster.df))
  })
  return(fit.list)
}

# Calculate KL between genes

calculateKL <- function(c1,c2){
  # average over all kl distances
  kl.av <- mean(unlist(sapply(1:nrow(c1), function(g){
    # Only include rows where both genes are represented in the cluster
    if (is.na(c1[g,'mu']) == FALSE && is.na(c2[g,'mu']) == FALSE){
      
      mu1 <- c1[g,'mu']
      r1 <- c1[g,'r']
      sigma1 <- sqrt(mu1 + (mu1^2/r1))
      
      mu2 <- c2[g,'mu']
      r2 <- c2[g,'r']
      sigma2 <- sqrt(mu2 + (mu2^2/r2))
      
      # calculate kl distance between 2 univariate gaussian distributions
      # https://tgmstat.wordpress.com/2013/07/10/kullback-leibler-divergence/ equation [4]
      kl <- log(sigma2/sigma1) + (sigma1^2 + (mu1-mu2)^2)/(2*(sigma2^2)) - 1/2
    }
  })))
  return(kl.av)
}

# Generate a distance matrix using caluclated KL values
computeClusterMatrix <- function(counts, cluster_assignments) {
  fits <- fit_clusters(counts,cluster_assignments)
  clusters <- sort(unique(cluster_assignments))
  kl.matrix <- data.frame()
  for (j in 1:length(clusters)){
    for (i in 1:length(clusters)){
      kl.matrix[i,j] <- calculateKL(fits[[i]],fits[[j]])
    }
  }
  colnames(kl.matrix) <- clusters
  max.kl <- max(rowSums(kl.matrix))
  scaled.kl <- matrix(unlist(kl.matrix/max.kl),nrow=length(clusters))
  t.kl <- t(scaled.kl) + scaled.kl
  return(t.kl)
}

# Make a tree topology from distance matrix 
# requires cell assignments in a vector
# tree type options are nj, hclust, upgma, fastme, network
# will save the distance matrix ($dist), as well as the tree items generated by each method

makeTree <- function(counts,cluster_assignments,tree_type){
  distMatrix <- as.dist(computeClusterMatrix(counts,cluster_assignments))
  # Neighbour joining
  nj.tree <- bionj(distMatrix) 
  if (tree_type == 'nj'){
    ggtree(nj.tree) + geom_tiplab(colour='grey20',size=5)
  }
  # Plain hierarchical clustering
  hclust.tree <- hclust(distMatrix) 
  if (tree_type == 'hclust'){
    plot(as.phylo(hclust.tree),type='phylogram',edge.width=2)
  }
  # UPGMA algorithm (type of hclust; branch lengths are halved)
  upgma.tree <- upgma(distMatrix)
  if (tree_type == 'upgma'){
    plot(upgma.tree,type='phylogram',edge.width=2)
  }
  # Balanced fastme (doesn't calculate BL)
  fastme.tree <- fastme.bal(distMatrix)
  if (tree_type == 'fastme'){
    ggtree(midpoint.root(fastme.tree)) + geom_tiplab(colour='grey20',size=5)
  }
  # Network
  if (tree_type == 'network'){
    library(qgraph)
    qgraph(distMatrix, layout='spring', vsize=5)
  }
  return(list(dist=distMatrix,nj=nj.tree,hclust=hclust.tree,upgma=upgma.tree,fastme=fastme.tree))
}

### FOR SIMULATED DATA ###
nc <- 1000
ng <- 1000
cluster_assigns <- sample(1:5,nc,replace=TRUE)
cluster_assigns <- as.numeric(as.character(cluster_assigns)) + 1
sim <- simulate_counts_from_ZINB(nc,ng)
sim.tree <- makeTree(sim,cluster_assigns,'nj')

#### REAL MOSQUITO DATA ####
mosq.data <- readRDS(file = "Mosquito_data/aggr81scat.rds")
mosq.counts <- round(counts(mosq.data)) # make counts round integers

# Save a vector of cell type assignments
mosq.clusters.0.6 <- colData(mosq.data)$res.0.6
mosq.clusters.man <- colData(mosq.data)$manual_clustering

# have to add one to cluster assignments
# or else they crash
# should probably be fixed...
mosq.clusters.0.6 <- as.numeric(as.character(mosq.clusters.0.6)) + 1
mosq.clusters.man <- as.numeric(as.character(mosq.clusters.man)) + 1

# NEIGHBOUR-JOINING TREE, 'refined' clusters
mosq.tree.0.6 <- makeTree(mosq.counts,mosq.clusters.0.6,'nj')
mosq.labels.0.6 <- c('Prohemocytes 1','Inactive',
                     'Mid-Active','Oenocytoids','Active',
                     'Regulatory','Anti-bacterial','Prohemocytes 2',
                     'Rapidly Cycling','Fat')
# NEIGHBOUR-JOINING TREE, manual clusters
mosq.tree.man <- makeTree(mosq.counts,mosq.clusters.man,'nj')
mosq.labels.man <- c('Prohemocytes 1','Prohemocytes 2','Mid-Active',
                     'Oenocytoid 1','Active','Regulatory','Anti-bacterial',
                     'Prohemocytes 3','Rapidly Cycling','Fat body',
                     'Secreting Hemocytes','Inactive','Oenocytoid 2',
                     'Light-sensors')

# Write trees in newick format
mosq.tree.0.6$nj$tip.label <- mosq.labels.0.6
write.tree(mosq.tree.0.6$nj, file = "Mosq0.6.bionj.nwk.tre", append = FALSE,
           digits = 10, tree.names = FALSE)
mosq.tree.man$nj$tip.label <- mosq.labels.man
write.tree(mosq.tree.man$nj, file = "MosqManualCluster.bionj.nwk.tre", append = FALSE,
           digits = 10, tree.names = FALSE)


# Plot trees in R - note, ggtree crashing RStudio for some reason??
# Best option is to save a nwk, look into why plot function has started crashing?
mosq.tree.0.6$nj$tip.label <- mosq.labels.0.6
t1 <- ggtree(mosq.tree.0.6$nj) + geom_tiplab(colour='grey20',size=5)
mosq.tree.0.6$fastme$tip.label <- mosq.labels.0.6
t2 <- ggtree(mosq.tree.0.6$fastme) + geom_tiplab(colour='grey20',size=5)
ggsave(filename = 'testtree.bionj.pdf',plot = t1,device = 'pdf')

mosq.tree.man$nj$tip.label <- mosq.labels.man
plot(mosq.tree.man$nj)
t1 <- ggtree(mosq.tree.man$nj) + geom_tiplab(colour='grey20',size=5)
mosq.tree.man$fastme$tip.label <- mosq.labels.man
t2 <- ggtree(mosq.tree.man$fastme) + geom_tiplab(colour='grey20',size=5)
ggsave(filename = 'testtree.bionj.pdf',plot = t1,device = 'pdf')
