
## This is a function that returns modularity-type scoring for individual clusters in a graph.
### Clusters with a higher ratio of edge weights inside the cluster vs connecting to other
### clusters will recieve a higher score. Clusters are scored on the scale of 0(worst) to 1(best).

### This function requires, at minimum, an igraph object with an edge attribute of 'pcor'.

## This function generates clusters through walktrap and cutreeDynamic. 
### Use parameter 't' to set walk length and 'minclusts' to set minimum cluster size.
### Previously generated clusters (through any method, not just walktrap) can be imputed as 
### a membership vector through the parameter 'clusters'.

## If you want to calculate cluster size, set 'sizes = TRUE'

## If you want to calculate p.value, set 'calculate.p = TRUE'
### This function uses a monte carlo simulation to calculate p values. If you already generated a null 
### set, pass it to the function through the 'nulldist' parameter. If you have not, the function will
### generate 1000 random clusters.

## This executes a system command
## Note that the environment path is incorrect and needs to be fixed for commands to work
## nbconvert needs this code

Sys.setenv(PATH = paste("/usr/lib64/qt-3.3/bin:/opt/applications/cytoscape/",
                        "3.3.0:/opt/applications/R/3.5.1/gnu/bin:/opt/applications/",
                        "gcc/4.9.4/bin:/opt/applications/java/jdk1.8.0_65/bin:/opt/",
                        "applications/ant/apache-ant-1.9.0/bin:/opt/applications/",
                        "python/3.6.3/gnu/bin:/usr/local/go/bin:/usr/local/sbin:/",
                        "usr/local/bin:/sbin:/bin:/usr/sbin:/usr/bin:/usr/local/",
                        "bin:/usr/lpp/mmfs/bin:/root/bin:/usr/local/bin:/usr/lpp/mmfs/bin",
                       sep = ''))

osc <- function (cmd) {
    scra <- '/gpfs/group/torkamani/devans/GTEx/Scratch.Area/' # Location of temporary files
    fullcmd <- paste(cmd, ' >& ', scra, 'temp.txt', sep = '')
    system(fullcmd)
    lines <- readLines(paste(scra, 'temp.txt', sep = ''))
    system(paste('rm ', scra, 'temp.txt', sep = ''))
    return(lines)
    }

cluster_score <- function(net, clusters = NULL, sizes = FALSE, max.str = FALSE, t = 10, minclusts = 30, calculate.p = FALSE, nulldist = NULL){
    
    #load packages
    library(igraph)
    library(corpcor)
    library(WGCNA)
    
    # generate clusters
    if(is.null(clusters)){
        wt <- cluster_walktrap(net, weights = E(net)$pcor, steps = t,
                            merges = TRUE, modularity = TRUE, membership = TRUE)
        hclust <- as.hclust(wt)
        dyntree <- cutreeDynamic(hclust, method = 'tree', minClusterSize = minclusts)
        clusters <- dyntree
    }
    
    # generate null distribution
    if(calculate.p & is.null(nulldist)){
        randperms <- vector(mode = 'numeric', length = 1000)
        for(i in 1:1000) {
            samplenamesrand <- sample(V(net)$name, size = minclusts)
            just1clustrand <- induced_subgraph(net, V(net)$name %in% samplenamesrand)
            clustedgesrand <- E(net)[E(net) %in% unlist(incident_edges(net, V(net)$name %in% samplenamesrand))]
            internalsumrand <- sum(abs(E(just1clustrand)$pcor))
            totalsumrand <- sum(abs(clustedgesrand$pcor))
            clusterscorerand <- internalsumrand/totalsumrand
            randperms[i] <- clusterscorerand
            nulldist <- randperms
        }
    }
    
    # generate cluster score 
    V(net)$clusters <- clusters
    scoresvec <- vector(mode = 'numeric', length = length(unique(clusters)))
    names(scoresvec) <- unique(clusters)[order(unique(clusters))]
    count <- 1
    for(j in unique(clusters)[order(unique(clusters))]) {
        clustnum <- j
        just1clust <- induced_subgraph(net, V(net)$clusters == j)
        clustedges <- E(net)[E(net) %in% unlist(incident_edges(net, V(net)$clusters ==j))]
        internalsum <- sum(abs(E(just1clust)$pcor))
        totalsum <- sum(abs(clustedges$pcor))
        clusterscore <- internalsum/totalsum
        scoresvec[count] <- clusterscore
        count = count+1
        
    }
    colnames <- 'score'
    returnthis <- as.data.frame(scoresvec)
    rownames(returnthis) <- names(scoresvec)
    colnames(returnthis) <- colnames
    
    # cluster size
    if(sizes){
        clustsizes <- as.data.frame((table(V(net)$clusters)))[,2]
        returnthis <- cbind(returnthis, clustsizes)
        colnames <- c(colnames, 'size')
        colnames(returnthis) <- colnames
    }
   
    # calculate pvals
    if(calculate.p){
        pvals <- vector(mode = 'numeric', length = length(scoresvec))
        names(pvals) <- names(scoresvec)
        for (k in length(scoresvec)){
            realscore <- scoresvec[k]
            nullishigher <- nulldist[nulldist > realscore]
            pvalue <- length(nullishigher)/length(nulldist)
            pvals[k] <- pvalue
        }
        colnames <- c(colnames, 'pvals')
        returnthis <- (cbind(returnthis, pvals))
        colnames(returnthis) <- colnames
    }
    
    # max strength
    if(max.str) {
        strength <- vector(mode = 'numeric', length = length(scoresvec))
        names(strength) <- names(scoresvec)
        count <- 1
        for(l in unique(clusters)[order(unique(clusters))]){
            strength[count] <- max(V(net)[V(net)$clusters == l]$strength)
            count <- count + 1
        }
        colnames <- c(colnames, 'max.str')
        returnthis <- (cbind(returnthis, strength))
        colnames(returnthis) <- colnames
    }
    # return a data frame
    return(returnthis)
}

cluster_score.help <- function(){
    cat("THANK YOU FOR CALLING THE HELP DESK :)")
    cat("# This is a function that returns modularity-type scoring for individual clusters in a graph.

    ## Clusters with a higher ratio of edge weights inside the cluster vs connecting to other
    ## clusters will recieve a higher score. Clusters are scored on the scale of 0(worst) to 1(best).

    # This function requires, at minimum, an igraph object with an edge attribute of 'pcor'.

    # This function generates clusters through walktrap and cutreeDynamic. 
    ## Use parameter 't' to set walk length and 'minclusts' to set minimum cluster size.
    ## Previously generated clusters (through any method, not just walktrap) can be imputed as 
    ## a membership vector through the parameter 'clusters'.

    # If you want to calculate cluster size, set 'sizes = TRUE'

    # If you want to calculate p.value, set 'calculate.p = TRUE'
    ## This function uses a monte carlo simulation to calculate p values. If you already generated a null 
    ## set, pass it to the function through the 'nulldist' parameter. If you have not, the function will
    ## generate 1000 random clusters.")
}

## This cell prevents the nbconvert code below from excuting when this module is loaded using "source" function
FirstRun <- TRUE

## Execute the next cell before running this cell to convert this file to an R script
if (!FirstRun) print(osc("jupyter nbconvert cluster_score.ipynb --to script"))

## Execute this cell, and then the cell above to convert this file to an R script
FirstRun <- FALSE
