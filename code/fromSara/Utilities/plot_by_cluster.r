# This is a function that plots only a subset of clusters from a larger igraph object

# it requires: an igraph object (graph), a membership vector (clusts), 
## and a table of scores generated by the cluster_score function (scores)
## when creating the scores table, be sure max.str is set to TRUE

# the user must specify the method of specifying clusters (mode)
## 'number' plots the top scoring clusters. The number of clusters to plot is set with (numclusts)
## 'cutoff' plots clusters scoring above the score set in (minscore)
## 'specify' plots given clusters irrespective of score. Clusters to plot are set in (specifyclusts)

# set a file path for the graph to ouptut too in (filepath)

# (crosses) is still a work in progress, do not use this variable

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

plot_by_cluster <- function(graph, clusts, scores, mode = c('number','cutoff','specify'), numclusts = 5, minscore = .3, specifyclusts = c(1,3,5,9), crosses = FALSE, filename = '/gpfs/group/torkamani/devans/GTEx/Data.Work.Area.2/xTissue/Unorganized/plot_by_clusts.pdf'){

    # get the subgraph, method depends on mode user selects
    if(identical(mode, 'number')){
        orderscore <- scores[order(scores$score, decreasing = TRUE),]
        topclusts <- as.numeric(sort(rownames(orderscore[1:numclusts,])))
        selectednames <- names(clusts[clusts %in% topclusts])
        subgraph <- induced_subgraph(graph, V(graph)$name %in% selectednames)
    }
    if(identical(mode, 'cutoff')){
        orderscore <- scores[order(scores$score, decreasing = TRUE),]
        topclusts <- as.numeric(sort(rownames(orderscore[orderscore$score >= minscore,])))
        selectednames <- names(clusts[clusts %in% topclusts])
        subgraph <- induced_subgraph(graph, V(graph)$name %in% selectednames)
    }
    if(identical(mode, 'specify')){
        orderscore <- scores[order(scores$score, decreasing = TRUE),]
        wantedclusts <- as.numeric(sort(specifyclusts))
        selectednames <- names(clusts[clusts %in% wantedclusts])
        subgraph <- induced_subgraph(graph, V(graph)$name %in% selectednames)
    }
 
    # generate comunities object
    rlcommunities <- create.communities(subgraph, membership = clusts[names(clusts) %in% selectednames])

    lt <- layout_with_dh
    clust_temp <- rlcommunities
    graph_temp <- delete.edges(subgraph, E(subgraph)[crossing(clust_temp, subgraph)])
    set.seed(103)
    graph_lo <- layout_components(graph_temp, layout = lt)
    graph_lo <- layout.norm(graph_lo, -1, 1, -1, 1)
    
    # Set up node color based on strength (partially transluscent)
    mxscale <- max(V(graph_temp)$strength) ##
    cl <- V(graph_temp)$strength / mxscale
    # V(subgraph)$color <- rgb(1-cl, 1-cl, 1, alpha = .75)
    # V(graph_temp)$color <- rgb(1-cl, 1-cl, 1, alpha = .75)

    # Set up more plotting parameters
    # Node border color 
    pcnodes <- which(V(graph_temp)$v27type == 'protein_coding')
    npcnodes <- which(V(graph_temp)$v27type != 'protein_coding')
    pccol <- rgb(1, 0, 0, alpha = .75)
    npccol <- rgb(0, 1, 0, alpha = .75)
    node_frames_cols <- vector(mode = "character", length = length(pcnodes))
    node_frames_cols[pcnodes] <- pccol
    node_frames_cols[npcnodes] <- npccol
    # Module label color (translucent blue)
    nodelabc <- rgb(0, 0, 1, alpha = .4)
    # Module background color (make it clear)
    polycol <- rgb(0, 1, 1, alpha = .0)
    
    ## Fixing an issue?
    node_frames_cols <- node_frames_cols[node_frames_cols!='']    
        
    # fixing names because igraph is cranky when clust nums aren't sequential
    clust_temp2 <- clust_temp
    numcl <- length(unique(membership(clust_temp)))
    clustname <- sort(unique(membership(clust_temp)))
    for (i in 1:numcl){
        clust_temp2$membership[clust_temp2$membership == clustname[i]] <- i
    }

    # cluster labels
    ctxtpos <- matrix(0, nrow = numcl, ncol = 2)
    colnames(ctxtpos) <- c('x', 'y')
    for (i in 1:numcl) {
        ctxtpos[i,'x'] <- max(graph_lo[membership(clust_temp2) == i, 1])
        ctxtpos[i,'y'] <- mean(graph_lo[membership(clust_temp2) == i, 2])
    }

    ## a fix for label size
    V(graph_temp)$label.cex = .1
    E(graph_temp)$label.cex = .1
    
    ## will be used later
    if(crosses){
        plotvar <- subgraph
    } else {
        plotvar <- graph_temp
    }

    # Set up plot file
    # Plot with legend and module ids
    pdf(filename, width = 30, height = 30)
    plot(plotvar, layout = graph_lo, mark.groups = clust_temp2, mark.shape = 1, mark.expand = 1.5, mark.col = polycol,
        vertex.frame.color = node_frames_cols, vertex.shape = 'circle',
        edge.width = 20 * abs(E(subgraph)$pcor), vertex.size = .6,
        main = "Partial Correlation Network\nRecursive Louvain")
    legend(x = -0.85, y = -.85, legend = c('Protein Coding biotype genes have Red Frames',
                                        'All Other biotype genes have Green Frames',
                                        'Nodes Shaded Blue Are Not Effected By Disease State',
                                        'Red Is Healthy, Green Is Diseased',
                                        'Strength Shown Below Gene Name in Each Node'), col = 'black')
    
    text(x = ctxtpos[,1], y = ctxtpos[,2], adj = c(0, 0), # adj = c(-.25, .5),
         paste('cluster = ', clustname, sep = ''),
         cex = .5, col = nodelabc, font = 4)
    
    text(x = ctxtpos[,1], y = ctxtpos[,2] - .006, adj = c(0, 0), # adj = c(-.25, .5),
         paste('score = ', round(scores[as.character(clustname),]$score, digits = 2), sep = ''),
         cex = .5, col = nodelabc, font = 4)
    
    text(x = ctxtpos[,1], y = ctxtpos[,2] - .012, adj = c(0, 0), # adj = c(-.25, .5),
         paste('max strength = ', round(scores[as.character(clustname),]$max.str, digits = 2), sep = ''),
         cex = .5, col = nodelabc, font = 4)
    
    dev.off()
}

plot_by_cluster.help <- function(){
    cat("THANK YOU FOR CALLING THE HELP DESK :)")
    cat("# This is a function that plots only a subset of clusters from a larger igraph object

    # it requires: an igraph object (graph), a membership vector (clusts),
    ## and a table of scores generated by the cluster_score function (scores)
    ## when creating the scores table, be sure max.str is set to TRUE

    # the user must specify the method of specifying clusters (mode)
    ## 'number' plots the top scoring clusters. The number of clusters to plot is set with (numclusts)
    ## 'cutoff' plots clusters scoring above the score set in (minscore)
    ## 'specify' plots given clusters irrespective of score. Clusters to plot are set in (specifyclusts)

    # set a file path for the graph to ouptut to in (filepath)

    # (crosses) is still a work in progress, do not use this variable")
}

## This cell prevents the nbconvert code below from excuting when this module is loaded using "source" function
FirstRun <- TRUE

## Execute the next cell before running this cell to convert this file to an R script
if (!FirstRun) print(osc("jupyter nbconvert plot_by_cluster.ipynb --to script"))

## Execute this cell, and then the cell above to convert this file to an R script
FirstRun <- FALSE


