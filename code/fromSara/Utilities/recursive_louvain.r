# This function recursively generates louvain clusters.

# the parameter (graph) takes the igraph object to be clustered

# It will continue to break down large clusters by calling
## itself until all of the clusters are smaller than the size
## set in the parameter (biggest)

# returns a membership vector

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

recursive_louvain <- function(graph, biggest = 150, count = 0, membership = NULL) {
    print('Hello')
    
    if(is.null(membership)){
        #Only runs on the first loop
        cl <- cluster_louvain(graph, weights = abs(E(graph)$pcor))
        # generates a vector from membership object
        membervec <- as.vector(membership(cl))
        names(membervec) <- names(membership(cl))
        membership <- membervec
        }
    
    for (i in unique(membership)){
         clustersize <- sum(membership == i)
         if(clustersize > biggest){
             whichtokeep <- membership == i
             smallerclust <- induced_subgraph(graph, V(graph)$name %in% names(membership[whichtokeep]))
             cl <- cluster_louvain(smallerclust, weights = abs(E(smallerclust)$pcor))
             membervec <- as.vector(membership(cl))
             names(membervec) <- names(membership(cl))
             membership[names(membervec)] = (membervec + count*100)
            }
        }
    if(sum(table(membership)>biggest) == 0){
        
        ## Added by Doug to return a communities object instead of a membership set
        selectednames <- names(membership)
        print(selectednames[1:5])
        # rlcommunities <<- create.communities(graph, membership = selectednames)
        rlcommunities <- make_clusters(graph, membership = membership, algorithm = 'multi level')

        return(rlcommunities) # membership)
    } else {
        return(recursive_louvain(graph, biggest = biggest, count = count+1, membership = membership))
    }
}

recursuve_louvain.help <- function() {
    cat("THANK YOU FOR CALLING THE HELP DESK :)")
    cat("# This function recursively generates louvain clusters.

    # the parameter (graph) takes the igraph object to be clustered

    # It will continue to break down large clusters by calling
    ## itself until all of the clusters are smaller than the size
    ## set in the parameter (biggest)

    # returns a membership vector")
}

## This cell prevents the nbconvert code below from excuting when this module is loaded using "source" function
FirstRun <- TRUE

## Execute the next cell before running this cell to convert this file to an R script
if (!FirstRun) print(osc("jupyter nbconvert recursive_louvain.ipynb --to script"))

## Execute this cell, and then the cell above to convert this file to an R script
FirstRun <- FALSE

test_f <- function(biggest = 150) {
    print("hello")
    print(biggest)
    return(biggest * 2)
    }


