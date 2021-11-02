# This function scores the matchup between two groups of clusters.

## 1) it gets cluster set 1's node names
## 2) it loops through cluster set 2 to find the best match ((set1 genes in set 2)/set1)
## 3) prints % overlap and overlapping genes between the best match cluster

# users need to input:
## (clusters1), the membership vector of cluster set 1
## (clusters2), the membership vector of cluster set 2
## (name1), the name of set 1 for the ouptput sentance
## (name2), the name of set 2 for the ouptput sentance

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

cluster_overlap <- function(clusters1, clusters2, name1 = 'Group 1', name2 = 'group2'){
    
    for (i in 1:length(unique(clusters1))){
        nodenames <- V(tgraph)$name[clusters1 == i]
        nodenum <- sum(clusters1 == i)
        bestmatch <- 1             #temporary, just to initialize variable
        whichmatch <- nodenames    #temporary, just to initialize variable
        percentmatch <- 0          #temporary, just to initialize variable
        for (j in 1:length(unique(clusters2))){
            othernames <- V(tgraph)$name[clusters2 == j]
            numinboth <- sum(nodenames %in% othernames)
            percent <- (numinboth/nodenum)*100
            if(percent > percentmatch){
                bestmatch <- j
                percentmatch <- percent
                whichmatch <- nodenames[nodenames %in% othernames]
            }
        }
        print(paste(name1, " cluster ", i , " matches ", name2, " cluster " , bestmatch , ". Percent match is " , percentmatch, sep = ''))
        print("The nodes in both are:")
        print(whichmatch)
    }
}

cluster_overlap.help <- function() {
    cat("THANK YOU FOR CALLING THE HELP DESK :)")
    cat("# This function scores the matchup between two groups of clusters.

    ## 1) it gets cluster set 1's node names
    ## 2) it loops through cluster set 2 to find the best match ((set1 genes in set 2)/set1)
    ## 3) prints % overlap and overlapping genes between the best match cluster

    # users need to input:
    ## (clusters1), the membership vector of cluster set 1
    ## (clusters2), the membership vector of cluster set 2
    ## (name1), the name of set 1 for the ouptput sentance
    ## (name2), the name of set 2 for the ouptput sentance")
}

## This cell prevents the nbconvert code below from excuting when this module is loaded using "source" function
FirstRun <- TRUE

## Execute the next cell before running this cell to convert this file to an R script
if (!FirstRun) print(osc("jupyter nbconvert cluster_overlap.ipynb --to script"))

## Execute this cell, and then the cell above to convert this file to an R script
FirstRun <- FALSE


