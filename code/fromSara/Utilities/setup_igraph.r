
## This is a library of functions that can be used
## to support network related activities, especially
## related to plotting and tabulations

## tested in 48.demo

## The functions defined in this notebook are:
##
## stat_parm <- function () # apply expression stats and strength to network
##
## annot_parm <- function () # apply biomaRt annotations to a network
##
## v19_parm <- function () # apply GENCODE v19 chromosomal information to a network
##
## plot_parm <- function () # apply standard plotting parameters to a network
## make_v19chrom <- function () # Create a file with just gene info from GENCODE 19 input file
##
## Notes added 5/22/2019
## Need to revise this to include the following
## -- Expand Network to Include Many More Biological Attributes, Add Global Network Details
## -- Provide Details on Inputs, Outputs, and Usage for Above Functions
## -- Create Examples of Simulated Data to Demonstrate Methodology



## This executes a system command
## Note that the environment path is incorrect and needs to be fixed for commands to work
## Two functions in this notebook need to perform command line functions, hence this code.

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

## Quiet version
qosc <- function (cmd) {
    scra <- '/gpfs/group/torkamani/devans/GTEx/Scratch.Area/' # Location of temporary files
    system(cmd)
    return(NULL)
    }

## This cell prevents the nbconvert code below from excuting when this module is loaded using "source" function
FirstRun <- TRUE

## Load packages needed by functions in this notebook
## load igraph
suppressMessages(library(igraph))
## load biomaRt
suppressMessages(library(biomaRt))

stat_parm <- function(gf, log10expdata){
    
    ## Remove expression values not in the network, get antilog
    expdata_filt <- 10 ^ log10expdata[,V(gf)$name]
     
    ## Compute the node expression stats, including coeff of variation, add to network
    V(gf)$mean <- apply(expdata_filt, 2, mean)
    V(gf)$sd <- apply(expdata_filt, 2, sd)
    V(gf)$cv <- V(gf)$sd/V(gf)$mean
    V(gf)$med <- apply(expdata_filt, 2, median)

    ## Compute the strength
    V(gf)$strength <- strength(gf)
    return(gf)
    }

annot_parm <- function(gf){
    
    ## Get the biomart object
    ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

    ## Get the annotation for this network
    
    ## Get the network list of genes to annotate
    r <- 1:length(V(gf)$gene)
    gene_sub <- substr(V(gf)$gene[r],1,15)
    
    ## Perform the biomaRt lookup (only lookup unique gene names)
    bm_annot <- getBM(attributes=c('external_gene_name', 'ensembl_gene_id', 'description',
                               'gene_biotype', 'chromosome_name', 'start_position', 'end_position'), 
                filters = 'ensembl_gene_id', 
                values = unique(gene_sub),
                mart = ensembl,
                uniqueRows = TRUE)
    
    ## Shorten the descriptions by removing source info (could not vectorize!)
    for (i in 1:dim(bm_annot)[1]) {
        item1 <- unlist(strsplit(bm_annot[i,3], split = 'Source'))[1]
        bm_annot$sdesc[i] <- substr(item1, 1, nchar(item1) - 2)
        }
    
    ## Set up row names for easy lookup
    rownames(bm_annot) <- bm_annot[,2]

    ## Put the values in the network
    V(gf)[r]$sdesc      <- bm_annot[gene_sub, 'sdesc']
    V(gf)[r]$ename      <- bm_annot[gene_sub, 'external_gene_name']
    V(gf)[r]$biotype    <- bm_annot[gene_sub, 'gene_biotype']
    V(gf)[r]$chr        <- bm_annot[gene_sub, 'chromosome_name']
    V(gf)[r]$strtpos    <- bm_annot[gene_sub, 'start_position']
    V(gf)[r]$endpos     <- bm_annot[gene_sub, 'end_position']
    
    ##Add attribute to show annotation date-stamp(this code causes a graph object to be undisplayable)
    ## graph_attr(gf, "biomaRt Annot Date") <- format(Sys.Date( ) , format="%B %d %Y")
    
    return(gf)
    }

v19_parm <- function(gf, v19cfile){
    
    ## Get the v19 chromosomal information
    v19chrom <- read.table(file = v19cfile, stringsAsFactors = FALSE, sep = '\t',
                           strip.white = TRUE)
    
    ## Set up the column and row names
    row.names(v19chrom) <- v19chrom[, 'V5']
    colnames(v19chrom) <- c('chr', 'desig', 'start', 'stop', 'ensembl_id')
    
    # return(gf)
    
    ## Add v19 chromosome info
    V(gf)$v19chr <- v19chrom[V(gf)$gene,'chr']
    V(gf)$v19start <- v19chrom[V(gf)$gene,'start']
    V(gf)$v19stop <- v19chrom[V(gf)$gene,'stop']
    
    return(gf)
    }

# Apply standard plotting parameters to a multi-tissue network. Currently 
# can differentiate up to 4 tissues

plot_parm <- function (gf) {

    tissues <- table(V(gf)$tissue) 
    tnames <- names(tissues)
    ntissues <- length(tnames)
    fcolchoice <- c('red', 'blue', 'gray', 'magenta')
    for (tid in 1:ntissues) {
        V(gf)[V(gf)$tissue == tnames[tid]]$color <- fcolchoice[tid]
    }
    
    ## Set the node color legend attribute (this code causes a graph object to be undisplayable)
    # tlegend <- paste(tnames, fcolchoice[1:ntissues], sep = ' = ', collapse = ', ')
    # graph_attr(gf, "tissue color legend") <- tlegend 
    
    E(gf)$color[E(gf)$signedw > 0] <- 'green'
    E(gf)$color[E(gf)$signedw < 0] <- 'orange'
    fedthick <- 100 * abs(E(gf)$pcor)
    V(gf)$label.cex = .1
    E(gf)$label.cex = .1
    V(gf)$label <- paste(substr(V(gf)$v19,1,10),
                                   round(V(gf)$strength,3), sep = '\n')
    E(gf)$width <- fedthick
    E(gf)$label <- round(E(gf)$signedw, 4)
    V(gf)$size <- 1
    
    # E(gf)$color[E(gf)$signedw > 0] <- 'green'
    # E(gf)$color[E(gf)$signedw < 0] <- 'orange'
    # fedthick <- 100 * abs(E(gf)$pcor)
    # V(gf)$label.cex = .1
    # E(gf)$label.cex = .1
    # V(gf)$vertex.label <- paste(substr(V(gf)$v19,1,10),
    #                                round(V(gf)$strength,3), sep = '\n')
    # E(gf)$edge.width <- fedthick
    # E(gf)$edge.label <- round(E(gf))$signedw
    # V(gf)$vertex.size <- 1
    
    return(gf)
    }

## Create a file with extracted v19 chromosomal info in it
make_v19chrom <- function(infile, outfile){

    t9 <- "cat "
    t10  <- " | grep -v '##' | grep -v 'exon' | cut -f1,4  -d ' '"
    t11 <- " | tr '\ ' '\t' | cut -f1,3,4,5,10 | grep 'gene' | tr ';' ' ' > "
    cmd <- paste(t9, infile, t10, t11, outfile, sep = '')
    
    ## Execute the command
    temp <- qosc(cmd) # Execute the command
    
    ## Return value is not really needed
    return(cmd)
}

## Execute the next cell before running this cell to convert this file to an R script
if (!FirstRun) print(osc("jupyter nbconvert setup_igraph.ipynb --to script"))

## Execute this cell, and then the cell above to convert this file to an R script
FirstRun <- FALSE

## removed old function save_net
