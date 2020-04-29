## The Regents of the University of California and The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2018) by the
## Regents of the University of California abd the 
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

# Load any packages used to in our code to interface with GenePattern.
# Note the use of suppressMessages and suppressWarnings here.  The package
# loading process is often noisy on stderr, which will (by default) cause
# GenePattern to flag the job as failing even when nothing went wrong. 
suppressMessages(suppressWarnings(library(getopt)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(monocle3)))

# Print the sessionInfo so that there is a listing of loaded packages, 
# the current version of R, and other environmental information in our
# stdout file.  This can be useful for reproducibility, troubleshooting
# and comparing between runs.
sessionInfo()

is.emptyString=function(a){return (trimws(a)=="")}

# Get the command line arguments.  We'll process these with optparse.
# https://cran.r-project.org/web/packages/optparse/index.html
arguments <- commandArgs(trailingOnly=TRUE)

# Declare an option list for optparse to use in parsing the command line.
option_list <- list(
  # Note: it's not necessary for the names to match here, it's just a convention
  # to keep things consistent.
  make_option("--input.seurat.rds.file", dest="input.seurat.rds.file", default="notARealFilename.impossible"),
  make_option("--output.file", dest="output.file"),
  make_option("--input.10x.file", dest="input.10x.file", default="notARealFilename.impossible"),
  make_option("--root.cells", dest="root.cells"),
  make_option("--root.nodes", dest="root.nodes"),
  make_option("--plot.color.by", dest="plot.colorbys"),  
  make_option("--plot.post.color.by", dest="plot.post.colorbys")
  
)

# Parse the command line arguments with the option list, printing the result
# to give a record as with sessionInfo.
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE, args=arguments)
print(opt)
opts <- opt$options

mat=NULL

if (file.exists(opts$input.seurat.rds.file)){
    filename = opts$input.seurat.rds.file

    if (endsWith(tolower(filename), ".rds")) {
        mat <- readRDS(opts$input.seurat.rds.file)
    } else { 
        # assumes just one object in the file
        env = new.env()
        matName <- load(filename, env, verbose=TRUE)[1]
        mat <- env[[matName]]
    }
   
    ############### Converting Seurat object to Monocle ###############

    # Extract expression data (raw UMI counts) from the Seurat object
    # We are extracting raw UMI counts since Monocle 3 is designed
    # for use with absolute transcript counts (e.g. from UMI experiments)
    data <- GetAssayData(object = mat, slot = "counts")

    # Extract phenotype data (cell metadata) from the Seurat object
    pData <- mat@meta.data

    # Extract feature data from the Seurat object
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

    # Construct a Monocle CDS (cell_data_set) object
    cds <- new_cell_data_set(data, cell_metadata = pData, gene_metadata = fData)

    # Construct and assign the made up partition
    recreate.partition <- c(rep(1, length(cds@colData@rownames)))
    names(recreate.partition) <- cds@colData@rownames
    recreate.partition <- as.factor(recreate.partition)
    cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

    # Assign the cluster info
    list_cluster <- Idents(object = mat)
    names(list_cluster) <- cds@colData@rownames
    cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

    # A space-holder to essentially fill out louvain parameters
    cds@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

    # Assign UMAP coordinate
    # cds@reducedDims@listData[["UMAP"]] <-mat@reductions[["umap"]]@cell.embeddings
    cds@reduce_dim_aux@listData[["UMAP"]] <-mat@reductions[["umap"]]@cell.embeddings
    cds@reduce_dim_aux@listData[["UMAP"]] <-mat@reductions[["umap"]]@cell.embeddings

    # Assign feature loading for downstream module analysis
    cds@preprocess_aux$gene_loadings <- mat@reductions[["pca"]]@feature.loadings

    # A space-holder needed for order_cells() function later on
    rownames(cds@principal_graph_aux[['UMAP']]$dp_mst) <- NULL
    colnames(cds@reduce_dim_aux@listData[["UMAP"]]) <- NULL
    cds@int_colData@listData$reducedDims@listData[["UMAP"]] <-mat@reductions[["umap"]]@cell.embeddings

    #  this other ways for colnames are for different versions of the packages
    #colnames(cds@reducedDims$UMAP) <- NULL

} else if (file.exists(opts$input.10x.file)){
    # instead of a seurat object, we are getting a zipped 10x genomics file
    
    # TODO:  Check if outs dir is duplicated and move files if necessary
    untar(opts$input.10x.file, exdir="./10x_data/outs")

    cds <- load_cellranger_data("./10x_data")    
    cds <- preprocess_cds(cds, num_dim = 100)
    cds <- align_cds(cds)
    cds <- reduce_dimension(cds)
    cds <- cluster_cells(cds) 

}

############### Trajectory analysis ###############

# Learn the trajectory graph
cds <- learn_graph(cds)

colorbys =  unlist(strsplit(opts$plot.colorbys, ","))
for (colorby in colorbys){
    col = trimws(colorby)
    pdf(paste(opts$output.file, "_", col, ".pdf", sep=""))
    print(plot_cells(cds, color_cells_by = col))
    dev.off()
}

# Order the cells in pseudotime, must have either root.cells or root.nodes specified
if  ( (is.null( opts$root.cells ) && is.null( opts$root.nodes ) ) )  { 
	 stop("No root cells or nodes specified as the start point for ordering. Skipping pseudotime ordering" )
} else {

    if ( ! is.null(opts$root.cells) ){

        cellList = unlist(strsplit(opts$root.cells,","))

        cds <- order_cells(cds, root_cell=cellList)

        
    } else if (!is.null(opts$root.nodes)){ 
        print(opts$root.nodes)
        nodeList = unlist(strsplit(opts$root.nodes,","))
        print(nodeList)
        print("A")
        cds <- order_cells(cds, root_pr_nodes=c(11))

    }

    post.colorbys =  unlist(strsplit(opts$plot.post.colorbys, ","))
    for (colorby in post.colorbys){
        col = trimws(colorby)
        pdf(paste(opts$output.file, "_", col, "_ordered.pdf", sep=""))
        print(plot_cells(cds, color_cells_by = col))
        dev.off()
    }

    # always include the pseudotime ordered
    pdf(paste(opts$output.file, "_pseudotime_ordered.pdf", sep=""))
    print(plot_cells(cds, color_cells_by = "pseudotime"))
    dev.off()

}
save(cds, file = paste(opts$output.file, "_", "monocle.robj", sep=""))
print("Done")




