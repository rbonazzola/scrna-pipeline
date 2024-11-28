library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(monocle3)
library(SeuratWrappers)
library(ff)
library(argparse)

# setwd("~/Postdoc/scRNA/pipeline_seurat_monocle/")
source(file="utils.R")

# Create an argument parser
parser <- ArgumentParser(description = "Single-cell RNA-seq pipeline with Seurat and Monocle3")

# Add arguments
parser$add_argument("--name_rep", default = "WT_rep123", help = "Name of the replicate set (default: WT_rep123)")
parser$add_argument("--data_dir", default = "data/", help = "Directory containing input data (default: current directory)")
parser$add_argument("--rds_filename", default = NULL, help = "Path to save the RDS file (default: auto-generated)")
parser$add_argument("--overwrite_rds", action = "store_true", help = "Overwrite existing RDS file (default: FALSE)")
parser$add_argument("--output_folder", default="output/", help = "Folder where output is stored")

parser$add_argument("--pdf_width", type = "double", default = 10, help = "Width of the output PDF plots (default: 10)")
parser$add_argument("--pdf_height", type = "double", default = 10, help = "Height of the output PDF plots (default: 10)")

parser$add_argument("--random_seed", type = "integer", default = 42, help = "Random seed for all the stochastic functions.")
parser$add_argument("--marker_genes_file", type = "character", default="", help = "Marker genes for different cell types.")

args <- parser$parse_args()

entropy <- function(count_vector) {
  freq_vector <- count_vector / sum(count_vector)
  ifelse(!is.na(-sum(log2(freq_vector) * freq_vector)), -sum(log2(freq_vector) * freq_vector), 0)
}

# name of the output pdf file and output rds file
NAME_REP      <- args$name_rep
DATA_DIR      <- args$data_dir

RDS_FILENAME <- ifelse(is.null(args$rds_filename), glue::glue("{args$output_folder}/sc_full{NAME_REP}__seed{args$random_seed}.rds"), args$rds_filename)

OVERWRITE_RDS <- args$overwrite_rds

pdf(glue::glue('{args$output_folder}/plots_{NAME_REP}.pdf'), width=args$pdf_width, height=args$pdf_height)

options(
  future.globals.maxSize = 3 * 1024^3
)

list_reps <- list("rep_1", "rep_2", "rep_3")

if (!file.exists(RDS_FILENAME) || OVERWRITE_RDS) {
  
  organize_files_in_folders()
  # import the data
  sc.data_WT_rep1 <- Read10X(data.dir = glue::glue("{DATA_DIR}/WT1"))
  sc.data_WT_rep2 <- Read10X(data.dir = glue::glue("{DATA_DIR}/WT2"))
  sc.data_WT_rep3 <- Read10X(data.dir = glue::glue("{DATA_DIR}/WT3"))
  
  # remember the columns to label each replicate  
  col_1 <- colnames(sc.data_WT_rep1)
  col_2 <- colnames(sc.data_WT_rep2)
  col_3 <- colnames(sc.data_WT_rep3)
  
  # Build sc_data as the concatenation of each dataset  
  sc.data = Reduce("cbind", list(sc.data_WT_rep1, sc.data_WT_rep2,sc.data_WT_rep3))  
  
  #>>> ESTO LO AGREGUE YO 
  # Remove repeated columns (TODO: check why they are repeated)
  sc.data_ = sc.data[, colnames(sc.data) %>% sort %>% .[duplicated(.) %>% !.]]
  
  # Create the Seurat Object  
  sc_full <- CreateSeuratObject(counts = sc.data_, project = "scRNA", min.cells = 3, min.features = 200)
  
  sc_full[['rep_label']] = 0
  # assigning to the Seurat object the replicate label for each cell
  sc_full[['rep_label']][col_1,] = 'rep_1'
  sc_full[['rep_label']][col_2,] = 'rep_2'
  sc_full[['rep_label']][col_3,] = 'rep_3'
  # Extract the pt and mt
  sc_full[["percent.pt"]] <- PercentageFeatureSet(sc_full, pattern="^ATCG")
  sc_full[["percent.mt"]] <- PercentageFeatureSet(sc_full, pattern="^ATMG")
  
  # Removing the useless cells.
  sc_full <- subset(sc_full, percent.pt < 5)
  sc_full <- subset(sc_full, percent.mt < 5)
  sc_full$rep_label <- c(sc_full[['rep_label']])
  
  sc_full <- Seurat::SCTransform(sc_full, vars.to.regress="percent.pt")
  sc_full <- Seurat::FindVariableFeatures(sc_full, selection.method="vst", nfeatures=2000)
  sc_full <- Seurat::RunPCA(sc_full, npcs = 50, features = Seurat::VariableFeatures(object = sc_full), seed.use = args$random_seed)
  sc_full <- Seurat::RunUMAP(sc_full, dims=1:30)
  # sc_full <- Seurat::FindVariableFeatures(sc_full, selection.method="vst", nfeatures=2000)
  sc_full <- Seurat::FindNeighbors(sc_full, reduction='pca', dims=1:30, verbose = TRUE)
  sc_full <- Seurat::FindClusters(sc_full, resolution=1, random.seed = args$random_seed, verbose = TRUE)
  saveRDS(sc_full, file=RDS_FILENAME)
} else {
  sc_full <- readRDS(RDS_FILENAME)
}

# ElbowPlot to determine the dimensionality of the dataset
clusters <- Seurat::Idents(sc_full)

# Visualization with UMAP, and ACP
ElbowPlot(sc_full)
Seurat::DimPlot(sc_full, reduction="umap", pt.size=0.4)
Seurat::DimPlot(sc_full, reduction='umap', group.by='rep_label', pt.size=0.2)

# Visualization of cell specific expression of selected genes in selected clusters
clusterChosen <- c(1, 2, 5, 8, 9, 10, 11, 13, 15, 16)
Seurat::DimPlot(subset(x=sc_full, idents=clusterChosen), reduction="umap", pt.size=0.2)

########################################################################################################################

# Next we build a dataframe that will contain the number of cells in each cluster for each rep
rep_label_frame  <- data.frame(sc_full[['rep_label']])
cluster_frame    <- data.frame(clusters)
df               <- data.frame(list(rep_label_frame, cluster_frame))

nb_clusters                  <- max(as.numeric(clusters))
# aggregate the dataframe so that we have the number of cells for each cluster for each rep
aggreg_df                    <- aggregate(df$rep_label, list(df$clusters, df$rep_label), FUN=length, drop=FALSE)
aggreg_df[is.na(aggreg_df)]  <- 0 # Replace the NA with zeros
abs_repart_cluster           <- data.frame( 0:(nb_clusters -1), list_reps, check.names=FALSE)
rownames(abs_repart_cluster) <- 0:(nb_clusters -1)
abs_repart_cluster[-1]       <- 0 # set all the entries to zero for clarity
abs_repart_cluster           <- abs_repart_cluster[-1]  # Remove the unwanted col
colnames(abs_repart_cluster) <- list_reps # set the right names of the columns

for (col in colnames(abs_repart_cluster)) {
  ind_rep_col <- (aggreg_df$Group.2 == col)
  abs_repart_cluster[col] <- aggreg_df$x[ind_rep_col]
}

rel_repart_cluster <- mapply("/", abs_repart_cluster, colSums(abs_repart_cluster))
plot_list_hist <- list() # vector('list', nrow(rel_repart_cluster))
group <- unlist(list_reps)

theme_bw(base_size = 20)

counts_by_cluster_and_rep_df <- matrix(ncol = 3)
clusters_with_replica_effect <- c()
clusters_without_replica_effect <- c()


for (i in 1:nb_clusters) {
  rel_repart <- rel_repart_cluster[i,]
  abs_repart <- abs_repart_cluster[i,]
  rel_value  <- as.numeric(rel_repart)
  abs_value  <- as.numeric(abs_repart)
  rel_df     <- data.frame(group, rel_value)
  abs_df     <- data.frame(group, abs_value)
  
  if (entropy(abs_df[,2]) < 0.5) {
    clusters_with_replica_effect <- c(clusters_with_replica_effect, i-1)
    next
  }
  
  clusters_without_replica_effect <- c(clusters_without_replica_effect, i-1) 
  # print(i)
  counts_by_cluster_and_rep_df <- rbind(counts_by_cluster_and_rep_df, abs_df[,2])
  
  bh <- ggplot(abs_df, aes(x=group, y=abs_value))
  bh <- bh + geom_bar(width=1, stat="identity" ) 
  bh <- bh + geom_text(aes(label=abs_value), vjust=1, size=3)
  bh <- bh + xlab(glue::glue("cluster {i-1}"))
  
  plot_list_hist <- c(plot_list_hist, list(bh))
  # plot_list_hist[[i]] <- bh
}

ggarrange(plotlist=plot_list_hist)

########################################################################################################################

# ¿QUÉ ES ESTO? NO ESTÁ GENERANDO NINGÚN PLOT
# Cluster repartition
pp <- ggplot(data.frame(clusters), aes(x=as.numeric(clusters)))
pp <- geom_bar(aes(x=clusters), position = "dodge", stat = "count")
print(pp)

dev.off()

########################################################################################################################
# Get the markers for each cluster
sc_full.markers <- Seurat::FindAllMarkers(sc_full, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, verbose = TRUE, random.seed = args$random_seed)
write.table(sc_full.markers, file=glue::glue("{args$output_folder}/gene_markers.txt"), sep="\t", quote=FALSE)

########################################################################################################################
###### MONOCLE
Seurat_obj <- readRDS(file=RDS_FILENAME)
nbClusters <- length(summary(Idents(Seurat_obj)))
print('The clusters are: ')
print(summary(Idents(Seurat_obj)))

########################################################################################################################
pdf(file=glue::glue('{args$output_folder}/first_plots_{basename(RDS_FILENAME)}{paste(as.character(clusterChosen), collapse="-")}.pdf'))

Seurat_obj[['clusters']]        <- Idents(Seurat_obj)
Seurat_obj                      <- subset(Seurat_obj, clusters %in% clusterChosen)
Seurat_obj@active.assay         <- 'RNA'

# From Seurat object to Monocle object
cds                             <- SeuratWrappers::as.cell_data_set(Seurat_obj)
cds                             <- monocle3::estimate_size_factors(cds)
rowData(cds)$gene_name          <- rownames(cds)
rowData(cds)$gene_short_name    <- rowData(cds)$gene_name
cds                             <- monocle3::preprocess_cds(cds, num_dim=100)
## Cluster the cells and compute the partitions. We will change the cluster with Seurat clusters.
cds                             <- monocle3::cluster_cells(cds)
cds@clusters[['UMAP']]$clusters <- Idents(Seurat_obj) ## Assign the clusters of the cds to the Seurat clusters
cds                             <- monocle3::learn_graph(cds) # Learn the graph, takes some time

## Plot the cells, according to their  clusters.
monocle3::plot_cells(cds, color_cells_by='cluster', label_groups_by_cluster=TRUE, labels_per_group=0)

get_earliest_principal_node <- function (cds, time_bin="130-170") {
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes  <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex))))]
  root_pr_nodes
}

# cds <- monocle3::order_cells(cds)
## Remove the line below if you want to choose the node yourself, and uncomment the line above.
cds <- monocle3::order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

# plot the cells by pseudotime
monocle3::plot_cells(cds, color_cells_by='pseudotime')

nbClusters <- length(summary(Idents(Seurat_obj)))

print(glue::glue('There are {nbClusters} clusters in total.'))
print(glue::glue('The clusters are: {clusterChosen}'))

nbClustersChosen <- length(summary(Idents(Seurat_obj)))

#### IF YOU WANT TO GET ALL THE CLUSTERS ALREADY CHOSEN BEFORE RUN THIS
subClusterChosen <- clusterChosen
# get only the wanted clusters the cell_data_set. It will be stored in sub_cds
BoolInChosenCluster <- as.vector(cds@clusters[['UMAP']][3][1]$clusters) %in% subClusterChosen

clusters                        <- as.vector(cds@clusters[['UMAP']][3][1]$clusters)
clusterIdOfCellsInChosenCluster <- clusters[clusters %in% subClusterChosen]
cellsInChosenCluster            <- row.names(subset(pData(cds),BoolInChosenCluster))
sub_cds                         <- cds[, cellsInChosenCluster]

## Takes a lot of time
sub_ciliated_cds_pr_test_res <- getCiliatedCds(sub_cds)

qValueChosen        <- 5e-10 # 0.0000000005
geneModuleAndAggMat <- getAggregatedMatrix(sub_cds, sub_ciliated_cds_pr_test_res, clusterIdOfCellsInChosenCluster, qValueChosen)
gene_module_df      <- geneModuleAndAggMat[[1]]
agg_mat             <- geneModuleAndAggMat[[2]]
showAggregatedMat(agg_mat)
dev.off()

########################################################################################################################
## If you want to select genes and plot their expression levels in the different cells
pdf(file=paste('output/first_plots_', rdsFile, paste(as.character(clusterChosen), collapse="-"), '.pdf', sep=""))

plot_cells(
  cds, 
  genes=c("AT2G23770","SC-AT1G78370","AT2G39518","MC2-AT3G46230","MC2-AT3G56880", "AT3G16570", "AT5G42980", "AT2G33380"),
  show_trajectory_graph=FALSE,
  label_cell_groups=FALSE,
  label_leaves=FALSE
)

## If you want to plot mean expression of genes (here Marker_cell) in  different clusters
Marker_cell <- read.csv("genes_to_plot.txt", header = FALSE)[,1]

plot_genes_by_group(
  cds,
  Marker_cell,
  group_cells_by="clusters",
  ordering_type="cluster_row_col",
  max.size=3
)

dev.off()

get_citations(cds)
