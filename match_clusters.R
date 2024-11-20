library(biomaRt)
library(dplyr)
library(readxl)

ensembl <- useMart("plants_mart", dataset = "athaliana_eg_gene", host = "https://plants.ensembl.org")

get_tair_ids <- function(gene_symbols) {
    tair_ids <- getBM(attributes = c("external_gene_name", "tair_locus"),
                  filters = "external_gene_name",
                  values = gene_symbols,
                  mart = ensembl)
    tair_ids
}

cluster_gene_df <- read.csv("gene_markers_2.txt", sep = "\t")
cluster_sheet_names <- c(
  "C0C4C10C11 healthy mesophyl", 
  "C1 responsive epidermal cells",
  "C2 vascular S cells",
  "C5",
  "C8 responsive mesophyl cells",
  "C9 healthy epidermal cells",
  "C13",
  "C16 guard cells"
)

clusters <- sapply(strsplit(cluster_sheet_names, " "), function(x) x[1])

cluster_names <- character()
clusters_paper <- list()

all_subcluster <- character()

for (cluster_idx  in seq_along(cluster_sheet_names)) {
    
    xls_file = "~/Postdoc/Papers/AT_PST_ScienceDirect_files_06Nov2024_13-40-25.581/1-s2.0-S2590346223002043-mmc5.xlsx"
    cluster_genes_paper_df <- readxl::read_xlsx(xls_file, sheet = cluster_sheet_names[cluster_idx])
    
    cluster <- clusters[cluster_idx]
    
    subclusters = na.omit(unique(cluster_genes_paper_df$cluster))
    # subclusters <- strsplit(cluster_sheet_names[[1]], "C")[[1]] %>% .[2:length(.)] %>% paste0("C", .) 
    # strsplit(rownames(overlap_matrix)[1], "C")[[1]] %>% .[2:length(.)] %>% paste0("C", .)
    all_subcluster <- c(all_subcluster, subclusters)
    
    print(unique(cluster_genes_paper_df$cluster))
    print(subclusters)
    
    for (subcluster in subclusters) {
        
        valid_rows <- !is.na(cluster_genes_paper_df['Gene model'])[,1]
        cluster_genes_paper_df <- cluster_genes_paper_df[valid_rows,]
        cluster_genes_paper_df <- cluster_genes_paper_df %>% filter(cluster == subcluster)
        genes <- cluster_genes_paper_df$`Gene model`
        
        clusters_paper <- c(clusters_paper, list(genes))
    
        # print(subcluster)
        # print(genes)#[1:sum(valid_rows)])
    }
}

print(glue::glue("All subclusters are {length(all_subcluster)}"))

clusters_mine = sapply(1:23, function(i) cluster_gene_df %>% filter(cluster == i) %>% rownames)

overlap_matrix <- matrix(0, nrow = length(clusters_paper), ncol = length(clusters_mine),
                         dimnames = list(names(clusters_paper), names(clusters_mine)))

for (i in seq_along(all_subcluster)) {
  print(i)
  for (j in seq_along(clusters_mine)) {
    
    genes_paper <- clusters_paper[[i]]
    genes_mine <- clusters_mine[[j]]
    
    tair_ids <- get_tair_ids(genes_mine)$tair_locus
    overlap <- length(intersect(genes_paper, tair_ids)) #  / length(tair_ids)
    overlap_matrix[i, j] <- overlap
  }
}

rownames(overlap_matrix) <- all_subcluster
colnames(overlap_matrix) <- paste0("C", 1:23)


suma_filas <- rowSums(overlap_matrix)
suma_columnas <- colSums(overlap_matrix)

# Reordenar filas y columnas en orden descendente
ordered_overlap_matrix <- overlap_matrix[order(-suma_filas), order(-suma_columnas)]