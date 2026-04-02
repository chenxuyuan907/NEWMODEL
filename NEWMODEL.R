### Bayesprism
library(BayesPrism)
library(Seurat)
library(dplyr)
load("bulk_data.RData")
load("sc_data.RData")

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = sc_data, 
                                 project = "scRNA_project", 
                                 min.cells = 3, 
                                 min.features = 200)

# Data preprocessing
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & 
                       percent.mt < 20)
seurat_obj <- NormalizeData(seurat_obj, 
                            normalization.method = "LogNormalize", 
                            scale.factor = 10000)

seurat_obj <- FindVariableFeatures(seurat_obj, 
                                   selection.method = "vst", 
                                   nfeatures = 2000)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15) 
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)
cell.type.labels <- Idents(seurat_obj)

filtered_cells <- colnames(seurat_obj)
all_cells <- colnames(sc_data)
filtered_out_cells <- setdiff(all_cells, filtered_cells)
sc_data <- sc_data[, !(colnames(sc_data) %in% filtered_out_cells)]
scRNA_data=t(sc_data)

#BayesPrism main function
myPrism <- new.prism(
  reference=scRNA_data, 
  mixture=bulk_data,
  input.type="GEP",
  cell.type.labels = cell.type.labels, 
  cell.state.labels = NULL,
  key=NULL,# 
  outlier.cut=0.01,
  outlier.fraction=0.1,
)
bp.res <- run.prism(prism = myPrism, n.cores=64)

### Run mr.mash
library(parallel)
library(mr.mash.alpha)
library(pbapply)  

# Get all gene names
gene_files <- list.files(zcellgene_dir, pattern = "\\.csv$")
gene_names <- sub("\\.csv$", "", gene_files)
# Take n genes for testing
test_genes <- gene_names[1:n]

# mr.mash main function
process_gene <- function(gene) {
  tryCatch({
    Y <- read.csv(file.path(zcellgene_dir, paste0(gene, ".csv")))
    X <- read.csv(file.path(genotype_dir, paste0(gene, "_sum.csv")))
    X <- X[, sapply(X, is.numeric)]
    X_std <- scale(X, center = TRUE, scale = TRUE)
    X <- as.matrix(X_std)
    Y <- as.matrix(Y)
    Y <- Y[, -1] 
    Y_num <- matrix(as.numeric(Y), nrow = nrow(Y), ncol = ncol(Y))
    colnames(Y_num) <- colnames(Y)
    rownames(Y_num) <- rownames(Y)
    Y <- scale(Y_num)
    Y[is.na(Y)] <- 0
    S0 <- compute_canonical_covs(ncol(Y))
    fit <- mr.mash(X, Y, S0, V = diag(ncol(Y)) * 0.00001 + cov(Y) * 0.99999, nthreads = 1)
    save(fit, file = file.path(output_dir, paste0(gene, "model_fit.RData")))
    return(paste(gene, "completed successfully"))
  }, error = function(e) {
    return(paste("Error processing", gene, ":", e$message))
  })
}

num_cores <- 64  
cl <- makeCluster(num_cores)
clusterExport(cl, c("zcellgene_dir", "genotype_dir", "output_dir", "compute_canonical_covs", "mr.mash"))
pboptions(type = "txt", style = 3, char = "=")
results <- pblapply(cl = cl, X = test_genes, FUN = process_gene)
stopCluster(cl)


### Calculate Z_score and use Cauchy combination to find P-value.
library(readr)
library(bigstatsr)
library(foreach)
library(doParallel)
cl <- makeCluster(64)
registerDoParallel(cl)
gene_chrom <- readRDS("genechrom_vector.rds")
gwas_data <- read_tsv(
  "expressionQCPC/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.updated.tsv", 
  col_select = c(CHROM, POS, BETA, SE)
)

analyze_gene <- function(i) {
  target_chrom <- as.character(gene_chrom[i])
  gene_name <- gsub(".RData", "", list.files("mrresult")[i])
  model_file <- paste0("mrresult/", gene_name, ".RData")
  load(model_file)
  beta <- fit$mu1
  SNPPOS <- gsub("^X", "", rownames(beta))
  target_data <- gwas_data[gwas_data$CHROM == target_chrom, ]
  target_pos <- as.numeric(target_data$POS) 
  snppos_num <- as.numeric(SNPPOS) 
  match_indices <- match(snppos_num, target_pos)
  if (anyNA(match_indices)) {
    match_indices[is.na(match_indices)] <- 1
  }
  Z <- matrix(
    target_data$BETA[match_indices] / target_data$SE[match_indices], 
    ncol = 1 
  )
  geno_file <- paste0("genotype/", gene_name, "_sum.csv")
  X <- read.csv(geno_file) 
  X_std <- scale(X[, sapply(X, is.numeric)])
  X_big <- as_FBM(as.matrix(X_std)) 
  V <- big_cor(X_big)[, ]
  betaVbeta <- diag(t(beta) %*% V %*% beta)
  betaVbeta <- abs(betaVbeta)  
  Zscore <- (t(beta) %*% Z) / sqrt(betaVbeta)
  robust_p <- function(Z, zmax = 20) {
    Z_trunc <- pmin(pmax(Z, -zmax), zmax)
    p <- 2 * pnorm(-abs(Z_trunc))
    0.5 - atan(mean(tan((0.5 - p) * pi))) / pi 
  }
  gene_p <- robust_p(Zscore)
  
  return(list(
    gene = gene_name,
    chrom = target_chrom,
    p_value = gene_p,
    z_score = as.numeric(Zscore)
  ))
}

results <- foreach(i = 1:n, .combine = 'rbind') %dopar% {
  tryCatch({
    analyze_gene(i)
  }, error = function(e) {
    list(gene = paste0("Gene_", i), chrom = NA, p_value = NA, z_score = NA)
  })
}

stopCluster(cl)
