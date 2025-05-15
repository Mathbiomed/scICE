library(Seurat)
library(ggplot2)
library(stringr)


if (.Platform$OS.type == "windows"){
  source('c:/Users/kimma/SynologyDrive/Code/R/chooseR/pipeline.R')
  parallel_flag = FALSE
  # tmp_dir = 'C:\\Users\\kimma\\SynologyDrive\\Data\\Project_with_Faeyza\\data'
  tmp_dir = 'C:\\Users\\kimma\\SynologyDrive\\Data\\gene_test_data\\harmony_ex\\balanced_new\\preprocessed'
  results_path <- "c:/Users/kimma/Downloads/"
} else {
  source('/home/khlab/Documents/Code/R/chooseR/pipeline.R')
  parallel_flag = TRUE
  # tmp_dir = '/home/khlab/Documents/Data/Project_with_Faeyza/data/'
  tmp_dir = '/home/khlab/Documents/Data/gene_test_data/harmony_ex/balanced_new/preprocessed'
  results_path <- "/home/khlab/다운로드/chooseR_results/"
}

read.tcsv = function(file, header=TRUE, sep=",", ...) {
  
  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)
  
  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }
  
  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep) 
  out = read.csv(text=x, sep=sep, header=header, ...)
  return(out)
  
}

run_pipeline <- function(obj, resolutions) {
  for (res in resolutions) {
    message(paste0("Clustering ", res, "..."))
    message("\tFinding ground truth...")
    
    # "Truths" will be stored at glue::glue("{reduction}.{assay}_res.{res}")
    obj <- find_clusters(
      obj,
      reduction = reduction,
      assay = assay,
      resolution = res
    )
    clusters <- obj[[glue::glue("{reduction}.{assay}_res.{res}")]]
    
    # Now perform iterative, sub-sampled clusters
    results <- multiple_cluster(
      obj,
      n = 100,
      size = 0.8,
      npcs = npcs,
      res = res,
      reduction = reduction,
      assay = assay
    )
    
    # Now calculate the co-clustering frequencies
    message(paste0("Tallying ", res, "..."))
    # This is the more time efficient vectorisation
    # However, it exhausts vector memory for (nearly) all datasets
    # matches <- purrr::map(columns, find_matches, df = results)
    # matches <- purrr::reduce(matches, `+`)
    columns <- colnames(dplyr::select(results, -cell))
    mtchs <- matrix(0, nrow = dim(results)[1], ncol = dim(results)[1])
    i <- 1 # Counter
    for (col in columns) {
      message(paste0("\tRound ", i, "..."))
      mtchs <- Reduce("+", list(
        mtchs,
        find_matches(col, df = results)
      ))
      i <- i + 1
    }
    
    message(paste0("Scoring ", res, "..."))
    mtchs <- dplyr::mutate_all(
      dplyr::as_tibble(mtchs),
      function(x) dplyr::if_else(Re(x) > 0, percent_match(x), 0)
    )
    
    # Now calculate silhouette scores
    message(paste0("Silhouette ", res, "..."))
    sil <- cluster::silhouette(
      x = as.numeric(as.character(unlist(clusters))),
      dmatrix = (1 - as.matrix(mtchs))
    )
    saveRDS(sil, paste0(results_path, "silhouette_", res, ".rds"))
    
    # Finally, calculate grouped metrics
    message(paste0("Grouping ", res, "..."))
    # grp <- group_scores(mtchs, unlist(clusters))
    # saveRDS(grp, paste0(results_path, "frequency_grouped_", res, ".rds"))
    sil <- group_sil(sil, res)
    saveRDS(sil, paste0(results_path, "silhouette_grouped_", res, ".rds"))
  }
  
  # Save original data, with ground truth labels
  saveRDS(obj, paste0(results_path, "clustered_data.rds"))
  
  # Create silhouette plot
  # Read in scores and calculate CIs
  scores <- purrr::map(
    paste0(results_path, "silhouette_grouped_", resolutions, ".rds"),
    readRDS
  )
  scores <- dplyr::bind_rows(scores) %>%
    dplyr::group_by(res) %>%
    dplyr::mutate("n_clusters" = dplyr::n()) %>%
    dplyr::ungroup()
  meds <- scores %>%
    dplyr::group_by(res) %>%
    dplyr::summarise(
      "boot" = list(boot_median(avg_sil)),
      "n_clusters" = mean(n_clusters)
    ) %>%
    tidyr::unnest_wider(boot)
  writexl::write_xlsx(meds, paste0(results_path, "median_ci.xlsx"))
  
  # Find thresholds
  threshold <- max(meds$low_med)
  choice <- as.character(
    meds %>%
      dplyr::filter(med >= threshold) %>%
      dplyr::arrange(n_clusters) %>%
      tail(n = 1) %>%
      dplyr::pull(res)
  )
  label_p <- obj@meta.data[str_glue("pca.SCT_res.",choice[[1]])]  
  return(label_p)
}


tmp_dir = commandArgs(trailingOnly=TRUE)

files = sort(Sys.glob(file.path(tmp_dir,"*.csv.gz")))
npcs <- 100
resolutions <- c(0.3, 0.5, 0.8, 1, 1.2, 1.6, 2, 4, 6, 8)
assay <- "SCT"
reduction <- "pca"
f = files[4]
t_df = data.frame()
for (f in files) {
  base_n = strsplit(basename(f),".csv.gz")[[1]]
  raw_d = as.sparse(read.tcsv(file = f, sep = ",",header = TRUE, row.names = 1))
  
  obj.rna <- raw_d
  obj.rna <- CollapseSpeciesExpressionMatrix(obj.rna)
  obj <- CreateSeuratObject(counts = obj.rna)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs=npcs,verbose = FALSE)
  
  s_time = Sys.time()
  result <- run_pipeline(obj,resolutions)
  e_time = Sys.time()
  
  t_diff = as.numeric(e_time - s_time,unit="secs")
  t_df = rbind(t_df,data.frame(filename = base_n,end_time = as.numeric(Sys.time()),d_time = t_diff))
}

outname3 = str_glue(dirname(dirname(f)),"/","chooseR_mem.csv")
write.csv(file=outname3,t_df)