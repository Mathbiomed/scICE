# install.packages("Seurat")
# install.packages("sigclust")
# library(devtools)
# options(download.file.method = "wget")
# install_github("siyao-liu/MultiK")

library(stringr)
library(Seurat)
library(MultiK)

findOptK = function(tog,ups) {
  
  tog.f <- tog[tog$Freq > ups | tog$Freq ==ups, ]
  hpts <- chull(tog.f[, c("one_minus_rpac", "Freq")]) # in clockwise order
  hpts <- c(hpts, hpts[1])
  ch.df <- tog.f[hpts, ]
  
  df <- ch.df[ , c("ks", "one_minus_rpac", "Freq")]
  colnames(df) <- c("k", "x", "y")
  b <- c()
  end_points <- c()
  for (i in 1: (nrow(df)-1)) {
    end_points[i] <- paste(as.character(df[i, ]$k), as.character(df[(i+1),]$k), sep="-")
    b[i] <- (df[(i+1), ]$y - df[i, ]$y)/(df[(i+1), ]$x - df[i, ]$x)
  }
  
  # put in data frame for each line segment
  lineseg.df <- data.frame("end_points"=end_points, "slope"=b)
  lineseg.df$p1 <- do.call("rbind", strsplit(lineseg.df$end_points, "-"))[, 1]
  lineseg.df$p2 <- do.call("rbind", strsplit(lineseg.df$end_points, "-"))[, 2]
  
  # step 1: find k with largest # of runs
  which.k <- as.character(ch.df[which.max(ch.df$Freq), ]$ks)
  
  # step 2: see if a line segment with negative slope coming out
  if ( all(lineseg.df[lineseg.df$p1==which.k | lineseg.df$p2==which.k, ]$slope > 0) ) {
    optK <- which.k
  }
  else {
    
    # follow the line segment with the negative slope
    tmp <- which(lineseg.df[lineseg.df$p1==which.k | lineseg.df$p2==which.k, ]$slope < 0)
    tmp <- lineseg.df[lineseg.df$p1==which.k | lineseg.df$p2==which.k, ][tmp, ]
    which.k2 <- as.character(c(tmp$p1, tmp$p2)[which(c(tmp$p1, tmp$p2)!=which.k)])
    
    # check if the slope becomes more negative
    lineseg.df.sub <- lineseg.df[lineseg.df$p1!=which.k & lineseg.df$p2 !=which.k, ]
    
    if ( #any(lineseg.df[lineseg.df$p1==which.k2 | lineseg.df$p2==which.k2, ]$slope > 0)
      lineseg.df.sub[lineseg.df.sub$p1==which.k2 | lineseg.df.sub$p2 == which.k2, ]$slope > tmp$slope ) {
      optK <- c(which.k, which.k2)
    }
    
    else {
      tmp <- which(lineseg.df.sub[lineseg.df.sub$p1==which.k2 | lineseg.df.sub$p2==which.k2, ]$slope < 0)
      tmp <- lineseg.df.sub[lineseg.df.sub$p1==which.k2 | lineseg.df.sub$p2==which.k2, ][tmp, ]
      which.k3 <- as.character(c(tmp$p1, tmp$p2)[ which(c(tmp$p1, tmp$p2)!=which.k & c(tmp$p1, tmp$p2)!=which.k2)])
      
      lineseg.df.sub <- lineseg.df[lineseg.df$p1!=which.k & lineseg.df$p2 !=which.k
                                   & lineseg.df$p1!=which.k2 & lineseg.df$p2 !=which.k2, ]
      
      if ( lineseg.df.sub[lineseg.df.sub$p1==which.k3 | lineseg.df.sub$p2 == which.k3, ]$slope > tmp$slope ) {
        optK <- c(which.k, which.k2, which.k3)
      }
      else {
        optK <- c(which.k, which.k2, which.k3)
      }
    }
  }
  return(optK)
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


if (.Platform$OS.type == "windows"){
  tmp_dir = 'C:\\Users\\kimma\\SynologyDrive\\Data\\gene_test_data\\harmony_ex\\balanced_new\\preprocessed'
} else {
  tmp_dir = '/home/khlab/Documents/Data/gene_test_data/harmony_ex/balanced_new/preprocessed'
}

multiK_perform <- function(cbmc) {
  multik <- MultiK(cbmc)
  ks = multik$k
  res = multik$consensus
  tog <- as.data.frame(table(ks)[table(ks) > 1])
  # calculate rpac
  pacobj <- CalcPAC(x1=0.1, x2=0.9, xvec = tog$ks, ml = res)
  tog$rpac <- pacobj$rPAC
  tog$one_minus_rpac  <- 1-tog$rpac
  
  k = tryCatch(findOptK(tog,10),
               error = function(e) findOptK(tog,5),
               warning = function(w) findOptK(tog,5),
               finally = 3)
  
  cluster_list <- getClusters(cbmc,k)
  # cluster_list <- getClusters(cbmc,c(k[1],k[2]))
  
  cluster_df <- cluster_list$clusters
  cluster_df <- cbind(true_l = colnames(cbmc@assays$RNA$data), cluster_df)
  return(cluster_df)
}

files = sort(Sys.glob(file.path(tmp_dir,"*.csv.gz")))
f = files[1]
t_df = data.frame()
for (f in c(files[5],files[1])) {
  base_n = strsplit(basename(f),".csv.gz")[[1]]
  raw_d = as.sparse(read.tcsv(file = f, sep = ",",header = TRUE, row.names = 1))

  cbmc.rna <- raw_d
  cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)
  cbmc <- CreateSeuratObject(counts = cbmc.rna)
  cbmc <- NormalizeData(cbmc)
  
  s_time = Sys.time()
  results <- multiK_perform(cbmc)
  e_time = Sys.time()

  
  t_diff = as.numeric(e_time - s_time,unit="secs")
  t_df = rbind(t_df,data.frame(filename = base_n,end_time = as.numeric(Sys.time()),d_time = t_diff))
  
}

outname3 = str_glue(dirname(dirname(f)),"/","multiK_mem2.csv")
write.csv(file=outname3,t_df)



