#' @import magrittr

rundeseq <- function(expr, md){
  expr %<>% as.data.frame()
  suppressWarnings({
    gene_col <- as.numeric(readline(prompt = "Please enter the column number of your gene id in your expression dataframe: "))
    while(is.na(gene_col) == TRUE || gene_col > ncol(expr)){
      message("Error detected, please try again :)")
      gene_col <- as.numeric(readline(prompt = "Please enter the column number of your gene id in your expression dataframe: "))
    }
  })
  colnames(expr)[gene_col] <- "geneid"
  expr %<>% dplyr::relocate(geneid)
  expr[, -1] <- lapply(expr[, -1], as.integer)
  View(md)
  suppressWarnings({
    sam_col <- as.numeric(readline(prompt = "Please enter the column number of your sample IDs in your metadata: "))
    while(is.na(sam_col) == TRUE || sam_col > ncol(md)){
      sam_col <- as.numeric(readline(prompt = "Please enter the column number of your sample IDs in your metadata: "))
    }
  })
  suppressWarnings({
    cdn_col <- as.numeric(readline(prompt = "Please enter the column number of your experimental condition in your metadata: "))
    while(is.na(cdn_col) == TRUE || cdn_col > ncol(md)){
      cdn_col <- as.numeric(readline(prompt = "Please enter the column number of your experimental condition in your metadata: "))
    }
  })
  colnames(md)[sam_col] <- "sampid"
  colnames(md)[cdn_col] <- "cond"
  md$cond <- as.factor(md$cond)
  expr %<>% dplyr::select(geneid, any_of(md$sampid))
  dds = DESeq2::DESeqDataSetFromMatrix(countData = expr,
                               colData = md,
                               design = ~cond, tidy = TRUE)
  ctrl <- readline(prompt = "Enter the control: ")
  dds$cond <- stats::relevel(dds$cond, ref = ctrl)
  n_groups <- as.numeric(readline(prompt = "Enter the number of experimental conditions excluding control: "))
  p = as.numeric(readline(prompt = "Enter your desired adjusted p-value: "))

  prefil <- readline(prompt = "Do you want to pre-filter the count input? (yes/no): ")
  if (tolower(prefil) == "yes" || tolower(prefil) == "y"){
    count <- as.numeric(readline(prompt = "Enter your desired filter value. Default = 10: "))
    print(forcats::fct_count(as.factor(md$cond)))
    spl <- as.numeric(readline(prompt = "Please enter the number of the smallest biological replicates among your experimental conditions: "))
    keep <- rowSums(counts(dds) >= count) >= spl
    dds <- dds[keep,]

    rundds = DESeq2::DESeq(dds, fitType = c("parametric", "local", "mean", "glmGamPoi"))

    if (n_groups == 1){
      print(head(DESeq2::results(rundds), tidy = TRUE))
      print(summary(DESeq2::results(rundds, alpha = p)))
      return(rundds)
    }
    else{
      nr <- as.numeric(DESeq2::resultsNames(rundds) %>% as.data.frame() %>% nrow())
      for (i in 2:nr){
        print(head(DESeq2::results(rundds, name = DESeq2::resultsNames(rundds)[i]), tidy = TRUE))
        print(summary(DESeq2::results(rundds, name = DESeq2::resultsNames(rundds)[i], alpha = p)))
      }
      return(rundds)
    }
  }
  else{
    rundds = DESeq2::DESeq(dds, fitType = c("parametric", "local", "mean", "glmGamPoi"))

    if (n_groups == 1){
      print(head(DESeq2::results(rundds), tidy = TRUE))
      print(summary(DESeq2::results(rundds, alpha = p)))
      return(rundds)
    }
    else{
      nr <- as.numeric(DESeq2::resultsNames(rundds) %>% as.data.frame() %>% nrow())
      for (i in 2:nr){
        print(head(DESeq2::results(rundds, name = DESeq2::resultsNames(rundds)[i]), tidy = TRUE))
        print(summary(DESeq2::results(rundds, name = DESeq2::resultsNames(rundds)[i], alpha = p)))
      }
      return(rundds)
    }
  }
  cat("Analysis complete :>")
}
