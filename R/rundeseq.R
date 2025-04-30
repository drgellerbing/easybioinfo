#' @import magrittr
#' @title Differential Expression Analysis Based on the Gamma-Poisson Distribution
#' @description
#' A streamlined pipeline that performs differential expression analysis through a series of prompts using the "DESeq2" package.
#' Produces DESeqDataSet with DESeq2::DESeqDataSetFromMatrix, establishes the control through your prompts, and runs the differential expression analysis using DESeq2::DESeq
#' @param expression Gene expression dataframe
#' @param md Clinical data dataframe
#' @returns a DESeqDataSet object with results stored as metadata columns. These results can be accessed by calling results function from "DESeq2" or from the getresults function
#' @examples
#' df <- easybioinfo::deseqexpr
#' md <- easybioinfo::deseqmd
#' dds <- rundeseq(df, md)

rundeseq <- function(expression, md){
  expression %<>% as.data.frame()
  View(expression)
  suppressWarnings({
    gene_col <- as.numeric(readline(prompt = "Please enter the column number of your gene id in your expression dataframe: "))
    while(is.na(gene_col) == TRUE || gene_col > ncol(expression)){
      message("Error detected, please try again :)")
      gene_col <- as.numeric(readline(prompt = "Please enter the column number of your gene id in your expression dataframe: "))
    }
  })
  colnames(expression)[gene_col] <- "geneid"

  #Checking for duplicates
  tryCatch({
    rownames(expression) <- expression$geneid
  }, error = function(e){
    cat("Error: Duplicates found in your gene id. Please make sure there are no duplicates \n")
    stop(e)
  })

  expression %<>% dplyr::relocate(geneid)
  expression[, -1] <- lapply(expression[, -1], as.integer)
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

  condlabel <- unique(md$cond)
  cat("Here are your provided experimental conditions: \n")
  for(cl in 1:length(condlabel)){
    cat(cl, ": ", paste(condlabel[cl]), "\n", sep = "")
  }
  n_ctrl <- as.numeric(readline(prompt = "Enter the number that corresponds to your control: "))
  ctrl <- condlabel[n_ctrl]

  md$cond <- as.factor(md$cond)
  expression %<>% dplyr::select(geneid, dplyr::any_of(md$sampid))
  dds = DESeq2::DESeqDataSetFromMatrix(countData = expression,
                               colData = md,
                               design = ~cond, tidy = TRUE)
  dds$cond <- stats::relevel(dds$cond, ref = ctrl)
  n_groups <- length(unique(dds$cond))

  prefil <- readline(prompt = "Do you want to pre-filter the count input? (yes/no): ")
  if (tolower(prefil) == "yes" || tolower(prefil) == "y"){
    count <- as.numeric(readline(prompt = "Enter your desired filter value. Default = 10: "))
    print(forcats::fct_count(as.factor(md$cond)))
    spl <- as.numeric(readline(prompt = "Please enter the number of the smallest biological replicates among your experimental conditions: "))
    keep <- rowSums(DESeq2::counts(dds) >= count) >= spl
    dds <- dds[keep,]
  }

  rundds = DESeq2::DESeq(dds, fitType = c("parametric", "local", "mean", "glmGamPoi"))

  if (n_groups == 2){
    print(utils::head(DESeq2::results(rundds), tidy = TRUE))
    print(DESeq2::summary(DESeq2::results(rundds, alpha = 0.05)))
    return(rundds)
  }
  else{
    nr <- as.numeric(DESeq2::resultsNames(rundds) %>% as.data.frame() %>% nrow())
    for (i in 2:nr){
      print(utils::head(DESeq2::results(rundds, name = DESeq2::resultsNames(rundds)[i]), tidy = TRUE))
      print(DESeq2::summary(DESeq2::results(rundds, name = DESeq2::resultsNames(rundds)[i], alpha = 0.05)))
    }
    return(rundds)
  }
  cat("Analysis complete :>")
}
