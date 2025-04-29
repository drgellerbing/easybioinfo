getresults <- function(deseq_dds){
  cs_dds <- class(deseq_dds) == "DESeqDataSet"
  if(cs_dds != "TRUE"){
    message("Please provide the DESeqDataSet from the rundeseq function")
    return(NULL)
  }

  ## Function for naming lists for DESeq results
  namelist <- function(x, deseq_dds){
    names(x) <- c(DESeq2::resultsNames(deseq_dds)[-1])
    names(x) <- sub("^[^_]*_", "", names(x))
    names(x) <- gsub("_", " ", names(x))
    names(x) <- gsub("\\.", " ", names(x))
    return(x)
  }

  message("This function will produce the significant results in excel files based on your provided cutoffs")
  nr <- as.numeric(length(DESeq2::resultsNames(deseq_dds)))
  p = as.numeric(readline(prompt = "Please enter the desired adjusted p-value between 0 and 1: "))
  message("The log2 fold change(LFC) is usually between 0.6 to 1")
  lfc = as.numeric(readline(prompt = "Please enter the desired LFC cutoff: "))

  if (nr > 2){
    x <- as.numeric(nr-1)
    l1 <- vector("list", length = x)
    for(i in 1:x){
      ddsn <- i + 1
      cat("Now processing", DESeq2::resultsNames(deseq_dds)[ddsn], "\n")
      res = DESeq2::results(deseq_dds, name = DESeq2::resultsNames(deseq_dds)[ddsn], alpha = p)
      res %>% as.data.frame() -> resdf
      resdf$genes <- rownames(resdf)
      sigresults <- filtersigres(resdf, p, lfc)
      l1[[i]] <- list(res, resdf, sigresults)
      names(l1[[i]]) <- c("deseqres", "fullres", "sigres")
    }
    l1 <- namelist(l1, deseq_dds)

    message("DESeq2 provides the function to compare the log fold change among the groups other than the control")
    com_res = readline(prompt = "Do you want the result comparisons among the query groups other than the control? (yes/no): ")
    if(tolower(com_res) == "yes" || tolower(com_res) == "y"){
      det_res = readline(prompt = "Among all your experimental control other than the control, do you want to compare all of them? (yes/no): ")

      if(det_res == "yes" || det_res == "y"){
        cdn <- as.character(unique(deseq_dds[[2]]))
        cat("The following are your provided experimental conditions: \n")
        for(cdni in 1:length(cdn)){
          cat(cdni, ": ", paste(cdn[cdni]), "\n", sep = "")
        }
        n_query_ctrl <- as.numeric(readline(prompt = "Please enter the number that corresponds to your control: "))
        query_ctrl <- as.character(cdn[n_query_ctrl])
        cdn <- stats::relevel(cdn, ref = query_ctrl)
        query <- levels(cdn)[-1]
        query_com <- utils::combn(query, 2)
        qcn <- query_com %>% as.data.frame() %>% ncol()

        l2 <- vector("list", length = qcn)

        for (b in 1:qcn){
          query_group <- query_com[,b]
          cat("Now processing", query_group[1], "vs", query_group[2], "\n")
          res <- DESeq2::results(deseq_dds, contrast = c(as.character(DESeq2::design(dds))[[2]], query_group), alpha = p)
          res %>% as.data.frame() -> resdf
          resdf$genes <- rownames(resdf)
          sigresults <- filtersigres(resdf, p, lfc)
          l2[[b]] <- list(res, resdf, sigresults)
          names(l2[[b]]) <- c("deseqres", "fullres", "sigres")

          y <- query_group[1]
          z <- query_group[2]
          names(l2)[b] <- paste(y, "vs", z, sep = " ")
        }
        #End of processing all results loop
      }

      if(det_res == "no" || det_res == "n"){
        cdn <- as.character(unique(deseq_dds[[2]]))
        cdn_len <- as.numeric(length(cdn))
        message("The following are your provided query groups: ")
        for (n in 1:cdn_len){cat(paste(n, cdn[n], sep = ": "), sep = "\n")}

        res_num <- as.numeric(readline(prompt = "How many results do you want to produce?: "))

        l2 <- vector("list", length = res_num)

        for (g in 1:res_num){
          cat("Now processing result", g, "\n")
          ng1 <- as.numeric(readline(prompt = "Please enter the number of your first query group: "))
          ng2 <- as.numeric(readline(prompt = "Please enter the number of your second query group: "))
          g1 <- cdn[ng1]
          g2 <- cdn[ng2]
          res = DESeq2::results(deseq_dds, contrast = c(as.character(DESeq2::design(dds))[[2]], g1, g2), alpha = p)
          res %>% as.data.frame() -> resdf
          resdf$genes <- rownames(resdf)
          sigresults <- filtersigres(resdf, p, lfc)
          l2[[g]] <- list(res, resdf, sigresults)
          names(l2[[g]]) <- c("deseq_dds", "fullres", "sigres")
          names(l2)[g] <- paste(g1, "vs", g2, sep = " ")
        }
      }
      #End of processing certain results loop
      l1 <- c(l1, l2)
    #End of comparing results loop
    }

    save_df <- readline(prompt = "Do you want to save your results into individual excel files? (yes/no): ")
    if(tolower(save_df) == "yes" || tolower(save_df) == "y"){
      message("Which results do you want to save? All results(both), unfiltered results only(unfil), or significant results only(sig)?")
      save_res <- readline(prompt = "Please enter the results you want to save (both/unfil/sig): ")
      dir <- readline(prompt = "Do you want to create a new folder for the saved results? (yes/no): ")
      if(tolower(dir) == "yes" || tolower(dir) == "y"){dir_n <- createdir()}
      nlist <- as.numeric(length(l1))

      for (nl in seq_along(l1)){
        filen <- names(l1)[nl]
        cat("Now saving", filen, "\n")
        res1 <- as.data.frame(l1[[nl]][[2]])
        res2 <- as.data.frame(l1[[nl]][[3]])

        if(tolower(dir) == "yes" || tolower(dir) == "y"){
          name1 <- paste("./", dir_n, "/", filen, "_fullres.xlsx", sep = "")
          name2 <- paste("./", dir_n, "/", filen, "_sigres.xlsx", sep = "")}

        if(tolower(dir) == "no" || tolower(dir) == "n"){
          name1 <- paste("./", filen, "_fullres.xlsx", sep = "")
          name2 <- paste("./", filen, "_sigres.xlsx", sep = "")}

        if(tolower(save_res) == "both"){
          writexl::write_xlsx(res1, name1)
          writexl::write_xlsx(res2, name2)
        }

        if(tolower(save_res) == "unfil"){writexl::write_xlsx(res1, name1)}
        if(tolower(save_res) == "sig"){writexl::write_xlsx(res2, name2)}
      }
      #End of saving result loop
    }
    #End of processing more than one results loop
    return(l1)
  }
  else{
    l1 <- vector("list", length = 1)
    res = DESeq2::results(deseq_dds, alpha = p)
    res %>% as.data.frame() -> resdf
    resdf$genes <- rownames(resdf)
    sigresults <- filtersigres(resdf, p, lfc)

    filen <- readline(prompt = "Please provide the name of your result: ")
    name1 <- paste("./", filen, "_fullres.xlsx", sep = "")
    name2 <- paste("./", filen, "_sigres.xlsx", sep = "")
    name3 <- paste("./", filen, "_deletedgenes.xlsx", sep = "")
    writexl::write_xlsx(resdf, name1)
    writexl::write_xlsx(sigresults, name2)
    l1[[1]] <- list(res, resdf, sigresults)
    names(l1[[1]]) <- c("deseqres", "fullres", "sigres")
    l1 <- namelist(l1, deseq_dds)

    return(l1)
  }
  cat("Done :>")
}
