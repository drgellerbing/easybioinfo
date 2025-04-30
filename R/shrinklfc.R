#' @title Shrink Log2 Fold Changes
#' @description
#' A simplified pipeline for the shrinkage of effect size (LFC estimates) which is useful for visualisation and ranking of genes.
#' Shrinks the DESeq2 results and filters the significant results based on your desired LFC and p-value that can be saved into excel files
#' Shrinkage of genes can be done using three different estimators: ashr/normal/apeglm(default).
#' Please visit http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.htmlhttp://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html for more information about the shrinkage estimators.
#' @param dds a DESeqDataSet object, either after running DESeq from DESeq2 package or rundeseq function
#' @returns A list that contains the DESeqResults object with shrunken LFC, the full shrunken result dataframe, and the significant shrunken result dataframe.
#' Both dataframes can be saved into excel files
#' @examples
#' expression = easybioinfo::deseqexpr
#' md = easybioinfo::deseqmd
#' 
#' exampledeseq = DESeq2::DESeqDataSetFromMatrix(countData = expression, colData = md, design = ~condition, tidy = TRUE)
#' dds <- DESeq2::DESeq(exampledeseq)
#' 
#' library(easybioinfo)
#' dds <- easybioinfo::rundeseq(expression, md)
#' lfc <- shrinklfc(dds)
#' 

shrinklfc <- function(dds){
  cs_dds <- class(dds) == "DESeqDataSet"
  if(cs_dds != "TRUE"){
    message("Please provide the DESeqDataSet from the rundeseq function")
    return(NULL)
  }

  tidyname <- function(x){
    x <- sub("^[^_]*_", "", x)
    x <- gsub("_", " ", x)
    x <- gsub("\\.", " ", x)
    return(x)
  }

  suppressWarnings({
    p <- as.numeric(readline(prompt = "Please state your desired p-value cutoff: "))
    while(is.na(p) == TRUE || p <= 0 || p > 1){
      message("The p-value provided should be within 0 and 1, please try again :>")
      p <- as.numeric(readline(prompt = "Please state your desired p-value cutoff: "))
    }
  })

  lfc <- as.numeric(readline(prompt = "Please enter your desired LFC cutoff: "))
  x = as.numeric(length(DESeq2::resultsNames(dds)))

  ## More than one result to shrink
  if (x > 2){
    ans1 <- readline(prompt = "Do you want to shrink all the results? (yes/no): ")
    message("Please note that the default estimator (apeglm) will only shrink results of your experimental conditions against the control")
    Sys.sleep(0.5)
    message("If you want to shrink the results comparing the exprimental conditions other than the control, \nplease change the shrinkage esimator to ashr/normal")
    message("However, do note that according to the author, apeglm and ashr tend to show less bias than normal")
    Sys.sleep(0.5)
    message("For more information about LFC shrinkage, please refer to the DESeq2 vignette")
    message("http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html")
    Sys.sleep(0.5)
    ans2 <- readline(prompt = "Do you want to change the shrinkage estimator? Default = apeglm (yes/no): ")

    ## Defining shrinking parameters
    if (tolower(ans2) == "no" || tolower(ans2) == "n"){shrinkingparameter = "apeglm"}

    if(tolower(ans2) == "yes" || tolower(ans2) == "y"){
      suppressWarnings({
        shrinkingparameter <- readline(prompt = "Please enter your desired shrinkage estimator(ashr/normal): ")
        while(is.na(shrinkingparameter) == TRUE || shrinkingparameter != "ashr" && shrinkingparameter != "normal"){
          shrinkingparameter <- readline(prompt = "Please enter your desired shrinkage estimator(ashr/normal): ")
        }
      })
    }

    ## Shrinking all results
    if (tolower(ans1) == "yes" || tolower(ans1) == "y"){
      y <- as.numeric(x - 1)

      l1 <- vector("list", length = y)
      for(i in 1:y){
        resn <- as.numeric(i + 1)
        resname <- c(DESeq2::resultsNames(dds)[-1])
        resname <- tidyname(resname)
        cat("Now shrinking result ", i, ": ", resname[i], "\n", sep = "")
        slr <- DESeq2::lfcShrink(dds, coef = resn, type = shrinkingparameter, res = DESeq2::results(dds, name = DESeq2::resultsNames(dds)[resn], alpha = p))
        slr %>% as.data.frame() -> slrdf
        slrdf$genes <- rownames(slrdf)
        rownames(slrdf) <- NULL
        sigresults <- filtersigres(slrdf, p, lfc)
        l1[[i]] <- list(slr, slrdf, sigresults)
        names(l1[[i]]) <- c("lfcdds", "shrinked_fullres", "shrinked_sigres")
        names(l1)[i] <- resname[i]
      }

      if(tolower(ans2) == "yes" || tolower(ans2) == "y"){
        com_res <- readline(prompt = "Do you want to shrink your results among the rest of your experimental conditions (yes/no): ")

        if(tolower(com_res) == "yes" || tolower(com_res) == "y"){
          cdn <- unique(dds[[2]])
          cat("The following are your provided experimental conditions: \n")
          for(cdni in 1:length(cdn)){
            cat(cdni, ": ", paste(cdn[cdni]), "\n", sep = "")
          }
          n_query_ctrl <- as.numeric(readline(prompt = "Please enter the number that corresponds to your control: "))
          query_ctrl <- as.character(cdn[n_query_ctrl])
          cdn <- stats::relevel(cdn, ref = query_ctrl)
          query <- levels(cdn)[-1]
          query_com <- utils::combn(query, 2)
          qcn <- as.numeric(query_com %>% as.data.frame() %>% ncol())

          l2 <- vector("list", length = qcn)

          for(g in 1:qcn){
            resn2 <- as.numeric(y + g)
            query_group <- query_com[,g]
            n1 <- query_group[1]
            n2 <- query_group[2]
            gslrname <- paste(n1, " vs ", n2, sep = "")
            cat("Now shrinking result ", resn2, ": ", gslrname, "\n", sep = "")
            slr2 <- DESeq2::lfcShrink(dds, contrast = c(as.character(DESeq2::design(dds))[[2]], query_group), type = shrinkingparameter, alpha = p)
            slr2 %>% as.data.frame() -> slrdf2
            slrdf2$genes <- rownames(slrdf2)
            rownames(slrdf2) <- NULL
            sigresults2 <- filtersigres(slrdf2, p, lfc)
            l2[[g]] <- list(slr2, slrdf2, sigresults2)
            names(l2[[g]]) <- c("lfcdds", "shrinked_fullres", "shrinked_sigres")
            names(l2)[g] <- gslrname
          }
          l1 <- c(l1,l2)
        }
      }
      #End of shrinking all results loop
  }

    ## Shrinking certain results
    if(tolower(ans1) == "no" || tolower(ans1) == "n"){
      res_num <- as.numeric(readline(prompt = "How many results do you want to shrink?: "))
      l1 <- vector("list", length = res_num)

      for(sresn in 1:res_num){
        if(tolower(ans2) == "no" || tolower(ans2) == "n"){
          cat("The following are the results you can shrink :) \n")
          dds_con <- DESeq2::resultsNames(dds)
          dds_con <- dds_con[-1]
          t_dds_con <- tidyname(dds_con)
          for (dc in 1:length(t_dds_con)) {cat(dc, ". ", t_dds_con[dc], "\n", sep = "")}
          c <- as.numeric(readline(prompt = "Enter the number of your desired query condition: "))
          resnlfc <- dds_con[c]
          slr <- DESeq2::lfcShrink(dds, coef = c, type = shrinkingparameter, res = DESeq2::results(dds, name = resnlfc, alpha = p))
        }

        if(tolower(ans2) == "yes" || tolower(ans2) == "y"){
          cdn <- as.character(unique(dds[[2]]))
          cat("The following are your provided experimental conditions: \n")
          for(cdni in 1:length(cdn)){cat(cdni, ": ", paste(cdn[cdni]), "\n", sep = "")}
          n1 <- as.numeric(readline(prompt = "Please enter the number of your first experimental condition query: "))
          n2 <- as.numeric(readline(prompt = "Please enter the number of your second experimental condition query: "))
          slr <- DESeq2::lfcShrink(dds, type = shrinkingparameter, contrast = c(as.character(DESeq2::design(dds))[[2]], cdn[n1], cdn[n2]), alpha = p)
          lname <- paste(cdn[n1], "vs", cdn[n2], sep = " ")
        }

        slr %>% as.data.frame() -> slrdf
        slrdf$genes <- rownames(slrdf)
        rownames(slrdf) <- NULL

        sigresults <- filtersigres(slrdf, p, lfc)
        l1[[sresn]] <- list(slr, slrdf, sigresults)
        names(l1[[sresn]]) <- c("lfcdds", "shrinked_fullres", "shrinked_sigres")

        if(tolower(ans2) == "no" || tolower(ans2) == "n"){names(l1)[sresn] <- t_dds_con[c]}
        if(tolower(ans2) == "yes" || tolower(ans2) =="y"){names(l1)[sresn] <- lname}
      }
    #End of shrinking results loop
    }

    ## Saving results
    save_df <- readline(prompt = "Do you want to save your results into individual excel files? (yes/no): ")
    if(tolower(save_df) == "yes" || tolower(save_df) == "y"){
      message("Which results do you want to save? All results(both), unfiltered results only(unfil), or significant results only(sig)?")
      save_res <- readline(prompt = "Please enter the results you want to save (both/unfil/sig): ")
      dir <- readline(prompt = "Do you want to create a new folder for the saved results? (yes/no): ")
      if (tolower(dir) == "yes" || tolower(dir) == "y"){dir_n <- createdir()}
      nlist <- as.numeric(length(l1))

      for (nl in seq_along(l1)){
        filen <- names(l1)[nl]
        cat("Now saving", filen, "\n")
        res1 <- as.data.frame(l1[[nl]][[2]])
        res2 <- as.data.frame(l1[[nl]][[3]])

        if (tolower(dir) == "yes" || tolower(dir) == "y"){
          name1 <- paste("./", dir_n, "/", filen, "_shrinked_fullres.xlsx", sep = "")
          name2 <- paste("./", dir_n, "/", filen, "_shrinked_sigres.xlsx", sep = "")
        }

        if(tolower(dir) == "no" || tolower(dir) == "n"){
          name1 <- paste("./", filen, "_shrinked_fullres.xlsx", sep = "")
          name2 <- paste("./", filen, "_shrinked_sigres.xlsx", sep = "")
        }

        if(tolower(save_res) == "both"){
          writexl::write_xlsx(res1, name1)
          writexl::write_xlsx(res2, name2)
        }

        if(tolower(save_res) == "unfil"){writexl::write_xlsx(res1, name1)}
        if(tolower(save_res) == "sig"){writexl::write_xlsx(res2, name2)}
      }
    }
    #End of saving results loop
  }

  ## Only one result to shrink
  if (x == 2){
    ans5 <- readline(prompt = "Do you want to change the shrinkage estimator? Default = apeglm (yes/no): ")
    l1 <- vector("list", length = 1)

    ## Defining Shrinking Parameter
    if(tolower(ans5) == "no" || tolower(ans5) == "n"){shrinkingparameter = "apeglm"}
    if(tolower(ans5) == "yes" || tolower(ans5) == "y") {
      suppressWarnings({
        shrinkingparameter <- readline(prompt = "Please enter your desired shrinkage estimator(ashr/normal): ")
        while(is.na(shrinkingparameter) == TRUE || shrinkingparameter != "ashr" && shrinkingparameter != "normal"){
          shrinkingparameter <- readline(prompt = "Please enter your desired shrinkage estimator(ashr/normal): ")
        }
      })
    }

    ## General shrinking functions
    slr <- DESeq2::lfcShrink(dds, coef = 2, type = shrinkingparameter, res = DESeq2::results(dds, name = DESeq2::resultsNames(dds)[2], alpha = p))
    slr %>% as.data.frame() -> slrdf
    slrdf$genes <- rownames(slrdf)
    rownames(slrdf) <- NULL

    sigresults <- filtersigres(slrdf, p, lfc)
    l1[[1]] <- list(slr, slrdf, sigresults)
    names(l1[[1]]) <- c("lfcdds", "shrinked_fullres", "shrinked_sigres")
    lname <- DESeq2::resultsNames(dds)[2]
    lname <- tidyname(lname)
    names(l1) <- lname
    res1 <- as.data.frame(l1[[1]][[2]])
    res2 <- as.data.frame(l1[[1]][[3]])
    filen <- readline(prompt = "Please provide a name for your files: ")
    name1 <- paste("./", filen, "_shrinked_fullres.xlsx", sep = "")
    name2 <- paste("./", filen, "_shrinked_sigres.xlsx", sep = "")
    writexl::write_xlsx(res1, name1)
    writexl::write_xlsx(res2, name2)
  }
  #End of one result shrink loop
  cat("Done :>")
  return(l1)
}
