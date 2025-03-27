#' @import magrittr

tidydata = function(){
  cat("Combining expression data with clinical data ....\n")

  md_success <- FALSE
  while(!md_success){
    n_md <- readline(prompt = "Please enter the name of your survival data: ")
    tryCatch({
      md <- get(n_md)
      md_success <- TRUE
    }, error = function(e){
      message("Error: Please provide a valid dataframe")
      md_success <- FALSE
    })
  }

  md <- as.data.frame(md)

  View(md)

  suppressWarnings({
    idcol <- as.numeric(readline(prompt = "Please provide the column number of the sample ID in your survival/clinical dataframe: "))
    while(is.na(idcol) == TRUE || idcol <= 0 || idcol > ncol(md)){
      message("Error detected, please try again :)")
      idcol <- as.numeric(readline(prompt = "Please provide the column number of the sample ID in your survival/clinical dataframe: "))
    }
  })
  colnames(md)[idcol] <- 'id'

  suppressWarnings({
    message("Please make sure the clinical data has all the necesssary data - Sample ID, survival times, and survival status")
    surv_type <- readline(prompt = "Are you looking to test overall survival, progression free survival, or both? (os/pfs/both): ")
    while(is.na(surv_type) == TRUE || surv_type != "os" && surv_type != "pfs" && surv_type != "both"){
      message("Error detected, please re-enter your answer :)")
      surv_type <- readline(prompt = "Are you looking into OS, PFS, or both? (OS/PFS/both): ")
    }
  })

  if(surv_type == "os"){
    suppressWarnings({
      os_col <- as.numeric(readline(prompt = "Please provide the column number of the overall survival in your survival/clinical dataframe: "))
      while(is.na(os_col) == TRUE || os_col <= 0 || os_col > ncol(md) || os_col == idcol){
        message("Error detected, please try again :)")
        os_col <- as.numeric(readline(prompt = "Please provide the column number of the overall survival in your survival/clinical dataframe: "))
      }
    })
    os_val <- readline(prompt = "Are your overall survival values in months? (yes/no): ")
    while(tolower(os_val) != "yes" && tolower(os_val) != "no"){
      message("Answer not given, please try again :)")
      os_val <- readline(prompt = "Are your overall survival values in months? (yes/no): ")
    }

    if(tolower(os_val) == "yes" || tolower(os_val) == "y"){
      cat("Converting overall survival to days and checking for NA values..\n")
      cat("Overall survival is converted to days by multiplying the values by 30\n")
      colnames(md)[os_col] <- "os_months"
      md <- reportingna_surv(md, os_col)
      days <- convertmonths(md$osmonths)
      colnames(days) <- "os"
      md <- cbind(md, days)
    }

    if(tolower(os_val) == "no" || tolower(os_val) == "n"){
      cat("Checking for NA values...\n")
      colnames(md)[os_col] <- "os"
      md <- reportingna_surv(md, os_col)
      md$os <- as.numeric(md$os)
    }

    q <- as.numeric(readline(prompt = "Please provide the column number of the overall survival status in your survival/clinical dataframe: "))
    suppressWarnings({
      while(is.na(q) == TRUE || q <= 0 || q > ncol(md) || q == idcol || q == os_col){
        message("Error detected, please try again :)")
        q <- as.numeric(readline(prompt = "Please provide the column number of the overall survival status in your survival/clinical dataframe: "))
      }
    })
    colnames(md)[q] <- "vitalos"

    md <- choosingeventos(md, q)

    md %<>% dplyr::select(id, os, vitalos, statusos)
  }

  if(surv_type == "pfs"){
    suppressWarnings({
      pfs_col <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival in your survival/clinical dataframe: "))
      while(is.na(pfs_col) == TRUE || pfs_col <= 0 || pfs_col > ncol(md) || pfs_col == idcol){
        message("Error detected, please try again :)")
        pfs_col <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival in your survival/clinical dataframe: "))
      }
    })
    pfs_val <- readline(prompt = "Are your progression free survival values in months? (yes/no): ")
    while(tolower(pfs_val) != "yes" && tolower(pfs_val) != "no"){
      message("Answer not given, please try again :)")
      pfs_val <- readline(prompt = "Are your progression free survival values in months? (yes/no): ")
    }

    if(tolower(pfs_val) == "yes" || tolower(pfs_val) == "y"){
      cat("Converting progression free survival to days and checking for NA values..\n")
      cat("PFS is converted to days by multiplying the values by 30\n")
      colnames(md)[pfs_col] <- "pfs_months"
      md <- reportingna_surv(md, pfs_col)
      days <- convertmonths(md$pfs_months)
      colnames(days) <- "pfs"
      md <- cbind(md, days)
    }

    if(tolower(pfs_val) == "no" || tolower(pfs_val) == "n"){
      cat("Checking for NA values...\n")
      colnames(md)[pfs_col] <- "pfs"
      md <- reportingna_surv(md, pfs_col)
      md$pfs <- as.numeric(md$pfs)
    }
    q <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival status in your survival/clinical dataframe: "))
    suppressWarnings({
      while(is.na(q) == TRUE || q <= 0 || q > ncol(md) || q == idcol || q == os_col){
        message("Error detected, please try again :)")
        q <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival status in your survival/clinical dataframe: "))
      }
    })
    colnames(md)[q] <- "vitalpfs"
    md <- choosingeventpfs(md, q)

    md %<>% dplyr::select(id, os, vitalpfs, statuspfs)
  }

  if(surv_type == "both"){
    suppressWarnings({
      os_col <- as.numeric(readline(prompt = "Please provide the column number of the overall survival in your survival/clinical dataframe: "))
      while(is.na(os_col) == TRUE || os_col <= 0 || os_col > ncol(md) || os_col == idcol){
        message("Error detected, please try again :)")
        os_col <- as.numeric(readline(prompt = "Please provide the column number of the overall survival in your survival/clinical dataframe: "))
      }
    })
    os_val <- readline(prompt = "Are your overall survival values in months? (yes/no): ")
    while(tolower(os_val) != "yes" && tolower(os_val) != "no"){
      message("Answer not given, please try again :)")
      os_val <- readline(prompt = "Are your overall survival values in months? (yes/no): ")
    }

    if(tolower(os_val) == "yes" || tolower(os_val) == "y"){
      cat("Converting overall survival to days and checking for NA values..\n")
      cat("Overall survival is converted to days by multiplying the values by 30\n")
      colnames(md)[os_col] <- "os_months"
      md <- reportingna_surv(md, os_col)
      days <- convertmonths(md$os_months)
      colnames(days) <- "os"
      md <- cbind(md, days)
    }

    if(tolower(os_val) == "no" || tolower(os_val) == "n"){
      cat("Checking for NA values...\n")
      colnames(md)[os_col] <- "os"
      md <- reportingna_surv(md, os_col)
      md$os <- as.numeric(md$os)
    }

    q <- as.numeric(readline(prompt = "Please provide the column number of the overall survival status in your survival/clinical dataframe: "))
    suppressWarnings({
      while(is.na(q) == TRUE || q <= 0 || q > ncol(md) || q == idcol || q == os_col){
        message("Error detected, please try again :)")
        q <- as.numeric(readline(prompt = "Please provide the column number of the overall survival status in your survival/clinical dataframe: "))
      }
    })
    colnames(md)[q] <- "vitalos"
    md <- choosingeventos(md, q)

    suppressWarnings({
      pfs_col <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival in your survival/clinical dataframe: "))
      while(is.na(pfs_col) == TRUE || pfs_col <= 0 || pfs_col > ncol(md) || pfs_col == idcol){
        message("Error detected, please try again :)")
        pfs_col <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival in your survival/clinical dataframe: "))
      }
    })
    pfs_val <- readline(prompt = "Are your progression free survival values in months? (yes/no): ")
    while(tolower(pfs_val) != "yes" && tolower(pfs_val) != "no"){
      message("Answer not given, please try again :)")
      pfs_val <- readline(prompt = "Are your progression free survival values in months? (yes/no): ")
    }

    if(tolower(pfs_val) == "yes" || tolower(pfs_val) == "y"){
      cat("Converting progression free survival to days and checking for NA values..\n")
      cat("PFS is converted to days by multiplying the values by 30\n")
      colnames(md)[pfs_col] <- "pfs_months"
      md <- reportingna_surv(md, pfs_col)
      days2 <- convertmonths(md$pfs_months)
      colnames(days2) <- "pfs"
      md <- cbind(md, days2)
    }

    if(tolower(pfs_val) == "no" || tolower(pfs_val) == "n"){
      cat("Checking for NA values...\n")
      colnames(md)[pfs_col] <- "pfs"
      md <- reportingna_surv(md, pfs_col)
      md$pfs <- as.numeric(md$pfs)
    }
    w <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival status in your survival/clinical dataframe: "))
    suppressWarnings({
      while(is.na(w) == TRUE || w <= 0 || w > ncol(md) || w == idcol || w == os_col){
        message("Error detected, please try again :)")
        w <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival status in your survival/clinical dataframe: "))
      }
    })
    colnames(md)[w] <- "vitalpfs"
    md <- choosingeventpfs(md, w)

    md %<>% dplyr::select(id, os, vitalos, statusos, pfs, vitalpfs, statuspfs)
  }

  if(all(md$os != 'NA')) {cat("Data is clean \n")}

  expr_success <- FALSE
  while(!expr_success){
    n_expr <- readline(prompt = "Please enter the name of your expression data: ")
    tryCatch({
      expressiondataframe <- get(n_expr)
      expr_success <- TRUE
    }, error = function(e){
      message("Error: Please provide a valid dataframe")
      expr_success <- FALSE
    })
  }

  expressiondataframe <- as.data.frame(expressiondataframe)
  gts <- getgenes()

  View(expressiondataframe)
  t_expr <- readline(prompt = "Does your expression data contain more than one column of gene identifiers? (yes/no): ")
  while(tolower(t_expr) != "yes" && tolower(t_expr) != "no"){
    message("Answer not given, please try again :)")
    t_expr <- readline(prompt = "Does your expression data contain more than one column of gene identifiers? (yes/no): ")
  }
  if(tolower(t_expr) == "yes" || tolower(t_expr) == "y"){
    message("If yes, columns that are not the official gene symbols have to be removed :)")
    suppressWarnings({
      hgnc_ncol <- as.numeric(readline(prompt = "Please state the column number that contains the official gene names: "))
      while(is.na(hgnc_ncol) == TRUE || hgnc_ncol <= 0 || hgnc_ncol > ncol(expressiondataframe)){
        message("Error detected, please try again :)")
        hgnc_ncol <- as.numeric(readline(prompt = "Please state the column number that contains the official gene names: "))
      }
    })
    colnames(expressiondataframe)[hgnc_ncol] <- "hgnc"
    expressiondataframe %<>% dplyr::filter(expressiondataframe$hgnc %in% gts)
    rm_genes <- gts[!gts %in% expressiondataframe$hgnc]
    gts <- expressiondataframe$hgnc
    message("If there are more than one column of gene identifiers to remove, please enter it in a comma separated format WITH spaces")
    message("Example: 1, 3, 5")
    Sys.sleep(0.5)
    message("If there is only one, just key in the specific column number")
    message("Example: 1")
    suppressWarnings({
      gene_ncol <- readline(prompt = "Please state the column number in a comma separated format with spaces: ")
      gene_ncol <- as.numeric(unlist(strsplit(gene_ncol, ", ")))
      while(is.na(gene_ncol) == TRUE){
        message("Error detected, please try again :)")
        gene_ncol <- readline(prompt = "Please state the column number in a comma separated format with spaces: ")
        gene_ncol <- as.numeric(unlist(strsplit(gene_ncol, ", ")))
      }
    })

    gene_ncol <- setdiff(gene_ncol, hgnc_ncol)  # Ensure the HGNC column isn't removed

    expressiondataframe <- expressiondataframe[, -c(gene_ncol)]
    if(length(rm_genes) > 0){
      for(i in 1:length(rm_genes)){cat(rm_genes[i], "was not found in your expression data \n")}
    }
  }

  if(tolower(t_expr) == "no" || tolower(t_expr) == "n"){
    suppressWarnings({
      gene_ncol <- as.numeric(readline(prompt = "Please state the column number that contains the official gene names: "))
      while(is.na(gene_ncol) == TRUE || gene_ncol <= 0 || gene_ncol > ncol(expressiondataframe)){
        message("Error detected, please try again :> ")
        gene_ncol <- as.numeric(readline(prompt = "Please state the column number that contains the official gene names: "))
      }
    })
    colnames(expressiondataframe)[gene_ncol] <- "hgnc"
    expressiondataframe %<>% dplyr::filter(expressiondataframe$hgnc %in% gts)
    rm_genes <- gts[!gts %in% expressiondataframe$hgnc]

    if(length(rm_genes) > 0){
      for(i in 1:length(rm_genes)){cat(rm_genes[i], "was not found in your expression data \n")}
    }
    gts <- expressiondataframe$hgnc
  }
  expressiondataframe %<>% dplyr::relocate(hgnc)
  expressiondataframe <- as.data.frame(t(expressiondataframe))
  colnames(expressiondataframe) <- expressiondataframe[1, ]
  expressiondataframe <- expressiondataframe[-1, ]
  expressiondataframe$id <- rownames(expressiondataframe)
  expressiondataframe %<>% dplyr::relocate(id)
  expressiondataframe[, -1] <- lapply(expressiondataframe[, -1], as.double)
  expressiondataframe <- as.data.frame(expressiondataframe)
  expr2 <- expressiondataframe[(expressiondataframe$id %in% md$id), ]
  rm_expr <- setdiff(expressiondataframe$id, expr2$id)

  tryCatch({
    if (!all(expr2$id %in% md$id) || !all(md$id %in% expr2$id)) {
      stop("Samples of the expression and metadata aren't equal")
    }
  }, error = function(e) {
    message(e$message)
  })

  if(length(rm_expr) <= 10){
    for(i in 1:length(rm_expr)){
      x = rm_expr[i]
      cat("Sample", x, "is not found in your clinical data \n")
    }
  } else{
    cat("A total of", length(rm_expr), "samples were not found in your clinical data \n")
  }


  if(all(md$id %in% expr2$id == FALSE)){
    rm_md <- setdiff(md$id, expr2$id)
    for(n in 1:length(rm_md)){
      mdl = rm_md[n]
      cat("Sample", mdl, "is not found in your expression data \n")
    }
    md %<>% dplyr::filter(md$id %in% expr2$id)
  }

  if(any(duplicated(colnames(expr2))) == TRUE){
    dupgenes <- colnames(expr2)[duplicated(colnames(expr2))]
    for(dupn in 1:length(dupgenes)){
      message("Please note ", dupgenes[dupn] , " is duplicated in your expression data")
    }
    cat("It is okay to proceed, but do double check your expression data \n")
  }

  md <- as.data.frame(md)
  md <- merge(md, expr2, by = "id")
  if(any(duplicated(colnames(expr2))) == FALSE){
    cat("Congrats :> It is safe to proceed :> \n")
  }
  rownames(md) <- NULL
  return(md)
}
