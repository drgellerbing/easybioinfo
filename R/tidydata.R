#' @import magrittr

#' @name tidykmdata
#' @title Data Wrangling for Downstream Survival Analysis
#' @description
#' A prompt-based function that produces a dataframe containing the overall survival or progression-free survival information for downstream survival analysis.
#' The function will prompt you for a list of genes and wrangle your expression dataframes and clinical data into one final dataframe.
#' 
#' @param expressiondataframe Dataframe containing expression values. The dataframe should have individual genes in each row and individual samples in each column.
#' @param clinicaldataframe Dataframe containing clinical data. The clinical data should contain the samples, overall survival/progression-free survival values, as well as the events for either the overall survival/progression-free survival.
#' @returns A dataframe containing information for the samples, the survival information and the expression of your query genes for downstream survival analysis using uniKM or multicox functions
#' @examples
#' df <- easybioinfo::kmexpr
#' md <- easybioinfo::kmclinical
#' 
#' kmdf <- tidykmdata(df, md)

tidykmdata = function(expressiondataframe, clinicaldataframe){
  cat("Combining expression data with clinical data ....\n")
  
  ## Repeating functions
  convertmonths = function(mdm){
    fd = function(x){x*30}
    days <- as.data.frame(sapply(as.double(mdm), fd))
    return(days)
  }
  
  ##Reporting NA in survival dataframe
  reportingna_surv = function(md, mdcol){
    mdna <- md[is.na(md[, mdcol]), ]
    if(length(mdna$id) > 0){
      for(poi in 1:length(mdna$id)){
        cat("NA values found in", mdna$id[poi], "\n")
      }
    }
    md <- md[!(is.na(md[, mdcol])), ]
    return(md)
  }
  
  ##Function to segregate patients based on desired OS events
  choosingevent = function(md, qcol, surv_type){
    if(any(is.na(md[, qcol])) == TRUE){
      vsna <- md$id[is.na(md[, qcol])]
      for(nai in 1:length(vsna)){
        cat(vsna[nai], "contains NA value\n")
      }
    }
    md <- md[!(is.na(md[, qcol])), ]
    
    vsvector <- unique(md[,qcol])
    message("Here are the events in the provided clinical/metadata:")
    for(vsvi in 1:length(vsvector)){
      cat(vsvi, ": ", vsvector[vsvi], "\n", sep = "")
    }
    vsp <- as.numeric(readline(prompt = "Please provide the number of your desired event: "))
    vs <- vsvector[vsp]
    if(surv_type == "os"){name = "statusos"}
    if(surv_type == "pfs"){name = "statuspfs"}
    
    md %>%
      dplyr::mutate("{name}" := as.numeric(
        ifelse(
          md[, qcol] == vs,
          '1',
          '0'
        ))) %>% as.data.frame() -> md
    return(md)
  }

  clinicaldataframe <- as.data.frame(clinicaldataframe)
  expressiondataframe <- as.data.frame(expressiondataframe)

  View(clinicaldataframe)

  suppressWarnings({
    idcol <- as.numeric(readline(prompt = "Please provide the column number of the sample ID in your survival/clinical dataframe: "))
    while(is.na(idcol) == TRUE || idcol <= 0 || idcol > ncol(clinicaldataframe)){
      message("Error detected, please try again :)")
      idcol <- as.numeric(readline(prompt = "Please provide the column number of the sample ID in your survival/clinical dataframe: "))
    }
  })
  colnames(clinicaldataframe)[idcol] <- 'id'

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
      while(is.na(os_col) == TRUE || os_col <= 0 || os_col > ncol(clinicaldataframe) || os_col == idcol){
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
      colnames(clinicaldataframe)[os_col] <- "os_months"
      clinicaldataframe <- reportingna_surv(clinicaldataframe, os_col)
      days <- convertmonths(clinicaldataframe$osmonths)
      colnames(days) <- "os"
      clinicaldataframe <- cbind(clinicaldataframe, days)
    }

    if(tolower(os_val) == "no" || tolower(os_val) == "n"){
      cat("Checking for NA values...\n")
      colnames(clinicaldataframe)[os_col] <- "os"
      clinicaldataframe <- reportingna_surv(clinicaldataframe, os_col)
      clinicaldataframe$os <- as.numeric(clinicaldataframe$os)
    }

    q <- as.numeric(readline(prompt = "Please provide the column number of the overall survival status in your survival/clinical dataframe: "))
    suppressWarnings({
      while(is.na(q) == TRUE || q <= 0 || q > ncol(clinicaldataframe) || q == idcol || q == os_col){
        message("Error detected, please try again :)")
        q <- as.numeric(readline(prompt = "Please provide the column number of the overall survival status in your survival/clinical dataframe: "))
      }
    })
    colnames(clinicaldataframe)[q] <- "vitalos"

    clinicaldataframe <- choosingevent(clinicaldataframe, q, surv_type)

    clinicaldataframe %<>% dplyr::select(id, os, vitalos, statusos)
  }

  if(surv_type == "pfs"){
    suppressWarnings({
      pfs_col <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival in your survival/clinical dataframe: "))
      while(is.na(pfs_col) == TRUE || pfs_col <= 0 || pfs_col > ncol(clinicaldataframe) || pfs_col == idcol){
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
      colnames(clinicaldataframe)[pfs_col] <- "pfs_months"
      clinicaldataframe <- reportingna_surv(clinicaldataframe, pfs_col)
      days <- convertmonths(clinicaldataframe$pfs_months)
      colnames(days) <- "pfs"
      clinicaldataframe <- cbind(clinicaldataframe, days)
    }

    if(tolower(pfs_val) == "no" || tolower(pfs_val) == "n"){
      cat("Checking for NA values...\n")
      colnames(clinicaldataframe)[pfs_col] <- "pfs"
      clinicaldataframe <- reportingna_surv(clinicaldataframe, pfs_col)
      clinicaldataframe$pfs <- as.numeric(clinicaldataframe$pfs)
    }
    q <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival status in your survival/clinical dataframe: "))
    suppressWarnings({
      while(is.na(q) == TRUE || q <= 0 || q > ncol(clinicaldataframe) || q == idcol || q == os_col){
        message("Error detected, please try again :)")
        q <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival status in your survival/clinical dataframe: "))
      }
    })
    colnames(clinicaldataframe)[q] <- "vitalpfs"
    clinicaldataframe <- choosingevent(clinicaldataframe, q, surv_type)

    clinicaldataframe %<>% dplyr::select(id, os, vitalpfs, statuspfs)
  }

  if(surv_type == "both"){
    suppressWarnings({
      os_col <- as.numeric(readline(prompt = "Please provide the column number of the overall survival in your survival/clinical dataframe: "))
      while(is.na(os_col) == TRUE || os_col <= 0 || os_col > ncol(clinicaldataframe) || os_col == idcol){
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
      colnames(clinicaldataframe)[os_col] <- "os_months"
      clinicaldataframe <- reportingna_surv(clinicaldataframe, os_col)
      days <- convertmonths(clinicaldataframe$os_months)
      colnames(days) <- "os"
      clinicaldataframe <- cbind(clinicaldataframe, days)
    }

    if(tolower(os_val) == "no" || tolower(os_val) == "n"){
      cat("Checking for NA values...\n")
      colnames(clinicaldataframe)[os_col] <- "os"
      clinicaldataframe <- reportingna_surv(clinicaldataframe, os_col)
      clinicaldataframe$os <- as.numeric(clinicaldataframe$os)
    }

    q <- as.numeric(readline(prompt = "Please provide the column number of the overall survival status in your survival/clinical dataframe: "))
    suppressWarnings({
      while(is.na(q) == TRUE || q <= 0 || q > ncol(clinicaldataframe) || q == idcol || q == os_col){
        message("Error detected, please try again :)")
        q <- as.numeric(readline(prompt = "Please provide the column number of the overall survival status in your survival/clinical dataframe: "))
      }
    })
    colnames(clinicaldataframe)[q] <- "vitalos"
    clinicaldataframe <- choosingevent(clinicaldataframe, q, "os")

    suppressWarnings({
      pfs_col <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival in your survival/clinical dataframe: "))
      while(is.na(pfs_col) == TRUE || pfs_col <= 0 || pfs_col > ncol(clinicaldataframe) || pfs_col == idcol){
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
      colnames(clinicaldataframe)[pfs_col] <- "pfs_months"
      clinicaldataframe <- reportingna_surv(clinicaldataframe, pfs_col)
      days2 <- convertmonths(clinicaldataframe$pfs_months)
      colnames(days2) <- "pfs"
      clinicaldataframe <- cbind(clinicaldataframe, days2)
    }

    if(tolower(pfs_val) == "no" || tolower(pfs_val) == "n"){
      cat("Checking for NA values...\n")
      colnames(clinicaldataframe)[pfs_col] <- "pfs"
      clinicaldataframe <- reportingna_surv(clinicaldataframe, pfs_col)
      clinicaldataframe$pfs <- as.numeric(clinicaldataframe$pfs)
    }
    w <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival status in your survival/clinical dataframe: "))
    suppressWarnings({
      while(is.na(w) == TRUE || w <= 0 || w > ncol(clinicaldataframe) || w == idcol || w == os_col){
        message("Error detected, please try again :)")
        w <- as.numeric(readline(prompt = "Please provide the column number of the progression free survival status in your survival/clinical dataframe: "))
      }
    })
    colnames(clinicaldataframe)[w] <- "vitalpfs"
    clinicaldataframe <- choosingevent(clinicaldataframe, w, "pfs")

    clinicaldataframe %<>% dplyr::select(id, os, vitalos, statusos, pfs, vitalpfs, statuspfs)
  }

  if(all(clinicaldataframe$os != 'NA')) {cat("Data is clean \n")}

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
  expr2 <- expressiondataframe[(expressiondataframe$id %in% clinicaldataframe$id), ]
  rm_expr <- setdiff(expressiondataframe$id, expr2$id)

  tryCatch({
    if (!all(expr2$id %in% clinicaldataframe$id) || !all(clinicaldataframe$id %in% expr2$id)) {
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


  if(all(clinicaldataframe$id %in% expr2$id == FALSE)){
    rm_md <- setdiff(clinicaldataframe$id, expr2$id)
    for(n in 1:length(rm_md)){
      mdl = rm_md[n]
      cat("Sample", mdl, "is not found in your expression data \n")
    }
    clinicaldataframe %<>% dplyr::filter(clinicaldataframe$id %in% expr2$id)
  }

  if(any(duplicated(colnames(expr2))) == TRUE){
    dupgenes <- colnames(expr2)[duplicated(colnames(expr2))]
    for(dupn in 1:length(dupgenes)){
      message("Please note ", dupgenes[dupn] , " is duplicated in your expression data")
    }
    cat("It is okay to proceed, but do double check your expression data \n")
  }

  clinicaldataframe <- as.data.frame(clinicaldataframe)
  clinicaldataframe <- merge(clinicaldataframe, expr2, by = "id")
  if(any(duplicated(colnames(expr2))) == FALSE){
    cat("Congrats :> It is safe to proceed :> \n")
  }
  rownames(clinicaldataframe) <- NULL
  return(clinicaldataframe)
}

#' @name transformexpr
#' @title Transforms RNA-Seq Expression into Raw Counts
#' @description
#' A function to convert log2 transformed RNA-Seq counts into raw counts for downstream DESeq2 differential expression analysis.
#' This function will also account for any psudocounts based on your input prompt.
#' @param expr An expression dataframe containing RNA-Seq log2 transformed counts
#' @returns A dataframe with transformed raw counts for downstream DESeq2 differential expression analysis.
#' @examples
#' df <- transformexpr(expr)
#' 

transformexpr <- function(expr){
  print("Transformation Ongoing...")
  suppressWarnings({
    gene_col <- as.numeric(readline(prompt = "Please enter the column number of your gene id in your expression dataframe: "))
    while(is.na(gene_col) == TRUE || gene_col > ncol(expr)){
      message("Error detected, please try again :)")
      gene_col <- as.numeric(readline(prompt = "Please enter the column number of your gene id in your expression dataframe: "))
    }
  })
  colnames(expr)[gene_col] <- "geneid"
  expr %<>% dplyr::relocate(geneid)
  fun <- function(x){2^x}
  expr[, -1] <- lapply(expr[, -1], fun)
  pseudo <- as.numeric(readline(prompt = "Pseudocount value: "))
  pseudofun <- function(x){as.integer(x-pseudo)}
  expr[, -1] <- lapply(expr[, -1], pseudofun)
  expr %<>% as.data.frame()
  return(expr)
  cat("Done :>")
}
