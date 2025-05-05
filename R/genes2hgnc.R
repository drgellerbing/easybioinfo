#' @name genes2hgnc
#' @title Converts Gene Identifiers to HGNC Symbols
#' @description
#' Using this function, you can easily convert different gene identifiers from your expression dataframe into official HGNC symbols. Conversion is done using the 'biomaRt' package.
#' 
#' @param expr Your expression dataframe
#' @param srcType Gene Identifier required by biomaRt. Defaults to "ensembl_gene_id". Can be changed in the function through user prompts
#' @returns A dataframe with converted gene identifiers
#' @examples
#' df <- easybioinfo::kmexpr
#' df <- genes2hgnc(df) 
#' @export

genes2hgnc <- function(expr, srcType = "ensembl_gene_id" )
{
  message("Please note that unmapped IDs are left as it is")
  ## Retrieve the EMSEMBL -> HUGO mapping
  expr <- as.data.frame(expr)
  View(expr)
  expr_ncol <- as.numeric(readline(prompt = "Please enter the column number of the gene IDs: "))


  #Loop to connect to Ensembl's database
  max_attempts <- 5
  attempt <- 1
  base_wait <- 2  # Initial wait time (seconds)
  ensembl <- NULL

  while (attempt <= max_attempts && is.null(ensembl)) {
    tryCatch({
      message("Attempt ", attempt, " to connect to Ensembl's database...")
      ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = "https://www.ensembl.org", dataset = "hsapiens_gene_ensembl")
      message("Success! :D")
    }, error = function(e) {
      wait_time <- base_wait * (2 ^ (attempt - 1))  # Exponential backoff
      message("Attempt ", attempt, " failed. Retrying in ", wait_time, " seconds...")
      Sys.sleep(wait_time)
      attempt <<- attempt + 1
    })
  }

  if (is.null(ensembl)) {
    stop("Failed after ", max_attempts, " attempts.")
  }

  #Grabbing specific attribute names
  src <- readline(prompt = "Are the gene names provided in the expression data Ensembl Gene IDs? (yes/no): ")
  if(tolower(src) == "yes" || tolower(src) == "y"){
    expr[,expr_ncol] <- gsub("\\.\\d+$", "", expr[[expr_ncol]])
    srcType = "ensembl_gene_id"
  }

  if(tolower(src) == "no" || tolower(src) == "n"){
    message("Please find the specific attribute in which to map your genes to")
    attrlist = biomaRt::listFilters(ensembl)
    View(attrlist)
    srcType = readline(prompt = "Please enter the specific attribute name in attrlist: ")

    #Loop to ensure the correct attribute name is provided
    while(!srcType %in% attrlist$name){
      message("Attribute name not found. Please ensure there are no additional SPACES/empty characters :)")
      srcType = readline(prompt = "Please enter the specific attribute name in attrlist: ")
    }
  }

  v <- expr[,expr_ncol]
  max_attempts2 <- 5
  attempt2 <- 1
  ID <- data.frame(matrix(nrow = 0, ncol = 2))

  ##Loop to make sure mapping is successful
  while(attempt2 <= max_attempts2 && nrow(ID) == 0){
    tryCatch({
      cat("Attempting to map genes....\n")
      ID <- biomaRt::getBM(attributes=c(srcType, "hgnc_symbol"), filters=srcType, values=v, mart=ensembl)
      if (nrow(ID) == 0) {
        stop("No IDs mapped successfully \n")
      }
    }, error = function(e){
      cat("No IDs mapped successfully \n")
      message("Attempt ", attempt2, " failed. Retrying....Please wait")
      attempt2 <<- attempt2 + 1
      ID <- data.frame(matrix(nrow = 0, ncol = 2))
    })
  }

  ## Drop empty duds
  j <- which( ID[,2] == "" )
  if( length(j) > 0 ) ID <- ID[-j,]
  idx <- match(expr[[expr_ncol]], ID[[1]])
  colnames(expr)[expr_ncol] <- "hgnc"
  expr$hgnc[!is.na(idx)] <- ID$hgnc_symbol[idx[!is.na(idx)]]
  cat("All done! :D \n")
  return(expr)
}

#' @name ddsgenes2hgnc
#' @title Gene Conversion of DESeq2 Results
#' @description
#' Converts your DESeq2 result into official HGNC symbols
#' Conversion is done using the 'biomaRt' package.
#' @param dds A DESeqDataSet input
#' @param srcType Gene Identifier required by biomaRt. Defaults to "ensembl_gene_id". Can be changed in the function through user prompts
#' @returns A DESeqDataSet with gene symbols converted into official HGNC symbols
#' @examples
#' library(DESeq2)
#' expression = easybioinfo::deseqexpr
#' md = easybioinfo::deseqmd
#' 
#' exampledeseq = DESeq2::DESeqDataSetFromMatrix(countData = expression, colData = md, design = ~condition, tidy = TRUE)
#' dds <- DESeq2::DESeq(exampledeseq)
#' 
#' library(easybioinfo)
#' dds <- easybioinfo::rundeseq(expressiondf, md)
#' ddsgenestohgnc(dds)
#' @export

ddsgenestohgnc <- function(dds, srcType = "ensembl_gene_id"){
  cs_dds <- class(dds) == "DESeqDataSet"
  if(cs_dds != "TRUE"){
    message("Please provide the DESeqDataSet from the rundeseq function")
    return(NULL)
  }

  ##Loop to connect to Ensembl's daatabase
  max_attempts <- 5
  attempt <- 1
  base_wait <- 2  # Initial wait time (seconds)
  ensembl <- NULL

  while (attempt <= max_attempts && is.null(ensembl)) {
    tryCatch({
      message("Attempt ", attempt, " to connect to Ensembl's database...")
      ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = "https://www.ensembl.org", dataset = "hsapiens_gene_ensembl")
      message("Success! :D")
    }, error = function(e) {
      wait_time <- base_wait * (2 ^ (attempt - 1))  # Exponential backoff
      message("Attempt ", attempt, " failed. Retrying in ", wait_time, " seconds...")
      Sys.sleep(wait_time)
      attempt <<- attempt + 1
    })
  }

  if (is.null(ensembl)) {
    stop("Failed after ", max_attempts, " attempts.")
  }

  ##Determining the correct ID to map
  src <- readline(prompt = "Are the gene names provided in the expression data Ensembl Gene IDs? (yes/no): ")
  if(tolower(src) == "yes" || tolower(src) == "y"){
    rownames(dds) <- gsub("\\.\\d+$", "", rownames(dds))
    srcType = "ensembl_gene_id"
  }

  if(tolower(src) == "no" || tolower(src) == "n"){
    message("Please find the specific attribute in which to map your genes to")
    attrlist = biomaRt::listFilters(ensembl)
    View(attrlist)
    srcType = readline(prompt = "Please enter the specific attribute name in attrlist: ")

    #Loop to ensure the correct attribute name is provided
    while(!srcType %in% attrlist$name){
      message("Attribute name not found. Please ensure there are no additional SPACES/empty characters :)")
      srcType = readline(prompt = "Please enter the specific attribute name in attrlist: ")
    }
  }

  v <- rownames(dds)
  max_attempts2 <- 5
  attempt2 <- 1
  ID <- data.frame(matrix(nrow = 0, ncol = 2))

  ##Loop to make sure mapping is successful
  while(attempt2 <= max_attempts2 && nrow(ID) == 0){
    tryCatch({
      cat("Attempting to map genes....\n")
      ID <- biomaRt::getBM(attributes=c(srcType, "hgnc_symbol"), filters=srcType, values=v, mart=ensembl)
      if (nrow(ID) == 0) {
        stop("No IDs mapped successfully \n")
      }
    }, error = function(e){
      cat("No IDs mapped successfully \n")
      message("Attempt ", attempt2, " failed. Retrying....Please wait")
      attempt2 <<- attempt2 + 1
      ID <- data.frame(matrix(nrow = 0, ncol = 2))
    })
  }

  ## Drop empty duds
  j <- which( ID[,2] == "" )
  if( length(j) > 0 ) ID <- ID[-j,]
  idx <- match(rownames(dds), ID[[1]])
  rownames(dds)[!is.na(idx)] <- ID$hgnc_symbol[idx[!is.na(idx)]]
  cat("All done! :D \n")
  return(dds)
}
