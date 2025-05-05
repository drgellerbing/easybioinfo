#' @name maplot
#' @title Customised MA-Plot 
#' @description
#' A function to produce customised scatter plot of log2 fold changes (on the y-axis) versus the mean of normalised counts (on the x-axis)
#' @param result Either a DESeqResults object or the list produced from the getresults/shrinklfc function
#' @returns A MA plot that can be saved into .png files
#' @examples
#' expression = easybioinfo::deseqexpr
#' md = easybioinfo::deseqmd
#' 
#' exampledeseq = DESeq2::DESeqDataSetFromMatrix(countData = expression, colData = md, design = ~condition, tidy = TRUE)
#' dds <- DESeq2::DESeq(exampledeseq)
#' result <- DESeq2::results(dds)
#' maplot(result)
#' 
#' df <- easybioinfo::deseqexpr
#' md <- easybioinfo::deseqmd
#' dds <- rundeseq(df, md)
#' result <- getresults(dds)
#' maplot(result)
#' 
#' lfc <- shrinklfc(dds)
#' maplot(lfc)
#' @export

maplot <- function(result){
  if(!inherits(result, c("DESeqResults", "list"))){
    message("Please provide the DESeqResults from the getresults/shrinklfc function")
    return(NULL)
  }

  if(class(result) == "list"){
    for (m in seq_along(result)){
      cs <- class(result[[m]][[1]]) == "DESeqResults"
      if (cs != "TRUE"){
        message("Please provide the DESeqResults from the getresults/shrinklfc function")
        return(NULL)
      }
    }
    nplot <- length(result)
  }

  suppressWarnings({
    pval <- as.numeric(readline(prompt = "Please state your desired p-value cutoff: "))
    while(is.na(pval) == TRUE || pval <= 0 || pval > 1){
      message("The p-value provided should be within 0 and 1, please try again :>")
      pval <- as.numeric(readline(prompt = "Please state your desired p-value cutoff: "))
    }
  })

  if(class(result) == "list" && nplot >= 2){
    all_res <- readline(prompt = "Do you want to plot all your results? (yes/no): ")
    ans_title <- readline(prompt = "Do you want to supply your own plot titles? (yes/no): ")
    dir <- readline(prompt = "Do you want to create a new directory for your MA plots? (yes/no): ")
    if(tolower(dir) == "yes" || tolower(dir) == "y"){dir_n <- createdir()}

    if(tolower(all_res) == "yes" || tolower(all_res) == "y"){
      plotseq <- as.numeric(seq(1, length(result)))
    }

    if(tolower(all_res) == "no" || tolower(all_res) == "n"){
      message("The following are your provided results: ")
      for (np in 1:length(result)){cat(np, ": ", names(result)[np], "\n", sep = "")}

      #Input loop
      input_success <- FALSE
      while(input_success == FALSE){
        message("Please specify the result number for the plots in a comma separated format WITH ONE SPACES: e.g. 1, 2, 4")
        plotseq <- readline(prompt = "Please enter the number of the results you want to plot: ")

        #TryCatch Loop
        tryCatch({
          plotseq <- gsub("\\s*,\\s*", ",", trimws(plotseq))
          plotseq <- gsub(",", ", ", plotseq)

          if (!grepl("^[1-9]\\d*(, [1-9]\\d*)*$", plotseq)) {
            message("Invalid format. Please use comma-separated numbers (e.g., '1, 2, 4').")
            stop()
          }

          plotseq <- as.numeric(unlist(strsplit(plotseq, ", ")))

          if(any(is.na(plotseq))){stop("Error found in the input.")}

          input_success <- TRUE
        }, error = function(e){
          input_success <- FALSE
          message("Error detected. Please try again :) \n")
        })
      }
      #End of input loop
    }

    ## Plotting Loop
    for(pl in 1:length(plotseq)){
      plot_number <- plotseq[pl]
      mares <- result[[plot_number]][[1]]
      res <- as.data.frame(mares)
      limval = max(abs(res$log2FoldChange), na.rm = TRUE)
      limval <- round(limval/5)*5

      if(tolower(ans_title) == "yes" || tolower(ans_title) == "y"){
        cat("Producing MA plot for", paste(names(result)[plot_number]), "\n")
        title <- readline(prompt = "Please provide the title for the current plot: ")
      }

      if(tolower(ans_title) == "no" || tolower(ans_title) == "n"){title <- names(result)[plot_number]}

      if(tolower(dir) == "yes" || tolower(dir) == "y"){plottitle <- paste(dir_n, "/", title, "_maplot.png", sep = "")}
      if(tolower(dir) == "no" || tolower(dir) == "n"){plottitle <- paste(title, "_maplot.png", sep = "")}

      maplotfunc(mares, pval, limval, title, plottitle)
    }
    ## End of plotting loop
  }

  if(class(result) == "DESeqResults" || class(result) == "list" && nplot == 1){
    title <- readline(prompt = "Please provide the title for the MA plot: ")

    if(class(result) == "DESeqResults"){
      res <- as.data.frame(result)
      limval = max(abs(res$log2FoldChange), na.rm = TRUE)
      limval <- round(limval/5)*5
    }

    if(class(result) == "list"){
      result <- result[[1]][[1]]
      res <- as.data.frame(result)
      limval = max(abs(res$log2FoldChange), na.rm = TRUE)
      limval <- round(limval/5)*5
    }

    maplotfunc(result, pval, limval, title)
  }
}
