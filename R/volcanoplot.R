volcanoplot <- function(result){
  message("Volcano Plots will be saved in your working directory :>")
  ##Checking provided result
  if(!inherits(result, c("DESeqResults", "list", "DESeqDataSet"))){
    message("Please provide the DESeqResults from the getresults/shrinklfc function")
    return(NULL)
  }

  if(class(result) == "list"){
    for (m in seq_along(result)){
      cs <- class(result[[m]][[1]]) == "DESeqResults"
      if (cs != "TRUE"){
        message("Please provide the results from either the getresults or shrinklfc function")
        return(NULL)
      }
    }
  }

  suppressWarnings({
    p <- as.numeric(readline(prompt = "Please state your desired p-value cutoff: "))
    while(is.na(p) == TRUE || p <= 0 || p > 1){
      message("The p-value provided should be within 0 and 1, please try again :>")
      p <- as.numeric(readline(prompt = "Please state your desired p-value cutoff: "))
    }
  })

  lfc = as.numeric(readline(prompt = "Enter the desired LFC cutoff: "))

  nplot <- as.numeric(seq(1, length(result)))

  ## Getting parameters for plot
  def_par <- readline(prompt = "Do you want the default plot settings? (yes/no): ")
  if(tolower(def_par) == "no" || tolower(def_par) == "n"){
    suppressWarnings({
      ps <- as.numeric(readline(prompt = "Please enter the value for the point size. Default = 2.5: "))
      while(is.numeric(ps) == FALSE){
        message("Input provided was not a number, please try again :>")
        ps <- as.numeric(readline(prompt = "Please enter the value for the point size. Default = 2.5: "))
      }
    })
    suppressWarnings({
      ls <- as.numeric(readline(prompt = "Please enter the value for the label size. Default = 3.5: "))
      while(is.numeric(ls) == FALSE){
        message("Input provided was not a number, please try again :>")
        ls <- as.numeric(readline(prompt = "Please enter the value for the label size. Default = 3.5: "))
      }
    })
    suppressWarnings({
      lim <- as.numeric(readline(prompt = "Please enter the value for the x-axis limit. Default = 10: "))
      while(is.numeric(lim) == FALSE){
        message("Input provided was not a number, please try again :>")
        lim <- as.numeric(readline(prompt = "Please enter the value for the x-axis limit. Default = 10: "))
      }
    })
  }
  else{
    ps <- as.numeric(2.5)
    ls <- as.numeric(3)
    lim <- as.numeric(10)
  }

  if(length(nplot) >= 2){
    all_res <- readline(prompt = "Do you want to plot all your results? (yes/no): ")
    if(tolower(all_res) == "yes" || tolower(all_res) == "y"){
      cat("Your volcano plots will be printed in the following order: \n")
      for(l in 1:length(nplot)){
        np = nplot[l]
        cat(paste(np, names(result)[np], sep = ". "), sep = "\n")
      }
    }

    if(tolower(all_res) == "no" || tolower(all_res) == "n"){
      message("The following are your provided results: ")
      for (np in 1:length(result)){cat(np, ". ", names(result)[np], "\n", sep = "")}

      #Input loop
      input_success <- FALSE
      while(input_success == FALSE){
        message("Please specify the result number for the plots in a comma separated format WITH ONE SPACE: e.g. 1, 2, 4")
        splot <- readline(prompt = "Please enter the result number of the plots required: ")

        #TryCatch Loop
        tryCatch({
          splot <- gsub("\\s*,\\s*", ",", trimws(splot))
          splot <- gsub(",", ", ", splot)

          if (!grepl("^[1-9]\\d*(, [1-9]\\d*)*$", splot)) {
            message("Invalid format. Please use comma-separated numbers (e.g., '1, 2, 4').")
            stop()
          }

          nplot <- as.numeric(unlist(strsplit(splot, ", ")))

          if(any(is.na(nplot))){stop("Error found in the input.")}

          input_success <- TRUE
        }, error = function(e){
          input_success <- FALSE
          message("Error detected. Please try again :) \n")
        })
      }
      #End of input loop
    }


    p_title <- readline(prompt = "Do you want to provide your own plot titles? (yes/no): ")

    slab <- readline(prompt = "Do you want to only label specific genes? (yes/no): ")
    if(tolower(slab) == "no" || tolower(slab) == "n"){slabel = NULL}
    if(tolower(slab) == "yes" || tolower(slab) == "y"){slab2 <- readline(prompt = "Do you have to specifically label different genes for each plot? (yes/no): ")}
    if(tolower(slab) == "yes" && tolower(slab2) == "no"){slabel = getgenes()}

    for (l in 1:length(nplot)){
      i = nplot[l]
      cat(paste("Now printing plot ", i, ". ", names(result)[i], "\n", sep = ""))

      if(tolower(p_title) == "yes" || tolower(p_title) == "y"){title1 = readline(prompt = "Please provide the title for plot: ")}
      if(tolower(p_title) == "no" || tolower(p_title) == "n"){title1 = names(result)[i]}

      if(tolower(slab) == "yes" && tolower(slab2) == "yes"){slabel = getgenes()}

      res <- as.data.frame(result[[i]][[2]])
      plot <- EnhancedVolcano::EnhancedVolcano(res,
                             lab = res$genes,
                             selectLab = slabel,
                             x = 'log2FoldChange',
                             y = 'padj',
                             title = title1,
                             pCutoff = p,
                             FCcutoff = lfc,
                             pointSize = ps,
                             labSize = ls,
                             subtitle = NULL,
                             xlim = c(lim,-lim),
                             legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
      )
      print(plot)
      ggsave(paste0(title1, "_volcanoplot.png", sep = ""), width = 8, height = 8, dpi = 500)
    }
  }

  else{
    if(class(result) == "list"){cat("This is your provided result: ", paste(names(result)), "\n")}
    if(class(result) == "DESeqResults"){
      resn <- DESeq2::resultsNames(result)
      resn <- sub("^[^_]*_", "", resn)
      resn <- gsub("_", " ", resn)
      resn <- gsub("\\.", " ", resn)
      cat("This is your provided result: ", resn, "\n")
    }
    title1 = readline(prompt = "Please provide a title for the plot: ")
    if(class(result) == "DESeqResults"){alllabel = rownames(res)}
    if(class(result) == "DESeqDataSet"){res = as.data.frame(DESeq2::results(result))
                                        alllabel = rownames(res)}
    if(class(result) == "list"){
      res = result[[1]][[1]]
      alllabel = rownames(res)
    }

    slab <- readline(prompt = "Do you want to only label specific genes? (yes/no): ")
    if(tolower(slab) == "no" || tolower(slab) == "n"){slabel = NULL}
    if(tolower(slab) == "yes" || tolower(slab) == "y"){slabel <- getgenes()}

    plot <- EnhancedVolcano::EnhancedVolcano(res,
                           lab = alllabel,
                           selectLab = slabel,
                           x = 'log2FoldChange',
                           y = 'padj',
                           title = title1,
                           pCutoff = p,
                           FCcutoff = lfc,
                           pointSize = ps,
                           labSize = ls,
                           subtitle = NULL,
                           xlim = c(lim,-(lim)),
                           legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", p, sep = " "), paste("padj <", p, "& |LFC| >", lfc, sep = " "))
    )
    print(plot)
    ggsave(paste0(title1, "_volcanoplot.png", sep = ""), width = 8, height = 8, dpi = 500)
  }
  cat("All done! :> \n")
}

