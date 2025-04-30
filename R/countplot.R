#' @title Normalised Count Plots for a Single Gene
#' 
#' @description
#' Produces count plots of your DESeq2 results
#' 
#'
#' @param dds DESeqDataSet Input
#'
#' @returns Count plots of your DESeq2 results in individual .png files
#'
#' @examples
#' expression = easybioinfo::deseqexpr
#' md = easybioinfo::deseqmd
#' 
#' exampledeseq = DESeq2::DESeqDataSetFromMatrix(countData = expression, colData = md, design = ~condition, tidy = TRUE)
#' dds <- DESeq2::DESeq(exampledeseq)
#' 
#' library(easybioinfo)
#' dds <- easybioinfo::rundeseq(expression, md)
#' countplot(dds)

countplot <- function(dds){
  cs_dds <- class(dds) == "DESeqDataSet"
  if(cs_dds != "TRUE"){
    message("Please provide the DESeqDataSet from the rundeseq function")
    return(NULL)
  }

  plot_all <- readline(prompt = "Do you want to plot all the genes in the provided dataset? (yes/no): ")
  if(tolower(plot_all) == "yes" || tolower(plot_all) == "y"){genelist = rownames(SummarizedExperiment::assay(dds))}
  if(tolower(plot_all) == "no" || tolower(plot_all) == "n"){genelist <- getgenes()}

  cat("Please note that detailed parameters cannot be altered and you'll have to save it individually if you want more than 1 plot per picture \n")
  plot_n <- readline(prompt = "Do you want more than 1 graph/plot per picture? (yes/no): ")
  if(tolower(plot_n) == "yes" || tolower(plot_n) == "y"){
    cat("Please state the number of rows and the number of plots per row\n")
    cat("For example, 2 rows and 3 plots per row will produce a 2x3 picture\n")
    row_n <- as.numeric(readline(prompt = "Please state the number of rows: "))
    column_n <- as.numeric(readline(prompt = "Please state the number of plot per row: "))
    par(mfrow = c(row_n, column_n))
    for (i in unique(1:length(genelist))){
      gp <- genelist[i]
      tryCatch({
        DESeq2::plotCounts(dds, gene = gp, intgroup = as.character(DESeq2::design(dds))[[2]])
      }, error = function(e) {
        cat(gp, "not found \n")
      })
    }
    par(mfrow = c(1, 1))
  }

  if(tolower(plot_n) == "no" || tolower(plot_n) == "n"){
    def_par <- readline(prompt = "Do you want to use the default parameters? (yes/no): ")
    if(def_par == "yes" || tolower(def_par) == "y"){
      box_w <- as.numeric(0.5)
      labx <- "Condition"
      laby <- "Normalised Counts"
      t_size <- as.numeric(18)
      a_size <- as.numeric(16)
      at_size <- as.numeric(15)
    }

    if(tolower(def_par) == "no" || tolower(def_par) == "n"){
      ox_w <- as.numeric(readline(prompt = "Please enter the width of the boxplot. Default = 0.5: "))
      labx = readline(prompt = "Please enter the x-axis label: ")
      laby = readline(prompt = "Please enter the y-axis label: ")
      t_size = as.numeric(readline(prompt = "Please enter the size of the title. Default = 18: "))
      a_size = as.numeric(readline(prompt = "Please enter the size of the axis labels. Default = 16: "))
      at_size = as.numeric(readline(prompt = "Please enter the size of the axis tick labels. Default = 15: "))
    }

    dir <- readline(prompt = "Do you want to create a new directory for your saved plots? (yes/no): ")
    if(tolower(dir) == "yes" || tolower(dir) == "y"){dir_n <- createdir()}

    for (i in unique(1:length(genelist))){
      gp <- genelist[i]
      tryCatch({
        gplot <- DESeq2::plotCounts(dds, gene= gp , intgroup = as.character(DESeq2::design(dds))[[2]], returnData = TRUE)
        colnames(gplot)[2] <- "cond"
        cat("Now printing the count plot for gene", gp, "\n")
        ggraph =
          ggplot2::ggplot(gplot, aes(x = cond, y = count, colour = cond)) +
          geom_boxplot(width = box_w) +
          ggbeeswarm::geom_beeswarm(cex = 1.0, size = 1.0) +
          ggplot2::scale_y_log10(labels = scales::comma) +
          xlab(labx)+
          ylab(laby)+
          ggplot2::scale_color_brewer(palette = "Set1") +
          ggtitle(gp) +
          ggpubr::theme_pubclean() +
          theme(plot.title = element_text(hjust = 0.5, size = t_size), legend.position = 'none',
                axis.line= element_line(colour='black', linewidth = 0.8, linetype = 'solid'),
                axis.title = element_text(size = a_size), axis.text = element_text(size = at_size))
        if(tolower(dir) == "yes" || tolower(dir) == "y"){ggsave(file = paste(dir_n, "/", gp, ".png", sep = ""), width = 8, height = 8, dpi = 400)}
        if(tolower(dir) == "no" || tolower(dir) == "n"){ggsave(file = paste(gp, ".png", sep = ""), width = 8, height = 8, dpi = 400)}
        print(ggraph)
      }, error = function(e) {
        cat(gp, "not found \n")
      })
    } #End of plotting loop
  }
}
