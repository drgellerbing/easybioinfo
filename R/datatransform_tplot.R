#' @import magrittr

#' @name datatransform
#' 
#' @title Transforms and Normalises DESeq2 Counts for Downstream Analysis
#' 
#' @description
#' Transforms your DESeq2 counts via:
#' 1. Variance Stabilising Transformation
#' 2. Regularised Log Transformation
#' 3. Normalised Count Transformation
#' 
#' Please check the DESeq2 vignette for further details on the transformations: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#' 
#' @param dds DESeqDataSet Input
#'
#' @returns A DESeqTransform object of transformed, normalised counts and saves them into excel files
#'
#' @examples
#' expression = easybioinfo::deseqexpr
#' md = easybioinfo::deseqmd
#' 
#' exampledeseq = DESeq2::DESeqDataSetFromMatrix(countData = expression, colData = md, design = ~condition, tidy = TRUE)
#' dds <- DESeq2::DESeq(exampledeseq)
#' 
#' library(easybioinfo)
#' dds <- easybioinfo::rundeseq(expressiondf, md)
#' datatransform(dds)
#' @export

datatransform <- function(dds){
  cs_dds <- class(dds) == "DESeqDataSet"
  if(cs_dds != "TRUE"){
    message("Please provide the DESeqDataSet from the rundeseq function")
    return(NULL)
  }

  message("There are three types of transformation for downstream visualisations:")
  cat("1. Variance Stabilising Transformation (vst)\n")
  cat("2. Normalised Count Transformation (nct)\n")
  cat("3. Regularised Log Transformation(rlog)\n")
  message("Please check the DESeq2 vignette below for further details on the transformations: ")
  cat("http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html\n")

  suppressWarnings({
    tmtd <- readline(prompt = "Please state the desired method of transformation (vst/nct/rlog): ")
    while(is.na(tmtd) == TRUE || tmtd != "vst" && tmtd != "nct" && tmtd != "rlog"){
      message("Error detected, please try again :>")
      tmtd <- readline(prompt = "Please state the desired method of transformation (vst/nct/rlog): ")
    }
  })

  if(tolower(tmtd) == "vst"){
    cat("Blinding should be false for downstream analysis. Blinding should be true when used for comparing samples in a manner unbiased by prior information on samples, for sample to perform sample QA")
    suppressWarnings({
      bt = readline(prompt = "Do you want the blind the transformation to the experimental design? (TRUE/FALSE): ")
      while(is.na(bt) == TRUE || bt != "TRUE" && bt != "FALSE"){
        bt = readline(prompt = "Do you want the blind the transformation to the experimental design? (TRUE/FALSE): ")
      }
    })
    bt %<>% as.logical()
    tdds <- DESeq2::vst(dds, blind = bt)
    xlsxn <- "_vst.xlsx"
  }

  if(tolower(tmtd) == "nct"){
    cat("Parameters for nct include function and pseudocount..\n")
    cat("Default function: log2\n")
    cat("Default pseudocount: 1\n")
    pcount <- as.numeric(readline(prompt = "Please state the desired pseudocount. Default = 1: "))
    tdds <- DESeq2::normTransform(dds, f = log2, pc = pcount)
    xlsxn <- "_nct.xlsx"
  }

  if(tolower(tmtd) == "rlog"){
    cat("Blinding should be false for downstream analysis. Blinding should be true when used for comparing samples in a manner unbiased by prior information on samples, for sample to perform sample QA")
    suppressWarnings({
      bt = readline(prompt = "Do you want the blind the transformation to the experimental design? (TRUE/FALSE): ")
      while(is.na(bt) == TRUE || bt != "TRUE" && bt != "FALSE"){
        bt = readline(prompt = "Do you want the blind the transformation to the experimental design? (TRUE/FALSE): ")
      }
    })
    bt %<>% as.logical()
    tdds <- DESeq2::rlog(dds, blind = bt)
    xlsxn <- "_rlog.xlsx"
  }

  tres <- SummarizedExperiment::assay(tdds) %>% as.data.frame()
  tres$genes <- rownames(tdds)
  tres %<>% dplyr::relocate(genes)
  rownames(tres) <- NULL
  filen <- readline(prompt = "Please provide the name for your file: ")
  writexl::write_xlsx(tres, paste("./", filen, xlsxn, sep = ""))

  return(tdds)
}

#' @name pcaplot
#' 
#' @title PCA Plotting Function
#' @description
#' Produces a PCA plot based on your DESeqTransform dataset. Plot parameters can be tweaked specifically by the user
#' 
#' 
#' @param tdds DESeqTransform Input
#'
#' @returns A customised PCA plot which is saved into a .png file
#'
#' @examples
#' library(DESeq2)
#' dds <- DESeq2::DESeq(exampledeseq)
#' tdds <- vst(dds)
#' tdds <- rlog(dds)
#' tdds <- normTransform(dds, f = log2, pc = pcount)
#' 
#' library(easybioinfo)
#' dds <- easybioinfo::rundeseq(expressiondf, md)
#' tdds <- easybioinfo::datatransform(dds)
#' 
#' pcaplot(dds)
#' @export

pcaplot <- function(tdds){
  t <- readline(prompt = "Please provide the title of the plot: ")
  defpar <- readline(prompt = "Do you want to use the default parameters? (yes/no): ")

  ## Default Parameters
  if(tolower(defpar) == "yes" || tolower(defpar) == "y"){
    tsize = as.numeric(18)
    ltextsize = as.numeric(12)
    ltitlesize = as.numeric(14)
    lposition = "bottom"
    atextsize = as.numeric(10)
    atitlesize = as.numeric(12)
  }

  ## User-generated Parameters
  if(tolower(defpar) == "no" || tolower(defpar) == "n"){
    tsize = as.numeric(readline(prompt = "Please provide the text size for your title (Default = 18): "))
    ltitlesize = as.numeric(readline(prompt = "Please provide the text size for your legend title (Default = 14): "))
    ltextsize = as.numeric(readline(prompt = "Please provide the text size for your legend texts (Default = 12): "))
    lposition = readline(prompt = "Please provide your desired position for the legend (top/right/bottom/left): ")
    atitlesize = as.numeric(readline(prompt = "Please provide the text size for your axis titles (Default = 10): "))
    atextsize = as.numeric(readline(prompt = "Please provide the text size for your axis texts (Default = 10): "))
  }

  intgroupn <- colnames(SummarizedExperiment::colData(tdds))[2]
  pcaplot <- DESeq2::plotPCA(tdds, intgroup = intgroupn, returnData = TRUE)
  percentVar <- round(100 * attr(pcaplot, "percentVar"))
  graph <- ggplot2::ggplot(pcaplot, aes(PC1, PC2, color=pcaplot[, 4])) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    theme_minimal() +
    scale_fill_brewer() +
    labs(color = "Condition") +
    ggtitle(t) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = tsize),
          legend.text = element_text(size = ltextsize),
          legend.title = element_text(size = ltitlesize),
          legend.position = lposition,
          axis.text = element_text(size = atextsize),
          axis.title = element_text(size = atitlesize))
  print(graph)
  ggsave(paste(t, ".png", sep = ""), plot = graph, dpi = 400, width = 10, height = 10)
}

distheatmapplot <- function(tdds){
  sampleDists = stats::dist(t(SummarizedExperiment::assay(tdds)))
  samplematrix = as.matrix(sampleDists)
  rownames(samplematrix) <- paste(tdds$cond, sep="-")
  colnames(samplematrix) <- NULL
  colors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9,"Blues"))) (255)
  graph <- pheatmap::pheatmap(samplematrix,
                              clustering_distance_rows = sampleDists,
                              clustering_distance_cols=sampleDists,
                              col=colors)
  print(graph)
}

#' @name countheatmapplot
#' @title Count Heatmap Plot of Specific Genes
#' @description
#' Produces a customised heatmap of transformed/normalised counts of genes of your choosing through a series of prompts
#' 
#' @param tdds DESeqTransform Input
#' @returns A customised count heatmap plot which is saved into a .png file
#' @examples
#' library(DESeq2)
#' dds <- DESeq2::DESeq(exampledeseq)
#' tdds <- DESeq2::vst(dds)
#' tdds <- DESeq2::rlog(dds)
#' tdds <- DESeq2::normTransform(dds, f = log2, pc = pcount)
#' 
#' library(easybioinfo)
#' dds <- easybioinfo::rundeseq(expressiondf, md)
#' tdds <- easybioinfo::datatransform(dds)
#' 
#' countheatmapplot(dds)
#' @export

countheatmapplot <- function(tdds){
  message("The heatmap will be saved as countheatmap.png in your working directory :>")
  cs_dds <- class(tdds) == "DESeqTransform"
  if(cs_dds != "TRUE"){
    message("Please provide the DESeqTransform from the datatransform function")
    return(NULL)
  }

  genelist <- getgenes()

  llab <- "Condition"

  suppressWarnings({
    scn <- as.logical(readline(prompt = "Do you want to show the sample names? (TRUE/FALSE): "))
    while(is.na(scn) == TRUE || scn != "TRUE" && scn != "FALSE"){
      scn <- as.logical(readline(prompt = "Do you want to show the sample names? (TRUE/FALSE): "))
    }
  })

  suppressWarnings({
    ccol <- as.logical(readline(prompt = "Do you want to show the clustering between the samples? (TRUE/FALSE): "))
    while(is.na(ccol) == TRUE || ccol != "TRUE" && ccol != "FALSE"){
      ccol <- as.logical(readline(prompt = "Do you want to show the clustering between the samples? (TRUE/FALSE): "))
    }
  })

  sample = SummarizedExperiment::assay(tdds)[genelist,]
  sample <- (sample - rowMeans(sample))/MatrixGenerics::rowSds(sample)
  collabel <- as.data.frame(SummarizedExperiment::colData(tdds)[2])
  colnames(collabel)[1] <- llab
  colour <- grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(50)
  paletteLength <- 50
  myBreaks <- c(seq(min(sample), 0, length.out=ceiling(paletteLength/2) + 1),
                seq(max(sample)/paletteLength, max(sample), length.out=floor(paletteLength/2)))
  graph <- pheatmap::pheatmap(sample, annotation_col = collabel, show_colnames = scn, cluster_cols = ccol, labels_row = genelist, color = colour, breaks = myBreaks)
  print(graph)
  ggsave("countheatmap.png", plot = graph$gtable, dpi = 400)
}

