#' @import magrittr
#' @import ggplot2

## General Function to grab gene list
getgenes = function(){
  message("For your query genes, do you want to enter the gene names directly (d) or provide a dataframe/vector of character strings (v)?")
  input = readline(prompt = "Enter D or V (d/v): ")
  while(tolower(input) != "d" && tolower(input) != "v"){
    message("Answer not given, please try again :)")
    input = readline(prompt = "Enter D or V (d/v): ")
  }
  if(tolower(input) == "d"){
    input_success <- FALSE
    while(input_success == FALSE){
      message("Please enter your specific genes in the following format with spaces in between: Gene1, Gene2, Gene3 ")
      genes <- readline(prompt = "Please enter a comma-separated list of gene names: ")

      #TryCatch Loop
      tryCatch({
        genes <- gsub("\\s*,\\s*", ",", trimws(genes))
        genes <- gsub(",", ", ", genes)

        gts <- unlist(strsplit(genes, ", "))

        if(any(is.na(genes))){stop("Error found in the input.")}

        input_success <- TRUE
      }, error = function(e){
        input_success <- FALSE
        message("Error detected. Please try again :) \n")
      })
    }
  }
  else if(tolower(input) == 'v'){
    message("Please provide a dataframe or a vector containing the gene names :)")
    repeat {
      gls <- readline(prompt = "Enter the name of the dataframe/vector: ")
      if (exists(gls)) {
        gts <- get(gls)
        break
      }
      message("Error: The object does not exist. Try again.")
    }

    cs_gts <- class(gts)
    cs_gts <- cs_gts[1]
    if (inherits(gts, c("data.frame", "data.table", "tbl_df"))) {
      repeat {
        genecolnum <- as.numeric(readline(prompt = "Enter the column number of the query genes: "))
        if (!is.na(genecolnum) && genecolnum > 0 && genecolnum <= ncol(gts)) break
        message("Invalid column number. Try again.")
      }
      gts <- gts[[genecolnum]]
    }
  }
  gts <- unique(stats::na.omit(gts))
  return(gts)
}

## General Function to create directories
createdir <- function(){
  message("Please make sure the name of your created folder does not exists in your current directory")
  success <- FALSE
  while(!success){
    dir_name <- readline(prompt = "Please enter the name for your newly created folder: ")
    if(dir.exists(dir_name)){
      message("Folder already exists. Please enter a different name :)")
    } else{
      dir.create(dir_name)
      success <- TRUE
    }
  }
  return(dir_name)
}

## General Function to create directories
createdir <- function(){
  message("Please make sure the name of your created folder does not exists in your current directory")
  success <- FALSE
  while(!success){
    dir_name <- readline(prompt = "Please enter the name for your newly created folder: ")
    if(dir.exists(dir_name)){
      message("Folder already exists. Please enter a different name :)")
    } else{
      dir.create(dir_name)
      success <- TRUE
    }
  }
  return(dir_name)
}
 
onegenecox = function(fsurv, gname, kmdf){
  unicox <- survival::coxph(fsurv, data = kmdf)
  ucoef = summary(unicox)$coef[1,1]
  uhaz = summary(unicox)$coef[1,2]
  upvalue = summary(unicox)$coef[1,5]
  ci = stats::confint(unicox, level = 0.95)
  ci <- exp(ci[1,])
  ci <- round(ci, digits = 2)
  uci <- paste(ci[1], '~', ci[2], sep = " ")
  unicoxdataframe = data.frame(gname, ucoef, uhaz, uci, upvalue)
  return(unicoxdataframe)
}

## Function for Plotting Survival Curves
kmplot <- function(fsurv, kmplotdf, kmgname, kmylabel = "Survival Probability", kmltitle = "Median Survival", survpvalue = NULL, kmpvalue = NULL, kmlpos = c(0.65, 0.9), kmmsize = 18, kmlsize = 16, kmpsize = 16, kmaxsize = 16, kmpxc = 0.45, kmpyc = 38){
  fit1 = survminer::surv_fit(fsurv, data = kmplotdf)
  slow = round(unname(summary(fit1)$table[,'median'][1]),0)
  shigh = round(unname(summary(fit1)$table[,'median'][2]),0)
  nlow = unname(summary(fit1)$table[,'records'][1])
  nhigh = unname(summary(fit1)$table[,'records'][2])
  pval <- survminer::surv_pvalue(fit1)
  pval <- pval[[4]]
  pval <- gsub("p = ", "", pval)
  pval <- as.numeric(pval)
  pval <- format(pval, scientific = FALSE)


  if(survpvalue < kmpvalue){
    pvalue <- bquote(bolditalic(p) == bold(.(pval)))
    custom_theme <- function() {
      survminer::theme_survminer() %+replace%
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5),legend.box = "vertical")+
        ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = kmpxc, vjust = -(kmpyc), face = "bold", size = kmpsize))}
  }

  if(survpvalue > kmpvalue || (is.null(survpvalue) && is.null(kmpvalue))){
    pvalue <- bquote(italic(p) == .(pval))
    custom_theme <- function() {
      survminer::theme_survminer() %+replace%
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5),legend.box = "vertical")+
        ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = kmpxc, vjust = -(kmpyc), size = kmpsize))}
  }


  graph = survminer::ggsurvplot(fit1, data=kmplotdf,
                                break.time.by = 500,
                                xlab = "Time (Days)", ylab = kmylabel,
                                palette = c( "steelblue2", "firebrick1"),
                                title = (paste(kmgname, " Survival Curve", sep = "")), ggtheme = custom_theme(),
                                legend.labs = c(paste("Low Expression (n= ", nlow, ")", " : ", slow, " days", sep = ""),
                                                paste("High Expression (n= ", nhigh, ")", " : ", shigh, " days", sep = "")),
                                legend.title = kmltitle, legend = kmlpos,
                                font.main = c(kmmsize,"bold"), font.legend = c(kmlsize),
                                font.x = c(kmaxsize, "bold"), font.y = c(kmaxsize, "bold"), font.tickslab = c(15)) + ggplot2::labs(subtitle = pvalue)
  print(graph)
  return(graph)
}

## Function to filter DESeq2 results
filtersigres <- function(df, pval, lfcval){
  results = df[which(df$padj < pval), ]
  resup = results[results$log2FoldChange > lfcval, ]
  resdown = results[results$log2FoldChange < -abs(lfcval), ]
  sigresults = rbind(resup, resdown)
  return(sigresults)
}

maplotfunc <- function(deseqres, pval = 0.05, limitvalue = 10, title = "MA Plot", plottitle = "./maplot.png"){
  DESeq2::plotMA(deseqres, ylim = c(-limitvalue, limitvalue),
                 xlab = "Mean of Normalised Counts",
                 ylab = "Log2 Fold Change",
                 main = title, alpha = pval)

  grDevices::png(plottitle, res = 400, units = "in", width = 8, height = 6)
  DESeq2::plotMA(deseqres, ylim = c(-limitvalue, limitvalue),
                 xlab = "Mean of Normalised Counts",
                 ylab = "Log2 Fold Change",
                 main = title, alpha = pval)
  grDevices::dev.off()
}

volcplotfunc <- function(res, alllabel, pvalue = 0.05, lfc = 1, slabel = NULL, title = "Volcano Plot", ps = 2.5, ls = 3, lim = 10){
  EnhancedVolcano::EnhancedVolcano(res,
                                   lab = alllabel,
                                   selectLab = slabel,
                                   x = 'log2FoldChange',
                                   y = 'padj',
                                   title = title,
                                   pCutoff = pvalue,
                                   FCcutoff = lfc,
                                   pointSize = ps,
                                   labSize = ls,
                                   subtitle = NULL,
                                   xlim = c(lim,-(lim)),
                                   xlab = bquote(~Log[2] ~ "Fold Change"),
                                   legendLabels = c('Insignificant', paste("|LFC| >", lfc, sep = " "), paste("padj <", pvalue, sep = " "), paste("padj <", pvalue, "& |LFC| >", lfc, sep = " "))
  )
}