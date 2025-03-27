#' @import magrittr
#' @import ggplot2

##Function to grab gene list
getgenes = function(){
  message("For your query genes, do you want to enter the gene names directly (d) or provide a dataframe/vector of character strings (v)?")
  input = readline(prompt = "Enter D or V (d/v): ")
  while(tolower(input) != "d" && tolower(input) != "v"){
    message("Answer not given, please try again :)")
    input = readline(prompt = "Enter D or V (d/v): ")
  }
  if(tolower(input) == "d"){
    message("Please enter your specific genes in the following format with spaces in between: Gene1, Gene2, Gene3 ")
    genes <- readline(prompt = "Please enter a comma-separated list of gene names: ")
    gts <- unlist(strsplit(genes, ", "))
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

convertmonths = function(mdm){
  fd = function(x){x*30}
  days <- as.data.frame(sapply(as.double(mdm), fd))
  return(days)
}

##Function to segregate patients based on desired OS events
choosingeventos = function(md, qcol){
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

  md %>%
    dplyr::mutate(statusos = as.numeric(
      ifelse(
        md[, qcol] == vs,
        '1',
        '0'
      ))) %>% as.data.frame() -> md
  return(md)
}

##Function to segregate patients based on desired PFS events
choosingeventpfs = function(md, qcol){
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

  md %>%
    dplyr::mutate(statuspfs = as.numeric(
      ifelse(
        md[, qcol] == vs,
        '1',
        '0'
      ))) %>% as.data.frame() -> md
  return(md)
}

## Function for Univariate Cox
onegenecox = function(fsurv, gname, ogcdf){
  unicox <- survival::coxph(fsurv, data = ogcdf)
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

kmplot <- function(fsurv, kmplotdf, kmgname, kmylabel, kmltitle, kmlpos, kmmsize, kmlsize, kmpsize, kmaxsize, kmpxc, kmpyc, survpvalue, kmpvalue){
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

  if(survpvalue > kmpvalue){
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
