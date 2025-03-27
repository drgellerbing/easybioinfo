#' @import magrittr

multicox = function(kmdf){
  cat("Executing Multivariate Cox Analysis ... \n")
  gts <- getgenes()
  gts <- gsub("-", "", gts)
  kmdf <- as.data.frame(kmdf)
  gts <- gts[gts %in% colnames(kmdf)]

  type <- readline(prompt = "Are you looking into the overall survival or progression free survival? (os/pfs): ")
  cutoff <- readline(prompt = "Please enter your cutoff method (median/quartile/custom/zscore/optimal/continuous): ")
  while(tolower(cutoff) != "median" && tolower(cutoff) != "quartile" && tolower(cutoff) != "custom" && tolower(cutoff) != "zscore" && tolower(cutoff) != "optimal" && tolower(cutoff) != "continuous"){
    message("Answer not given, please try again :)")
    cutoff <- readline(prompt = "Please enter your cutoff method (median/quartile/custom/zscore/optimal/continuous): ")
  }

  #Converting column names in dataframe into lowercase
  if(tolower(cutoff) == "median"){
    for(col_name in colnames(kmdf)) {
      if(tolower(col_name) %in% tolower(gts)) {
        colnames(kmdf)[colnames(kmdf) == col_name] <- tolower(col_name)
      }
    }

    #Segregating Genes based on chosen method
    rm_genes <- vector("character", length = 0)
    for (i in unique(1:length(gts))){
      z = gts[i]
      l = tolower(z)
      name = paste(z)
      if (exists(l, kmdf)){
        cat("Calculating multivariate survival for", z, "\n")
        kmdf %>%dplyr::select(dplyr::all_of(l)) %>%unlist() -> c
      } else {
        warning(paste("Column '", z, "' not found in kmdf. Skipping...", sep = ""))
      }
      mdn = stats::median(c)
      kmdf %<>% dplyr::mutate("{name}" := as.factor(
        ifelse(
          c >= mdn,
          'High',
          'Low'
        )))
      fctcheck <- forcats::fct_count(kmdf[[name]])
      if(length(rownames(fctcheck)) <2){
        rm_genes <- c(rm_genes, name)
        kmdf %<>% dplyr::select(-dplyr::all_of(name))
      }
      if(length(rownames(fctcheck))== 2){kmdf[[name]] <- stats::relevel(factor(kmdf[[name]]), ref = 'Low')}
    }

    #Reporting Errors
    if(length(rm_genes) > 0){
      for(m in 1:length(rm_genes)){cat("Gene", rm_genes[m], "was removed from the analysis. Unable to segregate expression into two groups\n")}
    }

    #Calculating Cox
    if(type == "os"){survdf <- survival::Surv(time = kmdf$os, event = kmdf$statusos)}
    if(type == "pfs"){survdf <- survival::Surv(time = kmdf$pfs, event = kmdf$statuspfs)}
    f1 <- stats::as.formula(paste('survdf ~', paste(gts, collapse = "+")))
    cox = survival::coxph(f1, data = kmdf)
    zph = survival::cox.zph(cox)
  }

  ## Quartile | Custom Cutoff
  if(tolower(cutoff) == "quartile" || tolower(cutoff) == "custom"){
    if(tolower(cutoff) == "custom"){
      qval_low <- as.numeric(readline(prompt = "Please enter the low % cutoff between 0 and 1: "))
      qval_high <- as.numeric(readline(prompt = "Please enter the high % cutoff between 0 and 1: "))
    }

    #Converting column names in dataframe into lowercase
    for(col_name in colnames(kmdf)) {
      if(tolower(col_name) %in% tolower(gts)) {
        colnames(kmdf)[colnames(kmdf) == col_name] <- tolower(col_name)
      }
    }

    #Segregating Genes based on chosen method
    rm_genes <- vector("character", length = 0)
    for (i in unique(1:length(gts))){
      z = gts[i]
      l = tolower(z)
      name = paste(z)
      if (exists(l, kmdf)){
        cat("Calculating multivariate survival for", z, "\n")
        kmdf %>%dplyr::select(dplyr::all_of(l)) %>%unlist() -> c
      } else {
        warning(paste("Column '", z, "' not found in kmdf. Skipping...", sep = ""))
      }
      if(tolower(cutoff) == "quartile"){quartile <- stats::quantile(c, probs = c(0.25, 0.75))}
      if(tolower(cutoff) == "custom"){quartile <- stats::quantile(c, probs = c(qval_low, qval_high))}

      kmdf %<>% dplyr::mutate("{name}" := as.factor(
        ifelse(
          c >= quartile[[2]],
          'High',
          ifelse(
            c <= quartile[[1]],
            "Low",
            "Medium"
          )
        )))
      fctcheck <- forcats::fct_count(kmdf[[name]])
      if(length(rownames(fctcheck)) <3){
        rm_genes <- c(rm_genes, name)
        kmdf %<>% dplyr::select(-dplyr::all_of(name))
      }
      if(length(rownames(fctcheck))== 3){kmdf[[name]] <- factor(kmdf[[name]], levels = c("Medium", "Low", "High"))}
    }

    #Reporting Errors
    if(length(rm_genes) > 0){
      for(m in 1:length(rm_genes)){cat("Gene", rm_genes[m], "was removed from the analysis. Unable to segregate expression into two groups\n")}
    }

    #Calculating Cox
    if(type == "os"){survdf <- survival::Surv(time = kmdf$os, event = kmdf$statusos)}
    if(type == "pfs"){survdf <- survival::Surv(time = kmdf$pfs, event = kmdf$statuspfs)}
    f1 <- stats::as.formula(paste('survdf ~', paste(gts, collapse = "+")))
    cox = survival::coxph(f1, data = kmdf)
    zph = survival::cox.zph(cox)
  }

  if(tolower(cutoff) == "zscore"){
    zcut <- as.numeric(readline(prompt = "Please enter your desired z-score cutoff: "))

    #Converting column names in dataframes to lowercase
    for(col_name in colnames(kmdf)) {
      if(tolower(col_name) %in% tolower(gts)) {
        colnames(kmdf)[colnames(kmdf) == col_name] <- tolower(col_name)
      }
    }

    #Segregating genes based on chosen method
    rm_genes <- vector("character", length = 0)
    for (i in unique(1:length(gts))){
      z = gts[i]
      l = tolower(z)
      name = paste(z)
      if (exists(l, kmdf)){
        cat("Calculating multivariate survival for", z, "\n")
        kmdf %>%dplyr::select(dplyr::all_of(l)) %>%unlist() -> c
      } else {
        warning(paste("Column '", z, "' not found in kmdf. Skipping...", sep = ""))
      }
      zsc = scale(c)

      kmdf %<>% dplyr::mutate("{name}" := as.factor(
        ifelse(
          zsc >= zcut,
          'High',
          ifelse(
            zsc <= -abs(zcut),
            "Low",
            "Medium"
          )
        )))
      fctcheck <- forcats::fct_count(kmdf[[name]])
      if(length(rownames(fctcheck)) <3){
        rm_genes <- c(rm_genes, name)
        kmdf %<>% dplyr::select(-dplyr::all_of(name))
      }
      if(length(rownames(fctcheck))== 3){kmdf[[name]] <- factor(kmdf[[name]], levels = c("Medium", "Low", "High"))}
    }

    #Reporting Errors
    if(length(rm_genes) > 0){
      for(m in 1:length(rm_genes)){cat("Gene", rm_genes[m], "was removed from the analysis. Unable to segregate expression into two groups\n")}
    }

    #Cox calculation
    if(type == "os"){survdf <- survival::Surv(time = kmdf$os, event = kmdf$statusos)}
    if(type == "pfs"){survdf <- survival::Surv(time = kmdf$pfs, event = kmdf$statuspfs)}
    f1 <- stats::as.formula(paste('survdf ~', paste(gts, collapse = "+")))
    cox = survival::coxph(f1, data = kmdf)
    zph = survival::cox.zph(cox)
  }

  ## Optimal Cutoff
  if(tolower(cutoff) == "optimal"){
    #Converting column names in dataframe to lowercase
    rm_genes <- vector("character", length = 0)
    for(col_name in colnames(kmdf)) {
      if(tolower(col_name) %in% tolower(gts)) {
        colnames(kmdf)[colnames(kmdf) == col_name] <- tolower(col_name)
      }
    }

    #Segregating genes based on chosen method
    rm_genes <- vector("character", length = 0)
    if(type == "os"){survdf <- survival::Surv(time = kmdf$os, event = kmdf$statusos)}
    if(type == "pfs"){survdf <- survival::Surv(time = kmdf$pfs, event = kmdf$statuspfs)}
    for (i in unique(1:length(gts))){
      z = gts[i]
      l = tolower(z)
      name = paste(z)
      if (exists(l, kmdf)){
        cat("Calculating multivariate survival for", z, "\n")
        kmdf %>%dplyr::select(dplyr::all_of(l)) %>%unlist() -> c
        fc <- stats::as.formula(paste('survdf ~', paste(l)))
        tryCatch({
          mscutoff <- maxstat::maxstat.test(fc, data = kmdf, smethod = "LogRank")
          o_cutoff <- mscutoff$estimate %>%as.double()

          kmdf %<>% dplyr::mutate("{name}" := as.factor(
            ifelse(
              c >= o_cutoff,
              'High',
              'Low'
            )))
        }, error = function(e) {
          # Handle the error
          # Print an error message or take any other necessary actions
          print(paste("Error:", z, "contains non-zero cell counts within the specified proportion range, removed from the analysis"))
          # You can also choose to continue the loop without any further action
        })
      } else {
        warning(paste("Column '", z, "' not found in kmdf. Skipping...", sep = ""))
      }

      # Checking Segregated genes for two factors before proceeding
      fctcheck <- forcats::fct_count(kmdf[[name]])
      if(length(rownames(fctcheck)) <2){
        rm_genes <- c(rm_genes, name)
        kmdf %<>% dplyr::select(-dplyr::all_of(name))
      }
      if(length(rownames(fctcheck))== 2){kmdf[[name]] <- stats::relevel(factor(kmdf[[name]]), ref = 'Low')}
    }

    #Reporting Errors
    if(length(rm_genes) > 0){
      for(m in 1:length(rm_genes)){cat("Gene", rm_genes[m], "was removed from the analysis. Unable to segregate expression into two groups\n")}
    }

    #Cox Calculation
    f1 <- stats::as.formula(paste('survdf ~', paste(gts, collapse = "+")))
    cox = survival::coxph(f1, data = kmdf)
    zph = survival::cox.zph(cox)
  }

  if(tolower(cutoff) == "continuous"){
    if(type == "os"){survdf <- survival::Surv(time = kmdf$os, event = kmdf$statusos)}
    if(type == "pfs"){survdf <- survival::Surv(time = kmdf$pfs, event = kmdf$statuspfs)}
    f1 <- stats::as.formula(paste('survdf ~', paste(gts, collapse = "+")))
    cox = survival::coxph(f1, data = kmdf)
  }

  unicolumn = c('Gene', 'Hazard Ratio', 'CI', 'P-value', 'PH Check')
  coxtable = data.frame(matrix(nrow=0, ncol = length(unicolumn)))
  suppressWarnings({
    p <- as.numeric(readline(prompt = "Please state your desired p-value cutoff: "))
    while(is.na(p) == TRUE || p <= 0 || p > 1){
      message("The p-value provided should be within 0 and 1, please try again :>")
      p <- as.numeric(readline(prompt = "Please state your desired p-value cutoff: "))
    }
  })
  suppressWarnings({
    civalue <- as.numeric(readline(prompt = "Please enter your desired confidence interval between 0 and 1: "))
    while(is.na(civalue) == TRUE || civalue <= 0 || civalue >= 1){
      message("The CI value provided should be within 0 and 1, please try again :)")
      civalue <- as.numeric(readline(prompt = "Please enter your desired confidence interval between 0 and 1: "))
    }
  })

  ucoef = summary(cox)$coef[,2]
  upvalue = summary(cox)$coef[,5]
  ci = stats::confint(cox, level = civalue)
  ci <- exp(ci)
  ci <- round(ci, digits = 2)
  uci <- paste(ci[,1], '~', ci[,2], sep = " ")
  phcheck <- zph$table[1:(length(rownames(zph$table)) - 1), 3]

  if(tolower(cutoff) == "median" || tolower(cutoff) == "optimal" || tolower(cutoff) == "continuous"){df = data.frame(gts, ucoef, uci, upvalue, phcheck)}
  if(tolower(cutoff) == "quartile" || tolower(cutoff) == "custom" || tolower(cutoff) == "zscore"){
    dfname = rownames(summary(cox)$coef)
    dfname <- gsub("High", " High", dfname)
    dfname <- gsub("Low", " Low", dfname)
    df = data.frame(dfname, ucoef, uci, upvalue, phcheck)
  }
  coxtable <- rbind(coxtable, df)
  coxtable[length(rownames(coxtable)) + 1, 5] <- paste("Global = ", zph$table[length(rownames(zph$table)), 3])
  colnames(coxtable) = unicolumn
  rownames(coxtable) = NULL
  coxtable %>%dplyr::filter(coxtable$`P-value` <= p) %>%as.data.frame() -> sigcox

  graph <- survminer::ggforest(cox, data = kmdf, refLabel = "Reference") %>%print()

  filen <- readline(prompt = "Please provide a name for your results: ")
  list1 <- (list(cox, kmdf, coxtable, sigcox))
  names(list1) <- c("cox", "kmdf", "full_cox", "sig_cox")
  res1 <- as.data.frame(list1[[3]])
  name1 <- paste("./", filen, "_full_multicox.xlsx", sep = "")

  writexl::write_xlsx(res1, name1)

  name3 <- paste("./", filen, "_forestplot.png", sep = "")
  ggplot2::ggsave(file = name3, width = 15, height = 15, dpi = 400)

  cat("Analysis Complete :> Hope you get something :> \n")
  return(list1)
}

plotforest = function(cox){
  coxdf <- cox[[1]]
  kmdf <- cox[[2]]

  title <- readline(prompt = "Please enter the title of the plot: ")
  fsize <- as.numeric(readline(prompt = "Please enter your desired font size (Default = 0.7): "))
  ref <- readline(prompt = "Please enter the label for reference levels of factor variables: ")
  ndigit <- as.numeric(readline(prompt = "Please enter your desired number of digits for estimates and p-values in the plot: "))

  survminer::ggforest(coxdf, kmdf, main = title, fontsize = fsize, refLabel = ref, noDigits = ndigit)

  filen <- readline(prompt = "Please provide a name for your forest plot: ")
  name <- paste("./", filen, "_forestplot.png", sep = "")

  nw <- as.numeric(readline(prompt = "Please enter the value of the width of the final figure (Default provided previously = 15): "))
  nh <- as.numeric(readline(prompt = "Please enter the value of the height of the final figure (Default provided previously = 15): "))

  ggplot2::ggsave(file = name, width = nw, height = nh, dpi = 500)
}
