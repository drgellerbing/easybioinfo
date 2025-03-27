#' @import magrittr

##uniKM Function
uniKM = function(kmdf){
  cat("Executing Kaplan Meier Analysis .... \n")
  kmdf <- data.table::as.data.table(kmdf)

  ## Backbone of Final Dataframe
  unicolumn = c('Gene', 'Coefficient', 'HR', '95%CI', 'Pvalue', 'PHCheck')
  kmtable = data.frame(matrix(nrow=0, ncol = length(unicolumn)))
  colnames(kmtable) = unicolumn

  ## Grabbing Query Genes
  gts <- getgenes()
  gts <- gsub("-", "", gts)
  gts <- gts[gts %in% colnames(kmdf)]

  ## Cutoffs and Plot Type Query
  message("Do you want to set the cutoff based off the median, quartile, custom %, z-score, or individual optimal cutoffs?")
  cutoff <- readline(prompt = "Please enter your cutoff method (median/quartile/custom/zscore/optimal): ")
  while(tolower(cutoff) != "median" && tolower(cutoff) != "quartile" && tolower(cutoff) != "custom" && tolower(cutoff) != "zscore" && tolower(cutoff) != "optimal"){
    message("Answer not given, please try again :)")
    cutoff <- readline(prompt = "Please enter your cutoff method (median/quartile/custom/zscore/optimal): ")
  }

  #P-value query
  suppressWarnings({
    p <- as.numeric(readline(prompt = "Please state your desired p-value cutoff: "))
    while(is.na(p) == TRUE || p <= 0 || p > 1){
      message("The p-value provided should be within 0 and 1, please try again :>")
      p <- as.numeric(readline(prompt = "Please state your desired p-value cutoff: "))
    }
  })

  #Plot Types
  plot_type <- readline(prompt = "Are you looking at overall survival or progression free survival? (osc/pfs): ")
  while(tolower(plot_type) != "osc" && tolower(plot_type) != "pfs"){
    message("Answer not given, please try again :>")
    plot_type <- readline(prompt = "Are you looking at overall survival curves or PFS curves? (osc/pfs): ")
  }

  plot <- readline(prompt = "Do you want to plot the survival curves? (yes/no): ")
  while(tolower(plot) != "yes" && tolower(plot) != "no"){
    message("Answer not given, please try again :>")
    plot <- readline(prompt = "Do you want to plot the survival curves? (yes/no): ")
  }
  if(tolower(plot) == "yes" || tolower(plot) == "y"){
    #Defining Plot Parameters
    lpos <- "0.65, 0.9"
    lpos <- unlist(strsplit(lpos, ", "))
    lpos <- as.numeric(lpos)
    msize <- as.numeric(18)
    lsize <- as.numeric(16)
    psize <- as.numeric(16)
    axsize <- as.numeric(16)
    pxc <- as.numeric(0.45)
    pyc <- as.numeric(38)

    if(length(gts) > 10){
      plot_sig <- readline(prompt = "Do you want to plot only the significant results? (yes/no): ")
      while(tolower(plot_sig) != "yes" && tolower(plot_sig) != "no" && tolower(plot_sig) != "y" && tolower(plot_sig) != "n"){
        message("Answer not given, please try again :>")
        plot_sig <- readline(prompt = "Do you want to plot only the significant results? (yes/no): ")
      }
    } else{
      plot_sig <- "no"
    }

    #Defining Label Title & Survival Calculations
    if(plot_type == "osc"){
      kmlabel <- "Survival Probability"
      ltitle <- "Median Survival"
    }

    if(plot_type == "pfs"){
      kmlabel <- "Progression-Free Survival Probability"
      ltitle <- "Median PFS"
    }

    #Folder Creation
    if(length(gts) > 10){
      dir <- readline(prompt = "Do you want to create a new folder for your plots? (yes/no): ")
      while(tolower(dir) != "yes" && tolower(dir) != "y" && tolower(dir) != "no" && tolower(dir) != "n"){
        message("Answer not given, please try again :>")
        dir <- readline(prompt = "Do you want to create a new folder for your plots & results? (yes/no): ")
      }

      if(tolower(dir) == "yes" || tolower(dir) == "y"){
        message("Please make sure the name of your created folder does not exists in your current directory")
        success <- FALSE
        while(!success){
          dir_n <- readline(prompt = "Please enter the name for your newly created folder: ")
          if(dir.exists(dir_n)){
            message("Folder already exists. Please enter a different name :)")
          } else{
            dir.create(dir_n)
            success <- TRUE
          }
        }
      }
    } else{
      dir <- "no"
    }
  }

  # Schoenfield Query
  plotsr <- readline(prompt = "Do you want the Schoenfeld residuals plot? (yes/no): ")
  while(tolower(plotsr) != "yes" && tolower(plotsr) != "no"){
    message("Answer not given, please try again :>")
    plotsr <- readline(prompt = "Do you want the Schoenfeld residuals plot? (yes/no): ")
  }
  if(tolower(plotsr) == "yes" || tolower(plotsr) == "y"){
    message("Schoenfield residuals plot will be placed under a folder named 'PH Check'")
    if(!file.exists("PH Check")){
      dir.create("PH Check")
    }
    if(length(gts) > 10){
      plotsr_sig <- readline(prompt = "Do you want to plot only the significant results? (yes/no): ")
      while(tolower(plotsr_sig) != "yes" && tolower(plotsr_sig) != "no" && tolower(plotsr_sig) != "y" && tolower(plotsr_sig) != "n"){
        message("Answer not given, please try again :>")
        plotsr_sig <- readline(prompt = "Do you want to plot only the significant results? (yes/no): ")
      }
    } else{
      plotsr_sig <- "no"
    }
  }

  ## Median Cutoff
  if(tolower(cutoff) == "median"){
    if(plot_type == "osc"){survdf <- survival::Surv(time = kmdf$os, event = kmdf$statusos)}
    if(plot_type == "pfs"){survdf <- survival::Surv(time = kmdf$pfs, event = kmdf$statuspfs)}

    # Segregating Expression Based on Median
    for (y in unique(1:length(gts)))
    {
      z = gts[y]
      name = paste(z, 'Expression', sep = "")
      if (exists(z, kmdf)) {kmdf %>% dplyr::select(dplyr::all_of(z)) %>% unlist() -> c}
      mdn = stats::median(c)
      kmdf %<>% dplyr::mutate("{name}" := as.factor(
        ifelse(
          c >= mdn,
          'High',
          'Low'
        )))
      }

    #Removing Genes That Cannot be Segregated
    rm_genes <- vector("character", length = 0)
    gtslist1 = grep("Expression", colnames(kmdf), value = TRUE)
    for (o in gtslist1){
      fctcheck <- forcats::fct_count(kmdf[[o]])
      if(length(rownames(fctcheck)) <2){
        rm_genes <- c(rm_genes, o)
        kmdf %<>% dplyr::select(-dplyr::all_of(o))
      }
      if(length(rownames(fctcheck))== 2){kmdf[[o]] <- stats::relevel(factor(kmdf[[o]]), ref = 'Low')}
    }
    rm_genes <- gsub("Expression", "", rm_genes)
    if(length(rm_genes) > 0){
      for(m in 1:length(rm_genes)){cat("Gene", rm_genes[m], "was removed from the analysis. Unable to segregate expression into two groups\n")}
    }
    Sys.sleep(0.5)

    gtslist = grep("Expression", colnames(kmdf), value = TRUE)
    gts <- gsub("Expression", "", gtslist)

    for (i in unique(1:length(gts))) {
      x = gts[i]
      y = gtslist[i]

      if (exists(y, kmdf)) {
        cat("Now processing survival analysis for", x, "\n" )
        f2 <- stats::as.formula(paste('survdf ~', paste(y)))
        if (!is.null(f2)){
          unidf2 <- onegenecox(f2, x, kmdf)
          zph = survival::cox.zph(survival::coxph(f2, data = kmdf))
          unidf2 <- data.frame(unidf2, zph$table[1,3])
          colnames(unidf2) <- unicolumn
          kmtable <- rbind(kmtable, unidf2)

          if(tolower(plot) == "yes" || tolower(plot) == "y"){
            pv <- survminer::surv_pvalue(survminer::surv_fit(f2, data = kmdf))[[2]]
            if(tolower(plot_sig) == "yes"){
              if(pv <= p){
                graph <- kmplot(f2, kmdf, x, kmlabel, ltitle, lpos, msize, lsize, psize, axsize, pxc, pyc, pv, p)

                if(tolower(dir) == "yes" || tolower(dir) == "y"){ggplot2::ggsave(paste(dir_n, "/", x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
                if(tolower(dir) == "no" || tolower(dir) =="n"){ggplot2::ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
              }
            }

            if(tolower(plot_sig) == "no" || tolower(plot_sig) == "n"){
              graph <- kmplot(f2, kmdf, x, kmlabel, ltitle, lpos, msize, lsize, psize, axsize, pxc, pyc, pv, p)

              if(tolower(dir) == "yes" || tolower(dir) == "y"){ggplot2::ggsave(paste(dir_n, "/", x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
              if(tolower(dir) == "no" || tolower(dir) =="n"){ggplot2::ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
            }
          }
          if(tolower(plotsr) == "yes" || tolower(plotsr) == "y"){
            if(tolower(plotsr_sig) == "yes" || tolower(plotsr_sig) == "y"){
              pv_sr <- survminer::surv_pvalue(survminer::surv_fit(f2, data = kmdf))[[2]]
              if(pv_sr <= p){
                grDevices::png(paste("./PH Check/", x, ".png", sep = ""), res = 300, units = "cm", width = 15, height = 15)
                plot(zph)
                grDevices::dev.off()
              }
            }
            if(tolower(plotsr_sig) == "no" || tolower(plotsr_sig) == "n"){
              grDevices::png(paste("./PH Check/", x, ".png", sep = ""), res = 300, units = "cm", width = 15, height = 15)
              plot(zph)
              grDevices::dev.off()
            }
          }
        }
        else {
          warning(paste("Gene", x, "not found in dataframe. Skipping..."))
        }
      }
    }
  }

  ## Quartile Cutoff
  if(tolower(cutoff) == "quartile" || tolower(cutoff) == "custom"){
    if(tolower(cutoff) == "custom"){
      valid_input <- FALSE

      while(!valid_input){
        qval_low <- suppressWarnings(as.numeric(readline(prompt = "Please enter the low % cutoff between 0 and 1: ")))
        qval_high <- suppressWarnings(as.numeric(readline(prompt = "Please enter the high % cutoff between 0 and 1: ")))

        if(!is.na(qval_low) && !is.na(qval_high) &&
           qval_low >= 0 && qval_low <= 1 &&
           qval_high >= 0 && qval_high <= 1 &&
           qval_low < qval_high){
          valid_input <- TRUE
        } else {
          cat("Invalid input. Please ensure that:\n")
          cat("1. Both values are between 0 and 1.\n")
          cat("2. The low cutoff is less than the high cutoff.\n")
        }
      }
    }

    # Segregating based on quartile
    for (y in unique(1:length(gts)))
    {
      z = gts[y]
      name = paste(z, 'Expression', sep = "")
      if (exists(z, kmdf)) {
        kmdf %>% dplyr::select(dplyr::all_of(z)) %>% unlist() -> c
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
    }

    # Removing Genes that failed to segregate in three different groups
    rm_genes <- vector("character", length = 0)
    gtslist1 = grep("Expression", colnames(kmdf), value = TRUE)
    for (o in gtslist1){
      fctcheck <- forcats::fct_count(kmdf[[o]])
      if(length(rownames(fctcheck)) <3){
        rm_genes <- c(rm_genes, o)
        kmdf %<>% dplyr::select(-dplyr::all_of(o))
      }
      if(length(rownames(fctcheck))== 3){kmdf[[o]] <- stats::relevel(factor(kmdf[[o]]), ref = 'Low')}
    }
    rm_genes <- gsub("Expression", "", rm_genes)
    if(length(rm_genes) > 0){
      for(m in 1:length(rm_genes)){cat("Gene", rm_genes[m], "was removed from the analysis. Unable to segregate expression into three groups\n")}
    }
    Sys.sleep(0.5)

    #Final Analysis
    gtslist = grep("Expression", colnames(kmdf), value = TRUE)
    gts <- gsub("Expression", "", gtslist)

    for (i in unique(1:length(gts))) {
      x = gts[i]
      y = gtslist[i]

      if (exists(y, kmdf)) {
        cat("Now processing survival analysis for", x, "\n" )
        f2 <- stats::as.formula(paste('survdf ~', paste(y)))
        if (!is.null(f2)){
          if(plot_type == "osc"){
            kmdf %>% dplyr::select(id, os, vitalos, statusos, dplyr::all_of(y)) %>% data.table::as.data.table() -> kmdf2
            kmdf2 %<>% dplyr::filter(kmdf2[,5] != "Medium")
            survdf <- survival::Surv(time = kmdf2$os, event = kmdf2$statusos)
          }
          if(plot_type == "pfs"){
            kmdf %>% dplyr::select(id, pfs, vitalpfs, statuspfs, dplyr::all_of(y)) %>% data.table::as.data.table() -> kmdf2
            kmdf2 %<>% dplyr::filter(kmdf2[,5] != "Medium")
            survdf <- survival::Surv(time = kmdf2$pfs, event = kmdf2$statuspfs)
          }
          unidf2 <- onegenecox(f2, x, kmdf2)
          zph = survival::cox.zph(survival::coxph(f2, data = kmdf2))
          unidf2 <- data.frame(unidf2, zph$table[1,3])
          colnames(unidf2) <- unicolumn
          kmtable <- rbind(kmtable, unidf2)

          if(tolower(plot) == "yes" || tolower(plot) == "y"){
            pv <- survminer::surv_pvalue(survminer::surv_fit(f2, data = kmdf2))[[2]]
            if(tolower(plot_sig) == "yes"){
              if(pv <= p){
                graph <- kmplot(f2, kmdf2, x, kmlabel, ltitle, lpos, msize, lsize, psize, axsize, pxc, pyc, pv, p)

                if(tolower(dir) == "yes" || tolower(dir) == "y"){ggplot2::ggsave(paste(dir_n, "/", x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
                if(tolower(dir) == "no" || tolower(dir) =="n"){ggplot2::ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
              }
            }

            if(tolower(plot_sig) == "no" || tolower(plot_sig) == "n"){
              graph <- kmplot(f2, kmdf2, x, kmlabel, ltitle, lpos, msize, lsize, psize, axsize, pxc, pyc, pv, p)

              if(tolower(dir) == "yes" || tolower(dir) == "y"){ggplot2::ggsave(paste(dir_n, "/", x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
              if(tolower(dir) == "no" || tolower(dir) =="n"){ggplot2::ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
            }
          }
          if(tolower(plotsr) == "yes" || tolower(plotsr) == "y"){
            if(tolower(plotsr_sig) == "yes" || tolower(plotsr_sig) == "y"){
              pv_sr <- survminer::surv_pvalue(survminer::surv_fit(f2, data = kmdf2))[[2]]
              if(pv_sr <= p){
                grDevices::png(paste("./PH Check/", x, ".png", sep = ""), res = 300, units = "cm", width = 15, height = 15)
                plot(zph)
                grDevices::dev.off()
              }
            }
            if(tolower(plotsr_sig) == "no" || tolower(plotsr_sig) == "n"){
              grDevices::png(paste("./PH Check/", x, ".png", sep = ""), res = 300, units = "cm", width = 15, height = 15)
              plot(zph)
              grDevices::dev.off()
            }
          }
        }
        else {
          warning(paste("Gene", x, "not found in dataframe. Skipping..."))
        }
      }
    }
  }

  ## Z-Score Cutoff
  if(tolower(cutoff) == "zscore"){
    zcut <- as.numeric(readline(prompt = "Please enter your desired z-score cutoff: "))

    for (y in unique(1:length(gts))){
      z = gts[y]
      name = paste(z, 'Expression', sep = "")
      if (exists(z, kmdf)) {
        kmdf %>% dplyr::select(dplyr::all_of(z)) %>% unlist() -> c
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
    }

    #Removing Genes That Cannot be Segregated
    rm_genes <- vector("character", length = 0)
    gtslist1 = grep("Expression", colnames(kmdf), value = TRUE)
    for (o in gtslist1){
      fctcheck <- forcats::fct_count(kmdf[[o]])
      if(length(rownames(fctcheck)) <3){
        rm_genes <- c(rm_genes, o)
        kmdf %<>% dplyr::select(-dplyr::all_of(o))
      }
      if(length(rownames(fctcheck))== 3){kmdf[[o]] <- stats::relevel(factor(kmdf[[o]]), ref = 'Low')}
    }
    rm_genes <- gsub("Expression", "", rm_genes)
    if(length(rm_genes) > 0){
      for(m in 1:length(rm_genes)){cat("Gene", rm_genes[m], "was removed from the analysis. Unable to segregate expression into three groups\n")}
    }
    Sys.sleep(0.5)

    #Final Analysis
    gtslist = grep("Expression", colnames(kmdf), value = TRUE)
    gts <- gsub("Expression", "", gtslist)

    for (i in unique(1:length(gts))) {
      x = gts[i]
      y = gtslist[i]

      if (exists(y, kmdf)) {
        cat("Now processing survival analysis for", x, "\n" )
        f2 <- stats::as.formula(paste('survdf ~', paste(y)))
        if (!is.null(f2)){
          if(plot_type == "osc"){
            kmdf %>% dplyr::select(id, os, vitalos, statusos, dplyr::all_of(y)) %>% data.table::as.data.table() -> kmdf2
            kmdf2 %<>% dplyr::filter(kmdf2[,5] != "Medium")
            survdf <- survival::Surv(time = kmdf2$os, event = kmdf2$statusos)
          }
          if(plot_type == "pfs"){
            kmdf %>% dplyr::select(id, pfs, vitalpfs, statuspfs, dplyr::all_of(y)) %>% data.table::as.data.table() -> kmdf2
            kmdf2 %<>% dplyr::filter(kmdf2[,5] != "Medium")
            survdf <- survival::Surv(time = kmdf2$pfs, event = kmdf2$statuspfs)
          }
          unidf2 <- onegenecox(f2, x, kmdf2)
          zph = survival::cox.zph(survival::coxph(f2, data = kmdf2))
          unidf2 <- data.frame(unidf2, zph$table[1,3])
          colnames(unidf2) <- unicolumn
          kmtable <- rbind(kmtable, unidf2)

          if(tolower(plot) == "yes" || tolower(plot) == "y"){
            pv <- survminer::surv_pvalue(survminer::surv_fit(f2, data = kmdf2))[[2]]
            if(tolower(plot_sig) == "yes"){
              if(pv <= p){
                graph <- kmplot(f2, kmdf2, x, kmlabel, ltitle, lpos, msize, lsize, psize, axsize, pxc, pyc, pv, p)

                if(tolower(dir) == "yes" || tolower(dir) == "y"){ggplot2::ggsave(paste(dir_n, "/", x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
                if(tolower(dir) == "no" || tolower(dir) =="n"){ggplot2::ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
              }
            }

            if(tolower(plot_sig) == "no" || tolower(plot_sig) == "n"){
              graph <- kmplot(f2, kmdf2, x, kmlabel, ltitle, lpos, msize, lsize, psize, axsize, pxc, pyc, pv, p)

              if(tolower(dir) == "yes" || tolower(dir) == "y"){ggplot2::ggsave(paste(dir_n, "/", x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
              if(tolower(dir) == "no" || tolower(dir) =="n"){ggplot2::ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
            }
          }
          if(tolower(plotsr) == "yes" || tolower(plotsr) == "y"){
            if(tolower(plotsr_sig) == "yes" || tolower(plotsr_sig) == "y"){
              pv_sr <- survminer::surv_pvalue(survminer::surv_fit(f2, data = kmdf2))[[2]]
              if(pv_sr <= p){
                grDevices::png(paste("./PH Check/", x, ".png", sep = ""), res = 300, units = "cm", width = 15, height = 15)
                plot(zph)
                grDevices::dev.off()
              }
            }
            if(tolower(plotsr_sig) == "no" || tolower(plotsr_sig) == "n"){
              grDevices::png(paste("./PH Check/", x, ".png", sep = ""), res = 300, units = "cm", width = 15, height = 15)
              plot(zph)
              grDevices::dev.off()
            }
          }
        }
        else {
          warning(paste("Gene", x, "not found in dataframe. Skipping..."))
        }
      }
    }
  }

  ## Segregating Genes based on optimal AOC
  if(tolower(cutoff) == "optimal"){
    if(plot_type == "osc"){survdf <- survival::Surv(time = kmdf$os, event = kmdf$statusos)}
    if(plot_type == "pfs"){survdf <- survival::Surv(time = kmdf$pfs, event = kmdf$statuspfs)}

    ## Backbone of final dataframes
    cutoffcol = c("Gene", "Cutoff")
    cuttable = data.frame(matrix(nrow=0, ncol = length(cutoffcol)))
    colnames(cuttable) = cutoffcol

    #Grouping Gene Expression into Groups
    for (y in unique(1:length(gts))){
      z = gts[y]
      name = paste(z, 'Expression', sep = "")
      if (exists(z, kmdf)) {
        kmdf %>% dplyr::select(dplyr::all_of(z)) %>% unlist() -> c
        fc <- stats::as.formula(paste('survdf ~', paste(z)))
        tryCatch({
          mscutoff <- maxstat::maxstat.test(fc, data = kmdf, smethod = "LogRank")
          o_cutoff <- mscutoff$estimate %>% as.double()
          cutdf <- data.frame(z, o_cutoff)
          cuttable <- rbind(cuttable, cutdf)

          kmdf %<>% dplyr::mutate("{name}" := as.factor(
            ifelse(
              c >= o_cutoff,
              'High',
              'Low'
            )))
        }, error = function(e) {
          cat(z, "contains non-zero cell counts within the specified proportion range, removed from the analysis \n")
        })
      }
    }

    #Removing Genes That Cannot be Segregated
    rm_genes <- vector("character", length = 0)
    gtslist1 = grep("Expression", colnames(kmdf), value = TRUE)
    for (o in gtslist1){
      fctcheck <- forcats::fct_count(kmdf[[o]])
      if(length(rownames(fctcheck)) <2){
        rm_genes <- c(rm_genes, o)
        kmdf %<>% dplyr::select(-dplyr::all_of(o))
      }
      if(length(rownames(fctcheck))== 2){kmdf[[o]] <- stats::relevel(factor(kmdf[[o]]), ref = 'Low')}
    }
    rm_genes <- gsub("Expression", "", rm_genes)
    if(length(rm_genes) > 0){
      for(m in 1:length(rm_genes)){cat("Gene", rm_genes[m], "was removed from the analysis. Unable to segregate expression into two groups\n")}
    }
    Sys.sleep(0.5)

    #Final Analysis
    gtslist = grep("Expression", colnames(kmdf), value = TRUE)
    gts <- gsub("Expression", "", gtslist)

    for (i in unique(1:length(gts))) {
      x = gts[i]
      y = gtslist[i]

      if (exists(y, kmdf)) {
        cat("Now processing survival analysis for", x, "\n" )
        f2 <- stats::as.formula(paste('survdf ~', paste(y)))
        if (!is.null(f2)){
          unidf2 <- onegenecox(f2, x, kmdf)
          zph = survival::cox.zph(survival::coxph(f2, data = kmdf))
          unidf2 <- data.frame(unidf2, zph$table[1,3])
          colnames(unidf2) <- unicolumn
          kmtable <- rbind(kmtable, unidf2)

          if(tolower(plot) == "yes" || tolower(plot) == "y"){
            pv <- survminer::surv_pvalue(survminer::surv_fit(f2, data = kmdf))[[2]]
            if(tolower(plot_sig) == "yes"){
              if(pv <= p){
                graph <- kmplot(f2, kmdf, x, kmlabel, ltitle, lpos, msize, lsize, psize, axsize, pxc, pyc, pv, p)

                if(tolower(dir) == "yes" || tolower(dir) == "y"){ggplot2::ggsave(paste(dir_n, "/", x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
                if(tolower(dir) == "no" || tolower(dir) =="n"){ggplot2::ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
              }
            }

            if(tolower(plot_sig) == "no" || tolower(plot_sig) == "n"){
              graph <- kmplot(f2, kmdf, x, kmlabel, ltitle, lpos, msize, lsize, psize, axsize, pxc, pyc, pv, p)

              if(tolower(dir) == "yes" || tolower(dir) == "y"){ggplot2::ggsave(paste(dir_n, "/", x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
              if(tolower(dir) == "no" || tolower(dir) =="n"){ggplot2::ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
            }
          }
          if(tolower(plotsr) == "yes" || tolower(plotsr) == "y"){
            if(tolower(plotsr_sig) == "yes" || tolower(plotsr_sig) == "y"){
              pv_sr <- survminer::surv_pvalue(survminer::surv_fit(f2, data = kmdf))[[2]]
              if(pv_sr <= p){
                grDevices::png(paste("./PH Check/", x, ".png", sep = ""), res = 300, units = "cm", width = 15, height = 15)
                plot(zph)
                grDevices::dev.off()
              }
            }
            if(tolower(plotsr_sig) == "no" || tolower(plotsr_sig) == "n"){
              grDevices::png(paste("./PH Check/", x, ".png", sep = ""), res = 300, units = "cm", width = 15, height = 15)
              plot(zph)
              grDevices::dev.off()
            }
          }
        }
        else {
          warning(paste("Gene", x, "not found in dataframe. Skipping..."))
        }
      }
    }
  }

  kmtable %>% dplyr::filter(kmtable$Pvalue < p) %>% as.data.frame() -> sigKM
  list1 <- (list(kmtable, sigKM, kmdf))
  names(list1) <- c("kmtable", "sigKM", "transformeddf")
  sig_return <- readline(prompt = "Do you only want to save all the results, or just the significant results? (all/sig/none): ")
  if(tolower(sig_return) == "all" || tolower(sig_return) == "sig"){filen <- readline(prompt = "Please provide a name for your results: ")}

  if(tolower(sig_return) == "all"){
    res1 <- as.data.frame(list1[[1]])
    res2 <- as.data.frame(list1[[2]])
    name1 <- paste("./", filen, "_full_uniKM.xlsx", sep = "")
    name2 <- paste("./", filen, "_sig_uniKM.xlsx", sep = "")
    if(tolower(cutoff) == "optimal"){namec <- paste("./", filen, "_gene_cutoff.xlsx", sep = "")}

    writexl::write_xlsx(res1, name1)
    writexl::write_xlsx(res2, name2)
    if(tolower(cutoff) == "optimal"){writexl::write_xlsx(cuttable, namec)}
  }
  if(tolower(sig_return) == "sig"){
    res <- as.data.frame(list1[[2]])
    name1 <- paste("./", filen, "_sig_uniKM.xlsx", sep = "")
    if(tolower(cutoff) == "optimal"){namec <- paste("./", filen, "_gene_cutoff.xlsx", sep = "")}
    writexl::write_xlsx(res, name1)
    if(tolower(cutoff) == "optimal"){writexl::write_xlsx(cuttable, namec)}
  }
  cat("All done :> \n")
  return(list1)
}

##Customisable Survival Curve Function
customisesurvivalplot = function(unidf){
  kmdf <- unidf[[3]]

  gts <- getgenes()
  gts <- gts[gts %in% colnames(kmdf)]

  plot_type <- readline(prompt = "Are you looking at overall survival or progression free survival? (osc/pfs): ")
  while(tolower(plot_type) != "osc" && tolower(plot_type) != "pfs"){
    message("Answer not given, please try again :>")
    plot_type <- readline(prompt = "Are you looking at overall survival curves or PFS curves? (osc/pfs): ")
  }

  if(plot_type == "osc"){
    kmlabel <- "Survival Probability"
    ltitle <- "Median Survival"
  }

  if(plot_type == "pfs"){
    kmlabel <- "Progression-Free Survival Probability"
    ltitle <- "Median PFS"
  }

  #P-value query
  suppressWarnings({
    p <- as.numeric(readline(prompt = "Please state your desired p-value cutoff: "))
    while(is.na(p) == TRUE || p <= 0 || p > 1){
      message("The p-value provided should be within 0 and 1, please try again :>")
      p <- as.numeric(readline(prompt = "Please state your desired p-value cutoff: "))
    }
  })

  ##Aspect Query
  aspect <- readline(prompt = "Which aspect(s) of the graph(s) do you want to alter? (all/legends and pvalue/title and axis/pvalue): ")
  suppressWarnings({
    while(is.na(aspect) == TRUE || tolower(aspect) != "all" && tolower(aspect) != "legends and pvalue" && tolower(aspect) != "title and axis" && tolower(aspect) != "pvalue"){
      message("Error detected, please try again :)")
      aspect <- readline(prompt = "Which aspect of the graph do you want to alter? (all/legends and pvalue/title and axis/pvalue): ")
    }
  })


  ##All Aspects
  if(tolower(aspect) == "all"){
    cat("You are required to enter the x and y coordinates for the position of the legends \n")
    cat("Please enter the comma-separated values WITH SPACES\n")
    cat("For example: 0.65, 0.75")
    suppressWarnings({
      while (TRUE) {
        lpos <- readline(prompt = "Please enter the x and y coordinates between 0 and 1 for the legend labels (Default Example: 0.65, 0.9): ")
        lpos <- unlist(strsplit(lpos, ", "))
        lpos <- as.numeric(lpos)

        # Check conditions
        if (length(lpos) != 2 || any(is.na(lpos))) {
          message("Please enter exactly two numeric values separated by a comma.")
        } else if (lpos[1] <= 0 || lpos[1] > 1 || lpos[2] <= 0 || lpos[2] > 1) {
          message("Coordinates must be between 0 and 1. Please try again.")
        } else {
          break  # Exit the loop if input is valid
        }
      }
    })
    lsize <- as.numeric(readline(prompt = "Please enter the size of the legends. (Default = 16): "))
    suppressWarnings({
      while(is.na(lsize) == TRUE){
        message("Error detected, please try again :)")
        lsize <- as.numeric(readline(prompt = "Please enter the size of the legends. (Default = 16): "))
      }
    })
    msize <- as.numeric(readline(prompt = "Please enter the size of the main title. (Default = 18): "))
    suppressWarnings({
      while(is.na(msize) == TRUE){
        message("Error detected, please try again :)")
        msize <- as.numeric(readline(prompt = "Please enter the size of the main title. (Default = 18): "))
      }
    })

    axsize <- as.numeric(readline(prompt = "Please enter the size of the axis labels. (Default = 16): "))
    suppressWarnings({
      while(is.na(axsize) == TRUE){
        message("Error detected, please try again :)")
        axsize <- as.numeric(readline(prompt = "Please enter the size of the axis labels. (Default = 16): "))
      }
    })

    psize <- as.numeric(readline(prompt = "Please enter the size of the p-value label. (Default = 16): "))
    suppressWarnings({
      while(is.na(psize) == TRUE){
        message("Error detected, please try again :)")
        psize <- as.numeric(readline(prompt = "Please enter the size of the p-value label. (Default = 16): "))
      }
    })
    pxc <- as.numeric(readline(prompt = "Please enter the x-coordinate between 0 and 1 for the p-value label: "))
    suppressWarnings({
      while(is.na(pxc) == TRUE){
        message("Error detected, please try again :)")
        pxc <- as.numeric(readline(prompt = "Please enter the x-coordinate between 0 and 1 for the p-value label. (Default = 0.45): "))
      }
    })
    cat("For the y-coordinate, the higher the value, the lower the position of the label.")
    pyc <- as.numeric(readline(prompt = "Please enter the y-coordinate for the p-value label. (Default = 38): "))
    suppressWarnings({
      while(is.na(pxc) == TRUE){
        message("Error detected, please try again :)")
        pyc <- as.numeric(readline(prompt = "Please enter the y-coordinate for the p-value label. (Default = 38): "))
      }
    })
  }

  ##Legends and P-value
  if(tolower(aspect) == "legends and pvalue"){
    cat("You are required to enter the x and y coordinates for the position of the legends \n")
    cat("Please enter the comma-separated values WITH SPACES\n")
    cat("For example: 0.65, 0.9")
    suppressWarnings({
      while (TRUE) {
        lpos <- readline(prompt = "Please enter the x and y coordinates between 0 and 1 for the legend labels (Default Example = 0.65, 0.9): ")
        lpos <- unlist(strsplit(lpos, ", "))
        lpos <- as.numeric(lpos)

        # Check conditions
        if (length(lpos) != 2 || any(is.na(lpos))) {
          message("Please enter exactly two numeric values separated by a comma.")
        } else if (lpos[1] <= 0 || lpos[1] > 1 || lpos[2] <= 0 || lpos[2] > 1) {
          message("Coordinates must be between 0 and 1. Please try again.")
        } else {
          break  # Exit the loop if input is valid
        }
      }
    })
    lsize <- as.numeric(readline(prompt = "Please enter the size of the legends. (Default = 16): "))
    suppressWarnings({
      while(is.na(lsize) == TRUE){
        message("Error detected, please try again :)")
        lsize <- as.numeric(readline(prompt = "Please enter the size of the legends. (Default = 16): "))
      }
    })
    msize <- as.numeric(18)
    axsize <- as.numeric(16)
    psize <- as.numeric(readline(prompt = "Please enter the size of the p-value label. (Default = 16): "))
    suppressWarnings({
      while(is.na(psize) == TRUE){
        message("Error detected, please try again :)")
        psize <- as.numeric(readline(prompt = "Please enter the size of the p-value label. (Default = 16): "))
      }
    })
    pxc <- as.numeric(readline(prompt = "Please enter the x-coordinate between 0 and 1 for the p-value label. (Default = 0.45): "))
    suppressWarnings({
      while(is.na(pxc) == TRUE){
        message("Error detected, please try again :)")
        pxc <- as.numeric(readline(prompt = "Please enter the x-coordinate between 0 and 1 for the p-value label: "))
      }
    })
    cat("For the y-coordinate, the higher the value, the lower the position of the label.")
    pyc <- as.numeric(readline(prompt = "Please enter the y-coordinate for the p-value label. (Default = 38): "))
    suppressWarnings({
      while(is.na(pxc) == TRUE){
        message("Error detected, please try again :)")
        pyc <- as.numeric(readline(prompt = "Please enter the y-coordinate for the p-value label. (Default = 38): "))
      }
    })
  }

  ## Title and Axis
  if(tolower(aspect) == "title and axis"){
    msize <- as.numeric(readline(prompt = "Please enter the size of the main title. (Default = 18): "))
    suppressWarnings({
      while(is.na(msize) == TRUE){
        message("Error detected, please try again :)")
        msize <- as.numeric(readline(prompt = "Please enter the size of the main title. (Default = 18): "))
      }
    })

    axsize <- as.numeric(readline(prompt = "Please enter the size of the axis labels. (Default = 16): "))
    suppressWarnings({
      while(is.na(axsize) == TRUE){
        message("Error detected, please try again :)")
        axsize <- as.numeric(readline(prompt = "Please enter the size of the axis labels. (Default = 16): "))
      }
    })

    psize <- as.numeric(16)
    pxc <- as.numeric(0.45)
    pyc <- as.numeric(38)

    lsize <- as.numeric(16)
    lpos <- "0.65, 0.9"
    lpos <- unlist(strsplit(lpos, ", "))
    lpos <- as.numeric(lpos)
  }

  ##P-value
  if(tolower(aspect) == "pvalue"){
    msize <- as.numeric(18)
    axsize <- as.numeric(16)
    lsize <- as.numeric(16)
    lpos <- "0.65, 0.9"
    lpos <- unlist(strsplit(lpos, ", "))
    lpos <- as.numeric(lpos)

    psize <- as.numeric(readline(prompt = "Please enter the size of the p-value label. (Default = 16): "))
    suppressWarnings({
      while(is.na(psize) == TRUE){
        message("Error detected, please try again :)")
        psize <- as.numeric(readline(prompt = "Please enter the size of the p-value label. (Default = 16): "))
      }
    })
    pxc <- as.numeric(readline(prompt = "Please enter the x-coordinate between 0 and 1 for the p-value label. (Default = 0.45): "))
    suppressWarnings({
      while(is.na(pxc) == TRUE){
        message("Error detected, please try again :)")
        pxc <- as.numeric(readline(prompt = "Please enter the x-coordinate between 0 and 1 for the p-value label: "))
      }
    })
    cat("For the y-coordinate, the higher the value, the lower the position of the label.")
    pyc <- as.numeric(readline(prompt = "Please enter the y-coordinate for the p-value label. (Default = 38): "))
    suppressWarnings({
      while(is.na(pxc) == TRUE){
        message("Error detected, please try again :)")
        pyc <- as.numeric(readline(prompt = "Please enter the y-coordinate for the p-value label. (Default = 38): "))
      }
    })
  }

  gtslist = grep("Expression", colnames(kmdf), value = TRUE)

  for(i in unique(1:length(gts))){
    x = gts[i]
    y = grep(x, gtslist, value = TRUE)

    fctcheck <- forcats::fct_count(kmdf[[y]])

    if(length(rownames(fctcheck))== 2){
      if(plot_type == "osc"){survdf <- survival::Surv(time = kmdf$os, event = kmdf$statusos)}
      if(plot_type == "pfs"){survdf <- survival::Surv(time = kmdf$pfs, event = kmdf$statuspfs)}

      f2 <- stats::as.formula(paste('survdf ~', paste(y)))
      pv <- survminer::surv_pvalue(survminer::surv_fit(f2, data = kmdf))[[2]]
      graph <- kmplot(f2, kmdf, x, kmlabel, ltitle, lpos, msize, lsize, psize, axsize, pxc, pyc, pv, p)
      ggplot2::ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)
    }

    if(length(rownames(fctcheck))== 3){
      if(plot_type == "osc"){
        kmdf %>% dplyr::select(id, os, statusos, dplyr::all_of(y)) %>% as.data.frame() -> kmdf2
        kmdf2 %>% dplyr::filter(kmdf2[,4] != "Medium") %>% as.data.frame() -> kmdf2
        survdf <- survival::Surv(time = kmdf2$os, event = kmdf2$statusos)
      }
      if(plot_type == "pfs"){
        kmdf %>% dplyr::select(id, pfs, statuspfs, dplyr::all_of(y)) %>% as.data.frame() -> kmdf2
        kmdf2 %>% dplyr::filter(kmdf2[,4] != "Medium") %>% as.data.frame() -> kmdf2
        survdf <- survival::Surv(time = kmdf2$pfs, event = kmdf2$statuspfs)
      }
      f2 <- stats::as.formula(paste('survdf ~', paste(y)))
      pv <- survminer::surv_pvalue(survminer::surv_fit(f2, data = kmdf2))[[2]]

      graph <- kmplot(f2, kmdf2, x, kmlabel, ltitle, lpos, msize, lsize, psize, axsize, pxc, pyc, pv, p)
      ggplot2::ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)
    }
  }
}
