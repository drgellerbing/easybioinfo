#' @import magrittr

##uniKM Function
uniKM = function(kmdf){
  cat("Executing Kaplan Meier Analysis .... \n")
  kmdf <- data.table::as.data.table(kmdf)

  ## Grabbing Query Genes
  gts <- getgenes()
  gts <- gsub("-", "", gts)
  gts <- gts[gts %in% colnames(kmdf)]

  ##P-value query
  suppressWarnings({
    p <- as.numeric(readline(prompt = "Please state your desired p-value cutoff: "))
    while(is.na(p) == TRUE || p <= 0 || p > 1){
      message("The p-value provided should be within 0 and 1, please try again :>")
      p <- as.numeric(readline(prompt = "Please state your desired p-value cutoff: "))
    }
  })

  ## Plot Types
  if(all(c("os", "statusos", "pfs", "statuspfs") %in% colnames(kmdf))){plot_para <- list(length = 2, type = c("osc", "pfs"))}
  else if (all(c("os", "statusos") %in% colnames(kmdf))){plot_para <- list(length = 1, type = "osc")}
  else if (all(c("pfs", "statuspfs") %in% colnames(kmdf))){plot_para <- list(length = 1, type = "pfs")}

  plot <- readline(prompt = "Do you want to plot the survival curves? (yes/no): ")
  while(tolower(plot) != "yes" && tolower(plot) != "no"){
    message("Answer not given, please try again :>")
    plot <- readline(prompt = "Do you want to plot the survival curves? (yes/no): ")
  }

  # Schoenfield Query
  message("Schoenfield residuals plots will be placed under a folder named 'PH Check'")
  plotsr <- readline(prompt = "Do you want the Schoenfeld residuals plot? (yes/no): ")
  while(tolower(plotsr) != "yes" && tolower(plotsr) != "no"){
    message("Answer not given, please try again :>")
    plotsr <- readline(prompt = "Do you want the Schoenfeld residuals plot? (yes/no): ")
  }

  if(tolower(plot) == "yes" || tolower(plot) == "y"){
    if(length(gts) > 10){
      plot_sig <- readline(prompt = "Do you want to plot only the significant results? (yes/no): ")
      while(tolower(plot_sig) != "yes" && tolower(plot_sig) != "no" && tolower(plot_sig) != "y" && tolower(plot_sig) != "n"){
        message("Answer not given, please try again :>")
        plot_sig <- readline(prompt = "Do you want to plot only the significant results? (yes/no): ")
      }
    } else{
      plot_sig <- "no"
    }

    ## Folder Creation
    if(length(gts) > 10){
      dir <- readline(prompt = "Do you want to create a new folder for your plots? (yes/no): ")
      while(tolower(dir) != "yes" && tolower(dir) != "y" && tolower(dir) != "no" && tolower(dir) != "n"){
        message("Answer not given, please try again :>")
        dir <- readline(prompt = "Do you want to create a new folder for your plots & results? (yes/no): ")
      }
      if(tolower(dir) == "yes" || tolower(dir) == "y"){dir_n <- createdir()}
    }

    if(plot_para$length == 2 && dir == "yes"){
      dir.create(paste(dir_n, "/", "OS", sep = ""))
      dir.create(paste(dir_n, "/", "PFS", sep = ""))
      dir_n <- c(paste(dir_n, "/", "OS", sep = ""), paste(dir_n, "/", "PFS", sep = ""))
    }
  } else{dir <- "no"}

  kmdflist <- vector("list", length = plot_para$length)

  ## Complete Analysis
  for(pl in 1:plot_para$length){
    plot_type = plot_para$type[pl]
    kmdf2 <- kmdf
    if(dir == "yes"){f_dir_n = dir_n[pl]}

    #Defining Label Title, Survival Calculations & Establishing Folders
    if(plot_type == "osc"){
      message("Now processing overall survival result :)")
      kmlabel <- "Survival Probability"
      ltitle <- "Median Survival"
      if(tolower(plotsr) == "yes" || tolower(plotsr) == "y"){
        if(tolower(dir) == "yes" || tolower(dir) == "y"){srfname <- paste(f_dir_n, "/", "PH_Check_OS", sep = "")}
        if(tolower(dir) == "no" || tolower(dir) == "n"){srfname <- "PH_Check_OS"}
        if(!file.exists(srfname)){dir.create(srfname)}
        }
    }

    if(plot_type == "pfs"){
      message("Now processing progression-free survival result :)")
      kmlabel <- "Progression-Free Survival Probability"
      ltitle <- "Median PFS"
      if(tolower(plotsr) == "yes" || tolower(plotsr) == "y"){
        if(tolower(dir) == "yes" || tolower(dir) == "y"){srfname <- paste(f_dir_n, "/", "PH_Check_PFS", sep = "")}
        if(tolower(dir) == "no" || tolower(dir) == "n"){srfname <- "PH_Check_PFS"}
        if(!file.exists(srfname)){dir.create(srfname)}
        }
    }

    message("Do you want to set the cutoff based off the median, quartile, custom %, z-score, or individual optimal cutoffs?")
    cutoff <- readline(prompt = "Please enter your cutoff method (median/quartile/custom/zscore/optimal): ")
    while(tolower(cutoff) != "median" && tolower(cutoff) != "quartile" && tolower(cutoff) != "custom" && tolower(cutoff) != "zscore" && tolower(cutoff) != "optimal"){
      message("Answer not given, please try again :)")
      cutoff <- readline(prompt = "Please enter your cutoff method (median/quartile/custom/zscore/optimal): ")
    }

    ## Establishing custom/z-score cutoff
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

    if(tolower(cutoff) == "zscore"){zcut <- as.numeric(readline(prompt = "Please enter your desired z-score cutoff: "))}



    ## Segregation loop
    for (y in unique(1:length(gts))){
      z = gts[y]
      name = paste(z, 'Expression', sep = "")
      tryCatch({
        kmdf2 %>% dplyr::select(dplyr::all_of(z)) %>% unlist() -> c
        #Median Cutoff
        if(tolower(cutoff) == "median"){
          mdn = stats::median(c)
          kmdf2 %<>% dplyr::mutate("{name}" := as.factor(
            ifelse(
              c >= mdn,
              'High',
              'Low'
            )))
        } #End of median cutoff loop

        # Quartile/custom cutoff
        if(tolower(cutoff) == "quartile" || tolower(cutoff) == "custom"){
          if(tolower(cutoff) == "quartile"){quartile <- stats::quantile(c, probs = c(0.25, 0.75))}
          if(tolower(cutoff) == "custom"){quartile <- stats::quantile(c, probs = c(qval_low, qval_high))}

          kmdf2 %<>% dplyr::mutate("{name}" := as.factor(
            ifelse(
              c >= quartile[[2]],
              'High',
              ifelse(
                c <= quartile[[1]],
                "Low",
                "Medium"
              )
            )))
        } #End of quartile/custom loop

        # Z-Score cutoff
        if(tolower(cutoff) == "zscore"){
          zsc = scale(c)
          kmdf2 %<>% dplyr::mutate("{name}" := as.factor(
            ifelse(
              zsc >= zcut,
              'High',
              ifelse(
                zsc <= -abs(zcut),
                "Low",
                "Medium"
              )
            )))
        } # End of zscore loop

        if(tolower(cutoff) == "optimal"){
          ## Backbone of final dataframes
          cutoffcol = c("Gene", "Cutoff")
          cuttable = data.frame(matrix(nrow=0, ncol = length(cutoffcol)))
          colnames(cuttable) = cutoffcol

          if(plot_type == "osc"){survcutoff <- survival::Surv(time = kmdf2$os, event = kmdf2$statusos)}
          if(plot_type == "pfs"){survcutoff <- survival::Surv(time = kmdf2$pfs, event = kmdf2$statuspfs)}
          fc <- stats::as.formula(paste('survcutoff ~', paste(z)))

          tryCatch({
            mscutoff <- maxstat::maxstat.test(fc, data = kmdf2, smethod = "LogRank")
            o_cutoff <- mscutoff$estimate %>% as.double()
            cutdf <- data.frame(z, o_cutoff)
            cuttable <- rbind(cuttable, cutdf)

            kmdf2 %<>% dplyr::mutate("{name}" := as.factor(
              ifelse(
                c >= o_cutoff,
                'High',
                'Low'
              )))
          }, error = function(e) {
            cat(z, "contains non-zero cell counts within the specified proportion range, removed from the analysis \n")
          })
        } # End of optimal cutoff loop
      }, error = function(e){
        cat(z, "was not found in your provided dataframe :) \n")
      }) # End of gene check loop
    } # End of segregation loop

    #Removing Genes That Cannot be Segregated
    rm_genes <- vector("character", length = 0)
    gtslist1 = grep("Expression", colnames(kmdf2), value = TRUE)

    if(tolower(cutoff) == "optimal" || tolower(cutoff) == "median"){
      for (o in gtslist1){
        fctcheck <- forcats::fct_count(kmdf2[[o]])
        if(length(rownames(fctcheck)) <2){
          rm_genes <- c(rm_genes, o)
          kmdf2 %<>% dplyr::select(-dplyr::all_of(o))
        }
        if(length(rownames(fctcheck))== 2){kmdf2[[o]] <- stats::relevel(factor(kmdf2[[o]]), ref = 'Low')}
      }
    }

    if(tolower(cutoff) == "quartile" || tolower(cutoff) == "custom" || tolower(cutoff) == "zscore"){
      for (o in gtslist1){
        fctcheck <- forcats::fct_count(kmdf2[[o]])
        if(length(rownames(fctcheck)) <3){
          rm_genes <- c(rm_genes, o)
          kmdf2 %<>% dplyr::select(-dplyr::all_of(o))
        }
        if(length(rownames(fctcheck))== 3){kmdf2[[o]] <- stats::relevel(factor(kmdf2[[o]]), ref = 'Low')}
      }
    }

    rm_genes <- gsub("Expression", "", rm_genes)
    if(length(rm_genes) > 0){
      for(m in 1:length(rm_genes)){cat("Gene", rm_genes[m], "was removed from the analysis. Unable to segregate expression into two groups\n")}
    }

    ##Final Analysis

    # Backbone of Final Dataframe
    unicolumn = c('Gene', 'Coefficient', 'HR', '95%CI', 'Pvalue', 'PHCheck')
    if(plot_type == "osc"){
      kmtableos = data.frame(matrix(nrow=0, ncol = length(unicolumn)))
      colnames(kmtableos) = unicolumn
    }
    if(plot_type == "pfs"){
      kmtablepfs = data.frame(matrix(nrow=0, ncol = length(unicolumn)))
      colnames(kmtablepfs) = unicolumn
    }

    gtslist = grep("Expression", colnames(kmdf2), value = TRUE)
    gts <- gsub("Expression", "", gtslist)

    # Final analysis loop
    for (i in unique(1:length(gts))) {
      x = gts[i]
      y = gtslist[i]

      #Trycatch Loop
      tryCatch({
        cat("Now processing survival analysis for", x, "\n" )
        f2 <- stats::as.formula(paste('survdf ~', paste(y)))
        if(plot_type == "osc"){
          kmdf2 %>% dplyr::select(id, os, vitalos, statusos, dplyr::all_of(y)) %>% data.table::as.data.table() -> kmdf3
          if(tolower(cutoff) == "quartile" || tolower(cutoff) == "custom" || tolower(cutoff) == "zscore"){kmdf3 %<>% dplyr::filter(kmdf3[,5] != "Medium")}
          survdf <- survival::Surv(time = kmdf3$os, event = kmdf3$statusos)
        }
        if(plot_type == "pfs"){
          kmdf2 %>% dplyr::select(id, pfs, vitalpfs, statuspfs, dplyr::all_of(y)) %>% data.table::as.data.table() -> kmdf3
          if(tolower(cutoff) == "quartile" || tolower(cutoff) == "custom" || tolower(cutoff) == "zscore"){kmdf3 %<>% dplyr::filter(kmdf3[,5] != "Medium")}
          survdf <- survival::Surv(time = kmdf3$pfs, event = kmdf3$statuspfs)
        }
        unidf2 <- onegenecox(f2, x, kmdf3)
        zph = survival::cox.zph(survival::coxph(f2, data = kmdf3))
        unidf2 <- data.frame(unidf2, zph$table[1,3])
        colnames(unidf2) <- unicolumn
        if(plot_type == "osc"){kmtableos <- rbind(kmtableos, unidf2)}
        if(plot_type == "pfs"){kmtablepfs <- rbind(kmtablepfs, unidf2)}

        ## KM Plot Loop
        if(tolower(plot) == "yes" || tolower(plot) == "y"){
          pv <- survminer::surv_pvalue(survminer::surv_fit(f2, data = kmdf3))[[2]]
          if(tolower(plot_sig) == "yes"){
            if(pv <= p){
              graph <- kmplot(f2, kmdf3, x, kmlabel, ltitle, pv, p)

              if(tolower(dir) == "yes" || tolower(dir) == "y"){ggplot2::ggsave(paste(f_dir_n, "/", x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
              if(tolower(dir) == "no" || tolower(dir) =="n"){ggplot2::ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
            }
          }

          if(tolower(plot_sig) == "no" || tolower(plot_sig) == "n"){
            graph <- kmplot(f2, kmdf3, x, kmlabel, ltitle, pv, p)

            if(tolower(dir) == "yes" || tolower(dir) == "y"){ggplot2::ggsave(paste(f_dir_n, "/", x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
            if(tolower(dir) == "no" || tolower(dir) =="n"){ggplot2::ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)}
          }
        } #End of KM plot loop

        #Schoenfield Plot loop
        if(tolower(plotsr) == "yes" || tolower(plotsr) == "y"){
          if(tolower(plot_sig) == "yes" || tolower(plot_sig) == "y"){
            pv_sr <- survminer::surv_pvalue(survminer::surv_fit(f2, data = kmdf3))[[2]]
            if(pv_sr <= p){
              grDevices::png(paste(srfname, "/", x, ".png", sep = ""), res = 300, units = "cm", width = 15, height = 15)
              plot(zph)
              grDevices::dev.off()
            }
          }
          if(tolower(plot_sig) == "no" || tolower(plot_sig) == "n"){
            grDevices::png(paste(srfname, "/", x, ".png", sep = ""), res = 300, units = "cm", width = 15, height = 15)
            plot(zph)
            grDevices::dev.off()
          }
        } #End of schoenfield plot loop
      }, error = function(e){
        warning(paste("Gene", x, "not found in dataframe. Skipping..."))
      })
    } #End of final analysis loop
    kmdflist[[pl]] <- kmdf2
  } #End of complete analysis loop

  ## Saving Results
  if(plot_para$length == 1){
    #OS Loop
    if(plot_para$type == "osc"){
      kmtableos %>% dplyr::filter(kmtableos$Pvalue < p) %>% as.data.frame() -> sigKMos
      list1 <- (list(kmtableos, sigKMos, kmdflist[[1]]))
      names(list1) <- c("kmtableos", "sigKMos", "transformeddf_os")
    }#End of os loop

    #PFS Loop
    if(plot_para$type == "pfs"){
      kmtablepfs %>% dplyr::filter(kmtablepfs$Pvalue < p) %>% as.data.frame() -> sigKMpfs
      list1 <- (list(kmtablepfs, sigKMpfs, kmdflist[[1]]))
      names(list1) <- c("kmtablepfs", "sigKMpfs", "transformeddf_pfs")
    }#End of pfs loop
  }

  if(plot_para$length == 2){
    kmtableos %>% dplyr::filter(kmtableos$Pvalue < p) %>% as.data.frame() -> sigKMos
    kmtablepfs %>% dplyr::filter(kmtablepfs$Pvalue < p) %>% as.data.frame() -> sigKMpfs
    list1 <- list(kmtableos, sigKMos, kmtablepfs, sigKMpfs, kmdflist[[1]], kmdflist[[2]])
    names(list1) <- c("kmtableos", "sigKMos", "kmtablepfs", "sigKMpfs", "transformeddf_os", "transformeddf_pfs")
  }

  res_return <- readline(prompt = "Do you only want to save the results? (yes/no): ")
  if(tolower(res_return) == "yes" || tolower(res_return) == "y"){filen <- readline(prompt = "Please provide a name for your results: ")}

  if(tolower(res_return) == "yes"){
    # Only OS or only PFS
    if(plot_para$length == 1){
      res1 <- as.data.frame(list1[[1]])
      res2 <- as.data.frame(list1[[2]])
      name1 <- paste("./", filen, "_full_uniKM.xlsx", sep = "")
      name2 <- paste("./", filen, "_sig_uniKM.xlsx", sep = "")

      writexl::write_xlsx(res1, name1)
      writexl::write_xlsx(res2, name2)
    } #End of loop

    # Both OS and PFS
    if(plot_para$length == 2){
      writexl::write_xlsx(as.data.frame(list1[[1]]), paste("./", filen, "_full_OS_uniKM.xlsx", sep = ""))
      writexl::write_xlsx(as.data.frame(list1[[2]]), paste("./", filen, "_sig_OS_uniKM.xlsx", sep = ""))
      writexl::write_xlsx(as.data.frame(list1[[3]]), paste("./", filen, "_full_PFS_uniKM.xlsx", sep = ""))
      writexl::write_xlsx(as.data.frame(list1[[4]]), paste("./", filen, "_sig_PFS_uniKM.xlsx", sep = ""))
    } #End of loop

    if(tolower(cutoff) == "optimal"){
      namec <- paste("./", filen, "_gene_cutoff.xlsx", sep = "")
      writexl::write_xlsx(cuttable, namec)}
  } #End of loop

  cat("All done :> \n")
  return(list1)
}

##Customisable Survival Curve Function
customisesurvivalplot = function(unidf){
  if(length(unidf) == 3){
    kmdf <- unidf[[3]]
    if(all(c("os", "statusos") %in% colnames(kmdf))){plot_type == "osc"}
    else if(all(c("pfs", "statuspfs") %in% colnames(kmdf))){plot_type == "pfs"}
  }

  if(length(unidf) == 6){
    plot_type <- readline(prompt = "Are you looking to customise your overall survival/progression-free survival plots? (osc/pfs): ")
    if(plot_type == "osc"){kmdf <- unidf[[5]]}
    if(plot_type == "pfs"){kmdf <- unidf[[6]]}
  }

  gts <- getgenes()
  gts <- gts[gts %in% colnames(kmdf)]

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

  ##Legends and P-value
  if(tolower(aspect) == "legends and pvalue" || tolower(aspect) == "all"){

    #lpos Input loop
    lpos_success <- FALSE
    cat("You are required to enter the x and y coordinates for the position of the legends \n")
    cat("Please enter the comma-separated values WITH SPACES\n")
    cat("For example: 0.65, 0.75")
    while(lpos_success == FALSE){
      lpos <- readline(prompt = "Please enter the x and y coordinates between 0 and 1 for the legend labels (Default Example: 0.65, 0.9): ")

      #TryCatch Loop
      tryCatch({
        lpos <- gsub("\\s*,\\s*", ",", trimws(lpos))
        lpos <- gsub(",", ", ", lpos)

        lpos <- unlist(strsplit(lpos, ", "))
        lpos <- as.numeric(lpos)

        if(any(is.na(lpos))){stop("Error found in the input.")}
        if (length(lpos) != 2 || any(is.na(lpos))) {
          message("Please enter exactly two numeric values separated by a comma.")
          stop()
        }
        if (lpos[1] <= 0 || lpos[1] > 1 || lpos[2] <= 0 || lpos[2] > 1) {
          message("Coordinates must be between 0 and 1. Please try again.")
          stop()
        }

        lpos_success <- TRUE

      }, error = function(e){
        lpos_success <- FALSE
        message("Error detected. Please try again :) \n")
      })
    }
    #End of lpos input loop

    lsize <- as.numeric(readline(prompt = "Please enter the size of the legends. (Default = 16): "))
    suppressWarnings({
      while(is.na(lsize) == TRUE){
        message("Error detected, please try again :)")
        lsize <- as.numeric(readline(prompt = "Please enter the size of the legends. (Default = 16): "))
      }
    })

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
  if(tolower(aspect) == "title and axis" || tolower(aspect) == "all"){
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
  }

  ##P-value
  if(tolower(aspect) == "pvalue"){
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
      graph <- kmplot(f2, kmdf, x, kmlabel, ltitle, pv, p, lpos, msize, lsize, psize, axsize, pxc, pyc)
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

      graph <- kmplot(f2, kmdf, x, kmlabel, ltitle, pv, p, lpos, msize, lsize, psize, axsize, pxc, pyc)
      ggplot2::ggsave(file = paste(x, ".png", sep = ""), width = 8, height = 9, dpi = 500)
    }
  }
}
