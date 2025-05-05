#' @import magrittr
#' @name multicox
#' @title Multivariate Cox Proprotional Hazard Analysis
#' @description
#' Performs multivariate cox proportional hazard analysis for overall survival or progression-free survival through a series of prompts.
#' The function will prompt you for a list of genes.
#' Survival can be calculated based on four different cutoffs: median, z-scores, quartiles, custom quartiles or optimal.
#' Optimal cutoffs are calculated by performing a test of independence of response using maximally selected rank statistics
#' @param kmdf Ideally, dataframe produced from the tidykmdata function. Dataframe should contain "os"/"pfs" columns for OS/PFS values, "statusos"/"statuspfs" with "1 or 0" as events, as well as expression values for individual genes in each column.
#' @returns A list containing the class coxph representing the fit, the full cox analysis result, the significant cox analysis result and the transformed input dataframe used for the analysis.
#' Also produces a forest plot that is saved into a .png file. Forest plot can be further customised using plotforest function.
#' @examples
#' df <- easybioinfo::kmexpr
#' md <- easybioinfo::kmclinical
#' kmdf <- tidykmdata(df, md)
#' unidf <- uniKM(kmdf)
#' cox <- multicox(kmdf)
#' @export

multicox = function(kmdf){
  cat("Executing Multivariate Cox Analysis ... \n")

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

  if(all(c("os", "statusos", "pfs", "statuspfs") %in% colnames(kmdf))){calc_para <- list(length = 2, type = c("osc", "pfs"))}
  else if (all(c("os", "statusos") %in% colnames(kmdf))){calc_para <- list(length = 1, type = "osc")}
  else if (all(c("pfs", "statuspfs") %in% colnames(kmdf))){calc_para <- list(length = 1, type = "pfs")}

  kmdflist <- vector("list", length = calc_para$length)

  #Converting column names in dataframe into lowercase
  for(pl in 1:calc_para$length){
    calc_type = calc_para$type[pl]
    if(calc_type == "osc"){message("Now processing overall survival result :)")}
    if(calc_type == "pfs"){message("Now processing progression-free survival result :)")}

    gts <- getgenes()
    gts <- gsub("-", "", gts)
    kmdf <- as.data.frame(kmdf)
    gts <- gts[gts %in% colnames(kmdf)]

    kmdf2 <- kmdf

    cutoff <- readline(prompt = "Please enter your cutoff method (median/quartile/custom/zscore/optimal/continuous): ")
    while(tolower(cutoff) != "median" && tolower(cutoff) != "quartile" && tolower(cutoff) != "custom" && tolower(cutoff) != "zscore" && tolower(cutoff) != "optimal" && tolower(cutoff) != "continuous"){
      message("Answer not given, please try again :)")
      cutoff <- readline(prompt = "Please enter your cutoff method (median/quartile/custom/zscore/optimal/continuous): ")
    }

    if(tolower(cutoff) != "continuous"){
      for(col_name in colnames(kmdf2)) {
        if(tolower(col_name) %in% tolower(gts)) {
          colnames(kmdf2)[colnames(kmdf2) == col_name] <- tolower(col_name)
        }
      }
    }

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

    #Segregating Genes based on chosen method
    rm_genes <- vector("character", length = 0)
    for (i in unique(1:length(gts))){
      z = gts[i]
      l = tolower(z)
      name = paste(z)

      #Loop to double check if genes are present in the query dataframe
      tryCatch({
        cat("Calculating multivariate survival for", z, "\n")
        kmdf2 %>% dplyr::select(dplyr::all_of(l)) %>% unlist() -> c

        #Median segregation loop
        if(tolower(cutoff) == "median"){
          mdn = stats::median(c)
          kmdf2 %<>% dplyr::mutate("{name}" := as.factor(
            ifelse(
              c >= mdn,
              'High',
              'Low'
            )))
        } # End of loop

        # Quartile or custom segregation loop
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
        } # End of quartile/custom loop

        #Z-score segregation loop
        if(tolower(cutoff) == "zscore"){
          zcut <- as.numeric(readline(prompt = "Please enter your desired z-score cutoff: "))
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
        } # End of z-score loop

        #Optimal cutoff segregation loop
        if(tolower(cutoff) == "optimal"){
          if(calc_type == "osc"){survdfop <- survival::Surv(time = kmdf2$os, event = kmdf2$statusos)}
          if(calc_type == "pfs"){survdfop <- survival::Surv(time = kmdf2$pfs, event = kmdf2$statuspfs)}

          fc <- stats::as.formula(paste('survdfop ~', paste(l)))
          tryCatch({
            mscutoff <- maxstat::maxstat.test(fc, data = kmdf2, smethod = "LogRank")
            o_cutoff <- mscutoff$estimate %>%as.double()

            kmdf2 %<>% dplyr::mutate("{name}" := as.factor(
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
        } # End of optimal loop
      }, error = function(e){
        warning(paste("Column '", z, "' not found in your provided dataframe. Skipping...", sep = ""))
      }) #End of loop

      fctcheck <- forcats::fct_count(kmdf2[[name]])
      # Fact Check loop
      if(tolower(cutoff) == "median" || tolower(cutoff) == "optimal"){
        if(length(rownames(fctcheck)) <2){
          rm_genes <- c(rm_genes, name)
          kmdf2 %<>% dplyr::select(-dplyr::all_of(name))
        }
        if(length(rownames(fctcheck))== 2){kmdf2[[name]] <- stats::relevel(factor(kmdf2[[name]]), ref = 'Low')}
      } # End of fact check loop

      # Fact check loop for custom, quartile and zscore cutoffs
      if(tolower(cutoff) == "quartile" || tolower(cutoff) == "custom" || tolower(cutoff) == "zscore"){
        if(length(rownames(fctcheck)) <3){
          rm_genes <- c(rm_genes, name)
          kmdf2 %<>% dplyr::select(-dplyr::all_of(name))
        }
        if(length(rownames(fctcheck))== 3){kmdf2[[name]] <- factor(kmdf2[[name]], levels = c("Medium", "Low", "High"))}
      } # End of fact check loop

    } # End of analysis loop

    #Reporting Errors
    if(length(rm_genes) > 0){
      for(m in 1:length(rm_genes)){cat("Gene", rm_genes[m], "was removed from the analysis. Unable to segregate expression into required groups\n")}
    }

    #Calculating Cox
    if(calc_type == "osc"){survdf <- survival::Surv(time = kmdf2$os, event = kmdf2$statusos)}
    if(calc_type == "pfs"){survdf <- survival::Surv(time = kmdf2$pfs, event = kmdf2$statuspfs)}
    f1 <- stats::as.formula(paste('survdf ~', paste(gts, collapse = "+")))
    cox = survival::coxph(f1, data = kmdf2)
    zph = survival::cox.zph(cox)

    unicolumn = c('Gene', 'Hazard Ratio', 'CI', 'P-value', 'PH Check')
    coxtable = data.frame(matrix(nrow=0, ncol = length(unicolumn)))

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
    if(calc_type == "osc"){
      coxos <- cox
      coxtableos <- coxtable
      coxtableos %>% dplyr::filter(coxtableos$`P-value` <= p) %>%as.data.frame() -> sigcoxos
    }

    if(calc_type == "pfs"){
      coxpfs <- cox
      coxtablepfs <- coxtable
      coxtablepfs %>% dplyr::filter(coxtablepfs$`P-value` <= p) %>% as.data.frame() -> sigcoxpfs
    }

    graph <- survminer::ggforest(cox, data = kmdf2, refLabel = "Reference") %>% print()
    kmdflist[[pl]] <- kmdf2
  } #End of entire analysis loop

  filen <- readline(prompt = "Please provide a name for your results: ")
  if(calc_para$length == 1){
    finalkmdf <- kmdflist[[1]]
    if(calc_para$type == "osc"){
      list1 <- (list(coxos, coxtableos, sigcoxos, finalkmdf))
      names(list1) <- c("cox_os", "full_cox", "sig_cox", "kmdf_os")
    }

    if(calc_para$type == "pfs"){
      list1 <- (list(coxpfs, coxtablepfs, sigcoxpfs, finalkmdf))
      names(list1) <- c("cox_pfs", "full_cox", "sig_cox", "kmdf_os")
    }
    res1 <- as.data.frame(list1[[2]])
    name1 <- paste("./", filen, "_multicox.xlsx", sep = "")

    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Data")
    openxlsx::writeData(wb, "Data", res1)
    openxlsx::conditionalFormatting(wb, "Data", rows = 1:nrow(res1), cols = 1, style = openxlsx::createStyle(textDecoration = "bold"), rule = 'D2 < 0.05')

    openxlsx::saveWorkbook(wb, name1, overwrite = TRUE)
  }

  if(calc_para$length == 2){
    finalkmdfos <- kmdflist[[1]]
    finalkmdfpfs <- kmdflist[[2]]
    list1 <- list(coxos, coxpfs, coxtableos, sigcoxos, coxtablepfs, sigcoxpfs, finalkmdfos, finalkmdfpfs)
    names(list1) <- c("cox_os", "cox_pfs","fullcox_os", "sigcox_os", "fullcox_pfs", "sigcox_pfs", "kmdf_os", "kmdf_pfs")

    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "OS")
    openxlsx::addWorksheet(wb, "PFS")
    openxlsx::writeData(wb, "OS", as.data.frame(list1[[3]]))
    openxlsx::writeData(wb, "PFS", as.data.frame(list1[[5]]))
    openxlsx::conditionalFormatting(wb, "OS", rows = 1:nrow(as.data.frame(list1[[3]])) + 1, cols = 1, style = openxlsx::createStyle(textDecoration = "bold"), rule = 'D2 < 0.05')
    openxlsx::conditionalFormatting(wb, "PFS", rows = 1:nrow(as.data.frame(list1[[5]])) + 1, cols = 1, style = openxlsx::createStyle(textDecoration = "bold"), rule = 'D2 < 0.05')
    openxlsx::saveWorkbook(wb, paste("./", filen, "_cox.xlsx", sep = ""), overwrite = TRUE)
  }
  cat("Analysis Complete :> Hope you get something :> \n")
  return(list1)
}

#' @name plotforest
#' @title Customise Forest Plot
#' @description
#' A helper function to allow to produce customised forest plots. 
#' The title, font size, reference label, and the number of digits for estimates and p-values in the plot can all be customised.
#' @param cox The list from the multicox function
#' @returns A forest plot that can be saved into a .png file
#' @examples
#' df <- easybioinfo::kmexpr
#' md <- easybioinfo::kmclinical
#' kmdf <- tidykmdata(df, md)
#' unidf <- uniKM(kmdf)
#' cox <- multicox(kmdf)
#' plotforest(cox)
#' @export

plotforest = function(cox){
  if(length(cox) == 4){
    coxdf <- cox[[1]]
    kmdf <- cox[[4]]
  }

  if(length(cox) == 8){
    type <- readline(readline(prompt = "Are you looking to customise the forest plot of your overall survival or progression-free survival result? (osc/pfs): "))
    if(type == "osc"){
      coxdf <- cox[[1]]
      kmdf <- cox[[7]]
    }
    if(type == "pfs"){
      coxdf <- cox[[2]]
      kmdf <- cox[[8]]
    }
  }


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
