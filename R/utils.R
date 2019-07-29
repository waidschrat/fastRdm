#' @title check data format
check_dat <- function(dat, wd_fastdm, res_prefix){
  if(class(dat) == "data.frame"){
    vars <- names(dat)
    if(sum(c("TIME","RESPONSE") %in% vars) != 2)
      stop("dat does not contain both variables 'TIME' and 'RESPONSE'")
    if(sum(dat$TIME > 10) > 0)
      warning("implausibly long response time(s) detected. consider rescaling to seconds")
    if(sum(unique(dat$RESPONSE) %in% c(0,1)) != 2)
      stop("erroneous response coding: only RESPONSE = 0/1 allowed")
    write.table(dat, paste(wd_fastdm,"/",res_prefix,".dat",sep=""), row.names = FALSE, col.names = FALSE)
    return(vars)
  }else if(class(dat) == "list" & !is.null(names(dat))){
    for(i in 1:length(dat)){
      vars <- names(dat[[i]])
      if(sum(c("TIME","RESPONSE") %in% vars) != 2)
        stop(paste0("dat['",names(dat)[i],"'] does not contain both variables 'TIME' and 'RESPONSE'"))
      if(sum(dat[[i]]$TIME > 10) > 0)
        warning(paste0("implausibly long response time(s) detected in dat['",names(dat)[i],"']. consider rescaling to seconds"))
      if(sum(unique(dat[[i]]$RESPONSE) %in% c(0,1)) != 2)
        stop(paste0("erroneous response coding in dat['",names(dat)[i],"']: only RESPONSE = 0/1 allowed"))
      write.table(dat[[i]], paste(wd_fastdm,"/",res_prefix,"_",names(dat)[i],".dat", sep=""), row.names = FALSE, col.names = FALSE)
      if(i == length(dat)) return(vars)
    }
  }else{stop("dat is no data.frame or named list")}
}


#' @title check and parse model formula
parse_frml <- function(formula){
  pterms <- strsplit(as.character(formula)[2], split = c(" + "), fixed=TRUE)[[1]] #extract all parameter terms
  res <- list("pterms"=pterms)
  
  if(sum(grepl("*",pterms, fixed=TRUE)) > 0){
    pfix <- strsplit(pterms[grepl("*",pterms, fixed=TRUE)], " * ", fixed=TRUE) #list of parameter constraints
    for(i in 1:length(pfix)){
      if(!(pfix[[i]][1] %in% c("a","zr","v","t0","d","szr","sv","st0","p"))) stop(paste("unknown parameter label in formula:", pfix[[i]][1]))
      if(suppressWarnings(is.na(as.numeric(pfix[[i]][2])))) stop(paste("non-numeric parameter constraint in formula:", pfix[[i]][2]))
    }
    res[["pfix"]] <- pfix
  }
  if(sum(grepl(":",pterms, fixed=TRUE)) > 0){
    pdep <- strsplit(pterms[grepl(":",pterms, fixed=TRUE)], ":", fixed=TRUE) #list of parameter dependencies
    for(i in 1:length(pdep)){
      if(!(pdep[[i]][1] %in% c("a","zr","v","t0","d","szr","sv","st0","p"))) stop(paste("unknown parameter label in formula:", pdep[[i]][1]))
      if(!(pdep[[i]][2] %in% vars)) stop(paste("unknown parameter dependency in formula:", pdep[[i]][2]))
    }
    res[["pdep"]] <- pdep
  }
  
  return(res)
}
