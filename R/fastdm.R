#' @title Fit model using fast-dm
#' @description Specify model and fit it using fast-dm (command line)
#' @param formula model formula.
#' @param dat data.frame or named list.
#' @param method character, optimization method: Kolmogorov-Smirnov (ks), Maximum Likelihood (ml), or Chi-Square (cs).
#' @param precision numeric, approximate precision of parameter estimates in decimals.
#' @param wd_fastdm character, working directory containing fast-dm executable (defaults to 'path').
#' @param res_prefix character, prefix of temporary data files.
#' @param removeTempFiles logical, should the temporary files used by fast-dm should be kept (FALSE) or removed (TRUE). 
#' @return data.frame, parameter estimates and fit statistics for each individual.
#'
#' @examples
#' \dontrun{
#' path <- "./inst/binaries_30_2"
#' dat <- rdm(100, a=1, zr=.5, v=1.25, t0=0, N = 10)
#' out <- fastdm(~ zr*.5 + szr*0 + sv*0 + st0*0 + d*0, dat, method = "ks", precision = 3, res_prefix = "ID")
#' }
#' 
fastdm <-
function(formula=NULL, dat=NULL, method="ks", precision=2.5, wd_fastdm=path, res_prefix="data", removeTempFiles=T){
  
  #check data format
  if(class(dat) == "data.frame"){
    vars <- names(dat)
    if(sum(dat$TIME > 10) > 0) warning("implausibly long response time(s) detected. consider rescaling to seconds")
    if(sum(unique(dat$RESPONSE) %in% c(0,1)) != 2) stop("erroneous response coding: only RESPONSE = 0/1 allowed")
    write.table(dat, paste(path,"/",res_prefix,".dat",sep=""), row.names = FALSE, col.names = FALSE)
  }else if(class(dat) == "list" & !is.null(names(dat))){
    vars <- names(dat[[1]])
    for(i in 1:length(dat)){
      if(sum(dat[[i]]$TIME > 10) > 0) warning("implausibly long response time(s) detected. consider rescaling to seconds")
      if(sum(unique(dat[[i]]$RESPONSE) %in% c(0,1)) != 2) stop("erroneous response coding: correct response -> RESPONSE = 1, incorrect response -> RESPONSE = 0")
      if(sum(names(dat[[i]]) != vars) > 0) stop(paste("inconsistent variable labels in list. individual ",names(dat)[i]))
      write.table(dat[[i]], paste(path,"/",res_prefix,"_",names(dat)[i],".dat", sep=""), row.names = FALSE, col.names = FALSE)
    }
  }else{stop("data is no data.frame or named list")}
  
  #check formula
  pterms <- strsplit(as.character(formula)[2], split = c(" + "), fixed=TRUE)[[1]] #extract all parameter terms
  if(sum(grepl("*",pterms, fixed=TRUE)) > 0){
    pfix <- strsplit(pterms[grepl("*",pterms, fixed=TRUE)], " * ", fixed=TRUE) #list of parameter constraints
    for(i in 1:length(pfix)){
      if(!(pfix[[i]][1] %in% c("a","zr","v","t0","d","szr","sv","st0","p"))) stop(paste("unknown parameter label in formula:", pfix[[i]][1]))
      if(suppressWarnings(is.na(as.numeric(pfix[[i]][2])))) stop(paste("non-numeric parameter constraint in formula:", pfix[[i]][2]))
    }
  }
  if(sum(grepl(":",pterms, fixed=TRUE)) > 0){
    pdep <- strsplit(pterms[grepl(":",pterms, fixed=TRUE)], ":", fixed=TRUE) #list of parameter dependencies
    for(i in 1:length(pdep)){
      if(!(pdep[[i]][1] %in% c("a","zr","v","t0","d","szr","sv","st0","p"))) stop(paste("unknown parameter label in formula:", pdep[[i]][1]))
      if(!(pdep[[i]][2] %in% vars)) stop(paste("unknown parameter dependency in formula:", pdep[[i]][2]))
    }
  }
  
  #generate design file
  wd_temp <- getwd(); setwd(wd_fastdm) #change working directory to fast-dm
  write(c(
    paste("method",method, sep=" "),
    paste("precision",precision, sep=" ")
  ), "experiment.ctl")
  
  if(sum(grepl("*",pterms, fixed=TRUE)) > 0){
    for(i in 1:length(pfix)) write(paste("set",pfix[[i]][1],pfix[[i]][2], sep=" "), "experiment.ctl", append = TRUE)  
  }
  if(sum(grepl(":",pterms, fixed=TRUE)) > 0){
    for(i in 1:length(pdep)) write(paste("depends",pdep[[i]][1],pdep[[i]][2], sep=" "), "experiment.ctl", append = TRUE)
  }
  
  write(c(
    paste("format",paste(vars, collapse = " "), sep=" "),
    paste("load","*.dat", sep=" "),
    paste("log", "results.dat", sep=" ")
  ), "experiment.ctl", append = TRUE)
  
  #run fast-dm and remove temporary files
  system2("fast-dm", stdout="history.txt")
  out <- read.table("results.dat", header = TRUE)
  
  if(removeTempFiles){
    if(sum(file.remove("experiment.ctl","history.txt",list.files(pattern = ".dat"))) > 2) message("analysis successfully completed")
  }
  setwd(wd_temp)
  return(out)
}
