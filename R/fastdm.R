#' @title Fit diffusion model
#' @description Specify diffusion model and fit it to response time (RT) data using the fast-dm binaries (v30.2)
#' 
#' @param formula model formula (see details).
#' @param dat data.frame or named list consisting of many data.frames. Each data.frame must contain the following variables:
#'  a numeric vector TIME, containing RT from a respective trial
#'  a boolean/integer vector RESPONSE, containing the given response FALSE/0 or TRUE/1 from the respective trial
#'  
#'  
#' @param method character, optimization method: Kolmogorov-Smirnov (ks), Maximum Likelihood (ml), or Chi-Square (cs).
#' @param precision numeric, approximate precision of parameter estimates in decimals.
#' @param wd_fastdm character, working directory containing fast-dm executable (defaults to 'path').
#' @param res_prefix character, prefix of temporary data files.
#' @param removeTempFiles logical, should the temporary files used by fast-dm should be kept (FALSE) or removed (TRUE). 
#' @return data.frame, parameter estimates and fit statistics for each fitted RT distribution.
#' @examples
#' \dontrun{
#' path <- "./inst/binaries_30_2"
#' dat <- rdm(100, a=1, zr=.5, v=1.25, t0=0, N = 10)
#' out <- fastdm(~ zr*.5 + szr*0 + sv*0 + st0*0 + d*0, dat, method = "ks", precision = 3, res_prefix = "ID")
#' }
#' @export
fastdm <- function(formula=NULL, dat=NULL, method="ks", precision=2.5, wd_fastdm=path, res_prefix="data", removeTempFiles=T){
  vars <- check_dat(dat, wd_fastdm, res_prefix)
  mterms <- parse_frml(formula)
  
  #generate design file
  wd_temp <- getwd(); setwd(wd_fastdm) #change working directory to fast-dm
  write(c(
    paste("method",method, sep=" "),
    paste("precision",precision, sep=" ")
  ), "experiment.ctl")
  
  if(sum(grepl("*",mterms$pterms, fixed=TRUE)) > 0){
    for(i in 1:length(mterms$pfix)) write(paste("set",mterms$pfix[[i]][1],mterms$pfix[[i]][2], sep=" "), "experiment.ctl", append = TRUE)  
  }
  if(sum(grepl(":",mterms$pterms, fixed=TRUE)) > 0){
    for(i in 1:length(meterms$pdep)) write(paste("depends",mterms$pdep[[i]][1],mterms$pdep[[i]][2], sep=" "), "experiment.ctl", append = TRUE)
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


#' @title Fit diffusion model
#' @description Specify diffusion model and fit it to response time (RT) data using the fast-dm library (v30.2)
#' 
#' @param formula model formula (see details).
#' @param dat data.frame or named list consisting of many data.frames. Each data.frame must contain the following variables:
#'  a numeric vector TIME, containing RT from a respective trial
#'  a boolean/integer vector RESPONSE, containing the given response FALSE/0 or TRUE/1 from the respective trial
#'  
#'  
#' @param method character, optimization method: Kolmogorov-Smirnov (ks), Maximum Likelihood (ml), or Chi-Square (cs).
#' @param precision numeric, approximate precision of parameter estimates in decimals.
#' @param wd_fastdm character, working directory containing fast-dm executable (defaults to 'path').
#' @param res_prefix character, prefix of temporary data files.
#' @param removeTempFiles logical, should the temporary files used by fast-dm should be kept (FALSE) or removed (TRUE). 
#' @return data.frame, parameter estimates and fit statistics for each fitted RT distribution.
#' @examples
#' \dontrun{
#' path <- "./inst/binaries_30_2"
#' dat <- rdm(100, a=1, zr=.5, v=1.25, t0=0, N = 10)
#' out <- fastdm(~ zr*.5 + szr*0 + sv*0 + st0*0 + d*0, dat, method = "ks", precision = 3, res_prefix = "ID")
#' }
#' @useDynLib "./inst/fastdm_30_2/fastdm.dll"
#' @export
fastRdm <- function() {
  
  trialcount = 30; # helper
  
  prec = 2; # precision
  meth = "ks"; # method (ml ks cs)
  cond = c('mycond1'); # condition headers
  
  # depends sets conditional dependancies, value sets constants,
  # NULL is dont override (not passing fields should also work...)
  oparams = list(
    list(name = "t0", depends = c('mycond1'), value = NULL),
    list(name = "d", depends = NULL, value = 1.2),
    list(name = "szr", depends = NULL, value = NULL)
  ); 
  ts = rnorm(trialcount); # time vector (numeric) (length = trial count)
  resps = rnorm(trialcount)>0; # response vector (logical 0=F, 1=T) (length = trial count)
  
  # condition values (string vector) (length = condition header count * trial count)
  conds = ifelse(rnorm(trialcount * length(cond)) > 0, "a", "b");
  
  res = .Call("process_r", prec, meth, cond, oparams, ts, resps, conds);
  print(res);
}