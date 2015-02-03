################################################
##
## @script: RNASeq_QC_lib.R
##
## @author: Joseph Tran (Joseph.Tran@versailles.inra.fr)
##
## @version: 0.0.1.0
##
## @date: 2015-01-23
##
## @description: This R library contains functions for running RNASeq Quality Control (QC) on gene count dataset and mapping statistics (BBRIC). 
##
## @copyright: Copyright (c) 2015 by Joseph Tran. This work is made available under the terms of the Creative Commons Attribution-ShareAlike 3.0 license, http://creativecommons.org/licenses/by-sa/3.0/.
##
###############################################

version <- "0.0.1.0"

copyright <- "Copyright (c) 2015 by Joseph Tran. This work is made available under the terms of the Creative Commons Attribution-ShareAlike 3.0 license, http://creativecommons.org/licenses/by-sa/3.0/."

###
### Dependencies ###
###

### check for installed package
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])

### testthat
if (is.installed('testthat'))
{
  suppressPackageStartupMessages(require('testthat'));
} else
{
  install.packages('testthat');
  suppressPackageStartupMessages(library('testthat'));
}

### log4r
if (is.installed('log4r'))
{
  suppressPackageStartupMessages(require('log4r'));
} else
{
  install.packages('log4r');
  suppressPackageStartupMessages(library('log4r'));
}

### getopt
if (is.installed('getopt'))
{
  suppressPackageStartupMessages(require('getopt'));
} else
{
  install.packages('getopt');
  suppressPackageStartupMessages(library('getopt'));
}

### optparse
if (is.installed('optparse'))
{
  suppressPackageStartupMessages(require('optparse'));
} else
{
  install.packages('optparse');
  suppressPackageStartupMessages(library('optparse'));
}

### knitr
if (is.installed('knitr'))
{
  suppressPackageStartupMessages(require('knitr'));
} else
{
  install.packages('knitr');
  suppressPackageStartupMessages(library('knitr'));
}

### markdown
if (is.installed('markdown'))
{
  suppressPackageStartupMessages(require('markdown'));
} else
{
  install.packages('markdown');
  suppressPackageStartupMessages(library('markdown'));
}

### rmarkdown
if (is.installed('rmarkdown'))
{
  suppressPackageStartupMessages(require('rmarkdown'));
} else
{
  install.packages('rmarkdown');
  suppressPackageStartupMessages(library('rmarkdown'));
}

###
### FUNCTIONS ###
###

###
### Capturing warnings/errors
###
tryCatch.W.E <- function(expr)
{
  W <- NULL
  w.handler <- function(w)
  { # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}

###
### Load count data
###
loadCountData <- function(file)
{
  # load count data
  warn_err <- tryCatch.W.E(read.table(opt$count, header=TRUE, sep="\t"))
  
  if (is.null(warn_err$warning) && is.null(warn_err$value$message))
  {
    count.df <- warn_err$value
    debug(logger, paste("input count data dimensions: ", dim(count.df)[1], " x ", dim(count.df)[2], sep=""))
    colNames <- names(count.df)
  } else {
    stop(paste(geterrmessage(), str(warn_err)))
  }
  count.df
}

###
### Switch/Case statement over vector
###
# BUG: not testing the case but perform all cases
# vswitch <- function(EXPR, ...) {
#   vars <- cbind(...)
#   vars[cbind(seq_along(EXPR), match(EXPR, names(list(...))))]
# }
# BUG: idem
# vswitch <- function(expr, ...) {
#   lookup <- list(...)
#   vec <- as.character(expr)
#   vec[is.na(vec)] <- "NA"
#   unname(do.call(c, lookup[vec]))
# }

###
### check BBRIC count data format
###
isCountDataBBRIC <- function(count.df)
{
  colNames <- names(count.df)
  if (length(colNames)>8) {
    debug(logger, "Count input file have more than 8 columns, count format seems to be BBRIC")
    # check for bbric format: check .count and .rpkm for 2 first columns pairs (a minimum to compare)
    if (length(colNames)>=12) {
      debug(logger, "Count input file have at least 12 columns, a minimum of 2 samples to compare")
      debug(logger, "Will check for BBRIC count format: check .count and .rpkm name extension for 2 first columns pairs")
      is_count <- grepl("\\.count$", colNames[c(9,11)])
      if (all(is_count)) { 
        debug(logger, "The 2 first columns pairs have .count columns") 
      } else
      {
        debug(logger, "Missing .count column in any of the 2 first columns pairs")
      }
      
      is_rpkm <- grepl("\\.rpkm$", colNames[c(10,12)])
      if (all(is_rpkm)) { 
        debug(logger, "The 2 first columns pairs have .rpkm columns") 
      } else
      {
        debug(logger, "Missing .rpkm column in any of the 2 first columns pairs")
      }
      
      if (all(is_count, is_rpkm)) {
        is_bbric_format = TRUE
        info(logger, "OK Count input file format is BBRIC")
        is_bbric_format
      } else
      {
        is_bbric_format = FALSE
        info(logger, "Count input file format is not BBRIC")
        is_bbric_format
      }
    } else
    {
      error(logger, "Count input file have less than 12 columns, need at least 2 samples to compare")
      stop("Count input file have less than 12 columns, need at least 2 samples to compare")
    }
  }
}

###
### check generic count data format
###
isCountDataGeneric <- function(count.df)
{
  colNames <- names(count.df)
  #debug(logger, "Count input file have less than 8 columns, count format seems to be generic")
  # check for generic format: at least 3 columns, 2 samples to compare
  if (length(colNames)>=3) {
    debug(logger, "Count input file have at least 3 columns, a minimum of 2 samples to compare")
    debug(logger, "Will check for generic count format: check that 2nd and 3rd columns values are numric vectors")
    #if (all(count.df[,2] == floor(count.df[,2])) && all(count.df[,3] == floor(count.df[,3]))) {
    if (all(is.numeric(count.df[,2])) && all(is.numeric(count.df[,3]))) {
      is_generic_format = TRUE
      debug(logger, "Generic format: all 2nd and 3rd columns values are numeric vectors")
      info(logger, "OK Count input file format is generic")
      is_generic_format
    } else {
      is_generic_format = FALSE
      debug(logger, "Not a generic format: all or any of the 2nd and 3rd columns values are not numeric vectors")
      info(logger, "Count input file format is not generic")
      is_generic_format
    }
    
  } else
  {
    error(logger, "Count input file have less than 3 columns, need at least 2 samples to compare")
    stop("Count input file have less than 3 columns, need at least 2 samples to compare")
  }
}
  
###
### Check count data format
###
isCountDataFormat <- function(count.df, format)
{
  formats=c("BBRIC", "generic")
  stopifnot(format %in% formats)
  
# BUG: not working, not testing the current case but performs all the cases  
#   vswitch(format,
#           BBRIC = isCountDataBBRIC(count.df),
#           isCountDataGeneric(count.df)
#           )
  
  if (format == "BBRIC")
  {
    isCountDataBBRIC(count.df)
    
  } else if (format == "generic")
  {
    isCountDataGeneric(count.df)
    
  }
  
}






