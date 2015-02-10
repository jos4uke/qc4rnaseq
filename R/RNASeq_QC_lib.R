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
## @license: GPL (>= 2)
##
###############################################

version <- "0.0.1.0"

copyright <- "GPL (>= 2)"

###
### FUNCTIONS ###
###

###
### check for installed package
###
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])

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
  warn_err <- tryCatch.W.E(read.table(file, header=TRUE, sep="\t", check.names=FALSE))
  
  if (is.null(warn_err$warning) && is.null(warn_err$value$message))
  {
    count.df <- warn_err$value
    #debug(logger, paste("input count data dimensions: ", dim(count.df)[1], " x ", dim(count.df)[2], sep=""))
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
  # check for bbric format
  if (length(colNames)>=12) {
    #debug(logger, "Count input file have at least 12 columns, a minimum of 2 samples to compare")
    #debug(logger, "Will check for BBRIC count format: check -count and -rpkm name extension for 2 first columns pairs")
    # check for bbric format: check -count and -rpkm for odd and even columns pairs 
    is_count <- grepl("\\-count$", colNames[seq(9, length(colNames), by=2)])
#     if (all(is_count)) { 
#       debug(logger, "The odd columns pairs have -count columns") 
#     } else
#     {
#       debug(logger, "Missing -count column in any of the odd columns pairs")
#     }    
    is_rpkm <- grepl("\\-rpkm$", colNames[seq(10, length(colNames), by=2)])
#     if (all(is_rpkm)) { 
#       debug(logger, "The even columns pairs have -rpkm columns") 
#     } else
#     {
#       debug(logger, "Missing -rpkm column in any of the even columns pairs")
#     }
    
    # check for bbric format: check -count and -rpkm columns have the same lib name
    is_idem <- gsub("-count", "", colNames[seq(9, length(colNames), by=2)]) == gsub("-rpkm", "", colNames[seq(10, length(colNames), by=2)])
#     if (all(is_idem)) { 
#       debug(logger, "The -count and -rpkm columns have the same lib name") 
#     } else
#     {
#       debug(logger, "All the -count and -rpkm columns don't have the same lib name")
#     }	  
    
    # check for bbric format: check -count and -rpkm columns contain only numeric values  
    is_num <- c()
    for (i in 9:length(colNames)) { is_num <- append(is.numeric(count.df[,i]), is_num)}
#     if (all(is_num)) { 
#       debug(logger, "-count and -rpkm columns contain only numeric values") 
#     } else
#     {
#       debug(logger, "Values not numeric are found in the -count and -rpkm columns")
#     }
    
    if (all(is_count, is_rpkm, is_idem, is_num)) {
      is_bbric_format = TRUE
#       info(logger, "OK Count input file format is BBRIC")
      is_bbric_format
    } else
    {
      is_bbric_format = FALSE
#       info(logger, "Count input file format is not BBRIC")
      is_bbric_format
    }
  } else
  {
    #print("Count input file have less than 12 columns, need at least 2 samples to compare", file=stderr())
    #stop("Count input file have less than 12 columns, need at least 2 samples to compare")
    is_bbric_format = FALSE
    is_bbric_format
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

    #if (all(is.numeric(count.df[,2])) && all(is.numeric(count.df[,3]))) {
	is_num <- c()
	for (i in 2:length(colNames)) { is_num <- append(is.numeric(count.df[,i]), is_num)}
	if (all(is_num)) {
		is_generic_format = TRUE
		debug(logger, "Generic format: 2nd and all following columns values are numeric vectors")
		info(logger, "OK Count input file format is generic")
		is_generic_format
    } else {
      is_generic_format = FALSE
      debug(logger, "Not a generic format: all or any of the 2nd and following columns values are not numeric vectors")
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


###
### Load Stats data
###
loadStatsData <- function(file)
{
  # load Stats data
  warn_err <- tryCatch.W.E(read.table(opt$stats, header=TRUE, sep="\t", check.names=FALSE))
  
  if (is.null(warn_err$warning) && is.null(warn_err$value$message))
  {
    stats.df <- warn_err$value
    debug(logger, paste("input stats data dimensions: ", dim(stats.df)[1], " x ", dim(stats.df)[2], sep=""))
    StatscolNames <- names(stats.df)
  } else {
    stop(paste(geterrmessage(), str(warn_err)))
  }
  stats.df
}

###
### Check Stats data format
###

isStatsDataFormat <- function(stats.df)
{
  StatscolNames <- names(stats.df)
  if (length(StatscolNames) == 7) {
    debug(logger, "Stats input file has 7 columns as expected")
    # check for stats format: check the StatscolNames
	TrueStatscolNames <- c("lib", "specific_hits", "mapping_hits",	"raw_reads/pairs_count", "mapped_reads/pairs_count", "feature_overlapping_hits", "feature_overlapping_hits/specific_hits_percent" )
	if (identical(StatscolNames, TrueStatscolNames)) {
		debug(logger, "Stats input file has the expected column names")
		info(logger, "OK Stats input file format is OK")
		is_stats_format = TRUE
		is_stats_format
		} else
		{
		debug(logger, "Stats input file doesnt' have the expected column names")
		info(logger, "Stats input file doesnt' have the expected column names")
		is_stats_format = FALSE
		is_stats_format
		}
	} else {
		error(logger, "Stats input file doesn't have 7 columns as expected")
		stop("Stats input file doesn't have 7 columns as expected")
	}
}
	

###
### Load Design data
###	
loadDesignData <- function(file)
{
  # load Design data
  warn_err <- tryCatch.W.E(read.table(opt$design, header=TRUE, sep="\t", check.names=FALSE))
  
  if (is.null(warn_err$warning) && is.null(warn_err$value$message))
  {
    design.df <- warn_err$value
    debug(logger, paste("input design data dimensions: ", dim(design.df)[1], " x ", dim(design.df)[2], sep=""))
  } else {
    stop(paste(geterrmessage(), str(warn_err)))
  }
  design.df
}	
	
	
###
### Check Design data format
###

isDesignDataFormat <- function(design.df)
{
  DesigncolNames <- names(design.df)
  if (length(DesigncolNames) == 2) {
    debug(logger, "Design input file has 2 columns as expected")
	is_design_format = TRUE
	is_design_format
  } else
  {
	is_design_format = FALSE
	is_design_format
	error(logger, "Design input file doesn't have 2 columns as expected")
	stop("Design input file doesn't have 2 columns as expected")
   }
}


##
## Check count data and design data consistency (same lib names)
##

isCountDesign <- function(count.df, design.df, format)
{
	colNames <- names(count.df)
	# Take the lib_names from count data
	if (format == "BBRIC") {
		countlibNames <- gsub("-count", "", colNames[seq(9, length(colNames), by=2)])
	} else if (format == "generic") {
		countlibNames <- colNames[2:length(colNames)]
	}
	# Take the lib_names from design data
	designlibNames <- design.df[,1]
	# Compare both count and design lib names
	if (all(countlibNames == designlibNames)){
		debug(logger, "Count and design data have the same lib names")
		info(logger, "OK Count and design data have the same lib names")
		is_count_design = TRUE
		is_count_design
		} else
		{
		debug(logger, "Count and design data don't have the same lib names")
		info(logger, "Count and design data don't have the same lib names")
		is_count_design = FALSE
		is_count_design
		}
}

##
## Check stats data and design data consistency (same lib names)
##
	
isStatsDesign <- function(stats.df, design.df)
{
	# Take the lib_names from stats data
	statslibNames <- stats.df[,1]
	# Take the lib_names from design data
	designlibNames <- design.df[,1]
	# Compare both stats and design lib names
	if (identical(statslibNames, designlibNames)) {
		debug(logger, "Stats and design data have the same lib names")
		info(logger, "OK Stats and design data have the same lib names")
		is_stats_design = TRUE
		is_stats_design
		} else
		{
		debug(logger, "Stats and design data don't have the same lib names")
		info(logger, "Stats and design data don't have the same lib names")
		is_stats_design = FALSE
		is_stats_design
		}	
}


