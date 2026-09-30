#######################################
# Decomposition of warped time series #
#######################################
#' @title FARM time warping decomposition
#' @name tw.decomp
#'
#' @author rch
#' @note   date: 11-17-2022


#' tw.decomp
#' @description Two time series present a warping scheme to each other. For comparison and advanced analysis
#' of warped time series, it requires de-warping time series according a given scheme. This function expands
#' the reference and query time series for resulting in a 1 to 1 data point alignment. The series are expanded
#' according to the farm warping solution by interpolation. There are no data points removed.
#'
#' @param rts   vector: original reference time series
#' @param qts   vector: original query time series
#' @param twp   data.frame: warped data point assignment. df requires an x (qry) and y (ref) variable providing assignment information.
#' e.g. a data point of query ts is assigned to two subsequent reference data points results to: x=\[.., x_i, x_i+1, x_i+1, ..\] y=\[.., y_i, y_i+1, y_i+2, ..\]
#' @param mode  keyword \{'uni', 'dual'\}: mode for decomposing. uni: decomposing only query time series. dual: decomposing query and reference time series. In the dual mode, decomposing only inserts data by interpolation but does not delete any information (no information loss). default: uni
#' @param abs   bool: true if series are amplitude values (absolute); false if series are delta values of subsequent samples (relative)
#'
#' @return     list(rts=c(rts1, rts2, ..), qts=c(qts1, qts2, ..)): de-warped reference (rts) and query (qts) time series.
#' @export
#' @examples   tw.decomp(c(0,1,0,0,0), c(0,0.5,0.5,0,0), data.frame(x=c(1,2,2,3,4), y=c(1,2,3,4,4)))
tw.decomp <- function(rts, qts, twp, mode = 'uni', abs = FALSE){
  stopifnot("mode must be either 'uni' or 'dual'" = (mode %in% c('uni','dual')))
  stopifnot('if the largest indices equals the number of data points in series, only dewarping on absolute values possible.' = ((abs == FALSE) == (max(length(rts), length(qts)) > max(twp$y, twp$x))))
  stopifnot('length of series must be equal for abs=TRUE or +1 the largest index of dewarping indices.' = ((length(rts) >= max(twp$y)) && (length(qts) >= max(twp$x))))

  if(abs){ # warping relying on data points like in dtw
    if(mode == 'uni'){ # adjust query ts only
      dec.qts <- qts[twp$x]   # inserts samples in query ts
      # deletions in query ts
      dec.qts <- sapply(unique(twp$y), function(x) mean(dec.qts[twp$y==x]))
      dec.rts <- rts          # keep reference ts unchanged
    } else if(mode == 'dual') { # adjust query and reference time series
      dec.qts <- qts[twp$x]   # apply query warping
      dec.rts <- rts[twp$y]   # apply reference warping
    }
  } else { # warping relying on differences between data points
    # insert samples in query ts (linear interpolation)
    dec.qts <- c(qts[1],unlist(sapply(unique(twp$x), function(x) seq(qts[x],qts[x+1],length.out = (length(twp$x[twp$x == x])+1))[-1])))
    if(mode == 'uni'){ # adjust query ts only
      dec.qts <- dec.qts[!c(1,diff(twp$y))==0] # delete samples in query ts instead inserts in reference ts
      dec.rts <- rts                           # keep reference ts untouched
    } else if(mode == 'dual'){ # adjust query and reference time series
      # insert samples in reference ts (linear interpolation)
      dec.rts <- c(rts[1],unlist(sapply(unique(twp$y), function(x) seq(rts[x],rts[x+1],length.out = (length(twp$y[twp$y == x])+1))[-1])))
    }
  }

  return(list(rts = dec.rts, qts = dec.qts))
}
