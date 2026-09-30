################################
# Farm distance of two vectors #
################################
#' @title FARM distance
#' @name farm.dist
#'
#' @author rch
#' @note   date: 11-17-2022


#' farm.dist
#' @description The returned distance is the sine of the absolute angle between vectors if
#' both have a positive or both have a negative argument. For contrasting vectors
#' (v1: +-x,+y and v2: +-x,-y) the resulting distance is: 1-exp(-phi*5), for non-metric space
#' and: . The vector argument 0
#' is treated as positive argument.
#'
#' @param dyr           double: dy reference normalized to dx = 1
#' @param dyq           double: dy query normalized to dx = 1. Not required if dyr also comprises query delta values in a data frame.
#' @param metric.space  boolean: triangle inequality compliance flag. If TRUE, the distance relies on the sine function only instead of the enhanced sine-exp function composition. - default: FALSE
#'
#' @return      double: distance value
#'
#' @export
#' @examples    farm.dist(0.23, -2.85)
farm.dist <- function(dyr, dyq=NA, metric.space = FALSE){
  vrc <- complex(0,1,0)
  v1c <- ifelse(is.na(dyq), complex(0,1,dyr[[1]]), complex(0,1,dyr))
  v2c <- ifelse(is.na(dyq), complex(0,1,dyr[[2]]), complex(0,1,dyq))

  phi.v1 <- Arg(v1c)
  phi.v2 <- Arg(v2c)

  # zero difference between data points implies no information
  # thus, we avoid an alignment to zero information
  if(phi.v1*phi.v2 == 0) return(Inf)

  # ignore orientation, only distance is relevant; thus we take the abs()
  phi.d <- abs(phi.v1 - phi.v2)

  # check if metric space required
  if(metric.space){
    dist <- ifelse(phi.d > pi/2, 1, sin(phi.d)/(1+sin(max(phi.v1, phi.v2))))           #v1
    # dist <- ifelse(phi.d > pi/2, 1, sin(phi.d)/(1+sin(max(abs(c(phi.v1, phi.v2))))))   #v2
    # dist <- sin(phi) if 0 <= phi <= pi/2 else 1                                        #v3
  } else {
    # set align neg if any is negative and none == 0 (0=0°)
    align <- sign(phi.v1) * sign(phi.v2)
    # set align positive if any phi has 0°
    if(align == 0) align <- 1
    # calculate dist based on the combined 1-exp and sine approach
    dist <- ifelse(align >=0, sin(phi.d), max(sin(phi.d), (1-exp(-phi.d*5))))
  }

  return(dist)
}


#' farm.mean.p
#' @description An absolute values based quality measure for feature alignment of two independent time series with the **same** length.
#' The measure returns the mean of the absolute product of two time series according the formula:
#' \deqn{\left |{1 \over n} \sum_{i=0}^{n} {ref_i \cdot qry_i} \right|}
#' Resulting from the product, the **higher** the values, the **better** the result.
#'
#' @param rts  double vector: reference time series of length n
#' @param qts  double vector: query time series of length n
#'
#' @return     double: quality measure *larger means better*
#'
#' @export
#' @examples   farm.mean.p(c(1,2,1.5,2.1,5,2.32,1,0.2), c(0.2,0.4,0.4,0.5,1,0.92,0.5,0.2))
farm.mean.p = function(rts, qts){
  if(length(qts) != length(rts)) warning('reference and query time sereis must be of same length!')
  else
    return(1/length(rts) * abs(sum(as.double(rts) * as.double(qts), na.rm = T)))
}


#' farm.slm.f
#' @description The slope-mean-factorial slm_f returns a quality measure for feature alignment of two independent time series with the **same** length.
#' The measure relies on the mean product of the slopes (i.e. differences) of subsequent data values in two time series, according the formula:
#' \deqn{{1 \over n} \sum_{i=2}^{n} {\Delta ref_{[i-1, i]} \cdot \Delta qry_{[i-1, i]}} \over {1 + { \left | \#smp_{orig} - \#smp_{method} \over \#smp_{orig} \right | }}}
#' Resulting from the inner product of the slopes, the **higher** the values, the **better** the result. For penalizing the number of changes in time series for feature matching,
#' the original and modified series length can be set by the parameters l.orig and l.adj respectively.
#'
#' @param rts     double vector: reference time series of length n
#' @param qts     double vector: query time series of length n
#' @param l.orig  int: length of original time series - default: 1
#' @param l.adj   int: length of adjusted time series - default: 1
#'
#' @return     double: quality measure *larger means better*
#'
#' @export
#' @examples   farm.slm.f(c(1,2,1.5,2.1,5,2.32,1,0.2), c(0.2,0.4,0.4,0.5,1,0.92,0.5,0.2))
farm.slm.f = function(rts, qts, l.orig = 1, l.adj = 1){
  if(length(qts) != length(rts)) warning('reference and query time sereis must be of same length!')
  else
    return(1/length(rts) * sum((diff(rts)*diff(qts)), na.rm = T) / (1 + abs((l.orig - l.adj)/l.orig)))
}
