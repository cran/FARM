##################
# farm algorithm #
##################
#' @title  Forward Angular Relevance Metric (FARM)
#' @name   farm
#'
#' @author rch
#' @note   date: 11-17-2022


#' farm
#' @description         farm identifies and compensates time warping between query and reference time series by sample and hold for an insertion
#'                      and a deletion for a remove. Subsequent, the algorithm calculates two relevance values: (1) the correlation of the query
#'                      to the reference time series for each point locally in the given time window for local relevance and (2) the l2 norm of
#'                      the mean of a proportion with highest local relevance values and the mean of all local relevance values. In case, the
#'                      variation of any local sequence (i.e. local query or reference sequence) equals 0 (var(xts) == 0), it adds a small white
#'                      noise with a std. variation of 0.000001 to normalized (mean = 0) sequences. This avoids NA in using local correlation of
#'                      constant series.
#'
#' @param refTS         double vector: reference time series
#' @param qryTS         double vector: query time series
#' @param lcwin         int: length of local correlation window - default: 5
#' @param rel.th        int: proportion of local relevance measure with highest values included in global measure - default: 15
#' @param ff.align      boolean: forward feature alignment. If set to FALSE, the algorithm calculates relevance on originally passed time series without prior feature alignment. - default: TRUE
#' @param reshape.fnc   function(x): applied function on calculated local relevance values that are multiplied with qryTS for reshaping (suppressing irrelevant features). For threshold reshaping (i.e. shap.series = orig IF rel >= th ELSE 0) set the function to: function(x) (x >= th)*1 - default: function(x) x
#' @param fuzzyc        double vector: fuzzy coefficients. Ideally, length(fuzzyc) = lcwin. To disable fuzzification, set fuzzyc = 1 - default: c(2,6,10,6,2)/10
#' @param metric.space  boolean: triangle inequality compliance flag. If TRUE, the alignment algorithm only applies the sine function for distance instead of the enhanced sine-exp function composition. - default: TRUE
#' @param reshape.mode  keyword \{'uni', 'dual'\}: mode for decomposing. uni: decomposing only query time series. dual: decomposing query and reference time series. In the dual mode, decomposing only inserts data by interpolation but does not delete any information (no information loss). default: uni
#'
#' @return              list(path, rts.decomp, qts.decomp, qts.shaped, rel.local, rel.local.fuzz, rel.global, qmeas.mean.p, qmeas.slm.f)
#'                      path            warping path
#'                      rts.decomp      time warping decomposed reference time series
#'                      qts.decomp      time warping decomposed query time series
#'                      qts.shaped      shaped dewarped (if ff.align = T) query time series according to the reshape function reshape.fnc. The reshaping relies on the fuzzified rel.local.fuzz.
#'                      rel.local       local relevance values (length(rel.local) = length(qryTS)) compliant with metric requirements (triang. ineq.)
#'                      rel.local.end   local relevance values (not fuzzyfied) with correlation coefficient spot to the end of the lcwin.
#'                      rel.local.fuzz  local relevance values fuzzyfied according to fuzzyc factors applied on rel.local. Set fuzzyc = 1 to turn off fuzzification.
#'                      rel.global      global relevance value (scalar) based on local relevance values rel.local (being compliant with metric space). The global relevance relies on the fuzzified rel.local.fuzz.
#'                      qmeas.mean.p    mean absolute product quality measure for alignment of query with respect to reference time series. **higher** values indicate **better** alignment results.
#'                      qmeas.slm.f     slope-mean-factorial slm_f quality measure (mean product of slopes) for alignment of query with respect to reference time series. **higher** values indicate **better** alignment results.
#'
#' @export
#' @importFrom          stats cor embed rnorm
#' @importFrom          dplyr slice %>% n
#' @importFrom          utils head tail
#' @examples            farm(c(1,2,2,1,3), c(1,1,2,1,4), lcwin = 3, rel.th = 10)
farm <- function(refTS, qryTS, lcwin = 5, rel.th = 15, ff.align = TRUE, reshape.fnc = function(x) x, fuzzyc = c(2,6,10,6,2)/10, metric.space = TRUE, reshape.mode = 'uni')
{
  dref <- diff(refTS)                                       # delta reference ts
  dqry <- diff(qryTS)                                       # delta query ts
  lref <- length(dref)                                      # length delta reference values
  lqry <- length(dqry)                                      # length delta query values
  path <- data.frame(y=1, x=1)                              # de-warped path for given ts
  path.last <- path %>% slice(n())                          # index of last de-warping path element
  dist.root <- farm.dist(dref[1], dqry[1], metric.space)    # summed farm distance up to the current root element


  # ensure local correlation window size is smaller than the length of given time series
  stopifnot('The length of the local correlation window (lcwin) exceeds the length of one of the time series.' = (lcwin <= min(length(refTS), length(qryTS))))
  # prove that lcwin equals odd number
  stopifnot('The local correlation window (lcwin) is not an odd number.' = (lcwin%%2 == 1))
  # ensure rel.th is a valid proportion value 1 <= th <= 100
  stopifnot('The local relevance threshold is not a valid percentage value (1 <= rel.th <= 100).' = (rel.th > 1 & rel.th < 100))
  # check mode
  stopifnot("The reshape.mode must be either 'uni' or 'dual'." = (reshape.mode %in% c('uni','dual')))
  # fuzzy array must have an odd number of coefficients
  stopifnot("Fuzzy array 'fuzzyc' must have an odd number of coefficients." = (length(fuzzyc)%%2 == 1))
  # number of fuzzy coefficients must be < 2*min(len(refTS, qryTS))
  stopifnot("Number of fuzzy coefficients must be < 2 times minimal length of refTS or qryTS." = (length(fuzzyc) < 2*min(length(refTS), length(qryTS))))


  # init misalignment weight
  missal.y <- 0

  # time series feature aligning
  if(ff.align){
    while(path.last$x+1 <= lqry || path.last$y+1 <= lref){
      # calculate diagonal, dx+2y+1 and dx+1y+2 distance
      dist.xy   <- ifelse(path.last$y+1 <= lref && path.last$x+1 <= lqry, farm.dist(dref[path.last$y+1],dqry[path.last$x+1],metric.space), NA)
      dist.xyy  <- ifelse(path.last$y+2 <= lref && path.last$x+1 <= lqry, farm.dist(dref[path.last$y+2],dqry[path.last$x+1],metric.space), NA)
      dist.xxy  <- ifelse(path.last$y+1 <= lref && path.last$x+2 <= lqry, ifelse(path.last$y == path.last$x, Inf, farm.dist(dref[path.last$y+1],dqry[path.last$x+2],metric.space)), NA)
      # dist.xxyy <- ifelse(path.last$y+2 <= lref && path.last$x+2 <= lqry, farm.dist(dref[path.last$y+2],dqry[path.last$x+2],metric.space), NA)

      next.step <- switch(which.min(c(dist.xy, dist.xxy, dist.xyy, Inf)),
                          '1' = list(root = dist.xy   + dist.root, path = data.frame(y=path.last$y+1     , x=path.last$x+1)),
                          '2' = list(root = dist.xxy  + dist.root, path = data.frame(y=path.last$y+c(0,1), x=path.last$x+c(1,2))),
                          '3' = list(root = dist.xyy  + dist.root, path = data.frame(y=path.last$y+c(1,2), x=path.last$x+c(0,1))),
                          # '4' = list(root = dist.xxyy + dist.root, path = data.frame(y=path.last$y+c(1,2), x=path.last$x+c(1,2))),
                          # all distances are NA: reference boundary reached -> path must increment x or y only to fill missing values.
                          '4' = list(root =            dist.root, path = data.frame(y=ifelse(path.last$y < lref, path.last$y+1, path.last$y),
                                                                                    x=ifelse(path.last$x < lqry, path.last$x+1, path.last$x))))

      dist.root <- next.step$root
      path.last <- next.step$path

      # update path and set path.last to the very last entry
      path <- rbind(path, path.last)
      path.last <- path %>% slice(n())
    }

    # decompose original time series
    ts.dec <- tw.decomp(refTS, qryTS, path, reshape.mode)
  } else {
    # no alignment: keep time series unchanged
    ts.dec <- list(rts = refTS, qts = qryTS)
  }

  # calculate quality meas
  qmeas.mean.p = farm.mean.p(ts.dec$rts, ts.dec$qts)
  qmeas.slm.f = farm.slm.f(ts.dec$rts, ts.dec$qts, l.orig = length(refTS), l.adj = length(ts.dec$qts))

  # local relevance: correlation and fuzzyfication
  {
    # scale data
    rel.rts <- scale(ts.dec$rts)
    rel.qts <- scale(ts.dec$qts)

    # add zero-padding at the tails to force equal series length for results.
    # lcwin-1 zero-padding at start: enables aligning correlation coeff. to the spot being at the end of the window.
    # CAUTION: this causes irregular boundary correlation effects
    rel.rts <- c(rep(0,lcwin-1), rel.rts, rep(0,floor(lcwin/2)))
    rel.qts <- c(rep(0,lcwin-1), rel.qts, rep(0,floor(lcwin/2)))

    # calculate farm.dist for all element pairings (used for last spot correlation weighting)
    w.rq.farm.dist <- mapply(farm.dist, rel.rts, rel.qts, metric.space = TRUE)

    # prepare ref and qry series for mapply by rolling series with length = lcwin
    # CAUTION: the correlation of constant series causes warnings and return NA
    rel.rts <- embed(rel.rts,lcwin)[,lcwin:1]
    rel.qts <- embed(rel.qts,lcwin)[,lcwin:1]
    w.rq.farm.dist <- embed(w.rq.farm.dist, lcwin)[,lcwin:1]


    # CORRELATION between rolling windows of query and reference time series results a series rel.local ϵ [0,1]
    # prove of triangle inequality -> see Van Dongen et al. https://arxiv.org/abs/1208.3145
    # replace NAs with 0: caused by constant values with standard deviation (sd) == 0
    rel.local <- suppressWarnings(mapply(function(x,y) 1-sqrt(0.5*(1-cor(x,y))),
                                         as.data.frame(t(rel.rts)), as.data.frame(t(rel.qts))))
    rel.local[is.na(rel.local)] <- 0

    # create exponential weight function
    sum.weights <- seq(0,1,length.out=ncol(rel.rts))^6
    sum.weights <- sum.weights/sum(sum.weights)

    # calculate standard deviation and normalize data
    # rel.rts.norm <- rel.rts * apply(rel.rts,1,function(x) sd(x)*(length(x)-1)/length(x))
    # rel.qts.norm <- rel.qts * apply(rel.qts,1,function(x) sd(x)*(length(x)-1)/length(x))

    # get sum of weighted farm similarity for all windows
    # matrix((1-mapply(farm.dist, abs(rel.rts.norm), abs(rel.qts.norm)))*sum.weights, nrow = nrow(rel.qts.norm))
    w.rq.farm.dist <- rowSums((1-w.rq.farm.dist)*(rep(1,nrow(w.rq.farm.dist)) %*% t(sum.weights)))

    # weight local relevance for last spot and
    # cut tail by floor(win/2) to set spot to the end of sliding window
    rel.local.end <- head(rel.local * w.rq.farm.dist, -floor(lcwin/2))
    # cut head by floor(win/2) to set spot centric of sliding window
    rel.local <- tail(rel.local, -floor(lcwin/2))


    # FUZZIFICATION - include the relevance of adjacent data points to the final relevance measure
    # prepare local relevance series
    fuzzy.l <- length(fuzzyc)
    # append l(filter)/2 times 0 (i.e. NAs) at start and end for avoiding tailing effects
    rel.local.fuzz <- c(rep(0,floor(fuzzy.l/2)), rel.local, rep(0,floor(fuzzy.l/2)))

    # FUZZIFY local relevance coefficients
    # create rolling matrix (ncol = fuzzy.l, nrow = length(fuzzyc)), multiply each element with fuzzyc and sum all
    rel.local.fuzz <- embed(rel.local.fuzz, fuzzy.l)[,fuzzy.l:1] %*% as.matrix(fuzzyc)
    # normalize output to 1
    rel.local.fuzz <- rel.local.fuzz/sum(fuzzyc)
  }

  # NOTE: to continue without fuzzification set fuzzyc = 1
  # global relevance
  minlen <- min(length(ts.dec$rts), length(ts.dec$qts))
  rel.global <- sqrt((mean(sort(rel.local.fuzz, T)[1:(rel.th/100*length(rel.local.fuzz))])^2 + mean(rel.local.fuzz)^2) /2)

  # reshape of query time series
  qts.shaped <- ts.dec$qts * sapply(rel.local.fuzz, reshape.fnc)

  return(list(path = path,
              rts.decomp = ts.dec$rts,
              qts.decomp = ts.dec$qts,
              qts.shaped = qts.shaped,
              rel.local = rel.local,
              rel.local.end = rel.local.end,
              rel.local.fuzz = rel.local.fuzz,
              rel.global = rel.global,
              qmeas.mean.p = qmeas.mean.p,
              qmeas.slm.f = qmeas.slm.f))
}





