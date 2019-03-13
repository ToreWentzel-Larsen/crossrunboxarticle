library(Rmpfr)
library(crossrun)

# Set parameters ----
nmin     <- 10
nmax     <- 100
shifts   <- seq(0, 3, by = 0.2)
prec.use <- 120
one      <- mpfr(1, prec.use)
two      <- mpfr(2, prec.use)
mone     <- mpfr(-1, prec.use)
shiftsc  <- format(shifts, nsmall = 1)
nshifts  <- length(shifts)

# bestbox function ----
## Function for box with lowest probability for the target shift, among boxes
## with probability >= target for shift 0.
bestbox <- function(pt0    = crs$pt_0.0,
                    pts    = crs$pt_0.8,
                    target = 0.925,
                    n1     = 100,
                    mult   = 2,
                    prec   = 120) {
  nill    <- mpfr(0, prec)
  one     <- mpfr(1, prec)
  two     <- mpfr(2, prec)
  multm   <- mpfr(mult, prec)
  targetm <- mpfr(target, prec)
  targt   <- targetm * (multm ^ (n1 - 1)) # target on "times" scale
  pt0n    <- pt0[[n1]]
  ptsn    <- pts[[n1]]
  bpt0    <- boxprobt(pt0n) # box probabilities for no shift
  bpttarg <- boxprobt(ptsn) # box probabilities for target shift
  boxprt  <-
    two * (multm ^ (n1 - 1)) # initialize to impossible high value
  
  for (cc in 0:(n1 - 1))
    for (ll in 1:n1) {
      if (pt0n[cc + 1, ll] > nill &
          bpt0[cc + 1, ll] >= targt &
          bpttarg[cc + 1, ll] < boxprt) {
        c1 <- cc
        l1 <- ll
        boxprt <- bpttarg[cc + 1, ll]
      }
    } # end search through (c, l) with positive no (, l)
  return(c(c1, l1))
} # end function bestbox

# cutbox function ----
## Function for cutting a box while keeping probability >= target for shift 0. No
## cutting if the corner cannot be removed. If the corner may be removed, it is
## attempted to remove parts of the border, starting from the corner, in the
## direction with highest point probability for the target shift:
cutbox <- function(pt0    = crs$pt_0.0,
                   pts    = crs$pt_0.8,
                   target = 0.925,
                   n1     = 100,
                   c1     = 41,
                   l1     = 10,
                   mult   = 2,
                   prec   = 120) {
  nill      <- mpfr(0, prec)
  one       <- mpfr(1, prec)
  two       <- mpfr(2, prec)
  multm     <- mpfr(mult, prec)
  targetm   <- mpfr(target, prec)
  targt     <- targetm * (multm ^ (n1 - 1)) # target on "times" scale
  pt0n      <- pt0[[n1]]
  ptsn      <- pts[[n1]]
  bpt0      <-
    boxprobt(pt0n) # box probabilities for no shift, pt scale
  boxpt0    <-
    bpt0[c1 + 1, l1] # no shift probability of actual box, pt scale
  cornerpt0 <-
    pt0n[c1 + 1, l1] # no shift corner probability, pt scale
  finished  <- FALSE
  cbord     <- NA
  lbord     <- NA
  
  if (boxpt0 - cornerpt0 >= targt) {
    cutboxpt0 <-
      boxpt0 - cornerpt0 # pt of cutted box after removed corner
    cbord     <- c1 + 1
    lbord     <- l1 - 1
    
    while (!finished) {
      pt0n.directionc <- pt0n[cbord + 1, l1]
      pt0n.directionl <- pt0n[c1 + 1, lbord]
      ptsn.directionc <- ptsn[cbord + 1, l1]
      ptsn.directionl <- ptsn[c1 + 1, lbord]
      if ((cutboxpt0 - pt0n.directionc < targt |
           pt0n.directionc == 0) &
          (cutboxpt0 - pt0n.directionl < targt |
           pt0n.directionl == 0)) {
        finished <- TRUE
      } else if (cutboxpt0 - pt0n.directionc < targt |
                 pt0n.directionc == 0) {
        lstrip    <- pt0n[c1 + 1, lbord:1]
        nlstrip   <- length(lstrip)
        maxlstrip <- max((1:nlstrip)[lstrip > 0])
        lstrip    <- lstrip[(1:nlstrip) <= maxlstrip]
        lstripcum <- cumsum(lstrip)
        if (cutboxpt0 - max(lstripcum) >= targt) {
          lbord <- 0
        } else {
          # 0 cannot occurr
          lbord <-
            lbord + 1 - min((1:nlstrip)[cutboxpt0 - lstripcum < targt])
        }
        finished <- TRUE
      } else if (cutboxpt0 - pt0n.directionl < targt |
                 pt0n.directionl == 0) {
        cstrip    <- pt0n[(cbord + 1):n1, l1]
        ncstrip   <- length(cstrip)
        maxcstrip <- max((1:ncstrip)[cstrip > 0])
        cstrip    <- cstrip[(1:ncstrip) <= maxcstrip]
        cstripcum <- cumsum(cstrip)
        if (cutboxpt0 - max(cstripcum) >= targt) {
          cbord <- n1
        } else {
          # n1 cannot occurr
          cbord <-
            cbord + min((1:ncstrip)[cutboxpt0 - cstripcum < targt]) - 1
        }
        finished <- TRUE
      } else if (ptsn.directionc >= ptsn.directionl) {
        cbord     <- cbord + 1
        cutboxpt0 <- cutboxpt0 - pt0n.directionc
      } else if (ptsn.directionc < ptsn.directionl) {
        lbord     <- lbord - 1
        cutboxpt0 <- cutboxpt0 - pt0n.directionl
      }
    } # end while loop
  } # end if corner may be removed
  return(c(cbord, lbord))
} # end function cutbox

system.time({
  
  # Compute joint distributions of C and L for n = 1, ..., nmax ----
  crs <- list()
  for (s in shifts) {
    r <- paste0('pt_', format(s, nsmall = 1))
    print(paste('Joint distribution:', r))
    crs[[r]] <- crossrunshift(nmax, s)$pt
  }
  
  # Table 1 in "Run charts revisited", PLOS ONE November 25, 2014 ----
  bounds <- data.frame(
    n  = nmin:nmax,
    ca = qbinom(0.05, nmin:nmax - 1, 0.5),
    la = round(log2(nmin:nmax) + 3)
  )
  row.names(bounds) <- bounds$n
  
  # find best boxes
  bounds$cb <- NA
  bounds$lb <- NA
  
  for (nn in nmin:nmax) {
    print(paste('bestbox:', nn))
    bounds[bounds$n == nn, c('cb', 'lb')] <- bestbox(n1 = nn)
  }
  
  # find cut  boxes
  bounds$cbord <- NA
  bounds$lbord <- NA
  
  for (nn in nmin:nmax) {
    print(paste('cutbox', nn))
    bounds[bounds$n == nn, c('cbord', 'lbord')] <-
      cutbox(n1 = nn,
             c1 = bounds$cb[bounds$n == nn],
             l1 = bounds$lb[bounds$n == nn])
  }
  
  # Find no signal probabilities ----
  ## initialize to impossible (negative) value:
  pat <- mpfr2array(rep(mone, (nmax - 9) * nshifts),
                    dim = c((nmax - 9), nshifts),
                    dimnames = list(nmin:nmax, NULL))
  
  pbt       <- pat
  pct       <- pat
  loglrposa <- pat[ , -1]
  loglrposb <- loglrposa
  loglrposc <- loglrposa
  loglrnega <- loglrposa
  loglrnegb <- loglrposa
  loglrnegc <- loglrposa
  
  colnames(pat)       <- paste0('pat_', shiftsc)
  colnames(pbt)       <- paste0('pbt_', shiftsc)
  colnames(pct)       <- paste0('pct_', shiftsc)
  colnames(loglrposa) <- paste0('loglrposa_', shiftsc[-1])
  colnames(loglrposb) <- paste0('loglrposb_', shiftsc[-1])
  colnames(loglrposc) <- paste0('loglrposc_', shiftsc[-1])
  colnames(loglrnega) <- paste0('loglrnega_', shiftsc[-1])
  colnames(loglrnegb) <- paste0('loglrnegb_', shiftsc[-1])
  colnames(loglrnegc) <- paste0('loglrnegc_', shiftsc[-1])
  
  ## calculations
  for (nn in nmin:nmax) {
    ca1    <- bounds$ca[bounds$n == nn]
    la1    <- bounds$la[bounds$n == nn]
    cb1    <- bounds$cb[bounds$n == nn]
    lb1    <- bounds$lb[bounds$n == nn]
    cbord1 <- bounds$cbord[bounds$n == nn]
    lbord1 <- bounds$lbord[bounds$n == nn]
    
    for (s in shifts) {
      i              <- match(s, shifts)
      p              <- format(s, nsmall = 1)
      p              <- paste0('pt_', p)
      pat[nn - 9, i] <- sum(crs[[p]][[nn]][(ca1 + 1):nn, 1:la1])
      pbt[nn - 9, i] <- sum(crs[[p]][[nn]][(cb1 + 1):nn, 1:lb1])
      pct[nn - 9, i] <- pbt[nn - 9, i]
      
      if (!is.na(cbord1)) {
        pct[nn - 9, i] <- 
          sum(crs[[p]][[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
          sum(crs[[p]][[nn]][(cbord1 + 1):nn, lb1]) +
          sum(crs[[p]][[nn]][cb1 + 1, 1:lbord1])
      }
    }
  }
  
  # Find likelihood ratios ----
  for (s in shiftsc[shifts > 0]) {
    pats <- paste0('pat_', s)
    pbts <- paste0('pbt_', s)
    pcts <- paste0('pct_', s)
    loglrposa1 <- paste0('loglrposa_', s)
    loglrposb1 <- paste0('loglrposb_', s)
    loglrposc1 <- paste0('loglrposc_', s)
    loglrnega1 <- paste0('loglrnega_', s)
    loglrnegb1 <- paste0('loglrnegb_', s)
    loglrnegc1 <- paste0('loglrnegc_', s)
    
    loglrposa[, loglrposa1] <-
      log(two ^ (nmin:nmax - 1) - pat[, pats]) - 
      log(two ^ (nmin:nmax - 1) - pat[, 'pat_0.0'])
    loglrposb[, loglrposb1] <-
      log(two ^ (nmin:nmax - 1) - pbt[, pbts]) - 
      log(two ^ (nmin:nmax - 1) - pbt[, 'pbt_0.0'])
    loglrposc[, loglrposc1] <-
      log(two ^ (nmin:nmax - 1) - pct[, pcts]) - 
      log(two ^ (nmin:nmax - 1) - pct[, 'pct_0.0'])
    
    loglrnega[, loglrnega1] <-
      log(pat[, pats]) - log(pat[, 'pat_0.0'])
    loglrnegb[, loglrnegb1] <-
      log(pbt[, pbts]) - log(pbt[, 'pbt_0.0'])
    loglrnegc[, loglrnegc1] <-
      log(pct[, pcts]) - log(pct[, 'pct_0.0'])
  }
  
  # Create bounds table including probability information ----
  pa <- pat / (two ^ (nmin:nmax - 1))
  pb <- pbt / (two ^ (nmin:nmax - 1))
  pc <- pct / (two ^ (nmin:nmax - 1))
  
  boundspll <- cbind(
    bounds,
    asNumeric(pa),        asNumeric(pb),        asNumeric(pc),
    asNumeric(loglrposa), asNumeric(loglrposb), asNumeric(loglrposc),
    asNumeric(loglrnega), asNumeric(loglrnegb), asNumeric(loglrnegc)
  )
  
  # Fix column names, pat -> pa etc.
  names(boundspll) <- sub("(^p.{1}).", "\\1\\", names(boundspll))
  
  # # Save objects ----
  # ## save boundspll with Rmpfr background arrays:
  # save(
  #   boundspll,
  #   pat, pbt, pct,
  #   pa, pb, pc,
  #   loglrposa, loglrposb, loglrposc,
  #   loglrnega, loglrnegb, loglrnegc,
  #   file = 'data/boundspll.Rdata'
  # )
  # 
  # ## save crs Rmpfr arrays
  # saveRDS(crs, 'data/crs_mpfr.rds')
  # 
  # ## save crs numeric arrays
  # crs_num <- lapply(crs, function(x) {
  #   lapply(x, asNumeric)
  # })
  # saveRDS(crs_num, 'data/crs_num.rds')
  # 
  # ## save cr distribution for shift = 0
  # saveRDS(cr_num$pt_0.0, 'data/cr_dist.rds')
  # 
  # ## save box limits and probabilities
  # saveRDS(boundspll, 'data/cr_bounds.rds')
})