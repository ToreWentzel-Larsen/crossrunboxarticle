# ###############################################################################
# This script reproduces the data objects and figures from the Smooth Operator
# article by Jacob Anhøj and Tore Wentzel-Larsen, R Journal 2019.
#
# The objects of interest are cr_dists and cr_bounds.
#
# cr_dists: a list with probabilities for the joint distribution of number of
# crossings (C) and longest run (L) in multiple precision (mpfr) format. To get
# the matrix of probabilities for, say, N = 11 and no shift (SD = 0), use
# cr_dists$pt_0.0[[11]].
#
# cr_bounds: a data frame with critical values for longest run and number of
# crossings together with probabilities and log-likelihood ratios for the Anhøj,
# best box and cut box rules.
# Variables:
#   ca = lower limit for number of crossings, Anhøj rules. 
#   la = upper limit for longest run, Anhøj rules. 
#   cb = lower limit for number of crossings, best box rules. 
#   lb = upper limit for longest run, best box rules. 
#   cbord/lbord = coordinates for the cut box adjustment.
#   pa_n.m/pb_n.m / pcbord_n.m = probality of no signal. n.m = size of the shift
#                              in standard deviation units.
#   loglrpos_n.m / loglrneg_n.m = positive and negative log-likehood ratios.
# 
# For the sake of speed and memory consumption the scrip is by default set to
# produce output for N = 10-40 and SD = 0-2. To reproduce all data from the
# article, change the parameters nmax and smax to 100 and 3 respectively.
# 
# Jacob Anhøj & Tore Wentzel-Larsen 22 Mar 2019
################################################################################

# Load libraries ----
library(Rmpfr)
library(crossrun)
library(tidyverse)

# Set parameters ----
nmax         <- 40      # Max N to include in computations.
smax         <- 2       # Max shift in SD units to include in computations.
target       <- 0.925   # Target specificity for best box and cut box.
target_shift <- 0.8     # Target shift for best box and cut box.

## Probably no need to change anything below this line
nmin     <- 10
shifts   <- seq(0, smax, by = 0.2)

# bestbox function ----
## Function for box with lowest probability for the target shift, among boxes
## with probability >= target for shift = 0.
bestbox <- function(pt0,
                    pts,
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
  bpt0    <- boxprobt(pt0n)  # box probabilities for no shift
  bpttarg <- boxprobt(ptsn)  # box probabilities for target shift
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
    }
  return(c(c1, l1))
} # end function bestbox

# cutbox function ----
## Function for cutting a box while keeping probability >= target for shift = 0.
## No cutting if the corner cannot be removed.
cutbox <- function(pt0,
                   pts,
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
    boxprobt(pt0n)   # box probabilities for no shift, pt scale
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
    }
  }
  return(c(cbord, lbord))
} # end function cutbox

# crs function to compute joint distributions of C and L ----
crs <- function(nmax = 12,
                shifts = seq(0, 2, by = 0.2)) {
  
  crs <- list()
  for (s in shifts) {
    r <- paste0('pt_', format(s, nsmall = 1))
    print(paste('Joint distribution:', r))
    crs[[r]] <- crossrunshift(nmax, s)$pt
  }
  return(crs)
} # End crs function

# bounds function to compute limits and diagnostics for runs rules ----
bounds <- function(crs, target = 0.925, target_shift = 0.8) {
  nmin     <- 10
  nmax     <- length(crs[[1]])
  shiftsc  <- format(shifts, nsmall = 1)
  nshifts  <- length(shifts)
  prec.use <- 120
  one      <- mpfr(1, prec.use)
  two      <- mpfr(2, prec.use)
  mone     <- mpfr(-1, prec.use)
  pt0      <- crs$pt_0.0
  # pts      <- crs$pt_0.8
  pts      <- crs[[paste0('pt_', format(target_shift, nsmall = 1))]]
  
  # Begin bounds table
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
    bounds[bounds$n == nn, c('cb', 'lb')] <- bestbox(pt0,
                                                     pts,
                                                     n1 = nn,
                                                     target = target)
  }
  
  # find cut  boxes
  bounds$cbord <- NA
  bounds$lbord <- NA
  
  for (nn in nmin:nmax) {
    print(paste('cutbox', nn))
    bounds[bounds$n == nn, c('cbord', 'lbord')] <-
      cutbox(pt0,
             pts,
             n1 = nn,
             target = target,
             c1 = bounds$cb[bounds$n == nn],
             l1 = bounds$lb[bounds$n == nn])
  }
  
  # Find no signal probabilities
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
  
  # Find likelihood ratios
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
  
  # Finish bounds table including probability information
  pa <- pat / (two ^ (nmin:nmax - 1))
  pb <- pbt / (two ^ (nmin:nmax - 1))
  pc <- pct / (two ^ (nmin:nmax - 1))
  
  bounds <- cbind(
    bounds,
    asNumeric(pa),        asNumeric(pb),        asNumeric(pc),
    asNumeric(loglrposa), asNumeric(loglrposb), asNumeric(loglrposc),
    asNumeric(loglrnega), asNumeric(loglrnegb), asNumeric(loglrnegc)
  )
  
  ## Fix column names, pat -> pa etc.
  names(bounds) <- sub("(^p.{1}).", "\\1\\", names(bounds))
  
  return(bounds)
} # End bounds function

# crplot function to plot joint CL probabilites and box bounds ----
crplot <- function(bounds, cr_dists, n = 11, labels = T) {
  ca    <- bounds$ca[bounds$n == n]
  la    <- bounds$la[bounds$n == n]
  pa    <- bounds$pa_0.0[bounds$n == n]
  cb    <- bounds$cb[bounds$n == n]
  lb    <- bounds$lb[bounds$n == n]
  pb    <- bounds$pb_0.0[bounds$n == n]
  cbord <- bounds$cbord[bounds$n == n]
  lbord <- bounds$lbord[bounds$n == n]
  pc    <- bounds$pc_0.0[bounds$n == n]
  pt    <- paste0('pt', n)
  # cr    <- crs
  
  d <- cr_dists[[pt]] %>% 
    as_tibble(rownames = NA) %>% 
    rownames_to_column('C') %>% 
    gather('L', 'times', -C) %>% 
    mutate(L = as.integer(L),
           C = as.integer(C),
           p = times / sum(times),
           y = times,
           col = (times / max(times)) < 0.5 & (times / max(times)) > 0)
  
  p <- ggplot(d, aes(L, C, 
                     fill = times > 0,
                     alpha = times / max(times))) +
    geom_raster() +
    geom_rect(aes(xmin = 0.5,          # Anhoej box
                  xmax = la + 0.5,
                  ymin = ca - 0.5,
                  ymax = max(C) + 0.5),
              colour = '#F8766D',
              size = 1,
              fill   = NA) +
    geom_rect(aes(xmin = 0.5,          # Best box
                  xmax = lb + 0.5,
                  ymin = cb - 0.5,
                  ymax = max(C) + 0.5),
              linetype = 2,
              size = 1,
              colour   = '#00BA38',
              fill = NA) +
    geom_rect(aes(xmin = lbord + 0.5,  # Cut box, horizontal
                  xmax = lb + 0.5,
                  ymin = cb - 0.5,
                  ymax = cb + 0.5),
              colour   = '#619CFF',
              fill     = NA,
              linetype = 4,
              size = 1,
              na.rm    = T) +
    geom_rect(aes(xmin = lb - 0.5,     # Cut box, vertical
                  xmax = lb + 0.5,
                  ymin = cb - 0.5,
                  ymax = cbord - 0.5), 
              colour   = '#619CFF',
              linetype = 4,
              fill     = NA,
              na.rm    = T)
  
  if(labels) {
    p <- p +
      geom_text(aes(label = y,
                    colour = col),
                alpha = 1,
                size = 3)
  }
  
  p <- p +
    scale_y_reverse(breaks = 0:max(d$C)) +
    scale_x_continuous(position = 'top',
                       breaks = 1:max(d$L)) +
    scale_fill_manual(values = c('white', 'black')) +
    scale_colour_manual(values = c('white', 'black')) +
    theme_minimal() +
    theme(axis.title.y = element_text(angle = 0, hjust = 0),
          axis.title.x = element_text(hjust = 0),
          panel.grid.major = element_blank(),
          aspect.ratio = 1,
          legend.position = 'none') +
    labs(caption = paste0('N =', n, '\n',
                          'Rules specificity: ',
                          'Anhøj = ', round(pa, 3), ', ',
                          'Best Box = ', round(pb, 3), ', ',
                          'Cut Box = ', round(pc, 3)))
  plot(p)
} # End crplot function

# Create data objects ----
## Joint distribution matrices
cr_dists    <- crs(nmax, shifts)

## Data frame with limits and diagnostics for runs rules
cr_bounds <- bounds(cr_dists, target, target_shift)

## Bounds data in tall format
cr_bounds_tall <- cr_bounds %>% 
  select(-(ca:lbord)) %>%
  gather('key', 'val', -n) %>% 
  separate(key, c('test', 'shift'), '_') %>% 
  mutate(rule = substring(test, nchar(test)),
         test = substring(test, 1, nchar(test) - 1),
         shift = as.numeric(shift)) %>% 
  mutate(rule = fct_recode(rule, 
                           anhoej = 'a',
                           `best box` = 'b',
                           `cut box` = 'c')) %>% 
  spread(test, val)

# Figures ----
## Plot joint distribution matrix
crplot(cr_bounds, map(cr_dists$pt_0.0, asNumeric), 11, labels = T)

## Plot power function
ggplot(cr_bounds_tall, aes(n, 1 - p, colour = rule)) +
  geom_line(size = 1) +
  facet_wrap(~ shift, ncol = 4) +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  theme_minimal() +
  labs(title = 'Power function',
       y = 'Probability of signal',
       x = 'N')

## Plot specificity
ggplot(filter(cr_bounds_tall, shift == 0), aes(n, p, colour = rule)) +
  geom_line(size = 1) +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  theme_minimal() +
  labs(title = 'Specificity',
       y = 'Probability of true negative',
       x = 'N')

## Plot LR+
ggplot(filter(cr_bounds_tall, !is.na(loglrpos)), 
       aes(n, exp(loglrpos), colour = rule)) +
  geom_line(size = 0.75) +
  geom_hline(yintercept = 10) +
  facet_wrap(~ shift, ncol = 5) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  theme_minimal() +
  labs(title = 'Positive likelihood ratios',
       y = 'LR+',
       x = 'N')

## Plot LR-
ggplot(filter(cr_bounds_tall, !is.na(loglrpos)),
       aes(n, exp(loglrneg), colour = rule)) +
  geom_line(size = 0.75) +
  geom_hline(yintercept = 0.1) +
  facet_wrap(~ shift, ncol = 5) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  theme_minimal() +
  labs(title = 'Negative likelihood ratios',
       y = 'LR-',
       x = 'N')

