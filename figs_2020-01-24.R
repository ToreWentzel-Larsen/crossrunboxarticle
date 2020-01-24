# Setup ----
library(tidyverse)
library(qicharts2)

save.plots <- T

x <- read_rds('data/cr100_sym.rds')
b <- read_rds('data/cr100_bounds.rds')

## Tall parameter value data
vals <- b %>% 
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

# Run chart ----
set.seed(4)
fig1_run <- qic(rnorm(20, 3),
                title = NULL, #'Run chart',
                ylab = 'Indicator value',
                xlab = 'Time/order')

# Plot power function ----
fig4_pwr <- ggplot(vals, aes(n, 1 - p, colour = rule)) +
  geom_line(size = 1) +
  facet_wrap(~ shift, ncol = 4) +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  theme_minimal() +
  labs(title = NULL, #'Power function',
       y = 'Probability of signal',
       x = 'N')

## Specificity
fig2_spec <- ggplot(filter(vals, shift == 0), aes(n, p, colour = rule)) +
  geom_line(size = 1) +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  theme_minimal() +
  labs(title = NULL, #'Specificity (no shif)',
       y = 'Probability of true negative',
       x = 'N')

## Sensitivity
fig5_sens <- ggplot(filter(vals, shift == 0.8), aes(n, 1 - p, colour = rule)) +
  geom_line(size = 1) +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  theme_minimal() +
  labs(title = NULL, #'Sensitivity (shift = 0.8 SD)',
       y = 'Probability of true positive',
       x = 'N')

# Plot likelihood ratios ----
## LR+
fig6_lrpos <- ggplot(filter(vals, !is.na(loglrpos)), 
                     aes(n, exp(loglrpos), colour = rule)) +
  geom_line(size = 0.75) +
  geom_hline(yintercept = 10) +
  facet_wrap(~ shift, ncol = 5) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  theme_minimal() +
  labs(title = NULL, #'Positive likelihood ratios',
       y = 'LR+',
       x = 'N')

## LR-
fig7_lrneg <- ggplot(filter(vals, !is.na(loglrpos)),
                     aes(n, exp(loglrneg), colour = rule)) +
  geom_line(size = 0.75) +
  geom_hline(yintercept = 0.1) +
  facet_wrap(~ shift, ncol = 5) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  theme_minimal() +
  labs(title = NULL, #'Negative likelihood ratios',
       y = 'LR-',
       x = 'N')

# Function to plot LC box figures ----
crplot <- function(n = 11, labels = T) {
  ca    <- b$ca[b$n == n]
  la    <- b$la[b$n == n]
  pa    <- b$pa_0.0[b$n == n]
  cb    <- b$cb[b$n == n]
  lb    <- b$lb[b$n == n]
  pb    <- b$pb_0.0[b$n == n]
  cbord <- b$cbord[b$n == n]
  lbord <- b$lbord[b$n == n]
  pc    <- b$pc_0.0[b$n == n]
  
  d <- x[[paste0('pt', n)]] %>% 
    as_tibble(rownames = NA) %>% 
    rownames_to_column('C') %>% 
    gather('L', 'times', -C) %>% 
    mutate(L = as.numeric(L),
           C = as.numeric(C),
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
  return(p)
}

fig3_box11 <- crplot(11, labels = T)

if(save.plots) {
  w <- 7
  ggsave('figures/fig1_run.tiff',    fig1_run,   width = w, height = w / 2, compression = 'lzw')
  ggsave('figures/fig2_spec.tiff',   fig2_spec,  width = w, height = w / 2, compression = 'lzw')
  ggsave('figures/fig3_box11.tiff',  fig3_box11, width = w, height = w, compression = 'lzw')
  ggsave('figures/fig4_pwr.tiff',    fig4_pwr,   width = w, height = w * 6 / 8, compression = 'lzw')
  ggsave('figures/fig5_sens.tiff',   fig5_sens,  width = w, height = w / 2, compression = 'lzw')
  ggsave('figures/fig6_lrpos.tiff',  fig6_lrpos, width = w, height = w * 5 / 8, compression = 'lzw')
  ggsave('figures/fig7_lrneg.tiff', fig7_lrneg,  width = w, height = w * 5 / 8, compression = 'lzw')
}
