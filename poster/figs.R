library(qicharts2)
library(gridExtra)
library(tidyverse)

b <- readRDS('data/cr100_bounds.rds')
x <- read_rds('data/cr100_sym.rds')

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
              fill = NA) #+
  # geom_rect(aes(xmin = lbord + 0.5,  # Cut box, horizontal
  #               xmax = lb + 0.5,
  #               ymin = cb - 0.5,
  #               ymax = cb + 0.5),
  #           colour   = '#619CFF',
  #           fill     = NA,
  #           linetype = 4,
  #           size = 1,
  #           na.rm    = T) +
  # geom_rect(aes(xmin = lb - 0.5,     # Cut box, vertical
  #               xmax = lb + 0.5,
  #               ymin = cb - 0.5,
  #               ymax = cbord - 0.5), 
  #           colour   = '#619CFF',
  #           linetype = 4,
  #           fill     = NA,
  #           na.rm    = T)
  
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
                          'Best Box = ', round(pb, 3)))
  plot(p)
}

f1 <- crplot(11, labels = T)


# f1 <- readRDS('box11.rds')

p1 <- f1 +
  annotate('segment', 3.2, 0.2, xend = 5.9, yend = 3.8, 
           colour = 'gray30',
           arrow = arrow()) +
  annotate('text', 3, 0, label = 'Observed Fretheim et al. data') +
  labs(subtitle = 'Figure 1: Joint distribution of number of crossings and longest run')
p1

f2 <- read.delim('poster/f3.txt', sep = ';', dec = '.')
f2 <- qic(time, yprop, 
          data = f2, 
          part = 11, 
          point.size = 1.5, 
          cl = median(yprop[1:11]), 
          y.expand = 0.68,
          show.labels = F,
          part.labels = c('Pre-intervention', 'Transition / Post-intervention'),
          title = NULL,
          subtitle = 'Figure 2: Patients receiving treatment goals',
          ylab = 'Proportion',
          xlab = 'Months')

x2 <- summary(f2)[c('n.useful',
                    'CL', 
                    'longest.run', 
                    'longest.run.max',
                    'n.crossings',
                    'n.crossings.min')]
names(x2) <- c('N', 'Median', 'L', 'la', 'C', 'ca')
x2a <- x2[1,]
x2b <- x2[2,]
x2a <- tableGrob(x2a, 
                 theme = ttheme_minimal(base_size = 9,
                                        colhead=list(fg_params=list(fontface=1L))),
                 rows = NULL)
x2b <- tableGrob(x2b, 
                 theme = ttheme_minimal(base_size = 9,
                                        colhead=list(fg_params=list(fontface=1L))),
                 rows = NULL)
p2 <- grid.arrange(f2, x2a, x2b, 
                   layout_matrix = rbind(c(1, 1), c(2, 3)),
                   heights = c(6, 1))

f3 <- read.csv2('poster/crfirstchange.txt')
f3 <- qic(nr, Result, 
          data = f3, 
          part = 78,
          cl = median(Result[1:78]),
          show.labels = F,
          part.labels = c('Pre-change (lot 1)', 'Post-change (lot 2)'),
          title = NULL,
          subtitle = 'Figure 3: Creatinine in control materials',
          ylab = expression(paste(mu, 'mol / L')),
          xlab = 'Sequence # pre and post lot change')
x3 <- summary(f3)[c('n.useful',
                    'CL', 
                    'longest.run', 
                    'longest.run.max',
                    'n.crossings',
                    'n.crossings.min')]
names(x3) <- c('N', 'Median', 'L', 'la', 'C', 'ca')
x3a <- x3[1,]
x3a <- tableGrob(x3a, 
                 theme = ttheme_minimal(base_size = 9,
                                        colhead=list(fg_params=list(fontface=1L))),
                 rows = NULL)
x3b <- x3[2,]
x3b <- tableGrob(x3b, 
                 theme = ttheme_minimal(base_size = 9,
                                        colhead=list(fg_params=list(fontface=1L))),
                 rows = NULL)
p3 <- grid.arrange(f3, x3a, x3b, 
                   layout_matrix = rbind(c(1, 1), c(2, 3)),
                   heights = c(6, 1))

ggsave('poster/p1.pdf', p1, width = 15, height = 15, units = 'cm')
ggsave('poster/p2.pdf', p2, width = 17, height = 7.5, units = 'cm')
ggsave('poster/p3.pdf', p3, width = 17, height = 7.5, units = 'cm')

system2('pdfjam', c("--papersize '{37cm, 37cm}'", "poster/p1.pdf", "--outfile poster/p1jam.pdf"))
system2('pdfjam', c("--papersize '{40cm, 17.5cm}'", "poster/p2.pdf", "--outfile poster/p2jam.pdf"))
system2('pdfjam', c("--papersize '{40cm, 17.5cm}'", "poster/p3.pdf", "--outfile poster/p3jam.pdf"))

