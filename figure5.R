############################
#LIBRARIES 
############################

library(cowplot)
library(readxl)
library(viridis)
library(ggsignif)
library(RColorBrewer)
library(parallel)
library(readxl)
library(reshape2)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggrepel)

############################
#RESOURCES 
############################

list.of.patients = c('CRUKP2986', 'CRUKP1842', 'CRUKP6746', 'CRUKP5107', 'CRUKP6170', 'CRUKP2567', 'CRUKP1614', 'CRUKP9359', 'CRUKP9097', 'CRUKP1047', 'CRUKP6216', 'CRUKP1599', 'CRUKP6553', 'CRUKP2378')
list.of.patients.new = c("PC001", "PC002", "PC003", "PC004", "PC005", "PC006", "PC007", "PU008", "PU009", "PU010", "PA011", "PA012", "PA013", "PM014")
cosmic = read.csv('./dat/cosmic_gene_census.csv')
patients.df = read.csv('./dat/patients.remapped.csv')
segtab3 = readRDS('./dat/segtab3.fig5.rds')
abs.summary.panel = readRDS('./dat/abs.summary.panel.fig5.rds')
res2 = readRDS('./dat/res2.rds')

############################
#FUNCTIONS 
############################

is.cosmic = function(x) x %in% cosmic$Gene.Symbol

rot.lab = function(size, angle){
  if(missing(size)) size = 10
  if(missing(angle)) angle = 40
  theme(axis.text.x = element_text(angle = angle, hjust = 1, size = size))
}

arrange.patient = function(df, patient.col = 'patient'){
  #for arranging data in which multiple patients are in one file
  #only for new patient ordering
  #args:
  #    df = data.frame to arrange
  #    patient.col = string; name of patient column
  library(dplyr)
  df[[patient.col]] = factor(df[[patient.col]], levels = list.of.patients.new)
  df.split = split(df, df[[patient.col]])
  bind_rows(df.split)
}

patient.gsub = function(x){
  x = gsub('CRUKP2986', 'PC001', x)
  x = gsub('CRUKP1842', 'PC002', x)
  x = gsub('CRUKP2567', 'PC003', x)
  x = gsub('CRUKP9097', 'PC004', x)
  x = gsub('CRUKP6216', 'PC005', x)
  x = gsub('CRUKP6746', 'PC006', x)
  x = gsub('CRUKP1599', 'PC007', x)
  x = gsub('CRUKP1614', 'PU008', x)
  x = gsub('CRUKP6553', 'PU009', x)
  x = gsub('CRUKP5107', 'PU010', x)
  x = gsub('CRUKP9359', 'PA011', x)
  x = gsub('CRUKP2378', 'PA012', x)
  x = gsub('CRUKP1047', 'PA013', x)
  x = gsub('CRUKP6170', 'PM014', x)
  x
}

multi.str.split = function(x, character.split, split.id=NA){
  #An implementation of strsplit() that allows selection of a particular element of the output value of strsplit(), 
  #for all elements of the input vector to which the split is applied.
  #--------------
  #args:
  #x = character vector
  #character.split = character specifying what character to split with 
  #split.id = integer specifying which element of the split you want
  
  if(!is.na(split.id)){
    return(unlist(lapply(x, function(y){
      g = strsplit(y, character.split)[[1]]
      paste(g[split.id], collapse = character.split)
    })))
  } else {
    g = strsplit(x, character.split)
    g2 = unlist(lapply(g, function(z) z[length(z)]))
    return(g2)
  }
}

make.key = function(chr, start.bp, end.bp){
  if(length(chr) == 0){
    return(as.character())
  } else {
    g = paste0(chr, ':', start.bp, '-', end.bp)
    g[grepl('NA', g)] = '0:1-1'
    return(g)
  }
}


calculate.wgii = function(chr, start, end, cn){
  #calculates the wgii of a sample from the weighted copy number states of a series of segments
  #args:
  #    chr: integer vector: chromosome for each segment
  #    start: integer vector, start positions for each segment
  #    end: integer vector, end positions for each segment
  #    cn: integer vector, total copy number for each segment
  x = data.frame(chr, start, end, cn)
  x$length = x$end - x$start
  x.split = split(x, x$cn)
  
  x.weights = lapply(x.split, function(x){
    data.frame(
      sum = sum(x$length),
      cn = unique(x$cn)
    )
  }) %>% bind_rows
  
  x.weights$weight = x.weights$sum / sum(x.weights$sum)
  x.weights
  x.wgii = x.weights[-which.max(x.weights$weight), ]
  wgii = sum(x.wgii$weight)
  wgii
}

calculate.ploidy = function(chr, start, end, cn){
  #calculates the ploidy of a sample from the weighted copy number states of a series of segments
  #args:
  #    chr: integer vector: chromosome for each segment
  #    start: integer vector, start positions for each segment
  #    end: integer vector, end positions for each segment
  #    cn: integer vector, total copy number for each segment
  y.ranges = GRanges(make.key(chr, start, end))
  tot.cov.y = sum(width(y.ranges))
  
  y.width = width(y.ranges)
  
  if(sum(y.width / tot.cov.y) != 1) stop('error, coverage calculation fail')
  
  ploidy = sum(cn * (y.width / tot.cov.y))
  return(ploidy)
}

############################
#calculate ploidy and wgii 
############################


segtab3.split = split(segtab3, segtab3$tumor.key)

wgii = lapply(segtab3.split, function(x){
  ploidy = calculate.ploidy(x$Chromosome, x$Start.bp, x$End.bp, x$modal_total_cn)
  wgii = calculate.wgii(x$Chromosome, x$Start.bp, x$End.bp, x$modal_total_cn)
  return(data.frame(ploidy = ploidy, wgii = wgii, tumor.key = unique(x$tumor.key)))
}) %>% bind_rows


wgii = left_join(wgii, abs.summary.panel, by = 'tumor.key')

wgii$Genome.doublings = factor(wgii$Genome.doublings)

wgii_wgd_plot = ggplot(wgii, aes(x = Genome.doublings, y = wgii)) + 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.1) +
  geom_signif(comparisons = list(c(1, 2), c(1, 3), c(2, 3)),
              map_signif_level = T,
              y_position = c(0.9, 1, 1.1)) +
  ylab('wGII') +
  xlab('Genome doublings') +
  scale_y_continuous(labels = c(0, 0.5, 1), breaks = c(0, 0.5, 1)) +
  coord_cartesian(ylim = c(0, 1.2))
#theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))

wgd.pie = data.frame(
  num = c(
    w.doub = length(which(abs.summary.panel$Genome.doublings > 0)),
    wo.doub = length(which(abs.summary.panel$Genome.doublings == 0))
  ),
  cat = c('w.doub', 'wo.doub')
)

wgd.pie2 = as.data.frame(table(abs.summary.panel$Genome.doublings))
colnames(wgd.pie2) = c('num.wgd', 'freq')
wgd.pie2$num.wgd = factor(wgd.pie2$num.wgd, levels = c(0, 1, 2))

wgd.pie.plot = ggplot(wgd.pie2, aes(x = "", y = freq, group = num.wgd, fill = num.wgd)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0, direction = 1) +
  ggtitle('WGD frequency') +
  theme_bw() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  guides(fill = guide_legend(title = 'Number of WGD events')) +
  scale_fill_brewer(palette = 'Blues', direction = 1) +
  theme(
    panel.border = element_blank(),
    legend.position = 'bottom',
    axis.ticks = element_blank()
  ) +
  ylab('') +
  xlab('')

abs.summary.panel$patient = multi.str.split(abs.summary.panel$tumor.key, '_', 1)

abs.split1 = split(abs.summary.panel, abs.summary.panel$patient)
abs.wgd.bar = bind_rows(lapply(abs.split1, function(x) data.frame(table(x$Genome.doublings), patient = unique(x$patient))))
abs.wgd.bar$Var1 = 
  factor(abs.wgd.bar$Var1, levels = c(2, 1, 0))

colnames(abs.wgd.bar) = c('WGD', 'freq', 'patient')

abs.wgd.bar$patient = patient.gsub(abs.wgd.bar$patient)
abs.wgd.bar = arrange.patient(abs.wgd.bar)

abs.wgd.bar$patient = patients.df$patient.w.sub[match(abs.wgd.bar$patient, patients.df$new_name_subtype)]
abs.wgd.bar$patient = factor(abs.wgd.bar$patient, levels = patients.df$patient.w.sub)

abs.wgd.bar.plot = ggplot(abs.wgd.bar, aes(x = patient, y = freq, fill = WGD)) +
  #geom_bar(stat = 'identity', position = position_stack(reverse = T)) +
  geom_bar(stat = 'identity') +
  scale_fill_brewer(palette = 'Blues', direction = -1) +  
  theme_bw() +
  theme(panel.grid = element_blank()) +
  rot.lab(size = 8) +
  xlab('Patient') +
  ylab('Frequency')

############################
#MORE JUN ANALYSIS 
############################

breaks_fun = function(x){
  if(abs(min(x) - max(x)) < 100){
    g1 = seq(
      plyr::round_any(min(x), 25),
      plyr::round_any(max(x), 25),
      50
    )
    return(g1)
  }
  seq(
    plyr::round_any(min(x), 25),
    plyr::round_any(max(x), 25),
    100
  )
}



res2$perm_pval = (res2$gs_ge + 1) / 900
res2$perm_pval[res2$perm_pval > 1] = 1

res3 = res2[which(res2$perm_pval < 0.05), ]
sig.genes = unique(res3$gene)
res3 = res2[res2$gene %in% sig.genes, ]
res3.split = split(res3, res3$gene)

res2.amp = res2[res2$type == 'Amp', ]

res2.amp = res2.amp[c('type', 'outcome', 'gene', 'chrom', 'start.pos', 'end.pos', 'perm_pval')]
res2.amp = res2.amp[res2.amp$outcome != 'all', ]
res2.amp$chrom = as.numeric(res2.amp$chrom)
res2.amp = arrange(res2.amp, chrom, start.pos)

custom.theme = theme(
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = 'white', color = 'white'),
  panel.spacing.x = unit(0, 'cm'),
  strip.background = element_rect(fill = 'white', color = 'white')
) +
  rot.lab(size = 8)

res2.amp$start.pos = res2.amp$start.pos / 1000000
res2.amp$end.pos = res2.amp$end.pos / 1000000

res2.amp.split = split(res2.amp, res2.amp$chrom)

panel.background.df = data.frame(
  chrom = unique(res2.amp$chrom),
  panel.color = rep(c(1, 2), length(unique(res2.amp$chrom)) / 2)
)

panel.background.df$panel.color = as.character(panel.background.df$panel.color)

res2.amp$panel.background = panel.background.df$panel.color[match(res2.amp$chrom, panel.background.df$chrom)]
res2.amp$panel.background = as.character(res2.amp$panel.background)

plot.amp = ggplot(res2.amp, aes(x = start.pos, y = -log(perm_pval), group = outcome, color = outcome)) + 
  geom_line() +
  geom_rect(
    inherit.aes = F,
    data = panel.background.df,
    aes(fill = panel.color),
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.2
  ) +
  scale_fill_manual(values = c('#FFFFFF', '#999999')) +
  facet_grid(cols = vars(chrom), scales = 'free_x', space = 'free_x') +
  geom_hline(yintercept = -log(0.05), linetype = 'dashed') +
  ggtitle('Amplifications') +
  xlab('Position (Mb)') +
  ylab('-log(p-value)') +
  scale_color_manual(values = c('#0062ff', '#ff0000'), labels = c('NR', 'R')) +
  scale_x_continuous(breaks = breaks_fun) +
  guides(color = guide_legend(title = 'Outcome')) +
  guides(fill = 'none') +
  custom.theme 

res2.amp.1.and.8 = res2.amp[res2.amp$chr %in% c(1, 8), ]
res2.amp.1.and.8.text1 =
  res2.amp.1.and.8[
    (res2.amp.1.and.8$outcome == 'nr' & res2.amp.1.and.8$chr == 8) |
      (res2.amp.1.and.8$outcome == 'resp' & res2.amp.1.and.8$chr == 1), ]
res2.amp.text = res2.amp.1.and.8.text1 %>% filter(perm_pval < 0.05)
res2.amp.text = res2.amp.text[is.cosmic(res2.amp.text$gene), ]

cosmic[cosmic$Gene.Symbol == 'PBX1', ]

head(res2.amp.1.and.8)
res2.del = res2[res2$type == 'Del', ]

res2.del = res2.del[c('type', 'outcome', 'gene', 'chrom', 'start.pos', 'end.pos', 'perm_pval')]
res2.del = res2.del[res2.del$outcome != 'all', ]
res2.del$chrom = as.numeric(res2.del$chrom)
res2.del = arrange(res2.del, chrom, start.pos)

res2.del$start.pos = res2.del$start.pos / 1000000
res2.del$end.pos = res2.del$end.pos / 1000000

panel.background.df = data.frame(
  chrom = unique(res2.del$chrom),
  panel.color = rep(c(1, 2), length(unique(res2.del$chrom)) / 2)
)

panel.background.df$panel.color = as.character(panel.background.df$panel.color)

res2.del$panel.background = panel.background.df$panel.color[match(res2.del$chrom, panel.background.df$chrom)]
res2.del$panel.background = as.character(res2.del$panel.background)

plot.del = ggplot(res2.del, aes(x = start.pos, y = -log(perm_pval), group = outcome, color = outcome)) + 
  geom_line() +
  geom_rect(
    inherit.aes = F,
    data = panel.background.df,
    aes(fill = panel.color),
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.2
  ) +
  scale_fill_manual(values = c('#FFFFFF', '#999999')) +
  facet_grid(cols = vars(chrom), scales = 'free_x', space = 'free_x') +
  geom_hline(yintercept = -log(0.05), linetype = 'dashed') +
  ggtitle('Deletions') +
  xlab('Position (Mb)') +
  ylab('-log(p-value)') +
  scale_color_manual(values = c('#0062ff', '#ff0000'), labels = c('NR', 'R')) +
  scale_x_continuous(breaks = breaks_fun) +
  guides(color = guide_legend(title = 'Outcome')) +
  guides(fill = 'none') +
  custom.theme 

res2.all = bind_rows(res2.amp, res2.del)

res2.all$type = gsub('Amp', 'Amplifications', res2.all$type)
res2.all$type = gsub('Del', 'Deletions', res2.all$type)

plot.all = ggplot(res2.all, aes(x = start.pos, y = -log(perm_pval), group = outcome, color = outcome)) + 
  geom_line() +
  geom_rect(
    inherit.aes = F,
    data = panel.background.df,
    aes(fill = panel.color),
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.2
  ) +
  scale_fill_manual(values = c('#FFFFFF', '#999999')) +
  facet_grid(cols = vars(chrom), rows = vars(type), scales = 'free_x', space = 'free_x') +
  geom_hline(yintercept = -log(0.05), linetype = 'dashed') +
  xlab('Position (Mb)') +
  ylab('-log(p-value)') +
  scale_color_manual(values = c('#0062ff', '#ff0000'), labels = c('NR', 'R')) +
  scale_x_continuous(breaks = breaks_fun) +
  guides(color = guide_legend(title = 'Outcome')) +
  guides(fill = 'none') +
  custom.theme +
  theme(legend.position = 'bottom')

res.del.15 = res2.del %>% filter(chrom == 15)
res.del.15[res.del.15$perm_pval < 0.05, ]

sig.genes.amp = res2.amp %>% filter(perm_pval < 0.05) %>% dplyr::select(gene)
sig.genes.amp = sig.genes.amp[[1]]
sig.genes.amp = sig.genes.amp[which(is.cosmic(sig.genes.amp))]

sig.genes.del = res2.del %>% filter(perm_pval < 0.05) %>% dplyr::select(gene)
sig.genes.del = sig.genes.del[[1]]
sig.genes.del = sig.genes.del[which(is.cosmic(sig.genes.del))]

sig.genes.amp = cosmic[cosmic$Gene.Symbol %in% sig.genes.amp, c('Gene.Symbol', 'Role.in.Cancer', 'Genome.Location', 'Chr.Band')]
sig.genes.amp$type = 'amp'
sig.genes.del = cosmic[cosmic$Gene.Symbol %in% sig.genes.del, c('Gene.Symbol', 'Role.in.Cancer', 'Genome.Location', 'Chr.Band')]
sig.genes.del$type = 'del'

sig.genes.all = bind_rows(sig.genes.amp, sig.genes.del)
sig.genes.all$chr = as.numeric(multi.str.split(sig.genes.all$Genome.Location, ":", 1))
sig.genes.all = sig.genes.all[sig.genes.all$chr %in% c(1, 8), ]

sig.genes.all$start = as.numeric(multi.str.split(multi.str.split(sig.genes.all$Genome.Location, ":", 2), "-", 1))
sig.genes.all$end = as.numeric(multi.str.split(multi.str.split(sig.genes.all$Genome.Location, ":", 2), "-", 2))
sig.genes.all = arrange(sig.genes.all, type, chr, start)

sig.genes.all.split = split(sig.genes.all, as.character(sig.genes.all$chr))

sig.genes.all.split = lapply(sig.genes.all.split, function(x){
  end1 = max(x$end)
  x$line2 = seq(100000, end1, end1 / nrow(x))
  x
})

sig.genes.all = bind_rows(sig.genes.all.split)
sig.genes.all$chr = factor(sig.genes.all$chr, levels = 1:22)

sig.genes.all.melt = melt(sig.genes.all, id = c('Gene.Symbol', 'Role.in.Cancer', 'Genome.Location',
                                                'Chr.Band', 'type', 'chr', 'end'))

sig.genes.all.melt$y[1:(nrow(sig.genes.all.melt) / 2)] = 1
sig.genes.all.melt$y[((nrow(sig.genes.all.melt) / 2) + 1):nrow(sig.genes.all.melt)] = 2

sig.genes.text.df = sig.genes.all.melt[sig.genes.all.melt$y == 2, ]

sig.genes.plot1 = ggplot(sig.genes.all, aes(x = start, y = 1)) +
  facet_grid(cols = vars(chr), scales = 'free_x') +
  theme_classic() +
  rot.lab() +
  geom_line(
    data = sig.genes.all.melt,
    inherit.aes = F,
    aes(x = value, y = y, group = Gene.Symbol),
    color = '#c2c2c2'
  ) +
  geom_point() +
  coord_cartesian(ylim = c(0.5, 4)) +
  geom_text(
    inherit.aes = F,
    data = sig.genes.text.df,
    aes(x = value, y = 2.03, label = Gene.Symbol), 
    size = 2,
    angle = 45,
    hjust = 0
  ) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    #plot.margin = unit(rep(-20, 4), 'cm'),
    panel.spacing.x = unit(0, 'cm')
  )

sig.genes.all.melt$value.mb = sig.genes.all.melt$value / 1000000
colnames(sig.genes.all.melt)[6] = 'chrom'

res2.temp = bind_rows(
  res2.amp.1.and.8[res2.amp.1.and.8$chrom == 1 & res2.amp.1.and.8$outcome == 'resp', ],
  res2.amp.1.and.8[res2.amp.1.and.8$chrom == 8 & res2.amp.1.and.8$outcome == 'nr', ]
)

sig.genes.all.melt$perm_pval = -log(res2.temp[match(sig.genes.all.melt$Gene.Symbol, res2.temp$gene), ]$perm_pval)

#space between line pointing to gene labels and data line
sig.genes.all.melt[sig.genes.all.melt$variable == 'start', ]$perm_pval =
  sig.genes.all.melt[sig.genes.all.melt$variable == 'start', ]$perm_pval + 0.05

#height of line pointing to gene labels
sig.genes.all.melt[sig.genes.all.melt$variable == 'line2', ]$perm_pval =
  sig.genes.all.melt[sig.genes.all.melt$variable == 'line2', ]$perm_pval + 0.3

s.test2 = sig.genes.all.melt %>%
  filter(chrom == 1, variable == 'start') 

s.test3 = s.test2[(which(s.test2$Gene.Symbol == 'PBX1') + 1):nrow(s.test2), ]
bumper.genes = s.test3$Gene.Symbol

#x-axis buffer for genes on the right-hand side of chromosome 1
sig.genes.all.melt$value.mb[sig.genes.all.melt$Gene.Symbol %in% bumper.genes & sig.genes.all.melt$variable == 'line2'] = 
  sig.genes.all.melt$value.mb[sig.genes.all.melt$Gene.Symbol %in% bumper.genes & sig.genes.all.melt$variable == 'line2'] + 50

#x-axis adjustment for all gene labels and lines on chromosome 8
sig.genes.all.melt$value.mb[sig.genes.all.melt$variable == 'line2' & sig.genes.all.melt$chrom == 8] =
  sig.genes.all.melt$value.mb[sig.genes.all.melt$variable == 'line2' & sig.genes.all.melt$chrom == 8] + 10

#x-axis buffer for right hand side gene on chromosome 8
sig.genes.all.melt$value.mb[sig.genes.all.melt$variable == 'line2' & sig.genes.all.melt$Gene.Symbol == 'FAM135B'] =
  sig.genes.all.melt$value.mb[sig.genes.all.melt$variable == 'line2' & sig.genes.all.melt$Gene.Symbol == 'FAM135B'] + 25

plot.amp.chr1and8 = ggplot(res2.amp.1.and.8, aes(x = start.pos, y = -log(perm_pval), group = outcome, color = outcome)) + 
  geom_hline(yintercept = -log(0.05), linetype = 'dashed') +
  custom.theme +
  geom_point(
    inherit.aes = F,
    aes(x = start.pos, y = -log(perm_pval)),
    data = res2.amp.text
  ) +
  geom_rect(
    inherit.aes = F,
    data = panel.background.df[panel.background.df$chrom %in% c(1, 8), ],
    aes(fill = panel.color),
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.2
  ) +
  geom_line(
    data = sig.genes.all.melt,
    inherit.aes = F,
    aes(x = value.mb, y = perm_pval, group = Gene.Symbol),
    color = '#c2c2c2',
    linetype = 'dashed'
  ) +
  geom_text(
    data = sig.genes.all.melt[sig.genes.all.melt$variable == 'line2', ],
    inherit.aes = F,
    #adjustment here of the buffer between line end and text label
    aes(x = value.mb, y = perm_pval + 0.05, label = Gene.Symbol, angle = 300, hjust = 1),
    color = '#000000',
    size = 2
  ) +
  geom_line() +
  scale_color_manual(values = c('#0062ff', '#ff0000'), labels = c('NR', 'R')) +
  scale_fill_manual(values = c('#FFFFFF', '#999999')) +
  guides(fill = 'none') +
  guides(color = guide_legend(title = 'Outcome')) +
  xlab('Position (Mb)') +
  ylab('-log(p-value)') +
  #facet_wrap(~ chrom, ncol = 2, scales = 'free_x', space = 'free')
  facet_grid(cols = vars(chrom), scales = 'free_x', space = 'free_x') +
  coord_cartesian(ylim = c(0, 6.5)) +
  theme(legend.position = 'none')

genes.chrom1 = res2.amp %>%
  filter(chrom == 1, perm_pval < 0.05) %>%
  dplyr::select(gene) %>% unique

chr1.cosmic = genes.chrom1[which(unlist(lapply(genes.chrom1, is.cosmic))), ]
cosmic[cosmic$Gene.Symbol %in% chr1.cosmic, c('Gene.Symbol', 'Role.in.Cancer')]

top.rowv3 = plot_grid(abs.wgd.bar.plot, wgii_wgd_plot, ncol = 2, labels = c('a', 'b'))

plot.test1 = plot_grid(
  top.rowv3,
  plot.all,
  plot.amp.chr1and8,
  align = 'v',
  axis = 'r',
  nrow = 3
)

pdf('./figure5.pdf', 8.8, 10)
print(plot.test1)
dev.off()

