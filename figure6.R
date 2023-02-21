############################
#LIBRARIES 
############################

library(ctc)
library(GenomicRanges)
library(ggrepel)
library(ecodist)
library(adephylo)
library(ggsci)
library(ggtree)
library(ape)
library(dplyr)
library(reshape2)
library(cowplot)

############################
#RESOURCES
############################

CRUKP2567.tree = readRDS('./dat/CRUKP2567.tree.rds')
low.ploidy = readRDS('./dat/low.ploidy.rds')
high.ploidy = readRDS('./dat/high.ploidy.rds')
CRUKP2567.exome = readRDS('./dat/CRUKP2567.exome.rds')
hum.chrom = read.csv('./dat/hg19.chrom.lengths.csv')
CRUKP2567.medicc.tree = read.tree('./dat/CRUKP2567-exome-MEDICC_input_final_tree.new')

############################
#FUNCTIONS 
############################

ulapply = function(x, y){
  unlist(lapply(x, y))
}

uulapply = function(x, y){
  unname(unlist(lapply(x, y)))
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

make.key = function(chr, start.bp, end.bp){
  if(length(chr) == 0){
    return(as.character())
  } else {
    g = paste0(chr, ':', start.bp, '-', end.bp)
    g[grepl('NA', g)] = '0:1-1'
    return(g)
  }
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

generate.min.consensus.segments = function(combined.seg, allelic=F, join.up=T, remove.missing=F, print.chr=T){
  #takes in segmentation data and makes minimum consensus segments
  #args:
  #    combined.seg: data.frame with columns:
  #        sample: character vector with names of samples
  #        chr: numeric vector specifying names of chromosomes
  #        start: numeric vector of segment start positions
  #        end: numeric vector of segment end positions
  #        cn: numeric vector of copy number of segments
  #    allelic: boolean; indicates whether allelic copy number columns
  #        are included in the combined.seg dataframe (cn.a; cn.b)
  #    remove.missing: for MEDICC2; remove start positions from all samples 
  #        where CN information is not present in some samples
  
  combined.seg.split = split(combined.seg, combined.seg$chr)
  
  count1 = 1
  min.con.seg.all = mclapply(combined.seg.split, function(x){
    if(print.chr) print(count1)
    all.cut = unique(sort(c(unique(sort(as.integer(x$start))), unique(sort(as.integer(x$end))))))
    start = c(min(all.cut), all.cut[-1])
    start = start[-length(start)]
    end = all.cut[-1]
    end[1:(length(end) - 1)] = end[1:(length(end) - 1)] - 1
    cuts = data.frame(start, end)
    cuts$chr = unique(x$chr)
    cuts.ranges = GRanges(cuts)
    
    x.split = split(x, x$sample)
    
    cn.vals = lapply(x.split, function(y){
      y.ranges = GRanges(y)
      overlaps1 = as.data.frame(findOverlaps(cuts.ranges, y.ranges, type = 'within'))
      vals = y$cn[overlaps1$subjectHits]
      overlaps1$vals = vals
      overlaps1
      cuts$cn = ""
      cuts$cn[overlaps1$queryHits] = overlaps1$vals
      cuts$sample = unique(y$sample)
      
      if(allelic){
        vals = y$cn.a[overlaps1$subjectHits]
        overlaps1$vals = vals
        overlaps1
        cuts$cn.a = ""
        cuts$cn.a[overlaps1$queryHits] = overlaps1$vals
        
        vals = y$cn.b[overlaps1$subjectHits]
        overlaps1$vals = vals
        overlaps1
        cuts$cn.b = ""
        cuts$cn.b[overlaps1$queryHits] = overlaps1$vals
      }
      
      cuts
    })
    
    
    min.con.seg = bind_rows(cn.vals)
    
    count1 <<- count1 + 1
    min.con.seg
  }, mc.cores = 12)
  
  
  
  
  min.con.seg.all2 = bind_rows(min.con.seg.all)
  min.con.seg.all2 = min.con.seg.all2[!min.con.seg.all2$chr == 'X', ]
  
  medicc = min.con.seg.all2
  #medicc = medicc[c('sample', 'chr', 'start', 'end', 'cn')]
  #colnames(medicc)[1:2] = c('sample_id', 'chrom')
  medicc$start = as.integer(medicc$start)
  medicc$end = as.integer(medicc$end)
  head(medicc)
  
  medicc$start = as.character(medicc$start)
  medicc$end = as.character(medicc$end)
  head(medicc)
  
  to.rm1 = unique(medicc[is.na(medicc$cn), ]$start)
  to.rm2 = unique(medicc[medicc$cn == "", ]$start)
  
  if(remove.missing){
    medicc = medicc[!medicc$start %in% to.rm1, ]
    medicc = medicc[!medicc$start %in% to.rm2, ]
  }
  
  
  head(medicc)
  medicc$key = paste0(medicc$chr, ':', medicc$start, ':', medicc$end)
  
  medicc.split = split(medicc, factor(medicc$key, levels = unique(medicc$key)))
  
  
  keys.all = uulapply(medicc.split, function(x) paste0('chr', unique(x$chr), ':', paste0(x$cn, collapse = ':'), '_', unique(x$start), ':', unique(x$end)))
  
  keys.all = lapply(medicc.split, function(x){
    data.frame(
      cn.profile = paste0(x$cn, collapse = ':'),
      chr = unique(x$chr),
      start = unique(x$start),
      end = unique(x$end)
    )
  })
  
  keys.all = bind_rows(keys.all)
  keys.all
  keys.all$start = as.numeric(keys.all$start)
  keys.all$end = as.numeric(keys.all$end)
  keys.split = split(keys.all, keys.all$chr)
  
  
  keys.all
  medicc.split2 = medicc.split
  
  to.rm1 = as.numeric()
  i.same = 1
  for(i in 2:(length(medicc.split2) - 1)){
    #print(i)
    x = medicc.split2[[i.same]]
    x2 = medicc.split2[[i]]
    
    x.df = data.frame(
      cn.profile = paste0(x$cn, collapse = ':'),
      chr = unique(x$chr),
      start = unique(x$start),
      end = unique(x$end)
    )
    
    x2.df = data.frame(
      cn.profile = paste0(x2$cn, collapse = ':'),
      chr = unique(x2$chr),
      start = unique(x2$start),
      end = unique(x2$end)
    )
    
    if(
      x.df$cn.profile == x2.df$cn.profile &
      x.df$chr == x2.df$chr &
      as.numeric(x.df$end) == (as.numeric(x2.df$start) - 1)
    ){
      #print(paste0('updating ', i.same, ' to ', i))
      medicc.split2[[i.same]]$end = medicc.split2[[i]]$end
      to.rm1 = c(to.rm1, (i + 1))
    } else {
      i.same = i
    }
  }
  
  medicc.split2 = medicc.split2[-to.rm1]
  medicc.split2 = bind_rows(medicc.split2)
  
  #current = keys.all[[1]]
  #count1 = 1
  #to.change = ulapply(keys.all, function(x){
  #x.cn = multi.str.split(x, '_', 1)
  #x.pos = multi.str.split(x, '_', 2)
  #current.cn = multi.str.split(current, '_', 1)
  #current.pos = multi.str.split(current, '_', 2)
  #x.start = multi.str.split(x.pos, ':', 1)
  #x.end = multi.str.split(x.pos, ':', 2)
  #current.start = multi.str.split(current.pos, ':', 1)
  #current.end = multi.str.split(current.pos, ':', 2)
  #if(x.cn == current.cn & count1 != 1 & as.numeric(x.start) == (as.numeric(current.end) + 1)){
  #count1 <<- count1 + 1
  #return(T)
  #} else {
  #current <<- x
  #count1 <<- count1 + 1
  #return(F)
  #}
  #})
  
  
  #to.change
  #keys.all
  #which(to.change)
  
  #to.change
  #for(i in 1:length(medicc.split)){
  #if(to.change[[i]] == T){
  #medicc.split[[i]]$end = medicc.split[[i + 1]]$end
  #}
  #}
  
  
  #medicc.split = medicc.split[-(which(to.change) + 1)]
  
  #df3.split.longer.seg = lapply(df3.split, function(x){
  #x
  #one = medicc.split[[x$start]]
  #two = medicc.split[[x$end]]
  
  #if(all(one$sample_id == two$sample_id)){
  #one$end = two$end
  #one$key = paste0(one$chr, ':', one$start, ':', one$end)
  #} else{
  #print('error, sample names dont match')
  #stop()
  #}
  #return(one)
  
  #})
  
  #length(df3.split.longer.seg)
  df4 = bind_rows(medicc.split2)
  df4$start = as.numeric(df4$start)
  df4$end = as.numeric(df4$end)
  df4$cn = as.numeric(df4$cn)
  df4 = df4[order(df4$sample), ]
  if(allelic) df4 = df4[c('sample', 'chr', 'start', 'end', 'cn', 'cn.a', 'cn.b', 'key')]
  if(!allelic) df4 = df4[c('sample', 'chr', 'start', 'end', 'cn', 'key')]
  if(join.up == T) return(df4)
  if(join.up == F) return(min.con.seg.all2)
}

melt.ss = function(ss.dat){
  colnames(ss.dat) = gsub('.bam.*', '', colnames(ss.dat))
  ss.dat = melt(ss.dat, id = c('seqnames', 'start', 'end', 'width', 'strand', 'remap_start', 'remap_end'))
  
  ss.dat = ss.dat[c('variable', 'seqnames', 'remap_start', 'remap_end', 'value')]
  
  colnames(ss.dat) = c('sample', 'chr', 'start', 'end', 'cn')
  
  ss.dat = ss.dat[!is.na(ss.dat$start), ]
  
  ss.dat.split = split(ss.dat, ss.dat$sample)
  
  #remove normal cells
  ss.normal.cells = lapply(ss.dat.split, function(x){
    data.frame(
      ploidy = calculate.ploidy(x$chr, x$start, x$end, x$cn),
      wgii = calculate.wgii(x$chr, x$start, x$end, x$cn)
    )
  })
  ss.normal.cells = bind_rows(ss.normal.cells)
  ss.normal.cells[!abs(2 - ss.normal.cells$ploidy) < 0.05, ]
  
  ss.dat.split = ss.dat.split[!ss.normal.cells$wgii < 0.1]
  ss.dat = bind_rows(ss.dat.split)
  
  ss.dat$sample = as.character(ss.dat$sample)
  ss.dat = ss.dat[!ss.dat$chr == "X", ]
  ss.dat$chr = as.numeric(ss.dat$chr)
  ss.dat
}


perform.processing = function(exome.file, ss.file){
  pea04.exome = exome.file
  ss.pea04 = melt.ss(ss.file)
  pea04.comb = bind_rows(ss.pea04, pea04.exome)
  pea04.comb
}

############################
#CLONETREE 
############################


ucsc.colours = c(
  "#000000",
  "#FFFF00",
  "#FFCC00",
  "#00FF00",
  "#0000CC",
  "#CC33FF",
  "#99991E",
  "#999999",
  "#FF00CC",
  "#CC0000",
  "#FFCCCC",
  "#FF9900",
  "#CCFF00",
  "#358000",
  "#6699FF",
  "#99CCFF",
  "#00FFFF",
  "#CCFFFF",
  "#9900CC",
  "#CC99FF",
  "#996600",
  "#666600",
  "#666666",
  "#CCCCCC",
  "#79CC3D",
  "#CCCC99"
) 

treeplot.CRUKP2567 = ggtree(CRUKP2567.tree) +
  geom_point(aes(color = as.character(label)), size = 5) +
  scale_color_manual(values = ucsc.colours, na.value = 'transparent') +
  guides(color=guide_legend(ncol=2)) +
  scale_linetype_manual(values = c('solid', 'dashed', 'solid')) +
  theme_tree2() +
  theme(legend.position = 'none') +
  xlab('Mutations')

clonetree.plot = treeplot.CRUKP2567 +
  geom_hilight(
    data = data.frame(node = 2, type = 'Whole-genome doubled clone'),
    aes(
      node = node,
      fill = type
    ),
    extend = 10
  ) +
  scale_fill_manual(values = '#BBBBBB') +
  xlab('Mutations') +
  ggtitle('Clonal phylogeny') +
  theme(
    legend.position = 'bottom',
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(color = 'none')

## Irene check error:
# Error in (function (node, mapping = NULL, fill = "steelblue", alpha = 0.5,  : 
#                       unused argument (data = list(2, "Whole-genome doubled clone"))

############################
#CRUKP2567 BULK TREE (MEDICC) 
############################


df.highlight1 = data.frame(node = c(17), type = c('Whole-genome doubled clade'))

CRUKP2567.tree.plot = ggtree(CRUKP2567.medicc.tree) + geom_tiplab(offset = 1, size = 3) +
  #geom_text(aes(label = node)) +
  geom_hilight(
    data = df.highlight1,
    aes(
      node = node,
      fill = type
    ),
    alpha = 0.2
  ) +
  #geom_hilight(node = 17, fill = 'steelblue', alpha = 0.2) +
  geom_point2(aes(subset = (node == 8)), shape = 23, size = 3, fill = 'black') +
  geom_point2(aes(subset = (node == 11)), shape = 23, size = 3, fill = 'black') +
  ggtitle('Copy-number phylogeny') +
  theme(
    plot.margin = unit(c(0, 2, 0, 0), 'cm'),
    legend.title = element_blank(),
    legend.position = 'bottom',
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_fill_manual(values = 'steelblue') +
  coord_cartesian(clip = 'off')
#geom_cladelabel(node = 17, label = 'WGD'),

## Irene check error:
# Error in (function (node, mapping = NULL, fill = "steelblue", alpha = 0.5,  : 
#                       unused argument (data = list(17, "Whole-genome doubled clade"))

############################
#CRUKP2567 
############################


#low.ploidy = ss.files.hg19.second.run[[1]]
#high.ploidy = ss.files.hg19[[5]]

CRUKP2567.exome2 = CRUKP2567.exome[c('sample', 'Chromosome', 'Start.bp', 'End.bp', 'modal_total_cn')]
colnames(CRUKP2567.exome2) = c('sample', 'chr', 'start', 'end', 'cn')

CRUKP2567.comb = perform.processing(CRUKP2567.exome2, high.ploidy)
CRUKP2567.comb2 = perform.processing(CRUKP2567.exome2, low.ploidy)

CRUKP2567.comb2 = CRUKP2567.comb2[!CRUKP2567.comb2$sample %in% CRUKP2567.comb$sample, ]

CRUKP2567.comb = bind_rows(CRUKP2567.comb, CRUKP2567.comb2)

CRUKP2567.mcc = generate.min.consensus.segments(CRUKP2567.comb)

head(CRUKP2567.mcc)
unique(CRUKP2567.mcc$sample)

CRUKP2567.mcc.bulk = CRUKP2567.mcc[grepl('^M|^TB', CRUKP2567.mcc$sample), ]
CRUKP2567.mcc.bulk

CRUKP2567.mcc$sample = gsub('STU210622', 'FH', CRUKP2567.mcc$sample)
CRUKP2567.mcc$sample = gsub('STU210712', 'FL', CRUKP2567.mcc$sample)

CRUKP2567.subset = CRUKP2567.mcc[CRUKP2567.mcc$sample %in%
                             c(
                               'STU210622_11_074',
                               'STU210622_11_081',
                               'STU210622_11_076',
                               'M40-MT-DI-3d1',
                               'STU210712_3_084',
                               'TB55d1',
                               'STU210622_11_077',
                               'STU210712_3_084',
                               'STU210712_3_073'
                             ), ]

#pg = split(CRUKP2567.subset, CRUKP2567.subset$sample)

head(CRUKP2567.mcc)
unique(CRUKP2567.mcc$sample)

CRUKP2567.mcc = CRUKP2567.mcc[!CRUKP2567.mcc$sample %in% c('FL_3_066', 'FH_11_049', 'FL_3_091', 'FH_11_049'), ]

pg = split(CRUKP2567.mcc, CRUKP2567.mcc$sample)

pg.ploidy = ulapply(pg, function(x){
  head(x)
  calculate.ploidy
  x2 = x[!is.na(x$cn), ]
  round(calculate.ploidy(x2$chr, x2$start, x2$end, x2$cn))
})

pg = Map(function(x, y){
  x$cn.diff = y - x$cn
  x
}, pg, pg.ploidy)


hum.chrom.split = split(hum.chrom, hum.chrom$chr)

sample.ranges = lapply(hum.chrom.split, function(y){
  to.test = seq(1, y$length, 10000)
  
  r2 = GRanges(paste0(
    y$chr,
    ":",
    to.test,
    "-",
    to.test
  ))
  
})

sample.ranges = sample.ranges[!names(sample.ranges) %in% c('X', 'Y')]

x1 = GRanges()
for(i in sample.ranges){
  x1 = c(x1, i)
}

sample.ranges = x1

dist.cn.raw = lapply(pg, function(z){
  p17.ranges = GRanges(z)
  z2 = z[as.data.frame(findOverlaps(sample.ranges, p17.ranges))$subjectHits, ]
  z2$cn
}) %>% bind_rows %>% t %>% dist

dist.cn.relative = lapply(pg, function(z){
  p17.ranges = GRanges(z)
  z2 = z[as.data.frame(findOverlaps(sample.ranges, p17.ranges))$subjectHits, ]
  z2$cn.diff
}) %>% bind_rows %>% t %>% dist


pco.df1 = as.data.frame(cmdscale(dist.cn.raw))
pco.df1$sample = rownames(pco.df1)


pco.df1.relative = as.data.frame(cmdscale(dist.cn.relative))
pco.df1.relative$sample = rownames(pco.df1.relative)


############################
#Non-10kb segmentation profiles 
############################

dist.cn.non10kb = lapply(pg, function(z){
  z$cn
}) %>% bind_rows %>% t %>% dist

#pl(plot(hclust(dist.cn.non10kb)))

pco.non10kb = as.data.frame(cmdscale(dist.cn.non10kb))
pco.non10kb$sample = rownames(pco.non10kb)

#pl(
#ggplot(pco.non10kb, aes(x = V1, y = V2)) +
#geom_point() +
#geom_text_repel(aes(label = sample), max.overlaps = 20)
#)

############################
#HCLUST ON BULK ONLY 
############################

head(CRUKP2567.mcc.bulk)
pg.bulk = split(CRUKP2567.mcc.bulk, CRUKP2567.mcc.bulk$sample)

dist.cn.raw.bulk = lapply(pg.bulk, function(z){
  p17.ranges = GRanges(z)
  z2 = z[as.data.frame(findOverlaps(sample.ranges, p17.ranges))$subjectHits, ]
  z2$cn
}) %>% bind_rows %>% t %>% dist

#pl(plot(hclust(dist.cn.raw.bulk)))

############################
#HCLUST PLOTS FOR PUBLICATION 
############################

df.highlight1 = data.frame(node = c(16), type = c('Whole-genome doubled cluster'))
#pl(
#ggtree(as.phylo(hclust(dist.cn.raw.bulk))) +
#geom_tiplab() +
#theme(
#plot.margin = unit(c(0, 3, 0, 0), 'cm'),
#legend.title = element_blank(),
#legend.position = 'bottom',
#plot.title = element_text(hjust = 0.5)
#) +
#geom_treescale(width = 100) +
#coord_cartesian(clip = 'off') +
##geom_text(aes(label = node)) +
#geom_hilight(
#data = df.highlight1,
#aes(
#node = node,
#fill = type
#)
#) +
#geom_tree() +
#scale_fill_manual(values = '#BBBBBB'),
#5, 5
#)

df.highlight2 = data.frame(
  node = c(
    81,
    80,
    75,
    21,
    15,
    70,
    54,
    66,
    18,
    40,
    68,
    35,
    20,
    58
  ),
  type = c(
    'SS (FH)',
    'Bulk (WGD)',
    'SS (FH)',
    'SS (FH)',
    'SS (FH)',
    'Bulk (WGD)',
    'Bulk (WGD)',
    'SS (FL)',
    'SS (FH)',
    'SS (FL)',
    'Bulk (non-WGD)',
    'SS (FL)',
    'SS (FH)',
    'SS (FL)'
  )
)

ss.bulk.tree.plot = ggtree(as.phylo(hclust(dist.cn.raw))) +
  #geom_tiplab() +
  theme(
    plot.margin = unit(c(0, 3, 0, 0), 'cm'),
    legend.title = element_blank(),
    legend.position = 'bottom',
    plot.title = element_text(hjust = 0.5)
  ) +
  geom_treescale(width = 100) +
  #geom_text(aes(label = node)) +
  geom_hilight(
    data = df.highlight2,
    aes(
      node = node,
      fill = type
    ) 
  ) +
  geom_tree() +
  scale_fill_igv() +
  coord_cartesian(clip = 'off')

## Irene check error:
# Error in (function (node, mapping = NULL, fill = "steelblue", alpha = 0.5,  : 
#                       unused argument (data = list(c(81, 80, 75, 21, 15, 70, 54, 66, 18, 40, 68, 35, 20, 58), c("SS (FH)", "Bulk (WGD)", "SS (FH)", "SS (FH)", "SS (FH)", "Bulk (WGD)", "Bulk (WGD)", "SS (FL)", "SS (FH)", "SS (FL)", "Bulk (non-WGD)", "SS (FL)", "SS (FH)", "SS (FL)")))

#pl(clonetree.plot)

############################
#COMBINE PANELS 
############################

plot.test1 = plot_grid(
  plot_grid(clonetree.plot, NULL, ncol = 2, labels = c('a', 'b')),
  plot_grid(CRUKP2567.tree.plot, ss.bulk.tree.plot, ncol = 2, labels = c('c', 'd')),
  ncol = 1,
  rel_heights = c(1, 1, 0.5, 1),
  labels = c('', '', '', 'e', 'f'),
  align = 'v',
  axis = 'lr'
)

## Irene check error:
# Error in plot_grid(clonetree.plot, NULL, ncol = 2, labels = c("a", "b")) : 
#   object 'clonetree.plot' not found

pdf('./figure6.pdf', 8.8, 10)
print(plot.test1)
dev.off()
