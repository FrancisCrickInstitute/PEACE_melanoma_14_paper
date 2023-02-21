############################
#LIBRARIES 
############################

library(ape)
library(dplyr)
library(grid)
library(ggplotify)
library(ggsignif)
library(ggtreesignatures) #https://github.com/alexcoulton/ggtreesignatures
library(ggsci)
library(ggrepel)
library(viridis)
library(cowplot)

############################
#FUNCTIONS 
############################

ulapply = function(x, y){
  unlist(lapply(x, y))
}

uulapply = function(x, y){
  unname(unlist(lapply(x, y)))
}

reset.rownames = function(x){
  rownames(x) = 1:nrow(x)
  return(x)
}

rot.lab = function(size, angle){
  if(missing(size)) size = 10
  if(missing(angle)) angle = 40
  theme(axis.text.x = element_text(angle = angle, hjust = 1, size = size))
}

read.tsv = function(x, header = T, ...) read.delim(x, header = header, sep = '\t', ...)

############################
#RESOURCES 
############################

clusters1 = read.tsv('./dat/clusters.num.mutations.tsv')
signatures = readRDS('./dat/signatures.rds')
patients.df = read.csv('./dat/patients.remapped.csv')
patients = patients.df
treefile = read.tree('./dat/example.newick.tree')
treefile.clad = read.tree('./dat/example.newick.cladogram.tree')

list.of.patients = c('CRUKP2986', 'CRUKP1842', 'CRUKP6746', 'CRUKP5107', 'CRUKP6170', 'CRUKP2567', 'CRUKP1614', 'CRUKP9359', 'CRUKP9097', 'CRUKP1047', 'CRUKP6216', 'CRUKP1599', 'CRUKP6553', 'CRUKP2378')
list.of.patients.new = c("PC001", "PC002", "PC003", "PC004", "PC005", "PC006", "PC007", "PU008", "PU009", "PU010", "PA011", "PA012", "PA013", "PM014")
remap.patient.coord = match(
  read.csv('./dat/patients.remapped.csv')$Patient,
  list.of.patients
)

############################
#SIGNATURE PROCESSING 
############################


signatures = lapply(signatures, function(x){
  colnames(x) = c('node', 'signature', 'percentage')
  x$percentage = x$percentage * 100
  x$node = paste0('c', x$node)
  x$signature = as.character(x$signature)
  #x$signature = as.numeric(x$signature)
  x
})

signatures.comb = bind_rows(signatures)
sig.split = split(signatures.comb, signatures.comb$signature)
signatures.newmap = signatures[remap.patient.coord]

############################
#OTHER 
############################


treefile2 = lapply(treefile, function(x){
  x$node.label[x$node.label == 'NA'] = NA
  x$tip.label[x$tip.label == 'NA'] = NA
  x
})

treefile2 = treefile2[-1]

#get all ucscgb scale colours as manual for manipulation
g10 = ggplot(data.frame(x = factor(as.character(1:26), levels = as.character(1:26)), y = 1:26), aes(x = x, y = y, color = x)) +
  geom_point(size = 8) +
  scale_color_ucscgb(na.value = 'transparent') +
  geom_text(label = 1:26, color = 'black', aes(y = 1:26 + 1))

ucsc.colours = substr(ggplot_build(g10)$data[[1]]$colour, 1, 7)
ucsc.colours.orig = ucsc.colours
ucsc.colours[[1]] = '#000000'
ucsc.colours[[5]] = ucsc.colours.orig[[15]]
ucsc.colours[[15]] = ucsc.colours.orig[[5]]
ucsc.colours[[2]] = ucsc.colours.orig[[12]]
ucsc.colours[[12]] = ucsc.colours.orig[[2]]

#plots without mutational signatures
treeplots = lapply(treefile2[remap.patient.coord], function(x){
  treeplot1 = ggtree(x) +
    geom_point(aes(color = as.character(label)), size = 8) +
    scale_color_manual(values = ucsc.colours, na.value = 'transparent') +
    guides(color=guide_legend(ncol=2)) +
    scale_linetype_manual(values = c('solid', 'dashed', 'solid')) +
    theme_tree2() +
    theme(legend.position = 'none') +
    xlab('Mutations')
  treeplot1
})

#plots without mutational signatures
treeplots_cladograms = lapply(treefile2[remap.patient.coord], function(x){
  treeplot1 = ggtree(x, branch.length = 'none') +
    geom_point(aes(color = as.character(label)), size = 8) +
    scale_color_manual(values = ucsc.colours, na.value = 'transparent') +
    guides(color=guide_legend(ncol=2)) +
    scale_linetype_manual(values = c('solid', 'dashed', 'solid')) +
    theme(legend.position = 'none') 
  treeplot1
})

count1 = 1
treeplots_no_point = lapply(treefile2[remap.patient.coord], function(x){
  treeplot1 = ggtree(x, branch.length = 'none') +
    scale_color_manual(values = ucsc.colours, na.value = 'transparent') +
    guides(color=guide_legend(ncol=2)) +
    scale_linetype_manual(values = c('solid', 'dashed', 'solid')) +
    theme_tree2() +
    theme(legend.position = 'none') +
    ggtitle(patients.df$new_name_subtype[count1])
  count1 <<- count1 + 1
  treeplot1
})

sig.sum = bind_rows(signatures)
sig.sum.split = split(sig.sum, sig.sum$signature)

#take top quartile of signatures and remove rest for simplicity of visualization
sig.to.keep = names(sort(ulapply(sig.sum.split, function(x) sum(x$percentage)), decreasing = T)[1:7])

signatures.new = lapply(signatures, function(x){
  x.split = split(x, x$node)
  
  count1 = 1
  x.new = lapply(x.split, function(y){
    y$signature
    y2 = y[!y$signature %in% sig.to.keep, ]
    
    if(nrow(y2) > 0){
      y2$percentage = sum(y2$percentage)
      y2 = y2[1, ]
      y2$signature = '?'
      y3 = bind_rows(y[y$signature %in% sig.to.keep, ], y2)
    } else {
      y3 = y
    }
    
    count1 <<- count1 + 1
    return(y3)
  })
  bind_rows(x.new)
})

sig.to.keep = unique(bind_rows(signatures.new)$signature)

hex.colours = scales::viridis_pal()(length(as.character(unique(bind_rows(signatures)$signature))))
sig.names = as.character(unique(bind_rows(signatures)$signature))

good.sig.order = c('?', 'Signature.12', 'Signature.R3', 'Signature.18', 'Signature.U1', 'Signature.6','Signature.U2', 'Signature.19', 'Signature.11', 'Signature.13', 'Signature.15', 'Signature.5', 'Signature.9', 'Signature.3', 'Signature.10', 'Signature.1A', 'Signature.8', 'Signature.1B', 'Signature.17', 'Signature.R1', 'Signature.14', 'Signature.2', 'Signature.R2', 'Signature.21', 'Signature.4', 'Signature.20', 'Signature.16', 'Signature.7')
names(hex.colours) = good.sig.order

signatures = signatures.new
hex.colours = hex.colours[names(hex.colours) %in% sig.to.keep]

names(hex.colours) = gsub('Signature.', '', names(hex.colours))

hex.colours[[1]] = '#b3b3b3'
hex.colours
sort(names(hex.colours))
hex.colours = hex.colours[match(c('?', '1A', '2', '3', '6', '7', '11', '15'), names(hex.colours))]

hex.colours[c(4, 5, 8)] = c(
  '#99ccff',
  '#07f900',
  '#cccc99'
)

colordf = data.frame(color = hex.colours, x = 1:length(hex.colours), Signature = names(hex.colours))
colordf$Signature = gsub('Signature.', '', colordf$Signature)
colordf = reset.rownames(colordf)

head(colordf)
colordf$Signature = factor(colordf$Signature, levels = unique(colordf$Signature))

colordf.plot = ggplot(colordf, aes(x = x, y = 1, fill = Signature)) +
  scale_fill_manual(values = hex.colours, labels = c(
    'Unknown',
    '1A (Ageing)',
    '2 (APOBEC)',
    '3 (BRCA1/2 mutations)',
    '6 (DNA MMR deficiency)',
    '7 (UV light)',
    '11 (Temozolomide)',
    '15 (DNA MMR deficiency)'
  )) +
  geom_bar(stat = 'identity') 

legend1 = get_legend(colordf.plot) 

signatures = lapply(signatures, function(x){
  x$signature = gsub("Signature.", "", x$signature)
  x$signature[x$signature == 'unknown'] = '?'
  x
})

count1 = 1
tt_revised = Map(function(x, y){
  hex.colours2 = hex.colours[names(hex.colours) %in% y$signature]
  treeplot5 = ggtree(x, layout = 'signatures', signature.df = y, size = 0.5) +
    scale_color_manual(values = ucsc.colours, na.value = 'transparent') +
    guides(color=guide_legend(ncol=2)) +
    scale_linetype_manual(values = c('solid', 'dashed', 'solid')) +
    scale_fill_manual(values = hex.colours2) +
    guides(fill = guide_legend('Signatures', ncol = 2)) +
    theme_tree2() +
    theme(
      legend.position = 'none',
      plot.title = element_text(size = 9),
      axis.text = element_text(size = 9)
    ) +
    ggtitle(patients.df$patient.w.sub[count1]) +
    rot.lab(size = 8)
  count1 <<- count1 + 1
  treeplot5
}, treefile2[remap.patient.coord], signatures[remap.patient.coord])


############################
#CLUSTER ANALYSIS 
############################

clusters1$cluster = as.numeric(gsub('c', '', clusters1$cluster))
clusters1 = arrange(clusters1, patient, cluster)

clusters2 = bind_cols(clusters1, patients[match(clusters1$patient, patients$Patient), c('Subtype', 'Chemotherapy')])
clusters2
clusters2$cluster_type[clusters2$cluster_type != 'trunk'] = 'branch'
clusters2 = reset.rownames(clusters2)

clusters2[clusters2$cluster_type == 'trunk', ]$num_snps
clusters2[clusters2$cluster_type != 'trunk', ]$num_snps

clus.split = split(clusters2, clusters2$patient)

count1 = 1
clus.log.ratio = lapply(clus.split, function(x){
  patient = names(clus.split)[[count1]]
  ratio1 = log2(x$num_snps[x$cluster_type == 'branch'] / x$num_snps[x$cluster_type == 'trunk'])
  df1 = data.frame(patient, ratio1)
  count1 <<- count1 + 1
  df1
})

sum.plot.df = bind_rows(clus.log.ratio)
#add data on subtype, chemo etc
sum.plot.df = bind_cols(sum.plot.df, patients[match(sum.plot.df$patient, patients$Patient), ])
sum.plot.df2 = sum.plot.df[!sum.plot.df$Chemotherapy, ]
sum.plot.df[which(sum.plot.df$Chemotherapy), ]


buffer = c(-9, 9)
sum.plot.df2$Subtype = gsub('Cutaneous \\(Anorectal\\)', 'Cutaneous (non-UV)', sum.plot.df2$Subtype)
sum.plot.df2$Subtype = gsub('Cutaneous', 'Cut', sum.plot.df2$Subtype)

sum.plot1 = ggplot(sum.plot.df2, aes(x = Subtype, y = ratio1)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.2) +
  geom_signif(comparisons = list(
    c(1, 2),
    c(1, 4),
    c(2, 3)
    #c(3, 4)
  ),
  y_position = c(1, 3, 5, 7), map_signif_level = T) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
  ) +
  ylab('log2(# branch mut. / # trunk mut.)') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  coord_cartesian(ylim = buffer) +
  rot.lab(size = 8)

sum.plot2 = ggplot(sum.plot.df, aes(x = Chemotherapy, y = ratio1)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.2) +
  geom_signif(
    comparisons = list(c(1, 2)),
    map_signif_level = T,
    y_position = 6
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 6)
  ) +
  ylab('log2(# branch mut. / # trunk mut.)') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  coord_cartesian(ylim = buffer) +
  rot.lab(size = 8)

############################
#cowplot grid plots for paper 
############################

tree.fig.all.no.cladograms = plot_grid(
  #mutation trees plot
  do.call(
    plot_grid,
    c(
      tt_revised,
      list(ncol = 2)
    )
  ),
  #boxplots
  plot_grid(
    plot_grid(NULL, legend1, rel_heights = c(0.2, 1), ncol = 1),
    sum.plot1,
    sum.plot2,
    align = 'h',
    axis = 'tb',
    ncol = 3,
    labels = c('', 'b'),
    hjust = 1
  ),
  #parameters
  rel_heights = c(20, 7),
  axis = 't',
  align = 'h',
  labels = c('a'),
  vjust = 3,
  nrow = 2
)


pdf('./figure2.pdf', 8.8, 10)
print(tree.fig.all.no.cladograms)
dev.off()

