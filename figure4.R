options(stringsAsFactors = F)

library(ggplot2)
library(ggpubr)
library(tidyverse)
library(reshape2)
library(ape)
library(ggtree)
library(GenomicRanges)
library(circlize)
library(ComplexHeatmap)
library(cowplot)



#### Panel a - KIT expression vs CN ####

# Import data
KIT_CN_TPM<-read.delim("./dat/KIT_totCN_TPM.txt",header=T)

# Order relevant patient on top
KIT_CN_TPM$Patient<-factor(KIT_CN_TPM$Patient,levels=c("CRUKP9359","Other"),
                           ordered=T)

# Plot TMP vs CN
KIT_corrplot<-ggplot(KIT_CN_TPM,aes(KIT_total_CN,KIT_TPM),alpha=0.8)+
  geom_point(aes(color=Patient),size=3)+
  theme_bw(base_size = 10)+
  stat_cor(p.accuracy = 0.0001,)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="black")+
  ylab("KIT expression (TPM)")+
  xlab("KIT copy number")+
  scale_colour_manual(values=c("#e31a1c","darkgrey"),name="Patient")+
  theme(axis.title.x = element_text(size=9),
        axis.title.y = element_text(size=9),
        legend.text = element_text(size=7),
        legend.title = element_text(size=8),
        legend.key.height= unit(2, 'mm'),
        legend.key.width= unit(4, 'mm'),
        # legend.box.spacing = unit(0, "pt"),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        legend.position = c(0.82,0.12))



#### Panel c - sc CN tree ####

# Import CN segments from patient's sc WGS

sc_cn <- read.delim("./dat/KIT_CRUKP9359_singlecell_CN.tsv")

#KIT position in hg19
kit.pos = data.frame(
  start = 55524124,
  end = 55606881,
  chr = 4,
  gene = 'KIT'
)


# Transform sample IDs to cell n
sc_cn$sample = gsub('.*?_.*?_.*?_(.*)', 'c\\1', sc_cn$sample)
sc_cn$sample = paste0('cell_', as.numeric(as.factor(sc_cn$sample)))

# Segment ID
sc_cn$seg = paste0(sc_cn$chr, '_', sc_cn$start, '_', sc_cn$end)
sc_cn.orig = sc_cn

#remove KIT segments so they don't drive the clustering
sc_cn = sc_cn[!(sc_cn$chr == 4 & sc_cn$start == 55182555), ]
sc_cn = sc_cn[!(sc_cn$chr == 4 & sc_cn$start == 56042707), ]


# Hierarchical clustering
hclust.test = dcast(sample ~ seg, value.var = 'cn', data = sc_cn)
hclust.test[1:5, 1:5]
rownames(hclust.test) = hclust.test$sample
hclust.test = hclust.test[, -1]
hclust.test = as.matrix(hclust.test)

hclust.test = hclust.test[!round(rowMeans(hclust.test)) == 2, ]

g1 = dist(hclust.test)

kit.cn = sc_cn.orig[findOverlaps(
  GRanges(kit.pos),
  GRanges(sc_cn.orig)
) %>%
  as.data.frame() %>%
  .$subjectHits, ]

# Plot as tree
tree1 = hclust(g1)

lab.order = unclass(tree1)$labels[
  unclass(tree1)$order
]

kit.cn = kit.cn[match(lab.order, kit.cn$sample), ]

phylo1 = as.phylo(tree1)

tree.plot = ggtree(phylo1) +
  geom_tiplab(size=3) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = unit(c(0, 2, 0, 0), "lines"))

tree.plot.test = ggplot_build(tree.plot)

label.order = arrange(tree.plot.test$data[[3]], y)$label

kit.cn$sample = factor(kit.cn$sample, levels = label.order)

#KIT copy number barplot
kit.barplot = ggplot(kit.cn, aes(y = sample, x = cn)) +
  geom_col(width = 0.5) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size=8),
    axis.title.x = element_text(size=9)
  ) +
  xlab('KIT copy number')


# Combine into a single plot
sc_plot = plot_grid(tree.plot, kit.barplot, ncol = 2, 
                      axis = 'tb', align = 'h',
                      rel_widths=c(0.65,0.35))





#### Panel e - APG heatmap ####

#### Import data

APG_mat <- read.delim("./dat/APG_mutation_matrix.txt")
APG_anno_vec<-read.delim("./dat/APG_anno_table.txt",header = T)


# Total alterations per sample for bottom anno

APG_df<-APG_mat %>% as.data.frame 
APG_df$Gene <- rownames(APG_df)

APG_long<-APG_df %>% gather(Sample,Alteration,-Gene)

APG_long$mut<-APG_long$Alteration %>% strsplit(.,"-") %>% sapply(.,"[",1)
APG_long$loh<-APG_long$Alteration %>% strsplit(.,"-") %>% sapply(.,"[",2)

APG_long$mut<-ifelse(APG_long$mut=="noMut",0,1)
APG_long$loh<-ifelse(APG_long$loh=="noLOH",0,1)


nAltsPerSample<-APG_long %>% group_by(Sample) %>% 
  summarise(Mut=sum(mut),
            LOH=sum(loh),
            Both=sum(mut==1 & loh==1)) %>%
  as.data.frame


# Total alterations per gene for right anno

nSamplesPerGene<-APG_long %>% group_by(Gene) %>% 
  summarise(Mut=sum(mut),
            LOH=sum(loh),
            Both=sum(mut==1 & loh==1)) %>%
  as.data.frame


# Colors

colors_WGD<-c("#FFFFFF","#9ECAE1","#3182BD")
names(colors_WGD) = c(0, 1, 2)

colors_ploidy = colorRamp2(c(0, 2, 4, 6), c('blue', 'white', 'red', 'yellow'))

colors <- c("#cb181d", "#807dba", "#225ea8", "#DDDDDD")
names(colors) <- c("Mut-LOH", "Mut-noLOH", "noMut-LOH", "noMut-noLOH")


vector_anno_top<-HeatmapAnnotation(
  #Vectors of variables to annotate
  Ploidy=APG_anno_vec$Ploidy,
  WGD=APG_anno_vec$WGD,
  # #List of named vectors for colours
  #Add white lines in between cells
  gp = gpar(col = "white", lwd = 0.01),
  col = list(WGD = colors_WGD, 
             Ploidy = colors_ploidy),
  #Legend font sizes
  annotation_legend_param = list(labels_gp = gpar(fontsize = 8), 
                                 title_gp = gpar(fontsize = 8)),
  show_legend = c(T,TRUE, TRUE))




barplot_anno_bottom <- HeatmapAnnotation(
  n = anno_barplot(
    nAltsPerSample[,c("Both","Mut","LOH")],
    #Assign colors with a named vector
    gp = gpar(fill = structure(c("#807dba","#cb181d","#225ea8"),
                               names=c("Both","Mut","LOH"))),
    bar_width = 0.8
  ),
  annotation_legend_param = list(labels_gp = gpar(fontsize = 7), 
                                 title_gp = gpar(fontsize = 6)),
  #Column annotation
  which="col"
)


barplot_anno_right <- HeatmapAnnotation(
  n = anno_barplot(
    nSamplesPerGene[,c("Both","Mut","LOH")],
    #Assign colors with a named vector
    gp = gpar(fill = structure(c("#807dba","#cb181d","#225ea8"),
                               names = c("Both", "Mut", "LOH"))),
    bar_width = 0.8
  ),
  annotation_legend_param = list(labels_gp = gpar(fontsize = 7), 
                                 title_gp = gpar(fontsize = 6)),
  #Column annotation
  which="row"
)


# Patient vector to split columns and plot in consistent order
patient_vec<- colnames(APG_mat) %>% strsplit(.,"\\.") %>% sapply(.,"[",1)
patient_vec<-factor(patient_vec, levels=unique(patient_vec),ordered=T)

ht_opt$ROW_ANNO_PADDING = unit(0, "cm")
ht_opt$COLUMN_ANNO_PADDING = unit(0, "cm")

loh_heatmap = Heatmap(
  APG_mat,
  na_col = "white",
  col=colors,
  rect_gp = gpar(col = "white", lwd = 0.001),
  #Add annotations
  top_annotation = vector_anno_top,
  right_annotation = barplot_anno_right,
  bottom_annotation = barplot_anno_bottom,
  #Split columns by variable
  column_split = patient_vec,
  #Move row names to the left
  row_names_side = "left",
  #Size of column names
  #column_names_gp = gpar(fontsize = 5),
  show_column_names = F,
  #Size of row names
  row_names_gp = gpar(fontsize = 8),
  #Column titles (in this case, Patient1 etc)
  column_title_gp = gpar(fontsize=9),
  column_title_rot=90,
  #Legend name and sizes
  heatmap_legend_param = list(
    title = " ", 
    legend_direction = "vertical", 
    legend_height = unit(10, "cm"), 
    legend_width = unit(6, "cm"), 
    labels_gp = gpar(fontsize = 7), 
    title_gp = gpar(fontsize = 8)
  )
)


loh_heatmap_grob<-grid.grabExpr(draw(loh_heatmap))

#### Put them together ####

labsize<-14

top_row<-plot_grid(KIT_corrplot,NULL,
                   nrow=1,
                   rel_widths = c(0.55,0.45),
                   labels=c("a","b"),
                   label_size = labsize,
                   scale = c(0.9,1))

mid_row<-plot_grid(sc_plot,NULL,
                   nrow=1,
                   rel_widths = c(0.65,0.35),
                   labels=c("c","d"),
                   label_size = labsize)

bottom_row<-plot_grid(loh_heatmap_grob,
                      nrow=1,
                      rel_widths = c(1),
                      labels=c("e"),
                      label_size = labsize)


plot_grid(top_row,
          mid_row,
          bottom_row,
          ncol=1,
          rel_heights = c(0.33,0.28,0.39))




