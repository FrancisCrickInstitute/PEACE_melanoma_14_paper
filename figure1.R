options(stringsAsFactors = F)

#panel a
library(ComplexHeatmap)
#panel b
library(ggplot2)
library(ggrepel)
library(tidyverse)
#panel c
library(RColorBrewer)
#combining
library(cowplot)


#### Panel a - oncoplot ####


### READ DATA ###

# Read oncoplot info 
oncoplot_main_mat <- read.delim("./dat/oncoplot_main_mat.txt")
top_annotation <- read.delim("./dat/oncoplot_top_anno.txt")
# Read colour files
colors_list<-readRDS("./dat/oncoplot_anno_colours.RDS")
var_colors<-readRDS("./dat/oncoplot_variant_colours.RDS")


### ORGANISE DATA ###

# Metastatic sites as ordered factor
top_annotation$Site<-factor(top_annotation$Site,
                            levels=c("PR","LI","LU",
                                     "LN","ST","BR",
                                     "LMS","PE","AD",
                                     "Other"), ordered = T)

# Patients as ordered factor
top_annotation$Patient<-factor(top_annotation$Patient, 
                               levels=unique(top_annotation$Patient),ordered=T)


### PLOT ###

# Set up ComplexHeatmap global vars to control legend width
ht_opt(legend_grid_height=unit(3, "mm"),
       legend_grid_width=unit(3, "mm"))

# Make top annotation object
anno_top<-HeatmapAnnotation(Subtype=top_annotation$Subtype,
                            TMB = anno_barplot(top_annotation$TMB,bar_width = 0.8,height = unit(5, "mm")),
                            Ploidy=anno_barplot(top_annotation$Ploidy,bar_width=0.8,height = unit(5, "mm")),
                            WGD=as.factor(top_annotation$WGD),
                            wGII=top_annotation$wGII,
                            Site=top_annotation$Site,
                            col = colors_list,
                            annotation_legend_param = list(labels_gp = gpar(fontsize = 6),
                                                           title_gp = gpar(fontsize = 7),
                                                           Subtype=list(nrow=1),
                                                           WGD=list(nrow=1),
                                                           wGII=list(direction = "horizontal"),
                                                           Site=list(nrow=1)),
                            which="col",
                            annotation_name_gp= gpar(fontsize = 7),
                            simple_anno_size = unit(5,"mm"),height = unit(28,"mm"))



# Make oncoplot heatmap with top annotation
heatmap<-Heatmap(as.matrix(oncoplot_main_mat),
                 name="Variant type",
                 na_col = "white",
                 col=var_colors,
                 #Add annotations
                 top_annotation = anno_top,
                 column_split=top_annotation$Patient,
                 #Move row names to the left
                 row_names_side = "right",
                 #Size of column names
                 #column_names_gp = gpar(fontsize = 3),
                 # Decided to remove sample names
                 show_column_names = F,
                 #Size of row names
                 row_names_gp = gpar(fontsize = 6),
                 #Column titles (in this case, Patient1 etc)
                 column_title_gp = gpar(fontsize=6.5),
                 column_title_rot=90,
                 #Legend name and sizes
                 heatmap_legend_param = list(
                                             labels_gp = gpar(fontsize = 6),
                                             title_gp = gpar(fontsize = 7)),
                 height = nrow(oncoplot_main_mat)*unit(2, "mm"),
                 show_heatmap_legend = FALSE)


# Heatmap legend as separate object to plot in single column
lgd = Legend(labels = names(var_colors), title = "Variant type", legend_gp = gpar(fill = var_colors),
             grid_height = unit(3, "mm"), grid_width = unit(3, "mm"),labels_gp = gpar(fontsize = 6),title_gp =gpar(fontsize = 7) )

# Draw heatmap and legend on top
heatmap
draw(lgd,x = unit(0.87, "npc"), y = unit(0.87, "npc"))


# Capture drawing for final composition
heatmap_grob<-grid.grab()


#### Panel b - CN freq per sample ####


### READ DATA ###

# Read CN call proportions data
CNcall_prop<-read.delim("./dat/CNcall_prop.txt")


### ORGANISE DATA ###

# Create palette for existing CN values (to plot as factor) - range from -6 to 6
maxCNcall<-6
palette<-colorRampPalette(brewer.pal(11,"RdBu"))(maxCNcall*2+1)
names(palette)<-seq(maxCNcall,-maxCNcall)
pal<-palette[which(names(palette)%in%as.numeric(unique(CNcall_prop$CN_call)))]

# CN calls as ordered factor
CN_call_levels<-CNcall_prop$CN_call %>% unique %>% as.numeric %>% sort %>% rev

CNcall_prop$CN_call<-factor(CNcall_prop$CN_call,
                         levels=CN_call_levels,
                         ordered=T)

# Samples as ordered factor
sample_levels<-unique(CNcall_prop$sample)
CNcall_prop$sample<-factor(CNcall_prop$sample,
                           levels=sample_levels,
                           ordered=T)


### PLOT ###

# Plot CN call proportions as barplot
CNcalls_WGD_plot<-ggplot(CNcall_prop)+
  geom_bar(aes(sample,Prop,fill=CN_call),stat="identity")+
  scale_fill_manual(values=pal,name="Copy number call")+
  theme_bw()+
  ylab("Proportion of the genome")+xlab("")+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))+
  facet_grid(.~WGD,scales="free_x",space="free")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=7),
        axis.title.y = element_text(size=8),
        legend.title = element_text(size=8),
        legend.text = element_text(size=7),
        legend.key.size = unit(0.35, 'cm'))+
  guides(fill=guide_legend(ncol=2))

# Give padding around the barplot so it fits better in the composition
CN_states<-plot_grid(NULL,CNcalls_WGD_plot,NULL,nrow=1,rel_widths = c(0.028,1,0.02))


#### Panel c - SCNA freq ####


### READ DATA ###

# Read CN rate data and TCGA GISTIC annotation info
CN_rates<-read.delim("./dat/CNfreq_rates.txt")
TCGA_GISTIC<-read.delim("./dat/CNfreq_anno_TCGA_GISTIC.txt")


### PLOTTING FUNCTIONS ###

# Plotting function: SCNA freq with LOH and AI rate bellow
plot_SCNA_LOH_AI_freq<-function(CN_rates_singleCoord, CN_rates, title){
  YMAX<-(-1.6)
  plot<-ggplot(CN_rates_singleCoord)+
    geom_linerange(aes(x=coord,ymin=0,ymax=gain.subclonal),color="#fc9272",alpha=0.5)+
    geom_linerange(aes(x=coord,ymin=-loss.subclonal,ymax=0),color="#6baed6",alpha=0.5)+
    #Plot clonal on top
    geom_linerange(aes(x=coord,ymin=0,ymax=gain.clonal),color="#d73027",alpha=0.5)+
    geom_linerange(aes(x=coord,ymin=-loss.clonal,ymax=0),color="#08519c",alpha=0.5)+
    geom_line(aes(x=coord,y=gain.overall),color="black")+
    geom_line(aes(x=coord,y=-loss.overall),color="black")+
    geom_hline(yintercept=0)+
    geom_rect(data=CN_rates,
              aes(xmin=start,xmax=end,ymin=YMAX-0.4,ymax=YMAX-0.2,fill=LOH.clonal))+
    geom_rect(data=CN_rates,
              aes(xmin=start,xmax=end,ymin=YMAX-0.6,ymax=YMAX-0.4,fill=LOH.subclonal))+
    geom_rect(data=CN_rates,
              aes(xmin=start,xmax=end,ymin=YMAX-0.9,ymax=YMAX-0.7,fill=AI.clonal))+
    geom_rect(data=CN_rates,
              aes(xmin=start,xmax=end,ymin=YMAX-1.1,ymax=YMAX-0.9,fill=AI.subclonal))+
    scale_fill_gradient(low="white",high="#016c59",name="Rate")+
    ylab("")+
    facet_grid(.~chr, space = "free",scales = "free_x")+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.spacing = unit(0, "lines"),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.background= element_blank(),
          strip.background = element_rect(fill = "white",colour = NA),
          panel.border = element_rect(color = "lightgrey", fill = NA, size = 0.3))+
    ggtitle(title)
  return(plot)
}

# Function to add TCGA GISTIC annotations to SCNA freq plot
addTCGA2<-function(plot){
  YMAX<-(-1.6)
  new_plot<-plot+
    geom_rect(data=TCGA_GISTIC[TCGA_GISTIC$Type%in%c("driverGain"),],
              aes(xmin=start,xmax=end,ymin=1.05,ymax=1.1),fill="#9B111E",color="#9B111E")+
    geom_rect(data=TCGA_GISTIC[TCGA_GISTIC$Type%in%c("driverLoss"),],
              aes(xmin=start,xmax=end,ymin=-1.1,ymax=-1.05),fill="#1034A6",color="#1034A6")+
    geom_rect(data=TCGA_GISTIC[TCGA_GISTIC$Type%in%c("driverGeneCNGain"),],
              aes(xmin=start,xmax=end,ymin=1.15,ymax=1.2),fill="#E60026",color="#E60026")+
    geom_rect(data=TCGA_GISTIC[TCGA_GISTIC$Type%in%c("driverGeneCNLoss"),],
              aes(xmin=start,xmax=end,ymin=-1.15,ymax=-1.2),fill="#007FFF",color="#007FFF")+
    geom_text_repel(data=TCGA_GISTIC[TCGA_GISTIC$Type%in%c("driverGeneCNGain"),],
                    aes(x=start,y=1.25,label=event),size=1.5,
                    #force_pull   = 0, # do not pull toward data points
                    # nudge_y      = 0.05,
                    direction    = "x",
                    angle        = 70,
                    hjust        = 0,
                    vjust= 0.7,
                    segment.color="lightgrey",
                    min.segment.length = 0.1,
                    # segment.size = 0.2,
                    segment.inflect = T)+
    geom_text_repel(data=TCGA_GISTIC[TCGA_GISTIC$Type%in%c("driverGeneCNLoss"),],
                    aes(x=start,y=-1.2,label=event),size=1.5,
                    # force_pull   = 0, # do not pull toward data points
                    nudge_y      = -0.15,
                    direction    = "x",
                    angle        = 110,
                    hjust        = 0,
                    vjust= -0.3,
                    segment.color="lightgrey",
                    min.segment.length = 0.1,
                    # segment.size = 0.2,
                    segment.inflect = T)+
    
    scale_y_continuous(limits=c(-2.8,1.6),breaks=c(YMAX-1,YMAX-0.8,YMAX-0.5,YMAX-0.3,
                                                   -1,-0.5,0.5,1),
                       labels=c("AI subclonal","AI clonal","LOH subclonal","LOH clonal",
                                "Loss 100%","Loss 50%","Gain 50%","Gain 100%"))
  
  
  
  return(new_plot)
}

# Function to extract legend
extract.legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}



### ORGANISE DATA ### 

#Expand rate df into single coordinate per row
start_rates<-CN_rates %>% select(-end) %>% dplyr::rename(coord=start) %>% mutate(type="start")
end_rates<-CN_rates %>% select(-start) %>% dplyr::rename(coord=end) %>% mutate(type="end")
CN_rates_singleCoord<-rbind.data.frame(start_rates,end_rates)


### PLOT ### 

# SCNA freq plot
CNfreq_plot<-plot_SCNA_LOH_AI_freq(CN_rates_singleCoord,CN_rates,"")+
  theme(axis.text.y = element_text(size=6),
        strip.text.x = element_text(size = 8))

# Remove legend
CNfreq_plot_noleg<-CNfreq_plot+theme(legend.position = "none")

# Annotate SCNA freq plot with TCGA GISTIC events
CNfreq_plot_anno<-addTCGA2(CNfreq_plot_noleg)


### LEGENDS ### 


## option A: create and extract ##

# Extract LOH and AI legend
LOH_AI_legend<-extract.legend(CNfreq_plot+
                                theme(legend.title=element_text(size=8),
                                      legend.text = element_text(size = 7),
                                      legend.key.size = unit(0.4, "cm")))

# Create dummy df fot TCGA legend
TCGA_GISTIC_forlegend<-TCGA_GISTIC

TCGA_GISTIC_forlegend$Type<-factor(TCGA_GISTIC_forlegend$Type,
                                   levels=unique(TCGA_GISTIC_forlegend$Type)[c(1,3,4,2)],
                                   ordered=T)

levels(TCGA_GISTIC_forlegend$Type)<-c("Driver gain","GISTIC gain",
                                      "GISTIC loss","Driver loss")

# Extract TCGA GISTIC legend
TCGA_legend<-extract.legend(
  ggplot(TCGA_GISTIC_forlegend)+geom_rect(aes(xmin=chr,xmax=chr,ymax=1,ymin=1,fill=Type))+
    scale_fill_manual(values = c("#E60026","#9B111E","#1034A6","#007FFF"),name="Annotations")+
    theme(legend.title=element_text(size=8),
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.5, "cm")))


#Create a dummy df for SCNA freq legend
rate_DF_forLegend<-cbind.data.frame(Rate=c("Clonal gain",
                                           "Subclonal gain",
                                           "Clonal loss",
                                           "Subclonal loss",
                                           "Any event"))

rate_DF_forLegend$Rate<-factor(rate_DF_forLegend$Rate,levels=rate_DF_forLegend$Rate,ordered = T)

# Extract SCNA freq legend
CNfreq_legend<-extract.legend(ggplot(rate_DF_forLegend)+
                                geom_linerange(aes(x=1,ymin=0,ymax=0,color=Rate),size=1.5)+
                                scale_color_manual(values=c("#d73027","#fc9272","#08519c","#6baed6","black"))+
                                theme_bw()+
                                theme(legend.title=element_blank(),
                                      legend.text = element_text(size = 8),
                                      legend.key.size = unit(0.4, "cm")))


### option B: load from RDS ##

# CNfreq_legend<-readRDS("./dat/CNfreq_freq_legend.RDS")
# TCGA_legend<-readRDS("./dat/CNfreq_TCGA_legend.RDS")
# LOH_AI_legend<-readRDS("./dat/CNfreq_LOH_AI_legend.RDS")



### PLOT with legends ### 

# Compose panel c
CNfreq_withLegend<-plot_grid(CNfreq_plot_anno,
                             plot_grid(CNfreq_legend,
                                       TCGA_legend,
                                       LOH_AI_legend,ncol=1),
                             nrow=1,
                             rel_widths = c(1,0.1))




#### Combine panels into fig1 ####


plot_grid(heatmap_grob,
          CN_states,
          CNfreq_withLegend,
          ncol=1,
          rel_heights = c(0.8,0.31,0.69),
          labels = c("a","b","c"),
          label_size = 13)




