options(stringsAsFactors = F)

library(ggplot2)
library(ggpubr)
library(tidyverse)
library(RColorBrewer)
library(cowplot)


#### Panel a - wgii site randomization ####


extract.legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


wgii_df <- read.delim("./dat/brain_wgii_all.txt",header=T)
site_colors <- read.delim("./dat/brain_site_colors.txt",header = T)


#list with wgii per patient

patients<-wgii_df %>% distinct(Patient) %>% pull %>% sort

wgii_list<-lapply(patients, function(PATIENT) 
  wgii_df %>% filter(Patient==PATIENT) %>% 
    select(wGII) %>% pull)

names(wgii_list)<-patients


#randomise site

sites<-wgii_df %>% distinct(Site) %>% pull %>% sort

n_iter=10000
means<-list()

for (SITE in sites){
  means[[SITE]]<-list()
  # Number of site samples per patient
  p_freq<-as.data.frame(table(wgii_df$Patient[wgii_df$Site==SITE]))
  for (i in 1:n_iter){
    rep_sample<-list()
    # One iteration per patient with samples from that site
    for (j in 1:nrow(p_freq)){
      PATIENT<-as.character(p_freq$Var1[j])
      # Sample wgii values from any site, sample size equal to the site n
      rep_sample[[j]]<-sample(wgii_list[[PATIENT]], p_freq$Freq[j])
    }
    # Combine wgiis from all patients
    rep_sample<-unlist(rep_sample)
    # Calculate mean 
    means[[SITE]][[i]]<-mean(rep_sample)
  }
  # Means from all iterations into a single column
  means[[SITE]]<-cbind.data.frame(Site=SITE, wGII=unlist(means[[SITE]]))
}

wgii_df_null<-do.call(rbind.data.frame,means)

# Tag as belonging to the null distribution
wgii_df_null$Group<-"null"

# Same columns for observed wgii
wgii_df_obs<-wgii_df %>% select(Site,wGII)
wgii_df_obs$Group<-"observed"

# Combine into single df
wgii_null<-rbind.data.frame(wgii_df_obs,
                            wgii_df_null)



## Vars to factor to plot
wgii_null$Site<-factor(wgii_null$Site,
                       levels=c("Brain","Leptom","Lung","Liver","Adrenal",
                                "Peritoneal","Lymph node"),
                       ordered = T)

wgii_null$Group<-factor(wgii_null$Group,
                        levels = sort(unique(wgii_null$Group),decreasing = T),
                        ordered = T)



# Calculate p values

pvalues<-sapply(levels(wgii_null$Site), function(SITE)
  if(wgii_df_obs %>% filter(Site==SITE) %>% summarise(mean(wGII)) %>% pull < 
     wgii_df_null %>% filter(Site==SITE) %>% summarise(mean(wGII)) %>% pull){
    sum(wgii_df_null %>% filter(Site==SITE) %>% pull(wGII) <
          wgii_df_obs %>% filter(Site==SITE) %>% summarise(mean(wGII)) %>% pull)/n_iter
  }
  else{
    sum(wgii_df_null %>% filter(Site==SITE) %>% pull(wGII) >
          wgii_df_obs %>% filter(Site==SITE) %>% summarise(mean(wGII)) %>% pull)/n_iter
  })

pval_df<-cbind.data.frame(Site=levels(wgii_null$Site),
                        Pval=pvalues)

#Add bars, remove pvals and only have **

pvals_ast<-pval_df[pval_df$Pval<=0.05,]
pvals_ast$asterisks<-ifelse(pvals_ast$Pval<=0.001,"***",
                            ifelse(pvals_ast$Pval<=0.01,"**",
                                   ifelse(pvals_ast$Pval<=0.05,"*","")))



wGII_site<-ggplot(wgii_null,aes(Site,wGII))+
  geom_boxplot(aes(fill=Site,alpha=Group))+
  scale_fill_manual(values=site_colors$HEX,name="Site")+
  scale_alpha_manual(values=c(0.8,0))+
  xlab("")+
  ylab("wGII")+
  ylim(0,1)+
  theme_bw(base_size =12)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,size=10),
        axis.text.y = element_text(size=8),
        axis.title.y = element_text(size=10))+
  #add pvals
  geom_text(data=pvals_ast,aes(Site,0.92,label=asterisks))+
  geom_segment(x=0.8,xend=1.2,y=0.9,yend=0.9)+
  geom_segment(x=3.8,xend=4.2,y=0.9,yend=0.9)+
  scale_color_manual(values=c("#969696","#cb181d"))


LP<-ggplot(wgii_null,aes(Site,wGII))+
  geom_boxplot(aes(fill="Site-specific null distribution"))+
  scale_fill_manual(values=c("white"),name="")+
  theme_bw(base_size = 12)+
  theme(legend.box.spacing = unit(0, "pt"),
        legend.position="bottom",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))



wGII_site_legend<-extract.legend(LP)




#### Panel b - CN dist vs brain latency ####

brain_dist_time_df <- read.delim("./dat/brain_dist_time.txt",header=T)
subtype_cols_df<-read.table("./dat/subtype_colours.txt",quote="\"", comment.char="")

subtype_cols<-subtype_cols_df$V2
names(subtype_cols)<-subtype_cols_df$V1

brain_dist_time_df$Type[brain_dist_time_df$Type=="Leptomeningeal"]<-"Leptomen."

# Arrange geom_text positions

# y axis
brain_dist_time_df$ylabpos<-brain_dist_time_df$Days_IV_to_brain+60

bellow_lm<-which(brain_dist_time_df$normDist>0 &
                   brain_dist_time_df$Days_IV_to_brain<500)

brain_dist_time_df$ylabpos[bellow_lm]<-brain_dist_time_df$Days_IV_to_brain[bellow_lm]-60

# x axis
brain_dist_time_df$xlabpos<-brain_dist_time_df$normDist

move_right<-which(brain_dist_time_df$normDist<0 |
                    brain_dist_time_df$Subtype=="Acral")

brain_dist_time_df$xlabpos[move_right]<-brain_dist_time_df$normDist[move_right]+0.35

move_left<-which(brain_dist_time_df$normDist>1.15)

brain_dist_time_df$xlabpos[move_left]<-brain_dist_time_df$normDist[move_left]-0.35


# plot

brain_dist_vs_time<-ggplot(brain_dist_time_df,
             aes(normDist,Days_IV_to_brain),alpha=0.99)+
  geom_smooth(method='lm',se=F,color="black")+
  geom_point(aes(color=Subtype, shape=Type),size=5)+
  scale_color_manual(values=subtype_cols)+
  xlab("Normalised CN distance BR to other sites")+
  ylab("Time from stage IV to brain metastasis")+
  geom_text(aes(xlabpos,ylabpos,label=Patient),size=2.5)+
  theme_bw(base_size = 9)+
  stat_cor(method = "pearson", label.x = -1.8, label.y = 925)+
  theme(legend.text = element_text(size=7),
        legend.key.height= unit(2, 'mm'),
        legend.key.width= unit(4, 'mm'),
        legend.box.spacing = unit(0, "pt"),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8))+
  guides(color = guide_legend(override.aes = list(size = 3.5)),
         shape = guide_legend(override.aes = list(size = 3.5)))






#### Panel c - CRUKP5107 radio dynamics #####


radio_info <- read.delim("./dat/brain_radio_sizes.txt",header=T)
treatment_info <- read.delim("./dat/brain_treatment_info.txt",header=T)

treatment_colors<-c("#80cdc1","#abd9e9","#355f8d","#EBC3DB")

lesion_colors<-colorRampPalette(brewer.pal(9,"Greys"))(12)
lesion_colors<-c(lesion_colors[5:12],"#fc8d59")


radio_info$Site<-factor(radio_info$Site,
                        levels=sort(unique(radio_info$Site))[c(2:9,1)],
                        ordered = T)

#Make treatments factor so colors are shown in order in legend
treatment_info$Treatment<-factor(treatment_info$Treatment,
                                 levels = unique(treatment_info$Treatment),
                                 ordered = T)


radio_plot<-ggplot(radio_info)+
  # geom_point(aes(MonthsElapsed,Size,color=Site),size=2)+
  geom_line(aes(MonthsElapsed,Size,color=Site),size=1.2)+
  scale_color_manual(values=lesion_colors)+
  ylab("Size (mm)")+
  geom_rect(data=treatment_info,aes(xmin=Start_inMonthsSinceDx,xmax=End_inMonthsSinceDx,ymin=-15,ymax=-5,fill=Treatment))+
  geom_text(data=treatment_info,aes(x=Start_inMonthsSinceDx+(End_inMonthsSinceDx-Start_inMonthsSinceDx)/2,y=-10,label=Best.response),
            size=2.5,fontface = "bold")+
  scale_fill_manual(values=treatment_colors)+
  scale_x_continuous(breaks = seq(0,round(max(treatment_info$End_inMonthsSinceDx,na.rm = T)),6))+
  ylab("Lesion size (mm)")+
  xlab("Time from diagnosis of stage IV disease (months)")+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=10),
        axis.text.x = element_text(size=8),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        legend.key.height= unit(3, 'mm'),
        legend.key.width= unit(3, 'mm'),
        legend.box.spacing = unit(4, "pt"))





#### Load panels d and e and edit ####

phylogeny<-readRDS("./dat/brain_phylogeny.RDS")
SBS_phylo<-readRDS("./dat/brain_sbs_phylogeny.RDS")

#Solve incompatibility between R versions
SBS_phylo$theme<-ggplot2::theme_void()

phylogeny$theme<-ggplot2::theme_bw()+theme(legend.position = "none",
                                           panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(),
                                           panel.border = element_blank(),
                                           axis.text.y = element_blank(),
                                           axis.ticks.y = element_blank(),
                                           axis.line.x = element_line())

phylogeny$layers[[3]]$aes_params$alpha=0.99
phylogeny$layers[[3]]$aes_params$size=4

SBS_phylo<-SBS_phylo+theme(legend.text = element_text(size=7),
                legend.title = element_text(size=0),
                legend.key.height= unit(3.5, 'mm'),
                legend.key.width= unit(3.5, 'mm'),
                legend.box.spacing = unit(4, "pt"),
                legend.position = c(0.25,0.6))




#### Put together ####

labsize<-14


wGII_site_shorter<-plot_grid(wGII_site+theme(plot.margin=grid::unit(c(2,2,0,2), "mm")),
                             wGII_site_legend,
                             NULL,
                             ncol=1,rel_heights = c(0.9,0.05,0.05))


distTime_shorter<-plot_grid(brain_dist_vs_time,NULL,
                            ncol=1,
                            rel_heights = c(1,0.15))

top_row<-plot_grid(wGII_site_shorter,NULL,distTime_shorter,
          nrow=1,
          rel_widths = c(0.47,0.02,0.5),
          labels=c("a","b",""),
          label_size = labsize)


####

mid_row<-plot_grid(radio_plot,
                   nrow=1,
                   rel_widths = c(1),
                   labels=c("c"),
                   label_size = labsize)


### 


phylo_panel<-plot_grid(phylogeny,
                       SBS_phylo,
                       nrow=2,
                       labels=c("d","e"),
                       label_size = labsize)


bottom_row<-plot_grid(phylo_panel,NULL,
                      nrow=1,
                      rel_widths = c(0.6,0.4),
                      labels=c("","f"),
                      label_size = labsize)



plot_grid(top_row,
          mid_row,
          bottom_row,
          ncol=1,
          rel_heights = c(0.8,0.6,1))




