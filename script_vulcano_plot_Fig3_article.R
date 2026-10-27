#script for vulcano plot 
#Figure 3

library(ggplot2)
library(ggrepel)
library(magrittr)
library(dplyr)
library(cowplot)
library(magick)
library(ggimage)
library(tidyverse)


Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs10.08.0/bin/gswin64c.exe")


#path of files
path<-file.path("A/PATH/FOR/FILES")


#file with the results calcuated with med fold change and pvalues over different comparisons
tab_res<-read.csv2(file.path(path,"tab_res_all_comp.csv"),header=TRUE,sep=";")


###########################################################################
#subpanel A) Peritumor area comparison over transp and non-transp patients
############################################################################



#define "Expression.PA" variable as up-reg and down-reg genes on peritumor area according to a positive/negative med fold change + pvalues below 0.05
# less than a pval >0.05  (lower than 1.30 (log10 of 0.05)) : consider as unchanged
#in order to visualize on plot the difference btw up or down regulated genes  
tab_res<- tab_res %>% 
  mutate(
    Expression.PA = case_when(fold_chang_PA.log2 >= 0 & raw.pvalues.PA <= 0.05 ~ "Up-regulated",
                               fold_chang_PA.log2 <= 0 & raw.pvalues.PA <= 0.05 ~ "Down-regulated",
                               TRUE ~ "Unchanged")
  )

tab_res$Expression.PA[which(is.na(tab_res$fold_chang_PA.log2)==TRUE | tab_recap$raw.pvalues.PA=="NaN")]<-NA

#list of up reg genes and down reg genes
top_genes_up.PA<-tab_res[which(tab_res$Expression.PA=="Up-regulated"),] 
top_genes_down.PA<-tab_res[which(tab_res$Expression.PA=="Down-regulated"),]

#keep only up and down , range by alphabetic order
#it will be used for naming points in vulcanoplot
top_genes_up_down.PA<-rbind(top_genes_up.PA,top_genes_down.PA)
top_genes_up_down.PA<-top_genes_up_down.PA[order(top_genes_up_down.PA$Expression.PA,top_genes_up_down.PA$genes),]




#load image logo created via ggplot and save as .eps to keep quality image in vulcanoplot to understand which comparison was plotted
#PA here
logo.ext <- image_read(file.path(path,"logocercle_ext.eps"))

#to define log pvalue in log10
tab_res$raw.pvalues.PA.log10<- -log(tab_res$raw.pvalues.PA,10)


#to plot vulcanoplot which is defined by x-axis :log2 of Fold of Change  and y-axis: -log10 of pval
#only over genes which have been catgorized as up/down or unchanged 
p2.PA <- ggplot(tab_res[which(is.na(tab_res$Expression.PA)==FALSE),], aes(fold_chang_PA.log2, raw.pvalues.PA.log10))+
  geom_point(aes(color = Expression.PA)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", size=1)+
  xlim(-3.5,1.25) + ylim(0,4)+ xlab(expression("log"[2]*"FC")) +   ylab(expression("-log"[10]*"pvalue")) +
  geom_label_repel(data = top_genes_up_down.PA,
                   mapping = aes(fold_chang_PA.log2, -log(raw.pvalues.PA,10), label = genes),
                   size = 2)+
  theme(legend.position = "none",plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")
  )+
  #red square
  geom_segment(x = -1, y = 1.2, xend = -0.25, yend = 1.2, colour = "red",linewidth=1.3)+
  geom_segment(x = -1, y = 3, xend = -0.25, yend = 3, colour = "red",linewidth=1.3)+
  geom_segment(x = -1, y = 1.2, xend = -1, yend = 3, colour = "red",linewidth=1.3)+
  geom_segment(x = -0.25, y = 1.2, xend = -0.25, yend = 3, colour = "red",linewidth=1.3)+
  #green square
  geom_segment(x = -0.2, y = 1.2, xend = 0.5, yend = 1.2, colour = "green",linewidth=1.3)+
  geom_segment(x = -0.2, y = 1.8, xend = 0.5, yend = 1.8, colour = "green",linewidth=1.3)+
  geom_segment(x = -0.2, y = 1.2, xend = -0.2, yend = 1.8, colour = "green",linewidth=1.3)+
  geom_segment(x = 0.5, y = 1.2, xend = 0.5, yend = 1.8, colour = "green",linewidth=1.3)+
  
  #yellow square
  geom_segment(x = 0, y = 2.65, xend = 1, yend = 2.65, colour = "yellow",linewidth=1.3)+
  geom_segment(x = 0, y = 3.3, xend = 1, yend = 3.3, colour = "yellow",linewidth=1.3)+
  geom_segment(x = 0, y = 2.65, xend = 0, yend = 3.3, colour = "yellow",linewidth=1.3)+
  geom_segment(x = 1, y = 2.65, xend = 1, yend = 3.3, colour = "yellow",linewidth=1.3)+

  #to visualize logo in order to understand which is compared btw transplant and non-transplant : peritumor area 
  annotate(geom="text",x=-3,y=2.8,label="Transplant",color="black",fontsize=8)+draw_image(logo.ext, x = -3.5, y = 3,scale=2.0)+
  annotate(geom="text",x=-2,y=2.8,label="Non-Transplant",color="black",fontsize=8)+draw_image(logo.ext, x = -2.6, y = 3,scale=2.0)



####################################################################
#subpanel B) Tumor comparison over transp and non-transp patients
####################################################################

#load image logo created via ggplot and save as .eps to keep quality image in vulcanoplot to understand which comparison was plotted
#Tumor here
logo.int <- image_read(file.path(path,"logocercle_int.eps"))



tab_res<- tab_res %>% 
  mutate(
    Expression.T = case_when(fold_chang_T.log2 >= 0 & raw.pvalues.T <= 0.05 ~ "Up-regulated",
                               fold_chang_T.log2 <= 0 & raw.pvalues.T <= 0.05 ~ "Down-regulated",
                               TRUE ~ "Unchanged")
  )

tab_res$Expression.T[which(is.na(tab_res$fold_chang_T.log2)==TRUE | tab_res$raw.pvalues.T=="NaN")]<-NA

top_genes_up.tum<-tab_res[which(tab_res$Expression.T=="Up-regulated"),] 
top_genes_down.tum<-tab_res[which(tab_res$Expression.T=="Down-regulated"),] 

top_genes_up_down.tum<-rbind(top_genes_up.tum,top_genes_down.tum)
top_genes_up_down.tum<-top_genes_up_down.tum[order(top_genes_up_down.tum$Expression.T,top_genes_up_down.tum$genes),]


tab_res$raw.pvalues.T.log10<- -log(tab_res$raw.pvalues.T,10)
p2.T <- ggplot(tab_res[which(is.na(tab_res$Expression.T)==FALSE),], aes(fold_chang_T.log2,raw.pvalues.T.log10 )) +
  geom_point(aes(color = Expression.T)) +geom_hline(yintercept=-log10(0.05), linetype="dashed", size=1)+
  ylim(0,4)+scale_x_continuous(breaks=seq(-2,1.5,by=0.5))+
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"pvalue")) +
 guides(colour = guide_legend(override.aes = list(size=1.5))) +  
  geom_label_repel(data = top_genes_up_down.tum,
                   mapping = aes(fold_chang_T.log2, -log(raw.pvalues.T,10), label = genes),
                   size = 2)+theme(legend.position = "bottom",legend.title = element_blank(),
                                   plot.background = element_rect(fill = "white"),
                                   panel.background = element_rect(fill = "white"),
                                   axis.line.x = element_line(color = "black"),
                                   axis.line.y = element_line(color = "black")
                   )+
  #red square
  geom_segment(x = -1, y = 1.6, xend = -0.17, yend = 1.6, colour = "red",linewidth=1.3)+
  geom_segment(x = -1, y = 3.95, xend = -0.17, yend = 3.95, colour = "red",linewidth=1.3)+
  geom_segment(x = -1, y = 1.6, xend = -1, yend = 3.95, colour = "red",linewidth=1.3)+
  geom_segment(x = -0.17, y = 1.6, xend = -0.17, yend = 3.95, colour = "red",linewidth=1.3)+
  #green square
  geom_segment(x = 0.5, y = 1.77, xend = 0.7, yend = 1.77, colour = "green",linewidth=1.3)+
  geom_segment(x = 0.5, y = 2, xend = 0.7, yend = 2, colour = "green",linewidth=1.3)+
  geom_segment(x = 0.5, y = 1.77, xend = 0.5, yend = 2, colour = "green",linewidth=1.3)+
  geom_segment(x = 0.7, y = 1.77, xend = 0.7, yend = 2, colour = "green",linewidth=1.3)+

  #to visualize logo in order to understand which is compared btw transplant and non-transplant : tumor 
  annotate(geom="text",x=0.7,y=2.8,label="Transplant",color="black",fontsize=8)+draw_image(logo.int, x = 0.2, y = 3,scale=1.2)+
  annotate(geom="text",x=1.2,y=2.8,label="Non-Transplant",color="black",fontsize=8)+draw_image(logo.int, x = 0.7, y = 3,scale=1.2)




###############################
#creation over one figure with A) and B) subpanels
###############################


prow <- plot_grid( p2.PA + theme(legend.position="none"),
                   p2.T + theme(legend.position="none"),
                   align = 'vh',
                   labels = c("  (A)  ", "  (B)  "),
                   hjust = -1,
                   nrow = 1
)

# extract the legend from one of the plots 
legend_b = cowplot::get_plot_component(p2.T, 'guide-box-bottom', return_all = TRUE)

#associated 2 subpanels with one legend
fig.3.comLeg <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .1))

#ggsave as pdf
ggsave(fig.3.comLeg,file=file.path(path2,"Fig3.pdf"),width=16,height=7,units=("in"),dpi=600)
