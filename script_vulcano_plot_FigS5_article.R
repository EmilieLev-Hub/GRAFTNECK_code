#script vulcano plot fig S5

library(ggplot2)
library(ggrepel)
library(magrittr)
library(dplyr)
library(cowplot)
library(magick)
library(ggimage)
library(tidyverse)


Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs10.08.0/bin/gswin64c.exe")


#path of file
path<-file.path("A/PATH/FOR/FILES")


#file with the results calcuated with med fold change and pvalues over different comparisons
tab_res<-read.csv2(file.path(path,"tab_res_all_comp.csv"),header=TRUE,sep=";")

#########################################################################
#Subpanel A : comparison peritumor are vs tumor in transplant patients
#########################################################################

#create expression.trans for defining coloured points on vulcano plot
tab_res<- tab_res %>%
  mutate(
    Expression.trans = case_when(fold_chang_Trans.log2 >= 0 & raw.pvalues.PA_T.transp <= 0.05 ~ "Up-regulated",
                                 fold_chang_Trans.log2 <= 0 & raw.pvalues.PA_T.transp <= 0.05 ~ "Down-regulated",
                                 TRUE ~ "Unchanged"))


tab_res$Expression.trans[which(is.na(tab_res$fold_chang_Trans.log2)==TRUE | tab_res$raw.pvalues.PA_T.transp=="Nan")]<-NA


top_genes_up.trans<-tab_res[which(tab_res$Expression.trans=="Up-regulated"),] 
top_genes_down.trans<-tab_res[which(tab_res$Expression.trans=="Down-regulated"),]


top_genes_down.trans<-top_genes_down.trans[order(top_genes_down.trans$raw.pvalues.PA_T.transp),]
top_genes_up_down.trans<-rbind(top_genes_up.trans,top_genes_down.trans)
top_genes_up_down.trans<-top_genes_up_down.trans[order(top_genes_up_down.trans$Expression.trans,top_genes_up_down.trans$genes),]



#load logo to visualize comparison
logo <- image_read("logocercle_dbRect.eps")


p2.Trans <- ggplot(tab_res[which(is.na(tab_res$Expression.trans)==FALSE),], aes(fold_chang_Trans.log2, -log(raw.pvalues.PA_T.transp,10))) +
  geom_point(aes(color = Expression.trans)) +geom_hline(yintercept=-log10(0.05), linetype="dashed", size=1)+
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"pvalue")) +
  guides(colour = guide_legend(override.aes = list(size=1.5)))+
  geom_label_repel(data = top_genes_up_down.trans,
                   mapping = aes(fold_chang_Trans.log2, -log(raw.pvalues.PA_T.transp,10), label = genes),max.overlaps = getOption("ggrepel.max.overlaps", default = 15),
                   size = 2)+theme(legend.position = "bottom",legend.title = element_blank(),
                                   plot.background = element_rect(fill = "white"),
                                   panel.background = element_rect(fill = "white"),
                                   axis.line.x = element_line(color = "black"),
                                   axis.line.y = element_line(color = "black")
                   )+
  #red square
  geom_segment(x = 1.0, y = 7.0, xend = 1.7, yend = 7.0, colour = "red",linewidth=1.3)+
  geom_segment(x = 1.0, y = 8.0, xend = 1.7, yend = 8.0, colour = "red",linewidth=1.3)+
  geom_segment(x = 1.0, y = 7.0, xend = 1.0, yend = 8.0, colour = "red",linewidth=1.3)+
  geom_segment(x = 1.7, y = 7.0, xend = 1.7, yend = 8.0, colour = "red",linewidth=1.3)+
  #green square
  geom_segment(x = 0.2, y = 8.3, xend = 1.3, yend = 8.3, colour = "green",linewidth=1.3)+
  geom_segment(x = 0.2, y = 9.5, xend = 1.3, yend = 9.5, colour = "green",linewidth=1.3)+
  geom_segment(x = 0.2, y = 8.3, xend = 0.2, yend = 9.5, colour = "green",linewidth=1.3)+
  geom_segment(x = 1.3, y = 8.3, xend = 1.3, yend = 9.5, colour = "green",linewidth=1.3)+
  
  #yellow square
  geom_segment(x = -0.5, y = 6.3, xend = -1.1, yend = 6.3, colour = "yellow",linewidth=1.3)+
  geom_segment(x = -0.5, y = 7.5, xend = -1.1, yend = 7.5, colour = "yellow",linewidth=1.3)+
  geom_segment(x = -0.5, y = 6.3, xend = -0.5, yend = 7.5, colour = "yellow",linewidth=1.3)+
  geom_segment(x = -1.1, y = 6.3, xend = -1.1, yend = 7.5, colour = "yellow",linewidth=1.3)+
  
  #blue square 
  geom_segment(x = 1, y = 5.3, xend = 2.35, yend = 5.3, colour = "darkblue",linewidth=1.3)+
  geom_segment(x = 1, y = 6.7, xend = 2.35, yend = 6.7, colour = "darkblue",linewidth=1.3)+
  geom_segment(x = 1, y = 5.3, xend = 1, yend = 6.7, colour = "darkblue",linewidth=1.3)+
  geom_segment(x = 2.35, y = 5.3, xend = 2.35, yend = 6.7, colour = "darkblue",linewidth=1.3)+
  
  #purple square
  geom_segment(x = -1.4, y = 7.5, xend = -2.1, yend = 7.5, colour = "purple",linewidth=1.3)+
  geom_segment(x = -1.4, y = 8.5, xend = -2.1, yend = 8.5, colour = "purple",linewidth=1.3)+
  geom_segment(x = -1.4, y = 7.5, xend = -1.4, yend = 8.5, colour = "purple",linewidth=1.3)+
  geom_segment(x = -2.1, y = 7.5, xend = -2.1, yend = 8.5, colour = "purple",linewidth=1.3)+

  #locate logo + defining population 
  annotate(geom="text",x=-3.0,y=7.8,label="Transplant",color="black",fontsize=8)+draw_image(logo, x = -3.5, y = 8.8,scale=2.5)




#########################################################################
#Subpanel A : comparison peritumor are vs tumor in non-transplant patients
#########################################################################
#define Expression variable
tab_res<- tab_res %>%
  mutate(Expression.NT = case_when(fold_chang_NT.log2 >= 0 & raw.pvalues.PA_T.cont <= 0.05 ~ "Up-regulated",
                                   fold_chang_NT.log2 <= 0 & raw.pvalues.PA_T.cont <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged"))

tab_res$Expression.NT[which(is.na(tab_res$fold_chang_NT.log2)==TRUE | tab_res$raw.pvalues.PA_T.contt=="Nan")]<-NA


top_genes_up<-tab_res[which(tab_res$Expression.NT=="Up-regulated"),]
top_genes_down<-tab_res[which(tab_res$Expression.NT=="Down-regulated"),]

top_genes_down<-top_genes_down[order(top_genes_down$raw.pvalues.PA_T.cont),]
top_genes_up_down<-rbind(top_genes_up,top_genes_down)
top_genes_up_down<-top_genes_up_down[order(top_genes_up_down$Expression.NT,top_genes_up_down$genes),]

#vulcano plot
p2.NT <- ggplot(tab_res[which(is.na(tab_res$Expression.NT)==FALSE),], aes(fold_chang_NT.log2, -log(raw.pvalues.PA_T.cont,10))) +
  geom_point(aes(color = Expression.NT)) +  geom_hline(yintercept=-log10(0.05), linetype="dashed", size=1)+
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"pvalue")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_label_repel(data = top_genes_up_down,max.overlaps = getOption("ggrepel.max.overlaps", default = 15),
                   mapping = aes(fold_chang_NT.log2, -log(raw.pvalues.PA_T.cont,10), label = genes),
                   size = 2)+theme(legend.position="none",legend.title = element_blank(),
                                   plot.background = element_rect(fill = "white"),
                                   panel.background = element_rect(fill = "white"),
                                   axis.line.x = element_line(color = "black"),
                                   axis.line.y = element_line(color = "black"))+
  #red square
  geom_segment(x = 1.2, y = 7.5, xend = 1.8, yend = 7.5, colour = "red",linewidth=1.3)+
  geom_segment(x = 1.2, y = 8.6, xend = 1.8, yend = 8.6, colour = "red",linewidth=1.3)+
  geom_segment(x = 1.2, y = 7.5, xend = 1.2, yend = 8.6, colour = "red",linewidth=1.3)+
  geom_segment(x = 1.8, y = 7.5, xend = 1.8, yend = 8.6, colour = "red",linewidth=1.3)+
  #green square
  geom_segment(x = 0.4, y = 9.3, xend = 1.2, yend = 9.3, colour = "green",linewidth=1.3)+
  geom_segment(x = 0.4, y = 11.5, xend = 1.2, yend = 11.5, colour = "green",linewidth=1.3)+
  geom_segment(x = 0.4, y = 9.3, xend = 0.4, yend = 11.5, colour = "green",linewidth=1.3)+
  geom_segment(x = 1.2, y = 9.3, xend = 1.2, yend = 11.5, colour = "green",linewidth=1.3)+
  
  #yellow square
  geom_segment(x = -0.8, y = 8.7, xend = -0.4, yend = 8.7, colour = "yellow",linewidth=1.3)+
  geom_segment(x = -0.8, y = 10, xend = -0.4, yend = 10, colour = "yellow",linewidth=1.3)+
  geom_segment(x = -0.8, y = 8.7, xend = -0.8, yend = 10, colour = "yellow",linewidth=1.3)+
  geom_segment(x = -0.4, y = 8.7, xend = -0.4, yend = 10, colour = "yellow",linewidth=1.3)+
  
  #blue square
  geom_segment(x = 0.85, y = 5.8, xend = 2.62, yend = 5.8, colour = "darkblue",linewidth=1.3)+
  geom_segment(x = 0.85, y = 7.1, xend = 2.62, yend = 7.1, colour = "darkblue",linewidth=1.3)+
  geom_segment(x = 0.85, y = 5.8, xend = 0.85, yend = 7.1, colour = "darkblue",linewidth=1.3)+
  geom_segment(x = 2.62, y = 5.8, xend = 2.62, yend = 7.1, colour = "darkblue",linewidth=1.3)+
  
  #purple square
  geom_segment(x = -1.35, y = 10.6, xend = -1.9, yend = 10.6, colour = "purple",linewidth=1.3)+
  geom_segment(x = -1.35, y = 11.5, xend = -1.9, yend = 11.5, colour = "purple",linewidth=1.3)+
  geom_segment(x = -1.35, y = 10.6, xend = -1.35, yend = 11.5, colour = "purple",linewidth=1.3)+
  geom_segment(x = -1.9, y = 10.6, xend = -1.9, yend = 11.5, colour = "purple",linewidth=1.3)+
  
  #load logo + legend
  annotate(geom="text",x=2.2,y=9.5,label="Non-Transplant",color="black",fontsize=8)+draw_image(logo, x = 1.68, y = 10.4,scale=2.5)


#create one figure with two subpanels

prow <- plot_grid( p2.Trans + theme(legend.position="none"),
                   p2.NT + theme(legend.position="none"),
                   align = 'vh',
                   labels = c("   (A)  ", "    (B)  "),
                   hjust = -1,
                   nrow = 1
)

# extract the legend from one of the plots
legend_b = cowplot::get_plot_component(p2.Trans, 'guide-box-bottom', return_all = TRUE)

#grid subpanel+ common legend
fig.S5.comLeg <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .1))

#save figure
ggsave(fig.S5.comLeg,file=file.path(path,"FigS5.pdf"),width=16,height=7,units=("in"),dpi=600)
