#script IHC for figure S2
#explain the database type and create example dataset for IHC analysis
library(cowplot)
library(ggplot2)
library(ggpubr)


path<-file.path("A/PATH/FOR/FILES")

#load IHC data :
#format : first column Id , second column Category (PANT, PAT, TNT, TT) and after CD8 et CD20 
#2 rows by unique ID  : 1 for PA , 1 for Tumor 
base_IHC<-read.csv2(file.path(path,"Compilation_IHC_like.csv"),header=TRUE,sep=";")


##############Subpanel A) CD8 
#define pariwise comparisons to put on boxplot
comp.CD8<-list(c("PANT","PAT"),c("TNT","TT"),c("PANT","TNT"),c("PAT","TT"))

#boxplot over categories
p.CD8<-ggplot(data=base_IHC,aes(x=Category,y=CD8,color=Category))+geom_boxplot(show.legend=FALSE)+
  labs(y="CD8",x="Category")+theme_classic()+
  theme(axis.title.x = element_text(size = 13),axis.title.y = element_text(size = 13),axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),axis.line.x = element_line(linewidth = 1),axis.line.y = element_line(linewidth = 1))+
  scale_x_discrete(labels=c("PANT","PAT","TNT","TT"))+
  scale_color_manual(values=c("blue","darkgoldenrod1","darkgray","darkred"))+ 
  stat_compare_means(comparisons=comp.CD8, method="wilcox.test",label="p.signif",bracket.size=1.0)


##############Subpanel A) CD20 
#define pariwise comparisons to put on boxplot
comp.CD20<-list(c("PANT","PAT"),c("TNT","TT"),c("PANT","TNT"),c("PAT","TT"))

#boxplot over categories
p.CD20<-ggplot(data=base_IHC,aes(x=Category,y=CD20,color=Category))+geom_boxplot(show.legend=FALSE)+
  labs(y="CD20",x="Category")+theme_classic()+
  theme(axis.title.x = element_text(size = 13),axis.title.y = element_text(size = 13),axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),axis.line.x = element_line(linewidth = 1),axis.line.y = element_line(linewidth = 1))+
  scale_x_discrete(labels=c("PANT","PAT","TNT","TT"))+
  scale_color_manual(values=c("blue","darkgoldenrod1","darkgray","darkred"))+ 
  stat_compare_means(comparisons=comp.CD20, method="wilcox.test",label="p.signif",bracket.size=1.0)


#combine subpanels over one plot
pFig.S2 <- plot_grid( p.CD8 + theme(legend.position="none"),
                      p.CD20 + theme(legend.position="none"),
                      align = 'vh',
                      labels = c("  (A)  ", "  (B)  "),
                      hjust = -1,
                      nrow = 1
)

#save fig as pdf
ggsave(pFig.S2,file=file.path(path,"FigS2.pdf"),width=16,height=7,units=("in"),dpi=600)
