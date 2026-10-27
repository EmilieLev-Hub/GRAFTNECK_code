#script to obtain results summarized in figure 2 and figure S4 

#wilcoxon + adjustement over pvalues

path<-file.path("A/PATH/FOR/FILES")


#load IHC data :
#format : first column Id , second column Category (PANT, PAT, TNT, TT) and after CD8 et CD20 / PDL1 and CD3
#2 rows by unique ID  : 1 for PA , 1 for Tumor 
base_IHC<-read.csv2(file.path(path,"Compilation_IHC_like.csv"),header=TRUE,sep=";")




#Figure 2 : comparison peritumor area over transplant vs non transp (left column of corresponding table in figure)
# comparison tumor  over transplant vs non transp (right column of corresponding table in figure)
#wilcoxon tests for independant data 

########in peritumor area (left column)
base_IHC_PA<-base_IHC[which(base_IHC$Category%in% c("PANT","PAT")),]

pvalue.PDL1.PA<-signif(wilcox.test(PDL1~Category, data=base_IHC_PA,paired=FALSE)$p.value,2)
pvalue.CD3.PA<-signif(wilcox.test(CD3~Category, data=base_IHC_PA,paired=FALSE)$p.value,2)
pvalue.CD8.PA<-signif(wilcox.test(CD8~Category, data=base_IHC_PA,paired=FALSE)$p.value,2)
pvalue.CD20.PA<-signif(wilcox.test(CD20~Category, data=base_IHC_PA,paired=FALSE)$p.value,2)

vec_pval.PA<-c(pvalue.PDL1.PA,pvalue.CD3.PA,pvalue.CD8.PA,pvalue.CD20.PA)


Tab_Pval_IHC_PA<-as.data.frame(cbind(c("PDL1","CD3","CD8","CD20"),vec_pval.PA))
colnames(Tab_Pval_IHC_PA)<-c(" ","p")


########in tumor (right column)
base_IHC_T<-base_IHC[which(base_IHC$Category%in% c("TNT","TT")),]

pvalue.PDL1.T<-signif(wilcox.test(PDL1~Category, data=base_IHC_T,paired=FALSE)$p.value,2)
pvalue.CD3.T<-signif(wilcox.test(CD3~Category, data=base_IHC_T,paired=FALSE)$p.value,2)
pvalue.CD8.T<-signif(wilcox.test(CD8~Category, data=base_IHC_T,paired=FALSE)$p.value,2)
pvalue.CD20.T<-signif(wilcox.test(CD20~Category, data=base_IHC_T,paired=FALSE)$p.value,2)

vec_pval.T<-c(pvalue.PDL1.T,pvalue.CD3.T,pvalue.CD8.T,pvalue.CD20.T)


Tab_Pval_IHC_T<-as.data.frame(cbind(c("PDL1","CD3","CD8","CD20"),vec_pval.T))
colnames(Tab_Pval_IHC_T)<-c(" ","p")




#Figure S4: comparison peritumor area vs tumor in non transplant patients (left column of corresponding table in figure)
# comparison peritumor area vs tumor in transplant patients (right column of corresponding table in figure)
#wilcoxon tests for paired data 


########in non transplant patients (left column)
base_IHC_NT<-base_IHC[which(base_IHC$Category%in% c("PANT","TNT")),]

pvalue.paired.PDL1.NT<-signif(wilcox.test(PDL1~Category, data=base_IHC_NT,paired=TRUE)$p.value,2)
pvalue.paired.CD3.NT<-signif(wilcox.test(CD3~Category, data=base_IHC_NT,paired=TRUE)$p.value,2)
pvalue.paired.CD8.NT<-signif(wilcox.test(CD8~Category, data=base_IHC_NT,paired=TRUE)$p.value,2)
pvalue.paired.CD20.NT<-signif(wilcox.test(CD20~Category, data=base_IHC_NT,paired=TRUE)$p.value,2)

vec_pval<-c(pvalue.paired.PDL1.NT,pvalue.paired.CD3.NT,pvalue.paired.CD8.NT,pvalue.paired.CD20.NT)
#adjuted pval with bonferroni correction
vec_pval_adj<-vec_pval*4
vec_pval_adj[which(vec_pval_adj>=1)]<-1

Tab_Pval_IHC_NT<-as.data.frame(cbind(c("PDL1","CD3","CD8","CD20"),vec_pval,vec_pval_adj))
colnames(Tab_Pval_IHC_NT)<-c(" ","p","adj.p")


########in transplant patients (right column)
base_IHC_Trans<-base_IHC[which(base_IHC$Category%in% c("PAT","TT")),]

pvalue.paired.PDL1.T<-signif(wilcox.test(PDL1~Category, data=base_IHC_Trans,paired=TRUE)$p.value,2)
pvalue.paired.CD3.T<-signif(wilcox.test(CD3~Category, data=base_IHC_Trans,paired=TRUE)$p.value,2)
pvalue.paired.CD8.T<-signif(wilcox.test(CD8~Category, data=base_IHC_Trans,paired=TRUE)$p.value,2)
pvalue.paired.CD20.T<-signif(wilcox.test(CD20~Category, data=base_IHC_Trans,paired=TRUE)$p.value,2)

vec_pval.trans<-c(pvalue.paired.PDL1.T,pvalue.paired.CD3.T,pvalue.paired.CD8.T,pvalue.paired.CD20.T)
#adjuted pval with bonferroni correction
vec_pval_adj.trans<-vec_pval.trans*4
vec_pval_adj.trans[which(vec_pval_adj.trans>=1)]<-1

Tab_Pval_IHC_Trans<-as.data.frame(cbind(c("PDL1","CD3","CD8","CD20"),vec_pval.trans,vec_pval_adj.trans))
colnames(Tab_Pval_IHC_Trans)<-c(" ","p","adj.p")

