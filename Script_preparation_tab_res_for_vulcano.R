#script for calculcation of fold change and pvalues which are needed to produce vulcano plots in corresponding R scripts

library(rstatix)
library(gtools)
library(tidyverse)
library(purrr)

#path where datasets and source files are
path<-file.path("A/PATH/FOR/FILES")


#Example_data.csv
data_GEP<-read.csv2(file.path(path,"Example_data.csv"),header=TRUE,sep=";",encoding = "latin1",dec=".")
data_GEP$interpretation<-as.factor(data_GEP$interpretation)

#subset peritumor area only and subset on tumor only
data_GEP_Env<-data_GEP[which(data_GEP$interpretation=="Environment"),]
data_GEP_Tum<-data_GEP[which(data_GEP$interpretation=="Tumor"),]


#creation of data.frame where each gene is on row in order to put on columns results of corresponding pvalues and  corresponding fold change 
#get list of all genes (111)
vec_names_genes<-colnames(data_GEP[,c(3,5:115)])[-1]
tab_res<-as.matrix(rep(NA,length(vec_names_genes)),ncol=1,nrow=length(vec_names_genes))
tab_res[,1]<-vec_names_genes
tab_res<-as.data.frame(tab_res)
colnames(tab_res)<-c("genes")


#################################################
#peritumor area and tumor in transplants patients
data_GEP.PA_T.transp<-data_GEP[which(data_GEP$pathology=="Transplanted"),c(3,5:115)]

#calculate pvalues from paired wilcoxon test btw Peritumor Area and Tumor
pvalues.PA_T.transp <- data_GEP.PA_T.transp%>% 
  select(-interpretation) %>% 
  map(~wilcox.test(.x~interpretation, data=data_GEP.PA_T.transp,paired=TRUE)) %>% 
  map_dbl("p.value") #to be at 0 when smaller than 10^-3


pvalues.PA_T.transp <- unname(signif(pvalues.PA_T.transp,2)) 

#to add in data frame "tab_res"
tab_res<-cbind(tab_res,pvalues.PA_T.transp)
colnames(tab_res)[length(tab_res[1,])]<-"raw.pvalues.PA_T.transp"
###################################################


###################################################
#peritumor area and tumor in non-transplants patients
data_GEP.PA_T.cont<-data_GEP[which(data_GEP$pathology=="Controls"),c(3,5:115)]

pvalues.PA_T.cont <- data_GEP.PA_T.cont%>% 
  select(-interpretation) %>% 
  map(~wilcox.test(.x~interpretation, data=data_GEP.PA_T.cont,paired=TRUE)) %>% 
  map_dbl("p.value") 


pvalues.PA_T.cont <- unname(signif(pvalues.PA_T.cont,2)) 

tab_res<-cbind(tab_res,pvalues.PA_T.cont)
colnames(tab_res)[length(tab_res[1,])]<-"raw.pvalues.PA_T.cont"
###################################################


###################################################
#peritumor area in non-transplants patients and transplants patients only on genes
data_GEP_PA<-data_GEP_Env[,c(2,5:115)]

#???pvalues on wilcoxon for independant data (non transp vs transp)
pvalues.PA <- data_GEP_PA%>% 
  select(-pathology) %>% 
  map(~wilcox.test(.x~pathology, data=data_GEP_PA,paired=FALSE)) %>% 
  map_dbl("p.value")



pvalues.PA <- unname(signif(pvalues.PA,2))


tab_res<-cbind(tab_res,pvalues.PA)
colnames(tab_res)[length(tab_res[1,])]<-"raw.pvalues.PA"


# 
tab_res$Signif05_PA<-rep(NA,length(tab_res[,1]))
tab_res$Signif05_PA[which(tab_res$raw.pvalues.PA<=0.05)]<-1
###################################################


###################################################
#tumor in non-transplants patients and transplants patients only on genes
data_GEP_T<-data_GEP_Tum[,c(2,5:115)]

pvalues.T <- data_GEP_T%>% 
  select(-pathology) %>% 
  map(~wilcox.test(.x~pathology, data=data_GEP_T,paired=FALSE)) %>% 
  map_dbl("p.value")

pvalues.T <- unname(signif(pvalues.T,2) )

tab_res<-cbind(tab_res,pvalues.T)
colnames(tab_res)[length(tab_res[1,])]<-"raw.pvalues.T"

tab_res$Signif05_T<-rep(NA,length(tab_res[,1]))
tab_res$Signif05_T[which(tab_res$raw.pvalues.T<=0.05)]<-1
########################################################




########################################################################################
#calculate med fold change
########################################################################################



########################################################################################
#on peritumor area in transp vs non-transp (controls) patients
foldchang_PA<-as.matrix(tab_res[,1:2],nrow=length(vec_names_genes)+1,ncol=2)
foldchang_PA[,1]<-colnames(data_GEP_Env)[5:115]

for(j in 5:115){
  foldchang_PA[j-4,2]<-round(median(data_GEP_Env[which(data_GEP_Env$pathology=="Transplanted"),j],na.rm=TRUE)/median(data_GEP_Env[which(data_GEP_Env$pathology=="Controls"),j],na.rm=TRUE),4)
                       #calculate median over transp and non transp as reference
  }
colnames(foldchang_PA)<-c("genes","fold change med")

foldchang_PA[which(foldchang_PA=="NaN")]<-NA
foldchang_PA[which(foldchang_PA=="Inf")]<-NA
foldchang_PA[,2]<-as.numeric(as.character(foldchang_PA[,2]))
#########################################################################################


#########################################################################################
#on tumor in transp vs non-transp (controls) patients
foldchang_T<-as.matrix(tab_res[,1:2],nrow=length(vec_names_genes)+1,ncol=2)
foldchang_T[,1]<-colnames(data_GEP_Tum)[5:115]

for(j in 5:115){
  foldchang_T[j-4,2]<-round(median(data_GEP_Tum[which(data_GEP_Tum$pathology=="Transplanted"),j],na.rm=TRUE)/median(data_GEP_Tum[which(data_GEP_Tum$pathology=="Controls"),j],na.rm=TRUE),4)
}

colnames(foldchang_T)<-c("genes","fold change med")

foldchang_T[which(foldchang_T=="NaN")]<-NA
foldchang_T[which(foldchang_T=="Inf")]<-NA

foldchang_T[,2]<-as.numeric(as.character(foldchang_T[,2]))
########################################################################################


#add to data.frame tab_res which contains pvalues of all comparisons
tab_res<-cbind(tab_res,foldchang_T[,2],foldchang_PA[,2])
colnames(tab_res)[c(8,9)]<-c("fold_chang_T","fold_chang_PA")

tab_res$fold_chang_T<-as.numeric(as.character(tab_res$fold_chang_T))
tab_res$fold_chang_PA<-as.numeric(as.character(tab_res$fold_chang_PA))


#add log2 for vulcano
tab_res$fold_chang_PA.log2<-log2(tab_res$fold_chang_PA)
tab_res$fold_chang_T.log2<-log2(tab_res$fold_chang_T)




##########
#comparison for PA and T in non transplant patients 
#ratio btw gene value PA over gene value Tumor for each pair and after calculate median for each gene 

data_GEP_NT<-data_GEP[which(data_GEP$pathology=="Controls"),]
foldchang_PA_T_NT<-as.matrix(data_GEP_NT[1:length(unique(data_GEP_NT$Id)),c(1,5:115)],nrow=length(unique(data_GEP_NT$Identif_paire)),ncol=length(vec_names_genes)+1)
foldchang_PA_T_NT[,1]<-unique(data_GEP_NT$Id)

for(i in 1:length(unique(data_GEP_NT$Id))){
  
  ind<-which(data_GEP_NT$Id==unique(data_GEP_NT$Id)[i])
  if(length(ind)>0){
    for(j in 5:115){
      foldchang_PA_T_NT[i,j-3]<-abs(data_GEP_NT[ind[1],j]/data_GEP_NT[ind[2],j])  #1=PA and 2: T
    }
  }
}

foldchang_PA_T_NT[which(foldchang_PA_T_NT=="NaN")]<-NA
foldchang_PA_T_NT[which(foldchang_PA_T_NT=="Inf")]<-NA


##########
#comparison for PA and T in transplant patients 
data_GEP_Trans<-data_GEP[which(data_GEP$pathology=="Transplanted"),]
foldchang_PA_T_T<-as.matrix(data_GEP_Trans[1:length(unique(data_GEP_Trans$Id)),c(1,5:115)],nrow=length(unique(data_GEP_Trans$Id)),ncol=length(vec_names_genes)+1)
foldchang_PA_T_T[,1]<-unique(data_GEP_Trans$Id)

for(i in 1:length(unique(data_GEP_Trans$Id))){
  
  ind<-which(data_GEP_Trans$Id==unique(data_GEP_Trans$Id)[i])
  if(length(ind)>0){
    for(j in 5:115){
      foldchang_PA_T_T[i,j-3]<-abs(data_GEP_Trans[ind[1],j]/data_GEP_Trans[ind[2],j])
    }
  }
}

foldchang_PA_T_T[which(foldchang_PA_T_T=="NaN")]<-NA
foldchang_PA_T_T[which(foldchang_PA_T_T=="Inf")]<-NA


#calculation median fold change for each gene
#in transplants patients
med_foldchange_PA_T_T<-rep(NA,length(vec_names_genes))
for(j in 2:(length(vec_names_genes)+1)){
  med_foldchange_PA_T_T[j-1]<- round(median(foldchang_PA_T_T[,j],na.rm=TRUE),4)
}

#in non-transplant patients
med_foldchange_PA_T_NT<-rep(NA,length(vec_names_genes))
for(j in 2:(length(vec_names_genes)+1)){
  med_foldchange_PA_T_NT[j-1]<- round(median(foldchang_PA_T_NT[,j],na.rm=TRUE),4)
}

#on log2 transformation for vulcano plot
med_fold_T_log2<-log(med_foldchange_PA_T_T,2)
med_fold_NT_log2<-log(med_foldchange_PA_T_NT,2)

med_fold_T_log2[which(med_fold_T_log2=="-Inf")]<-NA
med_fold_NT_log2[which(med_fold_NT_log2=="-Inf")]<-NA


#add to data.frame tab_res which contains pvalues of all comparisons
tab_res<-cbind(tab_res, med_foldchange_PA_T_T, med_foldchange_PA_T_NT)
colnames(tab_res)[c(12,13)]<-c("fold_chang_PA_T_Trans","fold_chang_PA_T_NT")

#add log2 for vulcano
tab_res$fold_chang_Trans.log2<-med_fold_T_log2
tab_res$fold_chang_NT.log2<-med_fold_NT_log2

#save csv file with all res pvalues and med fold change for vulcano plots
write.csv2(tab_res,file=file.path(path,"tab_res_all_comp.csv"))
