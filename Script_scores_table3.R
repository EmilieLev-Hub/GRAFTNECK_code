
library(rstatix)
library(gtools)


path<-file.path("A/PATH/FOR/FILES")
#path where your database was saved


#the main file where gene date is available for each patient, each environement and each tumor component
#there were 2 lines for each patient : 1 for its environment component and 1 for its tumor component
data_GEP<-read.csv2(file.path(path,"Example_data.csv"),header=TRUE,sep=";",encoding = "latin1",dec=".")


#status variable : Environment / Tumor
data_GEP$interpretation<-as.factor(data_GEP$interpretation)



#hypoxia score : 
#[VEGF-1, FLT1, KDR, VHL, HIF1A, CA9, NDUFA4L2]
data_GEP$SCORE_Hypoxia<-rep(NA,length(data_GEP$Id))
for(i in 1 :length(data_GEP$Id)){
  
  data_GEP$SCORE_Hypoxia[i]<-0
  if(data_GEP$VEGF.A[i]>= round(unname(quantile(data_GEP$VEGF.A[which(data_GEP$VEGF.A!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){ data_GEP$SCORE_Hypoxia[i]<-data_GEP$SCORE_Hypoxia[i]+1}
  if(data_GEP$FLT1[i]>= round(unname(quantile(data_GEP$FLT1[which(data_GEP$FLT1!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Hypoxia[i]<-data_GEP$SCORE_Hypoxia[i]+1}
  if(data_GEP$KDR[i]>= round(unname(quantile(data_GEP$KDR[which(data_GEP$KDR!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Hypoxia[i]<-data_GEP$SCORE_Hypoxia[i]+1}
  if(data_GEP$VHL[i]>= round(unname(quantile(data_GEP$VHL[which(data_GEP$VHL!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Hypoxia[i]<-data_GEP$SCORE_Hypoxia[i]+1}
  if(data_GEP$HIF1A[i]>= round(unname(quantile(data_GEP$HIF1A[which(data_GEP$HIF1A!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Hypoxia[i]<-data_GEP$SCORE_Hypoxia[i]+1}
  if(data_GEP$CA9[i]>=round(unname(quantile(data_GEP$CA9[which(data_GEP$CA9!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Hypoxia[i]<-data_GEP$SCORE_Hypoxia[i]+1}
  if(data_GEP$NDUFA4L2[i]>=round(unname(quantile(data_GEP$NDUFA4L2[which(data_GEP$NDUFA4L2!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Hypoxia[i]<-data_GEP$SCORE_Hypoxia[i]+1}
  
}


#apoptotic score :
#Bcl2, DIABLO, XIAP, BIRC7, BIRC5, BIRC6 were considered as apoptotic inhibitors :-1
#CASP3, CASP7, FasL, FADD, TRADD, BAX as apoptotic inducers:+1
data_GEP$SCORE_Apoptotic<-rep(NA,length(data_GEP$Id))
for(i in 1 :length(data_GEP$Id)){
  
  data_GEP$SCORE_Apoptotic[i]<-0
  if(data_GEP$CASP3[i]>=round(unname(quantile(data_GEP$CASP3[which(data_GEP$CASP3!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Apoptotic[i]<-data_GEP$SCORE_Apoptotic[i]+1}
  if(data_GEP$CASP7[i]>=round(unname(quantile(data_GEP$CASP7[which(data_GEP$CASP7!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Apoptotic[i]<-data_GEP$SCORE_Apoptotic[i]+1}
  if(data_GEP$FasL[i]>=round(unname(quantile(data_GEP$FasL[which(data_GEP$FasL!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Apoptotic[i]<-data_GEP$SCORE_Apoptotic[i]+1}
  if(data_GEP$FADD[i]>=round(unname(quantile(data_GEP$FADD[which(data_GEP$FADD!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Apoptotic[i]<-data_GEP$SCORE_Apoptotic[i]+1}
  if(data_GEP$TRADD[i]>=round(unname(quantile(data_GEP$TRADD[which(data_GEP$TRADD!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Apoptotic[i]<-data_GEP$SCORE_Apoptotic[i]+1}
  if(data_GEP$BAX[i]>=round(unname(quantile(data_GEP$BAX[which(data_GEP$BAX!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Apoptotic[i]<-data_GEP$SCORE_Apoptotic[i]+1}
  
  if(data_GEP$Bcl2[i]>=round(unname(quantile(data_GEP$Bcl2[which(data_GEP$Bcl2!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){ data_GEP$SCORE_Apoptotic[i]<-data_GEP$SCORE_Apoptotic[i]-1}
  if(data_GEP$DIABLO[i]>=round(unname(quantile(data_GEP$DIABLO[which(data_GEP$DIABLO!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Apoptotic[i]<-data_GEP$SCORE_Apoptotic[i]-1}
  if(data_GEP$XIAP[i]>=round(unname(quantile(data_GEP$XIAP[which(data_GEP$XIAP!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Apoptotic[i]<-data_GEP$SCORE_Apoptotic[i]-1}
  if(data_GEP$BIRC5[i]>=round(unname(quantile(data_GEP$BIRC7[which(data_GEP$BIRC7!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Apoptotic[i]<-data_GEP$SCORE_Apoptotic[i]-1}
  if(data_GEP$BIRC6[i]>=round(unname(quantile(data_GEP$BIRC5[which(data_GEP$BIRC5!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Apoptotic[i]<-data_GEP$SCORE_Apoptotic[i]-1}
  if(data_GEP$BIRC7[i]>=round(unname(quantile(data_GEP$BIRC6[which(data_GEP$BIRC6!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Apoptotic[i]<-data_GEP$SCORE_Apoptotic[i]-1}
  
}


#metabolic score :
#[GLUT1, GLUT3, GLUT4, GLUT5, GLUT6, SGLT1, LDHA, LDHB]
data_GEP$SCORE_Metabolic<-rep(NA,length(data_GEP$Id))
for(i in 1 :length(data_GEP$Id)){
  
  data_GEP$SCORE_Metabolic[i]<-0
  if(data_GEP$GLUT1[i]>=round(unname(quantile(data_GEP$GLUT1[which(data_GEP$GLUT1!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){ data_GEP$SCORE_Metabolic[i]<-data_GEP$SCORE_Metabolic[i]+1}
  if(data_GEP$GLUT3[i]>=round(unname(quantile(data_GEP$GLUT3[which(data_GEP$GLUT3!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Metabolic[i]<-data_GEP$SCORE_Metabolic[i]+1}
  if(data_GEP$GLUT4[i]>=round(unname(quantile(data_GEP$GLUT4[which(data_GEP$GLUT4!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Metabolic[i]<-data_GEP$SCORE_Metabolic[i]+1}
  if(data_GEP$GLUT5[i]>=round(unname(quantile(data_GEP$GLUT5[which(data_GEP$GLUT5!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Metabolic[i]<-data_GEP$SCORE_Metabolic[i]+1}
  if(data_GEP$GLUT6[i]>=round(unname(quantile(data_GEP$GLUT6[which(data_GEP$GLUT6!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Metabolic[i]<-data_GEP$SCORE_Metabolic[i]+1}
  if(data_GEP$SGLT1[i]>=round(unname(quantile(data_GEP$SGLT1[which(data_GEP$SGLT1!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Metabolic[i]<-data_GEP$SCORE_Metabolic[i]+1}
  if(data_GEP$LDHA[i]>=round(unname(quantile(data_GEP$LDHA[which(data_GEP$LDHA!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Metabolic[i]<-data_GEP$SCORE_Metabolic[i]+1}
  if(data_GEP$LDHB[i]>=round(unname(quantile(data_GEP$LDHB[which(data_GEP$LDHB!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Metabolic[i]<-data_GEP$SCORE_Metabolic[i]+1}
  
}


#proliferation score 
#[AURKA, Rb, EGFR, GARP, HGF, CCND1, CDK4, CDK6, CDKN2A, WEE1, MAPK, PI3K, FAT1, PLK1]
data_GEP$SCORE_proliferation<-rep(NA,length(data_GEP$Id))
for(i in 1 :length(data_GEP$Id)){
  
  data_GEP$SCORE_proliferation[i]<-0
  if(data_GEP$AURKA[i]>=round(unname(quantile(data_GEP$AURKA[which(data_GEP$AURKA!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){ data_GEP$SCORE_proliferation[i]<-data_GEP$SCORE_proliferation[i]+1}
  if(data_GEP$Rb[i]>=round(unname(quantile(data_GEP$Rb[which(data_GEP$Rb!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_proliferation[i]<-data_GEP$SCORE_proliferation[i]+1}
  if(data_GEP$EGFR[i]>=round(unname(quantile(data_GEP$EGFR[which(data_GEP$EGFR!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_proliferation[i]<-data_GEP$SCORE_proliferation[i]+1}
  if(data_GEP$GARP[i]>=round(unname(quantile(data_GEP$GARP[which(data_GEP$GARP!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_proliferation[i]<-data_GEP$SCORE_proliferation[i]+1}
  if(data_GEP$HGF[i]>=round(unname(quantile(data_GEP$HGF[which(data_GEP$HGF!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_proliferation[i]<-data_GEP$SCORE_proliferation[i]+1}
  if(data_GEP$CCND1[i]>=round(unname(quantile(data_GEP$CCND1[which(data_GEP$CCND1!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_proliferation[i]<-data_GEP$SCORE_proliferation[i]+1}
  if(data_GEP$CDK4[i]>=round(unname(quantile(data_GEP$CDK4[which(data_GEP$CDK4!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_proliferation[i]<-data_GEP$SCORE_proliferation[i]+1}
  if(data_GEP$CDK6[i]>=round(unname(quantile(data_GEP$CDK6[which(data_GEP$CDK6!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_proliferation[i]<-data_GEP$SCORE_proliferation[i]+1}
  if(data_GEP$CDKN2A[i]>=round(unname(quantile(data_GEP$CDKN2A[which(data_GEP$CDKN2A!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_proliferation[i]<-data_GEP$SCORE_proliferation[i]+1}
  if(data_GEP$WEE1[i]>=round(unname(quantile(data_GEP$WEE1[which(data_GEP$WEE1!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_proliferation[i]<-data_GEP$SCORE_proliferation[i]+1}
  if(data_GEP$MAPK[i]>=round(unname(quantile(data_GEP$MAPK[which(data_GEP$MAPK!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_proliferation[i]<-data_GEP$SCORE_proliferation[i]+1}
  if(data_GEP$PI3K[i]>=round(unname(quantile(data_GEP$PI3K[which(data_GEP$PI3K!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_proliferation[i]<-data_GEP$SCORE_proliferation[i]+1}
  if(data_GEP$FAT1[i]>=round(unname(quantile(data_GEP$FAT1[which(data_GEP$FAT1!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_proliferation[i]<-data_GEP$SCORE_proliferation[i]+1}
  if(data_GEP$PLK1[i]>=round(unname(quantile(data_GEP$PLK1[which(data_GEP$PLK1!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))),2)){data_GEP$SCORE_proliferation[i]<-data_GEP$SCORE_proliferation[i]+1}
  
}

#inflammatory score
#[CD3, CD4, TIA1, PRF, CD5, GRB, CD8, MS4A1, CR2, SDC1, CD68, Lysosyme, CD56, DNTT, CD25, FOXP3, CD34, CD31, ITGAE, ENTPD1, KIR3DL2]
data_GEP$SCORE_Infl<-rep(NA,length(data_GEP$Id))
for(i in 1 :length(data_GEP$Id)){
  
  data_GEP$SCORE_Infl[i]<-0
  if(data_GEP$CD3[i]>=round(unname(quantile(data_GEP$CD3[which(data_GEP$CD3!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){ data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$CD4[i]>=round(unname(quantile(data_GEP$CD4[which(data_GEP$CD4!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$TIA1[i]>=round(unname(quantile(data_GEP$TIA1[which(data_GEP$TIA1!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$PRF[i]>=round(unname(quantile(data_GEP$PRF[which(data_GEP$PRF!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$GRB[i]>=round(unname(quantile(data_GEP$GRB[which(data_GEP$GRB!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$CD5[i]>=round(unname(quantile(data_GEP$CD5[which(data_GEP$CD5!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$CD8[i]>=round(unname(quantile(data_GEP$CD8[which(data_GEP$CD8!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$MS4A1[i]>=round(unname(quantile(data_GEP$MS4A1[which(data_GEP$MS4A1!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$CR2[i]>=round(unname(quantile(data_GEP$CR2[which(data_GEP$CR2!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$SDC1[i]>=round(unname(quantile(data_GEP$SDC1[which(data_GEP$SDC1!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$CD68[i]>=round(unname(quantile(data_GEP$CD68[which(data_GEP$CD68!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$Lysosyme[i]>=round(unname(quantile(data_GEP$Lysosyme[which(data_GEP$Lysosyme!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$CD56[i]>=round(unname(quantile(data_GEP$CD56[which(data_GEP$CD56!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$DNTT[i]>=round(unname(quantile(data_GEP$DNTT[which(data_GEP$DNTT!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$CD25[i]>=round(unname(quantile(data_GEP$CD25[which(data_GEP$CD25!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$FOXP3[i]>=round(unname(quantile(data_GEP$FOXP3[which(data_GEP$FOXP3!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$CD34[i]>=round(unname(quantile(data_GEP$CD34[which(data_GEP$CD34!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$CD31[i]>=round(unname(quantile(data_GEP$CD31[which(data_GEP$CD31!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$ENTPD1[i]>=round(unname(quantile(data_GEP$ENTPD1[which(data_GEP$ENTPD1!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$ITGAE[i]>=round(unname(quantile(data_GEP$ITGAE[which(data_GEP$ITGAE!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
  if(data_GEP$KIR3DL2[i]>=round(unname(quantile(data_GEP$KIR3DL2[which(data_GEP$KIR3DL2!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Infl[i]<-data_GEP$SCORE_Infl[i]+1}
}




#immune synapse regulatory score
#CTLA4, PDCD1, PDCD1LG1, TIM3, CD86, TIGIT, SIRPA  as inhibitors :-1
#TNFRSF4, CD40, CD40L, ILDR2, CCR4, LAG3, ICOS, CD80, CD28, GITR as activators : +1
data_GEP$SCORE_Synapse<-rep(NA,length(data_GEP$Id))

for(i in 1 :length(data_GEP$Id)){
  #-1 for inhibitor 
  data_GEP$SCORE_Synapse[i]<-0
  if(data_GEP$CTLA4[i]>=round(unname(quantile(data_GEP$CTLA4[which(data_GEP$CTLA4!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){ data_GEP$SCORE_Synapse[i]<-data_GEP$SCORE_Synapse[i]-1}
  if(data_GEP$PD1[i]>=round(unname(quantile(data_GEP$PD1[which(data_GEP$PD1!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Synapse[i]<-data_GEP$SCORE_Synapse[i]-1}
  if(data_GEP$PDL1[i]>=round(unname(quantile(data_GEP$PDL1[which(data_GEP$PDL1!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Synapse[i]<-data_GEP$SCORE_Synapse[i]-1}
  if(data_GEP$TIM3[i]>=round(unname(quantile(data_GEP$TIM3[which(data_GEP$TIM3!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Synapse[i]<-data_GEP$SCORE_Synapse[i]-1}
  if(data_GEP$CD86[i]>=round(unname(quantile(data_GEP$CD86[which(data_GEP$CD86!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Synapse[i]<-data_GEP$SCORE_Synapse[i]-1}
  if(data_GEP$TIGIT[i]>=round(unname(quantile(data_GEP$TIGIT[which(data_GEP$TIGIT!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Synapse[i]<-data_GEP$SCORE_Synapse[i]-1}
  if(data_GEP$SIRPA[i]>=round(unname(quantile(data_GEP$SIRPA[which(data_GEP$SIRPA!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Synapse[i]<-data_GEP$SCORE_Synapse[i]-1}
  
  #+1 for activator
  if(data_GEP$TNFRSF4[i]>=round(unname(quantile(data_GEP$TNFRSF4[which(data_GEP$TNFRSF4!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){ data_GEP$SCORE_Synapse[i]<-data_GEP$SCORE_Synapse[i]+1}
  if(data_GEP$CD40[i]>=round(unname(quantile(data_GEP$CD40[which(data_GEP$CD40!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Synapse[i]<-data_GEP$SCORE_Synapse[i]+1}
  if(data_GEP$CD40LE2_3[i]>=round(unname(quantile(data_GEP$CD40LE2_3[which(data_GEP$CD40LE2_3!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Synapse[i]<-data_GEP$SCORE_Synapse[i]+1}
  if(data_GEP$LAG3[i]>=round(unname(quantile(data_GEP$LAG3[which(data_GEP$LAG3!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Synapse[i]<-data_GEP$SCORE_Synapse[i]+1}
  if(data_GEP$ICOS[i]>=round(unname(quantile(data_GEP$ICOS[which(data_GEP$ICOS!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Synapse[i]<-data_GEP$SCORE_Synapse[i]+1}
  if(data_GEP$CD80[i]>=round(unname(quantile(data_GEP$CD80[which(data_GEP$CD80!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Synapse[i]<-data_GEP$SCORE_Synapse[i]+1}
  if(data_GEP$CD28[i]>=round(unname(quantile(data_GEP$CD28[which(data_GEP$CD28!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Synapse[i]<-data_GEP$SCORE_Synapse[i]+1}
  if(data_GEP$TNFRSF18[i]>=round(unname(quantile(data_GEP$TNFRSF18[which(data_GEP$TNFRSF18!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Synapse[i]<-data_GEP$SCORE_Synapse[i]+1}
  if(data_GEP$ILDR2[i]>=round(unname(quantile(data_GEP$ILDR2[which(data_GEP$ILDR2!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Synapse[i]<-data_GEP$SCORE_Synapse[i]+1}
  if(data_GEP$CCR4[i]>=round(unname(quantile(data_GEP$CCR4[which(data_GEP$CCR4!=0)],probs=c(0.25,0.5,0.75,0.9,0.95))[3]),2)){data_GEP$SCORE_Synapse[i]<-data_GEP$SCORE_Synapse[i]+1}
  
  
}



#2 subsets : 1 for Controls only and 1 for Transplant only
data_GEP_Controls<-data_GEP[which(data_GEP$pathology=="Controls"),]
data_GEP_Trans<-data_GEP[which(data_GEP$pathology=="Transplanted"),]

#comparison tumour of non-transplant (TNT) vs tumor transplant (TT)
data_Tum<-rbind(data_GEP_Controls[which(data_GEP_Controls$interpretation=="Tumor"),],data_GEP_Trans[which(data_GEP_Trans$interpretation=="Tumor"),] ) 
data_Tum$pathology<-as.factor(data_Tum$pathology)

source(file.path(path,"report_functions_20210324.R"),encoding="latin1")

stat_base.score.hypoxia   <- StatDesc( data_Tum$SCORE_Hypoxia,"Hypoxia score", data_Tum$pathology, 
                                       stats="meansd med min max", digits = 2,test = "Wilcoxon")


stat_base.score.apopt  <- StatDesc( data_Tum$SCORE_Apoptotic,"Apoptotic score", data_Tum$pathology, 
                                         stats="meansd med min max", digits = 2,test = "Wilcoxon")

stat_base.score.metabo  <- StatDesc( data_Tum$SCORE_Metabolic,"Metabolic score", data_Tum$pathology, 
                                     stats="meansd med min max", digits = 2,test = "Wilcoxon")

stat_base.score.prolif  <- StatDesc( data_Tum$SCORE_proliferation,"Proliferation score", data_Tum$pathology, 
                                       stats="meansd med min max", digits = 2,test = "Wilcoxon")


stat_base.score.infl <- StatDesc( data_Tum$SCORE_Infl,"Inflammatory score", data_Tum$pathology, 
                                            stats="meansd med min max", digits = 2,test = "Wilcoxon")


stat_base.score.synapse  <- StatDesc( data_Tum$SCORE_Synapse,"Immune synapse regulatory score", data_Tum$pathology, 
                                        stats="meansd med min max", digits = 2,test = "Wilcoxon")


TableScore.tum<-rbind(stat_base.score.hypoxia,  stat_base.score.apopt,
                      stat_base.score.metabo, stat_base.score.prolif, 
                      stat_base.score.infl, stat_base.score.synapse)
TableScore.tum <- TableScore.tum[,1:4]
names(TableScore.tum)<-c(" ",levels(data_Tum$pathology),"pvalue")




#comparison peritumor area of non-transplant (PANT) vs tumor transplant (PAT)
data_Env<-rbind(data_GEP_Controls[which(data_GEP_Controls$interpretation=="Environment"),],data_GEP_Trans[which(data_GEP_Trans$interpretation=="Environment"),] ) 

data_Env$pathology<-as.factor(data_Env$pathology)
stat_base.score.hypoxia   <- StatDesc( data_Env$SCORE_Hypoxia,"Hypoxia score ", data_Env$pathology, 
                                       stats="meansd med min max", digits = 2,test = "Wilcoxon")

stat_base.score.apopt <- StatDesc( data_Env$SCORE_Apoptotic,"Apoptotic score", data_Env$pathology, 
                                         stats="meansd med min max", digits = 2,test = "Wilcoxon")

stat_base.score.metabo  <- StatDesc( data_Env$SCORE_Metabolic,"Metabolic score", data_Env$pathology, 
                                     stats="meansd med min max", digits = 2,test = "Wilcoxon")

stat_base.score.prolif  <- StatDesc( data_Env$SCORE_proliferation,"Proliferation score", data_Env$pathology, 
                                       stats="meansd med min max", digits = 2,test = "Wilcoxon")

stat_base.score.infl  <- StatDesc( data_Env$SCORE_Infl,"Inflammatory score", data_Env$pathology, 
                                            stats="meansd med min max", digits = 2,test = "Wilcoxon")

stat_base.score.synapse  <- StatDesc( data_Env$SCORE_Synapse,"Immune synapse regulatory score", data_Env$pathology, 
                                        stats="meansd med min max", digits = 2,test = "Wilcoxon")


TableScore.env<-rbind(stat_base.score.hypoxia, stat_base.score.apopt,
                      stat_base.score.metabo, stat_base.score.prolif, 
                      stat_base.score.infl,  stat_base.score.synapse)
TableScore.env <- TableScore.env[,1:4]
names(TableScore.env)<-c(" ",levels(data_Env$pathology),"pvalue")


#add adj.p : 1) transform to numeric data 
TableScore.env[,4]<-as.numeric(as.character(substr(TableScore.env[,4],1,5)))
TableScore.tum[,4]<-as.numeric(as.character(substr(TableScore.tum[,4],1,5)))

#2): calculation pvalue: *6 since there are 6 scores
pvalues_ajust_Bonf_env<-as.numeric(as.character(TableScore.env[,4]))*6
pvalues_ajust_Bonf_env[which(pvalues_ajust_Bonf_env>=1)]<-1
pvalues_ajust_Bonf_env[which(is.na(pvalues_ajust_Bonf_env)==TRUE)]<-" "

pvalues_ajust_Bonf_tum<-as.numeric(as.character(TableScore.tum[,4]))*6
pvalues_ajust_Bonf_tum[which(pvalues_ajust_Bonf_tum>=1)]<-1
pvalues_ajust_Bonf_tum[which(is.na(pvalues_ajust_Bonf_tum)==TRUE)]<-" "

#3) merge all results for comparisons
TableScore.comp.tum.env<-cbind(TableScore.env,pvalues_ajust_Bonf_env,TableScore.tum[,2:4],pvalues_ajust_Bonf_tum)
colnames(TableScore.comp.tum.env)<-c(" ", "PANT","PAT","p" ,"Adj.p",
                                     "TNT", "TT","p","Adj.p")


TableScore.comp.tum.env[which(is.na(TableScore.comp.tum.env[,4])==TRUE),4]<-" "
TableScore.comp.tum.env[which(is.na(TableScore.comp.tum.env[,8])==TRUE),8]<-" "

#4) write csv for Table3 
write.csv2(TableScore.comp.tum.env,file=file.path(path,"Results_Table3.csv"))











