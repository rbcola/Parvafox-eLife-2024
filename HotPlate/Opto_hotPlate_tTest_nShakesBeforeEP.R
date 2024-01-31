################################################
##### T-Test analysis of Opto hot plate shakes##
################################################

rm(list=ls())

library(reshape2)
library(openxlsx)
library(plyr)
library(pastecs)
library(car)
library(ggplot2)
library(dplyr)

#create outputdir
output_dir=paste(getwd(),"/nShakesBeforeEP/", sep="")
ifelse(!dir.exists(output_dir), dir.create(output_dir),FALSE)

#create empty DF to store results in
Results=data.frame(matrix(nrow=0, ncol=8))
colnames(Results)=c("Opsin","Test","t-value","df","p-value","Hedges-g","Genotype","Hedges-Magnitude")

#create functions to calculate Effect size
rFromTTest=function(myT, mydf){
  myR0= sqrt(myT^2/(myT^2+mydf))
  myR = round(myR0, 3)
}

rFromWilcox= function(wilcoxModel, N){
  z= qnorm(wilcoxModel$p.value/2)
  r= z/sqrt(N)
}

#read data
myData=read.xlsx("D:/Reto (old harddrive)/UNIFR/MSc.EBR/Masterarbeit/HotPlate/22-12-21_Opto_HotPlate_allAnimals/22-09-07_Opto_HotPlate_nShakesBeforeEP.xlsx")
#only keep avg columns
#myData <- myData[,c(1,9,10)]

#convert to long format
longData=melt(myData, id.vars=("Animal"), variable.name = "Condition", value.name = "Value")

#convert all columns to factor except the dependent Variable
for (ii in 1:ncol(longData)) {
  if (colnames(longData[ii]) != "Value") {
    longData[, ii] = as.factor(longData [, ii])
  }
  else{
    longData[, ii] = longData[, ii]
  }
}

#change colnames
colnames(myData)[c(2,3)] <- c("LEDoff","LEDon")
wideData <- myData

#define group allocation per animal
ChR2_list=c("61","75","76","77","106-21/A","106-21/B", 
            "106-21/C","106-21/D","106-21/10","34-21/7",
            "34-21/11","34-21/10","39-21/5","39-21/3",
            "39-21/2","39-21/6","39-21/4")
eArch3.0_list=c("4","10","11","15","16")
ConFonChR2=c("71","93")
eArch3.0_Control=c("52")
PVCre_list <- c("39-21/5","39-21/3",
                "39-21/2","39-21/4","39-21/6")
wideData$Opsin <- NA
wideData$Genotype <- NA

for(kk in 1:nrow(wideData)){
  if (is.element(wideData[kk,1],ChR2_list)){
    wideData[kk,4]="ChR2"
  }else if(is.element(wideData[kk,1],eArch3.0_list)){
    wideData[kk,4]="ArchT3.0"
  }else if(is.element(wideData[kk,1],eArch3.0_Control)){
    wideData[kk,4]="ArchT3.0_Control"
  }else if(is.element(wideData[kk,1],ConFonChR2)){
    wideData[kk,4]="Con_Fon_ChR2"
  }
}

for(oo in 1:nrow(wideData)){
  if (is.element(wideData[oo,1],PVCre_list)){
    wideData[oo,5]="PV-Cre"
  }else if(wideData[oo,4] == "Con_Fon_ChR2"){
    wideData[oo,5]="Foxb1-Cre/PV-Flp"
  }else{
    wideData[oo,5]="Foxb1-Cre"
  }
}

#creat DF for later plotting
plotData=melt(wideData, id.var=c("Animal","Opsin","Genotype"), variable.name = "Condition", value.name = "Value")
#convert all columns to factor except the dependent Variable
for (ll in 1:ncol(plotData)) {
  if (colnames(plotData[ll]) != "Value") {
    plotData[, ll] = as.factor(plotData [, ll])
  }
  else{
    plotData[, ll] = plotData[, ll]
  }
}

#filter out Opsin-Conditions with less than 3 animals per group
opsins <- unique(wideData$Opsin)
for (ee in opsins){
  if(sum(wideData$Opsin == ee) < 3){
    wideData <- wideData %>% dplyr::filter(Opsin != ee)
  }
}

#loop through the Genotype and Opsin groups to run the analysis
for (qq in unique(wideData$Genotype)){
  
  qqData <- wideData %>% dplyr::filter(Genotype == qq)
  
  for (jj in unique(qqData$Opsin)){
    print(jj)
    #filter qqData by DREADD_receptor
    groupDF=qqData[qqData$Opsin==jj,]
    #create difference vector
    diffVec=groupDF$LEDoff-groupDF$LEDon
    
    #test for normality
    normality=stat.desc(diffVec,basic=FALSE , norm= TRUE)
    #print(qqPlot(diffVec))
    
    #run the appropriate test
    if( normality["normtest.p"]> 0.05 & normality["skew.2SE"]<1 & normality["kurt.2SE"]<1){
      hypTest=t.test(groupDF$LEDoff, groupDF$LEDon, paired= TRUE)
      Test=paste(hypTest$method,hypTest$alternative,sep=" ")
      EffSize=cohen.d(groupDF$LEDoff, groupDF$LEDon, hedges.correction = TRUE, paired = TRUE)
      #Store the results from the hypothesis test in Results
      Results[nrow(Results)+1,1]=jj
      Results[nrow(Results),2]=Test
      Results[nrow(Results),3]=hypTest[[1]]
      Results[nrow(Results),4]=hypTest[[2]]
      Results[nrow(Results),5]=hypTest[[3]]
      Results[nrow(Results),6]=EffSize$estimate
      Results[nrow(Results),7]=qq
      Results[nrow(Results),8]=levels(EffSize$magnitude)[EffSize$magnitude]
    }else{
      hypTest=wilcox.test(groupDF$LEDoff, groupDF$LEDon, paired= TRUE, exact=TRUE)
      Test=paste(hypTest$method,hypTest$alternative,sep=" ")
      EffSize=cohen.d(groupDF$LEDoff, groupDF$LEDon, hedges.correction = TRUE, paired = TRUE)
      #Store the results from the hypothesis test in Results
      Results[nrow(Results)+1,1]=jj
      Results[nrow(Results),2]=Test
      Results[nrow(Results),3]=NA
      Results[nrow(Results),4]=NA
      Results[nrow(Results),5]=hypTest$p.value
      Results[nrow(Results),6]=EffSize$estimate
      Results[nrow(Results),7]=qq
      Results[nrow(Results),8]=levels(EffSize$magnitude)[EffSize$magnitude]
    }
    
    #save inividual hypothesis test output
    nam=paste(output_dir,"Opto_nShakesBeforeEP_",qq,"_",jj,"_paired_tTest-Wilcoxon",".doc",sep="")
    capture.output(print(hypTest),file=nam)
  }  
  
  #plot data
  long_qqData <- melt(qqData, id.var=c("Animal","Opsin","Genotype"), variable.name = "Condition", value.name = "Value")
  gg <-  ggplot(long_qqData, aes(x=`Condition`, y=`Value`, fill=Condition)) +
    #geom_boxplot(stat = "summary", fun = "mean", position = "dodge")+
    stat_boxplot(geom="errorbar", position = position_dodge(0.8), size = 1, width = 0.5)+
    geom_boxplot(position = position_dodge(0.8),lwd = 1)+
    geom_line(aes(group = Animal), alpha = 0.5)+
    geom_point(aes(group = Animal), size = 3)+
    facet_wrap(~ Opsin)+
    theme_bw()+
    #stat_summary(fun.y=median, geom="point", shape=18, size=3, color="blue", position=position_dodge(0.5))+
    labs(title = "Opto Hot Plate nShakes before EP", x=element_blank(), y="number of shakes")+
    theme(plot.title = element_text(
      hjust = 0.5,
      vjust = 0.5,
      face = 'bold'))
  ggsave(paste(output_dir,"Opto_HotPlate_nShakesBeforeEP_DotPlot-Boxplot_",qq,".tiff"), height=5, width=7)
  
  gg_color <-  ggplot(long_qqData, aes(x=`Condition`, y=`Value`, fill=Condition)) +
    #geom_boxplot(stat = "summary", fun = "mean", position = "dodge")+
    stat_boxplot(geom="errorbar", position = position_dodge(0.8), size = 1, width = 0.5)+
    geom_boxplot(position = position_dodge(0.8),lwd = 1)+
    geom_line(aes(group = Animal, color = Animal), lwd = 1)+
    geom_point(aes(group = Animal), size = 3)+
    facet_wrap(~ Opsin)+
    theme_bw()+
    #stat_summary(fun.y=median, geom="point", shape=18, size=3, color="blue", position=position_dodge(0.5))+
    labs(title = "Opto Hot Plate nShakes before EP", x=element_blank(), y="number of shakes")+
    theme(plot.title = element_text(
      hjust = 0.5,
      vjust = 0.5,
      face = 'bold'))
  ggsave(paste(output_dir,"ColorCoded_Opto_HotPlate_nShakesBeforeEP_DotPlot-Boxplot_",qq,".tiff"), height=5, width=7)
  
}
#reorder columns in Results DF and then save it as xlsx
Results <- Results[,c(1,7,2,3,4,5,6,8)]
ResNam=paste(output_dir,"Opto_HotPlate_nShakesBeforeEP_tTest_overview",".xlsx",sep="")
write.xlsx(Results,ResNam)





