################################################################
###################### plotting DREADD data ####################
################################################################

#Important: This script should be run with the working directory
#set to the folder "1.5IQR_Outlier_Removed_Files" that was created 
#with the script "01_DREADD_1.5IQR-outlier_removal.R".

rm(list=ls())

library(pastecs)
library(car)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(robustbase)
library(reshape2)
library(chron)
library(lubridate)
library(tidyr)

options(scipen=999) #Turns off scientific writing style
filedir= getwd()

#set all reference files and variables
DREADD_info=read.xlsx("DREADD_Additional_Information.xlsx")
AllAnimals=unique(DREADD_info[,3])

#create an empty DF for later storage of time filtered Data of ALL DREADD animals (BoxplotDF)
AllDREADD_Data=data.frame(matrix(nrow=0, ncol= 24))
colnames(AllDREADD_Data)= c("Animal","DREADD_receptor","Condition","ElapsedTime","MinutesPostInjection","PIF","PEF","TV","MV","BPM","IT","ET","TT","Pause","PEnh","TVadj","MVadj","IF50","EF50","E/IF50","TVadjPerGram(uL)","MVadjPerGram(mL)","PIFadj","PEFadj")

#create loop to run through all the DREADD Animal files
for (i in AllAnimals){
  print(i)
  PerAnimal=DREADD_info[DREADD_info[,3]==i,]

  #create variables for injection time in minutes
  Conditions=c("BL1","BL2","BL3","Exp1","Exp2","Exp3","Clo1","Clo2","Clo3")
  for (ii in Conditions){
    do.call("=", list(paste("InjMin", ii, sep = ""), PerAnimal[PerAnimal[, 4]==ii, 5]*60*24))
  }
  
  #load data and seperate BL and DF
  Muster=paste(i,".xlsx$",sep="")
  readthat=list.files(filedir, pattern=Muster)
  myDataOrig=read.xlsx(paste(filedir,readthat, sep=""))

  BL_DF= myDataOrig[myDataOrig$Condition=="BL1" | myDataOrig$Condition=="BL2" | myDataOrig$Condition=="BL3",]
  Exp_DF= myDataOrig[myDataOrig$Condition=="Exp1" | myDataOrig$Condition=="Exp2" | myDataOrig$Condition=="Exp3",]
  Clo_DF= myDataOrig[myDataOrig$Condition=="Clo1" | myDataOrig$Condition=="Clo2" | myDataOrig$Condition=="Clo3",]
  
  #convert realTime to only time of day
  a= strsplit(BL_DF$RealTime," ")
  b=1:length(a)
  for (s in 1:length(b)){
    b[s]=a[[s]][2]
  }
  BL_DF$RealTime= b
  
  c= strsplit(Exp_DF$RealTime," ")
  d=1:length(c)
  for (s in 1:length(d)){
    d[s]=c[[s]][2]
  }
  Exp_DF$RealTime= d
  
  e= strsplit(Clo_DF$RealTime," ")
  f=1:length(e)
  for (s in 1:length(f)){
    f[s]=e[[s]][2]
  }
  Clo_DF$RealTime= f
  
  #convert Realtime to minutes
  BL_DF$RealTime=sapply(strsplit(BL_DF$RealTime,":"),
                        function(Q) {
                          Q= as.numeric(Q)
                          Q[1]*60+Q[2]+Q[3]/60
                        }
  )
  Exp_DF$RealTime=sapply(strsplit(Exp_DF$RealTime,":"),
                         function(Q) {
                           Q=as.numeric(Q)
                           Q[1]*60+Q[2]+Q[3]/60
                         }
  )
  Clo_DF$RealTime=sapply(strsplit(Clo_DF$RealTime,":"),
                         function(Q) {
                           Q= as.numeric(Q)
                           Q[1]*60+Q[2]+Q[3]/60
                         }
  )
  
  
  #convert all columns to numeric class
  for(k in 5:ncol(BL_DF)){
    BL_DF[,k]=as.numeric(BL_DF[,k])
    Exp_DF[,k]=as.numeric(Exp_DF[,k])
    Clo_DF[,k]=as.numeric(Clo_DF[,k])
  }
  
  #convert RealTime colname to minutes post injection
  colnames(BL_DF)[5]="MinutesPostInjection"
  colnames(Exp_DF)[5]="MinutesPostInjection"
  colnames(Clo_DF)[5]="MinutesPostInjection"
  
  
  # combine BL_DF, Exp_DF and Clo_DF and change VALUES RealTime to Minutes post injection
  DF=bind_rows(BL_DF,Exp_DF,Clo_DF)
  
  for(jj in 1:nrow(DF)) {
    if (DF[jj, 3] == "BL1") {
      DF[jj, 5] = DF[jj, 5] - InjMinBL1
    } else if (DF[jj, 3] == "BL2") {
      DF[jj, 5] = DF[jj, 5] - InjMinBL2
    } else if (DF[jj, 3] == "BL3") {
      DF[jj, 5] = DF[jj, 5] - InjMinBL3
    } else if (DF[jj, 3] == "Exp1") {
      DF[jj, 5] = DF[jj, 5] - InjMinExp1
    } else if (DF[jj, 3] == "Exp2") {
      DF[jj, 5] = DF[jj, 5] - InjMinExp2
    } else if (DF[jj, 3] == "Exp3") {
      DF[jj, 5] = DF[jj, 5] - InjMinExp3
    } else if (DF[jj, 3] == "Clo1") {
      DF[jj, 5] = DF[jj, 5] - InjMinClo1
    } else if (DF[jj, 3] == "Clo2") {
      DF[jj, 5] = DF[jj, 5] - InjMinClo2
    } else if (DF[jj, 3] == "Clo3") {
      DF[jj, 5] = DF[jj, 5] - InjMinClo3
    }
  }
  
  #filter according to condition and only between 50' and 2h after injection for BL and CNO
  # for Clo: filter to keep data between 20' and 1h30'
  BL1_DF=DF[DF[,3]=="BL1" & DF[,5]>50 & DF[,5]<=120,]
  BL2_DF=DF[DF[,3]=="BL2" & DF[,5]>50 & DF[,5]<=120,]
  BL3_DF=DF[DF[,3]=="BL3" & DF[,5]>50 & DF[,5]<=120,]
  Exp1_DF=DF[DF[,3]=="Exp1" & DF[,5]>50 & DF[,5]<=120,]
  Exp2_DF=DF[DF[,3]=="Exp2" & DF[,5]>50 & DF[,5]<=120,]
  Exp3_DF=DF[DF[,3]=="Exp3" & DF[,5]>50 & DF[,5]<=120,]
  Clo1_DF=DF[DF[,3]=="Clo1" & DF[,5]>20 & DF[,5]<=90,]
  Clo2_DF=DF[DF[,3]=="Clo2" & DF[,5]>20 & DF[,5]<=90,]
  Clo3_DF=DF[DF[,3]=="Clo3" & DF[,5]>20 & DF[,5]<=90,]
  
  # join the condition separated DF into one DF again
  DF=rbind(BL1_DF,BL2_DF,BL3_DF,Exp1_DF,Exp2_DF,Exp3_DF,Clo1_DF,Clo2_DF,Clo3_DF)
  
  #write the BoxplotDF into the AllDREADD_data
  AllDREADD_Data=rbind(AllDREADD_Data,DF)
}

#save the AllDREADD_data to the working directory
write.xlsx(AllDREADD_Data,"AllDREADD_NoLog_1.5IQR-rem_timefilt.xlsx")



#######################################################################################
######## Plot lineplots ###############################################################
#######################################################################################

#read data
myData = AllDREADD_Data

# add a column with the super-conditions
myData$SupCon=NA

for (i in 1:nrow(myData)) {
  if (myData$Condition[i] == "BL1" |
      myData$Condition[i] == "BL2" | myData$Condition[i] == "BL3") {
    myData$SupCon[i] = "BL"
  } else if (myData$Condition[i] == "Exp1" |
             myData$Condition[i] == "Exp2" | myData$Condition[i] == "Exp3") {
    myData$SupCon[i] = "CNO"
  } else if(myData$Condition[i] == "Clo1" |
            myData$Condition[i] == "Clo2" | myData$Condition[i] == "Clo3") {
    myData$SupCon[i] = "Clo"
  }
}

# create a DF with 1' averages
#first convert DF columns 4:end to numeric
myData[,5:24]=lapply(5:24,
                     function(x)
                       as.numeric(myData[[x]]))
#convert to Matrix
DF_mat=as.matrix(myData[,5:24])

AvgMat=matrix(numeric(0),nrow(DF_mat)/12,ncol(DF_mat))
for (l in seq(1, nrow(myData), by = 12)) {
  AvgMat[ceiling(l / 12),] = colMeans(DF_mat[l:(l +11),],na.rm=TRUE)
}
AvgMat[,1]=ceiling(AvgMat[,1])

#add Condition column to Matrix
Con_DREADD=matrix(numeric(0),nrow(AvgMat),3)
for (ll in seq(1, nrow(myData), by = 12)) {
  Con_DREADD[ceiling(ll/12),1]=myData[ll,2]
  Con_DREADD[ceiling(ll/12),2]=myData[ll,3]
  Con_DREADD[ceiling(ll/12),3]=myData[ll,25]
}

DF=data.frame(Con_DREADD,AvgMat)
colnames(DF)=c("DREADD_receptor","Condition","SupCon","MinutesPostInjection","PIF","PEF","TV","MV","BPM","IT","ET","TT","Pause","PEnh","TVadj","MVadj","IF50","EF50","E_IF50","TVadjPerGram(uL)","MVadjPerGram(ml)","PIFadj","PEFadj")

#crete DF to plot the data relative to expected peak effect
shiftData=DF
shiftData[shiftData$SupCon=="Clo",4]=shiftData[shiftData$SupCon=="Clo",4]-30
shiftData[shiftData$SupCon=="CNO" | shiftData$SupCon=="BL" ,4]= shiftData[shiftData$SupCon=="CNO" | shiftData$SupCon=="BL" ,4]-60

#transform sortData to longformat
longDF=melt(DF, id.vars=c("DREADD_receptor","Condition","SupCon", "MinutesPostInjection"), variable.name = "Parameter", value.name = "Value")
longShift=melt(shiftData, id.vars=c("DREADD_receptor","Condition","SupCon", "MinutesPostInjection"), variable.name = "Parameter", value.name = "Value")

#calculate mean and SE per Opsin and Condition
SEDF=longDF %>%
  group_by(DREADD_receptor, Parameter, SupCon, MinutesPostInjection) %>%
  summarise(new=list(mean_se(Value))) %>%
  unnest(new)

SEShift=longShift %>%
  group_by(DREADD_receptor, Parameter, SupCon, MinutesPostInjection) %>%
  summarise(new=list(mean_se(Value))) %>%
  unnest(new)

# plot linePlot of sortData facetted by parameter
groupList=unique(longDF$DREADD_receptor)

for (jj in groupList){
  print(jj)
  groupDF=SEDF[SEDF$DREADD_receptor==jj,]
  groupShift=SEShift[SEShift$DREADD_receptor==jj,]
  
  #plot minutes post Injection
  g = ggplot(groupDF, aes(x = MinutesPostInjection, y = y),se=TRUE) +
    geom_line(aes(col = SupCon), size = 1, alpha=1) +
    geom_ribbon(aes(fill= SupCon, ymin=ymin, ymax=ymax), alpha=0.2, show.legend = FALSE)+
    scale_color_hue(l = 50, c = 150) +
    theme_bw() +
    labs(title = paste(jj,"condition averaged line plots",sep=" "), x = "Time after injection (min)") +
    scale_color_discrete(breaks=c("BL","Clo","CNO"),labels=c("Saline","Clozapine","CNO"))+
    scale_x_continuous(breaks = seq(20, 120, 10)) +
    facet_grid(rows=vars(Parameter),scales="free_y")+
    theme(plot.title = element_text(
      hjust = 0.5,
      vjust = 0.5,
      face = 'bold',
      size=30))+
    theme(axis.title.y=element_blank(),
          axis.title.x=element_text(size=20))+
    theme(legend.key.size=unit(3.5, "line"), legend.text=element_text(size=20), legend.title = element_blank())
  #print(g)
  
  ggsave(paste("DREADD_NoLog_ribbon-linePlot_",jj,".tiff",sep=""),width=15, height=25)
  
  #plot minutes relative to expected peak effect
  gg = ggplot(groupShift, aes(x = MinutesPostInjection, y = y),se=TRUE) +
    geom_line(aes(col = SupCon), size = 1, alpha=1) +
    geom_ribbon(aes(fill=SupCon, ymin=ymin, ymax=ymax), alpha=0.2, show.legend = FALSE)+
    scale_color_hue(l = 50, c = 150) +
    theme_bw() +
    labs(title = paste(jj,"condition averaged line plots",sep=" "), x = "Time relative to expected peak effect (min)") +
    scale_color_discrete(breaks=c("BL","Clo","CNO"),labels=c("Saline","Clozapine","CNO"))+
    scale_x_continuous(breaks = seq(-10, 60, 10)) +
    facet_grid(rows=vars(Parameter),scales="free_y")+
    theme(plot.title = element_text(
      hjust = 0.5,
      vjust = 0.5,
      face = 'bold',
      size=30))+
    theme(axis.title.y=element_blank(),
          axis.title.x=element_text(size=20))+
    theme(legend.key.size=unit(3.5, "line"), legend.text=element_text(size=20), legend.title = element_blank())
  #print(gg)
  
  ggsave(paste("DREADD_NoLog_ribbon-ShiftLinePlot_",jj,".tiff",sep=""),width=15, height=25)
  
  #plot short version of shiftLinePlot
  gFiltDF=groupShift[groupShift$Parameter== "BPM" | groupShift$Parameter== "IT" | groupShift$Parameter== "ET" | groupShift$Parameter== "TT" | groupShift$Parameter== "TVadjPerGram(uL)" | groupShift$Parameter== "MVadjPerGram(ml)" | groupShift$Parameter== "PIFadj" | groupShift$Parameter== "PEFadj",]
  gFilt = ggplot(gFiltDF, aes(x = MinutesPostInjection, y = y),se=TRUE) +
    geom_line(aes(col = SupCon), size = 1, alpha=1) +
    geom_ribbon(aes(fill=SupCon, ymin=ymin, ymax=ymax), alpha=0.2, show.legend = FALSE)+
    scale_color_hue(l = 50, c = 150) +
    theme_bw() +
    labs(title = paste(jj,"condition averaged line plots",sep=" "), x = "Time relative to expected peak effect (min)") +
    scale_color_discrete(breaks=c("BL","Clo","CNO"),labels=c("Saline","Clozapine","CNO"))+
    scale_x_continuous(breaks = seq(-10, 60, 10)) +
    facet_grid(rows=vars(Parameter),scales="free_y")+
    theme(plot.title = element_text(
      hjust = 0.5,
      vjust = 0.5,
      face = 'bold',
      size=30))+
    theme(axis.title.y=element_blank(),
          axis.title.x=element_text(size=20))+
    theme(legend.key.size=unit(3, "line"), legend.text=element_text(size=15), legend.title = element_blank())
  #print(gFilt)
  
  ggsave(paste("DREADD_NoLog_ribbon-shortShiftLinePlot_",jj,".tiff",sep=""),width=15, height=11)
}


#plot faceted short version of shiftLinePlot
gAggDF=SEShift[SEShift$Parameter== "BPM" | SEShift$Parameter== "IT" | SEShift$Parameter== "ET" | SEShift$Parameter== "TT" | SEShift$Parameter== "TVadjPerGram(uL)" | SEShift$Parameter== "MVadjPerGram(ml)" | SEShift$Parameter== "PIFadj" | SEShift$Parameter== "PEFadj",]
gAgg = ggplot(gAggDF, aes(x = MinutesPostInjection, y = y),se=TRUE) +
  geom_line(aes(col = SupCon), size = 1, alpha=1) +
  geom_ribbon(aes(fill=SupCon, ymin=ymin, ymax=ymax), alpha=0.2, show.legend = FALSE)+
  scale_color_hue(l = 50, c = 150) +
  theme_bw() +
  labs(title = "DREADD condition-averaged line plots", x = "Time relative to expected peak effect (min)") +
  scale_color_discrete(breaks=c("BL","Clo","CNO"),labels=c("Saline","Clozapine","CNO"),name="Condition")+
  scale_x_continuous(breaks = seq(-10, 60, 10)) +
  facet_grid(Parameter~DREADD_receptor,scales="free_y")+
  theme(plot.title = element_text(
    hjust = 0.5,
    vjust = 0.5,
    face = 'bold',
    size=30))+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_text(size=20))+
  theme(legend.key.size=unit(3, "line"), legend.text=element_text(size=15), legend.title = element_text(size=20, face="bold"))
#print(gAgg)
ggsave("DREADD_NoLog_ribbon-shortAggShiftLinePlot.tiff",width=15, height=11)


########################################################################
############ plot violin plots #########################################
########################################################################

myData = AllDREADD_Data

# name all Conditions either BL, Exp or Clo
for (i in 1:length(myData$Condition)){
  if (myData$Condition[i]=="BL1" | myData$Condition[i]=="BL2" | myData$Condition[i]=="BL3"){
    myData$Condition[i]="BL"
  }else if (myData$Condition[i]=="Exp1" |myData$Condition[i]=="Exp2" |myData$Condition[i]=="Exp3" |myData$Condition[i]=="Exp4" |myData$Condition[i]=="Exp5" |myData$Condition[i]=="Exp6"){
    myData$Condition[i]="Exp"
  }else{myData$Condition[i]="Clo"}
}

#transform DF into long format
longData=melt(myData, id.vars= c("Animal","DREADD_receptor", "Condition", "ElapsedTime", "MinutesPostInjection"), variable.name="Parameter", value.name = "Value") 


# violin-plot an filtered overview across all groups and conditions
gFiltDF=longData[longData$Parameter== "BPM" | longData$Parameter== "IT" | longData$Parameter== "ET" | longData$Parameter== "TT" | longData$Parameter== "TVadjPerGram(uL)" | longData$Parameter== "MVadjPerGram(mL)" | longData$Parameter== "PIFadj" | longData$Parameter== "PEFadj",]
gFilt= ggplot(gFiltDF, aes(y=Value, x=Condition, fill=Condition))+
  geom_violin() +
  stat_boxplot(geom="errorbar", width=0.08)+
  geom_boxplot(width=0.1,lwd=0.3,outlier.size=0) + theme_minimal()+
  theme_bw() +
  theme(axis.text.x=element_blank(),axis.title.y = element_blank(), axis.ticks=element_blank())+
  scale_fill_discrete(breaks=c("BL","Clo","Exp"),labels=c("Saline","Clozapine","CNO"))+
  labs(title = "DREADD raw data group comparison") +
  facet_grid(Parameter~DREADD_receptor,scales="free")+
  theme(plot.title = element_text(
    hjust = 0.5,
    vjust = 0.5,
    face = 'bold',
    size = 30))+
  theme(legend.key.size=unit(3, "line"), legend.text=element_text(size=15), legend.title = element_text(size=15, face = "bold"))

ggsave("unavg_DREADD_NoLog_shortOverview_Violin_allAnimals.tiff",width=10, height=11)
