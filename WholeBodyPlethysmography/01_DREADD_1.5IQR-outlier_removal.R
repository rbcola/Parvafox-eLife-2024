#####################################################
########  Outlier Removal (+/-1.5*IQR)   ############
#####################################################

#This script should be run in the directory contianing 
#the raw files (starting with "AllCond_Animal...")
#deposited in the companion file repositry.

rm(list=ls())
library(openxlsx)
library(rcompanion)

#create output directory to store the filtered data
output_dir=paste(getwd(),"1.5IQR_Outlier_Removed_Files/", sep="/")
ifelse(!dir.exists(output_dir), dir.create(output_dir),FALSE)

#create outlier removal function
removeOutlier =function(x){ # x has to be a column vector
  qnt = quantile(x, probs = c(.25, .75), na.rm= TRUE) #get 1st and 3rd quantile
  step = 1.5 * IQR(x, na.rm = TRUE) #calculate 1.5x IQR step
  y=x
  y[x < qnt[1]-step | x > qnt[2] + step] = NA  #convert outliers to NA
  return(y)
}

#list all DREADD files in which outliers should be removed
allFiles=list.files(getwd(), pattern="AllCond_Animal*")

#loop through all the files in allFiles
for (i in allFiles){
  print(i)
  myData=read.xlsx(i)
  data=myData[,c(6:13,15:18,20:22,38:39,41:42)]

  for (j in 1:ncol(data)){
    data[,j]=removeOutlier(data[,j])
    #print(boxplot(data[,j]))
  }

  #add the missing factor columns to the dataFilt
  dataFilt=cbind(myData[,c(44:45,43,1,2)],data)

  #save data file as xlsx
  write.xlsx(dataFilt, paste(output_dir,"1.5IQR-outRem_",i,sep=""))
}