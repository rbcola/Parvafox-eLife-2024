#####################################################################
# Script to analyze Deeplabcut outputs from open field data #########
# Reto B. Cola ######################################################
# 
# Run this script in the same directory, that contains the companion 
# raw tracking data downloaded from the data repository. Also make sure
# that the script containing the functions  
# ("DLCAnalyzer_Functions_final_RBCmodified.R")is in the same directory. 

rm(list=ls())


library(sp)         #tested with v1.3-2
library(imputeTS)   #tested with v3.0
library(ggplot2)    #tested with v3.1.0
library(ggmap)      #tested with v3.0.0
library(data.table) #tested with v1.12.8
library(cowplot)    #tested with v0.9.4
library(corrplot)   #tested with v0.84
library(readxl)
library(hms)
library(lubridate)
library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)
library(purrr)
library(effsize)
library(openxlsx)
library(gridExtra)
library(scico)


#load the DLCAnalyzer functions
source("DLCAnalyzer_Functions_final_RBCmodified.R")

##############################################################################
#define variables manually
ExpGroup <- "DREADD" #choose between "Opto" or "DREADD"
input_folder <- "H:/S3IT_Linux/DLC_Dec23/" #NEEDS TO END WITH "/"
HelperFilePath <- "H:/S3IT_Linux/DLC_Dec23/22-12-22_OpenField_Starting_times.xlsx" #helper file with FrameTimes of Stim etc.
OptoGroupFilePath <- "H:/S3IT_Linux/DLC_Dec23/22-12-22_OptoAnimalGroupAllocation_HelperFile.xlsx"
DREADDGroupFilePath <- "H:/S3IT_Linux/DLC_Dec23/22-12-22_DREADDAnimal_HelperFile.xlsx"
VideoFrameRate <- 25 #frames/s
Opto_BinData <- "H:/S3IT_Linux/DLC_Dec23/Opto_BinData.csv"
filtered <- TRUE

OFT_movement_cutoff <- 1
OFT_integration_period <- 13
Freezing_cutoff <- 0.01
Freezing_integration_period <- 10
##############################################################################

#Create timestamp
DateTime <- Sys.time() %>% unlist() %>% 
  str_replace(pattern= " ", replacement = "_") %>% 
  str_replace_all(pattern= ":", replacement = "-")

#create output folder name
output_folder <- paste0(input_folder,DateTime,"_",ExpGroup,"/")
if (!dir.exists(output_folder)){
  dir.create(output_folder)
}

# save the workspace to reconstruct the analysis if necessary
save.image(file = paste0(output_folder, DateTime,"_Workspace.RData"))



#load the helper file with starting times of stimulation
HelperFile <- read_excel(HelperFilePath, 
                         col_types = c("text", "date", "date", "text", "text", "text", "numeric")) %>% 
              drop_na(c(2,3)) %>% 
              select(1:4)

#The mm:ss format is read incorrectly. Therefore i need to modify
HelperFile <- HelperFile %>% 
  separate(StartHab, c(NA, "StartHab"), " ") %>% 
  mutate(StartHab = StartHab %>% 
           str_remove(":00") %>% 
           as.POSIXct(format="%M:%S"))
         
HelperFile <- HelperFile %>% 
  separate(StartRec, c(NA, "StartRec"), " ") %>% 
  mutate(StartRec = StartRec %>% 
           str_remove(":00") %>% 
           as.POSIXct(format="%M:%S"))

#Split Recording into AnimalNr & Condition %>% convert the Start Rec time into FrameCount
HelperFile <- HelperFile %>% separate(Recording, c("VideoName_AnimalID","Condition"), " ") %>% 
  mutate(StimStartFrame = minute(StartRec)*60*VideoFrameRate + second(StartRec)*VideoFrameRate)


#Load OptoGroupFile and get rid of the last 4 columns. Load DREADDGroupFile and remove WBPsession column
OptoGroupFile <- read_excel(OptoGroupFilePath) %>% 
  select(-c((ncol(.)-3):ncol(.)))

DREADDGroupFile <- read_excel(DREADDGroupFilePath) %>% 
  select(-WBPsession)


###############################################################################
#Perform DLCAnyalysis by groups
###############################################################################

print(ExpGroup)

HelperFileModality <- HelperFile %>% filter(Modality == ExpGroup) 


if( ExpGroup == "Opto"){
  # leftjoin the OptoGroupFile and HelperFileModality and remove NA containing rows
  DLCHelperFile <- left_join(HelperFileModality, OptoGroupFile, by = "VideoName_AnimalID") %>% 
    na.omit()
  # Only keep LED1 Condition for OptoTrials, because Marco and Diana only recorded animals once
  # Also only Keep Foxb1-Cre genotypes
  DLCHelperFile <- DLCHelperFile %>% 
    filter(Condition == "LED1") %>% 
    filter(genotype == "Foxb1-Cre")
  
  #create list of files for Optogenetics modality
  files1 <- list.files(input_folder, pattern= "^Maus_")
  files2 <- list.files(input_folder, pattern= "^OF_top_Opto_")
  if(filtered == TRUE){
    files3 <- list.files(input_folder, pattern= "_filtered.csv$")
  }else{
    files3 <- list.files(input_folder, pattern= "[0-9].csv$")
  }
  files <- intersect(union(files1,files2),files3)
  
  #add LED1 condition to Marco's and Diana's files 
  for ( s in 1:length(files)){
    if (grepl("MAH0", files[s], fixed = TRUE)){
      pre <- files[s] %>% str_split("(?<=.)(?=DLC)") %>% map(1) %>% unlist()# split string with lookahead regex
      post <- files[s] %>% str_split("(?<=.)(?=DLC)") %>% map(2) %>% unlist()# split string with lookahead regex
      files[s] <- paste0(pre, "_LED1", post)
    }
    if (startsWith(files[s],"Maus_35")){
      files[s] <- files[s] %>% str_replace("LED1_","")
    }
  }
  
  #Only keep files contained in the helperfile (i.e. only LED1 files)
  files <- files[files %>% grepl(paste0(unique(DLCHelperFile$Condition),"DLC"),.)] 
  
  StartPattern <- unique(DLCHelperFile$VideoName_AnimalID)
  files <- files[grepl(paste(StartPattern, collapse="|"), files)] #filter to keep only files of mice contained in HelperFile
  
  #remove "_LED1" again from list to being able to call the csv files further below
  for ( s in 1:length(files)){
    if (grepl("MAH0", files[s], fixed = TRUE)){
      files[s] <- files[s] %>% str_replace("_LED1","")
    }
    if(startsWith(files[s],"Maus_35")){
      pre <- files[s] %>% str_split("(?<=.)(?=MAH)") %>% map(1) %>% unlist()# split string with lookahead regex
      post <- files[s] %>% str_split("(?<=.)(?=MAH)") %>% map(2) %>% unlist()
      files[s] <- paste0(pre,"LED1_",post)
    }
  }
  

  #define your DLCAnalyzer processing pipelines
  OFT_pipeline <- function(path){
    
    if(path %>% str_replace(input_folder,"") %>% str_detect("Maus_\\d{2,3}[^\\w]")){
      AnimalID <- path %>% str_replace(input_folder,"") %>% str_split("_MAH") %>% map(1) %>% unlist()
    }else if(path %>% str_replace(input_folder,"") %>% str_detect("Maus_35\\w")){
      AnimalID <- path %>% str_replace(input_folder,"") %>% str_split("_LED1") %>% map(1) %>% unlist()
    }else if(path %>% str_replace(input_folder,"") %>% startsWith("OF_top")){
      AnimalID <- path %>% str_replace(input_folder,"") %>% str_split("-\\w{2,3}\\dDLC_") %>% map(1) %>% unlist()
    }
    
    StartStim <- DLCHelperFile %>% filter(VideoName_AnimalID == AnimalID) %>% 
      select(StimStartFrame) %>% pull()
    CropFrom <- DLCHelperFile %>% filter(VideoName_AnimalID == AnimalID) %>% 
      select(StimStartFrame) %>% pull() -(3*60*VideoFrameRate)-25 #3min BL +1s room for error
    CropTo <- DLCHelperFile %>% filter(VideoName_AnimalID == AnimalID) %>% 
      select(StimStartFrame) %>% pull() +(3*60*VideoFrameRate)+24#3min Stim +1s room for error (the first frame of Stim is inclusive, hence +24frames)
    
    #create bining info data frame
    bindat <- data.frame(bin = c("B1","B2"),
                         from = c(CropFrom, StartStim+25),
                         to = c(StartStim-26, CropTo))
    
    
    Tracking <- ReadDLCDataFromCSV(path, fps = VideoFrameRate)
    Tracking <- CutTrackingData(Tracking, keep.frames = c(CropFrom:CropTo))
    Tracking <- CleanTrackingData(Tracking, likelihoodcutoff = 0.95)
    Tracking <- CalibrateTrackingData(Tracking, "area",in.metric = 40*40, points = c("ArenaTL","ArenaTR","ArenaBR","ArenaBL"))
    Tracking <- AddBinData(Tracking, bindat = bindat, unit = "frame") # 3' BL, 2x1'' room for error, 3' Stim
    Tracking$seconds <- Tracking$seconds - Tracking$seconds[1]
    Tracking$point.info$PointType <- c("Arena","Arena","Arena","Arena","Mouse","Mouse","Mouse","Mouse","Mouse") #add PointType for FST analysis
    Tracking <- AddOFTZones(Tracking, scale_center = 0.5,scale_periphery  = 0.8 ,scale_corners = 0.4, points = c("ArenaTL","ArenaTR","ArenaBR","ArenaBL"))
    Tracking <- OFTAnalysis(Tracking, movement_cutoff = OFT_movement_cutoff, integration_period = OFT_integration_period, points = c("Occiput"))
    return(Tracking)
  }
  
  Freezing_pipeline <- function(path){
    
    if(path %>% str_replace(input_folder,"") %>% str_detect("Maus_\\d{2,3}[^\\w]")){
      AnimalID <- path %>% str_replace(input_folder,"") %>% str_split("_MAH") %>% map(1) %>% unlist()
    }else if(path %>% str_replace(input_folder,"") %>% str_detect("Maus_35\\w")){
      AnimalID <- path %>% str_replace(input_folder,"") %>% str_split("_LED1") %>% map(1) %>% unlist()
    }else if(path %>% str_replace(input_folder,"") %>% startsWith("OF_top")){
      AnimalID <- path %>% str_replace(input_folder,"") %>% str_split("-\\w{2,3}\\dDLC_") %>% map(1) %>% unlist()
    }
    
    StartStim <- DLCHelperFile %>% filter(VideoName_AnimalID == AnimalID) %>% 
      select(StimStartFrame) %>% pull()
    CropFrom <- DLCHelperFile %>% filter(VideoName_AnimalID == AnimalID) %>% 
      select(StimStartFrame) %>% pull() -(3*60*VideoFrameRate)-25 #3min BL +1s room for error
    CropTo <- DLCHelperFile %>% filter(VideoName_AnimalID == AnimalID) %>% 
      select(StimStartFrame) %>% pull() +(3*60*VideoFrameRate)+24#3min Stim +1s room for error (the first frame of Stim is inclusive, hence +24frames)
    
    #create bining info data frame
    bindat <- data.frame(bin = c("B1","B2"), 
                         from = c(CropFrom, StartStim+25), 
                         to = c(StartStim-26, CropTo))
    
    Tracking <- ReadDLCDataFromCSV(path, fps = VideoFrameRate)
    Tracking <- CutTrackingData(Tracking, keep.frames = c(CropFrom:CropTo))
    Tracking <- CleanTrackingData(Tracking, likelihoodcutoff = 0.95)
    Tracking <- CalibrateTrackingData(Tracking, "area",in.metric = 40*40, points = c("ArenaTL","ArenaTR","ArenaBR","ArenaBL"))
    Tracking <- AddBinData(Tracking, bindat = bindat, unit = "frame") # 3' BL, 2'' room for error, 3' Stim
    Tracking$point.info$PointType <- c("Arena","Arena","Arena","Arena","Mouse","Mouse","Mouse","Mouse","Mouse") #add PointType for FST analysis
    Tracking <- FreezingAnalysis(Tracking, cutoff_freezing = Freezing_cutoff, integration_period = Freezing_integration_period, Object = "Mouse", points = "Occiput")
    return(Tracking)
  }  
 
  
  
   
}




if(ExpGroup == "DREADD"){
  
  DLCHelperFile <- left_join(HelperFileModality, DREADDGroupFile, by = "VideoName_AnimalID")
  
  
  files1 <- list.files(input_folder, pattern= "^OF_top_DREADD_")
  if(filtered == TRUE){
    files2 <- list.files(input_folder, pattern= "_filtered.csv$")
    }else{
      files2 <- list.files(input_folder, pattern= "[0-9].csv$")
    }
  
  files <- intersect(files1, files2)

  #define your DLCAnalyzer processing pipeline
  #AnimalID <- files[1] %>% str_split("-\\w{2,3}\\dDLC_") %>% map(1) %>% unlist() #this will be overwritten. It's just needed to initialize the variable
  #Cond <- files[1] %>% str_split("DLC") %>% map(1) %>% str_split("-") %>% map(2) %>% unlist() 
  OFT_pipeline <- function(path){

    AnimalID <- path %>% str_replace(input_folder,"") %>% str_split("-\\w{2,3}\\dDLC_") %>% map(1) %>% unlist()
    Cond <- path %>% str_replace(input_folder,"") %>% str_split("DLC") %>% map(1) %>% str_split("-") %>% map(2) %>% unlist() 
    
    CropFrom <- DLCHelperFile %>%
      filter(VideoName_AnimalID == AnimalID & Condition == Cond) %>%
      select(StimStartFrame) %>% pull()
    CropTo <- DLCHelperFile %>%
      filter(VideoName_AnimalID == AnimalID & Condition == Cond) %>%
      select(StimStartFrame) %>% pull() +(5*60*VideoFrameRate)-1 #5min recording after 5min habituation

    Tracking <- ReadDLCDataFromCSV(path, fps = VideoFrameRate)
    Tracking <- CutTrackingData(Tracking, keep.frames = c(CropFrom:CropTo))
    Tracking <- CleanTrackingData(Tracking, likelihoodcutoff = 0.95)
    Tracking <- CalibrateTrackingData(Tracking, "area",in.metric = 40*40, points = c("ArenaTL","ArenaTR","ArenaBR","ArenaBL"))
    Tracking$seconds <- Tracking$seconds - Tracking$seconds[1]
    Tracking$point.info$PointType <- c("Arena","Arena","Arena","Arena","Mouse","Mouse","Mouse","Mouse","Mouse") #add PointType for FST analysis
    Tracking <- AddOFTZones(Tracking, scale_center = 0.5,scale_periphery  = 0.8 ,scale_corners = 0.4, points = c("ArenaTL","ArenaTR","ArenaBR","ArenaBL"))
    Tracking <- OFTAnalysis(Tracking, movement_cutoff = OFT_movement_cutoff, integration_period = OFT_integration_period, points = c("Occiput"))
    return(Tracking)
  }

  Freezing_pipeline <- function(path){
    
    AnimalID <- path %>% str_replace(input_folder,"") %>% str_split("-\\w{2,3}\\dDLC_") %>% map(1) %>% unlist()
    Cond <- path %>% str_replace(input_folder,"") %>% str_split("DLC") %>% map(1) %>% str_split("-") %>% map(2) %>% unlist() 
    
    CropFrom <- DLCHelperFile %>%
      filter(VideoName_AnimalID == AnimalID & Condition == Cond) %>%
      select(StimStartFrame) %>% pull()
    CropTo <- DLCHelperFile %>%
      filter(VideoName_AnimalID == AnimalID & Condition == Cond) %>%
      select(StimStartFrame) %>% pull() +(5*60*VideoFrameRate)-1 #5min recording after 5min habituation
    
    Tracking <- ReadDLCDataFromCSV(path, fps = VideoFrameRate)
    Tracking <- CutTrackingData(Tracking, keep.frames = c(CropFrom:CropTo))
    Tracking <- CleanTrackingData(Tracking, likelihoodcutoff = 0.95)
    Tracking <- CalibrateTrackingData(Tracking, "area",in.metric = 40*40, points = c("ArenaTL","ArenaTR","ArenaBR","ArenaBL"))
    Tracking$seconds <- Tracking$seconds - Tracking$seconds[1]
    Tracking$point.info$PointType <- c("Arena","Arena","Arena","Arena","Mouse","Mouse","Mouse","Mouse","Mouse") #add PointType for FST analysis
    Tracking <- FreezingAnalysis(Tracking, cutoff_freezing = Freezing_cutoff, integration_period = Freezing_integration_period, Object = "Mouse", points = "Occiput")
    return(Tracking)
  }    
  
}



#after the pipeline is defined, we can execute it for all files and 
#combine them into a list of Tracking objects
OFT_TrackingAll <- RunPipeline(files,input_folder,FUN = OFT_pipeline)
saveRDS(OFT_TrackingAll, file = paste0(output_folder, DateTime,"_",ExpGroup,"_OFT_TrackingAll_CutOff-",OFT_movement_cutoff,"_Integ-",OFT_integration_period,".RDS"), compress = TRUE)

Freezing_TrackingAll <- RunPipeline(files,input_folder,FUN = Freezing_pipeline)
saveRDS(Freezing_TrackingAll, file = paste0(output_folder, DateTime,"_",ExpGroup,"_Freezing_TrackingAll_CutOff-",Freezing_cutoff,"_Integ-",Freezing_integration_period,".RDS"), compress = TRUE)

#get a report across all files analysed in batch
Report_OFT_TrackingAll <- MultiFileReport(OFT_TrackingAll)
saveRDS(Report_OFT_TrackingAll, file = paste0(output_folder, DateTime,"_",ExpGroup,"_Report_OFT_TrackingAll.RDS"), compress = TRUE)

Report_Freezing_TrackingAll <- MultiFileReport(Freezing_TrackingAll) #note: This does not include binned information
saveRDS(Report_Freezing_TrackingAll, file = paste0(output_folder, DateTime,"_",ExpGroup,"_Report_Freezing_TrackingAll_CutOff-",Freezing_cutoff,"_Integ-",Freezing_integration_period,".RDS"), compress = TRUE)
#Report[,1:6] #here only first 6 columns are displayed

#create PDF reports for each file. They will be saved in the cwd
PlotDensityPathsMyColors.Multi.PDF(OFT_TrackingAll, filename = paste0(output_folder, DateTime,"_",ExpGroup,"_PlotDensityMaps_TrackingAll"),points = c("Occiput"), add_zones = TRUE)
PlotZoneVisitsWithMoving.Multi.PDF(OFT_TrackingAll, filename = paste0(output_folder, DateTime,"_",ExpGroup,"_PlotZoneVisits_TrackingAll"), points = c("Occiput"), width = 10, height = 4)
PlotLabels.Multi.PDF(Freezing_TrackingAll, filename = paste0(output_folder, DateTime,"_",ExpGroup,"_PlotLabels_TrackingAll"), width = 10, height = 2)



if(ExpGroup == "Opto"){
  OFTBinAnalysisAll <- MultiFileBinanalysis(OFT_TrackingAll, FUN = OFTAnalysis ,movement_cutoff = OFT_movement_cutoff, integration_period = OFT_integration_period, points = c("Occiput"))
  saveRDS(OFTBinAnalysisAll, file = paste0(output_folder, DateTime,"_",ExpGroup,"_OFTBinAnalysisAll_CutOff-",OFT_movement_cutoff,"_Integ-",OFT_integration_period,".RDS"), compress = TRUE)
  
  FreezingBinAnalysisAll <- MultiFileBinanalysis(Freezing_TrackingAll, FUN = FreezingAnalysis ,cutoff_freezing = Freezing_cutoff,integration_period = Freezing_integration_period, Object = "Mouse", points = "Occiput")
  saveRDS(FreezingBinAnalysisAll, file = paste0(output_folder, DateTime,"_",ExpGroup,"_FreezingBinAnalysisAll_CutOff-",Freezing_cutoff,"_Integ-",Freezing_integration_period,".RDS"), compress = TRUE)

  #plotting binned data by Opsin and FiberLocation
  dataPlot <- OFTBinAnalysisAll %>% separate(file, into = c("VideoName_AnimalID", NA), sep = "(?<=.)(?=_MAH)|(?<=.)(?=-LED1)")
  for (ww in 1:nrow(dataPlot)){
    if(startsWith(dataPlot[ww,1],"Maus_35")){
      dataPlot[ww,1] <- gsub("_LED1","",dataPlot[ww,1])
    } 
  }
  dataPlot <- dataPlot %>% left_join(OptoGroupFile, by = "VideoName_AnimalID") %>% filter(FiberLocGroup != "unknown")
  dataPlot$FiberLocGroup <- factor(dataPlot$FiberLocGroup, levels = c("OnTarget_antPAG","OffTarget"))
  dataPlot$`type of opsin` <- factor (dataPlot$`type of opsin`, levels = c("ChR2","ArchT3.0"))
  dataPlot$bin[dataPlot$bin == "B1"] <- "BL"
  dataPlot$bin[dataPlot$bin == "B2"] <- "Stim"
  
  PlotOFTparamOpto.Multi.PDF(dataPlot, OFTBinAnalysisAll, 
                             filename = paste0(output_folder, DateTime,"_OptoAll_OFTBinAnalysisAll_CutOff-",OFT_movement_cutoff,"_Integ-",OFT_integration_period), 
                             width = 10, height = 6)
  
  PlotDensityPaths_binned.Multi.PDF(OFT_TrackingAll, filename = paste0(output_folder, DateTime,"_OptoBinPlot_OFTBinAnalysisAll_CutOff-",OFT_movement_cutoff,"_Integ-",OFT_integration_period), points = "Occiput")
  
  
  #Statistical Analysis
  BLdata <- dataPlot %>% filter(bin == "BL") %>% arrange(VideoName_AnimalID)
  Stimdata <- dataPlot %>% filter(bin == "Stim")%>% arrange(VideoName_AnimalID)
  
  if(all(BLdata$VideoName_AnimalID == Stimdata$VideoName_AnimalID)){
    IncreaseData <- Stimdata[,3:37]/BLdata[,3:37]*100
    IncreaseData$VideoName_AnimalID <- BLdata$VideoName_AnimalID
    IncreaseData$TypeOfOpsin <- BLdata$`type of opsin` %>% as.character() 
    IncreaseData$FiberLocGroup <- BLdata$FiberLocGroup %>% as.character()
    IncreaseData <- IncreaseData %>% relocate(VideoName_AnimalID)
  }else{
    warning("Can not calculate ratio due to non-matching data")
  }
  
  for (FiberGroup in c("OnTarget_antPAG","OffTarget")){
    ChR2_BLvsStim_allAnimals <- Opto_BLvsStim_Ttest(BLdata[BLdata$FiberLocGroup == FiberGroup,], Stimdata[Stimdata$FiberLocGroup == FiberGroup,], 3:37, "ChR2" , Report = TRUE)
    write.xlsx(ChR2_BLvsStim_allAnimals$Report %>% extract(), 
               file = paste0(output_folder, DateTime,"_Opto_ChR2_",FiberGroup,"_OFTstats_CutOff-",OFT_movement_cutoff,"_Integ-",OFT_integration_period,".xlsx"))
    
    ArchT_BLvsStim_allAnimals <- Opto_BLvsStim_Ttest(BLdata[BLdata$FiberLocGroup == FiberGroup,], Stimdata[Stimdata$FiberLocGroup == FiberGroup,], 3:37, "ArchT3.0", Report = TRUE)
    write.xlsx(ArchT_BLvsStim_allAnimals$Report %>% extract(), 
               file = paste0(output_folder, DateTime,"_Opto_ArchT_",FiberGroup,"_OFTstats_CutOff-",OFT_movement_cutoff,"_Integ-",OFT_integration_period,".xlsx"))
    
  }
}

if(ExpGroup == "DREADD"){
  #plotting DREADD data by DREADDreceptor and Condition
  dataPlot <- Report_OFT_TrackingAll %>% separate(file, into = c("VideoName_AnimalID_Cond", NA), sep = "(?<=.)(?=\\dDLC_)") %>% 
    separate(VideoName_AnimalID_Cond, into = c("VideoName_AnimalID", "Condition"), sep = "-") %>% 
    group_by(VideoName_AnimalID, Condition) %>% 
    summarize_all(mean, na.rm = TRUE) %>%
    left_join(DREADDGroupFile, by = "VideoName_AnimalID")
    
  dataPlot$DREADDreceptor <- factor(dataPlot$DREADDreceptor, levels = c("DREADD_neg","hM3Dq","hM4Di"))

  PlotOFTparamDREADD.Multi.PDF(dataPlot, Report_OFT_TrackingAll, 
                               filename = paste0(output_folder, DateTime,"_DREADD_OFTAnalysisAll_Occiput-Time-Stationary_CutOff-",OFT_movement_cutoff,"_Integ-",OFT_integration_period,".tiff"), 
                               width = 8, height = 4)
  
  #Statistical Analysis
  BLdata <- dataPlot %>% filter(Condition == "BL") %>% arrange(VideoName_AnimalID)
  Stimdata <- dataPlot %>% filter(Condition == "Clo")%>% arrange(VideoName_AnimalID)
  
  if(all(BLdata$VideoName_AnimalID == Stimdata$VideoName_AnimalID)){
    IncreaseData <- Stimdata[,3:37]/BLdata[,3:37]*100
    IncreaseData$VideoName_AnimalID <- BLdata$VideoName_AnimalID
    IncreaseData$DREADDreceptor <- BLdata$DREADDreceptor %>% as.character() 
    IncreaseData <- IncreaseData %>% relocate(VideoName_AnimalID)
  }else{
    warning("Can not calculate ratio due to non-matching data")
  }
  

  hM3Dq_BLvsStim_allAnimals <- DREADD_BLvsStim_Ttest(BLdata, Stimdata, 3:37, "hM3Dq" , Report = TRUE)
  write.xlsx(hM3Dq_BLvsStim_allAnimals$Report %>% extract(), 
             file = paste0(output_folder, DateTime,"_DREADD_hM3Dq_OFTstats_CutOff-",OFT_movement_cutoff,"_Integ-",OFT_integration_period,".xlsx"))
  
  hM4Di_BLvsStim_allAnimals <- DREADD_BLvsStim_Ttest(BLdata, Stimdata, 3:37, "hM4Di", Report = TRUE)
  write.xlsx(hM4Di_BLvsStim_allAnimals$Report %>% extract(), 
             file = paste0(output_folder, DateTime,"_DREADD_hM4Di_OFTstats_CutOff-",OFT_movement_cutoff,"_Integ-",OFT_integration_period,".xlsx"))
  
  DREADDneg_BLvsStim_allAnimals <- DREADD_BLvsStim_Ttest(BLdata, Stimdata, 3:37, "DREADD_neg", Report = TRUE)
  write.xlsx(DREADDneg_BLvsStim_allAnimals$Report %>% extract(), 
             file = paste0(output_folder, DateTime,"_DREADD_DREADDneg_OFTstats_CutOff-",OFT_movement_cutoff,"_Integ-",OFT_integration_period,".xlsx"))
  
  hM3DqVsDREADDneg <- DREADD_BLvsStim_Ttest(DF1 = IncreaseData, DF2 = NULL, cols = 2:36, c("hM3Dq","DREADD_neg"), Report = TRUE)
  write.xlsx(hM3DqVsDREADDneg$Report %>% extract(), 
             file = paste0(output_folder, DateTime,"_DREADD_hM3DqVSDREADDneg_OFTstats_CutOff-",OFT_movement_cutoff,"_Integ-",OFT_integration_period,".xlsx"))
  
  hM4DiVsDREADDneg <- DREADD_BLvsStim_Ttest(DF1 = IncreaseData, DF2 = NULL, cols = 2:36, c("hM4Di","DREADD_neg"), Report = TRUE)
  write.xlsx(hM4DiVsDREADDneg$Report %>% extract(), 
             file = paste0(output_folder, DateTime,"_DREADD_hM4DiVSDREADDneg_OFTstats_CutOff-",OFT_movement_cutoff,"_Integ-",OFT_integration_period,".xlsx"))
  
  hM3DqVshM4Di <- DREADD_BLvsStim_Ttest(DF1 = IncreaseData, DF2 = NULL, cols = 2:36, c("hM3Dq","hM4Di"), Report = TRUE)
  write.xlsx(hM3DqVshM4Di$Report %>% extract(), 
             file = paste0(output_folder, DateTime,"_DREADD_hM3DqVShM4Di_OFTstats_CutOff-",OFT_movement_cutoff,"_Integ-",OFT_integration_period,".xlsx"))
  
  
}

saveRDS(dataPlot, file = paste0(output_folder, DateTime,"_",ExpGroup,"_OFT_PlotData_CutOff-",OFT_movement_cutoff,"_Integ-",OFT_integration_period,".RDS"), compress = TRUE)

 
# #checkout the results from a specific file
# Tracking <- TrackingAll$`OF_top_DREADD_25-BL2DLC_resnet50_OpenFieldNov2shuffle1_200000_filtered.csv`
# PlotZoneVisits(Tracking, point = c("Nose","Occiput"))









#  ###########################################################################
#  #Perspective Transform
#  library(reticulate)
#  use_condaenv("r-reticulate")
#  cv2 <- import("cv2")
#  np <- import("numpy")
#  
#  path <- "H:/S3IT_Linux/DLCcopy/Maus_34-21_10_MAH00866DLC_resnet50_OpenFieldNov2shuffle1_200000.csv"
#  Tracking <- ReadDLCDataFromCSV(path, fps = 25)
#  
#  MedianTLxy <- c(Tracking[["data"]][["ArenaTL"]][["x"]] %>% median(),
#                  Tracking[["data"]][["ArenaTL"]][["y"]] %>% median())
#  MedianTRxy <- c(Tracking[["data"]][["ArenaTR"]][["x"]] %>% median(),
#                  Tracking[["data"]][["ArenaTR"]][["y"]] %>% median())  
#  MedianBRxy <- c(Tracking[["data"]][["ArenaBR"]][["x"]] %>% median(),
#                  Tracking[["data"]][["ArenaBR"]][["y"]] %>% median())  
#  MedianBLxy <- c(Tracking[["data"]][["ArenaBL"]][["x"]] %>% median(),
#                  Tracking[["data"]][["ArenaBL"]][["y"]] %>% median())
#  MedianCorners <- rbind(MedianTLxy,MedianTRxy,MedianBRxy,MedianBLxy)
#  
#  VideoCorners = np_array(MedianCorners, dtype = "float32") # Order: TL,TR,BR,BL; col1 = x, col2 = y
# 
# 
#  #Define reference coordinates of arena
#  RealTL <- c(0,0)
#  RealTR <- c(4000,0)
#  RealBR <- c(4000,4000)
#  RealBL <- c(0,4000)
#  RealCorners <- rbind(RealTL, RealTR,RealBR,RealBL)
#  Square = np_array(RealCorners, dtype = "float32") # real dimension of arena
#  
#  #calculate the transformation matrix needed
#  #use this mixtrix and multiply it with all tracked deeplabcut coordinates to correct for the perspetive distortion
#  TransMatrix = cv2$getPerspectiveTransform(VideoCorners, Square) 
#  #FinalDim = c(400,400)
# 
#  #InputShape <- py_get_attr(VideoCorners,"shape")
#  #result = cv2$warpPerspective(VideoCorners, TransMatrix, FinalDim)
#  #OutMatrix <- np_array(np$zeros(InputShape), dtype = "float32")
#  #result <- cv2$perspectiveTransform(VideoCorners,TransMatrix)
#  
#  TL <- np$append(VideoCorners[0],0)
#  TR <- np$append(VideoCorners[1],0)
#  BR <- np$append(VideoCorners[2],0)
#  BL <- np$append(VideoCorners[3],0)
#  
#  result <- np$dot(TransMatrix, TL)
#  result2 <- np$dot(TransMatrix, TR)
#  result3 <- np$dot(TransMatrix, BR)
#  results4 <- np$dot(TransMatrix, BL)
#  
#  BindRes<- rbind(result, result2, result3, results4)
#  z <- as.data.frame(BindRes)
#  colnames(z) <- c("x","y","z")
#  z <- z %>% mutate(Id = 1:n())
#  
#  pl <- as.data.frame(MedianCorners)
#  colnames(pl) <- c("x","y")
#  pl <- pl %>% mutate(Id = 1:n())
# 
#  z %>% ggplot(aes(x=x, y=y))+
#    geom_point(aes(col = Id))
#  
#  
#  
# ############################################################################# 
