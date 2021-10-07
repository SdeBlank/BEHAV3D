set.seed(1234)
getwd()
library(plyr)
library(readr)
library(dplyr)
library(reshape2)
library(zoo)
library(spatstat)
library(sp)
library(stats)
library(yaml)

pars = yaml.load_file("/Users/samdeblank/surfdrive/Shared/T cell paper/Stats reports/t cells/2021-08-11_ror1_CART_n3/BEHAV3D_config.yml")

## set directories where the files are located
working_directory <- pars$data_dir
setwd(working_directory)
output_dir=paste0(pars$output_dir,"/")
model_path <- pars$randomforest

###############################
######### Data import #########
###############################

# Import file-specific metadata for all images used in this analysis.
pat = pars$metadata_csv
metadata=read.csv(pars$metadata_csv, sep="\t")

# Function to import organoid data specifically from Imaris generated csv files
read_ims_csv <- function(pattern, recursive=TRUE) {
  read_plus <- function(flnm) {
    read_csv(flnm, skip = 3, col_types = cols()) %>% 
      mutate(filename = flnm)
  }
  files <- list.files(path = working_directory, pattern = pattern,  recursive = recursive)
  ims_csv <- ldply(files, read_plus)
  return(ims_csv)
}

# import Displacement^2
pat = "*Displacement"
displacement <- read_ims_csv(pattern=pat)

# import Speed
pat = "*Speed"
speed <- read_ims_csv(pattern=pat)

# import mean dead dye intensity values
pat = "*Intensity_Mean_Ch=3_Img=1"
red_lym <- read_ims_csv(pattern=pat)

# import Minimal distance to organoids
pat = "*dist_org"
dist_org <- read_ims_csv(pattern=pat)

# import Position
pat = "*Position"
pos <- read_ims_csv(pattern=pat)

### Join all Imaris information
master <- cbind(
  displacement[,c("Displacement^2","Time","TrackID" ,"ID")], 
  speed[,c("Speed" )], 
  dist_org[,c("Intensity Min")], 
  red_lym[,c("Intensity Mean")], 
  pos[,c("Position X" ,"Position Y" ,"Position Z","filename")]
)

### Get the basename from the filename for combination with metadata
master$basename <- gsub("_Position.csv", "", master$filename, perl=TRUE)
master$basename=basename(master$basename)
colnames(master) <- c("displacement","Time","TrackID","ID","speed","dist_org","red_lym","X-pos","Y-pos","Z-pos", "filename", "basename")

### Join the information of metadata to master:
master<-left_join(master, metadata)

### Create a unique TRACKID. 
### Each file processes with Imaris has separate TRACKIDs and these must be made unique before merging
category <- as.factor(master$filename)
ranks <- rank(-table(category), ties.method="first")
ranks <- as.data.frame(ranks)
ranks$filename <- row.names(ranks)
master <- left_join(master, ranks) 
master$TrackID2 <- with(master, interaction(ranks, TrackID, sep="_"))

### Remove filename
master$filename<-NULL

### save RDS for later use (e.g. Backprojection of classified TrackIDs)
saveRDS(master, paste0(output_dir,"raw_tcell_track_data.rds"))

###############################
####### Data processing #######
###############################

detach("package:reshape2", unload=TRUE)
detach("package:plyr", unload=TRUE)

###!!! VARIABLE
master <- master[which(master$Time<=300), ] ##Make sure that all the time-series have the same length, in this case 10hours

### Perform check for duplciates, should be empty
data_dup <- master%>%group_by(Time)%>%
  count(TrackID2) %>% 
  filter(n > 1) %>% 
  select(-n)

if (dim(data_dup)[1]!=0){
  stop("There are duplicates in the data, which should not be the case. Stopping execution...")
}

### For each well separately:
### Estimate if a T cells interacts with another T cell interact. 
### This is done by calculating the minmal distance between to the nearest neighbor. 

List = list() ## create list for each timepoint
List2 = list() ## create list for each experiment

### For loop to look for the distance to the nearest neighbor at each timepoint and each experiment
for ( m in unique(master$well)){
  distance_1<-master[which(master$well==m),]
  List = list()
  for (i in unique(distance_1$Time)){
    distanceDFi <- distance_1[distance_1$Time==i,]
    distanceDF2<-as.data.frame(distanceDFi[,c("Time","X-pos","Y-pos","Z-pos")])
    coordin<-ppx(distanceDF2, coord.type=c("t","s","s", "s")) ## create a multidimensional space-time pattern
    dist<- nndist(coordin)
    distanceDFi_dist<-cbind(distanceDFi, dist)
    List[[length(List)+1]] <-distanceDFi_dist   ## store to list object
  }
  master_distance <- data.frame(do.call(rbind, List)) ## convert List to dataframe
  List2[[length(List2)+1]] <-master_distance 
}


master_dist<-do.call(rbind, List2)
colnames(master_dist)[which(names(master_dist) == "dist")] <- "nearest_Tcell"

### Create a binary variable for Tcell contact based on distance
# VARIABLE
master_dist$contact_lym<- ifelse(master_dist$nearest_Tcell<10,1,0)

### Remove the variable TrackID and only use unique TrackID2 (unique identifier instead)
master_dist$TrackID<-master_dist$TrackID2
master_dist$TrackID2<-NULL

library(reshape2)
library(zoo)

### Since not all the tracks are tracked at all timepoints,
### Interpolate missing values and fill the NA (for each variable)

### Select the variables for which we need to interpolate NAs (numeric)
column_names<-names(master_dist)
column_names <- c("displacement", "dist_org", "red_lym", "contact_lym")

### Create a first dataset with refilled values for speed:
time_series<-acast(master_dist, Time ~ TrackID, value.var='speed',fun.aggregate = mean)

### rownames timepoints:
row.names(time_series)<-unique(master_dist$Time)

### Get rid of NA by interpolation
time_series_zoo<-zoo(time_series, row.names(time_series))
time_series_zoo<-na.approx(time_series_zoo) ## replace by interpolated value
time_series<-as.matrix(time_series_zoo)
time_series2<-melt(time_series)
data<-time_series2[complete.cases(time_series2), ] 
colnames(data)<-c("Time", "TrackID", "speed")

### Store this data for calculating lagged speed later:
time_series2_speed<-data

for (i in column_names){
  time_series<-acast(master_dist, Time ~ TrackID, value.var=i,fun.aggregate = mean)
  row.names(time_series)<-unique(master_dist$Time)
  ### get rid of NA
  time_series_zoo<-zoo(time_series,row.names(time_series))
  time_series_zoo<-na.approx(time_series_zoo) ## replace by last value
  time_series<-as.matrix(time_series_zoo)
  time_series2<-melt(time_series)
  new<-time_series2[complete.cases(time_series2), ] 
  data[ , ncol(data) + 1] <- new[3]                  # Append new column
  colnames(data)[ncol(data)] <- paste0(i)
}

library(dplyr)
### For cell interaction we need to consider the following:
### When two cells interact it is often the one cell moves and interacts with another one that is static
### In this case one might consider that only the motile cell is actively interacting and the static cell is just passively interacting
### To determine when a cell is actively interacting we measure for each cell what was its mean speed over the last 10 timepoints (20 mins)
time_series2_meanspeed <-time_series2_speed %>% 
  group_by(TrackID) %>%mutate(meanspeed=rollapply(speed,10,mean,align='right',fill=NA))

### Refill all missing values with the last value
time_series<-acast(time_series2_meanspeed, Time ~ TrackID, value.var='meanspeed',fun.aggregate = mean)
time_series_zoo<-zoo(time_series)
time_series_zoo<-na.locf(time_series_zoo, fromLast=T) ## replace by last value
time_series<-as.matrix(time_series_zoo)
time_series2_meanspeed<-melt(time_series)
colnames(time_series2_meanspeed)<-c("Time", "TrackID", "meanspeed")

### Remove last NAs
time_series2_meanspeed<-na.omit(time_series2_meanspeed)

### Create a dataframe with all the variables with corrected missing values
master_corrected <- data

### Join the information on the cell line, experiment number and well number
master_temp<- master[c("TrackID2", colnames(metadata)[-1])]
master_temp<-master_temp[!duplicated(master_temp$TrackID2),]
master_corrected<- left_join(master_corrected, master_temp, by=c("TrackID"="TrackID2"))

### Merge the information for the mean speed over the last 20 mins
master_corrected1<- merge(master_corrected, time_series2_meanspeed, by = c("Time","TrackID"))

### Update the binary variable for contact with organoids
### It can vary between experiments depending on the intensity of the T cells or organoids. 
### Check the threshold of contact in the imaging data and update in the metadata csv
master_corrected1$contact <- ifelse(master_corrected1$dist_org>master_corrected1$contact_threshold, 0,1)

### Plot the number of touching vs. non-touching T cells
library(ggplot2)
ggplot(master_corrected1, aes(x=contact, color=as.factor(exp_nr))) +
  geom_histogram(fill="white", alpha=0.5, position="identity")+facet_grid(organoid_line~well, scales = "free")

### Remove contact threshold variable
master_corrected1$contact_threshold<-NULL

### For clustering it is necessary to compare T cell tracks that have a similar length. 
### For that we select cell track that have at least 100 timepoints. 
### Detach package  'plyr' as it can interfere with 'dplyr'
detach("package:reshape2", unload=TRUE)
detach("package:plyr", unload=TRUE)

master_corrected2<-master_corrected1 %>% 
  group_by(TrackID) %>%arrange(TrackID)%>% filter(Time>00&Time<300)%>% filter(n() > 99)
### Create a variable for the relative Time
master_corrected2<-master_corrected2 %>% 
  group_by(TrackID) %>%arrange(Time)%>%mutate(Time2 = Time - first(Time))
### For the Tracks that have more then 100 timepoints filter only the first 100.
master_corrected2<-master_corrected2 %>% 
  group_by(TrackID) %>%arrange(TrackID)%>% filter(Time2<100)
### To exclude noise due to dead cells remove the dead t cells from the beginning
master_corrected3deadT0 <-master_corrected2%>%group_by(TrackID)%>%filter((Time2==0) & red_lym<10 )
master_corrected3 <-master_corrected2%>%filter(TrackID %in% master_corrected3deadT0$TrackID )
### Create a binary variable for live or dead cells:
master_corrected3$death<- ifelse(master_corrected3$red_lym<10,0,1)
### Create a variable for cumulative interaction with organoids
master_corrected3<-master_corrected3 %>% 
  group_by(TrackID) %>%mutate(contact2=(ave(contact, cumsum(!contact), FUN = cumsum)))
### Create a variable for T cells interact with other T cells while in the environment
master_corrected3$contact_lym<- ifelse(master_corrected3$contact==1,0,master_corrected3$contact_lym)
### For T cells inteacting in the environment keep as "interacting" only cells that had a mean speed in the last 20 minutes that is in the upper quantile.
master_corrected3<-master_corrected3%>%group_by(exp_nr)%>%mutate(contact_lym=ifelse(meanspeed<quantile(meanspeed,p=0.75),0,contact_lym))

### Save processed data on the tcells for possible further analysis
saveRDS(master_corrected3, file = paste0(output_dir,"processed_tcell_track_data.rds"))

###############################
##### Behavior prediction #####
###############################

library(scales)
library(randomForest)

# Load in the previously trained Random Forest model
tryCatch(
  {
  load(model_path)
  },
  error=function(e) {
    message("Provided model path does not exist") 
    message("Please set the path to the provided Random Forest model or \ntrain it on your own reference set (see: 'train_randomforest')")
    message(paste("Provided model path:", model_path))
    message("\nHere's the original error message:")
    message(e)
    stop()
    # Choose a return value in case of error
  }
)
#### Import new dataset to predict behavior (e.g. import dataset called "master_corrected3_example")
master_test<-readRDS(paste0(output_dir, "processed_tcell_track_data.rds"))
### Normalize the data
master_test2<-master_test%>% ungroup()%>%
  group_by(exp_nr) %>% 
  mutate(z.disp = (displacement-mean(displacement))/sd(displacement),z.speed = (speed-mean(speed))/sd(speed), z.red = (red_lym-mean(red_lym))/sd(red_lym))%>%
  mutate(q.disp=ifelse(z.disp>(quantile(z.disp, p=0.75)),z.disp,min(z.disp)), q.speed=ifelse(z.speed>(quantile(z.speed, p=0.75)),z.speed,min(z.speed)),q.red=ifelse(z.red>(quantile(z.red, p=0.75)),z.red,min(z.red)))%>%
  mutate(q.disp=rescale(q.disp, to=c(0,100)),q.speed=rescale(q.speed, to=c(0,100)),q.red=rescale(q.red, to=c(0,100)),s.contact=rescale(contact, to=c(0,1)),s.contact_lym=rescale(contact_lym, to=c(0,1)))%>%
  mutate(q.disp=q.disp/mean(quantile(q.disp, p=0.9999)),q.speed=q.speed/mean(quantile(q.speed, p=0.9999)),q.red=q.red/mean(quantile(q.red, p=0.9999)))%>%ungroup()

### Calculate time-series descriptive statistics
test_dataset <- master_test2%>%group_by(TrackID)%>% arrange(Time)%>%
  summarise(mean_speed= mean(q.speed),median_speed= median(q.speed),speed_sd= sd(q.speed),q3_speed= quantile(q.speed,0.90),
            mean_displacement = mean(q.disp),median_displacement = median(q.disp),
            displacement_sd=sd(q.disp),q3_disp= quantile(q.disp,0.90),
            mean_red_lym = mean(q.red),red_lym_sd=sd(q.red),q3_red= quantile(q.red,0.90),
            contact=mean(s.contact),mean_contact2=mean(contact2),contact2=max(contact2))

### create a  matrix from the predictors
test1<-as.matrix(test_dataset[,2:15]) 
model_predict_test<-predict(model,test1,type="response")
test_dataset_predicted<-cbind(model_predict_test,test_dataset)
test_dataset_predicted$cluster<-test_dataset_predicted$model_predict

### Plot the proportion of cells for each behavioral signature
## Join information of cell ID
classification<-test_dataset_predicted[,c(2,17)]
master_classified<-left_join(master_test,classification, by="TrackID")
cell_ID<-master_classified[!duplicated(master_classified$TrackID),c("TrackID","organoid_line","tcell_line","exp_nr", "well")] 


classified_tracks<-test_dataset_predicted
classified_tracks$cluster2<-classified_tracks$cluster
classified_tracks<-left_join(classified_tracks,cell_ID)
classified_tracks<-classified_tracks%>%arrange(cluster2)

### Quantify the number of cells per well
Number_cell_exp<-classified_tracks%>%group_by(well)%>%
  summarise(total_cell = n())
Percentage_clus<-left_join(Number_cell_exp,classified_tracks)
Percentage_clus <- Percentage_clus%>%group_by(cluster2,tcell_line, well,exp_nr, organoid_line)%>%
  summarise(total_cell = mean(total_cell), num_cluster=n())%>%mutate(percentage=num_cluster*100/total_cell)%>%ungroup()

### Plot proportion per well and per cell type
Per<-ggplot(Percentage_clus, aes(fill=as.factor(cluster2), y=percentage, x="")) + 
  geom_bar( stat="identity", position="fill")+ coord_flip()+ scale_y_reverse()
Per <- Per + facet_grid(interaction(exp_nr,well,organoid_line)  ~ interaction(organoid_line,tcell_line))
Per<-Per+theme_void()+ scale_fill_manual(values=c("gold3",
                                                  "darkolivegreen3",
                                                  "seagreen3",
                                                  "forestgreen",
                                                  "dodgerblue",
                                                  "cyan1",
                                                  "indianred",
                                                  "firebrick",
                                                  "brown1"))+theme(aspect.ratio = 1,strip.text.x = element_text(angle = 90))
Per

ggsave(paste0(output_dir,"RF_ClassProp_WellvsCelltype.png"), device="png")


### Run a chisq.test to see if the distibution is different between conditions
library(reshape2)
table1<-acast(Percentage_clus, cluster2~tcell_line, value.var = "percentage")
table1[is.na(table1)]<-0  ## If any cluster are missing convert NA to 0

## Run Pearson Chi sq test to see if there is a different distribution of clusters between cell types
chisq.test(table1)


