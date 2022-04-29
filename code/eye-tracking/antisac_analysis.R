
#  This function runs an analysis pipeline for an antisaccade task recorded in .asc format (Eyelink 1000)
# Output files are a .csv with a summary (mean + sd) of the analyszed subject and a .csv with all valid trials (correct and errors)

#  Copyright (C) April 2020, last modified April 2022
#   J.Waldthaler, A. Sperlich, D. Pedrosa
#   University Hospital of Gießen and Marburg
#
#   This software may be used, copied, or redistributed as long as it is
#   not sold and this copyright notice is reproduced on each copy made.
#   This routine is provided as is without any express or implied
#   warranties whatsoever.

# load required packages
require(eyelinker)
require(dplyr)
require(ggplot2)
library(scales)
library(tidyr)
library(stringr)
library(naniar)

# choose subject and condition
subject='XX_XXX_XXXX'
condition = 'off'  # should be "off", "130", "60" or "180"

# select eye, should be "R" whenever right eye was recorded 
eye="R"

#select the three files with the eyelink 1000 data
dir <- paste("/Documents/",subject, condition, sep="/")
file1 <- paste (subject,condition,"1.asc", sep="_")
file2 <- paste (subject,condition,"2.asc", sep="_")
file3 <- paste (subject,condition,"3.asc", sep="_")

dat <- read.asc(file1)
raw <- dat$raw
msg<- dat$msg
msg$text <- gsub("green","ps",msg$text)
msg$text <- gsub("red","as",msg$text)
msg<- filter(msg, msg$text!="startFix")
msg<- filter(msg, msg$text!="startStim")
msg<- filter(msg, msg$text!="endStim")
#msg<- filter(msg, msg$text!="stop recording")
msg<- filter(msg, msg$text!="start recording")
msg<- filter(msg, msg$text!="Trial_correct")
msg<- filter(msg, msg$text!="BAD_EYE")

msg$text2 <- msg$text
msg<-replace_with_na_at(msg,.vars = "text2",condition = ~.x != "stop recording")

for(i in 2:nrow(msg)){
  if(!is.na(msg[(i),4])){
    msg[i,2]= ((msg[i,2]) - 1030)
  }
}
msg$text[msg$text == "stop recording"] <- "endFix"
msg$text2 <- NULL
msg$stime <-msg$time 
raw <-merge (dat$raw,msg, by="time")
raw$stime <- raw$time
#load sac data
sac <- dat$sac
sac$amplitude <- sac$exp-sac$sxp
#remove microsaccades
sac <- filter(sac, sac$ampl >5)
sac.R  <- filter(sac, sac$eye==(eye))
sac.R <- merge(sac.R, msg, by="stime", all=TRUE)

sac.R$type=sac.R$text
sac.R$typea=sac.R$text
sac.R$dir=sac.R$text
sac.R$dir2=sac.R$text
sac.R<-replace_with_na_at(sac.R,.vars = "text",condition = ~.x != "endFix")
sac.R<-replace_with_na_at(sac.R,.vars = "type",condition = ~.x != ("ps"))
sac.R<-replace_with_na_at(sac.R,.vars = "typea",condition = ~.x != ("as"))
sac.R<-replace_with_na_at(sac.R,.vars = "dir",condition = ~.x != ("right"))
sac.R<-replace_with_na_at(sac.R,.vars = "dir2",condition = ~.x != ("left"))
sac.R$typ<- sac.R$type
sac.R$typ[!is.na(sac.R$typea)] = sac.R$typea[!is.na(sac.R$typea)] 
sac.R$dir[!is.na(sac.R$dir2)] = sac.R$dir2[!is.na(sac.R$dir2)] 
sac.R$type <- NULL
sac.R$typea <- NULL
sac.R$dir2 <- NULL

sac.R<-fill(sac.R, typ)
sac.R<-fill(sac.R, dir)

for(i in 2:nrow(sac.R)){
  if(!is.na(sac.R[(i-1),15])){
    sac.R[i,18]= (sac.R[i,1] - sac.R[(i-1),1])
  }
}

sac.R<-sac.R %>% 
  rename(rt=`V18`)

sac.R <- fill(sac.R, block.y)

sac.R$block <- NULL
sac.R$dur <- NULL
sac.R$syp <- NULL
sac.R$eyp <- NULL
sac.R$eye <- NULL
sac.R$time <-NULL

# remove trials with latencies longer than stimulus presentation time, then remove based on 
# individual mean and sd, then remove anticipated sac also
sac.R <- filter (sac.R, sac.R$rt<1000)
sac.R <- filter (sac.R, sac.R$rt<(mean(sac.R$rt)+2*sd(sac.R$rt)))
sac.R <- filter (sac.R, sac.R$rt>89)


for(i in 1:nrow(sac.R)){
  if(sac.R[i,11]=="left"){
    sac.R[i,14]= -500
  }
}
for(i in 1:nrow(sac.R)){
  if(sac.R[i,11]=="right"){
    sac.R[i,14]= 500
  }
}

#  calculate gain - CAVE since we did not use gain in this study,
# this is just an approximation and shouldn't be considered the true value 
sac.R$gain <- sac.R$amplitude/sac.R$V14

# filter for prosac errors and correct antisac. separately
pro.R<-filter(sac.R, sac.R$gain>0)
as.R<-filter(sac.R, sac.R$gain<0)
pro.R1<-pro.R
as.R1 <- as.R
sac.R1 <- sac.R
as.R1$trial_nr <- as.R1$block.x
pro.R1$trial_nr <- pro.R1$block.x

# same for other two files
dat <- read.asc(file2)
raw <- dat$raw
msg<- dat$msg

msg$text <- gsub("green","ps",msg$text)
msg$text <- gsub("red","as",msg$text)
msg<- filter(msg, msg$text!="startFix")
msg<- filter(msg, msg$text!="startStim")
msg<- filter(msg, msg$text!="endStim")
#msg<- filter(msg, msg$text!="stop recording")
msg<- filter(msg, msg$text!="start recording")
msg<- filter(msg, msg$text!="Trial_correct")
msg<- filter(msg, msg$text!="BAD_EYE")

msg$text2 <- msg$text
msg<-replace_with_na_at(msg,.vars = "text2",condition = ~.x != "stop recording")

for(i in 2:nrow(msg)){
  if(!is.na(msg[(i),4])){
    msg[i,2]= ((msg[i,2]) - 1030)
  }
}
msg$text[msg$text == "stop recording"] <- "endFix"
msg$text2 <- NULL
msg$stime <-msg$time 
raw <-merge (dat$raw,msg, by="time")
raw$stime <- raw$time

sac <- dat$sac

sac$amplitude <- sac$exp-sac$sxp
sac <- filter(sac, sac$ampl >5)
sac.R  <- filter(sac, sac$eye==(eye))
sac.R <- merge(sac.R, msg, by="stime", all=TRUE)

sac.R$type=sac.R$text
sac.R$typea=sac.R$text
sac.R$dir=sac.R$text
sac.R$dir2=sac.R$text
sac.R<-replace_with_na_at(sac.R,.vars = "text",condition = ~.x != "endFix")
sac.R<-replace_with_na_at(sac.R,.vars = "type",condition = ~.x != ("ps"))
sac.R<-replace_with_na_at(sac.R,.vars = "typea",condition = ~.x != ("as"))
sac.R<-replace_with_na_at(sac.R,.vars = "dir",condition = ~.x != ("right"))
sac.R<-replace_with_na_at(sac.R,.vars = "dir2",condition = ~.x != ("left"))
sac.R$typ<- sac.R$type
sac.R$typ[!is.na(sac.R$typea)] = sac.R$typea[!is.na(sac.R$typea)] 
sac.R$dir[!is.na(sac.R$dir2)] = sac.R$dir2[!is.na(sac.R$dir2)] 
sac.R$type <- NULL
sac.R$typea <- NULL
sac.R$dir2 <- NULL

sac.R<-fill(sac.R, typ)
sac.R<-fill(sac.R, dir)
sac.R
for(i in 2:nrow(sac.R)){
  if(!is.na(sac.R[(i-1),15])){
    sac.R[i,18]= (sac.R[i,1] - sac.R[(i-1),1])
  }
}
sac.R<-sac.R %>% 
  rename(rt=`V18`)
sac.R <- fill(sac.R, block.y)

sac.R$block <- NULL
sac.R$dur <- NULL
sac.R$syp <- NULL
sac.R$eyp <- NULL
sac.R$eye <- NULL
sac.R$time <-NULL

sac.R <- filter (sac.R, sac.R$rt<1000)
sac.R <- filter (sac.R, sac.R$rt<(mean(sac.R$rt)+2*sd(sac.R$rt)))
sac.R <- filter (sac.R, sac.R$rt>89)

for(i in 1:nrow(sac.R)){
  if(sac.R[i,11]=="left"){
    sac.R[i,14]= -500
  }
}
for(i in 1:nrow(sac.R)){
  if(sac.R[i,11]=="right"){
    sac.R[i,14]= 500
  }
}
sac.R$gain <- sac.R$amplitude/sac.R$V14
sac.R2 <- sac.R
pro.R<-filter(sac.R, sac.R$gain>0)
as.R<-filter(sac.R, sac.R$gain<0)
pro.R2<-pro.R
as.R2 <- as.R
as.R2$trial_nr <- as.R2$block.x +50
pro.R2$trial_nr <- pro.R2$block.x +50


dat <- read.asc(file3)
raw <- dat$raw
msg<- dat$msg
msg$text <- gsub("green","ps",msg$text)
msg$text <- gsub("red","as",msg$text)
msg<- filter(msg, msg$text!="startFix")
msg<- filter(msg, msg$text!="startStim")
msg<- filter(msg, msg$text!="start recording")
msg<- filter(msg, msg$text!="Trial_correct")
msg<- filter(msg, msg$text!="BAD_EYE")

msg$text2 <- msg$text
msg<-replace_with_na_at(msg,.vars = "text2",condition = ~.x != "stop recording")

for(i in 2:nrow(msg)){
  if(!is.na(msg[(i),4])){
    msg[i,2]= ((msg[i,2]) - 1030)
  }
}
msg$text[msg$text == "stop recording"] <- "endFix"
msg$text2 <- NULL
msg$stime <-msg$time 
raw <-merge (dat$raw,msg, by="time")
raw$stime <- raw$time
raw.long <- dplyr::select(raw,time,xp,yp, block.x, text) 
raw.long <- mutate(raw.long,ts=(time-min(time))/1e3)

sac <- dat$sac

sac$amplitude <- sac$exp-sac$sxp
sac <- filter(sac, sac$ampl >5)
sac.R  <- filter(sac, sac$eye==(eye))
sac.R <- merge(sac.R, msg, by="stime", all=TRUE)
sac.R$type=sac.R$text
sac.R$typea=sac.R$text
sac.R$dir=sac.R$text
sac.R$dir2=sac.R$text
sac.R<-replace_with_na_at(sac.R,.vars = "text",condition = ~.x != "endFix")
sac.R<-replace_with_na_at(sac.R,.vars = "type",condition = ~.x != ("ps"))
sac.R<-replace_with_na_at(sac.R,.vars = "typea",condition = ~.x != ("as"))
sac.R<-replace_with_na_at(sac.R,.vars = "dir",condition = ~.x != ("right"))
sac.R<-replace_with_na_at(sac.R,.vars = "dir2",condition = ~.x != ("left"))
sac.R$typ<- sac.R$type
sac.R$typ[!is.na(sac.R$typea)] = sac.R$typea[!is.na(sac.R$typea)] 
sac.R$dir[!is.na(sac.R$dir2)] = sac.R$dir2[!is.na(sac.R$dir2)] 
sac.R$type <- NULL
sac.R$typea <- NULL
sac.R$dir2 <- NULL

sac.R<-fill(sac.R, typ)
sac.R<-fill(sac.R, dir)

for(i in 2:nrow(sac.R)){
  if(!is.na(sac.R[(i-1),15])){
    sac.R[i,18]= (sac.R[i,1] - sac.R[(i-1),1])
  }
}
sac.R<-sac.R %>% 
  rename(rt=`V18`)
sac.R <- fill(sac.R, block.y)

sac.R$block <- NULL
sac.R$dur <- NULL
sac.R$syp <- NULL
sac.R$eyp <- NULL
sac.R$eye <- NULL
sac.R$time <-NULL
sac.R <- filter (sac.R, sac.R$rt<1000)
sac.R <- filter (sac.R, sac.R$rt<(mean(sac.R$rt)+2*sd(sac.R$rt)))
sac.R <- filter (sac.R, sac.R$rt>89)

for(i in 1:nrow(sac.R)){
  if(sac.R[i,11]=="left"){
    sac.R[i,14]= -500
  }
}
for(i in 1:nrow(sac.R)){
  if(sac.R[i,11]=="right"){
    sac.R[i,14]= 500
  }
}
sac.R$gain <- sac.R$amplitude/sac.R$V14
sac.R3 <- sac.R
pro.R<-filter(sac.R, sac.R$gain>0)
as.R<-filter(sac.R, sac.R$gain<0)

pro.R3<-pro.R
as.R3 <- as.R
as.R3$trial_nr <- as.R3$block.x +100
pro.R3$trial_nr <- pro.R3$block.x +100

# combine data from all files
as.R <- bind_rows(as.R1, as.R2, as.R3)
pro.R <- bind_rows(pro.R1, pro.R2, pro.R3)
sac.R <-bind_rows(sac.R1, sac.R2, sac.R3)

#  summarize results
as_lat=mean(as.R$rt)
as_lat_sd=sd(as.R$rt)
error_lat=mean(pro.R$rt)
error_lat_sd=sd(pro.R$rt)
er_rate = count(pro.R)/count(sac.R)
results <- data.frame(count(sac.R), as_lat, as_lat_sd, error_lat, error_lat_sd, er_rate)

# some more cleaning for the csv with all saccades
as.R[,1:8] <- NULL
as.R$text <- NULL
as.R$typ <- NULL
as.R$V14 <- NULL

pro.R[,1:8] <- NULL
pro.R$text <- NULL
pro.R$typ <- NULL
pro.R$V14 <- NULL

# 1 = correct, 0 = error
as.R$response <-'1'
as.R$cond <- cond
pro.R$response <- '0'
pro.R$cond <- cond

all <- rbind(as.R, pro.R)

# save summary and all trials 
summary <- paste(subject,condition,"results_with-ex.csv", sep ="_")
all_trials <- paste(subject,condition,"all_with-ex.csv", sep="_")
write.csv(results,file=summary)
write.csv(all,file=all_trials)
