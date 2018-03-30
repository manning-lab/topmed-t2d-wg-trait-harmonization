

args = commandArgs(trailingOnly=TRUE)
f.dir <- args[1]
out.pref <- args[2]
id.col <- args[3]
ped.file <- args[4]
trait <- args[5]

fulldata <- read.table(paste(f.dir,"/",ped.file,sep=""),header=T,sep=",")
print(dim(fulldata))
head(fulldata)

## Restrict to people with age < 25
print("Before removing people with age<25")
table(fulldata$t2d,AGE_test=fulldata$t2d_age >= 25 | is.na(fulldata$t2d_age),useNA='always')

# March 22

## APPLY EXCLUSIONS
summary(fulldata$t2d_age)
fulldata = fulldata[(fulldata$t2d_age >= 25 | is.na(fulldata$t2d_age)),]

summary(fulldata$last_exam_age)
fulldata = subset(fulldata, last_exam_age >= 25) 

print("After removing people with age<25")
table(fulldata$t2d,AGE_test=fulldata$t2d_age >= 25 | is.na(fulldata$t2d_age),useNA='always')

##  !!  NEED TO CREATE THESE VARIABLES !!  ##
## t2d_ctrl t2d_superctrl
# create another t2d status classification for pre-DM 

table(fulldata$t2d,useNA='always')
fulldata$t2d_ctrl[fulldata$t2d == 2] = 1
fulldata$t2d_ctrl[fulldata$t2d == 1] = 0
fulldata$t2d_ctrl[fulldata$t2d == 0] = 0
fulldata$t2d_ctrl[is.na(fulldata$t2d)] = NA
with(fulldata,table(t2d,t2d_ctrl,useNA='always'))

table(fulldata$t2d,useNA='always')
fulldata$t2d_nopre.ctrl[fulldata$t2d == 2] = 1
fulldata$t2d_nopre.ctrl[fulldata$t2d == 1] = NA
fulldata$t2d_nopre.ctrl[fulldata$t2d == 0] = 0
fulldata$t2d_nopre.ctrl[is.na(fulldata$t2d)] = NA
with(fulldata,table(t2d,t2d_nopre.ctrl,useNA='always'))

table(fulldata$t2d,useNA='always')
fulldata$t2dpre_ctrl[fulldata$t2d == 2] = 1
fulldata$t2dpre_ctrl[fulldata$t2d == 1] = 1
fulldata$t2dpre_ctrl[fulldata$t2d == 0] = 0
fulldata$t2dpre_ctrl[is.na(fulldata$t2d)] = NA
with(fulldata,table(t2d,t2dpre_ctrl,useNA='always'))
table(fulldata$t2dpre_ctrl,useNA='always')


write.table(fulldata, paste(f.dir,"/",out.pref,".removed.IDs.post-processed.csv",sep=""), row.names=F, col.names=T, quote=F, sep=',')
