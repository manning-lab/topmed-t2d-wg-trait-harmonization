

args = commandArgs(trailingOnly=TRUE)
f.dir <- args[1]
out.pref <- args[2]
id.col <- args[3]
ped.file <- args[4]
trait <- args[5]
clusterfile <- args[6]

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

## add clusters

clusters <- read.table(paste(f.dir,"/",clusterfile,sep=""),sep=",",as.is=T,header=T)

clusters$clustered.ancestry.notPCscaled <- NA

clusters$clustered.ancestry.notPCscaled[which(clusters$cluster == 1)] <- "cAF" 
clusters$clustered.ancestry.notPCscaled[which(clusters$cluster %in% c(2,3,7,8))] <- "cEU"
clusters$clustered.ancestry.notPCscaled[which(clusters$cluster == 4)] <- "cAS"
clusters$clustered.ancestry.notPCscaled[which(clusters$cluster == 5)] <- "cHS" 
clusters$clustered.ancestry.notPCscaled[which(clusters$cluster == 6)] <- "cSAS" 
clusters$clustered.ancestry.notPCscaled[which(clusters$cluster %in% c(9, 10))] <- "cAmish"
head(clusters)

fulldata <- merge(fulldata,clusters[,c("sample.id","cluster","clustered.ancestry.notPCscaled")],by="sample.id",all.x=TRUE)
table(fulldata$ancestry, fulldata$clustered.ancestry.notPCscaled,useNA = "always")

write.table(fulldata, paste(f.dir,"/",out.pref,".for_analysis.csv",sep=""), row.names=F, col.names=T, quote=F, sep=',')
