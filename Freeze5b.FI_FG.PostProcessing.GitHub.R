

args = commandArgs(trailingOnly=TRUE)
f.dir <- args[1]
out.pref <- args[2]
id.col <- args[3]
ped.file <- args[4]
trait <- args[5]
clusterfile <- args[6]

fulldata <- read.table(paste(f.dir,"/",ped.file,sep=""),header=T,sep=",",as.is=T)
print(dim(fulldata))
head(fulldata)

## add clusters
load(paste(f.dir,"/",clusterfile,sep=""))
i <- 7

table(fulldata$ancestry,clusters.list.sqrt[[i]]$clust.list.to.return[fulldata[,id.col]],useNA = "always")

##Assign ancestry labels
anc.labels.sqrt <- c("c.AF","c.EU","c.AF","c.AS","c.HS","c.SAS","c.EU")

# add cluster assignment to dataset
fulldata.sqrt <- cbind(fulldata, cluster.ancestry.sqrt=anc.labels.sqrt[clusters.list.sqrt[[i]]$clust.list.to.return[fulldata[,id.col]]])

table(fulldata.sqrt$ancestry,fulldata.sqrt$cluster.ancestry.sqrt,useNA = "always")

print(dim(fulldata))

###
#Add AgeSq and log(FI) recode sex as M/F
fulldata.sqrt$age_FG_sq = fulldata.sqrt$age_FG**2
fulldata.sqrt$age_FI_sq = fulldata.sqrt$age_FI**2
fulldata.sqrt$logFI = log10(fulldata.sqrt$FastingInsulin)
for(i in 1:length(fulldata.sqrt$sex)){
  if(fulldata.sqrt$sex[i]==1){
    fulldata.sqrt$sex[i]='M'
  }
  if(fulldata.sqrt$sex[i]==2){
    fulldata.sqrt$sex[i]='F'
  }
  
}
#######################
## exclude sparse cell new t2d definition PC defined ancestry
table(fulldata.sqrt$STUDY_ANCESTRY, fulldata.sqrt$cluster.ancestry.sqrt)

fulldata.sqrt$FastingGlucose.PCancestry <- fulldata.sqrt$FastingGlucose

fulldata.sqrt$FastingGlucose.PCancestry[which(fulldata.sqrt$cluster.ancestry.sqrt=="c.AF" & fulldata.sqrt$STUDY_ANCESTRY %in% c("ARIC_EU","MESA_EU","SAFS_HS","WHI_EU","SAS"))] <- NA


write.table(fulldata.sqrt, paste(f.dir,"/",out.pref,".for_analysis.csv",sep=""), row.names=F, col.names=T, quote=F, sep=',')

### Write out files for analysis

## PC-clustered ancestries
for(anc in c("c.AF","c.EU","c.AS","c.HS","c.SAS")) {
  write.table(fulldata.sqrt[which(fulldata.sqrt$cluster.ancestry.sqrt==anc),], paste(f.dir,"/",out.pref,".",anc,".for_analysis.csv",sep=""), row.names=F, col.names=T, quote=F, sep=',')
  
  
}


