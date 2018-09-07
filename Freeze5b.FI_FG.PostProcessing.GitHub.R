

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

## add clusters
load(paste(f.dir,"/",clusterfile,sep=""))
i <- 7

table(fulldata$population,clusters.list.sqrt[[i]]$clust.list.to.return[fulldata[,id.col]],useNA = "always")

##Assign ancestry labels
anc.labels.sqrt <- c("c.AF","c.EU","c.AF","c.AS","c.HS","c.SAS","c.EU")

# add cluster assignment to dataset
fulldata.sqrt <- cbind(fulldata, cluster.ancestry.sqrt=anc.labels.sqrt[clusters.list.sqrt[[i]]$clust.list.to.return[fulldata[,id.col]]])

table(fulldata.sqrt$ancestry,fulldata.sqrt$cluster.ancestry.sqrt,useNA = "always")

print(dim(fulldata))

write.table(fulldata, paste(f.dir,"/",out.pref,".for_analysis.csv",sep=""), row.names=F, col.names=T, quote=F, sep=',')
