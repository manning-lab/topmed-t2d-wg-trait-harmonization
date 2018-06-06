

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

clusters <- read.table(paste(f.dir,"/",clusterfile,sep=""),sep=",",as.is=T,header=T)

clusters$clustered.ancestry.notPCscaled <- NA

clusters$clustered.ancestry.notPCscaled[which(clusters$cluster == 1)] <- "cAF" 
clusters$clustered.ancestry.notPCscaled[which(clusters$cluster %in% c(2,3,7,8))] <- "cEU"
clusters$clustered.ancestry.notPCscaled[which(clusters$cluster == 4)] <- "cAS"
clusters$clustered.ancestry.notPCscaled[which(clusters$cluster == 5)] <- "cHS" 
clusters$clustered.ancestry.notPCscaled[which(clusters$cluster == 6)] <- "cSAS" 
clusters$clustered.ancestry.notPCscaled[which(clusters$cluster %in% c(9, 10))] <- "cAmish"
head(clusters)

print(dim(fulldata))
fulldata <- merge(fulldata,clusters[,c("sample.id","cluster","clustered.ancestry.notPCscaled")],by.x=id.col,by.y="sample.id",all.x=TRUE)

print(dim(fulldata))
table(fulldata$ancestry, fulldata$clustered.ancestry.notPCscaled,useNA = "always")

write.table(fulldata, paste(f.dir,"/",out.pref,".for_analysis.csv",sep=""), row.names=F, col.names=T, quote=F, sep=',')
