

args = commandArgs(trailingOnly=TRUE)
f.dir <- args[1]
out.pref <- args[2]
id.col <- args[3]
ped.file <- args[4]
trait <- args[5]
clusterfile <- args[6]

fulldata <- read.table(paste(f.dir,"/",ped.file,sep=""),header=T,sep=",",as.is=T) #44783
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

### define 'populations' and subset to groups that are well represented

fulldata$population[fulldata$ancestry == c("EU")] = 'EU'
fulldata$population[fulldata$ancestry == c("AF")] = 'AF'
fulldata$population[fulldata$ancestry == c("AS")] = 'AS'
fulldata$population[fulldata$ancestry == c("HS")] = 'HS'
fulldata$population[fulldata$ancestry == c("AMR")] = 'AMR'
fulldata$population[fulldata$ancestry == c("MIXED")] = 'MIXED'
fulldata$population[fulldata$ancestry == c("OTHER")] = 'OTHER'
#fulldata$population[fulldata$study == c("Amish")] = 'AMISH'
fulldata$population[fulldata$ancestry == c("SAS")] = 'SAS'
table(fulldata$population, useNA = 'always')
table(fulldata$population,fulldata$ancestry, useNA = 'always')

#subset populations
#subset to 5 major ancestry groups represented
fulldata_sub = subset(fulldata, subset = population %in% c("AF","EU",
                                                           "HS", "AS", "SAS")) #N=44769

table(fulldata_sub$population, useNA = 'always')

load(paste(f.dir,"/",clusterfile,sep=""))
i <- 7

table(fulldata_sub$population,clusters.list.sqrt[[i]]$clust.list.to.return[fulldata_sub$sample.id],useNA = "always")

##Assign ancestry labels
anc.labels.sqrt <- c("c.AF","c.EU","c.AF","c.AS","c.HS","c.SAS","c.EU")

# add cluster assignment to dataset
fulldata_sub_ancestry.sqrt <- cbind(fulldata_sub,
                                    cluster.ancestry.sqrt=anc.labels.sqrt[clusters.list.sqrt[[i]]$clust.list.to.return[fulldata_sub$sample.id]])

table(fulldata_sub_ancestry.sqrt$population,fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt,useNA = "always")

#gender vs ancestry & project or study
table(fulldata_sub_ancestry.sqrt$study,fulldata_sub_ancestry.sqrt$t2d,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt)
table(fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$t2d,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt)



###############################
#exclude sparse cells - T2D all
fulldata_sub_ancestry.sqrt$T2D.PCancestry <- fulldata_sub_ancestry.sqrt$t2d_ctrl

##
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.AF" & fulldata_sub_ancestry.sqrt$topmed_project=="SAFS")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.AF" & fulldata_sub_ancestry.sqrt$topmed_project=="SAS")] <- NA

fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.AS" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.AS" & fulldata_sub_ancestry.sqrt$topmed_project=="HyperGEN_GENOA")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.AS" & fulldata_sub_ancestry.sqrt$topmed_project=="JHS")] <- NA

fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.EU" & fulldata_sub_ancestry.sqrt$topmed_project=="HyperGEN_GENOA")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.EU" & fulldata_sub_ancestry.sqrt$topmed_project=="SAFS")] <- NA

fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="CFS")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="COPD")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="FHS")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="GeneSTAR")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="GOLDN")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="JHS")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="SAS")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="VTE")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.SAS" & fulldata_sub_ancestry.sqrt$topmed_project=="WHI")] <- NA

table(fulldata_sub_ancestry.sqrt$T2D.PCancestry,fulldata_sub_ancestry.sqrt$t2d_ctrl,useNA = 'always')
table(fulldata_sub_ancestry.sqrt$T2D.PCancestry,fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt,useNA = 'always')
table(fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$T2D.PCancestry,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt,useNA = 'always')

#exclude sparse cells - no preT2D
table(fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$t2d_nopre.ctrl,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt)

fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry <- fulldata_sub_ancestry.sqrt$t2d_nopre.ctrl

##
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.AF" & fulldata_sub_ancestry.sqrt$topmed_project=="SAFS")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.AF" & fulldata_sub_ancestry.sqrt$topmed_project=="SAS")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.AS" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.AS" & fulldata_sub_ancestry.sqrt$topmed_project=="HyperGEN_GENOA")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.AS" & fulldata_sub_ancestry.sqrt$topmed_project=="JHS")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.EU" & fulldata_sub_ancestry.sqrt$topmed_project=="HyperGEN_GENOA")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.EU" & fulldata_sub_ancestry.sqrt$topmed_project=="SAFS")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="CFS")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="COPD")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="FHS")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="GeneSTAR")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="GOLDN")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="JHS")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="SAS")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="VTE")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.SAS" & fulldata_sub_ancestry.sqrt$topmed_project=="WHI")] <- NA

table(fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt,useNA = 'always')

########################################
#########################################################
## recode high pret2d to t2d based on fg & hba1c
summary(as.numeric(fulldata_sub_ancestry.sqrt$last_exam_fg)) #n=16833
miss <- fulldata_sub_ancestry.sqrt[is.na(fulldata_sub_ancestry.sqrt$last_exam_fg),]
table(miss$study,useNA = 'always')
table(fulldata_sub_ancestry.sqrt$study,useNA = 'always')
summary(as.numeric(fulldata_sub_ancestry.sqrt$last_exam_hba1c)) #n=30989
miss <- fulldata_sub_ancestry.sqrt[is.na(fulldata_sub_ancestry.sqrt$last_exam_hba1c),]
table(miss$study,useNA = 'always')


fulldata_sub_ancestry.sqrt$t2dnew <- NA
fulldata_sub_ancestry.sqrt$t2dnew[fulldata_sub_ancestry.sqrt$t2d == 2 ] = 3
fulldata_sub_ancestry.sqrt$t2dnew[fulldata_sub_ancestry.sqrt$t2d == 0 ] = 4
fulldata_sub_ancestry.sqrt$t2dnew[fulldata_sub_ancestry.sqrt$t2d == 0 & is.na(fulldata_sub_ancestry.sqrt$last_exam_hba1c)] = 5
fulldata_sub_ancestry.sqrt$t2dnew[fulldata_sub_ancestry.sqrt$t2d == 2 & is.na(fulldata_sub_ancestry.sqrt$last_exam_hba1c)] = 6
fulldata_sub_ancestry.sqrt$t2dnew[fulldata_sub_ancestry.sqrt$t2d == 0 & is.na(fulldata_sub_ancestry.sqrt$last_exam_fg)] = 7
fulldata_sub_ancestry.sqrt$t2dnew[fulldata_sub_ancestry.sqrt$t2d == 2 & is.na(fulldata_sub_ancestry.sqrt$last_exam_fg)] = 8
fulldata_sub_ancestry.sqrt$t2dnew[fulldata_sub_ancestry.sqrt$last_exam_fg >= 6.105 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_fg)] = 9
fulldata_sub_ancestry.sqrt$t2dnew[fulldata_sub_ancestry.sqrt$last_exam_fg < 5.5 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_fg) & fulldata_sub_ancestry.sqrt$t2d != 2] = 10
fulldata_sub_ancestry.sqrt$t2dnew[fulldata_sub_ancestry.sqrt$last_exam_fg < 5.5 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_fg) & fulldata_sub_ancestry.sqrt$t2d == 1] = NA
fulldata_sub_ancestry.sqrt$t2dnew[fulldata_sub_ancestry.sqrt$last_exam_hba1c >= 6.1 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_hba1c)] = 11
fulldata_sub_ancestry.sqrt$t2dnew[fulldata_sub_ancestry.sqrt$last_exam_hba1c < 5.4 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_hba1c)& fulldata_sub_ancestry.sqrt$t2d != 2 ] = 12
fulldata_sub_ancestry.sqrt$t2dnew[fulldata_sub_ancestry.sqrt$last_exam_hba1c < 5.4 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_hba1c)& fulldata_sub_ancestry.sqrt$t2d == 1 ] = NA

fulldata_sub_ancestry.sqrt$t2dnew[fulldata_sub_ancestry.sqrt$t2dnew %in% c(4,5,7,10,12)] <- 0
fulldata_sub_ancestry.sqrt$t2dnew[fulldata_sub_ancestry.sqrt$t2dnew %in% c(3,6,8,9,11)] <- 1


table(fulldata_sub_ancestry.sqrt$t2dnew,fulldata_sub_ancestry.sqrt$t2d,useNA = 'always')
table(fulldata_sub_ancestry.sqrt$t2dnew,useNA = 'always')

##exclude sparse cells - T2D / no preT2D
table(fulldata_sub_ancestry.sqrt$t2dnew,fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt,useNA = 'always')

fulldata_sub_ancestry.sqrt$t2dnew_excl_pcancestry <- fulldata_sub_ancestry.sqrt$t2dnew
##
fulldata_sub_ancestry.sqrt$t2dnew_excl_pcancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.AF" & fulldata_sub_ancestry.sqrt$topmed_project=="SAFS")] <- NA
fulldata_sub_ancestry.sqrt$t2dnew_excl_pcancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.AF" & fulldata_sub_ancestry.sqrt$topmed_project=="SAS")] <- NA
fulldata_sub_ancestry.sqrt$t2dnew_excl_pcancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.AS" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen")] <- NA
fulldata_sub_ancestry.sqrt$t2dnew_excl_pcancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.AS" & fulldata_sub_ancestry.sqrt$topmed_project=="HyperGEN_GENOA")] <- NA
fulldata_sub_ancestry.sqrt$t2dnew_excl_pcancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.AS" & fulldata_sub_ancestry.sqrt$topmed_project=="JHS")] <- NA
fulldata_sub_ancestry.sqrt$t2dnew_excl_pcancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.EU" & fulldata_sub_ancestry.sqrt$topmed_project=="HyperGEN_GENOA")] <- NA
fulldata_sub_ancestry.sqrt$t2dnew_excl_pcancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.EU" & fulldata_sub_ancestry.sqrt$topmed_project=="SAFS")] <- NA
fulldata_sub_ancestry.sqrt$t2dnew_excl_pcancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen")] <- NA
fulldata_sub_ancestry.sqrt$t2dnew_excl_pcancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="CFS")] <- NA
fulldata_sub_ancestry.sqrt$t2dnew_excl_pcancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="COPD")] <- NA
fulldata_sub_ancestry.sqrt$t2dnew_excl_pcancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="FHS")] <- NA
fulldata_sub_ancestry.sqrt$t2dnew_excl_pcancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="GOLDN")] <- NA
fulldata_sub_ancestry.sqrt$t2dnew_excl_pcancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="GeneSTAR")] <- NA
fulldata_sub_ancestry.sqrt$t2dnew_excl_pcancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="JHS")] <- NA
fulldata_sub_ancestry.sqrt$t2dnew_excl_pcancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="SAS")] <- NA
fulldata_sub_ancestry.sqrt$t2dnew_excl_pcancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.HS" & fulldata_sub_ancestry.sqrt$topmed_project=="VTE")] <- NA
fulldata_sub_ancestry.sqrt$t2dnew_excl_pcancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.SAS" & fulldata_sub_ancestry.sqrt$topmed_project=="WHI")] <- NA
table(fulldata_sub_ancestry.sqrt$t2dnew_excl_pcancestry,fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt,useNA = 'always')



table(fulldata_sub_ancestry.sqrt$t2dnew,fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$population,useNA = 'always')

fulldata_sub_ancestry.sqrt$t2dnew_excl <- fulldata_sub_ancestry.sqrt$t2dnew

fulldata_sub_ancestry.sqrt$t2dnew_excl[which(fulldata_sub_ancestry.sqrt$population=="AF" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen")] <- NA
fulldata_sub_ancestry.sqrt$t2dnew_excl[which(fulldata_sub_ancestry.sqrt$population=="AS" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen")] <- NA
fulldata_sub_ancestry.sqrt$t2dnew_excl[which(fulldata_sub_ancestry.sqrt$population=="HS" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen")] <- NA
table(fulldata_sub_ancestry.sqrt$t2dnew_excl,fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$population,useNA = 'always')


### Write out files for analysis

## self-report ancestries
for(anc in c("AF","EU","AS","HS","SAS")) {
  write.table(fulldata_sub_ancestry.sqrt[which(fulldata_sub_ancestry.sqrt$ancestry==anc),], paste(f.dir,"/",out.pref,".",anc,".for_analysis.csv",sep=""), row.names=F, col.names=T, quote=F, sep=',')

}
## PC-clustered ancestries
for(anc in c("c.AF","c.EU","c.AS","c.HS","c.SAS")) {
  write.table(fulldata_sub_ancestry.sqrt[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt==anc),], paste(f.dir,"/",out.pref,".",anc,".for_analysis.csv",sep=""), row.names=F, col.names=T, quote=F, sep=',')

  
}

## Full analysis set
write.table(fulldata_sub_ancestry.sqrt, paste(f.dir,"/",out.pref,".for_analysis.csv",sep=""), row.names=F, col.names=T, quote=F, sep=',')

