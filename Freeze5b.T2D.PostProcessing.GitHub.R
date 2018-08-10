

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

#gender vs ancestry & project 
table(fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$t2d,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt)


###############################
#exclude sparse cells - T2D all - PC cluster ancestry
table(fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$t2d_ctrl,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt)

fulldata_sub_ancestry.sqrt$T2D.PCancestry <- fulldata_sub_ancestry.sqrt$t2d_ctrl

##
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtAF" & fulldata_sub_ancestry.sqrt$topmed_project=="SAFS")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtAF" & fulldata_sub_ancestry.sqrt$topmed_project=="SAS")] <- NA

fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtAS" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtAS" & fulldata_sub_ancestry.sqrt$topmed_project=="HyperGEN_GENOA")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtAS" & fulldata_sub_ancestry.sqrt$topmed_project=="JHS")] <- NA

fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtEU" & fulldata_sub_ancestry.sqrt$topmed_project=="HyperGEN_GENOA")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtEU" & fulldata_sub_ancestry.sqrt$topmed_project=="SAFS")] <- NA

fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project=="CFS")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project=="COPD")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project=="FHS")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project=="GeneSTAR")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project=="GOLDN")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project=="JHS")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project=="SAS")] <- NA
fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project=="VTE")] <- NA

fulldata_sub_ancestry.sqrt$T2D.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtSAS" & fulldata_sub_ancestry.sqrt$topmed_project=="WHI")] <- NA


table(fulldata_sub_ancestry.sqrt$T2D.PCancestry,fulldata_sub_ancestry.sqrt$t2d_ctrl,useNA = 'always')
table(fulldata_sub_ancestry.sqrt$T2D.PCancestry,fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt,useNA = 'always')
table(fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$T2D.PCancestry,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt,useNA = 'always')


###############################
#exclude sparse cells - T2D all - self-report ancestry
table(fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$t2d_ctrl,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$population)

fulldata_sub_ancestry.sqrt$t2d_ctrl_excl <- fulldata_sub_ancestry.sqrt$t2d_ctrl

##
fulldata_sub_ancestry.sqrt$t2d_ctrl_excl[which(fulldata_sub_ancestry.sqrt$population=="AF" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen")] <- NA
fulldata_sub_ancestry.sqrt$t2d_ctrl_excl[which(fulldata_sub_ancestry.sqrt$population=="AS" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen")] <- NA
fulldata_sub_ancestry.sqrt$t2d_ctrl_excl[which(fulldata_sub_ancestry.sqrt$population=="HS" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen")] <- NA

table(fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$t2d_ctrl_excl,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$population)

###############################
#exclude sparse cells - no preT2D - PC cluster ancestry 

table(fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$t2d_nopre.ctrl,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt)

fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry <- fulldata_sub_ancestry.sqrt$t2d_nopre.ctrl

##
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtAF" & fulldata_sub_ancestry.sqrt$topmed_project=="SAFS")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtAF" & fulldata_sub_ancestry.sqrt$topmed_project=="SAS")] <- NA

fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtAS" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtAS" & fulldata_sub_ancestry.sqrt$topmed_project=="HyperGEN_GENOA")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtAS" & fulldata_sub_ancestry.sqrt$topmed_project=="JHS")] <- NA

fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtEU" & fulldata_sub_ancestry.sqrt$topmed_project=="HyperGEN_GENOA")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtEU" & fulldata_sub_ancestry.sqrt$topmed_project=="SAFS")] <- NA

fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project=="CFS")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project=="COPD")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project=="FHS")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project=="GeneSTAR")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project=="GOLDN")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project=="JHS")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project=="SAS")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project=="VTE")] <- NA
fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtSAS" & fulldata_sub_ancestry.sqrt$topmed_project=="WHI")] <- NA

table(fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$T2Dnopre.PCancestry,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt,useNA = 'always')

###############################
#exclude sparse cells - no preT2D - self report ancestry
table(fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$t2d_nopre.ctrl,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$population)

fulldata_sub_ancestry.sqrt$t2d_nopre.ctrl_excl <- fulldata_sub_ancestry.sqrt$t2d_nopre.ctrl

##
fulldata_sub_ancestry.sqrt$t2d_nopre.ctrl_excl[which(fulldata_sub_ancestry.sqrt$population=="AF" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen")] <- NA

fulldata_sub_ancestry.sqrt$t2d_nopre.ctrl_excl[which(fulldata_sub_ancestry.sqrt$population=="AS" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen")] <- NA

fulldata_sub_ancestry.sqrt$t2d_nopre.ctrl_excl[which(fulldata_sub_ancestry.sqrt$population=="HS" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen")] <- NA

table(fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$t2d_nopre.ctrl_excl,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$population,useNA = 'always')


########################################
#########################################################
## recode high pret2d to t2d based on fg & hba1c
summary(as.numeric(fulldata_sub_ancestry.sqrt$last_exam_fg)) #NA=16889
miss <- fulldata_sub_ancestry.sqrt[is.na(fulldata_sub_ancestry.sqrt$last_exam_fg),]
table(miss$study,useNA = 'always') #which studies have missing data? can we pull in data from the FG/FI/a1c files???
## ARIC could be added from glycemic trait files, but when was it measured?? harmonization algorithm has fg set to na when t2d=2. should we consider updating harmonization document ???
table(fulldata_sub_ancestry.sqrt$study,useNA = 'always')
summary(as.numeric(fulldata_sub_ancestry.sqrt$last_exam_hba1c)) #NA=31028
miss <- fulldata_sub_ancestry.sqrt[is.na(fulldata_sub_ancestry.sqrt$last_exam_hba1c),]
table(miss$study,useNA = 'always')

table(fulldata_sub_ancestry.sqrt$t2d, useNA='always')

## Reclassified based on FG alone
table(fulldata_sub_ancestry.sqrt$t2d, fulldata_sub_ancestry.sqrt$last_exam_fg >= 6.0, useNA='always')
table(fulldata_sub_ancestry.sqrt$study, fulldata_sub_ancestry.sqrt$last_exam_fg >= 6.0,fulldata_sub_ancestry.sqrt$t2d, useNA='always')

## Reclassified based on HbA1c alone
table(fulldata_sub_ancestry.sqrt$t2d, fulldata_sub_ancestry.sqrt$last_exam_hba1c >= 6.0, useNA='always')
table(fulldata_sub_ancestry.sqrt$study, fulldata_sub_ancestry.sqrt$last_exam_hba1c >= 6.0,fulldata_sub_ancestry.sqrt$t2d, useNA='always')

## Reclassified based on FG or HbA1c
table(fulldata_sub_ancestry.sqrt$t2d, fulldata_sub_ancestry.sqrt$last_exam_fg >= 6.0 | fulldata_sub_ancestry.sqrt$last_exam_hba1c >= 6.0, useNA='always')
table(fulldata_sub_ancestry.sqrt$study, fulldata_sub_ancestry.sqrt$last_exam_fg >= 6.0 | fulldata_sub_ancestry.sqrt$last_exam_hba1c >= 6.0,fulldata_sub_ancestry.sqrt$t2d, useNA='always')

fulldata_sub_ancestry.sqrt$t2d.who.pret2d <- NA

fulldata_sub_ancestry.sqrt$t2d.who.pret2d[(fulldata_sub_ancestry.sqrt$t2d==0)] = 0 # baseline all previous controls will be controls
fulldata_sub_ancestry.sqrt$t2d.who.pret2d[(fulldata_sub_ancestry.sqrt$t2d==1)] = 0 # baseline all previous prediabetes will be controls
fulldata_sub_ancestry.sqrt$t2d.who.pret2d[(fulldata_sub_ancestry.sqrt$t2d==2)] = 1 # all previous cases will be controls

# Reclassify controls as cases

fulldata_sub_ancestry.sqrt$t2d.who.pret2d[(fulldata_sub_ancestry.sqrt$last_exam_fg >= 6.0 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_fg)) 
                                  | (fulldata_sub_ancestry.sqrt$last_exam_hba1c >= 6.0 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_hba1c)) 
                                  & (fulldata_sub_ancestry.sqrt$t2d==0)] = 1

# Reclassify prediabetes as cases
fulldata_sub_ancestry.sqrt$t2d.who.pret2d[(fulldata_sub_ancestry.sqrt$last_exam_fg >= 6.0 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_fg)) 
                                  | (fulldata_sub_ancestry.sqrt$last_exam_hba1c >= 6.0 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_hba1c)) 
                                  & (fulldata_sub_ancestry.sqrt$t2d==1)] = 1

table(t2d.who.pret2d=fulldata_sub_ancestry.sqrt$t2d.who.pret2d,t2d=fulldata_sub_ancestry.sqrt$t2d,useNA="always")

###############################
#exclude sparse cells - t2d.who.pret2d - cluster ancestry


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

