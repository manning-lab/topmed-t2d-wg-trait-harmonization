

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

fulldata_sub_ancestry.sqrt$t2dnew <- NA
fulldata_sub_ancestry.sqrt$t2dnew[(fulldata_sub_ancestry.sqrt$t2d==0)] = 0 #when defining controls this line should be first
fulldata_sub_ancestry.sqrt$t2dnew[(fulldata_sub_ancestry.sqrt$last_exam_fg < 5.5 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_fg)) | (fulldata_sub_ancestry.sqrt$last_exam_hba1c < 5.4 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_hba1c)) & (fulldata_sub_ancestry.sqrt$t2d != 2)] = 0
fulldata_sub_ancestry.sqrt$t2dnew[(fulldata_sub_ancestry.sqrt$last_exam_fg >= 6.0 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_fg)) | (fulldata_sub_ancestry.sqrt$last_exam_hba1c >= 6.1 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_hba1c)) | (fulldata_sub_ancestry.sqrt$t2d==2)] = 1
fulldata_sub_ancestry.sqrt$t2dnew[fulldata_sub_ancestry.sqrt$last_exam_fg < 5.5 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_fg) & fulldata_sub_ancestry.sqrt$t2d == 1] = NA
fulldata_sub_ancestry.sqrt$t2dnew[fulldata_sub_ancestry.sqrt$last_exam_hba1c < 5.4 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_hba1c) & fulldata_sub_ancestry.sqrt$t2d == 1 ] = NA

# table to show how newdefn affects traditional definition
table(fulldata_sub_ancestry.sqrt$t2dnew,fulldata_sub_ancestry.sqrt$t2d,useNA = 'always')


## exclude excess controls and exclude sparse cells - new t2d definition PC cluster ancestry

table(fulldata_sub_ancestry.sqrt$t2dnew,fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt)
table(fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$t2dnew,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt)

# keep cases where n>=5 by ancestry & project
fulldata_sub_ancestry.sqrt$t2dnew_pccluster <- NA

fulldata_sub_ancestry.sqrt$t2dnew_pccluster[which((fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtAF") & (fulldata_sub_ancestry.sqrt$topmed_project %in% c("AA_CAC","AFGen","CFS","COPD","GeneSTAR","HyperGEN_GENOA","JHS","MESA","VTE","WHI")) & (fulldata_sub_ancestry.sqrt$t2dnew==1))] <- 1
fulldata_sub_ancestry.sqrt$t2dnew_pccluster[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtAS" & fulldata_sub_ancestry.sqrt$topmed_project %in% c("GenSalt","MESA","WHI") & fulldata_sub_ancestry.sqrt$t2dnew==1)] <- 1
fulldata_sub_ancestry.sqrt$t2dnew_pccluster[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtEU" & fulldata_sub_ancestry.sqrt$topmed_project %in% c("AFGen","Amish","CFS","COPD","FHS","GeneSTAR","GOLDN","MESA","VTE","WHI") & fulldata_sub_ancestry.sqrt$t2dnew==1)]  <- 1
fulldata_sub_ancestry.sqrt$t2dnew_pccluster[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project %in% c("MESA","SAFS","WHI") & fulldata_sub_ancestry.sqrt$t2dnew==1)]  <- 1
fulldata_sub_ancestry.sqrt$t2dnew_pccluster[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtSAS" & fulldata_sub_ancestry.sqrt$topmed_project=="SAS" & fulldata_sub_ancestry.sqrt$t2dnew==1)] <- 1

#exclude controls based on fg and/or a1c if available otherwise use age
fulldata_sub_ancestry.sqrt$t2dnew_pccluster[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtAF" & fulldata_sub_ancestry.sqrt$topmed_project=="COPD"  & fulldata_sub_ancestry.sqrt$t2dnew==0 &
 fulldata_sub_ancestry.sqrt$last_exam_age >= 50)] <- 0
fulldata_sub_ancestry.sqrt$t2dnew_pccluster[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtEU" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen"  & fulldata_sub_ancestry.sqrt$t2dnew==0 & fulldata_sub_ancestry.sqrt$last_exam_age >= 54)] <- 0 #n=1325
fulldata_sub_ancestry.sqrt$t2dnew_pccluster[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtEU" & fulldata_sub_ancestry.sqrt$topmed_project=="Amish"  & fulldata_sub_ancestry.sqrt$t2dnew==0 &
fulldata_sub_ancestry.sqrt$last_exam_fg < 4.45 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_fg))] <- 0 #n=217

fulldata_sub_ancestry.sqrt$t2dnew_pccluster[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtEU" & fulldata_sub_ancestry.sqrt$topmed_project=="COPD"  & fulldata_sub_ancestry.sqrt$t2dnew==0 &
fulldata_sub_ancestry.sqrt$last_exam_age >= 61)] <- 0 #n=2720
fulldata_sub_ancestry.sqrt$t2dnew_pccluster[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtEU" & fulldata_sub_ancestry.sqrt$topmed_project=="GeneSTAR"  & fulldata_sub_ancestry.sqrt$t2dnew==0 & fulldata_sub_ancestry.sqrt$last_exam_fg < 5.2 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_fg))] <- 0 #n=440

# controls from the rest of the studies
fulldata_sub_ancestry.sqrt$t2dnew_pccluster[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtAF" & fulldata_sub_ancestry.sqrt$topmed_project %in% c("AA_CAC","AFGen","CFS","GeneSTAR","HyperGEN_GENOA","JHS","MESA","VTE","WHI") & fulldata_sub_ancestry.sqrt$t2dnew==0)] <- 0
fulldata_sub_ancestry.sqrt$t2dnew_pccluster[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtAS" & fulldata_sub_ancestry.sqrt$topmed_project %in% c("GenSalt","MESA","WHI") & fulldata_sub_ancestry.sqrt$t2dnew==0)] <- 0
fulldata_sub_ancestry.sqrt$t2dnew_pccluster[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtEU" & fulldata_sub_ancestry.sqrt$topmed_project %in% c("CFS","FHS","GOLDN","MESA","VTE","WHI") & fulldata_sub_ancestry.sqrt$t2dnew==0)] <- 0
fulldata_sub_ancestry.sqrt$t2dnew_pccluster[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtHS" & fulldata_sub_ancestry.sqrt$topmed_project %in% c("MESA","SAFS","WHI") & fulldata_sub_ancestry.sqrt$t2dnew==0)] <- 0
fulldata_sub_ancestry.sqrt$t2dnew_pccluster[which(fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt=="c.sqrtSAS" & fulldata_sub_ancestry.sqrt$topmed_project %in% c("SAS") & fulldata_sub_ancestry.sqrt$t2dnew==0)] <- 0

table(fulldata_sub_ancestry.sqrt$t2dnew_pccluster,useNA='always')

table(fulldata_sub_ancestry.sqrt$t2dnew_pccluster,fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt)

table(fulldata_sub_ancestry.sqrt$t2dnew_pccluster,fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$cluster.ancestry.sqrt)

################ exclude excess controls and sparse cells - new t2d definition self report ancestry   ################
table(fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$t2dnew,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$population)
table(fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$t2dnew,fulldata_sub_ancestry.sqrt$population)

fulldata_sub_ancestry.sqrt$t2dnew_sr <- NA

fulldata_sub_ancestry.sqrt$t2dnew_sr[which(fulldata_sub_ancestry.sqrt$population=="AF" & fulldata_sub_ancestry.sqrt$topmed_project %in% c("AA_CAC","CFS","COPD","GeneSTAR","HyperGEN_GENOA","JHS","MESA","VTE","WHI") & fulldata_sub_ancestry.sqrt$t2dnew==1)] <-1
fulldata_sub_ancestry.sqrt$t2dnew_sr[which(fulldata_sub_ancestry.sqrt$population=="AS" & fulldata_sub_ancestry.sqrt$topmed_project %in% c("GenSalt","MESA","WHI") & fulldata_sub_ancestry.sqrt$t2dnew==1)] <-1
fulldata_sub_ancestry.sqrt$t2dnew_sr[which(fulldata_sub_ancestry.sqrt$population=="EU" & fulldata_sub_ancestry.sqrt$topmed_project %in% c("AFGen","Amish","CFS","COPD","FHS","GeneSTAR","GOLDN","MESA","VTE","WHI") & fulldata_sub_ancestry.sqrt$t2dnew==1)] <-1
fulldata_sub_ancestry.sqrt$t2dnew_sr[which(fulldata_sub_ancestry.sqrt$population=="HS" & fulldata_sub_ancestry.sqrt$topmed_project %in% c("MESA","SAFS","WHI") & fulldata_sub_ancestry.sqrt$t2dnew==1)] <-1
fulldata_sub_ancestry.sqrt$t2dnew_sr[which(fulldata_sub_ancestry.sqrt$population=="SAS" & fulldata_sub_ancestry.sqrt$topmed_project=="SAS"  &fulldata_sub_ancestry.sqrt$t2dnew==1)] <-1

fulldata_sub_ancestry.sqrt$t2dnew_sr[which(fulldata_sub_ancestry.sqrt$population=="AF" & fulldata_sub_ancestry.sqrt$topmed_project=="COPD"  & fulldata_sub_ancestry.sqrt$t2dnew==0 & fulldata_sub_ancestry.sqrt$last_exam_age >= 50)]  <-0 #n=1639
fulldata_sub_ancestry.sqrt$t2dnew_sr[which(fulldata_sub_ancestry.sqrt$population=="EU" & fulldata_sub_ancestry.sqrt$topmed_project=="AFGen"  & fulldata_sub_ancestry.sqrt$t2dnew==0 & fulldata_sub_ancestry.sqrt$last_exam_age >= 54)] <- 0 #n=1316
fulldata_sub_ancestry.sqrt$t2dnew_sr[which(fulldata_sub_ancestry.sqrt$population=="EU" & fulldata_sub_ancestry.sqrt$topmed_project=="Amish"  & fulldata_sub_ancestry.sqrt$t2dnew==0 & fulldata_sub_ancestry.sqrt$last_exam_fg < 4.45 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_fg))] <- 0 #n=217
fulldata_sub_ancestry.sqrt$t2dnew_sr[which(fulldata_sub_ancestry.sqrt$population=="EU" & fulldata_sub_ancestry.sqrt$topmed_project=="COPD"  & fulldata_sub_ancestry.sqrt$t2dnew==0 & fulldata_sub_ancestry.sqrt$last_exam_age >= 61)] <-0 #n=2727
fulldata_sub_ancestry.sqrt$t2dnew_sr[which(fulldata_sub_ancestry.sqrt$population=="EU" & fulldata_sub_ancestry.sqrt$topmed_project=="GeneSTAR"  & fulldata_sub_ancestry.sqrt$t2dnew==0 & fulldata_sub_ancestry.sqrt$last_exam_fg < 5.2 & !is.na(fulldata_sub_ancestry.sqrt$last_exam_fg))] <- 0 #n=445

fulldata_sub_ancestry.sqrt$t2dnew_sr[which(fulldata_sub_ancestry.sqrt$population=="AF" & fulldata_sub_ancestry.sqrt$topmed_project %in% c("AA_CAC","CFS","GeneSTAR","HyperGEN_GENOA","JHS","MESA","VTE","WHI") & fulldata_sub_ancestry.sqrt$t2dnew==0)] <- 0
fulldata_sub_ancestry.sqrt$t2dnew_sr[which(fulldata_sub_ancestry.sqrt$population=="AS" & fulldata_sub_ancestry.sqrt$topmed_project %in% c("GenSalt","MESA","WHI") & fulldata_sub_ancestry.sqrt$t2dnew==0)] <- 0
fulldata_sub_ancestry.sqrt$t2dnew_sr[which(fulldata_sub_ancestry.sqrt$population=="EU" & fulldata_sub_ancestry.sqrt$topmed_project %in% c("CFS","FHS","GOLDN","MESA","VTE","WHI") & fulldata_sub_ancestry.sqrt$t2dnew==0)] <- 0
fulldata_sub_ancestry.sqrt$t2dnew_sr[which(fulldata_sub_ancestry.sqrt$population=="HS" & fulldata_sub_ancestry.sqrt$topmed_project %in% c("MESA","SAFS","WHI") & fulldata_sub_ancestry.sqrt$t2dnew==0)] <- 0
fulldata_sub_ancestry.sqrt$t2dnew_sr[which(fulldata_sub_ancestry.sqrt$population=="SAS" & fulldata_sub_ancestry.sqrt$topmed_project %in% c("SAS") & fulldata_sub_ancestry.sqrt$t2dnew==0)] <- 0


table(fulldata_sub_ancestry.sqrt$t2dnew_sr,useNA='always')

table(fulldata_sub_ancestry.sqrt$t2dnew_sr,fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$population)

table(fulldata_sub_ancestry.sqrt$t2dnew_sr,fulldata_sub_ancestry.sqrt$topmed_project,fulldata_sub_ancestry.sqrt$sex,fulldata_sub_ancestry.sqrt$population)


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

