
## T2D Harmonization  
## GENESIS needs sex as M,F
# FamilyID = SUBJECT_ID or individual_id depending on dataset
## PaternalID = 0 if not provided
## MaternalID = 0 if not provided
# create BMI variable
# BMI = 703·weight(lb)/height2(in2)
##

args <- commandArgs(trailingOnly=T)
f.dir <- args[1]
out.pref <- args[2]

# ### testing inputs ####
#  f.dir <- "/N/dc2/scratch/wesselj/OriginalFiles"
# # source.file <- "/Users/tmajaria/Documents/projects/topmed/code/topmed-traitHarmonization/Harmonization.19JAN2017.GitHub.SourceFiles.R"
#  out.pref <- "/N/dc2/scratch/wesselj/OriginalFiles"
# # ########################

# load all of the phenotype info through the source file
source("Freeze5b.T2D.inputfiles.R")
dat <- get_pheno_data(f.dir)

map <- dat$map

afp <- dat$afp
afvub <- dat$afvub
afccaf <- dat$afccaf
amish <- dat$amish
aric <- dat$aric
cfs <- dat$cfs
chs <- dat$chs
copd <- dat$copd
dhs <- dat$dhs
dhsPed <- dat$dhsPed
fhs <- dat$fhs
genestar <- dat$genestar
genoa <- dat$genoa
gensalt <- dat$gensalt
goldn <- dat$goldn
hypergen <- dat$hypergen
jhs <- dat$jhs
mesa <- dat$mesa
mesafam <- dat$mesafam
raw.HVH <- dat$raw.HVH
raw.MGH <- dat$raw.MGH
raw.MGHPed <- dat$raw.MGHPed
raw.VUdaw <- dat$raw.VUdaw
raw.VUdawPed <- dat$raw.VUdawPed
safs <- dat$safs
#safs.ids <- dat$safs.ids
sas <- dat$sas
whi <- dat$whi

rm(dat)

####### CREATE STANDARDIZED PHENOTYPE DATASETS TO POOL TOGETHER AND THEN MERGE WITH THE MAP FILE ###########


####################### AFib CCAF ##############################################
####################### AFib CCAF ##############################################

#
## AFib CCAF
# downlaoded new dataset 5AUG2017
# 7DEC2017 checked for new files - none
#names(afccaf)

# recode & check variable names & distributions
# rename column for diabetes
names(afccaf)[names(afccaf) == "diabetes"] <- "t2d"
# recode t2d status
afccaf$t2d[afccaf$t2d == 'yes'] = 2
afccaf$t2d[afccaf$t2d == 'no'] = 0
table(afccaf$t2d,useNA='always')
#afccaf = subset(afccaf, subset = t2d %in% c(0,1,2))
table(afccaf$race,useNA='always')
table(afccaf$ethnicity,useNA='always')
afccaf$ancestry[afccaf$race == 'white' & afccaf$ethnicity == 'no'] = "EU"
afccaf$ancestry[afccaf$race == 'white' & afccaf$ethnicity == 'yes'] = "HS"
table(afccaf$ancestry,useNA='always')
table(afccaf$sex,useNA='always')
afccaf$origsex = afccaf$sex
afccaf$sex[afccaf$sex == 'male'] = 'M'
afccaf$sex[afccaf$sex == 'female'] = 'F'
with(afccaf,table(origsex,sex,useNA='always'))
summary(afccaf$age)
afccaf$last_exam_age = afccaf$age
#afccaf = subset(afccaf, age >= 25) #n=352, lose 3 individuals
summary(afccaf$height)
summary(afccaf$weight)
afccaf$last_exam_bmi = ((703*afccaf$weight)/(afccaf$height*afccaf$height))
summary(afccaf$last_exam_bmi) # ! low BMI
afccaf$last_exam_fg = NA
afccaf$last_exam_hba1c = NA
afccaf$last_exam_t2d_treatment = NA
afccaf$last_exam_visit = NA
afccaf$t2d_age = NA
afccaf$t2d_bmi = NA 
afccaf$FamilyID = afccaf$SUBJECT_ID
afccaf$PaternalID = 0
afccaf$MaternalID = 0
afccaf$sequenced = NA
afccaf$study = "CCAF"
afccaf$unique_id <-  paste(afccaf$study,afccaf$SUBJECT_ID, sep = "_")
afccaf$individual_id <- afccaf$SUBJECT_ID
afccaf$study_ancestry <- paste(afccaf$study,afccaf$ancestry, sep = "_")
afccaf$JWsource = "dbGaP"


afccaf <- afccaf[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
                    'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
                    'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]
#n=363
# 20 missing T2D

table(duplicated(afccaf$unique_id)) #0
####################### AFib CCAF ##############################################
####################### AFib CCAF ##############################################


####################### AFib Partners ##############################################
####################### AFib Partners ##############################################
#
## AFib Partners
# 7DEC2017 using downloaded newly released data 6OCT2017
names(afp)
afp$sequenced = NA

# recode & check variable names & distributions
# rename column for diabetes
names(afp)[names(afp) == "diabetes"] <- "t2d"
# recode t2d status
afp$t2d[afp$t2d == 'yes'] = 2
afp$t2d[afp$t2d == 'no'] = 0
#afp = subset(afp, subset = t2d %in% c(0,1,2))
table(afp$t2d,useNA='always') # no missing t2d

# recode ancestry
table(afp$race,useNA='always')
# asian    black hispanic    other    white     <NA> 
#   1        4        1        1      120        1 
table(afp$ethnicity,useNA='always')
# no  yes <NA> 
#   120    1    7 
afp$ancestry[afp$race == 'black' & afp$ethnicity == 'no'] = "AF"
afp$ancestry[afp$race == 'white' & afp$ethnicity == 'no'] = "EU"
afp$ancestry[afp$race == 'white'] = "EU"
afp$ancestry[afp$race == 'asian' & afp$ethnicity == 'no'] = "AS"
afp$ancestry[afp$race == 'hispanic' & afp$ethnicity == 'yes'] = "HS"
afp$ancestry[afp$race == 'other' & afp$ethnicity == 'no'] = "OTHER"
afp$ancestry[is.na(afp$race) & is.na(afp$ethnicity)] = "UNKNOWN"
table(afp$ancestry,useNA='always')
with(afp,table(ancestry,ethnicity,race,useNA='always')) # 1  missing race & ancestry data

# recode sex
table(afp$sex,useNA='always')
afp$origsex = afp$sex
afp$sex[afp$sex == 'male'] = 'M'
afp$sex[afp$sex == 'female'] = 'F'
with(afp,table(origsex,sex,useNA='always'))

# recode age
summary(afp$age)
afp$last_exam_age = afp$age
#afp = subset(afp, last_exam_age >= 25) #drop 2 individuals
afp = subset(afp, weight <= 401) #drop 1 individual
## exclude 1 individual with weight = 400+ ##

# recode weight
afp$weight <- as.numeric(afp$weight)
class(afp$weight)
summary(afp$height)
summary(afp$weight)
afp$last_exam_bmi = ((703*afp$weight)/(afp$height*afp$height))
summary(afp$last_exam_bmi)
afp$last_exam_fg = NA
afp$last_exam_hba1c = NA
afp$last_exam_t2d_treatment = NA
afp$last_exam_visit = NA
afp$t2d_age = NA
afp$t2d_bmi = NA
afp$FamilyID = afp$SUBJECT_ID
afp$PaternalID = 0
afp$MaternalID = 0
afp$study = "Partners"
afp$unique_id <-  paste(afp$study,afp$SUBJECT_ID, sep = "_")
afp$individual_id <- afp$SUBJECT_ID
afp$study_ancestry <- paste(afp$study,afp$ancestry, sep = "_")
afp$JWsource = "dbGaP"

afp <- afp[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
              'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
              'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]
#n=127
table(duplicated(afp$unique_id)) #0
####################### AFib Partners ##############################################
####################### AFib Partners ##############################################

####################### VU_Ben ##############################################
####################### VU_Ben ##############################################

#
## AFib VU_Ben
# downlaoded new dataset 5AUG2017
# 7DEC2017 checked for new files - none
names(afvub)
afvub$sequenced = NA


# recode & check variable names & distributions
# rename column for diabetes
names(afvub)[names(afvub) == "diabetes"] <- "t2d"
# recode t2d status
afvub$t2d[afvub$t2d == 'yes'] = 2
afvub$t2d[afvub$t2d == 'no'] = 0
table(afvub$t2d,useNA='always')
#afvub = subset(afvub, subset = t2d %in% c(0,1,2))

table(afvub$race,useNA='always')
table(afvub$ethnicity,useNA='always')
afvub$ancestry[afvub$race == 'white' & afvub$ethnicity == 'no'] = "EU"
afvub$ancestry[afvub$ethnicity == 'yes'] = "HS"
table(afvub$ancestry,useNA='always')
table(afvub$sex,useNA='always')
afvub$origsex = afvub$sex
afvub$sex[afvub$sex == 'male'] = 'M'
afvub$sex[afvub$sex == 'female'] = 'F'
with(afvub,table(origsex,sex,useNA='always'))
afvub$last_exam_bmi = ((703*afvub$weight)/(afvub$height*afvub$height))
summary(afvub$last_exam_bmi)
summary(afvub$age)
afvub$last_exam_age = afvub$age
#afvub = subset(afvub, age >= 25) #n=117, no drops
summary(afvub$height)
summary(afvub$weight)
afvub$last_exam_bmi = ((703*afvub$weight)/(afvub$height*afvub$height))
summary(afvub$last_exam_bmi)
afvub$last_exam_fg = NA
afvub$last_exam_hba1c = NA
afvub$last_exam_t2d_treatment = NA
afvub$last_exam_visit = NA
afvub$t2d_age = NA
afvub$t2d_bmi = NA
afvub$FamilyID = afvub$SUBJECT_ID
afvub$PaternalID = 0
afvub$MaternalID = 0
afvub$study = "VAFAR"
afvub$unique_id <-  paste(afvub$study,afvub$SUBJECT_ID, sep = "_")
afvub$individual_id <- afvub$SUBJECT_ID
afvub$study_ancestry <- paste(afvub$study,afvub$ancestry, sep = "_")
afvub$JWsource = "dbGaP"

afvub <- afvub[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
                  'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
                  'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]
#n=171
table(duplicated(afvub$unique_id)) #0
####################### VU_Ben ##############################################
####################### VU_Ben ##############################################


####################### Amish ##############################################
####################### Amish ##############################################

## Amish
## Amish has 2 pairs of twins. they left the co-twin out "were removed on our end as
#twins caused some issues with some of our phenotypes in our in house analysis (but not lipids)"
# 15AUG2017

names(amish)

# recode & check variable names & distributions
table(amish$sequenced,useNA='always')

# rename column for diabetes
names(amish)[names(amish) == "Diabetes_Status"] <- "t2d"
# recode t2d status
amish$t2d[amish$t2d == 1] = 2
amish$t2d[amish$t2d == 0] = 0
table(amish$t2d,useNA='always')
#amish = subset(amish, subset = t2d %in% c(0,1,2))

amish$JWsource = "dbGaP_Ex"
amish$ancestry = 'EU'
table(amish$ancestry,useNA='always')
table(amish$sex,useNA='always')
amish$origsex = amish$sex
amish$sex[amish$sex == 1] = 'M'
amish$sex[amish$sex == 2] = 'F'
with(amish,table(origsex,sex,useNA='always'))

summary(amish$T2D_age) #exclude by t2d age first then see if we need to exclude by age 
summary(amish$last_exam_age)
#amish = subset(amish, last_exam_age >= 25) #n=962

summary(amish$last_exam_BMI)
summary(amish$T2D_BMI)
amish$last_exam_bmi = amish$last_exam_BMI
amish$last_exam_fg = amish$last_exam_FG
amish$last_exam_hba1c = amish$last_exam_HbA1c
amish$last_exam_t2d_treatment = amish$last_exam_T2D_treatment
amish$t2d_age = amish$T2D_age
amish$t2d_bmi = amish$T2D_BMI
amish$last_exam_visit = NA
amish$FamilyID = amish$FAMID
amish$PaternalID = amish$FATHER
amish$MaternalID = amish$MOTHER
amish$study = "Amish"
amish$unique_id <-  paste(amish$study,amish$SUBJID, sep = "_")
amish$individual_id <- amish$SUBJID
amish$study_ancestry <- paste(amish$study,amish$ancestry, sep = "_")


amish <- amish[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
                  'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
                  'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]
#n=1013
table(duplicated(amish$unique_id)) #0
####################### Amish ##############################################
####################### Amish ##############################################

####################### ARIC ##############################################
####################### ARIC ##############################################

## ARIC
# NEED TO FIX CODING OF last exam t2d treatment (Y/N/U) if used in analyses
names(aric)

# recode & check variable names & distributions
# rename column for diabetes
names(aric)[names(aric) == "T2D"] <- "t2d"
table(aric$t2d,useNA='always')
#aric = subset(aric, subset = t2d %in% c(0,1,2))

table(aric$ancestry,useNA='always')
table(aric$sex,useNA='always')
aric$origsex = aric$sex
aric$sex[aric$sex == 1] = 'M'
aric$sex[aric$sex == 2] = 'F'
with(aric,table(origsex,sex,useNA='always'))
summary(aric$last_exam_age)
aric$last_exam_bmi = aric$last_exam_BMI
summary(aric$last_exam_bmi) #!low BMI
aric$last_exam_fg = NA
aric$last_exam_t2d_treatment = aric$last_exam_T2D_treatment
aric$last_exam_hba1c = aric$last_exam_HbA1c
aric$t2d_age = aric$T2D_age
aric$t2d_bmi = aric$T2D_BMI
summary(aric$last_exam_hba1c)
summary(aric$t2d_age)
summary(aric$t2d_bmi)
table(aric$last_exam_visit)
table(aric$last_exam_t2d_treatment)
aric$FamilyID = aric$Family_ID
aric$PaternalID = aric$Father_ID
aric$MaternalID = aric$Mother_ID
aric$study = "ARIC"
aric$unique_id <-  paste(aric$study,aric$Individual_ID, sep = "_")
aric$individual_id <- aric$Individual_ID
aric$study_ancestry <- paste(aric$study,aric$ancestry, sep = "_")
aric$JWsource = "dbGaP_Ex"


aric <- aric[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
                'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
                'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]
#n=4230 49 missing T2D
table(duplicated(aric$unique_id)) #0


####################### ARIC ##############################################
####################### ARIC ##############################################


####################### CFS ##############################################
####################### CFS ##############################################

## CFS
# imported data 5AUG2017
# CFS has 1 MZ twin and both are included in pheno, although they end up being excluded because of age
# 15AUG2017

#T2D n=880 NA
names(cfs)
#add a check at each study to check that ID is in map file 
table(cfs$Sequenced, cfs$IID %in% map$unique_subject_key)

# recode & check variable names & distributions
table(cfs$Sequenced,useNA='always')
# 0    1 <NA> 
#  1537  994    0 
names(cfs)[names(cfs) == "Sequenced"] <- "sequenced"

# rename column for diabetes
names(cfs)[names(cfs) == "T2D"] <- "t2d"
table(cfs$t2d,useNA='always')
#cfs = subset(cfs, subset = t2d %in% c(0,1,2))

table(cfs$Population,useNA='always')
cfs$ancestry[cfs$Population == 'CFS-blacks'] = "AF"
cfs$ancestry[cfs$Population == 'CFS-whites'] = "EU"
cfs$ancestry[is.na(cfs$Population)] = "UNKNOWN"
with(cfs,table(Population,ancestry,useNA='always'))
table(cfs$ancestry,useNA = 'always')
table(cfs$sex,useNA='always')
cfs$origsex = cfs$sex
cfs$sex[cfs$sex == 1] = 'M'
cfs$sex[cfs$sex == 2] = 'F'
with(cfs,table(origsex,sex,useNA='always'))

summary(cfs$T2D_age)
#cfs = cfs[(cfs$T2D_age >= 25 | is.na(cfs$T2D_age)),] #n=723 drop 4
summary(cfs$last_exam_age)
#cfs = subset(cfs, last_exam_age >= 25) #n=727 drop 267

summary(cfs$last_exam_BMI)
cfs$last_exam_bmi = cfs$last_exam_BMI
cfs$last_exam_fg = cfs$FastingGlucose
cfs$last_exam_hba1c = NA
cfs$last_exam_t2d_treatment = cfs$last_exam_T2D_treatment
cfs$t2d_age = cfs$T2D_age
cfs$t2d_bmi = cfs$T2D_BMI
cfs$last_exam_visit = NA
cfs$FamilyID = cfs$FamID
cfs$PaternalID = cfs$FID
cfs$MaternalID = cfs$MID
cfs$study = "CFS"
cfs$unique_id <-  paste(cfs$study,cfs$IID, sep = "_") 
cfs$individual_id <- cfs$IID
cfs$study_ancestry <- paste(cfs$study,cfs$ancestry, sep = "_")
cfs$JWsource = "dbGaP_Ex"


cfs <- cfs[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
              'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
              'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]
#n=2531 880 missing T2D 
table(duplicated(cfs$unique_id)) #5 probably not sequenced??? make sure we get rid of them at merge
table(cfs$sequenced, cfs$unique_id %in% map$unique_subject_key)

# remove rows that are absolute duplicates
cfs.dups <- cfs$unique_id[duplicated(cfs$unique_id)]
cfs.final <- cfs[!(cfs$unique_id %in% cfs.dups),]
for (id in unique(cfs.dups)){
  my.rows <- cfs[cfs$unique_id == id,]
  is_equal = T
  for (col.id in seq(1,NCOL(my.rows))){
    this.col <- my.rows[,col.id]
    # change nas to some string
    this.col[is.na(this.col)] <- "This is na"
    
    # see how many unique elements we ahve
    nu <- length(unique(this.col)) > 1
    
    if (nu){
      is_equal = F
      break
    }
  }
  if(is_equal){
    cfs.final <- rbind(cfs.final,my.rows[1,])
  } else {
    cfs.final <- rbind(cfs.final,my.rows)
  }
}

#n=2531 still; individuals are not eliminated!!
table(duplicated(cfs.final$unique_id)) #5

cfs.final.dup <- cfs.final[duplicated(cfs.final$unique_id),]

####################### CFS ##############################################
####################### CFS ##############################################


####################### CHS ##############################################
####################### CHS ##############################################


#
# CHS
# 
names(chs)
# recode & check variable names & distributions
table(chs$sequenced,useNA='always')

# rename column for diabetes
names(chs)[names(chs) == "T2D"] <- "t2d"
table(chs$t2d,useNA='always')
#chs = subset(chs, subset = t2d %in% c(0,1,2))

table(chs$race01,useNA='always')
# 1    2    3    4    5 <NA> 
#   3123  785    5    3   13    0 
# Per email from J Floyd (27AUG2017) 3 is American Indian/Alaskan and 4 is Asian/Pacific Islander
chs$ancestry[chs$race01 == 2] = 'AF'
chs$ancestry[chs$race01 == 1] = 'EU'
chs$ancestry[chs$race01 == 5] = 'OTHER'
chs$ancestry[chs$race01 == 3] = 'AMR'
chs$ancestry[chs$race01 == 4] = 'AS'
with(chs,table(race01,ancestry,useNA='always'))
table(chs$ancestry,useNA='always')
table(chs$sex,useNA='always')
chs$sex[chs$sex == 1] = 'M'
chs$sex[chs$sex == 2] = 'F'
table(chs$sex,useNA='always')
summary(chs$last_exam_age)
summary(chs$last_exam_BMI)
summary(chs$last_exam_FG)
summary(chs$T2D_age)
chs$last_exam_bmi = chs$last_exam_BMI
chs$last_exam_fg = chs$last_exam_FG
chs$last_exam_hba1c = NA
chs$t2d_age = chs$T2D_age
chs$t2d_bmi = chs$T2D_BMI
chs$last_exam_t2d_treatment = chs$last_exam_T2D_treatment
chs$last_exam_visit = NA
chs$study = "CHS"
chs$unique_id <-  paste(chs$study,chs$SampleID, sep = "_")
chs$individual_id <- chs$SampleID
chs$study_ancestry <- paste(chs$study,chs$ancestry, sep = "_")
chs$JWsource = "dbGaP_Ex"


chs <- chs[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
              'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
              'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]
#n=3929
table(duplicated(chs$unique_id)) #0
####################### CHS ##############################################
####################### CHS ##############################################


####################### COPD Gene Sample & C1 & C2 ##############################################
####################### COPD Gene Sample & C1 & C2 ##############################################


##COPD Gene Sample & C1 & C2
names(copd)
copd$sequenced = NA
# recode & check variable names & distributions
# rename column for diabetes
names(copd)[names(copd) == "Diabetes"] <- "t2d"
# recode t2d status
copd$t2d[copd$t2d == 1] = 2
copd$t2d[copd$t2d == 0] = 0
table(copd$t2d,useNA='always')
#copd = subset(copd, subset = t2d %in% c(0,1,2))

table(copd$race,useNA='always')
table(copd$ethnic,useNA='always')
copd$ancestry[copd$race == 2 & copd$ethnic == 2] = "AF"
copd$ancestry[copd$race == 1 & copd$ethnic == 2] = "EU"
table(copd$ancestry,useNA='always')
table(copd$gender,useNA='always')
copd$sex = copd$gender
copd$sex[copd$sex == 1] = 'M'
copd$sex[copd$sex == 2] = 'F'
with(copd,table(gender,sex,useNA='always'))
summary(copd$Age_Enroll)
copd$last_exam_age = copd$Age_Enroll
copd$last_exam_bmi = copd$BMI
summary(copd$last_exam_bmi) # low BMI - exclude?
copd$last_exam_fg = NA
copd$last_exam_hba1c = NA
copd$last_exam_t2d_treatment = NA
copd$last_exam_visit = NA
copd$t2d_age = NA
copd$t2d_bmi = NA
copd$FamilyID = copd$SUBJECT_ID
copd$PaternalID = 0
copd$MaternalID = 0
copd$study = "COPDGene"
copd$unique_id <-  paste(copd$study,copd$SUBJECT_ID, sep = "_")
copd$individual_id <- copd$SUBJECT_ID
copd$study_ancestry <- paste(copd$study,copd$ancestry, sep = "_")
copd$JWsource = "dbGaP"

copd <- copd[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
                'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
                'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]
#n=10371 
#n=102 with missing phenotype data
table(duplicated(copd$unique_id)) #0
####################### COPD Gene Sample & C1 & C2 ##############################################
####################### COPD Gene Sample & C1 & C2 ##############################################


####################### DHS ##############################################
####################### DHS ##############################################

## DHS download
names(dhs)
dhsPed <- dhsPed[,c("FAMILY_ID",  "SUBJECT_ID", "FATHER", "MOTHER")]
dhs <- merge(dhs,dhsPed,by.x='SUBJECT_ID',by.y='SUBJECT_ID',all.x=TRUE) #n=405 & 33 variables
rm(dhsPed)

dhs$sequenced = NA
# recode & check variable names & distributions
# rename column for diabetes
names(dhs)[names(dhs) == "DIABETES"] <- "t2d"
# recode t2d status
dhs$t2d[dhs$t2d == 1] = 2
dhs$t2d[dhs$t2d == 0] = 0
table(dhs$t2d,useNA='always')
#dhs = subset(dhs, subset = t2d %in% c(0,1,2))
dhs$JWsource = "dbGaP_Ex_General"
table(dhs$RACE,useNA='always')
dhs$ancestry[dhs$RACE == 'AA'] = 'AF'
table(dhs$ancestry,useNA='always')
table(dhs$SEX,useNA='always')
dhs$sex = dhs$SEX
table(dhs$sex,useNA='always')
dhs$AGE<- as.numeric(dhs$AGE)
summary(dhs$AGE)
summary(dhs$BMI)
summary(dhs$HBA1C)
summary(dhs$GLUCOSE)
dhs$last_exam_fg = dhs$GLUCOSE*0.0555  # !! convert !!
summary(dhs$last_exam_fg)
dhs$t2d_age = NA
dhs$t2d_bmi = NA
dhs$last_exam_age = dhs$AGE
dhs$last_exam_bmi = dhs$BMI
dhs$last_exam_hba1c = dhs$HBA1C
dhs$last_exam_t2d_treatment = NA
dhs$FamilyID = dhs$FAMILY_ID
dhs$PaternalID = dhs$FATHER
dhs$MaternalID = dhs$MOTHER

dhs$study = "DHS"
dhs$unique_id <-  paste(dhs$study,dhs$SUBJECT_ID, sep = "_")
dhs$individual_id <- dhs$SUBJECT_ID
dhs$study_ancestry <- paste(dhs$study,dhs$ancestry, sep = "_")

dhs <- dhs[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
              'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
              'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]
#n=405
table(duplicated(dhs$unique_id)) #0
####################### DHS ##############################################
####################### DHS ##############################################


####################### FHS ##############################################
####################### FHS ##############################################
##FHS
# #13AUG2017 currently n=4024 since 1 individual has not been WGSed
# 15AUG2017 FHS has 6 MZ twin pairs. 
#
###############################################


# recode & check variable names & distributions
table(fhs$sequenced,useNA='always') 

# rename column for diabetes
names(fhs)[names(fhs) == "T2D"] <- "t2d"
table(fhs$t2d,useNA='always')
#fhs = subset(fhs, subset = t2d %in% c(0,1,2))

table(fhs$sex,useNA='always')
fhs$origsex = fhs$sex
fhs$sex[fhs$sex == 1] = 'M'
fhs$sex[fhs$sex == 2] = 'F'
with(fhs,table(sex,origsex,useNA='always'))
summary(fhs$last_exam_age)
fhs$last_exam_bmi = fhs$last_exam_BMI
summary(fhs$last_exam_bmi)
fhs$last_exam_fg = fhs$last_exam_FG
summary(fhs$last_exam_fg)
fhs$last_exam_hba1c = fhs$last_exam_HbA1c
summary(fhs$last_exam_hba1c)
fhs$last_exam_t2d_treatment = fhs$last_exam_T2D_treatment
table(fhs$last_exam_t2d_treatment,useNA='always')
fhs$t2d_age = fhs$T2D_age
summary(fhs$t2d_age)
fhs$t2d_bmi = fhs$T2D_BMI
summary(fhs$t2d_bmi)
table(fhs$last_exam_visit,useNA = 'always')
fhs$FamilyID = fhs$shareid
fhs$PaternalID = 0
fhs$MaternalID = 0
fhs$ancestry = 'EU'
table(fhs$ancestry,useNA='always')

fhs$study = "FHS"
fhs$unique_id <-  paste(fhs$study,fhs$shareid, sep = "_")
fhs$individual_id <- fhs$shareid


table(fhs$idtype,useNA = 'always')
fhs$studygen = ifelse(is.na(fhs$idtype),"fhsNA",fhs$idtype)
fhs$studygen[fhs$idtype == 0] = "fhs0"
fhs$studygen[fhs$idtype == 1] = "fhs1"
fhs$studygen[fhs$idtype == 3] = "fhs3"
table(fhs$studygen,useNA = 'always')
fhs$study_ancestry <- paste(fhs$studygen,fhs$ancestry, sep = "_")
table(fhs$study_ancestry,useNA = 'always')

fhs$JWsource = "dbGaP_Ex"

fhs <- fhs[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
              'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
              'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]

#n=14329
table(duplicated(fhs$unique_id)) #0
####################### FHS ##############################################
####################### FHS ##############################################





####################### GenSalt  ##############################################
####################### GenSalt ##############################################


# GenSalt 
#

# recode & check variable names & distributions
table(gensalt$Sequenced,useNA='always')
names(gensalt)[names(gensalt) == "Sequenced"] <- "sequenced"

# rename column for diabetes
names(gensalt)[names(gensalt) == "T2D"] <- "t2d"
table(gensalt$t2d,useNA='always')
#gensalt = subset(gensalt, subset = t2d %in% c(0,1,2))

gensalt$ancestry = 'AS'
table(gensalt$ancestry,useNA='always')
table(gensalt$SEX,useNA='always')
gensalt$sex[gensalt$SEX == 1] = 'M'
gensalt$sex[gensalt$SEX == 2] = 'F'
with(gensalt,table(SEX,sex,useNA='always'))
#describeBy(gensalt$Last_exam_age, gensalt$Sequenced, mat = TRUE)
summary(gensalt$Last_exam_age,useNA='always') #! younger than 18
#gensalt = subset(gensalt, Last_exam_age >= 18) 
gensalt$Last_exam_BMI <- as.numeric(as.character(gensalt$Last_exam_BMI)) #convert class character to numeric
summary(gensalt$Last_exam_BMI)
summary(gensalt$T2D_age) #drop 0 individuals
gensalt$last_exam_age = gensalt$Last_exam_age
gensalt$last_exam_bmi = gensalt$Last_exam_BMI
gensalt$last_exam_fg = gensalt$Last_exam_FG
gensalt$last_exam_hba1c = gensalt$Last_exam_HbA1c
gensalt$last_exam_t2d_treatment = gensalt$Last_exam_T2D_treatment
gensalt$t2d_age = gensalt$T2D_age
gensalt$t2d_bmi = gensalt$T2D_BMI
table(gensalt$Last_exam_visit,useNA = 'always')
gensalt$FamilyID = gensalt$Family_ID
gensalt$PaternalID = gensalt$Father_ID
gensalt$MaternalID = gensalt$Mother_ID
gensalt$individual_id <- gensalt$Individual_ID

gensalt$study = "GenSalt"
gensalt$unique_id <-  paste(gensalt$study,gensalt$Individual_ID, sep = "_")
gensalt$individual_id <- gensalt$Individual_ID
gensalt$study_ancestry <- paste(gensalt$study,gensalt$ancestry, sep = "_")
gensalt$JWsource = "dbGaP_Ex"

gensalt <- gensalt[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
                      'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
                      'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]
#n=1906 93 missing T2D
table(duplicated(gensalt$unique_id)) #0

####################### GenSalt Phase 2 ##############################################
####################### GenSalt Phase 2 ##############################################


####################### GeneSTAR ##############################################
####################### GeneSTAR ##############################################

#
# GeneSTAR
#

# recode & check variable names & distributions
table(genestar$sequenced,useNA='always')

# rename column for diabetes
names(genestar)[names(genestar) == "T2D"] <- "t2d"
table(genestar$t2d,useNA='always')
#genestar = subset(genestar, subset = t2d %in% c(0,1,2))

table(genestar$ancestry,useNA='always')
table(genestar$sex,useNA='always')
genestar$origsex = genestar$sex
genestar$sex[genestar$origsex == 1] = 'M'
genestar$sex[genestar$origsex == 2] = 'F'
with(genestar,table(origsex,sex,useNA='always'))
summary(genestar$last_exam_age,useNA='always')
summary(genestar$T2D_age)
#genestar = subset(genestar, T2D_age >= 25) #n=, drop XX individuals
summary(genestar$last_exam_BMI) # ! low BMI & missing
genestar$last_exam_bmi = genestar$last_exam_BMI
genestar$last_exam_fg = genestar$last_exam_FG
genestar$last_exam_hba1c = genestar$last_exam_HbA1C
genestar$last_exam_t2d_treatment = genestar$last_exam_T2D_treatment
genestar$t2d_age = genestar$T2D_age
genestar$t2d_bmi = genestar$T2D_BMI
table(genestar$last_exam_visit,useNA = 'always')
genestar$FamilyID = genestar$FAMILY_ID
genestar$PaternalID = genestar$FATHER
genestar$MaternalID = genestar$MOTHER

genestar$study = "GeneSTAR"
genestar$unique_id <-  paste(genestar$study,genestar$SUBJECT_ID, sep = "_")
genestar$individual_id <- genestar$SUBJECT_ID
genestar$study_ancestry <- paste(genestar$study,genestar$ancestry, sep = "_")
genestar$JWsource = "dbGaP_Ex"

genestar <- genestar[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
                        'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
                        'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]

#n=3207 1113 missing t2d ??
table(duplicated(genestar$unique_id)) #0
####################### GeneSTAR ##############################################
####################### GeneSTAR ##############################################


####################### GENOA ##############################################
####################### GENOA ##############################################


## GENOA download 31OCT2017

# recode & check variable names & distributions
table(genoa$sequenced,useNA='always')

# rename column for diabetes
names(genoa)[names(genoa) == "T2D"] <- "t2d"
table(genoa$t2d,useNA='always')
#genoa = subset(genoa, subset = t2d %in% c(0,1,2))

genoa$ancestry = 'AF'
table(genoa$ancestry,useNA='always')
table(genoa$sex,useNA='always')
genoa$sex[genoa$sex == 1] = 'M'
genoa$sex[genoa$sex == 2] = 'F'
table(genoa$sex,useNA='always')
summary(genoa$last_exam_age,useNA='always') 
genoa$last_exam_age[genoa$last_exam_age == -9] = NA
summary(genoa$last_exam_age,useNA='always') 
summary(genoa$last_exam_BMI)
genoa$last_exam_bmi = genoa$last_exam_BMI
genoa$last_exam_bmi[genoa$last_exam_bmi == -9] = NA
summary(genoa$last_exam_bmi)
summary(genoa$T2D_AGE)
summary(genoa$T2D_BMI)
genoa$t2d_age = genoa$T2D_AGE
genoa$t2d_bmi = genoa$T2D_BMI
genoa$last_exam_fg = genoa$last_exam_FG
genoa$last_exam_fg[genoa$last_exam_fg == -9] = NA
summary(genoa$last_exam_fg)
genoa$last_exam_hba1c = genoa$last_exam_HbA1c
genoa$last_exam_t2d_treatment = genoa$last_exam_T2D_treatment
genoa$FamilyID = genoa$FAMILY_ID
genoa$PaternalID = genoa$FATHER
genoa$MaternalID = genoa$MOTHER

genoa$study = "GENOA"
genoa$unique_id <-  paste(genoa$study,genoa$SUBJECT_ID, sep = "_")
genoa$individual_id <- genoa$SUBJECT_ID
genoa$study_ancestry <- paste(genoa$study,genoa$ancestry, sep = "_")
genoa$JWsource = "dbGaP_Ex"

genoa <- genoa[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
                  'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
                  'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]
#n=1854 58 missing t2d
table(duplicated(genoa$unique_id)) #0

####################### GENOA ##############################################
####################### GENOA ##############################################

####################### goldn ##############################################
####################### goldn ##############################################

## goldn download
names(goldn)
goldn$sequenced = NA
# recode & check variable names & distributions
# rename column for diabetes
names(goldn)[names(goldn) == "DIAB_Status_Final"] <- "t2d"
table(goldn$t2d,useNA='always')
#goldn = subset(goldn, subset = t2d %in% c(0,1,2))


# all individuals are EU per investigators correspondence
goldn$ancestry = 'EU'
table(goldn$ancestry,useNA='always')
table(goldn$Sex,useNA='always')
goldn$sex = goldn$Sex
goldn$sex[goldn$Sex == 1] = 'M'
goldn$sex[goldn$Sex == 2] = 'F'
with(goldn,table(Sex,sex,useNA='always'))

goldn$last_exam_age = goldn$age.last_exam_age
summary(goldn$last_exam_age)
goldn$last_exam_bmi = goldn$BMI.last_exam_BMI.DIAB_BMI
summary(goldn$last_exam_bmi)
goldn$last_exam_fg = as.numeric(goldn$FastingGlucose..mg.dL.)*0.0555  # !! convert !!
summary(goldn$last_exam_fg)
goldn$t2d_age = as.numeric(goldn$DIAB_age)
summary(goldn$t2d_age)
table(goldn$last_exam_visit,useNA='always')
goldn$t2d_bmi = NA
goldn$last_exam_hba1c = NA
goldn$last_exam_t2d_treatment = NA
goldn$FamilyID = NA
goldn$PaternalID = NA
goldn$MaternalID = NA
goldn$individual_id <- goldn$subject_id

goldn$study = "GOLDN"
goldn$unique_id <-  paste(goldn$study,goldn$subject_id, sep = "_")
goldn$study_ancestry <- paste(goldn$study,goldn$ancestry, sep = "_")
goldn$JWsource = "dbGaP_Ex"

goldn <- goldn[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
                  'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
                  'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')] 
#n=1100
table(duplicated(goldn$unique_id)) #31

# remove rows that are absolute duplicates
goldn.dups <- goldn$unique_id[duplicated(goldn$unique_id)]
#unique(goldn.dups)
goldn.final <- goldn[!(goldn$unique_id %in% goldn.dups),]
for (id in unique(goldn.dups)){
  my.rows <- goldn[goldn$unique_id == id,]
  is_equal = T
  for (col.id in seq(1,NCOL(my.rows))){
    this.col <- my.rows[,col.id]
    # change nas to some string
    this.col[is.na(this.col)] <- "This is na"
    
    # see how many unique elements we ahve
    nu <- length(unique(this.col)) > 1
    
    if (nu){
      is_equal = F
      break
    }
  }
  if(is_equal){
    goldn.final <- rbind(goldn.final,my.rows[1,])
  } else {
    goldn.final <- rbind(goldn.final,my.rows)
  }
}

#n=1069
table(duplicated(goldn.final$unique_id)) #0

####################### goldn ##############################################
####################### goldn ##############################################

####################### hypergen ##############################################
####################### hypergen ##############################################

## hypergen download
names(hypergen)
hypergen$sequenced = NA
# recode & check variable names & distributions
# rename column for diabetes
names(hypergen)[names(hypergen) == "X_DBFINAL"] <- "t2d"
table(hypergen$t2d,useNA='always')
#hypergen = subset(hypergen, subset = t2d %in% c(0,1,2))

#HyperGEN LABID (M2394) and NUMID (7736) self-reported ethnicity is asian. 
#However the asian sample was seen as an outlier in the PCA plots by investigators.
hypergen$ancestry = 'AF'
table(hypergen$ancestry,useNA='always')
hypergen$sex = hypergen$SEX
hypergen$sex[hypergen$SEX == 1] = 'M'
hypergen$sex[hypergen$SEX == 2] = 'F'
with(hypergen,table(SEX,sex,useNA='always'))

hypergen$last_exam_age = hypergen$AGE
summary(hypergen$last_exam_age)
hypergen$last_exam_bmi = as.numeric(hypergen$X_BMI)
summary(hypergen$last_exam_bmi)
hypergen$last_exam_fg = as.numeric(hypergen$GLUCCLEAN)*0.0555  # !! convert !!
summary(hypergen$last_exam_fg)
hypergen$t2d_age = as.numeric(hypergen$DBDXAGE)
summary(hypergen$t2d_age)
hypergen$t2d_bmi = NA
hypergen$last_exam_hba1c = NA
hypergen$last_exam_t2d_treatment = NA
hypergen$FamilyID = NA
hypergen$PaternalID = NA
hypergen$MaternalID = NA
hypergen$individual_id <- hypergen$NUMID

hypergen$study = "HyperGEN"
hypergen$unique_id <-  paste(hypergen$study,hypergen$NUMID, sep = "_")
hypergen$study_ancestry <- paste(hypergen$study,hypergen$ancestry, sep = "_")
hypergen$JWsource = "dbGaP_Ex"

hypergen <- hypergen[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
                        'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
                        'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]
#n=1898
table(duplicated(hypergen$unique_id)) #0
####################### hypergen ##############################################
####################### hypergen ##############################################

####################### JHS ##############################################
####################### JHS ##############################################

###############################################
##JHS
# 15AUG2017 JHS has 2 MZ twins, currently co-twins are being excluded
# 22AUG2017 now all co-twins are being inlcuded

# recode & check variable names & distributions
table(jhs$sequenced,useNA='always')

# rename column for diabetes
names(jhs)[names(jhs) == "T2D"] <- "t2d"
table(jhs$t2d,useNA='always')
#jhs = subset(jhs, subset = t2d %in% c(0,1,2))


jhs$ancestry = 'AF'
table(jhs$ancestry,useNA='always')
table(jhs$sex,useNA='always')
jhs$origsex = jhs$sex
jhs$sex[jhs$sex == 1] = 'M'
jhs$sex[jhs$sex == 2] = 'F'
with(jhs,table(origsex,sex,useNA='always'))
jhs$last_exam_age[jhs$last_exam_age == 'x'] = 'NA'
jhs$last_exam_age <- as.numeric(jhs$last_exam_age)
summary(jhs$last_exam_age)
jhs$last_exam_BMI[jhs$last_exam_BMI == 'x'] = 'NA'
jhs$last_exam_bmi <- as.numeric(jhs$last_exam_BMI)
summary(jhs$last_exam_bmi)
jhs$last_exam_FG[jhs$last_exam_FG == 'x'] = 'NA'
jhs$last_exam_fg <- as.numeric(jhs$last_exam_FG)
summary(jhs$last_exam_fg)
jhs$last_exam_HbA1c[jhs$last_exam_HbA1c == 'x'] = 'NA'
jhs$last_exam_hba1c <- as.numeric(jhs$last_exam_HbA1c)
summary(jhs$last_exam_hba1c)
jhs$last_exam_T2D_treatment[jhs$last_exam_T2D_treatment == 'x'] = 'NA'
jhs$last_exam_t2d_treatment <- as.numeric(jhs$last_exam_T2D_treatment)
table(jhs$last_exam_t2d_treatment)
jhs$T2D_age[jhs$T2D_age == 'x'] = 'NA'
jhs$t2d_age <- as.numeric(jhs$T2D_age)
summary(jhs$t2d_age)
jhs$T2D_BMI[jhs$T2D_BMI == 'x'] = 'NA'
jhs$t2d_bmi <- as.numeric(jhs$T2D_BMI)
summary(jhs$t2d_bmi)
jhs$last_exam_visit = NA
jhs$FamilyID = jhs$Family_ID
jhs$PaternalID = jhs$Father_ID
jhs$MaternalID = jhs$Mother_ID
jhs$study = "JHS"
jhs$unique_id <-  paste(jhs$study,jhs$Individual_ID, sep = "_")
jhs$individual_id <- jhs$Individual_ID
jhs$study_ancestry <- paste(jhs$study,jhs$ancestry, sep = "_")
jhs$JWsource = "dbGaP_Ex"

jhs <- jhs[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
              'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
              'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]
#n=3406 10 missng t2d
table(duplicated(jhs$unique_id)) #0
####################### JHS ##############################################
####################### JHS ##############################################


####################### MESA ##############################################
####################### MESA ##############################################

# MESA 
# download 1NOV2017

# recode & check variable names & distributions
table(mesa$sequenced,useNA='always')
# rename column for diabetes
names(mesa)[names(mesa) == "T2D"] <- "t2d"
table(mesa$t2d,useNA='always')
#mesa = subset(mesa, subset = t2d %in% c(0,1,2))

table(mesa$ancestry,useNA='always')
table(mesa$Sex,useNA='always')
mesa$sex[mesa$Sex == 1] = 'M'
mesa$sex[mesa$Sex == 2] = 'F'
with(mesa,table(Sex,sex,useNA='always'))
summary(mesa$last_exam_age,useNA='always')
summary(mesa$last_exam_bmi) # ! low BMI
mesa$last_exam_fg = mesa$last_exam_FG
summary(mesa$last_exam_fg)
summary(mesa$last_exam_HBa1c)
mesa$last_exam_hba1c = mesa$last_exam_HBa1c
mesa$last_exam_t2d_treatment = mesa$last_exam_T2D_treatment
mesa$t2d_age <- as.numeric(as.character(mesa$T2D_Age))
summary(mesa$t2d_age)
mesa$t2d_bmi <- as.numeric(as.character(mesa$T2D_BMI))
summary(mesa$t2d_bmi)
table(mesa$last_exam_visit,useNA = 'always')
mesa$FamilyID = mesa$Family_ID
mesa$PaternalID = mesa$Father_ID
mesa$MaternalID = mesa$Mother_ID
mesa$individual_id <- mesa$Individual_ID

mesa$study = "MESA"
mesa$unique_id <-  paste(mesa$study,mesa$Individual_ID, sep = "_")
mesa$individual_id <- mesa$Individual_ID
mesa$study_ancestry <- paste(mesa$study,mesa$ancestry, sep = "_")


mesa$JWsource = "dbGaP_Ex"

mesa <- mesa[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
                'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
                'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]

#n=6429 9 missing t2d
table(duplicated(mesa$unique_id)) #0
####################### MESA ##############################################
####################### MESA ##############################################

####################### MESA Family  ##############################################
####################### MESA Family  ##############################################
#
# MESA Family 
## download 1NOV2017

# recode & check variable names & distributions
table(mesafam$sequenced,useNA='always')
#mesafam = subset(mesafam, sequenced == 1) #n=694 drop 8 individuals

# rename column for diabetes
names(mesafam)[names(mesafam) == "T2D"] <- "t2d"
table(mesafam$t2d,useNA='always')
#mesafam = subset(mesafam, subset = t2d %in% c(0,1,2))

mesafam$ancestry = 'AF'
table(mesafam$ancestry,useNA='always')
table(mesafam$SEX,useNA='always')
mesafam$sex[mesafam$SEX == 1] = 'M'
mesafam$sex[mesafam$SEX == 2] = 'F'
with(mesafam,table(SEX,sex,useNA='always'))
summary(mesafam$last_exam_age,useNA='always')

summary(mesafam$t2d_age,useNA='always') #! less than 25

summary(mesafam$last_exam_BMI) # ! low BMI
mesafam$last_exam_bmi = mesafam$last_exam_BMI
mesafam$last_exam_fg = mesafam$last_exam_FG
mesafam$last_exam_hba1c = NA
mesafam$last_exam_t2d_treatment = mesafam$last_exam_T2D_treatment
table(mesafam$last_exam_visit,useNA = 'always')
mesafam$FamilyID = mesafam$FAMID
mesafam$PaternalID = mesafam$FAID
mesafam$MaternalID = mesafam$MOID
mesafam$individual_id <- mesafam$INDID

mesafam$study = "MESA"
mesafam$unique_id <-  paste(mesafam$study,mesafam$INDID, sep = "_")
mesafam$individual_id <- mesafam$INDID
mesafam$study_ancestry <- paste(mesafam$study,mesafam$ancestry, sep = "_")
mesafam$JWsource = "dbGaP_Ex"


mesafam <- mesafam[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
                      'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
                      'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]
#n=1134 14 missing t2d 
table(duplicated(mesafam$unique_id)) #0
mesafam <- mesafam[!(mesafam$unique_id %in% mesa$unique_id),]
####################### MESA Family  ##############################################
####################### MESA Family  ##############################################


####################### AFib MGH ##############################################
####################### AFib MGH ##############################################
## new file downloaded from dbGaP 5AUG2017
# 7DEC2017 checked for new files - none
names(raw.MGH)

# combine pedigree with pheno
raw.MGHPed <- raw.MGHPed[,c("FAMILY_ID",  "SUBJECT_ID", "FATHER", "MOTHER")]
raw.MGH <- merge(raw.MGH,raw.MGHPed,by.x='SUBJECT_ID',by.y='SUBJECT_ID',all.x=TRUE) #n=1025 & 25 variables

# remove the old variable
rm(raw.MGHPed)

raw.MGH$sequenced = NA

# note where this file came from
raw.MGH$JWsource = "dbGaP"

# recode & check variable names & distributions
# rename column for diabetes
names(raw.MGH)[names(raw.MGH) == "diabetes"] <- "t2d"
# recode t2d status
raw.MGH$t2d[raw.MGH$t2d == 'yes'] = 2
raw.MGH$t2d[raw.MGH$t2d == 'no'] = 0
#raw.MGH = subset(raw.MGH, subset = t2d %in% c(0,1,2))
table(raw.MGH$t2d,useNA='always')
# no  yes <NA>
#  934   59    0

table(raw.MGH$race,useNA='always')
# african_american       american_indian american_indian_black                 asian                 other     pacific_islanders                 white 
# 3                     9                     1                     3                     2                     1                  1006 
# <NA> 
#   0 
table(raw.MGH$ethnicity,useNA='always')
#  no  yes <NA>
# 984   15   26

# recode ancestry
## AMR = american indian
raw.MGH$ancestry[raw.MGH$race == 'american_indian' & raw.MGH$ethnicity == 'no'] = "AMR"
raw.MGH$ancestry[raw.MGH$race == 'american_indian_black' & raw.MGH$ethnicity == 'no'] = "MIXED"
raw.MGH$ancestry[raw.MGH$race == 'white' & raw.MGH$ethnicity == 'no'] = "EU"
raw.MGH$ancestry[raw.MGH$race == 'white' & is.na(raw.MGH$ethnicity)] = "EU"
raw.MGH$ancestry[raw.MGH$race == 'white' & raw.MGH$ethnicity == 'yes'] = "HS"
raw.MGH$ancestry[raw.MGH$race == 'african_american' & raw.MGH$ethnicity == 'no'] = "AF"
raw.MGH$ancestry[raw.MGH$race == 'african_american' & is.na(raw.MGH$ethnicity)] = "AF"
raw.MGH$ancestry[raw.MGH$race == 'asian' & is.na(raw.MGH$ethnicity)] = "AS"
raw.MGH$ancestry[raw.MGH$race == 'pacific_islanders' & is.na(raw.MGH$ethnicity)] = "AS"
raw.MGH$ancestry[raw.MGH$race == 'pacific_islanders' & raw.MGH$ethnicity == 'yes'] = "MIXED"
raw.MGH$ancestry[raw.MGH$race == 'other' & is.na(raw.MGH$ethnicity)] = "OTHER"

table(raw.MGH$ancestry,useNA='always')
# AF AMR  AS  EU    HS MIXED  OTHER  <NA>
# 3   9   4   991    15     1   2   0

table(raw.MGH$sex,useNA='always')
# female   male   <NA>
#    205    788      0

# recode sex
raw.MGH$origsex = raw.MGH$sex
raw.MGH$sex[raw.MGH$sex == 'male'] = 'M'
raw.MGH$sex[raw.MGH$sex == 'female'] = 'F'
with(raw.MGH,table(origsex,sex,useNA='always'))

summary(raw.MGH$age)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 18.00   48.00   56.00   54.41   61.00   82.00

# recode age
raw.MGH$last_exam_age = raw.MGH$age

#raw.MGH = subset(raw.MGH, age >= 25) #drop 11 individuals
summary(raw.MGH$height)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 54.00   68.00   71.00   70.28   73.00   82.00       4
summary(raw.MGH$weight)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 100.0   175.0   198.0   201.6   225.0   480.0       4 

# calc bmi
raw.MGH$last_exam_bmi = ((703*raw.MGH$weight)/(raw.MGH$height*raw.MGH$height))

summary(raw.MGH$last_exam_bmi) # low BMI - exclude?

# assign missing values for things we dont have data for
raw.MGH$last_exam_fg = NA
raw.MGH$last_exam_hba1c = NA
raw.MGH$last_exam_t2d_treatment = NA
raw.MGH$last_exam_visit = NA
raw.MGH$t2d_age = NA
raw.MGH$t2d_bmi = NA
raw.MGH$FamilyID = raw.MGH$FAMILY_ID
raw.MGH$PaternalID = raw.MGH$FATHER
raw.MGH$MaternalID = raw.MGH$MOTHER
raw.MGH$study = 'MGH_AF'
raw.MGH$unique_id <-  paste(raw.MGH$study,raw.MGH$SUBJECT_ID, sep = "_")
raw.MGH$individual_id <- raw.MGH$SUBJECT_ID
raw.MGH$study_ancestry <- paste(raw.MGH$study,raw.MGH$ancestry, sep = "_")

raw.MGH <- raw.MGH[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
                      'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
                      'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]


#n=1025 32 missing t2d 
table(duplicated(raw.MGH$unique_id)) #0
####################### AFib MGH ##############################################
####################### AFib MGH ##############################################

####################### AFib VU ##############################################
####################### AFib VU ##############################################

#
## AFib VU
# downloaded new dataset 7DEC2017
names(raw.VUdaw)
raw.VUdawPed <- raw.VUdawPed[,c("FAMILY_ID",  "SUBJECT_ID", "FATHER", "MOTHER")]
raw.VUdaw <- merge(raw.VUdaw,raw.VUdawPed,by.x='SUBJECT_ID',by.y='SUBJECT_ID',all.x=TRUE) #n=1134
rm(raw.VUdawPed)
raw.VUdaw$JWsource = "dbGaP"
raw.VUdaw$sequenced = NA

# recode & check variable names & distributions
# rename column for diabetes
names(raw.VUdaw)[names(raw.VUdaw) == "diabetes"] <- "t2d"
# recode t2d status
raw.VUdaw$t2d[raw.VUdaw$t2d == 'yes'] = 2
raw.VUdaw$t2d[raw.VUdaw$t2d == 'no'] = 0
table(raw.VUdaw$t2d,useNA='always')
raw.VUdaw = subset(raw.VUdaw, subset = t2d %in% c(0,1,2))

table(raw.VUdaw$race,useNA='always')
table(raw.VUdaw$ethnicity,useNA='always')
raw.VUdaw$ancestry[raw.VUdaw$race == 'native_american' & raw.VUdaw$ethnicity == 'no'] = "AMR"
raw.VUdaw$ancestry[raw.VUdaw$race == 'black' & raw.VUdaw$ethnicity == 'no'] = "AF"
raw.VUdaw$ancestry[raw.VUdaw$race == 'white' & raw.VUdaw$ethnicity == 'no'] = "EU"
raw.VUdaw$ancestry[raw.VUdaw$race == 'asian' & raw.VUdaw$ethnicity == 'no'] = "AS"
raw.VUdaw$ancestry[raw.VUdaw$race == 'white' & raw.VUdaw$ethnicity == 'yes'] = "HS"
table(raw.VUdaw$ancestry,useNA='always')
table(raw.VUdaw$sex,useNA='always')
raw.VUdaw$origsex = raw.VUdaw$sex
raw.VUdaw$sex[raw.VUdaw$sex == 'male'] = 'M'
raw.VUdaw$sex[raw.VUdaw$sex == 'female'] = 'F'
with(raw.VUdaw,table(origsex,sex,useNA='always'))
summary(raw.VUdaw$age)  ## Contains some < 18 yrs
raw.VUdaw$last_exam_age = raw.VUdaw$age
#raw.VUdaw = subset(raw.VUdaw, age >= 25) #n=1085, drop 24 individuals
summary(raw.VUdaw$height)
summary(raw.VUdaw$weight)
raw.VUdaw$last_exam_bmi = ((703*raw.VUdaw$weight)/(raw.VUdaw$height*raw.VUdaw$height))
summary(raw.VUdaw$last_exam_bmi) # low BMI - exclude?
raw.VUdaw$last_exam_fg = NA
raw.VUdaw$last_exam_hba1c = NA
raw.VUdaw$last_exam_t2d_treatment = NA
raw.VUdaw$last_exam_visit = NA
raw.VUdaw$t2d_age = NA
raw.VUdaw$t2d_bmi = NA
raw.VUdaw$FamilyID = raw.VUdaw$FAMILY_ID
raw.VUdaw$PaternalID = raw.VUdaw$FATHER
raw.VUdaw$MaternalID = raw.VUdaw$MOTHER
raw.VUdaw$study = "VU_AF"
raw.VUdaw$unique_id <-  paste(raw.VUdaw$study,raw.VUdaw$SUBJECT_ID, sep = "_")
raw.VUdaw$individual_id <- raw.VUdaw$SUBJECT_ID
raw.VUdaw$study_ancestry <- paste(raw.VUdaw$study,raw.VUdaw$ancestry, sep = "_")

raw.VUdaw <- raw.VUdaw[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
                          'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
                          'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]
#n=1134
table(duplicated(raw.VUdaw$unique_id)) #0
####################### AFib VU ##############################################
####################### AFib VU ##############################################
####################### HVH part of AFib and VTE projects ##############################################
####################### HVH part of AFib and VTE projects ##############################################


## HVH part of AFib and VTE projects
# HVH updated 18AUG2017

names(raw.HVH)

# recode & check variable names & distributions
# rename column for diabetes
names(raw.HVH)[names(raw.HVH) == "T2D"] <- "t2d"
# recode t2d status
raw.HVH$t2d[raw.HVH$t2d == 'yes'] = 2
raw.HVH$t2d[raw.HVH$t2d == 'no'] = 0
table(raw.HVH$t2d,useNA='always')
#raw.HVH = subset(raw.HVH, subset = t2d %in% c(0,1,2))

table(raw.HVH$RaceCalc,useNA='always')
raw.HVH$ancestry = raw.HVH$RaceCalc
raw.HVH$ancestry[raw.HVH$ancestry == 0] = 'AF'
raw.HVH$ancestry[raw.HVH$ancestry == 1] = 'EU'
raw.HVH$ancestry[raw.HVH$ancestry == 2] = 'OTHER'
table(raw.HVH$ancestry,useNA='always')
table(raw.HVH$sex,useNA='always')
raw.HVH$origsex = raw.HVH$sex
raw.HVH$sex[raw.HVH$sex == 1] = 'M'
raw.HVH$sex[raw.HVH$sex == 2] = 'F'
with(raw.HVH,table(origsex,sex,useNA='always'))
summary(raw.HVH$last_exam_age)
raw.HVH$last_exam_bmi = raw.HVH$last_exam_BMI
summary(raw.HVH$last_exam_bmi)
raw.HVH$last_exam_fg = NA
raw.HVH$last_exam_hba1c = NA
raw.HVH$last_exam_t2d_treatment = NA
raw.HVH$last_exam_visit = NA
raw.HVH$t2d_age = raw.HVH$T2D_age
summary(raw.HVH$T2D_age)
raw.HVH$t2d_bmi = NA
raw.HVH$study = "HVH"
raw.HVH$unique_id <-  paste(raw.HVH$study,raw.HVH$SampleID, sep = "_")
raw.HVH$individual_id <- raw.HVH$SampleID
raw.HVH$study_ancestry <- paste(raw.HVH$study,raw.HVH$ancestry, sep = "_")
raw.HVH$JWsource = "dbGaP_Ex"

raw.HVH <- raw.HVH[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
                      'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
                      'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]
#n=708 3 missing t2d 
table(duplicated(raw.HVH$unique_id)) #0
####################### HVH part of AFib and VTE projects ##############################################
####################### HVH part of AFib and VTE projects ##############################################

####################### SAFSCVD ##############################################
####################### SAFSCVD ##############################################


# ###### ADDED THIS RECODE OF STUDIES  13FEB2018
# ## BEGIN

# # #
# # SAFS CVD
# ## 9 SEP 2017 DOWNLOADED new dataset
# # 7DEC2017 checked for new files - none
# 

# 1: Get rows of map file that correspond to NWD ids from the id file
# safs.map <- map[map$sample.id %in% safs.ids$NWD.ID,] # 1509

#dim(safs.map)

#safs.map <- safs.map[which(!is.na(safs.map$submitted_subject_id)),]
#dim(safs.map)

## 2: Get the unique deidentified subject IDs from the id file
##uniq.deident <- unique(safs.ids$Deidentified.Subject) # 2624

## 3: Check that we have the only deidentified subject ids in the map file
#length(unique(safs.map$submitted_subject_id)) # 1474, so we have 35 duplicated NWD ids

# 4: Get the NWD ids for our unique deidentified subject ids that are in the map file
#safs.nwd <- safs.map[safs.map$submitted_subject_id %in% uniq.deident,] # 1502

# 5: Subset the phenotype file by these nwd ids
# restrict to NWD_ID's that are mapped to a subject_id in DCC map file

###### TDM edits Mar142018
# safs = 2457 rows
# safs.ids = 2713 rows

# -1: Subset map to remove any rows that do not have a submitted subject ID (these have most likely been NA'ed out by DCC during their QC and should not be included in analysis)
# safs.map <- map[!is.na(map$submitted_subject_id),]
# 
# # 1: Subset deidentified subject ids from SAFS by those deidentified subject ids that are in the map file (DCC)
# safs.ids <- safs.ids[safs.ids$Deidentified.Subject %in% safs.map$submitted_subject_id,] # 1531 rows
# 
# # 2: Merge id file with ped file by deidentified subject id, keeping only those deidentified subject ids that are in the map file
# safs <- merge(safs.ids, safs, by.x = "NWD.ID", by.y = "Individual_ID", all.x=T) #1531 rows
# 
# # 2.1: Check that we that we do not have duplicated NWD ids (this would signal an NWD id mapping to multiple deidentified subject ids, which would be bad.)
# safs[duplicated(safs$NWD.ID),]
# 
# # [1] NWD.ID                Sample.ID             Deidentified.Subject  Mexam.Freeze.Number   Mexam.Revision.Number Family_ID             Mother_ID             Father_ID            
# # [9] T2D                   Sex                   last_exam_age         last_exam_BMI         last_exam_FG          last_exam_HbA1c       last_exam_treatment   last_exam_visit      
# # [17] T2D_age               T2D_BMI               Sequenced            
# # <0 rows> (or 0-length row.names)
# ########
# 
# #head(safs.ped)
# CHANGED safs.ped to safs starting here

names(safs)
names(safs)[names(safs) == "Sequenced"] <- "sequenced"
safs$JWsource = "dbGaP_Ex"

# # recode & check variable names & distributions
table(safs$sequenced,useNA='always') 

# rename column for diabetes
names(safs)[names(safs) == "T2D"] <- "t2d"
table(safs$t2d,useNA='always')
#safs = subset(safs, subset = t2d %in% c(0,1,2))

safs$ancestry = 'HS'
table(safs$ancestry,useNA='always')
table(safs$Sex,useNA='always')
safs$sex[safs$Sex == 1] = 'M'
safs$sex[safs$Sex == 2] = 'F'
with(safs,table(Sex,sex,useNA='always'))
safs$last_exam_age = as.numeric(safs$last_exam_age)
summary(safs$last_exam_age,useNA='always')
#safs = subset(safs, last_exam_age >= 25) 
safs$last_exam_bmi = as.numeric(safs$last_exam_BMI)
summary(safs$last_exam_bmi) # ! low BMI
safs$last_exam_fg = as.numeric(safs$last_exam_FG)
summary(safs$last_exam_fg)
safs$last_exam_hba1c = as.numeric(safs$last_exam_HbA1c)
summary(safs$last_exam_hba1c)
table(safs$last_exam_treatment,useNA = 'always')
safs$last_exam_t2d_treatment = safs$last_exam_treatment
safs$t2d_age = as.numeric(safs$T2D_age)
summary(safs$t2d_age)
safs$t2d_bmi = as.numeric(safs$T2D_BMI)
summary(safs$t2d_bmi)
table(safs$t2d_bmi,useNA = 'always')
table(safs$last_exam_visit,useNA = 'always')
safs$FamilyID = safs$Family.ID
safs$PaternalID = safs$Father_ID
safs$MaternalID = safs$Mother_ID

safs$study = "SAFS"
safs$unique_id <-  paste(safs$study,safs$Subject.ID..deID., sep = "_")
safs$individual_id <- safs$Subject.ID..deID.
safs$study_ancestry <- paste(safs$study,safs$ancestry, sep = "_")

safs <- safs[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
                'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
                'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]

# n=4734
table(duplicated(safs$unique_id)) #5

# remove rows that are absolute duplicates
safs.dups <- safs$unique_id[duplicated(safs$unique_id)]
safs.final <- safs[!(safs$unique_id %in% safs.dups),]
for (id in unique(safs.dups)){
  my.rows <- safs[safs$unique_id == id,]
  is_equal = T
  for (col.id in seq(1,NCOL(my.rows))){
    this.col <- my.rows[,col.id]
    # change nas to some string
    this.col[is.na(this.col)] <- "This is na"
    
    # see how many unique elements we ahve
    nu <- length(unique(this.col)) > 1
    
    if (nu){
      is_equal = F
      break
    }
  }
  if(is_equal){
    safs.final <- rbind(safs.final,my.rows[1,])
  } else {
    safs.final <- rbind(safs.final,my.rows)
  }
}

#n=4730; 1 individual is not eliminated!!
table(duplicated(safs.final$unique_id)) #1

####################### SAS ##############################################
####################### SAS ##############################################

## SAS

#T2D NA=0 & -9=497

#sas2 <- sas[is.na(sas$T2D),]

#n=8 with WGS 7 no phenotype
names(sas)
# recode & check variable names & distributions
table(sas$sequenced,useNA='always')
#   0    1 <NA> 
#  3090  380    0 

# rename column for diabetes
names(sas)[names(sas) == "T2D"] <- "t2d"
table(sas$t2d,useNA='always')
#sas = subset(sas, subset = t2d %in% c(0,1,2))

sas$ancestry = 'AS'
table(sas$ancestry,useNA='always')
table(sas$sex,useNA='always')
sas$origsex = sas$sex
sas$sex[sas$sex == 1] = 'M'
sas$sex[sas$sex == 2] = 'F'
with(sas,table(origsex,sex,useNA='always'))
sas$last_exam_age = as.numeric(sas$last_exam_age)
summary(sas$last_exam_age)
summary(sas$last_exam_BMI)
sas$last_exam_bmi = sas$last_exam_BMI
sas$last_exam_fg = sas$last_exam_FG
summary(sas$last_exam_fg)
sas$last_exam_hba1c = NA
sas$last_exam_t2d_treatment = sas$last_exam_T2D_treatment
sas$t2d_age = sas$T2D_age
summary(sas$T2D_age)
summary(sas$T2D_BMI)
sas$t2d_bmi = as.numeric(sas$T2D_BMI)
summary(sas$t2d_bmi)
table(sas$last_exam_visit,useNA = 'always')
sas$FamilyID = sas$Family_ID
sas$PaternalID = sas$Father_ID
sas$MaternalID = sas$Mother_ID
sas$study = "SAS"
sas$unique_id <-  paste(sas$study,sas$SUBJECT_ID, sep = "_")
sas$individual_id <- sas$SUBJECT_ID
sas$study_ancestry <- paste(sas$study,sas$ancestry, sep = "_")
sas$JWsource = "dbGaP_Ex"


sas <- sas[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
              'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
              'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]

#n=3470
table(duplicated(sas$unique_id)) #0
####################### SAS ##############################################
####################### SAS ##############################################


####################### WHI ##############################################
####################### WHI ##############################################

## WHI
names(whi)
# recode & check variable names & distributions
table(whi$t2d,useNA='always') 
#whi = subset(whi, subset = t2d %in% c(0,1,2))

table(whi$ancestry,useNA='always')
table(whi$sex,useNA='always')
whi$origsex = whi$sex
whi$sex[whi$sex == 2] = 'F'
with(whi,table(origsex,sex,useNA='always'))
summary(whi$last_exam_age)
summary(whi$t2d_age)
summary(whi$last_exam_bmi) #!low BMI
whi$study = "WHI"
whi$unique_id <-  paste(whi$study,whi$individual_id, sep = "_")
whi$individual_id <- whi$individual_id
whi$study_ancestry <- paste(whi$study,whi$ancestry, sep = "_")
whi$JWsource = "dbGaP_Ex"

whi <- whi[,c('unique_id','individual_id','FamilyID','MaternalID','PaternalID',
              'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg', 'sequenced',
              'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]
#n=10938 2245 missing t2d ??? 
table(duplicated(whi$unique_id)) #0
####################### WHI ##############################################
####################### WHI ##############################################

# ##END OF RECODE TO STANDARDIZE VARIABLES ACROSS STUDIES FOR POOLING

########### CREATE POOLED DATASET ###########################
# 
# studies to add in future thrv, , 
pooled <- rbind(afccaf,afp,afvub,amish,aric,cfs.final,chs,copd,dhs,fhs,genestar,genoa,
                gensalt,goldn.final,hypergen,jhs,mesa,mesafam,raw.HVH,raw.MGH,raw.VUdaw,safs.final,sas,whi)
#n=80377

# # remove rows that are absolute duplicates
# pooled.dups <- pooled$unique_id[duplicated(pooled$unique_id)]
# unique(pooled.dups)
# pooled.final <- pooled[!(pooled$unique_id %in% pooled.dups),]
# for (id in unique(pooled.dups)){
#   my.rows <- pooled[pooled$unique_id == id,]
#   is_equal = T
#   for (col.id in seq(1,NCOL(my.rows))){
#     this.col <- my.rows[,col.id]
#     # change nas to some string
#     this.col[is.na(this.col)] <- "This is na"
#     
#     # see how many unique elements we ahve
#     nu <- length(unique(this.col)) > 1
#     
#     if (nu){
#       is_equal = F
#       break
#     }
#   }
#   if(is_equal){
#     pooled.final <- rbind(pooled.final,my.rows[1,])
#   } else {
#     pooled.final <- rbind(pooled.final,my.rows)
#   }
# }
# #n=80377 none removed

table(pooled$t2d,useNA='always') 
pooled$t2d[pooled$t2d == -9] = NA
pooled$t2d[pooled$t2d == 9] = NA
pooled$t2d[pooled$t2d == 'NA'] = NA
table(pooled$t2d,useNA='always') #n=7423 missing t2d
table(pooled$sex,useNA='always') #n=0 missing sex
table(pooled$ancestry,useNA='always') #n=0 NA
table(pooled$sequenced,useNA='always') #n=16563 NA
pooled$unique_id2 = pooled$unique_id

### JW added this bit to check for duplicates at the phenotype level 20MAR2018
table(duplicated(pooled$unique_id))
# FALSE  TRUE 
# 80200   177
table(is.na(pooled$unique_id)) #0
table(duplicated(pooled$individual_id))
# FALSE  TRUE 
# 70601  9776
table(is.na(pooled$individual_id)) #0

pooled.dup <- pooled[duplicated(pooled$unique_id),]
table(duplicated(pooled.dup$unique_id))
table(duplicated(pooled.dup$individual_id))
table(pooled.dup$study_ancestry, useNA = 'always')
# CFS_AF  CFS_EU MESA_AF SAFS_HS    <NA> 
#   2       3     171       1       0 
# TOTAL N=177

########################################

## SET-UP MAP FILE FOR MERGING WITH POOLED PHENOTYPE DATA

names(map) #n=54499
#  [1] "sample.id"            "unique_subject_key"   "submitted_subject_id"
#  [4] "consent"              "sex"                  "sexchr.kary"
#  [7] "topmed_phs"           "study"                "topmed_project"
# [10] "CENTER"               "geno.cntl"            "TRIO.dups"
# [13] "MZtwinID"             "keep"                 "unique.geno"
# [16] "unique.subj"


# per Josee Dupuis, keep trios DEC 2017 !! 


table(map$study,useNA='always')
#       Amish        ARIC        BAGS        CCAF         CFS         CHS
#        1030        3616         968         329         923          70
#    COPDGene         CRA         DHS      EOCOPD         FHS      GALAII
#        8742        1043         339          66        3758         914
#    GeneSTAR       GENOA     GenSalt       GOLDN         HVH    HyperGEN
#        1639        1144        1695         904          66        1777
#         JHS    Mayo_VTE        MESA      MGH_AF    Partners        SAFS
#        3136        1251        4880         918         111        1509
#        SAGE Sarcoidosis         SAS       VAFAR       VU_AF        WGHS
#         452         608        1208         157        1024          99
#         WHI   NA
#       10047   76
#

## checked freeze5b_duplicates_2018-01-10 for studies that overlap with the studies we use for T2D & traits; 
## excluding the studies not used and dont have overlapping studies
map <- subset(map, study %in% c('Amish','ARIC','CCAF','CFS','CHS','COPDGene','DHS','FHS','GeneSTAR','GENOA','GenSalt','GOLDN','HVH','HyperGEN',
                               'Mayo_VTE','JHS','MESA','MGH_AF','Partners','SAFS','SAS','VAFAR','VU_AF','WGHS','WHI')) #n=50372, drop 4127

## Create variables for final pooled file
map$topmedid = map$sample.id

# CHECK ####
table(is.na(map$unique_subject_key)) 
# FALSE  TRUE 
# 50100   272
table(is.na(map$sample.id)) #0
table(is.na(map$submitted_subject_id)) 
# FALSE  TRUE 
# 50100   272

map <- map[which(!is.na(map$unique_subject_key) & !is.na(map$sample.id)),] #n=50100 drop 272

table(is.na(map$unique_subject_key))
table(is.na(map$sample.id))

table(duplicated(map$unique_subject_key))
# FALSE  TRUE 
# 50035   65
table(duplicated(map$submitted_subject_id))
# FALSE  TRUE 
# 46063  4037 
table(duplicated(map$sample.id))
# FALSE 
# 50100 

map.dup <- map[duplicated(map$unique_subject_key),]
table(is.na(map.dup$unique_subject_key)) #0
table(duplicated(map.dup$submitted_subject_id))
# FALSE  TRUE 
# 61       4
table(map.dup$study,useNA = 'always')
# COPDGene      FHS  GenSalt      HVH     SAFS      SAS     <NA> 
#   10       10        4        2       29       10        0 
table(map.dup$study,map.dup$unique_subject_key,useNA = 'always')
map.dup2 <- map[duplicated(map$submitted_subject_id),]
table(duplicated(map.dup2$unique_subject_key))
# FALSE  TRUE 
# 4031   6


#########################################################
# merge this map data with Pooled data
fulldata <- merge(map,pooled,by.x=c('unique_subject_key'),by.y=c('unique_id'),all.x=TRUE) #n=50245; add 145 >!>!>

table(is.na(fulldata$unique_subject_key)) #0
table(is.na(fulldata$sample.id)) #0
table(is.na(fulldata$unique_id2))
# FALSE  TRUE 
# 48552  1693 # delete these individuals??
table(is.na(fulldata$individual_id))
# FALSE  TRUE 
# 48552  1693 

## DELETE INDIVIDUALS WITH ABSOLUTELY NO PHENOTYPE DATA??
#fulldata <- fulldata[!is.na(fulldata$unique_id2),] #n=48552
# fulldata.na <- fulldata[is.na(fulldata$unique_id2),]
# table(fulldata.na$study,useNA = 'always')
# table(fulldata.na$t2d,useNA = 'always')
#table(fulldata.na$JWsource,useNA = 'always')
#summary(fulldata.na,useNA = 'always')

### JW added this bit to check for duplicates at the merge level 20MAR2018
table(duplicated(fulldata$unique_subject_key))
# FALSE  TRUE 
# 50035   210
table(duplicated(fulldata$sample.id))
# FALSE  TRUE 
# 50100   145 ## we are creating duplicated NWD ids...
table(duplicated(fulldata$individual_id))
# FALSE  TRUE 
# 44407  5838
fulldata.dup <- fulldata[duplicated(fulldata$unique_subject_key),]
table(duplicated(fulldata.dup$sample.id))
# FALSE  
# 210     
table(duplicated(fulldata.dup$individual_id))
# FALSE  TRUE 
# 206     4 
table(fulldata.dup$study, useNA = 'always')

# COPDGene      FHS  GenSalt      HVH     MESA     SAFS      SAS     <NA> 
#   10       10        4        2      145       29       10        0 
# WTF?
########################################


# check for sex discordance
print(table(fulldata$sex.x,fulldata$sex.y,useNA='always')) #n=1 discordant for sex
wtfsex <- fulldata[which(fulldata$sex.x == 'F' & fulldata$sex.y == 'M'),]
table(wtfsex$study) #SAFS
# remove discordant sex here
fulldata <- subset(fulldata, sex.x == sex.y | is.na(sex.x) | is.na(sex.y)) #n=50244
table(fulldata$t2d) # no missing
table(fulldata$ancestry,useNA = 'always') #n=1693

#table(MAP=fulldata$sex.x,POOLED=fulldata$sex.y,T2D=fulldata$t2d,useNA = "always")

# choose a sex variable
names(fulldata)[names(fulldata) == "sex.x"] <- "sex"
names(fulldata)

fulldata <- fulldata[,c('unique_subject_key',"sample.id","submitted_subject_id","consent",
                        "sex","sexchr.kary","topmed_phs","study","topmed_project","CENTER","geno.cntl",
                        "TRIO.dups","MZtwinID","keep","unique.geno","unique.subj","topmedid",
                        'individual_id','FamilyID','MaternalID','PaternalID','t2d', 'sequenced',
                        'last_exam_age','last_exam_bmi','last_exam_fg','last_exam_hba1c',
                        'last_exam_t2d_treatment','t2d_age','t2d_bmi','JWsource','ancestry', 'study_ancestry')]


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

#collapse AFib studies
table(fulldata$study_ancestry,useNA='always')
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "ARIC_EU"] = "Afib_EU"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "CCAF_EU"] = "Afib_EU"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "HVH_EU"] = "Afib_EU"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "MGH_AF_EU"] = "Afib_EU"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "Partners_EU"] = "Afib_EU"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "VAFAR_EU"] = "Afib_EU"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "VU_AF_EU"] = "Afib_EU"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "Partners_AF"] = "Afib_AF"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "VU_AF_AF"] = "Afib_AF"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "CCAF_HS"] = "Afib_HS"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "MGH_AF_HS"] = "Afib_HS"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "VU_AF_HS"] = "Afib_HS"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "HVH_OTHER"] = "Afib_Other"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "Partners_OTHER"] = "Afib_Other"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "MGH_AF_AMR"] = "Afib_AMR"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "VU_AF_AMR"] = "Afib_AMR"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "MGH_AF_MIXED"] = "Afib_Mixed"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "Partners_AS"] = "Afib_AS"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "VU_AF_AS"] = "Afib_AS"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "Partners_NA"] = "Afib_NA"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "Amish_EU"] = "Amish_EU"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "CFS_AF"] = "CFS_AF"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "CFS_EU"] = "CFS_EU"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "COPDGene_AF"] = "COPDGene_AF"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "COPDGene_EU"] = "COPDGene_EU"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "FHS_EU"] = "FHS_EU"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "JHS_AF"] = "JHS_AF"
fulldata$study_ancestry_AFib[fulldata$study_ancestry == "SAS_AS"] = "SAS_AS"
with(fulldata,table(study_ancestry,study_ancestry_AFib,useNA='always'))
table(fulldata$t2d,useNA = 'always') #n=4754

# write.csv(fulldata,row.names=F,quote=F,file=paste(out.pref,"/",'FULLDATA_Test_T2D_15MAR2018.csv',sep=""))
write.csv(fulldata,row.names=F,quote=F,file=paste(f.dir,"/",out.pref,'.csv',sep=""))

# save(pooled, file = paste(f.dir,"Pooled_MIXED_WesselJ_27AUG2017_T2D.RData",sep="/"))


# ### AWAITING GENOTYPE DATA
# ##
# # AFib Australia
# ## new file downloaded from dbGaP 30OCT2017
# Aus =read.table(paste(f.dir,'Fatkin_Australia_dbGaP_SubjectPhenotypesDS_v1.txt',sep="/"),
#                 header=T,sep='\t',as.is=T,fill = TRUE) #n=120 & 17 variables
# Aus$JWsource = "dbGaP"
# mapAus <- map[map$study == "",] #n=
# Aus2 <- merge(mapAus,Aus,by.x='individual_id',by.y='SUBJECT_ID',all.x=TRUE) #n=789
# rm(mapAus)

# # recode & check variable names & distributions
# table(Aus$diabetes,useNA='always')
# Aus$t2d = ifelse(is.na(Aus$diabetes),-9,Aus$diabetes)
# Aus$t2d[Aus$t2d == 'yes'] = 2
# Aus$t2d[Aus$t2d == 'no'] = 0
# with(Aus,table(diabetes,t2d,useNA='always')) #n=5 unknown T2D
# table(Aus$race,useNA='always')
# table(Aus$ethnicity,useNA='always')
# Aus$ancestry[Aus$race == 'white' & Aus$ethnicity == 'no'] = "EU"
# table(Aus$ancestry,useNA='always')
# table(Aus$sex,useNA='always')
# Aus$origsex = Aus$sex
# Aus$sex[Aus$sex == 'male'] = 'M'
# Aus$sex[Aus$sex == 'female'] = 'F'
# with(Aus,table(origsex,sex,useNA='always'))
# summary(Aus$age)
# Aus = subset(Aus, age >= 25) #drop 1 individuals
# Aus$last_exam_age = Aus$age
# summary(Aus$height)
# summary(Aus$weight)
# Aus$last_exam_bmi = ((703*Aus$weight)/(Aus$height*Aus$height))
# summary(Aus$last_exam_bmi) 
# Aus$last_exam_fg = NA
# Aus$last_exam_hba1c = NA
# Aus$last_exam_t2d_treatment = NA
# Aus$last_exam_visit = NA
# Aus$t2d_age = NA
# Aus$t2d_bmi = NA
# Aus$FamilyID = Aus$individual_id
# Aus$PaternalID = 0
# Aus$MaternalID = 0
# Aus$study_ancestry <- paste(Aus$study,Aus$ancestry, sep = "_")

# Aus <- Aus[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
#               'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
#               'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
#               'study_topmedid','study_ancestry','JWsource','ancestry',
#               "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
#               "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
#               "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
# #n=65
# write.csv(Aus,row.names=F,quote=F,file=paste(f.dir,'Harmonized.AFib.Aus.XXXX2017.T2D.csv',sep="/"))


# ##
# # AFib UMass
# ## new file downloaded from dbGaP 30OCT2017
# UMass =read.table(paste(f.dir,'UMASS_dbGaP_SubjectPhenotypesDS_v1.txt',sep="/"),
#                   header=T,sep='\t',as.is=T,fill = TRUE) #n=65 & 17 variables
# UMass$JWsource = "dbGaP"
# mapUMass <- map[map$study == "",] #n=
# UMass2 <- merge(mapUMass,UMass,by.x='individual_id',by.y='SUBJECT_ID',all.x=TRUE) #n=789
# rm(mapUMass)

# # recode & check variable names & distributions
# table(UMass$diabetes,useNA='always')
# UMass$t2d = ifelse(is.na(UMass$diabetes),-9,UMass$diabetes)
# UMass$t2d[UMass$t2d == 'yes'] = 2
# UMass$t2d[UMass$t2d == 'no'] = 0
# with(UMass,table(diabetes,t2d,useNA='always')) #n=5 unknown T2D
# table(UMass$race,useNA='always')
# table(UMass$ethnicity,useNA='always')
# UMass$ancestry[UMass$race == 'white' & UMass$ethnicity == 'no'] = "EU"
# table(UMass$ancestry,useNA='always')
# table(UMass$sex,useNA='always')
# UMass$origsex = UMass$sex
# UMass$sex[UMass$sex == 'male'] = 'M'
# UMass$sex[UMass$sex == 'female'] = 'F'
# with(UMass,table(origsex,sex,useNA='always'))
# summary(UMass$age)
# UMass$last_exam_age = UMass$age
# summary(UMass$height)
# summary(UMass$weight)
# UMass$last_exam_bmi = ((703*UMass$weight)/(UMass$height*UMass$height))
# summary(UMass$last_exam_bmi) 
# UMass$last_exam_fg = NA
# UMass$last_exam_hba1c = NA
# UMass$last_exam_t2d_treatment = NA
# UMass$last_exam_visit = NA
# UMass$t2d_age = NA
# UMass$t2d_bmi = NA
# UMass$FamilyID = UMass$individual_id
# UMass$PaternalID = 0
# UMass$MaternalID = 0
# UMass$study_ancestry <- paste(UMass$study,UMass$ancestry, sep = "_")


# UMass <- UMass[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
#                   'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
#                   'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
#                   'study_topmedid','study_ancestry','JWsource','ancestry',
#                   "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
#                   "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
#                   "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
# #n=65
# write.csv(UMass,row.names=F,quote=F,file=paste(f.dir,'Harmonized.AFib.UMass.XXXX2017.T2D.csv',sep="/"))


# ## THRV download
# thrv = read.table(paste(f.dir,'THRV_Chinese_KS_09142017_T2D.ped',sep="/"), header=T,sep='\t',as.is=T) #n=2346
# thrv$JWsource = "dbGaP_Ex"
# # mapthrv <- map[map$study == "THRV",] #n=1144
# # thrv <- merge(mapthrv,thrv,by.x='individual_id',by.y='SUBJECT_ID',all.x=TRUE) #n=1144
# thrv <- merge(thrv,map,by.x='subject_id',by.y='individual_id',all.x=TRUE) #n=1144
# rm(mapthrv)

# # recode & check variable names & distributions
# table(thrv$sequenced,useNA='always')
# table(thrv$T2D,useNA='always')
# thrv$t2d = ifelse(is.na(thrv$T2D),-9,thrv$T2D)
# with(thrv,table(T2D,t2d,useNA='always'))
# thrv$ancestry = 'AS'
# table(thrv$ancestry,useNA='always')
# table(thrv$sex,useNA='always')
# thrv$sex[thrv$sex == 1] = 'M'
# thrv$sex[thrv$sex == 2] = 'F'
# table(thrv$sex,useNA='always')
# summary(thrv$last_exam_age,useNA='always')
# summary(thrv$last_exam_BMI)
# summary(thrv$T2D_age)
# summary(thrv$T2D_BMI)
# summary(thrv$last_exam_HbA1c)
# thrv$t2d_age = thrv$T2D_AGE
# thrv$t2d_bmi = thrv$T2D_BMI
# thrv$last_exam_bmi = thrv$last_exam_BMI
# thrv$last_exam_fg = thrv$last_exam_FG
# thrv$last_exam_hba1c = thrv$last_exam_HbA1c
# thrv$last_exam_t2d_treatment = thrv$last_exam_T2D_treatment
# thrv$FamilyID = thrv$family_ID
# thrv$PaternalID = thrv$father_id
# thrv$MaternalID = thrv$mother_id
# thrv$study_ancestry <- paste(thrv$study,thrv$ancestry, sep = "_")

# thrv <- thrv[,c('topmedid','individual_id','FamilyID','MaternalID','PaternalID',
#                 'sex','t2d','last_exam_age','last_exam_bmi','last_exam_fg',
#                 'last_exam_hba1c','last_exam_t2d_treatment','t2d_age','t2d_bmi','study',
#                 'study_topmedid','study_ancestry','JWsource','ancestry',
#                 "sample.id", "unique_subject_key", "submitted_subject_id", "consent", 
#                 "sexchr.kary",  "topmed_phs", "topmed_project", "CENTER", 
#                 "geno.cntl",  "TRIO.dups", "MZtwinID",  "keep", "unique.geno",  "unique.subj")]
# #n=
