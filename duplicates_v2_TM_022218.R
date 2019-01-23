args = commandArgs(trailingOnly=TRUE)
f.dir <- args[1]
out.pref <- args[2]

#### testing inputs ####
# f.dir <- "/N/dc2/scratch/wesselj/OriginalFiles"
# id.col <- "topmedid"
# ped.file <- "FULLDATA_Test_T2D_13MAR2018.csv"
# trait <- "t2d"
# out.pref <- "/N/dc2/scratch/wesselj/OriginalFiles/"
### testing inputs ####


if(length(args) < 1) {
 stop ("Not enough input args")
}


# if (length(args) == 1){
  ## load duplicates 
  dups <- read.table(paste(f.dir,"freeze5b_duplicates_2018-01-10.txt",sep="/"),header=T)
  dim(dups) #[1] 1528    5
  
  dups <- subset(dups, study1 %in% c('Amish','ARIC','CCAF','CFS','CHS','COPDGene','DHS','FHS','GeneSTAR','GENOA','GenSalt','GOLDN','HVH','HyperGEN',
                                  'Mayo_VTE','JHS','MESA','MGH_AF','Partners','SAFS','SAS','VAFAR','VU_AF','WGHS','WHI') &
                   study2 %in% c('Amish','ARIC','CCAF','CFS','CHS','COPDGene','DHS','FHS','GeneSTAR','GENOA','GenSalt','GOLDN','HVH','HyperGEN',
                                 'Mayo_VTE','JHS','MESA','MGH_AF','Partners','SAFS','SAS','VAFAR','VU_AF','WGHS','WHI')) #n=50372, drop 4127
  dim(dups) #  [1] 376   5

  ## load map file
  map <- read.table(paste(f.dir,"freeze5b_sample_annot_2018-07-13.txt",sep="/"),header=T)
  map <- subset(map, study %in% c('Amish','ARIC','CCAF','CFS','CHS','COPDGene','DHS','FHS','GeneSTAR','GENOA','GenSalt','GOLDN','HVH','HyperGEN',
                                  'Mayo_VTE','JHS','MESA','MGH_AF','Partners','SAFS','SAS','VAFAR','VU_AF','WGHS','WHI')) #n=50372, drop 4127
  
  # linker <- read.table("freeze5b_sample_annot_2017-12-01.txt",header=T)
  
  ## make some table?
  t <- table(map$study,map$CENTER)
  t <- prop.table(t,1)
  
  ###add center
  map <- map[,c("sample.id","CENTER")]
  # linker<-linker[,c("sample.id","CENTER")]
  
  names(map)<-c("sample.id","center1")
  # names(linker)<-c("sample.id","center1")
  
  dups = merge(dups, map, by.x = "ID1", by.y = "sample.id", all.x = T)
  # d<-merge(d,linker,by.x="ID1",by.y="sample.id",all.x=T)
  
  names(map)<-c("sample.id","center2")
  # names(linker)<-c("sample.id","center2")
  
  dups = merge(dups, map, by.x = "ID2", by.y = "sample.id", all.x = T)
  # d<-merge(d,linker,by.x="ID2",by.y="sample.id",all.x=T)
  
  ####add call rate
  qc <- read.table(paste(f.dir,"freeze_5b_pass_minDP10_qc.tsv",sep="/"),header=T)
  # q<-read.table("freeze_5b_pass_minDP10_qc.tsv",header=T)
  
  qc = qc[,c("s","sa.qc.callRate")]
  # q<-q[,c("s","sa.qc.callRate")]
  
  names(qc) <- c("s","cr1")
  # names(q)<-c("s","cr1")
  
  dups = merge(dups,qc,by.x="ID1",by.y="s",all.x=T)
  # d<-merge(d,q,by.x="ID1",by.y="s",all.x=T)
  
  names(qc) <- c("s","cr2")
  dups = merge(dups,qc,by.x="ID2",by.y="s",all.x=T)
  
  ###add center percentage
    
  for (j in 1:nrow(dups)) {
      dups$p1[j]<-ifelse(dups$study1[j]%in%row.names(t),t[as.character(dups$study1[j]),as.character(dups$center1[j])],NA)
      dups$p2[j]<-ifelse(dups$study2[j]%in%row.names(t),t[as.character(dups$study2[j]),as.character(dups$center2[j])],NA)
      
  }
    
    ##remove controls####
  dups <- dups[which(!(dups$study1 %in% c("CONTROL")) & !(dups$study2%in%c("CONTROL"))),]
    
  # Within study duplicates
  ## Choose the duplicate from the sequencing center that has the highest sequencing percentage.
  
  ## print table of duplicate pairs in the duplicates file
  print(table(dups$study1,dups$study2))
  
  #####FOR WITHIN THE SAME STUDY: Keeps monozygotic twins and chooses duplicates with highest sequencing call rate
  for (j in 1:nrow(dups)) {
      if (as.character(dups$study1[j])==as.character(dups$study2[j])) {
        if (!is.na(dups$MZtwinID[j])) { ###keep MZtwin
          dups$keep1[j]<-1;dups$keep2[j]<-1
        } else if (as.character(dups$center1[j])!=as.character(dups$center2[j]) & !is.na(dups$p1[j])& !is.na(dups$p2[j])) {
          dups$keep1[j]<-ifelse(dups$p1[j]>=dups$p2[j],1,0)
          dups$keep2[j]<-1-dups$keep1[j]
        } else if (as.character(dups$center1[j])==as.character(dups$center2[j])){
          dups$keep1[j]<-ifelse(dups$cr1[j]>=dups$cr2[j],1,0)
          dups$keep2[j]<-1-dups$keep1[j]
        }###Choose duplicate with highest sequencing call rate
      }
      
      ###FOR ACROSS STUDIES: ranks duplicates according to cohort type, such as population-based or with longest follow up#####
      #ARIC>  DHS, GOLDN, GENOA, HyperGEN, Mayo_VTE, MESA, WHI, JHS, GeneSTAR
      else if ((as.character(dups$study1[j])=='ARIC'&as.character(dups$study2[j])%in%c("DHS","GOLDN","GENOA","HyperGEN","Mayo_VTE", "MESA", "WHI", "JHS", "GeneSTAR")) |
               (as.character(dups$study2[j])=='ARIC'&as.character(dups$study1[j])%in%c("DHS","GOLDN","GENOA","HyperGEN","Mayo_VTE", "MESA", "WHI", "JHS", "GeneSTAR")))
      {
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='ARIC',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
      }#FHS>   COPDGene, MGH_AF, WHI
      else if ((as.character(dups$study1[j])=='FHS'&as.character(dups$study2[j])%in%c("COPDGene","MGH_AF","WHI"))|
               (as.character(dups$study2[j])=='FHS'&as.character(dups$study1[j])%in%c("COPDGene","MGH_AF","WHI")))
      {
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='FHS',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
      }#GeneSTAR>  COPDGene
      else if ((as.character(dups$study1[j])=='GeneSTAR'&as.character(dups$study2[j])%in%c("COPDGene"))|
               (as.character(dups$study2[j])=='GeneSTAR'&as.character(dups$study1[j])%in%c("COPDGene")))
      {
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='GeneSTAR',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
      }#HyperGEN>          COPDGene, DHS
      else if ((as.character(dups$study1[j])=='HyperGEN'&as.character(dups$study2[j])%in%c("COPDGene", "DHS"))|
               (as.character(dups$study2[j])=='HyperGEN'&as.character(dups$study1[j])%in%c("COPDGene", "DHS")))
      {
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='HyperGEN',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
      }
      ##JHS>       GENOA
      else if ((as.character(dups$study1[j])=='JHS'&as.character(dups$study2[j])%in%c("GENOA"))|
               (as.character(dups$study2[j])=='JHS'&as.character(dups$study1[j])%in%c("GENOA")))
      {
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='JHS',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
      }
      ##MESA>   COPDGene, DHS, GENOA, HyperGEN, GeneSTAR
      else if ((as.character(dups$study1[j])=='MESA'&as.character(dups$study2[j])%in%c("COPDGene", "DHS", "GENOA", "HyperGEN", "GeneSTAR"))|
               (as.character(dups$study2[j])=='MESA'&as.character(dups$study1[j])%in%c("COPDGene", "DHS", "GENOA", "HyperGEN", "GeneSTAR")))
      {
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='MESA',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
      }
      ###VU_AF> WGHS
      else if ((as.character(dups$study1[j])=='VU_AF'&as.character(dups$study2[j])%in%c("WGHS"))|
               (as.character(dups$study2[j])=='VU_AF'&as.character(dups$study1[j])%in%c("WGHS")))
      {
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='VU_AF',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
      }
      ###WHI> COPDGene, HyperGEN, Mayo_VTE, MESA
      else if ((as.character(dups$study1[j])=='WHI'&as.character(dups$study2[j])%in%c("COPDGene", "HyperGEN", "Mayo_VTE", "MESA"))|
               (as.character(dups$study2[j])=='WHI'&as.character(dups$study1[j])%in%c("COPDGene", "HyperGEN", "Mayo_VTE", "MESA")))
      {
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='WHI',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
      } 
      ####CHS > MESA, WHI 
      else if ((as.character(dups$study1[j])=='CHS'&as.character(dups$study2[j])%in%c("MESA", "WHI"))|
               (as.character(dups$study2[j])=='CHS'&as.character(dups$study1[j])%in%c("MESA", "WHI")))
      {
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='CHS',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
      } 
      ####MGH_AF > VU_AF, CCAF, Partners
      else if ((as.character(dups$study1[j])=='MGH_AF'&as.character(dups$study2[j])%in%c("CCAF", "VU_AF", "Partners"))|
               (as.character(dups$study2[j])=='MGH_AF'&as.character(dups$study1[j])%in%c("CCAF", "VU_AF", "Partners")))
      {
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='MGH_AF',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
      } 
      else { dups$keep1[j]<-NA; dups$keep2[j]<-NA}
    }
    
    
    dups <- dups[c("ID1", "ID2", "study1", "study2", "MZtwinID", "center1", "center2", "cr1", "cr2", "p1", "p2", "keep1", "keep2")]
    ### for reordering ID1 and ID2. Used column names for clarity but can use column numbers for simplicity. 
    
    
    write.table(dups,paste(f.dir,"/",out.pref,"duplicates.txt",sep=""),row.names=F,col.names=T,quote=F,sep='\t')
    
    
  # }


if(length(args) == 5) {
     id.col <- args[3]
     ped.file <- args[4]
     trait <- args[5]
  
  
  ped <- read.table(paste(f.dir,ped.file,sep="/"),header=T,sep=",",as.is=T) #n=50168
  table(duplicated(ped$unique_subject_key),useNA = "always")
  
  # dups <- read.table("duplicates.txt",header=T)
  # dups <- read.table(paste(f.dir,"freeze5b_duplicates.txt",sep="/"),header=T)
  
  if(! trait %in% colnames(ped)) {
   print(trait)
    stop("Trait entered must match column name in phenotype file")
  }
  
    ####ACROSS STUDIES: indicates which duplicate to remove based on missing trait data and cohort type 
    
    ####TRAIT
    p3 <- ped[,c(id.col, trait)]
    trait1 = paste(trait,"1",sep="_")
    trait2 = paste(trait,"2",sep="_")
    names(p3) <- c(id.col, "TRAIT1")
    dups<-merge(dups,p3,by.x = 'ID1',by.y = id.col,all.x=T)
    
    p3 <- ped[,c(id.col, trait)]
    names(p3) <- c(id.col, "TRAIT2")
    dups<-merge(dups,p3,by.x = 'ID2',by.y = id.col,all.x=T)
    
    dups$keep1_TRAIT<-NA
    dups$keep2_TRAIT<-NA
    
    dups$reason_TRAIT <- NA
    
    table(is.na(dups$TRAIT1),is.na(dups$TRAIT2),useNA = "always")
    #      FALSE TRUE <NA> # 7/19/2018 AKM
    #FALSE   327   17    0
    #TRUE     12   20    0
    #<NA>      0    0    0
    
    # old.dups <- dups  # for debugging
    for (j in 1:nrow(dups)){
      
      # For situations when both TRAIT1 and TRAIT2 are not missing
      if(!is.na(dups$TRAIT1[j]) & !is.na(dups$TRAIT2[j])) {
        
        # use same decision as before
        dups$keep1_TRAIT[j] <- dups$keep1[j]
        dups$keep2_TRAIT[j] <- dups$keep2[j]
        
        ##for additional output files
        if(!is.na(dups$keep1[j]) & !is.na(dups$keep2[j])){
          
          if(dups$keep1[j]==1 & as.character(dups$study1[j]) != as.character(dups$study2[j])) {
            dups$reason_TRAIT[j] <- "Both TRAIT1 and TRAIT2 not NA; used duplicate from population-based and/or older study"
          }
          else if(dups$keep2[j]==1 & as.character(dups$study1[j]) != as.character(dups$study2[j])) {
            dups$reason_TRAIT[j] <- "Both TRAIT1 and TRAIT2 not NA; used duplicate from population-based and/or older study"
          }
          else if(dups$keep1[j]==1 & dups$keep2[j]==1){
            dups$reason_TRAIT[j] <- "Both TRAIT1 and TRAIT2 not NA; kept both; monozygotic twins" 
          }
          else if(dups$keep1[j]==1 & dups$keep2[j]==0 & as.character(dups$study1[j])==as.character(dups$study2[j])){
            dups$reason_TRAIT[j] <- "Both TRAIT1 and TRAIT2 not NA; same cohort; kept duplicate with higher callrate"
          }
          else if(dups$keep2[j]==1&dups$keep1[j]==0 & as.character(dups$study1[j])==as.character(dups$study2[j])){
            dups$reason_TRAIT[j] <- "Both TRAIT1 and TRAIT2 not NA; same cohort; kept duplicate with higher callrate"
          }
        }
      } else if(is.na(dups$TRAIT1[j]) & !is.na(dups$TRAIT2[j])) # for situations where TRAIT1 is missing and TRAIT2 is not missing
      {
         dups$keep1_TRAIT[j]<-0
         dups$keep2_TRAIT[j]<-1
         dups$reason_TRAIT[j] <- "TRAIT1 individual missing trait data"
       } else if(!is.na(dups$TRAIT1[j]) & is.na(dups$TRAIT2[j])) # for situations where TRAIT1 is not missing and TRAIT2 is missing
       {
         dups$keep1_TRAIT[j]<-1
         dups$keep2_TRAIT[j]<-0
         dups$reason_TRAIT[j] <- "TRAIT2 individual missing trait data"
      }
      # changes so that we keep individuals with NA for trait in phenotype file
      else if (is.na(dups$TRAIT1[j]) & is.na(dups$TRAIT2[j]))  # for situations where TRAIT1 is missing and TRAIT2 is missing 
       {
         dups$keep1_TRAIT[j]<-dups$keep1[j]
         dups$keep2_TRAIT[j]<-dups$keep2[j]
         dups$reason_TRAIT[j] <- "TRAIT1 and TRAIT2 individuals both missing trait data"
      } else { 
         dups$keep1_TRAIT[j]<-NA
         dups$keep2_TRAIT[j]<-NA
         dups$reason_TRAIT[j] <- "something else; should be fixed"
      }
    }
    
    print(table(keep1_TRAIT=dups$keep1_TRAIT,keep2_TRAIT=dups$keep2_TRAIT,dups$reason_TRAIT,useNA='always'))
  #  print("Error if NAs occur")
    
    
    ####OUTPUT: to subset for TRAIT and to remove duplicates as determined by algorithm####
    notkeep1 <- subset(dups[,c("ID1","keep1_TRAIT")], keep1_TRAIT==0)
    notkeep2 <- subset(dups[,c("ID2","keep2_TRAIT")], keep2_TRAIT==0)
    names(notkeep1) <- c("ID","notkeep")
    names(notkeep2) <- c("ID","notkeep")
    
    notkeep <- rbind(notkeep1,notkeep2)
    notkeep.ids <- as.character(notkeep$ID)
  
    ped$keep_trait <- rep(1,NROW(ped))
    ped$keep_trait[ped[,id.col] %in% notkeep.ids] <- 0
    
    table(ped$keep_trait,ped$study,useNA="always")
    
    # does this take care of duplicate IDs from trait harmonization step?
    fulldata.dup.all <- ped[which(ped$unique_subject_key %in% ped$unique_subject_key[which(duplicated(ped$unique_subject_key))]),]
    table(duplicated(fulldata.dup.all$unique_subject_key),useNA = "always")
    
    fulldata.dup.all.dups.removed <- fulldata.dup.all[which(fulldata.dup.all$keep_trait==1),]
    # ped <-merge(ped, not_keep,by.x=id.col,by.y="ID",all.x=T)
    table(duplicated(fulldata.dup.all.dups.removed$unique_subject_key),useNA = "always")

    ped <- ped[!is.na(ped[,trait]),]
    
    table(duplicated(ped$unique_subject_key[ped$keep_trait == 1]),useNA = "always")
    
    
    write.table(ped[ped$keep_trait == 1,],paste(f.dir,"/",out.pref,".no.duplicates.csv",sep=""),row.names=F,col.names=T,quote=F,sep=',')
    
    ## NOTE: 13MAR2018 sex.x and sex.y from fulldata (harmonization) file is dropped
    
    ## for pooled phenotypes file with additional columns (for inspection purposes)
    
    # write.table(p,"./TRAIT/all_TRAIT.ped",row.names=F,col.names=T,quote=F,sep='\t')
    
    ###for the duplicates.txt file with additional trait-specific keep columns, if desired: 
    # dups <- dups[c("ID1", "ID2", "study1", "study2", "MZtwinID", "center1", "center2", "cr1", "cr2", "p1", "p2", "keep1", "keep2", "TRAIT1", "TRAIT2", "keep1_TRAIT", "keep2_TRAIT", "reason_TRAIT")]
    
    # write.table(dups,"TRAIT/duplicates_TRAIT.txt",row.names=F,col.names=T,quote=F,sep='\t')
    
    #### for lists of duplicates that were removed####
    
    names <- c("ID", "reason")
    
    removed1_TRAIT <- subset(dups, keep1_TRAIT==0)
    removed1_TRAIT <- removed1_TRAIT[c("ID1", "reason_TRAIT")]
    colnames(removed1_TRAIT) <- names
    
    removed2_TRAIT <-subset(dups, keep2_TRAIT==0)
    removed2_TRAIT <- removed2_TRAIT[c("ID2", "reason_TRAIT")]
    colnames(removed2_TRAIT) <- names
    
    removed_TRAIT <- rbind(removed1_TRAIT, removed2_TRAIT)

    write.table(dups,paste(f.dir,"/",out.pref,"duplicates.with.trait.",trait,".txt",sep=""),row.names=F,col.names=T,quote=F,sep='\t')
    
    write.table(removed_TRAIT, paste(f.dir,"/",out.pref,".removed.IDs.csv",sep=""), row.names=F, col.names=T, quote=F, sep=',')
    
}

    