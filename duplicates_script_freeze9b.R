process_duplicates <- function( f.dir, f9_dup, f9_sample, id.col=NULL, ped=ped, qc=qc, trait=NULL, out.pref, t2d.trait=NULL) {
  
  ## load duplicates 
  dups <- read.table(paste(f.dir,f9_dup,sep="/"),header=T,stringsAsFactors = F)
  dups <- subset(dups, study1 %in% unique(ped$study) &
                   study2 %in% unique(ped$study)) 
  print(dim(dups))

  ## load map file
  map <- read.table(paste(f.dir,f9_sample,sep="/"),header=T,stringsAsFactors = F)
  map <- subset(map, study %in% unique(ped$study)) 
  
  t<-table(map$study,map$seq_center)
  t<-prop.table(t,1)
  t
    
  if(! trait %in% colnames(ped)) {
    print(trait)
    stop("Trait entered must match column name in phenotype file")
  }
  
  ### add center
  map <- map[,c("sample.id","seq_center","topmed_project")]
  names(map)<-c("sample.id","center1","topmed_project1")
  dups <- merge(dups, map, by.x = "ID1", by.y = "sample.id", all.x = T)
  names(map)<-c("sample.id","center2","topmed_project2")
  dups <- merge(dups, map, by.x = "ID2", by.y = "sample.id", all.x = T)

  #### add call rate
  qc = qc[,c("s","sa.qc.callRate")]
  names(qc) <- c("s","cr1")
  dups = merge(dups,qc,by.x="ID1",by.y="s",all.x=T)
  names(qc) <- c("s","cr2")
  dups = merge(dups,qc,by.x="ID2",by.y="s",all.x=T)
  
    
  ###add center percentage
  for (j in 1:nrow(dups)) {
    dups$p1[j]<-ifelse(dups$study1[j]%in%row.names(t),t[as.character(dups$study1[j]),as.character(dups$center1[j])],NA)
    dups$p2[j]<-ifelse(dups$study2[j]%in%row.names(t),t[as.character(dups$study2[j]),as.character(dups$center2[j])],NA)
  }
  
  ### add trait column
  ####ACROSS STUDIES: indicates which duplicate to remove based on missing trait data and cohort type 
  ####TRAIT
  print(table(dups[,'ID1'] %in% ped[,id.col] ))
  print(table(dups[,'ID2'] %in% ped[,id.col] ))
  print(table(dups[,'ID1'] %in% ped[,id.col],dups[,'ID2'] %in% ped[,id.col] ))
    
  #only deal with the duplicates that are in the current ped file
  dups <- dups[which(dups[,'ID1'] %in% ped[,id.col] & dups[,'ID2'] %in% ped[,id.col]),]
    
  p3 <- ped[,c(id.col, trait)]
  names(p3) <- c(id.col, "TRAIT1")
  dups<-merge(dups,p3,by.x = 'ID1',by.y = id.col,all.x=T)
  p3 <- ped[,c(id.col, trait)]
  names(p3) <- c(id.col, "TRAIT2")
  dups<-merge(dups,p3,by.x = 'ID2',by.y = id.col,all.x=T)
    
  print(dim(dups))
  
  dups$keep1_TRAIT<-NA;    
  dups$keep2_TRAIT<-NA
  dups$reason_TRAIT <- NA
  
  table(is.na(dups$TRAIT1),is.na(dups$TRAIT2),useNA = "always")
  
  # Within study duplicates
  ## Choose the duplicate from the sequencing center that has the highest sequencing percentage.
  ## print table of duplicate pairs in the duplicates file
  print(table(dups$study1,dups$study2))
  dups$keep1 <- NA; 
  dups$keep2 <- NA
  dups$REASON_keep <- NA
  
  #####FOR WITHIN THE SAME STUDY: Keeps monozygotic twins and chooses duplicates with highest sequencing call rate
  for (j in 1:nrow(dups)) {
    if (as.character(dups$study1[j])==as.character(dups$study2[j])) {
      if (!is.na(dups$MZtwinID[j])) { ###keep MZtwin
        dups$keep1[j]<-1;
        dups$keep2[j]<-1
        dups$REASON_keep[j] <- "keep MZtwin"
      } else if (as.character(dups$center1[j])!=as.character(dups$center2[j]) & !is.na(dups$p1[j])& !is.na(dups$p2[j])) {
        dups$keep1[j]<-ifelse(dups$p1[j]>=dups$p2[j],1,0)
        dups$keep2[j]<-1-dups$keep1[j]
        dups$REASON_keep[j] <- "Samples are from different centers, choose sample with higher proportion of cohort"
      } else if (as.character(dups$center1[j])==as.character(dups$center2[j])){
        dups$keep1[j]<-ifelse( (dups$cr1[j]>=dups$cr2[j]) | (!is.na(dups$cr1[j])&is.na(dups$cr2[j])),1,0)
        dups$keep1[j]<-ifelse( (is.na(dups$cr1[j])&!is.na(dups$cr2[j])),0,dups$keep1[j])
        dups$keep2[j]<-1-dups$keep1[j]
        dups$REASON_keep[j] <- "Samples are from same centers, Choose duplicate with highest sequencing call rate / choose sample with meeting callrate"
      }
   } else if ((as.character(dups$study1[j])=='ARIC' & as.character(dups$study2[j])%in%c("DHS","GOLDN","GENOA","HyperGEN","Mayo_VTE", "MESA", "WHI", "JHS", "GeneSTAR", "CHS", "COPDGene")) |
                 (as.character(dups$study2[j])=='ARIC' & as.character(dups$study1[j])%in%c("DHS","GOLDN","GENOA","HyperGEN","Mayo_VTE", "MESA", "WHI", "JHS", "GeneSTAR", "CHS", "COPDGene")))
   {
        ###FOR ACROSS STUDIES: ranks duplicates according to cohort type, such as population-based or with longest follow up#####
        #ARIC>  DHS, GOLDN, GENOA, HyperGEN, Mayo_VTE, MESA, WHI, JHS, GeneSTAR, CHS, COPDGene
        
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='ARIC',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
        dups$REASON_keep[j] <- "ARIC>  DHS, GOLDN, GENOA, HyperGEN, Mayo_VTE, MESA, WHI, JHS, GeneSTAR, CHS, COPDGene"
    } else if ((as.character(dups$study1[j])=='FHS'&as.character(dups$study2[j])%in%c("COPDGene","MGH_AF","WHI"))|
                 (as.character(dups$study2[j])=='FHS'&as.character(dups$study1[j])%in%c("COPDGene","MGH_AF","WHI")))
    {
        #FHS>   COPDGene, MGH_AF, WHI
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='FHS',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
        dups$REASON_keep[j] <- "#FHS>   COPDGene, MGH_AF, WHI"
    } else if ((as.character(dups$study1[j])=='GeneSTAR'&as.character(dups$study2[j])%in%c("COPDGene"))|
                 (as.character(dups$study2[j])=='GeneSTAR'&as.character(dups$study1[j])%in%c("COPDGene")))
    {
        #GeneSTAR>  COPDGene
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='GeneSTAR',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
        dups$REASON_keep[j] <- "#GeneSTAR>  COPDGene"
    } else if ((as.character(dups$study1[j])=='HyperGEN'&as.character(dups$study2[j])%in%c("COPDGene", "DHS"))|
                 (as.character(dups$study2[j])=='HyperGEN'&as.character(dups$study1[j])%in%c("COPDGene", "DHS")))
    {
        #HyperGEN>          COPDGene, DHS
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='HyperGEN',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
        dups$REASON_keep[j] <- "#HyperGEN>          COPDGene, DHS"
    } else if ((as.character(dups$study1[j])=='JHS'&as.character(dups$study2[j])%in%c("GENOA"))|
                 (as.character(dups$study2[j])=='JHS'&as.character(dups$study1[j])%in%c("GENOA")))
    {
        ##JHS>       GENOA
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='JHS',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
        dups$REASON_keep[j] <- "#JHS>       GENOA"
    } else if ((as.character(dups$study1[j])=='MESA'&as.character(dups$study2[j])%in%c("COPDGene", "DHS", "GENOA", "HyperGEN", "GeneSTAR", "BioMe"))|
                 (as.character(dups$study2[j])=='MESA'&as.character(dups$study1[j])%in%c("COPDGene", "DHS", "GENOA", "HyperGEN", "GeneSTAR", "BioMe")))
    {
        ##MESA>   COPDGene, DHS, GENOA, HyperGEN, GeneSTAR, BioMe
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='MESA',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
        dups$REASON_keep[j] <- "#MESA>   COPDGene, DHS, GENOA, HyperGEN, GeneSTAR, BioMe"
    } else if ((as.character(dups$study1[j])=='VU_AF'&as.character(dups$study2[j])%in%c("WGHS"))|
                 (as.character(dups$study2[j])=='VU_AF'&as.character(dups$study1[j])%in%c("WGHS")))
    {
        ###VU_AF> WGHS
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='VU_AF',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
        dups$REASON_keep[j] <- "#VU_AF> WGHS"
    } else if ((as.character(dups$study1[j])=='WHI'&as.character(dups$study2[j])%in%c("COPDGene", "HyperGEN", "Mayo_VTE", "MESA", "THRV", 'HVH'))|
                 (as.character(dups$study2[j])=='WHI'&as.character(dups$study1[j])%in%c("COPDGene", "HyperGEN", "Mayo_VTE", "MESA", "THRV", 'HVH')))
    {
        ###WHI> COPDGene, HyperGEN, Mayo_VTE, MESA, THRV, HVH
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='WHI',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
        dups$REASON_keep[j] <- "#WHI> COPDGene, HyperGEN, Mayo_VTE, MESA, THRV, HVH"
    } 
      ####CHS > MESA, WHI, DHS, HyperGEN
      else if ((as.character(dups$study1[j])=='CHS'&as.character(dups$study2[j])%in%c("MESA", "WHI", 'DHS', "HyperGEN"))|
               (as.character(dups$study2[j])=='CHS'&as.character(dups$study1[j])%in%c("MESA", "WHI", 'DHS', "HyperGEN")))
    {
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='CHS',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
        dups$REASON_keep[j] <- "#CHS > MESA, WHI, DHS, HyperGEN"
    } else if ((as.character(dups$study1[j])=='MGH_AF'&as.character(dups$study2[j])%in%c("CCAF", "VU_AF", "Partners"))|
                 (as.character(dups$study2[j])=='MGH_AF'&as.character(dups$study1[j])%in%c("CCAF", "VU_AF", "Partners")))
    {
        ####MGH_AF > VU_AF, CCAF, Partners
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='MGH_AF',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
        dups$REASON_keep[j] <- "#MGH_AF > VU_AF, CCAF, Partners"
    } else if ((as.character(dups$study1[j])=='CARDIA'&as.character(dups$study2[j])%in%c("COPDGene", "HyperGEN"))|
                 (as.character(dups$study2[j])=='CARDIA'&as.character(dups$study1[j])%in%c("COPDGene", "HyperGEN")))
    { 
        ####CARDIA > COPDGene, HyperGEN
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='CARDIA',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
        dups$REASON_keep[j] <- "#CARDIA > COPDGene, HyperGEN"
    } else if ((as.character(dups$study1[j])=='BioMe'&as.character(dups$study2[j])%in%c("COPDGene"))|
                 (as.character(dups$study2[j])=='BioMe'&as.character(dups$study1[j])%in%c("COPDGene")))
    {
        ####BioMe > COPDGene
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='BioMe',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
        dups$REASON_keep[j] <- "#BioMe > COPDGene"
    } else if ((as.character(dups$study1[j])=='HCHS_SOL'&as.character(dups$study2[j])%in%c("BioMe"))|
                 (as.character(dups$study2[j])=='HCHS_SOL'&as.character(dups$study1[j])%in%c("BioMe")))
    {
        ####HCHS_Sol > BioMe
        dups$keep1[j]<-ifelse(as.character(dups$study1[j])=='HCHS_SOL',1,0)
        dups$keep2[j]<-1-dups$keep1[j]
        dups$REASON_keep[j] <- "#HCHS_Sol > BioMe"
    } 
  }
    
    #dups <- dups[c("ID1", "ID2", "study1", "study2", "MZtwinID", "center1", "center2", "cr1", "cr2", "p1", "p2", "keep1", "keep2","topmed_project1","topmed_project2","REASON_keep")]
    print("After Step 1, if there are any rows with missing keep1 and keep2, this needs to be fixed")
    print(table(keep1=dups$keep1,keep2=dups$keep2,dups$REASON_keep,useNA='always'))
    
    for (j in 1:nrow(dups)) {
      if((!is.na(dups$topmed_project1) & !is.na(dups$topmed_project2)) & dups$topmed_project1[j] == "Control" & dups$topmed_project2[j] == "Control") 
      {     # finally, make sure that "control samples" are not chosen
        dups$keep1_TRAIT[j] <- 0
        dups$keep2_TRAIT[j] <- 0  
        dups$reason_TRAIT[j] <- "Both ID1 and ID2 from CONTROL samples"
      } else if((!is.na(dups$topmed_project1) & !is.na(dups$topmed_project2)) & dups$topmed_project1[j] == "Control" & dups$topmed_project2[j] != "Control") {
        dups$keep1_TRAIT[j] <- 0
        dups$keep2_TRAIT[j] <- 1  
        dups$reason_TRAIT[j] <- "ID1 from CONTROL samples"
        
      } else if((!is.na(dups$topmed_project1) & !is.na(dups$topmed_project2)) & dups$topmed_project1[j] != "Control" & dups$topmed_project2[j] == "Control") {
        dups$keep1_TRAIT[j] <- 1
        dups$keep2_TRAIT[j] <- 0  
        dups$reason_TRAIT[j] <- "ID2 from CONTROL samples"
      } else if (is.na(dups$TRAIT1[j]) & is.na(dups$TRAIT2[j])) {
        # use same decision as before
        dups$keep1_TRAIT[j] <- dups$keep1[j]
        dups$keep2_TRAIT[j] <- dups$keep2[j]
        dups$reason_TRAIT[j] <- "Both Trait1 and Trait2 are missing; Using Same decision as keep1 and keep2"
      } else if(!is.na(dups$TRAIT1[j]) & !is.na(dups$TRAIT2[j])) { # For situations when both TRAIT1 and TRAIT2 are not missing
        if(is.null(t2d.trait)) {
          # use same decision as before
          dups$keep1_TRAIT[j] <- dups$keep1[j]
          dups$keep2_TRAIT[j] <- dups$keep2[j]
          dups$reason_TRAIT[j] <- "Quantitative trait; Using Same decision as keep1 and keep2"
        } else if(!is.null(t2d.trait)) {
          ## Added for incident T2D AKM 2/4/2022
          if(dups$TRAIT1[j] == 1 & dups$TRAIT2[j] == 0) {
            dups$keep1_TRAIT[j] <- 1
            dups$keep2_TRAIT[j] <- 0
            dups$reason_TRAIT[j] <- "T2D trait; Prioritizing incident T2D cases over control"
            
          } else if (dups$TRAIT1[j] == 0 & dups$TRAIT2[j] == 1) {
            dups$keep1_TRAIT[j] <- 0
            dups$keep2_TRAIT[j] <- 1
            dups$reason_TRAIT[j] <- "T2D trait; Prioritizing incident T2D cases over control"
          } else {
            # use same decision as before
            dups$keep1_TRAIT[j] <- dups$keep1[j]
            dups$keep2_TRAIT[j] <- dups$keep2[j]
            dups$reason_TRAIT[j] <- "T2D trait; Using Same decision as keep1 and keep2"
          }
        }
      } else if(is.na(dups$TRAIT1[j]) & !is.na(dups$TRAIT2[j])) { # for situations where TRAIT1 is missing and TRAIT2 is not missing
        dups$keep1_TRAIT[j]<-0
        dups$keep2_TRAIT[j]<-1
        dups$reason_TRAIT[j] <- "TRAIT1 individual missing trait data"
      } else if(!is.na(dups$TRAIT1[j]) & is.na(dups$TRAIT2[j])) { # for situations where TRAIT1 is not missing and TRAIT2 is missing
        dups$keep1_TRAIT[j]<-1
        dups$keep2_TRAIT[j]<-0
        dups$reason_TRAIT[j] <- "TRAIT2 individual missing trait data"
      } else { 
        dups$keep1_TRAIT[j]<-NA
        dups$keep2_TRAIT[j]<-NA
        dups$reason_TRAIT[j] <- "something else; should be fixed"
      }
    }
      print("After Step 2, if there are any rows with missing keep1_TRAIT and keep2_TRAIT, this needs to be fixed")
      print(table(keep1_TRAIT=dups$keep1_TRAIT,keep2_TRAIT=dups$keep2_TRAIT,dups$reason_TRAIT,useNA='always'))

      return(dups)
}
  
do.subsetting <- function(dups) {
    ####OUTPUT: to subset for TRAIT and to remove duplicates as determined by algorithm####
    notkeep1 <- subset(dups[,c("ID1","keep1_TRAIT")], keep1_TRAIT==0)
    notkeep2 <- subset(dups[,c("ID2","keep2_TRAIT")], keep2_TRAIT==0)
    names(notkeep1) <- c("ID","notkeep");    
    names(notkeep2) <- c("ID","notkeep")
    notkeep <- rbind(notkeep1,notkeep2)
    notkeep.ids <- as.character(notkeep$ID)
    
    return(notkeep.ids)
}
