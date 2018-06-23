input_args <- commandArgs(trailingOnly=T)
print(input_args)
ped.file <- input_args[1]
outcome <- input_args[2]
covars.continuous <- unlist(strsplit(input_args[3],","))
covars.categorical <- unlist(strsplit(input_args[4],","))
label <- input_args[5]
cohort_column <- input_args[6]
ancestry_column <- input_args[7]
f.dir <- input_args[8]

# Load phenotype data
library(data.table)
ped.data <- fread(paste(f.dir,"/",ped.file,sep=""),header=T,stringsAsFactors=FALSE,showProgress=TRUE,data.table=FALSE)
ped.data = na.omit(as.data.frame(ped.data[,unique(c(cohort_column,outcome,covars.continuous,covars.categorical,ancestry_column)),drop=F]))

# Determine outcome type
is.continuous <- ifelse(length(unique(ped.data[,outcome])) < 4 ,F,T)

# Get number of possible outcomes
if (!(is.continuous)){
  outcome.vals <- sort(unique(ped.data[,outcome]))
} else {
  covars.continuous = c(outcome,covars.continuous)
}

# Remove cohort column and any study based columns if its in covars
#covars = covars[!(covars %in% c(cohort_column,"STUDY_ANCESTRY","study_ancestry","topmed_project"))]

# Determine types of covars
#val.type <- apply(ped.data[,covars], 2, function(x) ifelse(length(unique(x)) <= 2 ,"categorical","continuous"))

# Gather each cohorts stats in a list to be combined later
stats <- list()

# Two different loops for catgorical vs continuous
if (!(is.continuous)){
  # Loop through cohorts
  for (study in unique(ped.data[,cohort_column])){
    print(study)
    for(ancestry in unique(ped.data[which(ped.data[,cohort_column]==study),ancestry_column])) {
    print(ancestry)
    # Generate the row names that we'll use
    row_names <- c(paste(study,ancestry,"total",sep=" "))
    for (v in outcome.vals){
      row_names <- c(row_names,paste(study,ancestry,":",outcome,"=",v,sep=" "))
    }
    
    # Subset the phenotype data down to only the study we're currently using
    ped.cur = ped.data[ped.data[,cohort_column] == study & ped.data[,ancestry_column] == ancestry,]
    
    # Get the total number of samples in this cohort
    total <- length(ped.cur[,1])
    
    # Start collecting numbers of samples per outcome condition
    samp <- c(total)
    for (v in outcome.vals){
      samp <- c(samp, length(ped.cur[ped.cur[,outcome] == v,1]))
    }
    
    # This is where we will store all stats for this cohort
    all_dat <- data.frame(V1 = row_names, Samples=samp)
    
    # Loop through the continuous covariates to calculate mean, median, std, min, max
    for (c in covars.continuous){
      if (c == cohort_column){
        next
      }
      
      # Store the values for the whole cohort first
      means <- c(mean(ped.cur[,c]))
      medians <- c(median(ped.cur[,c]))
      sds <- c(sd(ped.cur[,c]))
      mins <- c(min(ped.cur[,c]))
      maxes <- c(max(ped.cur[,c]))
      
      # Loop through each outcome condition and calculate values
      for (v in outcome.vals){
        print(v)
        ped.new <- ped.cur[ped.cur[,outcome] == v,c]
        means = c(means,mean(ped.new))
        medians = c(medians, median(ped.new))
        sds = c(sds, sd(ped.new))
        mins = c(mins, min(ped.new))
      }
      
      # Store the values back in the whole data frame
      all_dat <- cbind(all_dat,means,medians,sds,mins,maxes)
      
    }
    
    # Loop through the categorical covariates
    for (c in covars.categorical) {
      if (c == cohort_column){
        next
      }
      
      # Get all possible values of the current covar
      curvals <- unique(ped.data[,c])
      
      # Loop through the possible values of this covar
      for (cv in curvals){
        
        # Subset the phenotypes to only cases with this covar value
        pd.c <- ped.cur[ped.cur[,c] == cv,]
        
        # Total percent with this value for the covar
        percents <- c(length(pd.c[,1])/total)
        
        # Calculate percents of samples with this value over total with given condition
        for (v in outcome.vals){
          ped.new <- pd.c[pd.c[,outcome] == v,c]
          percents <- c(percents,length(ped.new)/length(ped.cur[ped.cur[,outcome] == v,1]))
        }
        
        # Store back in the total data frame
        all_dat <- cbind(all_dat,percents)  
      }
    }
    
    # Store the cohort specific stats in the list
    stats[[length(stats)+1]] <- all_dat
  }
  }
  # Get totals for all cohorts
  print("All cohorts")
  row_names <- c("All cohorts")
  for (v in outcome.vals){
    print(v)
    print(paste("All cohorts",":",outcome,"=",v,sep=" "))
    row_names <- c(row_names,paste("All cohorts",":",outcome,"=",v,sep=" "))
  }
  
  # total samples
  total <- length(ped.data[,1])
  samp <- c(total)
  for (v in outcome.vals){
    samp <- c(samp, length(ped.data[ped.data[,outcome] == v,1]))
  }
  print(row_names)
  all_dat <- data.frame(V1 = row_names, Samples=samp)
  print("continuous covariates")
  # Loop through the continuous covariates to calculate mean, median, std, min, max
  for (c in covars.continuous){
    if (c == cohort_column){
      next
    }
    
    # Store the values for the whole cohort first
    means <- c(mean(ped.data[,c]))
    medians <- c(median(ped.data[,c]))
    sds <- c(sd(ped.data[,c]))
    mins <- c(min(ped.data[,c]))
    maxes <- c(max(ped.data[,c]))
    
    # Loop through each outcome condition and calculate values
    for (v in outcome.vals){
      ped.new <- ped.data[ped.data[,outcome] == v,c]
      means = c(means,mean(ped.new))
      medians = c(medians, median(ped.new))
      sds = c(sds, sd(ped.new))
      mins = c(mins, min(ped.new))
    }
    
    # Store the values back in the whole data frame
    print(cbind(means,medians,sds,mins,maxes))
    all_dat <- cbind(all_dat,means,medians,sds,mins,maxes)
    
  }
  print("categorical covariates")
  # Loop through the categorical covariates
  for (c in covars.categorical){
    if (c == cohort_column){
      next
    }
    
    # Get all possible values of the current covar
    curvals <- unique(ped.data[,c])
    
    # Loop through the possible values of this covar
    for (cv in curvals){
      
      # Subset the phenotypes to only cases with this covar value
      pd.c <- ped.data[ped.data[,c] == cv,]
      
      # Total percent with this value for the covar
      percents <- c(length(pd.c[,1])/total)
      
      # Calculate percents of samples with this value over total with given condition
      for (v in outcome.vals){
        ped.new <- pd.c[pd.c[,outcome] == v,c]
        percents <- c(percents,length(ped.new)/length(ped.data[ped.data[,outcome] == v,1]))
      }
      
      # Store back in the total data frame
      all_dat <- cbind(all_dat,percents)  
    }
    stats[[length(stats)+1]] <- all_dat
  
  } 
  print("Done with categorical")
# If we have a continuous outcome
} else {
  for (study in unique(ped.data[,cohort_column])){
    
    # Subset the phenotype data down to only the study we're currently using
    ped.cur = ped.data[ped.data[,cohort_column] == study,]
    
    # Get the total number of samples in this cohort
    total <- length(ped.cur[,1])
    
    # This is where we will store all stats for this cohort
    all_dat <- data.frame(V1 = study, Samples=total)
    
    # Loop through the continuous covariates to calculate mean, median, std, min, max
    for (c in covars.continuous){
      if (c == cohort_column){
        next
      }
      
      # Store the values for the whole cohort first
      means <- c(mean(ped.cur[,c]))
      medians <- c(median(ped.cur[,c]))
      sds <- c(sd(ped.cur[,c]))
      mins <- c(min(ped.cur[,c]))
      maxes <- c(max(ped.cur[,c]))
      
      # Store the values back in the whole data frame
      all_dat <- cbind(all_dat,means,medians,sds,mins,maxes)
      
    }
    
    # Loop through the categorical covariates
    for (c in covars.categorical){
      if (c == cohort_column){
        next
      }
      
      # Get all possible values of the current covar
      curvals <- unique(ped.data[,c])
      
      # Loop through the possible values of this covar
      for (cv in curvals){
        
        # Subset the phenotypes to only cases with this covar value
        pd.c <- ped.cur[ped.cur[,c] == cv,]
        
        # Total percent with this value for the covar
        percents <- c(length(pd.c[,1])/total)
        
        # Store back in the total data frame
        all_dat <- cbind(all_dat,percents)  
      }
    }
    
    # Store the cohort specific stats in the list
    stats[[length(stats)+1]] <- all_dat
  }
  
  # Get the total number of samples in this cohort
  total <- length(ped.data[,1])
  
  # This is where we will store all stats for this cohort
  all_dat <- data.frame(V1 = "All cohorts", Samples=total)
  
  # Loop through the continuous covariates to calculate mean, median, std, min, max
  for (c in covars.continuous){
    if (c == cohort_column){
      next
    }
    
    # Store the values for the whole cohort first
    means <- c(mean(ped.data[,c]))
    medians <- c(median(ped.data[,c]))
    sds <- c(sd(ped.data[,c]))
    mins <- c(min(ped.data[,c]))
    maxes <- c(max(ped.data[,c]))
    
    # Store the values back in the whole data frame
    all_dat <- cbind(all_dat,means,medians,sds,mins,maxes)
    
  }
  
  # Loop through the categorical covariates
  for (c in covars.categorical){
    if (c == cohort_column){
      next
    }
    
    # Get all possible values of the current covar
    curvals <- unique(ped.data[,c])
    
    # Loop through the possible values of this covar
    for (cv in curvals){
      
      # Subset the phenotypes to only cases with this covar value
      pd.c <- ped.data[ped.data[,c] == cv,]
      
      # Total percent with this value for the covar
      percents <- c(length(pd.c[,1])/total)
      
      # Store back in the total data frame
      all_dat <- cbind(all_dat,percents)  
    }
  }
  
  # Store the cohort specific stats in the list
  stats[[length(stats)+1]] <- all_dat
}

# Determine the column and row names for the data frame
# Continuous covars go first
cont_vals <- covars.continuous
cont_cols <- c()
for (c in cont_vals){
  if (c == cohort_column){
    next
  }
  cont_cols <- c(cont_cols, paste(c,"mean",sep=" "), paste(c,"median",sep=" "), paste(c,"STD",sep=" "), paste(c,"min",sep=" "), paste(c,"max",sep=" "))
}

# Next are categorical covars
dich_vals <- covars.categorical
dich_cols <- c()
for (c in dich_vals){
  if (c == outcome){
    next
  }
  uvals <- unique(ped.cur[,c])
  dich_cols <- c(dich_cols, paste(c,"%",uvals[1],sep=" "),paste(c,"%",uvals[2],sep=" "))
}

# Add samples to headers and combine
col_names <- c("Samples",cont_cols,dich_cols)

# Combine list to data frame containing all cohorts
stats.df <- as.data.frame(do.call(rbind,stats))#, row.names = row_names)

# Move the row names to row.names of data frame
row.names(stats.df) <- stats.df$V1

# Remove row names column
stats.df <- stats.df[,2:length(stats.df[1,])]

# Assign the correct column names
colnames(stats.df) <- col_names

# Write the final stats out to a csv file
fwrite(stats.df,file=paste(f.dir,"/",label,"_stats.csv",sep=""),sep=",",row.names = T)

# Make some plots
# Load ggplot for plots
library(ggplot2)

# Gather continuous covariates
quant_covars <- covars.continuous

# Make sure we dont include the cohort column
quant_covars = quant_covars[quant_covars != cohort_column]

# Initialize the pdf
pdf(paste(f.dir,"/",label,"_plots.pdf",sep=""),width=20)
# if (length(quant_covars) > 1){
#   layout(matrix(seq(1,6*(length(quant_covars)-1)),nrow=length(quant_covars)-1,ncol=6,byrow=T))
# }else{
#   layout(matrix(seq(1,6*(length(quant_covars))),nrow=length(quant_covars),ncol=6,byrow=T))
# }
ped.data$ancestry_column <- ped.data[,ancestry_column]
# For each continuous covariate
for (i in covars.continuous){
  print(i)
  
  # Plot all samples
  plot <- ggplot(ped.data, aes_string(cohort_column, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) 
  print(plot + geom_boxplot(aes(fill=factor(ancestry_column))) + labs(title = "All samples",x="Cohort",y=i))
  
  # Plot samples devided by sex
  plot <- ggplot(ped.data, aes_string(cohort_column, i)) + theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) + facet_grid(. ~ ancestry_column) 
  print(plot + geom_boxplot(aes(fill=factor(sex))) + labs(title = "All samples by sex",x="Cohort",y=i,fill="sex"))
}


dev.off()

