##########################################################
############## Prognostic miRs modelling #################
##########################################################

library(dplyr)
library(tidyr)
library(gdata)
library(edgeR)
library(caret)
library(pROC)
library(ROCR)
library(TCGAbiolinks)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(gridExtra)

source('~/Documents/Scripts/scoreSignature.R') # By Ed curry


#############################################
################ FUNCTIONS ##################
#############################################


## 1) Compute ROC
makeROC <- function(x,y){
  data.frame(y=c(0,cumsum(1-y[order(x,decreasing=T)]))/sum(1-y),x=c(0,cumsum(y[order(x,decreasing=T)])/sum(y)))
}

## 2) Get AUC
getAUC <- function(ROC){
  sum((1-ROC$x[1:(nrow(ROC)-1)])*(ROC$y[2:nrow(ROC)]-ROC$y[1:(nrow(ROC)-1)]))
}

## 3) Plot ROC
plotROC <- function(ROC, col, roc.title=''){
  # Get AUC
  AUC <- sum((1-ROC$x[1:(nrow(ROC)-1)])*(ROC$y[2:nrow(ROC)]-ROC$y[1:(nrow(ROC)-1)]))
  # graph
  plot(x=ROC$x,y=ROC$y,type="l",xlab="1-specificity",ylab="sensitivity", lwd=2, col=col, main=roc.title, cex.main=1.4, cex.lab=1.4, cex.axis=1.4)
  abline(a=0, b= 1)
  legend('bottomright', legend=paste0('AUC = ', signif(AUC, 3)), bty="n", text.font=2, cex=1.4)
}


## 4) Permutation based p-value for AUC
AUC.pval <- function(ROC,nperm=1000){
  AUC <- getAUC(ROC)
  null.AUCs <- rep(NA,nperm)
  for(i in 1:nperm){
    null.AUCs[i] <- getAUC(data.frame(x=ROC$x,y=sample(ROC$y)))
  }
  (1+sum(null.AUCs>=AUC))/(1+nperm)
}


## 5) Elastic net regression
E.net.reg <- function(Train.df=counts.norm.Train$part1, Test.df=counts.norm.Test$part1, weight=c('Y','N'), clinical=c('Y','N','only')){
  
  temp.train <- Train.df %>% select_if(is.numeric)
  temp.test <- Test.df %>% select_if(is.numeric)
  remove.ind <- c(names(which(!is.finite(colSums(temp.train)))), names(which(!is.finite(colSums(temp.test)))))
  # remove these columns
  Train.df <- Train.df[,!colnames(Train.df) %in% remove.ind]
  Test.df <- Test.df[,!colnames(Test.df) %in% remove.ind]
  
  # Option to remove clinical variables
  if (clinical=='N'){
    
    Train.df <- Train.df %>%
      select(-c(ISUP.gleason,pT.grp,stage_event_psa))
    Test.df <- Test.df %>%
      select(-c(ISUP.gleason,pT.grp,stage_event_psa))
    print('Models with only isomiRs as predictors. No CLINICAL VARIABLES AS PREDICTORS')
  } else if (clinical=='Y'){
    print('Model with isomiRs + clinical variables also as predictors')
  } else if (clinical=='only'){
    
    Train.df <- Train.df %>%
      select(ISUP.gleason,pT.grp,stage_event_psa,BCR.event)
    Test.df <- Test.df %>%
      select(ISUP.gleason,pT.grp,stage_event_psa,BCR.event)
  }
  
  # Set Control values
  fitControl <- trainControl(method='cv', number=10, summaryFunction=twoClassSummary, classProbs=TRUE)
  
  # Select hyperparamters with CV + run model on traindf
  if (weight=='N'){
    
    print('No weightings applied')
    set.seed(123)
    fit <- train(BCR.event ~., data=Train.df, method='glmnet', trControl=fitControl, tuneLength=10, metric='ROC')
    
  } else if (weight=='Y') {
    
    print('proportional weightings applied')
    # weights as proportion
    wgts <- ifelse(Train.df$BCR.event=='BCR', (1-(sum(Train.df$BCR.event=='BCR')/length(Train.df$BCR.event))), (1-(sum(Train.df$BCR.event=='no_BCR')/length(Train.df$BCR.event))))
    
    set.seed(123)
    fit <- caret::train(BCR.event ~., data=Train.df, method='glmnet', trControl=fitControl, tuneLength=10, metric='ROC', weights=wgts)
  }
  
  # Predictions on train dataset
  train <- list()
  train$prob <- fit %>% predict(Train.df, type='prob')
  train$class <- fit %>% predict(Train.df, type='raw')
  train$con.mat <- confusionMatrix(data=train$class, Train.df$BCR.event)
  # Make predictions on test
  test <- list()
  test$prob <- fit %>% predict(Test.df, type='prob')
  test$class <- fit %>% predict(Test.df, type='raw')
  test$con.mat <- confusionMatrix(data=test$class, Test.df$BCR.event)
  ## Get AUC + ROC ED
  test$ROC.ED <- makeROC(x=test$prob[,2], y=ifelse(Test.df$BCR.event=='BCR', 1, 0))
  test$AUC.ED <- getAUC(test$ROC.ED)
  ## GET AUC + ROC pROC
  ## pROC
  test$ROC.pROC <- roc(Test.df$BCR.event, test$prob[,1],  levels=c('BCR', 'no_BCR'), direction='>', plot=F)
  test$AUC.pROC <- summary(auc(test$ROC.pROC))['Median']
  ## GET AUC + ROC ROCR
  test$ROC.ROCR <- prediction(test$prob[,2], Test.df$BCR.event)
  y <- performance(test$ROC.ROCR, measure='auc'); test$AUC.ROCR <- y@y.values[[1]]
  
  # Return results
  res <- list(model=fit, train=train, test=test)
  return(res)
}

## 6) Function to apply normalisation
normalise.TMM <- function(train.df=TrainSet.raw[[1]], test.df=TestSet.raw[[1]]){
  
  # Train
  dgList <- DGEList(counts=train.df, genes=rownames(train.df))
  dgList <- calcNormFactors(dgList, method='TMM')
  Train.cts <- cpm(dgList, normalized.lib.sizes = T, log=T)
  # Normalise as later on we are applying reglarised regression whose assumption is normally distributed data
  temp <- rownames(Train.cts)
  Train.cts <- data.frame(t(Train.cts))
  colnames(Train.cts) <- temp
  
  # test
  dgList <- DGEList(counts=test.df, genes=rownames(test.df))
  dgList <- calcNormFactors(dgList, method='TMM')
  Test.cts <- cpm(dgList, normalized.lib.sizes = T, log=T)
  temp <- rownames(Test.cts)
  Test.cts <- data.frame(t(Test.cts))
  colnames(Test.cts) <- temp
  
  return(list(Train.cts, Test.cts))
}

## 7) Conducting permutation Elastic Net regressions
null.perm <- function(Train.df=counts.norm.Train$part1, Test.df=counts.norm.Test$part1, weight='Y', clinical=c('only'), iterations=3){
  
  # 1) DF to save results  
  Perf.DF <- data.frame(matrix(NA, nrow=9, ncol=iterations))
  rownames(Perf.DF)[1:2] <- paste0('train:', c('ACC-random.model', 'ACC-correct.label'))
  rownames(Perf.DF)[3:9] <- paste0('test:', c('AUC.ED','AUC.pROC','AUC.ROCR', 'ACC','Sensitivity.cor', 'Specificity.cor', 'Precision.cor')) 
  
  # 2) Save correct labels
  corr.labels <- Train.df$BCR.event
  
  # 3) perform permutations in a loop
  set.seed(123)
  for (i in 1:iterations){
    # permutation
    Train.df$BCR.event <- sample(Train.df$BCR.event, replace=F)
    
    # run model
    model <- E.net.reg(Train.df=Train.df, Test.df=Test.df, weight=weight, clinical=clinical)
    
    # Populate Perf DF
    Perf.DF[1:2,i] <- c(mean(model$train$class==Train.df$BCR.event), mean(model$train$class==corr.labels))
    
    #
    Perf.DF[3:9,i] <- c(model$test$AUC.ED, model$test$AUC.pROC, model$test$AUC.ROCR, model$test$con.mat$overall['Accuracy'], model$test$con.mat$byClass['Specificity'], model$test$con.mat$byClass['Sensitivity'], model$test$con.mat$byClass['Neg Pred Value'])
    print(i)
  }
  
  return(Perf.DF)
}

#############################################
############ TCGA clinical data #############
#############################################

## 1) Get clinical info using TCGAbiolinks
clin.query <- GDCquery(project = 'TCGA-PRAD',
                       data.category = 'Clinical',
                       file.type='xml')
GDCdownload(clin.query)
c_patient <- GDCprepare_clinic(clin.query, clinical.info = 'patient')
c_follow_up <- GDCprepare_clinic(clin.query, clinical.info = 'follow_up')


## 2) Prepare clinical data for TCGA
c_patient <- c_patient %>%
  # select clinicopathological features
  select(bcr_patient_barcode, stage_event_psa, stage_event_tnm_categories, stage_event_gleason_grading, age_at_initial_pathologic_diagnosis, vital_status, race_list, ethnicity, tissue_source_site) %>%
  # delete duplicated entries which are identical
  filter(!duplicated(bcr_patient_barcode)) %>%
  # append sum gleason
  mutate(gleason=case_when(substring(stage_event_gleason_grading,1,1)=='1'~'10',
                           TRUE~substring(stage_event_gleason_grading,1,1))) %>%
  # separate TNM staging to separate categories
  mutate(stage_N=gsub('.*(N[0-9]).*', "\\1", stage_event_tnm_categories)) %>%
  mutate(stage_M=gsub('.*(M[0-9]).*', "\\1", stage_event_tnm_categories)) %>%
  mutate(stage_combT=gsub('M[0-9]|N[0-9]', "\\1", stage_event_tnm_categories)) %>%
  # Get the second/ more extreme pT scoring
  mutate(stage_T=gsub('^T.*(T.*)', '\\1', stage_combT)) #500 10
# If T N or M dont have associated scores, replace them with NA
c_patient$stage_N[grep('N', c_patient$stage_N, invert=T)] <- NA
c_patient$stage_M[grep('M', c_patient$stage_M, invert=T)] <- NA
c_patient$stage_T[grep('T', c_patient$stage_T, invert=T)] <- NA


## 3) Get follow up data ready so it can be merged with c_patient
c_follow_up <- c_follow_up %>% 
  # select relevant coluns
  select(bcr_patient_barcode, lost_follow_up, days_to_last_followup, days_to_death, days_to_first_biochemical_recurrence, year_of_form_completion) %>%
  # BCR experienced? Yes: 2, No: 1
  mutate(BCR = case_when(!is.na(days_to_first_biochemical_recurrence) ~ 2,
                         TRUE ~ 1)) %>%
  # if days to death available, combine it with days to last follow up
  mutate(last_followup = case_when(!is.na(days_to_death) ~ days_to_death,
                                   TRUE ~ days_to_last_followup)) %>%
  # merge days to BCR and last follow-up
  mutate(event_time=case_when(!is.na(days_to_first_biochemical_recurrence)~c_follow_up$days_to_first_biochemical_recurrence,
                              TRUE ~ last_followup)) %>% 
  # remove samples without follow-up info as they cannot be right censored
  filter(!is.na(last_followup)) #504 4


## 4) Some samples have multiple follow-up entries. If they experience BCR, select the first recorded BCR entry. If they do not experience BCR, select the latest follow up entry.
ind <- c_follow_up$bcr_patient_barcode[duplicated(c_follow_up$bcr_patient_barcode)]
unique.df <- c_follow_up[(!c_follow_up$bcr_patient_barcode %in% ind),] #372
## create temporary dfs: with BCR and without BCR event
multiple.df.BCR <- c_follow_up %>% 
  # select non unique samples
  filter(bcr_patient_barcode %in% ind) %>%
  # select samples that experience BCR
  filter(BCR==2) %>% #17
  # sort according to patient ID and event time
  arrange(bcr_patient_barcode, event_time) %>%
  # select first entry: as its the earliest time point when BCR experience recorded
  filter(!duplicated(bcr_patient_barcode)) #12
### w/o BCR
multiple.df.no.BCR <- c_follow_up %>% 
  filter(bcr_patient_barcode %in% ind) %>%
  filter(BCR==1) %>% #115
  arrange(bcr_patient_barcode, desc(event_time)) %>%
  # select first entry: as its the latest follow-up entry
  filter(!duplicated(bcr_patient_barcode)) %>%
  # remove entries with BCR info in multiple.df.no.BCR
  filter(!bcr_patient_barcode %in% multiple.df.BCR$bcr_patient_barcode) #54
## merge unique.df, multiple.df.BCR and multiple.df.no.BCR
c_follow_up <- bind_rows(unique.df, multiple.df.BCR, multiple.df.no.BCR) #438


## 5) Add follow up info to c_patient
c_patient <- c_patient %>%
  # append BCR column
  mutate(BCR=c_follow_up$BCR[match(bcr_patient_barcode, c_follow_up$bcr_patient_barcode)]) %>% 
  # append time to BCR column
  mutate(event_time=c_follow_up$event_time[match(bcr_patient_barcode, c_follow_up$bcr_patient_barcode)])


## 6) Upload corresponding pheno metadata of maturemirnacounts
pheno_manifest <- read.csv('/home/sharmila/Documents/github_projects/ch3_Review_MetaAnalysis/data/TCGA/mir_metadata.csv', check.names=F, stringsAsFactors=F)
# Select and rename useful columns
pheno_manifest <- pheno_manifest[,c(23,32,90,97,93)]
colnames(pheno_manifest) <- c('file_id.isomir', 'file_id.maturemirna', 'submitter_id', 'sample_type','is_ffpe')
rownames(pheno_manifest) <- pheno_manifest$submitter_id 
## Add file_ID to c_patient as later useful to merge with miR counts later
temp <- pheno_manifest %>%
  #  c_patient only contains tumour samples so remove normal samples
  filter(sample_type =='Primary Tumor') %>%
  # select only first 3 terms of submitter id
  mutate(bcr_patient_barcode=str_sub(submitter_id, 1, 12))
c_patient$file.id.maturemiR <- temp$file_id[match(c_patient$bcr_patient_barcode, temp$bcr_patient_barcode)] 
c_patient$file.id.isomir <-temp$file_id.isomir[match(c_patient$bcr_patient_barcode, temp$bcr_patient_barcode)] 


########################################################
########### Data partitioning + normalisation ##########
########################################################

## Reads were annotated using command line tool miraligner and isomiR counts were quantified using R package IsomiRs. In addition to raw counts (object: counts) also gives features information (object: features.iso)
## Load raw counts
# load('raw_isomiR_counts.Rdata')

## 1) Subset samples with follow up info in the counts
recur.samp.df <- c_patient %>%
  # only select samples with followup
  filter(!is.na(BCR)) %>% 
  # select samples present in counts
  filter(file.id.isomir %in% colnames(counts))
intersect(recur.samp.df$file.id.isomir, colnames(counts)) #433(met+4ffpe samples removed?)
counts <- counts[,recur.samp.df$file.id.isomir]


## 2) Split data into train:test set (85:15)
set.seed(1000)
outcome <- as.character(recur.samp.df$BCR) #433
samples.part <- createDataPartition(y=outcome, times=1, p=0.85)
table(outcome[samples.part$Resample1]) #17.14% recurrent samples in train set
table(outcome[setdiff(seq(1, 433, by=1), samples.part$Resample1)]) #16.36% recurrent samples in test set

## 3) Create list and populate with train and test dfs
TrainSet <- lapply(samples.part, function(x) counts[,x])
TestSet <- lapply(samples.part, function(x) counts[,setdiff(seq(1, 433, by=1), x)])


## 4) Filtering: write function to select isomiRs with >= 1 cpm in atleast 80% samples. Filtering on only the train set reduces bias introduced from test set.
# function
filter.df <- function(TrainSet=TrainSet[[1]]){
  
  temp <- cpm(TrainSet)
  # Calculate %samples with counts above 1 for each miRs
  Train_prop_samp_1 <- vector(mode='numeric', length=0)
  
  for (i in 1:nrow(temp)){
    prop <- sum(temp[i,]>=1)/ncol(temp) #1cpm,
    Train_prop_samp_1 <- c(Train_prop_samp_1, prop)
    #print(paste( 'cpm >= 1', i))
  }
  names(Train_prop_samp_1) <- rownames(temp)
  return(Train_prop_samp_1)
}

## 4) Run function
prop_samp_1 <- filter.df(TrainSet = TrainSet[[1]])
# Select isomiRs with RPM >= 1 in 80% samples
sel.mirs1 <- names(prop_samp_1)[which(prop_samp_1>=0.8)] 

## 5) Subset train and test sets
TrainSet.raw <- list()
TrainSet.raw$part1 <- TrainSet[[1]][sel.mirs1,]
TestSet.raw <- list()
TestSet.raw$part1 <- TestSet[[1]][sel.mirs1,]


## 6) Normalisation using TMM
counts.norm <- list()
counts.norm$part1 <- normalise.TMM(train.df=TrainSet.raw$part1, test.df=TestSet.raw$part1)

## 7) Prepare clinical df for merging
# Group gleason score into ISUP grade groups
recur.samp.df <- recur.samp.df %>% 
  # firstly separate composite gleason scores
  mutate(GSS.gp=case_when(str_sub(stage_event_gleason_grading, 2, 3)=='05'~'55',
                          TRUE ~ str_sub(stage_event_gleason_grading, 2, 3))) %>%
  mutate(GSS.gp = unlist(lapply(strsplit(GSS.gp, split=''), function(x) paste0(x[[1]],'+',x[[2]])))) %>%
  # convert composites into groups
  mutate(ISUP.gleason=case_when(GSS.gp=='4+2'|GSS.gp=='3+3'|GSS.gp=='2+4'~ 'gp1',
                                GSS.gp=='3+4'~ 'gp2',
                                GSS.gp=='4+3'~ 'gp3',
                                GSS.gp=='4+4'|GSS.gp=='3+5'|GSS.gp=='5+3'~ 'gp4',
                                GSS.gp=='5+5'|GSS.gp=='5+4'|GSS.gp=='4+5'~ 'gp5',
                                TRUE ~ GSS.gp)) %>% 
  # convert T1, T2, T3, T4 into 2 groups or not? No leave as is
  mutate(pT.grp=gsub('a|b|c', '', stage_T)) %>%
  # Convert psa to numeric
  mutate_at(.funs=as.numeric, vars(stage_event_psa)) %>%
  # BCR event
  mutate(BCR.event = case_when(BCR==2 ~ 'BCR', TRUE ~ 'no_BCR'))
rownames(recur.samp.df) <- recur.samp.df$file.id.isomir

# Asign control level
recur.samp.df$ISUP.gleason <- factor(recur.samp.df$ISUP.gleason, levels=c('gp1','gp2','gp3','gp4','gp5'))
recur.samp.df$pT.grp <- factor(recur.samp.df$pT.grp, levels=c('T1','T2','T3','T4'))
recur.samp.df$BCR.event <- factor(recur.samp.df$BCR.event, levels=c('BCR','no_BCR'))


## 8) Add outcome info + categorical variables in the Train and Test set
recur.samp.df.tmp <- recur.samp.df %>% 
  select(BCR.event,pT.grp,ISUP.gleason,stage_event_psa,file.id.isomir) %>%
  filter(!is.na(pT.grp)) %>%
  filter(!is.na(stage_event_psa)) 
rownames(recur.samp.df.tmp) <- recur.samp.df.tmp$file.id.isomir
intersect(rownames(recur.samp.df.tmp), rownames(counts.norm$part1[[1]]))
intersect(rownames(recur.samp.df.tmp), rownames(counts.norm$part1[[2]]))

# Merge the Train.dfs
counts.norm.Train <- lapply(counts.norm, function(x) merge(x[[1]], recur.samp.df.tmp, by='row.names'))
counts.norm.Train <- lapply(counts.norm.Train,
                            function(x){rownames(x) <- x$Row.names;
                            x$Row.names <- NULL;
                            x$file.id.isomir <- NULL;
                            return(x)})
# Merge the Test.dfs
counts.norm.Test <- lapply(counts.norm, function(x) merge(x[[2]], recur.samp.df.tmp, by='row.names'))
counts.norm.Test <- lapply(counts.norm.Test,
                           function(x){rownames(x) <- x$Row.names;
                           x$Row.names <- NULL;
                           x$file.id.isomir <- NULL;
                           return(x)})
# Recalculate the proportion of BCR samples in each set
lapply(counts.norm.Test, function(x) table(x$BCR.event))


#############################################################
##### Group samples into biologically meaningful groups #####
#############################################################

## 1) Adding additional feature information
# Load miR cluster database and add cluster info
miR.cluster.db <- read.csv('~/miRBase_miRs.csv', stringsAsFactors = F, na.strings = c('')) # Contains information on miRNA primary transcript and mature miRNA transcript. hsa.gff3 was edited to create this file, which was downloaded from mirbase (ftp://mirbase.org/pub/mirbase/21/genomes/).


miR.cluster.db  <- miR.cluster.db %>% 
  filter(miRNA_type == 'miRNA_primary_transcript') %>%
  mutate_at(.vars='ID', .funs=gsub, pattern='.*\\=', replacement='') %>%
  mutate_at(.vars='Alias', .funs=gsub, pattern='.*\\=', replacement='') %>%
  mutate_at(.vars='Name', .funs=gsub, pattern='.*\\=hsa-', replacement='') %>%
  mutate_at(.vars='Derives_from', .funs=gsub, pattern='.*\\=', replacement='') %>% 
  mutate_at(.vars='Name', .funs=gsub, pattern='r', replacement='R')
# Assign miRs to a cluster? Cluster if inter-miR distance < 1kb (TargetScan)
temp <- sapply(miR.cluster.db$ch_start, function(x) ifelse(abs(x-miR.cluster.db$ch_start) < 1000, miR.cluster.db$Name, NA))
miR.cluster.db$cluster <- apply(temp, 2, function(x) paste(x[!is.na(x) & x != "No"], collapse = ";"))

## Adding additional feature information
features.iso <- features.iso %>%
  # append seed 7mer
  mutate(seed.7mer=substr(seq, 2, 8)) %>%
  # append seed 6mer
  mutate(seed.6mer=substr(seq, 2, 7)) %>%
  # append cluster info
  mutate(pre.miR=gsub('-[0-9]p$', '', mir)) %>%
  mutate(cluster=miR.cluster.db$cluster[match(pre.miR, miR.cluster.db$Name)]) %>%
  # append isotype info
  mutate(isotype=case_when((t5=='0' & t3=='0' & mism=='0' & add=='0') ~ 'canonical',
                           (t5!='0' & t3=='0' & mism=='0' & add=='0') ~ 'iso.5p',
                           (t5=='0' & t3!='0' & mism=='0' & add=='0') ~ 'iso.3p',
                           (t5=='0' & t3=='0' & mism!='0' & add=='0') ~ 'mism',
                           (t5=='0' & t3=='0' & mism=='0' & add!='0') ~ 'add.3p',
                           (t5!='0' & t3!='0' & mism=='0' & add=='0') ~ 'iso.5p+iso.3p',
                           (t5!='0' & t3=='0' & mism!='0' & add=='0') ~ 'iso.5p+mism',
                           (t5!='0' & t3=='0' & mism=='0' & add!='0') ~ 'iso.5p+add.3p',
                           (t5=='0' & t3!='0' & mism!='0' & add=='0') ~ 'iso.3p+mism',
                           (t5=='0' & t3!='0' & mism=='0' & add!='0') ~ 'iso.3p+add.3p',
                           (t5=='0' & t3=='0' & mism!='0' & add!='0') ~ 'mism+add.3p',
                           (t5!='0' & t3!='0' & mism!='0' & add=='0') ~ 'iso.5p+iso.3p+mism',
                           (t5!='0' & t3!='0' & mism=='0' & add!='0') ~ 'iso.5p+iso.3p+add.3p',
                           (t5!='0' & t3=='0' & mism!='0' & add!='0') ~ 'iso.5p+mism+add.3p',
                           (t5=='0' & t3!='0' & mism!='0' & add!='0') ~ 'iso.3p+mism+add.3p',
                           (t5!='0' & t3!='0' & mism!='0' & add!='0') ~ 'all'))

## Subset to include only features that passed filtering criteria
features.part1 <- features.iso[features.iso$isomir %in% colnames(counts.norm$part1[[1]]),]
## Also add isotypes of miR-148a-3p
features.part1$isotype_148a <- ifelse(features.part1$mir == 'miR-148a-3p',
                                      features.part1$isotype, 'not_miR-148a')

## 3) Create module sets and assign features into the sets
# i) Parent miR
mod.parent.mir <- lapply(unique(features.part1$mir), function(x) features.part1$isomir[which(features.part1$mir %in% x)]); names(mod.parent.mir) <- unique(features.part1$mir)
# ii) Cluster
mod.cluster <- lapply(unique(features.part1$cluster), function(x) features.part1$isomir[which(features.part1$cluster %in% x)])
names(mod.cluster) <- unique(features.part1$cluster)
mod.cluster <- mod.cluster[-which(is.na(names(mod.cluster)))]
# iii) isotype
mod.isotype <- lapply(unique(features.part1$isotype), function(x) features.part1$isomir[which(features.part1$isotype %in% x)])
names(mod.isotype) <- unique(features.part1$isotype)
# iv) t3 size
mod.t3size <- lapply(unique(features.part1$t3.size), function(x) features.part1$isomir[which(features.part1$t3.size %in% x)])
names(mod.t3size) <- unique(features.part1$t3.size)
# v) t5 size
mod.t5size <- lapply(unique(features.part1$t5.size), function(x) features.part1$isomir[which(features.part1$t5.size %in% x)])
names(mod.t5size) <- unique(features.part1$t5.size)
# vi) 7mer-m8
mod.seed7mer <- lapply(unique(features.part1$seed.7mer), function(x) features.part1$isomir[which(features.part1$seed.7mer %in% x)])
names(mod.seed7mer) <- unique(features.part1$seed.7mer)
# vii) 6mer
mod.seed6mer <- lapply(unique(features.part1$seed.6mer), function(x) features.part1$isomir[which(features.part1$seed.6mer %in% x)])
names(mod.seed6mer) <- unique(features.part1$seed.6mer)
# viii) miR-148a-isoytpes
mod.148a.isotype <- lapply(unique(features.part1$isotype_148a), function(x) features.part1$isomir[which(features.part1$isotype_148a %in% x)]);
names(mod.148a.isotype) <- unique(features.part1$isotype_148a)
# Remove set which are not miR-148a-3p
mod.148a.isotype <- mod.148a.isotype[-1]


## 3) Transpose counts norm
counts.norm.t <- list()
counts.norm.t$part1 <- lapply(counts.norm$part1, function(x) data.frame(t(x), check.names=F))


#########################################################################
######################## mean signature/ modules ########################
#########################################################################

## 4) Get the mean of the isomiRs in each modules sets
mean.mod <- list()
# parent miR
mean.mod$parent.train <- sapply(mod.parent.mir,scoreSignature,exprs=counts.norm.t$part1[[1]],method="mean")
mean.mod$parent.test <- sapply(mod.parent.mir,scoreSignature,exprs=counts.norm.t$part1[[2]],method="mean")
# cluster
mean.mod$cluster.train <- sapply(mod.cluster,scoreSignature,exprs=counts.norm.t$part1[[1]],method="mean")
mean.mod$cluster.test <- sapply(mod.cluster,scoreSignature,exprs=counts.norm.t$part1[[2]],method="mean")
# isotype
mean.mod$isotype.train <- sapply(mod.isotype,scoreSignature,exprs=counts.norm.t$part1[[1]],method="mean")
mean.mod$isotype.test <- sapply(mod.isotype,scoreSignature,exprs=counts.norm.t$part1[[2]],method="mean")
# t3size
mean.mod$t3size.train <- sapply(mod.t3size,scoreSignature,exprs=counts.norm.t$part1[[1]],method="mean")
mean.mod$t3size.test <- sapply(mod.t3size,scoreSignature,exprs=counts.norm.t$part1[[2]],method="mean")
# t5size
mean.mod$t5size.train <- sapply(mod.t5size,scoreSignature,exprs=counts.norm.t$part1[[1]],method="mean")
mean.mod$t5size.test <- sapply(mod.t5size,scoreSignature,exprs=counts.norm.t$part1[[2]],method="mean")
# seed 7merm8
mean.mod$seed7mer.train <- sapply(mod.seed7mer,scoreSignature,exprs=counts.norm.t$part1[[1]],method="mean")
mean.mod$seed7mer.test <- sapply(mod.seed7mer,scoreSignature,exprs=counts.norm.t$part1[[2]],method="mean")
# seed 6mer
mean.mod$seed6mer.train <- sapply(mod.seed6mer,scoreSignature,exprs=counts.norm.t$part1[[1]],method="mean")
mean.mod$seed6mer.test <- sapply(mod.seed6mer,scoreSignature,exprs=counts.norm.t$part1[[2]],method="mean")
# miR-148a isotype
mean.mod$isotype.148a.train <- sapply(mod.148a.isotype,scoreSignature,exprs=counts.norm.t$part1[[1]],method="mean")
mean.mod$isotype.148a.test <- sapply(mod.148a.isotype,scoreSignature,exprs=counts.norm.t$part1[[2]],method="mean")

## 5) Merge the mean.mod and clinical df
mean.mod.clin <- lapply(mean.mod, function(x) merge(x, recur.samp.df.tmp, by='row.names'))
mean.mod.clin <- lapply(mean.mod.clin,
                        function(x){rownames(x) <- x$Row.names;
                        x$Row.names <- NULL;
                        x$file.id.isomir <- NULL;
                        return(x)})


######################################
########### Perform ENR ##############
######################################


######## WITHOUT WEIGHTS 
## 1) EN model without weight for main dataset (partition 1)
EN.res <- list()
EN.res$CVonly <- E.net.reg(Train.df=counts.norm.Train$part1, Test.df=counts.norm.Test$part1, weight='N', clinical=c('only'))
# coef(EN.res$CVonly$model$finalModel, EN.res$CVonly$model$bestTune$lambda)
EN.res$indIso.only <- E.net.reg(Train.df=counts.norm.Train$part1, Test.df=counts.norm.Test$part1, weight='N', clinical=c('N'))
EN.res$indIso <- E.net.reg(Train.df=counts.norm.Train$part1, Test.df=counts.norm.Test$part1, weight='N', clinical=c('Y'))


####### WITH WEIGHTS
## 2) EN model with weight for main dataset (partition 1)
EN.res.wt.mean <- list()
# Individual isomiRs models
EN.res.wt.mean$CVonly <- E.net.reg(Train.df=counts.norm.Train$part1, Test.df=counts.norm.Test$part1, weight='Y', clinical=c('only'))
EN.res.wt.mean$indIso.only <- E.net.reg(Train.df=counts.norm.Train$part1, Test.df=counts.norm.Test$part1, weight='Y', clinical=c('N'))
EN.res.wt.mean$indIso <- E.net.reg(Train.df=counts.norm.Train$part1, Test.df=counts.norm.Test$part1, weight='Y', clinical=c('Y'))
## Model with only 148a isomiRs
train.temp.148a <- counts.norm.Train$part1 %>%
  select(c(matches('148a-3p'),'BCR.event','pT.grp','ISUP.gleason','stage_event_psa'))
test.temp.148a <- counts.norm.Test$part1 %>%
  select(c(matches('148a-3p'),'BCR.event','pT.grp','ISUP.gleason','stage_event_psa'))
EN.res.wt$iso_148a <- E.net.reg(Train.df=train.temp.148a, Test.df=test.temp.148a, weight='Y', clinical=c('Y'))
## Model with only 582-5p isomiRs
train.temp.582 <- counts.norm.Train$part1 %>%
  select(c(matches('582-5p'),'BCR.event','pT.grp','ISUP.gleason','stage_event_psa'))
test.temp.582 <- counts.norm.Test$part1 %>%
  select(c(matches('582-5p'),'BCR.event','pT.grp','ISUP.gleason','stage_event_psa'))
EN.res.wt$iso_582 <- E.net.reg(Train.df=train.temp.582, Test.df=test.temp.582, weight='Y', clinical=c('Y'))
## Model with only 148a + 582-5p isomiRs
train.temp.148a_582 <- counts.norm.Train$part1 %>%
  select(c(matches('148a-3p|582-5p'),'BCR.event','pT.grp','ISUP.gleason','stage_event_psa'))
test.temp.148a_582 <- counts.norm.Test$part1 %>%
  select(c(matches('148a-3p|582-5p'),'BCR.event','pT.grp','ISUP.gleason','stage_event_psa'))
EN.res.wt$iso_148a_582 <- E.net.reg(Train.df=train.temp.148a_582, Test.df=test.temp.148a_582, weight='Y', clinical=c('Y'))
# Signature isomiRs models
EN.res.wt.mean$parentmiR <- E.net.reg(Train.df=mean.mod.clin$parent.train, Test.df=mean.mod.clin$parent.test, weight='Y', clinical=c('Y'))
EN.res.wt.mean$cluster <- E.net.reg(Train.df=mean.mod.clin$cluster.train, Test.df=mean.mod.clin$cluster.test, weight='Y', clinical=c('Y'))
EN.res.wt.mean$isotype <- E.net.reg(Train.df=mean.mod.clin$isotype.train, Test.df=mean.mod.clin$isotype.test, weight='Y', clinical=c('Y'))
EN.res.wt.mean$t3size <- E.net.reg(Train.df=mean.mod.clin$t3size.train, Test.df=mean.mod.clin$t3size.test, weight='Y', clinical=c('Y'))
EN.res.wt.mean$t5size <- E.net.reg(Train.df=mean.mod.clin$t5size.train, Test.df=mean.mod.clin$t5size.test, weight='Y', clinical=c('Y'))
EN.res.wt.mean$seed7mer <- E.net.reg(Train.df=mean.mod.clin$seed7mer.train, Test.df=mean.mod.clin$seed7mer.test, weight='Y', clinical=c('Y'))
EN.res.wt.mean$seed6mer <- E.net.reg(Train.df=mean.mod.clin$seed6mer.train, Test.df=mean.mod.clin$seed6mer.test, weight='Y', clinical=c('Y'))
EN.res.wt.mean$isotype.148a <- E.net.reg(Train.df=mean.mod.clin$isotype.148a.train, Test.df=mean.mod.clin$isotype.148a.test, weight='Y', clinical=c('Y'))


## 3) Extract the performance metrics in a df
EN.wt.df.mean <- data.frame(models=names(EN.res.wt.mean), stringsAsFactors=F) 
EN.wt.df.mean <-  EN.wt.df.mean %>%
  mutate(AUC.ED=signif(sapply(EN.res.wt.mean, function(x) x$test$AUC.ED), 3)) %>%
  mutate(AUC.pROC=signif(sapply(EN.res.wt.mean, function(x) x$test$AUC.pROC), 3)) %>%
  #mutate(AUC.ROCR=signif(sapply(EN.res.wt, function(x) x$test$AUC.ROCR), 3)) %>%
  mutate(p.value=NA) %>%
  mutate(accuracy=signif(as.numeric(sapply(EN.res.wt.mean, function(x) 
    coords(roc=x$test$ROC.pROC, x=0.6, input='specificity', ret='precision', transpose=F))), 3)) %>%
  mutate(sensitivity=signif(as.numeric(sapply(EN.res.wt.mean, function(x) 
    coords(roc=x$test$ROC.pROC, x=0.6, input='specificity', ret='sensitivity', transpose=F))), 3)) %>%
  mutate(specificity=signif(as.numeric(sapply(EN.res.wt.mean, function(x) 
    coords(roc=x$test$ROC.pROC, x=0.6, input='specificity', ret='specificity', transpose=F))), 3)) %>% 
  mutate(precision=signif(as.numeric(sapply(EN.res.wt.mean, function(x) 
    coords(roc=x$test$ROC.pROC, x=0.6, input='specificity', ret='precision', transpose=F))), 3)) %>% 
  mutate(NPV=signif(as.numeric(sapply(EN.res.wt.mean, function(x) 
    coords(roc=x$test$ROC.pROC, x=0.6, input='specificity', ret='npv', transpose=F))), 3))



## 4) Plot ROC curves
# colours
n=length(EN.res.wt.mean); qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

## Plot ROC curves pROC
pdf('Documents/github_projects_with_data/ch4_PrognosticModels/results/EN_wt_pAUC.pdf', height=5, width=6)

for (i in 1:length(EN.res.wt.mean)){
  plot(EN.res.wt.mean[[i]]$test$ROC.pROC, print.thres='no',  print.auc=T, print.auc.y=0.1, print.auc.x=0.3, print.auc.cex=1.4, type='S', col=col_vector[6], lwd=2, xlim=c(1,0), legacy.axes=T, cex.main=1.4,cex.lab=1.4,cex.axis=1.4) #  main=names(EN.res.wt.mean)[[i]],
}
dev.off()


################################################################
######################## Permutate model #######################
################################################################

## 1) Run permutation models
Perf.DF.CVonly <- null.perm(Train.df=counts.norm.Train$part1, Test.df=counts.norm.Test$part1, weight='Y', clinical='only', iterations=1000)
Perf.DF.indIso.only <- null.perm(Train.df=counts.norm.Train$part1, Test.df=counts.norm.Test$part1, weight='Y', clinical='N', iterations=1000)
Perf.DF.indIso <- null.perm(Train.df=counts.norm.Train$part1, Test.df=counts.norm.Test$part1, weight='Y', clinical=c('Y'), iterations=1000)
## Model with only 148a isomiRs
train.temp <- counts.norm.Train$part1 %>%
  select(c(matches('148a-3p'),'BCR.event','pT.grp','ISUP.gleason','stage_event_psa'))
test.temp <- counts.norm.Test$part1 %>%
  select(c(matches('148a-3p'),'BCR.event','pT.grp','ISUP.gleason','stage_event_psa'))
Perf.DF.iso148a <- null.perm(Train.df=train.temp, Test.df=test.temp, weight='Y', clinical=c('Y'), iterations=1000)
## Model with only 582-5p isomiRs
train.temp <- counts.norm.Train$part1 %>%
  select(c(matches('582-5p'),'BCR.event','pT.grp','ISUP.gleason','stage_event_psa'))
test.temp <- counts.norm.Test$part1 %>%
  select(c(matches('582-5p'),'BCR.event','pT.grp','ISUP.gleason','stage_event_psa'))
Perf.DF.iso582 <- null.perm(Train.df=train.temp, Test.df=test.temp, weight='Y', clinical=c('Y'), iterations=1000)
## Model with only both 148 582-5p isomiRs
train.temp <- counts.norm.Train$part1 %>%
  select(c(matches('148a-3p|582-5p'),'BCR.event','pT.grp','ISUP.gleason','stage_event_psa'))
test.temp <- counts.norm.Test$part1 %>%
  select(c(matches('148a-3p|582-5p'),'BCR.event','pT.grp','ISUP.gleason','stage_event_psa'))
Perf.DF.iso148a_582 <- null.perm(Train.df=train.temp, Test.df=test.temp, weight='Y', clinical=c('Y'), iterations=1000)
Perf.DF.parentmiR <- null.perm(Train.df=segal.mod.clin$parent.train, Test.df=segal.mod.clin$parent.test, weight='Y', clinical=c('Y'), iterations=1000)
Perf.DF.cluster <- null.perm(Train.df=segal.mod.clin$cluster.train, Test.df=segal.mod.clin$cluster.test, weight='Y', clinical=c('Y'), iterations=1000)
Perf.DF.isotype <- null.perm(Train.df=segal.mod.clin$isotype.train, Test.df=segal.mod.clin$isotype.test, weight='Y', clinical=c('Y'), iterations=1000)
Perf.DF.t3size <- null.perm(Train.df=segal.mod.clin$t3size.train, Test.df=segal.mod.clin$t3size.test, weight='Y', clinical=c('Y'), iterations=1000)
Perf.DF.t5size <- null.perm(Train.df=segal.mod.clin$t5size.train, Test.df=segal.mod.clin$t5size.test, weight='Y', clinical=c('Y'), iterations=1000)
Perf.DF.seed7mer <- null.perm(Train.df=segal.mod.clin$seed7mer.train, Test.df=segal.mod.clin$seed7mer.test, weight='Y', clinical=c('Y'), iterations=1000)
Perf.DF.seed6mer <- null.perm(Train.df=segal.mod.clin$seed6mer.train, Test.df=segal.mod.clin$seed6mer.test, weight='Y', clinical=c('Y'), iterations=1000)


## 2) Save the AUC in dfs
perm.AUC.df <- data.frame(CVs=as.numeric(Perf.DF.CVonly['test:AUC.pROC',]), indIso.only=as.numeric(Perf.DF.indIso.only['test:AUC.pROC',]), indIso=as.numeric(Perf.DF.indIso['test:AUC.pROC',]), iso148=as.numeric(Perf.DF.iso148a['test:AUC.pROC',]), iso582=as.numeric(Perf.DF.iso582['test:AUC.pROC',]), iso148_582=as.numeric(Perf.DF.iso148a_582['test:AUC.pROC',]), parentmiR=as.numeric(Perf.DF.parentmiR['test:AUC.pROC',]), cluster=as.numeric(Perf.DF.cluster['test:AUC.pROC',]), isotype=as.numeric(Perf.DF.isotype['test:AUC.pROC',]), t3size=as.numeric(Perf.DF.t3size['test:AUC.pROC',]), t5size=as.numeric(Perf.DF.t5size['test:AUC.pROC',]), seed7=as.numeric(Perf.DF.seed7mer['test:AUC.pROC',]), seed6=as.numeric(Perf.DF.seed6mer['test:AUC.pROC',]), isotype.148a=as.numeric(Perf.DF.148aisotype['test:AUC.pROC',]))

actual.AUC.df <- c(EN.res.wt.mean$CVonly$test$AUC.pROC, EN.res.wt.mean$indIso.only$test$AUC.pROC, EN.res.wt.mean$indIso$test$AUC.pROC, EN.res.wt.mean$iso_148a$test$AUC.pROC, EN.res.wt.mean$iso_582$test$AUC.pROC, EN.res.wt.mean$iso_148a_582$test$AUC.pROC, EN.res.wt.mean$parentmiR$test$AUC.pROC, EN.res.wt.mean$cluster$test$AUC.pROC, EN.res.wt.mean$isotype$test$AUC.pROC, EN.res.wt.mean$t3size$test$AUC.pROC, EN.res.wt.mean$t5size$test$AUC.pROC, EN.res.wt.mean$seed7mer$test$AUC.pROC, EN.res.wt.mean$seed6mer$test$AUC.pROC, EN.res.wt.mean$isotype.148a$test$AUC.pROC)

## 3) Compute p-value for permutation test
p.val <-vector()
for(i in 1:ncol(perm.AUC.df)){
  iterations <- 1000
  p.val <- round(c(p.val, (1+sum(perm.AUC.df[,i]>=actual.AUC.df[i]))/(1+iterations)), 3) # significance test for the AUC we got from actual test is no different than expected by chance
}
p.val
