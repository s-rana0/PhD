library(Biobase)
library(GEOquery)
library(dplyr)
library(tidyr)
library(edgeR)
library(stringr)
library(TCGAbiolinks)
library(miRBaseConverter)

setwd('~/Documents/github_projects/ch3_Review_MetaAnalysis/')
source('data/TCGA/mergeTCGAmaturemirna.R') # Code to get miR counts out of GDC

#######################################
############## Functions ##############
#######################################

# 1) Function to get the median value if a miR has multiple probes
calc.df.median <- function(df=Suer_gse88958_exprs, miR='miR-27a'){
  df <- df[, grep(miR, colnames(df))]
  df.median <- apply(df, 1, median)
  return(df.median)
}

#####################################################
################# Obtain datasets ###################
#####################################################

##########################
########## TCGA ##########
##########################

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
  select(bcr_patient_barcode, stage_event_psa, stage_event_tnm_categories, stage_event_gleason_grading, age_at_initial_pathologic_diagnosis) %>%
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


## 6) Now get TCGA expression directly from GDC data portal (requires level 1 data access)
mir_manifest <- read.table('data/TCGA/gdc_download_20170213_163009/MANIFEST.txt', sep='\t', header=T, stringsAsFactors=F)
mirbasehuman <- read.csv('data/TCGA/mirbasehuman.csv', stringsAsFactors = F)
maturemirnacounts <- mergeTCGAmaturemirna(path='data/TCGA/gdc_download_20170213_163009/', mirfileIDs=mir_manifest$id, mirbasenames=mirbasehuman)
# Set rownames
rownames(maturemirnacounts) <- maturemirnacounts[,1]
maturemirnacounts <- maturemirnacounts[,-1]
# FYI: maturemirnacounts original datafile downloaded from GDC is same as miRNA isomiR counts downloaded from TCGAbiolinks 


## 7) Upload corresponding pheno metadata of maturemirnacounts
pheno_manifest <- read.csv('/home/sharmila/Documents/github_projects/ch3_Review_MetaAnalysis/data/TCGA/mir_metadata.csv', check.names=F, stringsAsFactors=F)
# Select and rename useful columns
pheno_manifest <- pheno_manifest[,c(32,90,97,93)]
colnames(pheno_manifest) <- c('file_id', 'submitter_id', 'sample_type','is_ffpe')
rownames(pheno_manifest) <- pheno_manifest$submitter_id 
# rename maturemirnacounts colnames with submitter_id
setdiff(colnames(maturemirnacounts), pheno_manifest$file_id)
colnames(maturemirnacounts) <- pheno_manifest$submitter_id[match(colnames(maturemirnacounts), pheno_manifest$file_id)]
maturemirnacounts <- maturemirnacounts[,order(colnames(maturemirnacounts))]
## Add file_ID to c_patient as later useful to merge with maturemirnacounts
temp <- pheno_manifest %>%
  #  c_patient only contains tumour samples so remove normal samples
  filter(sample_type =='Primary Tumor') %>%
  # select only first 3 terms of submitter id
  mutate(bcr_patient_barcode=str_sub(submitter_id, 1, 12))
c_patient$file_id <-temp$file_id[match(c_patient$bcr_patient_barcode, temp$bcr_patient_barcode)] 
# subset c_patient to include entries for samples with follow up info only
c_patient_fp <- c_patient %>% filter(!is.na(c_patient$BCR)) #438


## 8) Remove unwanted samples
# Remove metastatic sample from expression info as we focus on priamry samples only
met <- rownames(pheno_manifest)[which(pheno_manifest$sample_type == 'Metastatic')]
maturemirnacounts <- maturemirnacounts[,!colnames(maturemirnacounts) %in% met]
# remove ffpe samples http://gdac.broadinstitute.org/runs/stddata__2016_01_28/samples_report/UVM_FFPE_Cases.html, https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/163
pheno_manifest <- pheno_manifest %>% filter(is_ffpe !='true')
maturemirnacounts <- maturemirnacounts[,colnames(maturemirnacounts) %in% pheno_manifest$submitter_id]

## 9) Now, filter miRs to include samples with cpm counts >= 1cpm in 80% of the samples
# Normalise the raw counts with TNM 
dgList <- DGEList(counts= maturemirnacounts, genes=rownames(maturemirnacounts))
dgList_norm <- calcNormFactors(dgList, method='TMM') 
norm <- cpm(dgList_norm, normalized.lib.sizes = T, log=F)
norm <- t(norm) #Transpose so miRs as columns
# Calculate % samples with counts above >=1cpm for each miR
prop_samp <- vector(mode = 'numeric', length = 0) 
# start loop to populate prop_sample
for (i in 1:ncol(norm)){
  prop <- sum(norm[,i]>=1)/nrow(norm)
  prop_samp <- c(prop_samp, prop)
}
# 80% samples with miRNA counts >= 1cpm
mirs <- colnames(norm)[which(prop_samp>=0.8)] #328
# Get log2 normalised data and include these miRs 
norm <- cpm(dgList_norm, normalized.lib.sizes = T, log=T)
norm <- norm[rownames(norm) %in% mirs,] #Transpose so miRs as columns


## 10) subset and remove normal samples for Cox PH
miR_normal <- pheno_manifest$submitter_id[pheno_manifest$sample_type=='Solid Tissue Normal']
surv_counts <- norm[,!colnames(norm) %in% miR_normal] #328 494
colnames(surv_counts) <- str_sub(colnames(surv_counts),1,12)
# Remove samples without follow-up data for Cox PH
surv_counts <- surv_counts[ ,colnames(surv_counts) %in% c_patient_fp$bcr_patient_barcode] #328 433


## 11) Calculate z-score scale for the 328 miRs
surv_counts <- (surv_counts- apply(surv_counts, 1, median, na.rm=T))/apply(surv_counts, 1, sd, na.rm=T)
surv_counts <- data.frame(t(surv_counts)) #433 328
colnames(surv_counts) <- gsub('\\.', '-', colnames(surv_counts))
rownames(surv_counts) <- gsub('\\.', '-', rownames(surv_counts))


## 12) Merge pheno + exprs datasets for Cox PH
intersect(c_patient_fp$bcr_patient_barcode, rownames(surv_counts))
# merge 
surv_counts$bcr_patient_barcode <- as.vector(rownames(surv_counts))
TCGA.cox.df <- inner_join(surv_counts, c_patient_fp, by="bcr_patient_barcode", keep=T) #433 341



###########################
########## MSKCC ##########
###########################

## 1) Load data from GEO
Taylor_gse21036 <- getGEO("GSE21036", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(Taylor_gse21036) > 1) idx <- grep("GPL8227", attr(Taylor_gse21036, "names")) else idx <- 1
Taylor_gse21036 <- Taylor_gse21036[[idx]]

# If not log2 transformed, log2 it
ex <- exprs(Taylor_gse21036)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(Taylor_gse21036) <- log2(ex) }

## 2) Pheno, features and expression datasets
Taylor_gse21036_feat <- fData(Taylor_gse21036) 
Taylor_gse21036_exprs <- exprs(Taylor_gse21036)
Taylor_gse21036_pheno <- pData(Taylor_gse21036) 

## 3) Editing pheno data
Taylor_gse21036_pheno$characteristics_ch1 <- substr(Taylor_gse21036_pheno$characteristics_ch1, 12,18)
rownames(Taylor_gse21036_pheno) <- Taylor_gse21036_pheno$characteristics_ch1
## As this phenodata does not contain enough info, get other clinical info directly from https://cbio.mskcc.org/cancergenomics/prostate/data/
clinicaldata <- read.csv('data/MSKCC_PCa_Clinical_Annotation.csv', stringsAsFactors = F) 

## 4) Merge clinicaldata + Taylor_gse21036_pheno
clinicaldata <- clinicaldata %>% 
  # add geo_accession
  mutate(GEO_accession=Taylor_gse21036_pheno$geo_accession[match(Sample.ID, row.names(Taylor_gse21036_pheno))]) %>%
  # Only select samples with expression info
  filter(GEO_accession %in% colnames(Taylor_gse21036_exprs)) %>%
  # floor age
  mutate_at(.funs=floor, vars(DxAge)) %>%
  # Convert time to BCR From month to days
  mutate(event_time = floor(BCR_FreeTime*30.4167)) %>%
  # Convert BCR occurence to binary, 1: no BCR, 2: BCR
  mutate(BCR = case_when(BCR_Event == 'NO' ~ 1,
                         TRUE ~ 2)) %>%
  mutate_at(.funs=as.numeric, vars(PathGGS)) %>%
  # subset useful columns
  select(Sample.ID, GEO_accession, Type, DxAge, PreTxPSA, BCR, event_time, PreDxBxPSA, PathStage, PathGGS) %>%
  # select primary samples only
  filter(Type=='PRIMARY')


## 5) Subset exprs to include primary samples only 
Taylor.exprs <- Taylor_gse21036_exprs[, clinicaldata$GEO_accession]
# Select only HSA samples
grep('hsa', rownames(Taylor.exprs), ignore.case=T, value=T, invert=T)
# Calculate z-score 
Taylor.exprs <- (Taylor.exprs- apply(Taylor.exprs, 1, median, na.rm=T))/apply(Taylor.exprs, 1, sd, na.rm=T)
Taylor.exprs <- t(Taylor.exprs)


## 6) Merge pheno + exprs datasets for Cox PH
clinicaldata$GEO_accession == row.names(Taylor.exprs)
Taylor.cox.df <- cbind(Taylor.exprs, clinicaldata) 

  

##########################
####### Leite 2015 #######
##########################

## 1) Load data from GEO
Leite_gse46738 <- getGEO("GSE46738", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(Leite_gse46738) > 1) idx <- grep("GPL8786", attr(Leite_gse46738, "names")) else idx <- 1
Leite_gse46738 <- Leite_gse46738[[idx]]

# If not log2 transformed, log2 it
ex <- exprs(Leite_gse46738)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(Leite_gse46738) <- log2(ex) }


## 2) Pheno, features and expression datasets
Leite_gse46738_feat <- fData(Leite_gse46738) # 7815 11
Leite_gse46738_exprs <- exprs(Leite_gse46738) # 7815 57
Leite_gse46738_pheno <- pData(Leite_gse46738) # 57 43
colnames(Leite_gse46738_pheno) <- gsub(' ', '\\.', gsub('\\:ch1', '', colnames(Leite_gse46738_pheno)))

## 3) Getting pheno data ready 
Leite_gse46738_pheno  <- Leite_gse46738_pheno %>% 
  # remove control samples
  filter(tumor.stage!='N/A') %>%
  # case ID
  mutate(cases = as.numeric(gsub('([0-9])_.*', '\\1', title)))
# Change colnames for PSA
colnames(Leite_gse46738_pheno)[grep('psa.level.\\(ng\\/ml\\)', colnames(Leite_gse46738_pheno))] <- 'Pre.surgery.PSA'
Leite_gse46738_pheno$Pre.surgery.PSA <- as.numeric(Leite_gse46738_pheno$Pre.surgery.PSA)


## 4) Leite_gse46738_pheno does not contain event (BCR) + time to event info.  Upload info sent by original author
Leite_pheno2 <- read.csv('data/GSE46738_PSA_data.csv', stringsAsFactors = F, check.names = T, na.strings = c('', ' '))
Leite_pheno2  <- Leite_pheno2 %>%
  #calculate time to last follow up/ BCR (days)
  mutate(last.followup = as.Date(as.character(last.visit), format="%d/%m/%Y") - as.Date(as.character(date), format="%d/%m/%Y")) %>%
  mutate(time.to.BCR = as.Date(as.character(PSA.recurrence.date), format="%d/%m/%Y")- as.Date(as.character(date), format="%d/%m/%Y"))  %>%
  # Calculate time to event
  mutate(event_time = case_when(is.na(time.to.BCR) ~ last.followup,
                                TRUE ~ time.to.BCR)) %>%
  # Code recurrence binary, 1: no BCR, 2: BCR
  mutate(BCR = case_when(is.na(PSA.recurrence.date)~1,
                         TRUE~2))


## 5) Merge Leite_pheno2 and Leite_gse46738_pheno
Leite_gse46738_pheno <- left_join(Leite_gse46738_pheno, Leite_pheno2, by=c('cases', 'Pre.surgery.PSA'))
# Select columns of interest + clean pheno data
Leite_gse46738_pheno <- Leite_gse46738_pheno %>%
  select(geo_accession, cases, age, gleason.score, Pre.surgery.PSA, tumor.stage, last.PSA, PSA.level.at.recurrence, BCR, event_time) %>% 
  # clean gleason and pT columns
  mutate(gleason.score = as.numeric(gleason.score)) %>%
  mutate(tumor.stage = gsub('^p|(.*)N.', '\\1', tumor.stage)) %>% 
  # filter samples without follow up info
  filter(!is.na(event_time))
rownames(Leite_gse46738_pheno) <- Leite_gse46738_pheno$geo_accession


## 6) Get expression data ready
# Subset expression profile of human miRs only
ind <- grep('hsa', rownames(Leite_gse46738_exprs))
Leite_gse46738_exprs <- Leite_gse46738_exprs[ind,]
# Replace star in miR nomenclature with * symbol
rownames(Leite_gse46738_exprs) <- gsub('hsa-', '', gsub('_st$', '', rownames(Leite_gse46738_exprs)))
rownames(Leite_gse46738_exprs) <- gsub('-star', '\\*', rownames(Leite_gse46738_exprs))
## Remove control samples 
Leite_gse46738_exprs <- Leite_gse46738_exprs[,Leite_gse46738_pheno$geo_accession]
## Scale/ z-normalise data here
Leite_gse46738_exprs <- (Leite_gse46738_exprs- apply(Leite_gse46738_exprs, 1, median, na.rm=T))/apply(Leite_gse46738_exprs, 1, sd, na.rm=T)
Leite_gse46738_exprs <- t(Leite_gse46738_exprs)


## 7) Merge pheno + exprs datasets for Cox PH
rownames(Leite_gse46738_pheno) == rownames(Leite_gse46738_exprs)
Leite.cox.df <- cbind(Leite_gse46738_exprs, Leite_gse46738_pheno)



############################
########### Suer ###########
############################

## 1) Load data from GEO
Suer_gse88958 <- getGEO("GSE88958", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(Suer_gse88958) > 1) idx <- grep("GPL14767", attr(Suer_gse88958, "names")) else idx <- 1
Suer_gse88958 <- Suer_gse88958[[idx]]

# log2 transform if not already so
ex <- exprs(Suer_gse88958)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(Suer_gse88958) <- log2(ex) }

## 2) Pheno, features and expression datasets
Suer_gse88958_feat <- fData(Suer_gse88958) # 15739 40
Suer_gse88958_exprs <- exprs(Suer_gse88958) # 15739 40
Suer_gse88958_pheno <- pData(Suer_gse88958) # 57 43
Suer_gse88958_pheno$sample_ID <- gsub('.*#', '', Suer_gse88958_pheno$title)
Suer_gse88958_pheno$sample_ID <- ifelse(Suer_gse88958_pheno$`recurrence status:ch1` == 'recurrent', paste0('R', Suer_gse88958_pheno$sample_ID), paste0('NR', Suer_gse88958_pheno$sample_ID))


## 3) Suer_gse88958_pheno does not contain time to event info or clinical info. Upload info sent by original author and use that dataset
Suer_pheno <- read.csv('data/GSE88958_Suer.csv', stringsAsFactors = F)
setdiff(Suer_gse88958_pheno$sample_ID, Suer_pheno$LAB.ID) # "NR31_2" "R42" "R43"
setdiff(Suer_pheno$LAB.ID, Suer_gse88958_pheno$sample_ID) # "R29"  "R30"  "NR32"
# Convert NR31_2 in Suer_gse88958_pheno to NR32
Suer_gse88958_pheno$sample_ID[Suer_gse88958_pheno$sample_ID == 'NR31_2'] <- 'NR32'


## 4) Prepare Suer_pheno
Suer_pheno <- Suer_pheno %>%
  # Gleason
  mutate(gleason = case_when(Gleason.Score=='2+3' ~ 5,
                             Gleason.Score=='3+3' ~ 6,
                             Gleason.Score=='3+4'| Gleason.Score=='4+3' ~ 7,
                             Gleason.Score=='4+4' ~ 8,
                             Gleason.Score=='4+5' ~ 9,
                             Gleason.Score=='5+5' ~ 10,
                             TRUE~0)) %>%
  # Age
  mutate(AGE= floor(AGE)) %>%
  # Get recurrence column
  mutate(BCR = Suer_gse88958_pheno$`recurrence status:ch1`[match(LAB.ID, Suer_gse88958_pheno$sample_ID)]) %>% 
  # Convert it into binary, 1:no BCR, 2: no BCR
  mutate(BCR=case_when(BCR=='recurrent' ~ 2,
                       BCR=='non-recurrent' ~ 1,
                       TRUE ~ 0)) %>% 
  # Get time to BCR
  mutate(days_to_BCR=floor(as.numeric(DMOS) * 30.4167)) %>%
  # include patients with expression info only
  mutate(GSM.ID = Suer_gse88958_pheno$geo_accession[match(LAB.ID, Suer_gse88958_pheno$sample_ID)]) %>% 
  filter(!is.na(GSM.ID)) %>%
  # remove patients without follow-up data
  filter(!is.na(days_to_BCR))
rownames(Suer_pheno) <- Suer_pheno$GSM.ID  


## 4) Get expression data ready
# Subset expression profile of hsa miRs + samples with time to BCR info
ind <- grep('hsa', Suer_gse88958_feat$miRNA_ID) # hsa miRs
ind.GSM <- Suer_pheno$GSM.ID #samples with time to BCR info
Suer_gse88958_feat <- Suer_gse88958_feat[ind,]
Suer_gse88958_exprs <- Suer_gse88958_exprs[ind, ind.GSM]
rownames(Suer_gse88958_exprs) <- Suer_gse88958_feat$miRNA_ID
# As multiple probes for miRs, calculate median of each miR and select that probe
# Calculate median for all miRs
Suer_gse88958_exprs <- t(Suer_gse88958_exprs)
miRs <- unique(colnames(Suer_gse88958_exprs))
# Initiate vector to save results from the loop
Suer_exprs <- vector()
for(i in 1:length(miRs)){
  temp <- calc.df.median(Suer_gse88958_exprs, miR=miRs[i])
  Suer_exprs <- cbind(Suer_exprs, temp)
  print(i)
}
colnames(Suer_exprs) <- miRs
# Scale/ z-normalise data here
Suer_exprs <- t(Suer_exprs)
Suer_exprs <- (Suer_exprs- apply(Suer_exprs, 1, median, na.rm=T))/apply(Suer_exprs, 1, sd, na.rm=T)
Suer_exprs <- t(Suer_exprs)

## 5) Merge pheno + exprs datasets for Cox PH
rownames(Suer_exprs) == rownames(Suer_pheno)
Suer.cox.df <- cbind(Suer_exprs, Suer_pheno)

# Remove columns wiht NA values
NAs <- apply(Suer.cox.df[,1:851],2, function(x) length(which(is.na(x))))
-which(NAs>28) # if more than 28 samples have NA value for the miR remove them, as cannot calculate stuff
Suer.cox.df <- Suer.cox.df[,-which(NAs>28)] # 30 856

###################################
########## Long gse26247 ##########
###################################

## 1) Load data from GEO
Long_gse26247 <- getGEO("GSE26247", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(Long_gse26247) > 1) idx <- grep("GPL11350", attr(Long_gse26247, "names")) else idx <- 1
Long_gse26247 <- Long_gse26247[[idx]]


# log2 transform as the values are not in log spcace
ex <- exprs(Long_gse26247)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(Long_gse26247) <- log2(ex) }


## 2) Pheno, features and expression datasets
Long_gse26247_exprs <- exprs(Long_gse26247) # 1145 82
Long_gse26247_pheno <- pData(Long_gse26247) # 82 514
colnames(Long_gse26247_pheno) <- gsub('\\/| ', '\\.', gsub('\\:ch1', '', colnames(Long_gse26247_pheno)))
Long_gse26247_feat <- fData(Long_gse26247) # 1145 10

## 3) Dealing with replicate samples
# Find replicate samples and average their exprs info
rep.samples <- names(which(table(Long_gse26247_pheno$patient.number) > 1))
Long_gse26247_pheno$rep <- NA
Long_gse26247_pheno$rep <- ifelse(Long_gse26247_pheno$patient.number %in% rep.samples, 'rep', 'non-rep')
# For the replicates get the mean counts
Long_gse26247_exprs <- data.frame(t(Long_gse26247_exprs))
Long_gse26247_exprs$patient.numbers <- Long_gse26247_pheno$patient.number[match(rownames(Long_gse26247_exprs), Long_gse26247_pheno$geo_accession)]
Long_gse26247_exprs <- Long_gse26247_exprs %>% group_by(patient.numbers) %>% summarise_all(.funs=mean)
# Set rownames for exprs data
Long_gse26247_exprs <- data.frame(Long_gse26247_exprs)
rownames(Long_gse26247_exprs) <- Long_gse26247_exprs$patient.numbers
Long_gse26247_exprs$patient.numbers <- NULL
# Subset expression profile of human miRs only
colnames(Long_gse26247_exprs) <- gsub('hsa-', '', Long_gse26247_feat$TargetID[match(colnames(Long_gse26247_exprs), Long_gse26247_feat$PROBE_ID)])


## 4) Getting pheno data ready
Long_gse26247_pheno <- Long_gse26247_pheno %>%
  # Subset useful columns
  select(geo_accession, biochemical.recurrence, gleason.score, patient.number, psa, time.bcr, time.clin.rec, time.f.u, time.met.rec, time.no.bcr, tumor.stage) %>% 
  # remvoe replicates from pheno data
  group_by(patient.number) %>% distinct(patient.number, .keep_all=T) %>% 
  # Set as numeric certain columns with continuous variables
  mutate_at(.funs=as.numeric, vars(gleason.score, psa, time.bcr, time.clin.rec, time.f.u, time.met.rec, time.no.bcr)) %>% 
  # Calculate time to event
  mutate(time.to.event = case_when(biochemical.recurrence=='Recurrence'~time.bcr,
                                   TRUE ~ time.f.u)) %>%
  mutate(time.to.event = time.to.event * 30.417) %>%
  # Assign binary
  mutate(BCR=case_when(biochemical.recurrence == 'Non-Recurrence' ~ 1,
                       TRUE ~ 2)) %>%
  mutate(tumor.stage =  paste0('T', tumor.stage)) %>% 
  # remove samples without follow up info
  filter(!is.na(time.to.event))
# Set rownames 
Long_gse26247_pheno <- data.frame(Long_gse26247_pheno)
rownames(Long_gse26247_pheno) <- Long_gse26247_pheno$patient.number


## 5) Long_gse26247_pheno does not contain Age information. Upload the supplementary file with clinical info from the publication.
temp_pheno <- read.csv('data/GSE26247_clinical.csv', stringsAsFactors = F, check.names = F)
temp_pheno <- temp_pheno[, c('Time Clin Rec', 'Pre-PSA', 'Gleason Score', 'Age', 'Stage')]
# change colnames of columns used to merge the two datasets
colnames(temp_pheno)[1:2] <- c('time.clin.rec', 'psa')
# Merge and delete duplicated columns used for double checking
Long_gse26247_pheno <- left_join(Long_gse26247_pheno, temp_pheno, by=c('time.clin.rec', 'psa'))
rownames(Long_gse26247_pheno) <- Long_gse26247_pheno$patient.number

## 6) Scale/ z-normalise exprs data here
Long_gse26247_exprs <- t(Long_gse26247_exprs[rownames(Long_gse26247_pheno),]) 
Long_gse26247_exprs <- (Long_gse26247_exprs- apply(Long_gse26247_exprs, 1, median, na.rm=T))/apply(Long_gse26247_exprs, 1, sd, na.rm=T)
Long_gse26247_exprs <- t(Long_gse26247_exprs)

## 7) Merge pheno + exprs datasets for Cox PH
rownames(Long_gse26247_exprs) == rownames(Long_gse26247_pheno)
Long.test.cox.df <- cbind(Long_gse26247_exprs, Long_gse26247_pheno)


###################################
########## Long_GSE26245 ##########
###################################

## 1) Load data from GEO
Long_GSE26245 <- getGEO("GSE26245", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(Long_GSE26245) > 1) idx <- grep("GPL11350", attr(gset, "names")) else idx <- 1
Long_GSE26245 <- Long_GSE26245[[idx]]

# log2 transform as the values are not in log spcace 
ex <- exprs(Long_GSE26245)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(Long_GSE26245) <- log2(ex) }

## 2) Pheno, features and expression datasets
Long_gse26245_exprs <- exprs(Long_GSE26245) # 733 82
Long_gse26245_feat <- fData(Long_GSE26245) # 733 10
Long_gse26245_pheno <- pData(Long_GSE26245) # 82 51
colnames(Long_gse26245_pheno) <- gsub('\\/| ', '\\.', gsub('\\:ch1', '', colnames(Long_gse26245_pheno)))
Long_gse26245_pheno$title2 <- gsub('N_|Pca_', '', Long_gse26245_pheno$title)

## 3) Dealing with replicate samples
# Find replicate samples and remove them in Long_gse26245_pheno
rep.samples <- names(which(table(Long_gse26245_pheno$patient.number) > 1)) # PT131, PT5, Unknown
# For the normal samples with Unknown ID name them with GEO_accession
Long_gse26245_pheno$patient.number[Long_gse26245_pheno$patient.number == 'Unknown'] <- Long_gse26245_pheno$geo_accession[Long_gse26245_pheno$patient.number == 'Unknown']
# PT5 is actually not a replicate but 2 samples : 1N, 1T
ind <- which(Long_gse26245_pheno$patient.number == 'PT5' & Long_gse26245_pheno$group == 'Normal')
Long_gse26245_pheno$patient.number[ind] <- Long_gse26245_pheno$geo_accession[ind]
## Add rep columns
rep.samples <- names(which(table(Long_gse26245_pheno$patient.number) > 1))
Long_gse26245_pheno$rep <- NA
Long_gse26245_pheno$rep <- ifelse(Long_gse26245_pheno$patient.number %in% rep.samples, 'rep', 'non-rep')
# Mean counts for the replicates
Long_gse26245_exprs <- data.frame(t(Long_gse26245_exprs))
Long_gse26245_exprs$patient.numbers <- Long_gse26245_pheno$patient.number[match(rownames(Long_gse26245_exprs), Long_gse26245_pheno$geo_accession)]
Long_gse26245_exprs <- Long_gse26245_exprs %>% group_by(patient.numbers) %>% summarise_all(.funs=mean)
# Set rownames
Long_gse26245_exprs <- data.frame(Long_gse26245_exprs)
rownames(Long_gse26245_exprs) <- Long_gse26245_exprs$patient.numbers
Long_gse26245_exprs$patient.numbers <- NULL
colnames(Long_gse26245_exprs) <- gsub('hsa-', '', Long_gse26245_feat$TargetID[match(colnames(Long_gse26245_exprs), Long_gse26245_feat$PROBE_ID)])


## 4) Long_gse26247_pheno does not contain Age information. Upload the supplementary file with clinical info from the publication
temp_pheno <- read.csv('data/GSE26245_clinical.csv', stringsAsFactors = F, check.names = F) # 70 21
temp_pheno <- temp_pheno[,c('our id', 'Sample ID', 'Time Clin Rec', 'Pre-PSA', 'Gleason Score', 'Age', 'Stage')]
# Rename "PCS-000017, PCS-000064" as "PCS-000017"
temp_pheno$`Sample ID`[temp_pheno$`Sample ID` == "PCS-000017, PCS-000064"] <- 'PCS-000017'


## 5) Prepare pheno data
Long_gse26245_pheno[Long_gse26245_pheno == 'Unknown'] <- NA
Long_gse26245_pheno <- Long_gse26245_pheno %>%
  # remove samples without follow-up
  filter(!is.na(biochemical.recurrence)) %>%
  # remove rep samples
  group_by(patient.number) %>%
  distinct(patient.number, .keep_all=T) %>%
  # add age info 
  mutate(age = floor(temp_pheno$Age[match(title2, temp_pheno$`Sample ID`)]))%>% 
  # seect useful columns
  select(title, title2, geo_accession, biochemical.recurrence, gleason.score, group, patient.number, psa, time.bcr, time.clin.rec, time.f.u, time.met.rec, time.no.bcr, tumor.stage, age) %>%
  # set as numeric
  mutate_at(.funs=as.numeric, vars(gleason.score, psa, time.bcr, time.clin.rec, time.f.u, time.met.rec, time.no.bcr, age)) %>% 
  # get time to event info
  mutate(time.to.event=case_when(biochemical.recurrence == 'Recurrence'~time.bcr,
                                 TRUE ~ time.f.u)) %>%
  mutate(time.to.event=floor(time.to.event * 30.417)) %>%
  # Get binary BCR
  mutate(BCR = case_when(biochemical.recurrence == 'Recurrence' ~ 2,
                         biochemical.recurrence == 'Non-Recurrence' ~ 1,
                         TRUE ~ 0)) %>% 
  # Remove the normal, without follow-up
  filter(group == 'Prostate Cancer') %>%
  filter(!is.na(tumor.stage))
# Set rownames 
Long_gse26245_pheno <- as.data.frame(Long_gse26245_pheno)
rownames(Long_gse26245_pheno) <- Long_gse26245_pheno$patient.number


## 6) Scale/ z-normalise data here
Long_gse26245_exprs <- t(Long_gse26245_exprs[rownames(Long_gse26245_pheno),])
Long_gse26245_exprs <- (Long_gse26245_exprs- apply(Long_gse26245_exprs, 1, median, na.rm=T))/apply(Long_gse26245_exprs, 1, sd, na.rm=T)
Long_gse26245_exprs <- t(Long_gse26245_exprs)


## 7) Merge pheno + exprs datasets for Cox PH
rownames(Long_gse26245_pheno) == rownames(Long_gse26245_exprs)
Long.train.cox.df <- cbind(Long_gse26245_exprs, Long_gse26245_pheno)


################################
###### MIRNAME CONVERSION ######
################################

## Convert miR names for the databases to the most recent mirbase version name

## 1) Get miR version 22 names
checkMiRNAVersion(colnames(Suer.cox.df)[1:847], verbose=T) #v21
# Get accession ID 
temp.TCGA <- miRNA_NameToAccession(paste('hsa-',colnames(TCGA.cox.df)[1:328]), version='v21')
temp.Leite <- miRNA_NameToAccession(paste('hsa-',  colnames(Leite.cox.df)[1:847]), version='v11')
temp.Long.test <- miRNA_NameToAccession(paste('hsa-', colnames(Long.test.cox.df)[1:1145]), version='v12')
temp.Long.train <- miRNA_NameToAccession(paste('hsa-', colnames(Long.train.cox.df)[1:733]), version='v14')
temp.Taylor <- miRNA_NameToAccession(colnames(Taylor.cox.df)[1:373], version='v10_1')
temp.Suer <- miRNA_NameToAccession(colnames(Suer.cox.df)[1:847], version='v17')


## 2) Union of the accession ID for all the studies
miRnames.df <- data.frame(Accession=Reduce(union, list(temp.TCGA$Accession, temp.Leite$Accession, temp.Long.train$Accession, temp.Long.test$Accession, temp.Taylor$Accession, temp.Suer$Accession)), stringsAsFactors=F)
miRnames.df$V22 <- miRNA_AccessionToName(miRnames.df$Accession, targetVersion='v22')
miRnames.df$V22.mirname <- gsub('hsa-', '', miRnames.df$V22$TargetName)

## 3) Add database names
miRnames.df$TCGA <- temp.TCGA$miRNAName_v21[match(miRnames.df$V22$Accession, temp.TCGA$Accession)] 
miRnames.df$Taylor <- temp.Taylor$miRNAName_v10_1[match(miRnames.df$V22$Accession, temp.Taylor$Accession)] 
miRnames.df$Leite <- temp.Leite$miRNAName_v11[match(miRnames.df$V22$Accession, temp.Leite$Accession)] 
miRnames.df$Long.train <- temp.Long.train$miRNAName_v14[match(miRnames.df$V22$Accession, temp.Long.train$Accession)] 
miRnames.df$Long.test <- temp.Long.test$miRNAName_v12[match(miRnames.df$V22$Accession, temp.Long.test$Accession)] 
miRnames.df$Suer <- temp.Suer$miRNAName_v17[match(miRnames.df$V22$Accession, temp.Suer$Accession)] 

## 4) Get miR version 22 names
# TCGA
temp.TCGA <- miRnames.df$V22$TargetName[match(colnames(TCGA.cox.df)[1:328], gsub('hsa-', '', miRnames.df$TCGA))]
# For the miRs without the name copy their name manually
ind <- which(is.na(temp.TCGA))
temp.TCGA[ind] <- colnames(TCGA.cox.df)[ind]
# Leite
temp.Leite <- miRnames.df$V22$TargetName[match(colnames(Leite.cox.df)[1:847], gsub('hsa-', '', miRnames.df$Leite))]
ind <- which(is.na(temp.Leite))
temp.Leite[ind] <- colnames(Leite.cox.df)[ind]
# Long.test
temp.Long.test <- miRnames.df$V22$TargetName[match(colnames(Long.test.cox.df)[1:1145], gsub('hsa-', '', miRnames.df$Long.test))]
ind <- which(is.na(temp.Long.test))
temp.Long.test[ind] <- colnames(Long.test.cox.df)[ind]
# Long_GSE26245
temp.Long.train <- miRnames.df$V22$TargetName[match(colnames(Long.train.cox.df)[1:733], gsub('hsa-', '', miRnames.df$Long.train))]
ind <- which(is.na(temp.Long.train))
temp.Long.train[ind] <- colnames(Long.train.cox.df)[ind]
# Taylor
temp.Taylor <- miRnames.df$V22$TargetName[match(colnames(Taylor.cox.df)[1:373], miRnames.df$Taylor)]
ind <- which(is.na(temp.Taylor))
temp.Taylor[ind] <- colnames(Taylor.cox.df)[ind]
# Suer
temp.Suer <- miRnames.df$V22$TargetName[match(colnames(Suer.cox.df)[1:847], miRnames.df$Suer)]
ind <- which(is.na(temp.Suer))
temp.Suer[ind] <- colnames(Suer.cox.df)[ind]


## 5) Reasign colnames and remove hsa from the name too
colnames(TCGA.cox.df)[1:328] <- gsub('hsa-', '', temp.TCGA)
colnames(Leite.cox.df)[1:847] <- gsub('hsa-', '', temp.Leite)
colnames(Long.test.cox.df)[1:1145] <- gsub('hsa-', '', temp.Long.test)
colnames(Long.train.cox.df)[1:733] <- gsub('hsa-', '', temp.Long.train)
colnames(Taylor.cox.df)[1:373] <- gsub('hsa-', '', temp.Taylor)
colnames(Suer.cox.df)[1:847] <- gsub('hsa-', '', temp.Suer)

common.miRs <- Reduce(intersect, list(colnames(TCGA.cox.df)[1:328],
                                      colnames(Leite.cox.df)[1:847],
                                      colnames(Long.test.cox.df)[1:1145],
                                      colnames(Long.train.cox.df)[1:733],
                                      colnames(Taylor.cox.df)[1:373],
                                      colnames(Suer.cox.df)[1:847])) #162


save.image('code/MetaAnalysis_db.Rdata')
