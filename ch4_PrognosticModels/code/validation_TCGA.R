
########################################################
###### Whole TCGA-dataset for internal validation ######
########################################################

# load('/home/sharmila/Documents/github_projects_with_data/ch4_PrognosticModels/code/PrognosticMirModels_200630.Rdata'); keep(counts, sure=T)
# load('/home/sharmila/Documents/github_projects_with_data/ch4_PrognosticModels/code/PrognosticMirModels_200722.Rdata'); keep(EN.res.wt.mean, counts.norm.Train, counts.norm.Test, recur.samp.df.tmp, sure=T)
# load('/home/sharmila/Documents/github_projects_with_data/ch4_PrognosticModels/code/PrognosticMirModels_200708.Rdata'); keep(features.iso, mod.cluster, mod.isotype, mod.parent.mir, mod.seed6mer, mod.seed7mer, mod.t3size, mod.t5size, sure=T) 
# save.image('/home/sharmila/Documents/github_projects_with_data/ch4_PrognosticModels/code/validation/allTCGA_validation.Rdata')
# load('/home/sharmila/Documents/github_projects_with_data/ch4_PrognosticModels/code/validation/allTCGA_validation.Rdata')

## 1) Subset raw data for the counts_norm
valid.counts <- counts[grep('hsa-', colnames(counts.norm.Train$part1), value=T), row.names(counts.norm.Train$part1)]

## 2) Normalise
dgList <- DGEList(counts=valid.counts, genes=rownames(valid.counts))
dgList <- calcNormFactors(dgList, method='TMM')
valid.counts.norm <- cpm(dgList, normalized.lib.sizes = T, log=T)
# Log as later on we are applying reglarised regression whose assumption is normally distributed data
temp <- rownames(valid.counts.norm)
valid.counts.norm <- data.frame(t(valid.counts.norm))
colnames(valid.counts.norm) <- temp

## 5) Firstly create module sets and assign features into the sets
## 2) Firstly create module sets and assign features into the sets
features.part1 <- features.iso[features.iso$isomir %in% grep('hsa-', colnames(counts.norm.Train$part1), value=T),]
features.part1$isotype_148a <- ifelse(features.part1$mir == 'miR-148a-3p', features.part1$isotype, 'not_miR-148a')
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
# viii) isotypes + 148a
mod.148a.isotype <- lapply(unique(features.part1$isotype_148a), function(x) features.part1$isomir[which(features.part1$isotype_148a %in% x)]);
names(mod.148a.isotype) <- unique(features.part1$isotype_148a)
# Remove not_miR-148a
mod.148a.isotype <- mod.148a.isotype[-1]


## 3) Transpose counts norm
#counts.norm.t <- list()
valid.counts.norm.t <- data.frame(t(valid.counts.norm), check.names=F)

## 3.5)
scoreSignature <- function(signature,exprs){
  fcs <- exprs-apply(exprs,MARGIN=1,median) # centre miR expression around its median
  outScores <- sapply(signature, function(x) colMeans(fcs[x,]))
  outScores
}

## 4) Calculate z-score modules
mean.mod <- list()
# parent miR
mean.mod$parent <- scoreSignature(signature=mod.parent.mir,exprs=valid.counts.norm.t)
# cluster
mean.mod$cluster <- scoreSignature(signature=mod.cluster,exprs=valid.counts.norm.t)
# isotype
mean.mod$isotype <- scoreSignature(signature=mod.isotype,exprs=valid.counts.norm.t)
# t3size
mean.mod$t3size <- scoreSignature(signature=mod.t3size,exprs=valid.counts.norm.t)
# t5size
mean.mod$t5size <- scoreSignature(signature=mod.t5size,exprs=valid.counts.norm.t)
# seed 7merm8
mean.mod$seed7mer <- scoreSignature(signature=mod.seed7mer,exprs=valid.counts.norm.t)
# seed 6mer
mean.mod$seed6mer <- scoreSignature(signature=mod.seed6mer,exprs=valid.counts.norm.t)
# isotypes 148a
mean.mod$isotype.148a <- scoreSignature(signature=mod.148a.isotype,exprs=valid.counts.norm.t)

## 5) Append clinical info on to the datasets
mean.mod.clin <- lapply(mean.mod, function(x) merge(x, recur.samp.df.tmp, by='row.names'))
mean.mod.clin <- lapply(mean.mod.clin,
                        function(x){rownames(x) <- x$Row.names;
                        x$Row.names <- NULL;
                        x$file.id.isomir <- NULL;
                        return(x)})

valid.counts.norm <- data.frame(t(valid.counts.norm.t), check.names=F)
valid.counts.norm <- merge(valid.counts.norm, recur.samp.df.tmp, by='row.names')
rownames(valid.counts.norm) <- valid.counts.norm$Row.names
valid.counts.norm$Row.names <- NULL; valid.counts.norm$file.id.isomir <- NULL

## 5) Perform EN validation
## Function to do internal validation
ind.val <- function(MODEL = EN.res.wt.mean$CVonly, df=Test.df) {
  
  # Make predictions on test
  test <- list()
  test$prob <- MODEL$model %>% predict(df, type='prob')
  test$class <- MODEL$model %>% predict(df, type='raw')
  test$con.mat <- confusionMatrix(data=test$class, df$BCR.event)
  ## Get AUC + ROC ED
  test$ROC.ED <- makeROC(x=test$prob[,2], y=ifelse(df$BCR.event=='BCR', 1, 0))
  test$AUC.ED <- getAUC(test$ROC.ED)
  ## GET AUC + ROC pROC
  ## pROC
  test$ROC.pROC <- roc(df$BCR.event, test$prob[,1],  levels=c('BCR', 'no_BCR'), direction='>', plot=F)
  test$AUC.pROC <- summary(auc(test$ROC.pROC))['Median']
  ## GET AUC + ROC ROCR
  test$ROC.ROCR <- prediction(test$prob[,2], df$BCR.event)
  y <- performance(test$ROC.ROCR, measure='auc'); test$AUC.ROCR <- y@y.values[[1]]
  
  return(test)
}

## Validation
Valid.res <- list()
Valid.res$CVonly <- ind.val(MODEL = EN.res.wt.mean$CVonly, df=valid.counts.norm)
Valid.res$indIso.only <- ind.val(MODEL = EN.res.wt.mean$indIso.only, df=valid.counts.norm)
Valid.res$indIso <- ind.val(MODEL = EN.res.wt.mean$indIso, df=valid.counts.norm)
temp.148a <- valid.counts.norm %>%
  select(c(matches('148a-3p'),'BCR.event','pT.grp','ISUP.gleason','stage_event_psa'))
Valid.res$iso_148a <- ind.val(MODEL = EN.res.wt.mean$iso_148a, df=temp.148a)
## Model with only 582-5p isomiRs
temp.582 <- valid.counts.norm %>%
  select(c(matches('582-5p'),'BCR.event','pT.grp','ISUP.gleason','stage_event_psa'))
Valid.res$iso_582 <- ind.val(MODEL = EN.res.wt.mean$iso_582, df=temp.582)
## Model with only 148a + 582-5p isomiRs
temp.148a_582 <- valid.counts.norm %>%
  select(c(matches('148a-3p|582-5p'),'BCR.event','pT.grp','ISUP.gleason','stage_event_psa'))
Valid.res$iso_148a_582 <- ind.val(MODEL = EN.res.wt.mean$iso_148a_582, df=temp.148a_582)
Valid.res$parentmiR <- ind.val(MODEL = EN.res.wt.mean$parentmiR, df=mean.mod.clin$parent)
Valid.res$cluster <- ind.val(MODEL = EN.res.wt.mean$cluster, df=mean.mod.clin$cluster)
Valid.res$isotype <- ind.val(MODEL = EN.res.wt.mean$isotype, df=mean.mod.clin$isotype)
Valid.res$t3size <- ind.val(MODEL = EN.res.wt.mean$t3size, df=mean.mod.clin$t3size)
Valid.res$t5size <- ind.val(MODEL = EN.res.wt.mean$t5size, df=mean.mod.clin$t5size)
Valid.res$seed7mer <- ind.val(MODEL = EN.res.wt.mean$seed7mer, df=mean.mod.clin$seed7mer)
Valid.res$seed6mer <- ind.val(MODEL = EN.res.wt.mean$seed6mer, df=mean.mod.clin$seed6mer)
Valid.res$isotype.148a <- ind.val(MODEL = EN.res.wt.mean$isotype.148a, df=mean.mod.clin$isotype.148a) #mhere

### Assess results
Valid.res.df <- data.frame(models=names(Valid.res), stringsAsFactors=F) #add: 't5size', 148a+582 isomiRs
Valid.res.df <- Valid.res.df %>%
  #mutate(AUC.ED=signif(sapply(Valid.res, function(x) x$AUC.ED), 3)) %>%
  mutate(AUC.pROC=signif(sapply(Valid.res, function(x) x$AUC.pROC), 3)) %>%
  #mutate(AUC.ROCR=signif(sapply(Valid.res, function(x) x$AUC.ROCR), 3)) %>%
  #mutate(p.value=NA) %>%
  mutate(accuracy=signif(as.numeric(sapply(Valid.res, function(x) 
    coords(roc=x$ROC.pROC, x=0.6, input='specificity', ret='precision', transpose=F))), 3)) %>%
  mutate(sensitivity=signif(as.numeric(sapply(Valid.res, function(x) 
    coords(roc=x$ROC.pROC, x=0.6, input='specificity', ret='sensitivity', transpose=F))), 3)) %>%
  mutate(specificity=signif(as.numeric(sapply(Valid.res, function(x) 
    coords(roc=x$ROC.pROC, x=0.6, input='specificity', ret='specificity', transpose=F))), 3)) %>%
  mutate(precision=signif(as.numeric(sapply(Valid.res, function(x) 
    coords(roc=x$ROC.pROC, x=0.6, input='specificity', ret='precision', transpose=F))), 3)) %>% 
  mutate(NPV=signif(as.numeric(sapply(Valid.res, function(x) 
    coords(roc=x$ROC.pROC, x=0.6, input='specificity', ret='npv', transpose=F))), 3))




## 1) Plot ROC curves
# colours
n=length(Valid.res); qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# # plot
# pdf('Documents/github_projects_with_data/ch4_PrognosticModels/results/EN_wt_Ed_AUC1_validation.pdf', height=5, width=5.5)
# for (i in 1:length(Valid.res)){
#   plotROC(Valid.res[[i]]$ROC.ED, col=col_vector[6], roc.title=names(Valid.res)[i])}
# dev.off()


pdf('Documents/github_projects_with_data/ch4_PrognosticModels/results/EN_wt_pAUC_validation.pdf', height=5, width=6)

for (i in 1:length(Valid.res)){
  plot(Valid.res[[i]]$ROC.pROC, print.thres='no',  print.auc=T, print.auc.y=0.1, print.auc.x=0.3, print.auc.cex=1.4, type='S', col=col_vector[6], lwd=2, xlim=c(1,0), legacy.axes=T, cex.main=1.4,cex.lab=1.4,cex.axis=1.4) #main=names(Valid.res)[[i]], 
}
dev.off()


####
# save.image('/home/sharmila/Documents/github_projects_with_data/ch4_PrognosticModels/code/validation/allTCGA_validation.Rdata')
