
### Merge TCGA Mature miRNA data ###

# Function to generate a table of mature miRNA data from TCGA miRseq data downloaded from GDC
# Adapted from code on github/rptashkin/TCGA_miRNASeq_matrix

## INPUT : 	list of mirfileIDs - folder names in GDC download folder (each containing isoforms.quantification.txt)
##			table of mirbase names -  table containing mature miRNA names and accession numbers, from miRbase (mature.fa)

## OUTPUT : data frame containing counts for each mature miRNA in each sample

mergeTCGAmaturemirna <- function(path='/home/sharmila/Documents/github_projects/ch3_Review_MetaAnalysis/data/TCGA/gdc_download_20170213_163009/', mirfileIDs=mir_manifest$id, mirbasenames=mirbasehuman){
  
  testcounts<- lapply(mirfileIDs, function(x){
    # read datafile with isoforms
    tempfile <- read.table(file=paste(path, x ,"/isoforms.quantification.txt", sep=""), header=TRUE, stringsAsFactors=FALSE)
    
    # split out mirbase accession numbers (MIMAT)
    tempfile <- splitstackshape::cSplit(tempfile, "miRNA_region", sep=",")
    
    # create column with mature miRNA name matched to mirbase accession	(MIMAT)
    tempfile$mirname <- qdapTools::lookup(tempfile$miRNA_region_2, mirbasenames$accession, mirbasenames$maturemirna)
    
    # take mature mir names and counts
    tempcounts <- data.frame(miRNA = tempfile$mirname, count = tempfile$read_count)
    
    # remove rows not matched to mature mirna names
    tempcounts <- tempcounts[!(is.na(tempcounts$miRNA)),]
    tempcounts <- tempcounts[order(tempcounts$miRNA),]
    
    # sum counts for rows containing the same mature mirna
    tempcounts <- aggregate(data=tempcounts, count ~ miRNA, FUN=sum)
    
    # rename count column with file ID
    colnames(tempcounts) <- c("miRNA", paste(x))
    return(tempcounts)
  })
  
  
  datamatrix <- plyr::join(testcounts[[1]], testcounts[[2]], by="miRNA", type="full")
  
  if (length(mirfileIDs) > 2){
    for (i in 3:length(mirfileIDs)){
      datamatrix <- plyr::join(datamatrix, testcounts[[i]], by="miRNA", type="full")
    }
  }
  
  datamatrix <- datamatrix[order(datamatrix$miRNA),]
  
  datamatrix[is.na(datamatrix)] <- 0
  
  return(datamatrix)
}
