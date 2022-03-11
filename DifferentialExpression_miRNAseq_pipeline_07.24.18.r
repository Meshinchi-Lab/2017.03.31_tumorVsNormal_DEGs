#Jenny Smith

#Aug 16, 2017

#Differential Expression Analysis pipeline

#Purpose: Given matrix of RNA-seq raw counts and clinical annoations, create a series of differential expression Analyses for two group comparisions. 
setwd(file.path(TARGET,"RNA/miRNAseq/2017.03.31_tumorVsNormal_DEGs"))

library(ggplot2)
library(pryr)

theme_numX %<a-% { theme(plot.title = element_text(hjust = 0.5, size = 18),
                         panel.background = element_rect(fill="white"),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.border = element_rect(color = "black", fill=NA),
                         axis.text = element_text(color = "black"),
                         axis.text.x = element_text(angle = 0,hjust=0.5,vjust = 0.5, size = 16),
                         axis.text.y = element_text(size = 16),
                         axis.title = element_text(size = 18))
}

GroupIDs <- function(clinData, col){
  #clindata has patient IDs as rownames. 
  #col is a chracter string of the column name in the clinical data with the factor/variable information. 
  list <- list()
  grps <- unique(clinData[,col])
  N <- length(grps)
  
  for (i in 1:length(grps)){
    if (grepl("[^a-zA-Z0-9 :]", grps[i])){
      grps[i] <- gsub("[^a-zA-Z0-9 :]", "\\.", grps[i]) #remove special characters and replace with a "."
    }
    IDs <- rownames(clinData[grepl(grps[i], clinData[,col]), ])
    list[[i]] <- IDs
  }
  
  names(list) <- grps
  return(list)
}

phenoVectors <- function(groupA, groupB){
  library(magrittr)
  #groupA and GroupB are character vectors with the patients IDs in each group
  g1 <- as.character(substitute(groupA))
  g2 <- as.character(substitute(groupB)) 
  
  vector <- c(rep(g1, length(groupA)), rep(g2, length(groupB)))
  names(vector) <- c(groupA, groupB)
  
  return(vector)
}

merge_CDE_Expn <- function(clinData, expnMatrix,geneList, phenoVector=NULL){
  #expnData is a data frame with patient IDs as the column names, genes as rownames
  #Clindata is a data frame with patient IDs in rownames.  
  
  #subset for genes of interest
  expnMatrix <- expnMatrix[rownames(expnMatrix) %in% geneList, ] 
  
  #ensure # of rows == # cols of expression data for later merging
  expnMatrix <- expnMatrix[,intersect(colnames(expnMatrix), rownames(clinData))]
  
  if (is.null(phenoVector)){
    tmp <- data.frame(t(expnMatrix))
  }else{
    phenoVector <- phenoVector[intersect(names(phenoVector), rownames(clinData))]
    expnMatrix <- expnMatrix[, names(phenoVector)]
    tmp <- data.frame(t(expnMatrix),
                      Status=phenoVector)
  }
  
  #merge the clinical and expn data
  srt_clinData <- merge(clinData, tmp, by.y=0, by.x=0)
  rownames(srt_clinData) <- srt_clinData$Row.names
  
  return(srt_clinData)
}


calcDEMir <- function(expnData, pheno,ref,RPM=FALSE, nsamp=0.05) {
  # expnData is a matrix or data frame with the raw counts. Patient IDs as colnames, genes as rownames
  # pheno is a phenotype vector, a named character vector, like c("patient1"="AML", "patient2"="NBM")
  #ref is which of the groups in teh phenovector is the reference group, eg ref="NBM"
  #logCPM is whether to log2 CPM normalize the raw counts.
  
  library(limma)
  library(edgeR)
  
  #ensure correct order
  expnData <- expnData[,match(names(pheno), colnames(expnData))]
  
  if (!all(complete.cases(expnData))){
    print("Names DO NOT match in between phenovector and colnames of expression matrix")
    return(list("expnData"=expnData,"pheno"=pheno))
  }
  
  #order and relevel the phenovector 
  groups <- unique(pheno)
  groups <- c(groups[groups != ref], ref) #order so that reference is second 
  pheno.f <- factor(pheno, levels=groups)
  
  #Define the DGE object and select AML samples for filtering
  dge <- DGEList(counts = expnData, group=pheno.f) 
  AML <- ! grepl("^BM|^RO", colnames(expnData))
  AMLsamples <- sum(AML)
  
  #Create a design matrix and contrasts. 
  design <- model.matrix(~0 + pheno.f, data=dge$samples)#~0 means no intercept. 
  colnames(design) <- levels(pheno.f)
  cont.matrix <- makeContrasts(contrasts = paste(groups, collapse = "-"), 
                               levels = design) #contrast is approx. log2(mean(Pos)) - log2(mean(Neg)) per gene. 
  
  if(nsamp < 1){
    n <- nsamp*ncol(expnData)
  }else if( nsamp > 1){
    n <- nsamp 
  }
  
  if(RPM){
    #at least 1 read per million (RPM) in at least 2 or greater AML samples
    keep.dge <- rowSums(dge$counts[,AML] >= 1) >= max(2,n)
    dge <- dge[keep.dge,] #subset for those genes with RPM >= 1 per gene
    
    #log2 transform the RPMs only. 
    dge.norm <- apply(dge$counts, 2, function(x) log2(x + 1))
    NormFactors <- "Log2"
    
  }else if (!RPM){
    #at least 5 raw counts in at least 2 or greater AML samples
    keep.dge <- rowSums(dge$counts[,AML] >= 5) >= max(2,n)
    dge <- dge[keep.dge,]
    
    #yes, ok to keep. 
    dge <- calcNormFactors(dge) #calculate the TMM normalization factors
    
    dge.norm <- voom(dge, design, plot = FALSE)  
    NormFactors <- "Voom"  #voom transformed counts for sample to sample comparisons.
  }

  print(NormFactors)
  
  #fit the linear model and the contrasts fit.
  fit <- lmFit(dge.norm, design)
  fit <- contrasts.fit(fit, contrasts = cont.matrix)

  #compute moderated t-statistics using empirical bayes moderation.
  fit2 <- eBayes(fit)


  #Find the DE genes
  DE <- topTable(fit2,adjust.method="BH",sort.by="P",
               number=20000,p.value=0.05, lfc=1)

  list <- list( dge.norm, fit2, DE)
  names(list) <- c( NormFactors,"eBayesFit", "DE")

  
  return(list)
}


dendrograms <- function(df, pheno, genelist, method,log2=TRUE){
  #df with TPM or RPKM normalized data, patient IDs are column names and rows are genes. 
  #pheno is a character vector with patient IDs as names, and the status for each in each group (eg pos,neg)
  #genelist is a character vector with the genes of interest
  #log2 is to whether the input df is logged already. 
  require(edgeR)
  library(dendextend)
  
  #ensure correct order, drop rows with nas just in case
  df <- df[, intersect(names(pheno), colnames(df))] 
  df <- df[which(rownames(df) %in% genelist), ] #subset the matrix to genes of interest
  
  print(dim(df))
  
  
  if(!log2){
    #updated the prior count to be 1 rather than 0.01. Found 0.01 can cause outliers with
    #very negative log2 scale expression values bc they were so close to zero
    df <- as.data.frame(apply(df, 2, function(x) log2(x + 1))) #log2 transform counts 
    # TMMCPM10 <- as.data.frame(apply(TMMCPM, 2, function(x) log10(x + 0.01))) #log10 transform counts
  }
  
  
  d1 <- dist(t(df), method = "euclidean", diag = FALSE, 
             upper = FALSE) #sample distances WITHOUT SCALING 
  d2 <- dist(df, method = "euclidean", diag = FALSE,
             upper = TRUE) #gene distances WITHOUT SCaling
  
  c1 <- hclust(d1, method = method, members = NULL) #sample clustering
  c2 <- hclust(d2, method = method, members = NULL) #gene clustering
  
  
  list <- list(df, d1,d2,c1,c2)
  names(list) <- c("expnData", "d1", "d2", "c1", "c2")
  
  return(list)
}

basicHeatmap <- function(ExpnMatrix, geneDend, sampleDend, colors,title){
  require(gplots)
  require(pryr)
  library(colorspace)
  library(dendextend)
  
  #ExpnMatrix is the genes as rownames, patient IDs as colnames
  #genedend is from hclust object
  #sample dend is from hclust objest
  #rowlabels is the rownames of the initial expn matrix
  #colors is a character vector of colors of equal length of samples to illustrate the different groups
  
  
  # colors = c(seq(-3,-2,length=100),seq(-2,0.5,length=100),seq(0.5,6,length=100))
  # breaks = seq(-5,5,length.out = 300), #this MUST be checked between datasets or it will be inaccurate
  # colorPal <- colorRampPalette(c("blue", "white", "red"))(n=299)
  # colorPal <- colorRampPalette(c("darkgreen", "forestgreen", "green3", "green2", "black", "firebrick1", "red3", "red4", "darkred"))(n=299)
  colorPal <- colorRampPalette(c("deepskyblue4", "deepskyblue3", "deepskyblue2", "deepskyblue1","white","red1", "red2", "red3", "red4"))(n=299)
  N <- ncol(ExpnMatrix)
  N.genes <- nrow(ExpnMatrix)
  rowLabels <- rownames(ExpnMatrix)
  
  if ( N.genes < 100){
    cex=0.8
  }else if(length(N.genes) > 100 & length(N.genes) < 175){
    cex=0.5
  }else{
    cex=0.25
  }
  
  z.pryr %<a-% {
    par(cex.main=1.5, cex=0.75, font=2, font.axis=1, lend=1) 
    heatmap.2(as.matrix(ExpnMatrix), 
              Colv=dendextend::rotate(as.dendrogram(sampleDend), 
                                      order = c(N:1)), 
              Rowv=as.dendrogram(geneDend), 
              labRow=rowLabels,
              labCol = "",
              ColSideColors = colors,
              density.info="density", #density.info="density",
              trace = "none",
              scale="row",
              col = colorPal, 
              cexRow=cex,
              margins=c(2,10), 
              lwid=c(.8,3), 
              lhei=c(.8,3), 
              srtCol=75, 
              adjCol=c(1,1),
              keysize=0.75, 
              key.title="",
              key.ylab ="",
              key.par = list(cex=0.75),
              main=title)
  }
  return(z.pryr)
}


plotPCoA <- function(expnData,phenovector, title="",colorCodes=NULL, geneList=NULL){
  #factor is the name of the factor column 
  library(vegan)
  library(ggplot2)
  
  #Ensure correct order of patients in both datasets
  expnData <- expnData[ ,intersect(names(phenovector), colnames(expnData))] 
  phenovector <- phenovector[intersect(names(phenovector),colnames(expnData))]
  
  if (! is.null(geneList)){
    expnData <- t(expnData[geneList, ])
  }else{
    expnData <- t(expnData) #note: must remove all zero count genes or will  fail on an error
  }
  
  PCoA <- capscale(expnData ~ 1, distance = "bray", add=TRUE)
  scores <- data.frame(scores(PCoA, display="sites"), 
                       Group=phenovector) 
  
  
  p <- ggplot(scores, aes(x=MDS1, MDS2)) +
    geom_point(aes(color=scores[,"Group"]), size=3) +
    theme_numX +
    labs(title=title) 
  
  if(!is.null(colorCodes)){
    p <- p + 
      scale_color_manual(values=colorCodes)
  }
  
  
  # aov <- aov(scores$MDS1 ~ df[,factor]) #NOTE: This is only valid for balanced experimental designs! Equal # of obs in each factor level.
  
  list <- list(PCoA,scores, p)
  names(list) <- c("PCoA","scores","plot")
  
  return(list)
}

#from https://github.com/mikelove/DESeq2/blob/master/R/plots.R
#Want to return the whole scores matrix so can examine 3d pca plots. 
plotPCA.DESeq.mod <- function(object, intgroup="condition", ntop=500, returnData=FALSE)
{
  library(ggplot2)
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  # d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, name=colnames(object))
  d <- data.frame(as.data.frame(pca$x)[,1:10], group=group, intgroup.df, name=colnames(object))
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:10]
    return(d)
  }
  
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed()
}


#Updated on 6/9/17 to use variance stabilized transformed data as input (not center scaled log2, like in princomp)
PCA <- function(expnData,phenovector,title="",round=TRUE,colorCodes=NULL){
  suppressPackageStartupMessages(library(DESeq2))
  library(ggplot2)
  #expnData has patient IDs as colnames and genes as rownames. 
  
  #Order the mir counts and sampes
  samples <- intersect(names(phenovector), colnames(expnData))
  countData <- expnData[,samples]
  phenovector <- phenovector[samples]
  
  countData <- round(countData, digits = 0)
  colData <- as.data.frame(phenovector)
  
  #Create a deseq object
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ phenovector)
  
  #Filter by 5 counts for each mir
  dds <- dds[ rowSums(counts(dds)) > 5, ]
  
  #variance stabalizing transformation because miRs have so few rows. vst() works on a subset of rows. 
  varianceStab <- varianceStabilizingTransformation(dds, blind = TRUE)
  
  plot <- plotPCA.DESeq.mod(varianceStab, intgroup = "phenovector") + 
    theme_numX + 
    labs(title=title) 
  
  if (!is.null(colorCodes)){
    plot <- plot + 
      scale_color_manual(values=colorCodes)
  }
  
  
  pca.dat <- plotPCA.DESeq.mod(varianceStab, intgroup = "phenovector", returnData=TRUE)
  
  list <- list(dds, varianceStab, plot, pca.dat)
  names(list) <- c("dds", "vst", "pca_plot","pca_data")
  
  return(list)
}


twoGroups_DEMirs <- function(expnData, clinData, col, ref,RPM=FALSE,BM=FALSE,nsamp=0.05){
  # expnData is a matrix or data frame with the raw counts. Patient IDs as colnames, genes as rownames
  #clindata has patient IDs as rownames. 
  #col is a character string of the factor column of interest
  #ref is the character strign of the reference group level (eg BM, Neg, or control)
  library(magrittr)
  library(genefilter)
  library(gtools)
  
  #remove unknown categories from the datasets since only want yes/no or 0/1 groups
  rmUnknowns <- function(clinData, cols){
    removeUnknowns <- clinData
    
    for (i in 1:length(cols)){
      removeUnknowns <- removeUnknowns[!grepl("Unknown",removeUnknowns[, cols[i]]), ] 
    }
    return(removeUnknowns)
  }
  
  dataName <- substitute(expnData)
  variantName <- col
  print(variantName)
  
  #remove unknowns from the analysis
  clinData <- rmUnknowns(clinData, col)
  groups <- GroupIDs(clinData, col) #list of patient IDs, one for each group
  
  #Define Groups to compare based on group IDs from clinical data. Intersect with expression matrix to subset. 
  if (BM){
    BM <- grep("^BM|^RO", colnames(expnData), value = TRUE)
    GroupB <- BM #select the reference group 
    
    GroupA <- groups[[which(names(groups) != ref)]] %>%
      intersect(. , colnames(expnData)) #the second group (mutant, AML, treated, etc)
      
  }else{
    GroupB <- groups[[ref]] %>% 
      intersect(. , colnames(expnData)) #select the reference group (eg No, normal, wt, control, etc.) Must be a character(level) from the column of clinData selected. 
    
    GroupA <- groups[[which(names(groups) != ref)]] %>% 
      intersect(. , colnames(expnData)) #the second group (mutant, AML, treated, etc)
  }
  
  #Only analyze at least 3x3 comparisons. 
  if (any(lapply(list(GroupA,GroupB), length) < 3)){
    list <- list(expnData, clinData, GroupA,GroupB)
    names(list) <- c("InputExpnData", "InputClinData", "CompGroup", "RefGroup")
    return(list)
  }
  
  #Define the pheno vector 
  phenoVector <- phenoVectors(GroupA, GroupB)
  
  #update clinical data 
  if (identical(GroupB,BM)){
    clinData <- as.data.frame(phenoVector) %>%
      set_colnames(., value=col)
    
  }else{
    clinData <- clinData[intersect(c(GroupA,GroupB),rownames(clinData)), ]
  }
  
  #subset and order the dataframes.
  samps <- c(GroupA, GroupB)
  expnData <- expnData[ ,match(samps, colnames(expnData))] #mutant (groupA), then WT (groupB)

  
  if (any(is.na(expnData))){print("NAs Introduced. Check Rownames and colnames of inputs")}
  
  
  # Calculate Differential Expression
  DE <- calcDEMir(expnData, 
                  pheno = phenoVector,
                  ref = "GroupB", 
                  RPM = RPM, 
                  nsamp=nsamp) #mutant - wild type.
  DE$DE$FoldChange <- logratio2foldchange(DE$DE$logFC)

  
  if(!RPM){
    expn.norm <- DE$Voom$E #voom() function reports log2 + 0.5 CPM expression
  }else{
    expn.norm <- DE$Log2 #anti-log2 to put back on linear scale
  }
  
  # return(list(phenoVector, expnData, clinData, DE))
 
  if (nrow(DE$DE) <= 9){
    cc <- c("GroupB"="black", "GroupA"="firebrick")
    
    #MDS with non-log, normalized expression values
    MDS <- plotPCoA(expnData = 2^expn.norm,
                    phenovector = phenoVector, 
                    geneList = rownames(expn.norm), 
                    colorCode=cc, 
                    title = variantName)

    list <- list( phenoVector, DE, MDS)
    names(list) <- c( "phenovector","DE","PCA", "MDS")
    return(list)

  }else{

    #Unsupervised Heirachrach clustering
    DEMirs <- rownames(DE$DE)
    dends_DE <- dendrograms(expn.norm, phenoVector, DEMirs, method="ward.D2", log2=TRUE) #dendrograms based on all differentially expressed genes.
  
    #Create a heatmap of all DEGs
    colorBar <- ifelse(phenoVector == "GroupB", "black", "firebrick")
    title <- variantName
    heatmap <- basicHeatmap(dends_DE$expnData, dends_DE$c2, dends_DE$c1,colorBar, title=title)

    #Unconstrained Cluster Analysis
    # PCoA on non-log, normalized expression values of DEGs
    cc <- c("GroupB"="black", "GroupA"="firebrick")
    MDS <- plotPCoA(expnData = 2^expn.norm,
                    phenovector = phenoVector, 
                    geneList = DEMirs, 
                    colorCode=cc, 
                    title = variantName)
    eigen <- PCA(expnData = expnData,
                 phenovector =  phenoVector, 
                 colorCodes=cc, 
                 title=variantName)
    

    #return the objects
    list <- list( phenoVector,DE,dends_DE, heatmap, MDS, eigen)
    names(list) <- c("phenovector", "DE", "dends_DE", "Heatmap", "MDS","PCA")


    return(list)
  }
}


####### Extraction methods to get items of interest. ############

extract_DEMirs <- function(twoGroups_DEMirs.res, filter=FALSE,goi=NULL){
  library(dplyr)
  
  if(length(nrow(twoGroups_DEMirs.res$DE$DE)) < 1){
    return("No DEGs")
  }else{
    
    DE <- twoGroups_DEMirs.res$DE$DE %>%
        mutate(gene=rownames(.)) %>%
        arrange(dplyr::desc(logFC)) %>%
        dplyr::select(gene, everything())
    
    if (filter){
      DE <- DE %>%
        filter(gene %in% goi)
    }
    
    return(DE)
  }
}


extract_MDS <- function(twoGroups_DEGs.res){
  twoGroups_DEGs.res$MDS$plot
}


extract_PCA <- function(twoGroups_DEGs.res){
  twoGroups_DEGs.res$PCA$pca_plot
}







############# Old #########################


# plotPCoA <- function(expnData,clinData, geneList=NULL, factor){
#   #factor is the name of the factor column 
#   library(vegan)
#   
#   #Ensure correct order of patients in both datasets
#   expnData <- expnData[ ,intersect(rownames(clinData), colnames(expnData))] 
#   clinData <- clinData[intersect(rownames(clinData),colnames(expnData)), ]
#   
#   if (! is.null(geneList)){
#     expnData <- t(expnData[geneList, ])
#   }else{
#     expnData <- t(expnData) #note: must remove all zero count genes or will  fail on an error
#   }
#   
#   PCoA <- capscale(expnData ~ 1, distance = "bray", add=TRUE)
#   scores <- as.data.frame(scores(PCoA, display="sites"))
#   
#   p <- ggplot(scores, aes(x=MDS1, MDS2)) +
#     geom_point(aes(color=clinData[,factor])) +
#     theme_numX +
#     labs(title="")
#   
#   # aov <- aov(scores$MDS1 ~ df[,factor]) #NOTE: This is only valid for balanced experimental designs! Equal # of obs in each factor level.
#   
#   list <- list(df, PCoA,scores, p)
#   names(list) <- c("MDS_df","PCoA","scores","plot")
#   
#   return(list)
# }







# calcDEMir <- function(expnData, g1, g2, logCPM=NULL,Trend=TRUE) {
#   # expnData is a matrix or data frame with the raw counts. Patient IDs as colnames, genes as rownames
#   # g1,g2 are the patient IDs for each group
#   #logCPM is whether to log2 CPM normalize the raw counts.
#   
#   library(limma)
#   library(edgeR)
#   
#   dge <- DGEList(counts = expnData) #miRNA input is already in RPM normalized (reads per million)
#   
#   keep.dge <- rowSums(dge$counts >= 1) >= (0.05*ncol(expnData))
#   dge <- dge[keep.dge,] #subset for those genes with RPM >= 1 per gene
#   
#   if (is.null(logCPM)){
#     dge <- calcNormFactors(dge) #calculate  the TMM normalization factors
#     return(dge)
#     
#   }else if (logCPM==TRUE){
#     dge <- calcNormFactors(dge) #calculate the normalization factors
#     NormFactors <- dge$samples
#     dge <- cpm(dge, log=TRUE, prior.count = 1) #log2 CPM transformation.
#     
#   }else if (logCPM==FALSE){
#     dge <- apply(dge$counts, 2, function(x) log2(x + 1))
#     NormFactors <- "None"
#   }
#   
#   
#   designMatrix <- matrix(0, nrow=dim(dge)[2] , ncol=2)
#   colnames(designMatrix) <- c("g1","g2")
#   rownames(designMatrix) <- colnames(dge)
#   designMatrix[g1, 1] <- 1
#   designMatrix[g2, 2] <- 1
#   
#   fit <- lmFit(dge,designMatrix)
#   tmp <- paste("g1","g2",sep="-") #contrast is ~ log2(mean(g1)) - log2(mean(g2)) per gene
#   
#   cont.matrix <- makeContrasts(contrasts=tmp,
#                                levels=designMatrix)
#   
#   fit2<-contrasts.fit(fit, cont.matrix)
#   fit2<-eBayes(fit2, trend=Trend) #added on 6/4/2017 by JS
#   DE<-topTable(fit2,adjust.method="BH",sort.by="P",
#                number=20000,p.value=0.05, lfc=1)
#   
#   
#   
#   list <- list(dge,NormFactors, designMatrix,fit2, DE)
#   names(list) <- c("dge","NormFactors", "design","eBayesFit", "DE")
#   
#   return(list)
# }


# plotPCA <- function(expnData,clinData, factor,log2=TRUE){
#   #ExpnData is the expression data with patients as columns, genes as rows. 
#   #clindata  has the factor data, such as mutation status pos,neg. Patient USI as rownames
#   #cols specifies the numeric columns with expn values. 
#   #factor is the name of the factor column 
#   #log2 indicates wether the input expn data matrix is already log2 or not. 
#   library(genefilter)
#   
#   Ngenes <- length(which(colnames(expnData) %in% rownames(clinData))) - 1 #ensure that number of genes is less than #samples
#   
#   topVarGenes <- rownames(expnData[order(rowVars(as.matrix(expnData)),decreasing=TRUE), ] %>%
#                             .[1:Ngenes, ])
#   
#   pca_df <- merge_CDE_Expn(clinData,expnData, topVarGenes) #merge causes loss of rownames. 
#   
#   topVarGenes <- gsub("\\-", "\\.", topVarGenes) #remove hyphens to match colnames
#   topVarGenes <- gsub("(^[0-9].+)", "X\\1", topVarGenes) #add an X to genes names that begin with numbers
#   
#   
#   if (log2 == FALSE){
#     expn <- log2(pca_df[,intersect(colnames(pca_df),topVarGenes)] + 0.01)
#   }else if (log2 == TRUE){
#     expn = pca_df[,intersect(colnames(pca_df),topVarGenes)]
#   }
#   
#   pca.scaled <- scale(expn)
#   
#   pca <- princomp(pca.scaled,cor = T, scores = T)
#   scores <- as.data.frame(unclass(pca$scores))
#   percVar <- round((pca$sdev^2)/sum(pca$sdev^2)*100, digits = 2)
#   
#   pca_plot <- ggplot(scores, aes(x=scores$Comp.1, y=scores$Comp.2))
#   
#   pca_plot <- pca_plot + geom_point(aes(color=factor(pca_df[,factor]))) +
#     theme_bw() +
#     labs(title = "PCA of TARGET AML: Most Varied Genes",
#          x = paste("principal Component 1:", percVar[1], "% Variance", sep=" "),
#          y = paste("principal Component 2:  ", percVar[2], "% Variance", sep=" "))
#   
#   # scree_plot <-  plot(percVar[1:10], xlab = "Principal Component",
#   #                       ylab = "Percentage of Variance Explained",
#   #                       type = "b", pch = 16,
#   #                       main = "Percent Variance Exlpained in PCA analysis",
#   #                       col = "dark blue")
#   
#   list <- list(pca_df,pca, scores, pca_plot, topVarGenes)
#   names(list) <- c("pca_df", "pca", "scores", "pca_plot", "topVarGenes")
#   
#   return(list)
# }



