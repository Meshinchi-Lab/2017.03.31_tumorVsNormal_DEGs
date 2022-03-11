#Jenny Smith

#Aug 16, 2017

#Differential Expression Analysis pipeline

#Purpose: Given matrix of RNA-seq raw counts and clinical annoations, create a series of differential expression Analyses for two group comparisions. 
setwd(file.path(TARGET,"RNA/miRNAseq/2017.03.31_tumorVsNormal_DEGs"))

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

calcDEMir <- function(expData, g1, g2, logCPM=NULL,Trend=TRUE) {
  # expnData is a matrix or data frame with the raw counts. Patient IDs as colnames, genes as rownames
  # g1,g2 are the patient IDs for each group
  #logCPM is whether to log2 CPM normalize the raw counts.
  
  library(limma)
  library(edgeR)
  
  dge <- DGEList(counts = expData) #miRNA input is already in RPM normalized (reads per million)
  keep.dge <- rowSums(dge$counts >= 1) >= (0.05*ncol(expData))#5% of dataset has RPM of at least 1 for a gene
  dge <- dge[keep.dge,] #subset for those genes with RPM >= 1 per gene
  
  if (is.null(logCPM)){
    dge <- calcNormFactors(dge) #calculate  the TMM normalization factors
    return(dge)
    
  }else if (logCPM==TRUE){
    dge <- calcNormFactors(dge) #calculate the normalization factors
    NormFactors <- dge$samples
    dge <- cpm(dge, log=TRUE, prior.count = 1) #log2 CPM transformation.
    
  }else if (logCPM==FALSE){
    dge <- apply(dge$counts, 2, function(x) log2(x + 1))
    NormFactors <- "None"
  }
  
  
  designMatrix <- matrix(0, nrow=dim(dge)[2] , ncol=2)
  colnames(designMatrix) <- c("g1","g2")
  rownames(designMatrix) <- colnames(dge)
  designMatrix[g1, 1] <- 1
  designMatrix[g2, 2] <- 1
  
  fit <- lmFit(dge,designMatrix)
  tmp <- paste("g1","g2",sep="-") #contrast is ~ log2(mean(g1)) - log2(mean(g2)) per gene
  
  cont.matrix <- makeContrasts(contrasts=tmp,
                               levels=designMatrix)
  
  fit2<-contrasts.fit(fit, cont.matrix)
  fit2<-eBayes(fit2, trend=Trend) #added on 6/4/2017 by JS
  DE<-topTable(fit2,adjust.method="BH",sort.by="P",
               number=20000,p.value=0.05, lfc=1)
  
  
  
  list <- list(dge,NormFactors, designMatrix,fit2, DE)
  names(list) <- c("dge","NormFactors", "design","eBayesFit", "DE")
  
  return(list)
}

dendrograms <- function(df, pheno, genelist, method){
  #df with TPM or RPKM normalized data, patient IDs are column names and rows are genes. 
  #pheno is a character vector with patient IDs as names, and the status for each in each group (eg pos,neg)
  #genelist is a character vector with the genes of interest
  require(edgeR)
  library(dendextend)
  
  #ensure correct order, drop rows with nas just in case
  df <- df[, intersect(names(pheno), colnames(df))] 
  
  df <- df[which(rownames(df) %in% genelist), ] #subset the matrix to genes of interest
  
  # names <- rownames(df)
  df <- as.data.frame(apply(df, 2, function(x) log2(x + 0.01))) #log2 transform counts
  # TMMCPM10 <- as.data.frame(apply(TMMCPM, 2, function(x) log10(x + 0.01))) #log10 transform counts
  
  
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


plotPCoA <- function(expnData,clinData, geneList=NULL, factor){
  #factor is the name of the factor column 
  library(vegan)
  
  #Ensure correct order of patients in both datasets
  expnData <- expnData[ ,intersect(rownames(clinData), colnames(expnData))] 
  clinData <- clinData[intersect(rownames(clinData),colnames(expnData)), ]
  
  if (! is.null(geneList)){
    expnData <- t(expnData[geneList, ])
  }else{
    expnData <- t(expnData) #note: must remove all zero count genes or will  fail on an error
  }
  
  PCoA <- capscale(expnData ~ 1, distance = "bray", add=TRUE)
  scores <- as.data.frame(scores(PCoA, display="sites"))
  
  p <- ggplot(scores, aes(x=MDS1, MDS2)) +
    geom_point(aes(color=clinData[,factor])) +
    theme_numX +
    labs(title="")
  
  # aov <- aov(scores$MDS1 ~ df[,factor]) #NOTE: This is only valid for balanced experimental designs! Equal # of obs in each factor level.
  
  list <- list(df, PCoA,scores, p)
  names(list) <- c("MDS_df","PCoA","scores","plot")
  
  return(list)
}






twoGroups_DEMirs <- function(expnData, clinData, col, ref,logCPM=FALSE,BM=FALSE){
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
  # print(name)
  clinData <- rmUnknowns(clinData, col)
  groups <- GroupIDs(clinData, col) #list of patient IDs, one for each group
  
  #Define Groups to compare based on group IDs from clinical data. Intersect with expression matrix to subset. 
  if (BM == TRUE){
    BM <- grep("^BM", colnames(expnData), value = TRUE)
    GroupB <- BM #select the reference group 
    GroupA <- groups[[which(names(groups) != ref)]] %>% intersect(. , colnames(expnData)) #the second group (mutant, AML, treated, etc)
  }else{
    GroupB <- groups[[ref]] %>% intersect(. , colnames(expnData)) #select the reference group (eg No, normal, wt, control, etc.) Must be a character(level) from the column of clinData selected. 
    GroupA <- groups[[which(names(groups) != ref)]] %>% intersect(. , colnames(expnData)) #the second group (mutant, AML, treated, etc)
  }
  
  if (any(lapply(list(GroupA,GroupB), length) < 3)){
    list <- list(expnData, clinData, GroupA,GroupB)
    names(list) <- c("InputExpnData", "InputClinData", "CompGroup", "RefGroup")
    return(list)
  }
  
  phenoVector <- phenoVectors(GroupA, GroupB)
  
  if (identical(GroupB,BM)){
    clinData <- as.data.frame(phenoVector) %>% set_colnames(., "Group")
    col <- "Group"
  }else{
    clinData = clinData
  }
  
  #subset and order the expression values dataframe.
  expnData <- expnData[,match(c(GroupA, GroupB), colnames(expnData))] #mutant, then WT
  if (any(is.na(expnData))){print("NAs Introduced. Check Rownames and colnames of inputs")}

  # Calculate Differential Expression
  print(c("logCPM", logCPM))
  DE <- calcDEMir(expnData,GroupA, GroupB, logCPM=logCPM, Trend = TRUE) #mutant - wild type. logCPM the counts
  DE$DE$FoldChange <- logratio2foldchange(DE$DE$logFC)
  # NOTE: I included a more stringent filter here, so 5% of samples must have logCPM of greater than 1 for inclusion in analysis
  # this usually results in ~18,000 genes included in each analysis.

  if (nrow(DE$DE) < 5){
    MDS <- plotPCoA(2^DE$dge, clinData, rownames(expnData),col)

    list <- list(clinData, phenoVector, expnData, DE, PCA, MDS)
    names(list) <- c("InputClinData", "phenovector", "InputExpnMatrix", "DE","PCA", "MDS")
    return(list)

  }else{


    #Unsupervised Heirachrach clustering
    dends_DE <- dendrograms(expnData, phenoVector, rownames(DE$DE), method="ward.D2") #dendrograms based on all differentially expressed genes.
    colorBar <- ifelse(phenoVector == "GroupB", "black", "firebrick")
    # title <- paste(variantName, dataName, sep=" ")
    title <- variantName
    heatmap <- basicHeatmap(dends_DE$expnData, dends_DE$c2, dends_DE$c1,colorBar, title=title)

    #Unconstrained Cluster Analysis
    # PCoA
    genes <- rownames(DE$dge)
    MDS <- plotPCoA(2^DE$dge,clinData,genes,col) #non-log CPM normalized values

    #return the objects
    list <- list(clinData, phenoVector, expnData, DE, genes, dends_DE, heatmap, MDS)
    names(list) <- c("InputClinData", "phenovector", "InputExpnMatrix", "DE", "topVargenes", "dends_DE", "Heatmap", "MDS")


    return(list)
  }
}



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



