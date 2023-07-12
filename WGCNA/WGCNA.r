### Display the current working directory
getwd();

# #Load the package
# if (!require("WGCNA")) {
#   source("http://bioconductor.org/biocLite.R") 
#   biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
#   install.packages("WGCNA", dependencies = TRUE,repos="http://cran.rstudio.com/")
#   library(WGCNA)
# }
# if (!require("flashClust")) {
#   install.packages("flashClust", dependencies = TRUE,repos="http://cran.rstudio.com/")
#   library(flashClust)
# }
library(WGCNA)
library(flashClust)
library(limma)
library(zip)

###################################
# --- Parameters ------------------
###################################

args <- commandArgs()

IN_expression <- args[6]
IN_trait <- args[7]

if(args[8] == "") {
  my_minModuleSize <- as.numeric(50)
} else {
  my_minModuleSize = as.numeric(args[8])
}

OutputDir <- args[9]

temp_sft_TEST = args[10]
if(temp_sft_TEST == 'NULL'){
  sft_TEST = NULL
}else{
  if( !is.na(suppressWarnings(as.numeric(temp_sft_TEST))) ){ # if is is able to convert to numeric
    sft_TEST <- as.numeric(temp_sft_TEST)
  }else{
    print("Error: scale_free_model_fit must be NULL or number ranging from 0 to 1")
    stop()
  }
}

mod_trait_corr <- as.character(args[11])
quantile_or_raw <- as.character(args[12])
signed_or_unsigned <- as.character(args[13])

if(args[14] == "") {
  my_deepSplit <- as.numeric(1) # Default
} else {
  my_deepSplit = as.numeric(args[14])
}

if(args[15] == 'T'){
  my_pamRespecDendro = TRUE
}else if(args[15] == 'F'){
  my_pamRespecDendro = FALSE
}else{
  print("Error: pamRespecDendro must be 'T' or 'F'")
  stop()
}

print(paste0("min_mod_size=", my_minModuleSize))
print(paste0("OutputDir=", OutputDir))
print(paste0("scale-free_topology_fit=", sft_TEST))
print(paste0("mod_trait_corr = ", mod_trait_corr))
print(paste0("quantile_or_raw = ", quantile_or_raw))
print(paste0("signed_or_unsigned = ", signed_or_unsigned))
print(paste0("deepSplit = ", my_deepSplit))
print(paste0("pamRespecDendro = ", my_pamRespecDendro))

if(file.exists(OutputDir)==FALSE){
    dir.create(OutputDir)
} else{
    old_files <- list.files(OutputDir)
    print(old_files)
    if(length(old_files) != 0){
        for(i in 1:length(old_files)){
            file.remove(paste(OutputDir, old_files[i], sep="/"))
        }
    }
}

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

print("#####################################################")


################################################################################
# --- Data Import --------------------------------------------------------------
################################################################################

# Read Expression Dataset
femData = read.table(file=IN_expression, sep="\t", 
                     header=T, row.names=1, check.names=F);

if(quantile_or_raw == "quantile"){
    print("quantile_or_raw == quantile")
    femData <- normalizeBetweenArrays(as.matrix(femData))
} else if(quantile_or_raw == "raw"){
    print("quantile_or_raw == raw")    
    # do nothing
} else{
    print("Error. 'quantile' or 'raw' must be chosen as 'quantile_or_raw'")
    stop()
}

write.table(femData, 
            file=(paste(OutputDir, "expression_all.txt", sep="/")),
            col.names=T,sep="\t",row.names=T)

# Take a quick look at what is in the data set:
print(paste0('Expression data - ', 
             'nrow: ', nrow(femData), ', ncol: ', ncol(femData)))
print('Expression data colnames:')
colnames(femData);

datExpr0 <- femData

# We first check for genes and samples with too many missing values:
gsg = goodSamplesGenes(datExpr0, verbose = 3);

if (gsg$allOK) {
  print('goodSamplesGene: TRUE')
}else{
  print('goodSamplesGene: FALSE')
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0){
    printFlush(paste("Removing genes:", 
                     paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  }
  if (sum(!gsg$goodSamples)>0){
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  }
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

#' Next we cluster the samples (in contrast to clustering genes that will come later) 
#' to see if there are any obvious outliers. We use the function flashClust that 
#' provides faster hierarchical clustering than the standard function hclust.
sampleTree = flashClust(dist(t(datExpr0)), method = "average");

#' Plot the sample tree: Open a graphic output window of size 12 by 9 inches
#' The user should change the dimensions if the window is too large or too small.
print('Making sampleClustering_1.pdf')
pdf(paste(OutputDir, "sampleClustering_1.pdf", sep="/"), width = 25, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

print("############## done (1) Import expression, flashClust by Sample ########")

#' It appears there is one outlier (sample F2_221, see Fig. 1). One can remove it by hand, or use an automatic approach.
#' Choose a height cut that will remove the oending sample, say 15 (the red line in the plot), and use a branch cut at
#' that height.(In the example analyses, we did not remove outliers.)

# Chunk of Comment Out --------------------------------------------------
# Plot a line to show the cut
#abline(h = 15, col = "red");
# Determine cluster under the line
#clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
#table(clust)
# clust 1 contains the samples we want to keep.
#keepSamples = (clust==1)
#datExpr = datExpr0[keepSamples, ]
#nGenes = ncol(datExpr)
#nSamples = nrow(datExpr)
# -----------------------------------------------------------------------
datExpr = datExpr0

# We now read in the trait data and match the samples for which they were measured to the expression samples.
traitData = read.table(file = IN_trait, sep="\t",
                       row.names=1, header=T, check.names=F);
#dim(traitData)
#names(traitData)

# if(ncol(traitData) >= 2) {
#   noteda <- "The number of traits should not be greater than 1."
#   write.table(noteda,"../Output/Error.txt",row.names=F,sep="",col.names=F)
#   q("no")
# }
#print (traitData)

# Remove columns that hold information we do not need.
allTraits = traitData
print(paste0('Trait data - ', 
             'nrow: ', nrow(allTraits), ', ncol: ', ncol(allTraits)))
print('Trait data colnames:')
names(allTraits);

datTraits <- allTraits

print("############## done (2) Import trait ##################################")


#' We now have the expression data in the variable datExpr, and the corresponding clinical traits in the variable
#' datTraits. Before we continue with network construction and module detection, we visualize how the clinical traits
#' relate to the sample dendrogram.
sampleTree2 = flashClust(dist(t(datExpr)), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
print('Making sampleClustering_2.pdf')
pdf(paste(OutputDir, "sampleClustering_2.pdf", sep="/"), width = 25, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()
# save(datExpr, datTraits, file = "01-dataInput.RData")

print("############## done (3) Plot Sample Dendro & Trait Heatmap ############")


################################################################################
# --- Network Construction -----------------------------------------------------
################################################################################

# The variable lnames contains the names of loaded variables.
# lnames = load(file = "01-dataInput.RData");
# lnames

#' Choosing the soft-thresholding power: analysis of network topology
#' Choose a set of soft-thresholding powers
powers = c(c(1:20), seq(from = 25, to=100, by=5))

# Call the network topology analysis function

if(signed_or_unsigned == "signed"){
    print("signed_or_unsigned == signed")
    similarity <- cor(t(datExpr))
    similarity <- (1+similarity)/2
} else if(signed_or_unsigned == "unsigned"){
    print("signed_or_unsigned == unsigned")    
    similarity <- abs(cor(t(datExpr)))    
} else{
    print("Error. 'signed' or 'unsigned' must be chosen as 'signed_or_unsigned'")
    stop()
}

# write.table(as.data.frame(similarity), file="test_cor.txt", sep="\t",
#             quote=F, row.names=F, col.names=colnames(similarity))

sft = pickSoftThreshold.fromSimilarity(similarity, 
                                       powerVector = powers, 
                                       verbose = 5)
#sft = pickSoftThreshold(t(datExpr), powerVector = powers, verbose = 5)

# Plot the results
print('Making SF_cut.pdf')
pdf(paste(OutputDir, "SF_cut.pdf", sep="/"), width = 15, height = 7.5);
par(mfrow = c(1,2));
cex1 = 0.8;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",
     type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red");

# This line corresponds to using an R^2 cut-off of h
if(!is.null(sft_TEST)){ 
  abline(h=sft_TEST, col="red")
}

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], 
     labels=powers, cex=cex1, col="red")
dev.off()

if(!is.null(sft_TEST)){ # find the best soft power
  numSBIsftInfo <- length(sft$fitIndices[,1])
  SBIsftRR2 <- 0
  SBIsftPWR <- 0
  for(i in 1:numSBIsftInfo) {
    if(-sign(sft$fitIndices[i,3])*sft$fitIndices[i,2] < sft_TEST) {
      
    } else {
      print(-sign(sft$fitIndices[i,3])*sft$fitIndices[i,2])
      if(SBIsftRR2 > 0) {
        
      } else {
        SBIsftRR2 <- -sign(sft$fitIndices[i,3])*sft$fitIndices[i,2]
        SBIsftPWR <- sft$fitIndices[i,1]
      }
    }
  }
  print(paste0("soft power =", SBIsftPWR))
}else{  # use the fixed soft power
  if(ncol(femData) < 20){
    if(signed_or_unsigned == "signed"){
      SBIsftPWR = 18
    }else{
      SBIsftPWR = 9
    }
  }else if(ncol(femData) > 20 && ncol(femData) <= 30){
    if(signed_or_unsigned == "signed"){
      SBIsftPWR = 16
    }else{
      SBIsftPWR = 8
    }
  }else if(ncol(femData) > 30 && ncol(femData) <= 40){
    if(signed_or_unsigned == "signed"){
      SBIsftPWR = 14
    }else{
      SBIsftPWR = 7
    }
  }else if(ncol(femData) > 40){
    if(signed_or_unsigned == "signed"){
      SBIsftPWR = 12
    }else{
      SBIsftPWR = 6
    }
  }
  print(paste0("scale-free_topology_fit == NULL -> HARD-coded soft power: ", SBIsftPWR, 
               ' (', signed_or_unsigned, ' & sample num: ', ncol(femData), ')'))
}

if(SBIsftPWR == 0) {
  SBIsftPWR <- 10
  print("soft power == 0 -> HARD-coded soft power: 10")
}

print("############## done (4) Calc Similarity & pickSoftThreshold ###########")		


# We now calculate the adjacencies, using the soft thresholding power SBIsftPWR:
softPower = SBIsftPWR;
similarity = adjacency.fromSimilarity(similarity, power = softPower);

##--##write.table(rownames(similarity),file="ADJ.txt",sep="\t")
##--##write.table(colnames(similarity),file="ADJC.txt",sep="\t")

# Turn adjacency into topological overlap
similarity = TOMsimilarity(similarity);

#out_similarity = similarity
##--##write.table(out_similarity,file="similarity.txt",sep="\t",row.names=T,col.names=T)

similarity = 1-similarity

print("############## done (5) Calc Adjacency & TOM ##########################")		


#' Clustering using TOM
#' Call the hierarchical clustering function
geneTree = flashClust(as.dist(similarity), method = "average");

# Plot the resulting clustering tree (dendrogram)
print('Making Gene_Clustering_1.pdf')
pdf(paste(OutputDir, "Gene_Clustering_1.pdf", sep="/"), width = 15, height = 7.5);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(geneTree, xlab="", sub="", 
     main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()

print("############## done (6) Clustering using TOM ##########################")		

print(paste0("deepSplit = ", my_deepSplit))
print(paste0("pamRespecDendro = ", my_pamRespecDendro))

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = my_minModuleSize;

# Module identification using dynamic tree cut:

dynamicMods = 
  cutreeDynamic(dendro = geneTree, distM = similarity,
                minClusterSize = minModuleSize,
                deepSplit = my_deepSplit,
                pamRespectsDendro = my_pamRespecDendro);
table(dynamicMods)


print("############## done (7) Cut-tree ######################################")		


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
print('Making Gene_Clustering_2.pdf')
pdf(paste(OutputDir, "Gene_Clustering_2.pdf", sep="/"), width = 15, height = 10.0);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

#' Merging of modules whose expression profiles are very similar
#' The Dynamic Tree Cut may identify modules whose expression profiles are very similar. It may be prudent to merge
#' such modules since their genes are highly co-expressed. To quantify co-expression similarity of entire modules, we
#' calculate their eigengenes and cluster them on their correlation:

# Calculate eigengenes
firstPCeigen <- moduleEigengenes(t(datExpr), dynamicColors, returnValidOnly=FALSE)
write.table(firstPCeigen$eigengenes, 
            file=paste(OutputDir, "MEs_1stPC.txt", sep="/"), 
            sep="\t", row.names=T)

MElist <- chooseTopHubInEachModule(t(datExpr), dynamicColors)
MEvalue <- t(datExpr[MElist, , drop=FALSE])
#MEvalue <- cbind(MEvalue,t(t((1:ncol(datExpr)))))
#colnames(MEvalue) <- append(paste(names(MElist),sep=""),
#                            "grey", after=length(MElist))
colnames(MEvalue) <- paste(names(MElist), sep="")
# ##--##write.table(MEvalue,file="testME.txt",sep="\t")

MEs = MEvalue
print('Eigengene Expression')
print(head(MEs)) 


# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
#print(MEDiss)

if (ncol(MEvalue) > 2) {
  # Cluster module eigengenes
  METree = flashClust(as.dist(MEDiss), method = "average");
  
  # Plot the result
  print('Making Module_Clustering_1.pdf')
  pdf(paste(OutputDir, "Module_Clustering_1.pdf", sep="/"), width = 5, height = 3.0);
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
  dev.off()
  
  #We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge
  #MEDissThres = 0.25
  #Plot the cut line into the dendrogram
  #abline(h=MEDissThres, col = "red")
}

print("############## done (8) Calc Eigengene and Module Similarity###########")		


# Call an automatic merging function
# print(t(datExpr))
# MEDissThres = 0.5
# merge = mergeCloseModules(t(datExpr), dynamicColors, 
#                           MEs = MEs, cutHeight = MEDissThres, verbose = 3)
# merge = mergeCloseModules(t(datExpr), dynamicColors, verbose = 3)

# The merged module colors
# mergedColors = merge$colors;
mergedColors = dynamicColors;

# Eigengenes of the new merged modules:
#mergedMEs = merge$newMEs;
mergedMEs = MEs;

print('Making Plots_geneDendrogram.pdf')
pdf(file = paste(OutputDir, "Plots_geneDendrogram.pdf", sep="/"), wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"), 
                    dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, 
                    guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

#Save module colors and labels for use in subsequent parts
#save(MEs, moduleLabels, moduleColors, geneTree, 
#     file = "02-networkConstruction-stepByStep.RData")

print("############## done (9) ###############################################")		


##################################################################################
# ---- Relating modules to external information and identifying important genes --
##################################################################################

#lnames = load(file = "02-networkConstruction-stepByStep.RData");
#lnames

#' Quantifying module{trait associations
#' In this analysis we would like to identify modules that are signficantly associated with the measured clinical traits.
#' Since we already have a summary profile (eigengene) for each module, we simply correlate eigengenes with external
#' traits and look for the most significant associations:

# Define numbers of genes and samples
nGenes = ncol(t(datExpr));
nSamples = nrow(t(datExpr));

# Recalculate MEs with color labels
MElist0 <- chooseTopHubInEachModule(t(datExpr), moduleColors)
MEvalue0 <- t(datExpr[MElist0,,drop=FALSE])
# MEvalue0 <- cbind(MEvalue0,t(t((1:ncol(datExpr)))))
# colnames(MEvalue0) <- append(paste(names(MElist0),sep=""),"grey",after=length(MElist0))
colnames(MEvalue0) <- paste(names(MElist0),sep="")
MEs0 <- MEvalue0

##--##write.table(MEs0,file="testtesttest.txt",sep="\t")

MEs = orderMEs(MEs0)

print("############## !! ###########")		
#!!!!!my-----
print(MEs)
#!!!!!my-----
print("############## !!! ###########")		

moduleTraitCor = cor(MEs, datTraits, use = "p", method = mod_trait_corr);
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#Pthresfold = moduleTraitPvalue
#Pthresfold[Pthresfold[,] > 0.05] <- 2
#Pthresfold[Pthresfold[,] <= 0.05] <- 1
#Pthresfold[Pthresfold[,] == 2] <- 0
#Pthresfold <- Pthresfold*moduleTraitCor

moduleCorPva <- cbind(moduleTraitCor, moduleTraitPvalue)
hdCorPva = c(paste("Correlation", t(colnames(traitData))),
             paste("pvalue", t(colnames(traitData))))
hdCorPva = c("module",hdCorPva)

write.table(t(hdCorPva), 
            file= paste(OutputDir, "Module_trait_relationships.txt", sep="/"),
            row.names=F, col.names=F, sep="\t", quote=F)
write.table(moduleCorPva,
            file= paste(OutputDir, "Module_trait_relationships.txt", sep="/"),
            row.names=T, append=TRUE, col.names=F, sep="\t",quote=F)
write.table(moduleTraitPvalue,
            file= paste(OutputDir, "ModPvalue.txt", sep="/"), sep="\t")
write.table(moduleTraitCor,
            file=paste(OutputDir, "ModCol.txt", sep="/"), sep="\t")
#write.table(Pthresfold,file="../Output/ModCThres.txt",sep="\t")

# Chunk of Comment Out ---------------------------------------------------------
# #sizeGrWindow(10,6)
# # Will display correlations and their p-values
# textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
# signif(moduleTraitPvalue, 1), ")", sep = "");
# dim(textMatrix) = dim(moduleTraitCor)
# number_traits <- ncol(traitData)
# number_modules <- nrow(moduleTraitCor)
# height_mapda <- number_modules*1.5 + 2
# width_mapda <- number_traits*2.0 + 2
# 
# pdf(file = "../Output/Module_trait_relationships_visualization.pdf", wi = width_mapda, he = height_mapda)
# par(mar = c(6, 8.5, 3, 3));
# # Display the correlation values within a heatmap plot
# labeledHeatmap(Matrix = moduleTraitCor,
# xLabels = names(datTraits),
# #yLabels = names(MEs),
# yLabels = paste(names(MEs),"\n","Vipul",sep=""),
# ySymbols = names(MEs),
# colorLabels = FALSE,
# colors = blueWhiteRed(50),
# textMatrix = textMatrix,
# setStdMargins = FALSE,
# cex.text = 0.5,
# zlim = c(-1,1),
# main = paste("Module-trait relationships"))
# dev.off()
# ----------------------------------------------------------------------------
print("############## done (10) ###########")		


#' Gene relationship to trait and important modules: Gene Significance and Module Membership
#' We quantify associations of individual genes with our trait of interest (weight) by dedining Gene Signdicance GS as
#' (the absolute value of) the correlation between the gene and the trait. For each module, we also ddine a quantitative
#' measure of module membership MM as the correlation of the module eigengene and the gene expression prdile. This
#' allows us to quantify the similarity of all genes on the array to every module.

#' Define variable weight containing the weight column of datTrait
#' In this line, please specify the phenotype of interest (in example, we focused on "Braak Severity")######
weight = as.data.frame(datTraits[,1]); 
names(weight) = "weight"

# Names (colors) of the modules
#modNames = substring(names(MEs), 3)
modNames = names(MEs)
print(modNames)

geneModuleMembership = as.data.frame(cor(t(datExpr), MEs, use = "p"));
names(geneModuleMembership) = paste("MM", modNames, sep="");

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(t(datExpr), weight, use = "p"));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");

GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

TraitSignificance = as.data.frame(cor(t(datExpr),datTraits,use="p"));
names(TraitSignificance) = paste("GS.",names(datTraits),sep="")

TSPvalue = as.data.frame(corPvalueStudent(as.matrix(TraitSignificance), nSamples))
names(TSPvalue) = paste("p.GS.",names(datTraits),sep="");

#' Intramodular analysis: identifying genes with high GS and MM
#' Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module
#' membership in interesting modules. As an example, we look at the brown module that has the highest association
#' with weight. We plot a scatterplot of Gene Significance vs. Module Membership in the brown module:

sigmodsIDs <- which(moduleTraitPvalue[,1] < 1.00)
modNames <- modNames[sigmodsIDs] 

numSBImodules <- length(modNames)

for(i in 1:numSBImodules) {
	if(modNames[i] == "grey") {
	  
	} else {
	  #####in this line, please specify module of interest, in example, we focused "BLUE" module #####
	  module = modNames[i]
	  column = match(module, modNames);
  	moduleGenes = moduleColors==module;
	  #sizeGrWindow(7, 7);
	  outputSBImod <- paste(OutputDir, "/GS_vs_MM_",module,".pdf",sep="")
	  
	  pdf(file = outputSBImod, wi = 6, he = 6)
	  par(mfrow = c(1,1));
	  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
	                     abs(geneTraitSignificance[moduleGenes, 1]),
	                     xlab = paste("Module Membership in", module, "module"),
	                     ylab = "Gene significance for phenotype",
	                     main = paste("Module membership vs. gene significance\n"),
	                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
	  dev.off()
	}
}

print("############## done (11) ###########")		


#' Summary output of network analysis results
#' We have found modules with high association with our trait of interest, and have identigied their central players by
#' the Module Membership measure. We now merge this statistical information with gene annotation and write out agile that summarizes the most important results and can be inspected in standard spreadsheet software such as MS
#' Excel or Open Oce Calc. Our expression data are only annotated by probe ID names: the command

analysesrt <- cbind(moduleColors, TraitSignificance, TSPvalue, 
                    geneModuleMembership, MMPvalue)

# Chunk of Comment Out --------------------------------------------------
#print(analysesrt)
#analysesrtOR <- order(analysesrt$moduleColors, -abs(analysesrt$GS.weight));
#analysesrt = analysesrt[analysesrtOR,]

#names(t(datExpr))
#names(t(datExpr))[moduleColors=="blue"]
#annot = read.csv(file = "GeneAnnotation.csv");
#dim(annot)
#names(annot)
#probes = names(t(datExpr))
#probes2annot = match(probes, annot$substanceBXH)

## The following is the number or probes without annotation:
#sum(is.na(probes2annot))
## Should return 0.
## Create the starting data frame
#geneInfo0 = data.frame(substanceBXH = probes,
#geneSymbol = annot$gene_symbol[probes2annot],
#LocusLinkID = annot$LocusLinkID[probes2annot],
#moduleColor = moduleColors,
#geneTraitSignificance,
#GSPvalue)
## Order modules by their significance for weight
#modOrder = order(-abs(cor(MEs, weight, use = "p")));
## Add module membership information in the chosen order
#for (mod in 1:ncol(geneModuleMembership))
#{
#oldNames = names(geneInfo0)
#geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
#MMPvalue[, modOrder[mod]]);
#names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
#paste("p.MM.", modNames[modOrder[mod]], sep=""))
#}
## Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
#geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
#geneInfo = geneInfo0[geneOrder, ]
#write.csv(geneInfo, file = "geneInfo.csv")
# -------------------------------------------------------------------------

# print("OK test")
# print("OK test2")
genesT <- rownames(datExpr)

expressionrt1 <- datExpr0 

rownames(analysesrt) <- chartr(".","-",rownames(analysesrt))

write.table(t(c("",colnames(analysesrt))), 
            file=paste(OutputDir, "Analyses_Result.txt", sep="/"),
            sep="\t", row.names=F, col.names=F, quote=F)

write.table(analysesrt,
            file=paste(OutputDir, "Analyses_Result.txt", sep="/"),
            sep="\t", row.names=T, col.names=F, quote=F, append = T)

numSBImodules <- length(modNames)

ExprFilenames <- c()
for(i in 1:numSBImodules) {
  temp_module <- modNames[i]
  moduleG <- which(moduleColors == temp_module)
  valueP <- paste(OutputDir, "/p.MM",temp_module,sep="")
  
  outputFa <- paste(OutputDir, "/expression_data_",temp_module,".txt",sep="")
  outputF2 <- paste(OutputDir, "/Analyses_result_",temp_module,".txt",sep="")
  ExprFilenames <- append(ExprFilenames ,c(outputFa))
  
  temp_expression1 <- expressionrt1[moduleG,] 
  temp_analyses_result <- analysesrt[moduleG,]
  
  rownames(temp_expression1) <- chartr(".","-",rownames(temp_expression1))
  header_MODMOD <- colnames(temp_expression1)
  header_MODMOD <- c("genes",header_MODMOD)
  write.table(t(header_MODMOD), file=outputFa,
              row.names=F,sep="\t",col.names=F,quote=F)
  write.table(temp_expression1, file=outputFa,
              sep="\t",col.names=F,append=TRUE,quote=F)
  
  rownames(temp_analyses_result) <- chartr(".", "-", rownames(temp_analyses_result))
  
  write.table(t(c("", colnames(temp_analyses_result))),
              file=outputF2, sep="\t", row.names=F, col.names=F, quote=F)
  
  write.table(temp_analyses_result,
              file=outputF2, sep="\t", row.names=T, col.names=F, quote=F, append=T)
  
  # Top # of genes is 20 genes#
  genesLTS1 <- rownames(temp_analyses_result)
  print(genesLTS1)
  print(moduleG)
  
  # Chunk of Comment Out -----------------------------------
  #geneLSTexp1 <- temp_expression1[genesLTS1,]
  #write.table(geneLSTexp1,file=outputF5a,sep="\t",col.names=F)
  #geneLSTexp2 <- temp_expression2[genesLTS1,]
  #write.table(geneLSTexp2,file=outputF5b,sep="\t",col.names=F)
  #genesTOP50 <- rownames(geneLSTexp1)
  #write(genesTOP50,file=outputF7)
  #print("ckckck0")
  #rownames(adjacency) <- chartr(".","-",rownames(femdata1))
  #colnames(adjacency) <- chartr(".","-",rownames(femdata1))
  #print("ckckck0-1")
  #modTOM <- adjacency[genesLTS1,genesLTS1]
  #print("ckckck1")
  #diag(TOM) <- 0
  #print("ckckck01")
  #write.table(modTOM,file=outputF6,sep="\t",row.names=F)
  #print("ckckck12")
  # ---------------------------------------------------------
}

zipr(paste(OutputDir, "expression_files.zip", sep="/"), ExprFilenames)

print("############## done (12) ###########")	


textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
number_traits <- ncol(traitData)
number_modules <- nrow(moduleTraitCor)
height_mapda <- number_modules*1.5 + 2
width_mapda <- number_traits*2.0 + 2

num_genes <- array()
for(i in 1:length(names(MEs))){
  num_genes[i] <- nrow(subset(analysesrt , analysesrt$moduleColors %in% names(MEs)[i]))
}

# Display the correlation values within a heatmap plot
pdf(file = paste(OutputDir, "Module_trait_relationships_visualization.pdf", sep="/"), 
    wi = width_mapda, he = height_mapda)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               #yLabels = names(MEs),
               yLabels = paste(names(MEs),"\n#genes=",(num_genes),sep=""),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

print("############## done (13) ###########")	


if(F){
  #file.remove("Rplots.pdf")
  #file.remove("Rplots1.pdf")
  
  # #---------------------------
  # ######################################
  # ###Network output for visualization###
  # ######################################
  # library(WGCNA) # Load the WGCNA package
  # options(stringsAsFactors = FALSE);# The following setting is important, do not omit.
  # lnames = load(file = "01-dataInput.RData");# Load the expression and trait data saved in the first part
  # lnames = load(file = "02-networkConstruction-stepByStep.RData");# Load network data saved in the second part.
  # #softPower = 30
  # TOM = TOMsimilarityFromExpr(t(datExpr), power = softPower);# Recalculate topological overlap if needed
  # modules = modNames;
  # numSBImodules <- length(modules)
  # for(i in 1:numSBImodules) {
  #       module <- modules[i];
  #       inModule = is.finite(match(moduleColors, module));
  #       genes <- rownames(datExpr);
  #       genesModule <- genes[inModule];
  #       genesModule <- chartr(".","-",genesModule);
  #       modTOM = TOM[inModule, inModule];
  #       dimnames(modTOM) <- list(genesModule,genesModule);
  #       modTOMtemp <- modTOM;
  #       diag(modTOMtemp) <- 0;
  #       outputNET <- paste(module,"_WGCNA",".txt",sep="");
  #       write.table(modTOMtemp,file=outputNET,row.names=F,sep="\t");
  
  #       # Export the network into edge and node list files Cytoscape can read
  #       cyt = exportNetworkToCytoscape(modTOM,edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),weighted = TRUE,threshold = 0.02,nodeNames = genesModule,nodeAttr = moduleColors[inModule]);}
  # #---------------------------
}

