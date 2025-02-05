#!/usr/bin/env Rscript

pdf(NULL)

library("optparse")

###################################
#######  PARSING OPTIONS #########
###################################
option_list = list(
    make_option(c("-p", "--phenodata"), type="character", default=NULL, help="Phenodata file path", metavar="character"),
    make_option(c("-o", "--outdir"), type="character", default="out/", help="Output directory (must exist at launch time) [default= %default]", metavar="character"),
    make_option(c("-i", "--gene_counts"), type="character", default="gene_count_matrix.csv", help="input DESeq2 count file [default= %default]", metavar="character"),
    make_option(c("-v", "--var_to_test"), type="character", help="the variable/column to test against", metavar="character"),
    make_option(c("-n", "--control_name"), type="character", help="the name of controls inside the variable/column to test against", metavar="character"),
    make_option(c("-q", "--padjusted"), type="double", default=0.05, help="p-adjusted threshold value [default= %default]", metavar="QVALUE"),
    make_option(c("-f", "--log_fold_change"), type="double", default=1, help="fold-change threshold value [default= %default]", metavar="FOLD_CHANGE_THRESHOLD"),
    make_option(c("-t", "--separator"), type="character", default="\t", help="input field separator [default= %default]", metavar="SEP"),
    make_option(c("-s", "--subset"), type="character", help="Boolean formula to subset the input data", metavar="character"),
    make_option(c("-m", "--model_organism"), type="character", help="the model organism to use during Gene Ontology step, eg Mus_musculus", metavar="character"),
    make_option(c("-d", "--organism_library"), type="character", help="the R DB/library of the model organism to use during Gene Ontology step, e.g. org.Mm.eg.db for Ensembl version of Mus musculus", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$phenodata))
    stop("Phenodata parameter must be provided (-p). See script usage (--help)")

if (is.null(opt$outdir))
    stop("Output directory parameter must be provided (-o). See script usage (--help)")

if (is.null(opt$var_to_test))
    stop("Variable to test parameter must be provided (-v). See script usage (--help)")

if (is.null(opt$control_name)){
    print("Using lexicographic order of values in column VAR_TO_TEST instead of custom label (e.g., 'control') for differential analysis")
} else {
    print(paste("Using", opt$control_name, "as control value for differential analysis"))
}
    
PHENO_DATA = opt$phenodata
PATH = opt$outdir
if (is.null(PATH)) PATH = "./"
GENE_COUNTS = opt$gene_counts
VAR_TO_TEST = opt$var_to_test
SUBS = opt$subset
CONTROL_NAME = opt$control_name
PADJ = opt$padjusted                 #Set the PADJ value cutoff
FC_THR = opt$log_fold_change            #Set the log2FC value cutoff
ORGANISM=opt$model_organism
ORGANISM_LIBRARY = opt$organism_library    # "org.Mm.eg.db"

#PHENO_DATA = "phenodata.csv"
#GENE_COUNTS =  "gene_count_matrix.Mus_musculus.csv"
#SUBS = "Sex=='f' & Region=='pfc'"
#VAR_TO_TEST ="Stress.protocol"
#CONTROL_NAME="Control"
#ORGANISM_LIBRARY = "org.Mm.eg.db"
#ORGANISM="Mus_musculus"
#FC_THR=0.38
#PADJ=0.05
#PATH="output_pfc_F"

# Saving parameters to file
print("Saving parameters to file")
write.table(t(as.data.frame(opt)), file=paste(PATH, "parameters.txt", sep=""))


###################################
#######  DATA IMPORT ##############
###################################
print("#######  DATA IMPORT ##############")
print(paste("Reading gene counts from ", GENE_COUNTS))
countData <- as.matrix(read.csv(GENE_COUNTS, row.names="gene_id"))
countData <- countData[complete.cases(countData), ]
print("Input matrix dimension is: dim(countData)")

print(paste("Reading phenodata form ", PHENO_DATA))
colData <- read.csv(PHENO_DATA, sep="\t", row.names=1)
print("Metadata dimension is: dim(colData)")

#qualche check preliminare. Devo avere due TRUE per poter continuare:
print("Do phenodata and gene count matrix refer to the same set of samples?")
print("Check n.1:")
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]

print("Check n.2:")
all(rownames(colData) == colnames(countData))


library("DESeq2")
library("limma")
library("genefilter")

###################################
#######  DATA SUBSET AND PRE-FILTERING ##############
###################################

# ## costruisco a monte la stringa per il subset che abbia un aspetto come questa:
# ## SUBS <- quote(condition=="drought" | condition=="watered")
# if (! is.null(SUBS) )
# {
# #Posso così subsettare phenodata e matrice di conta:
# colData_filt <- subset(colData, eval(parse(text = SUBS)))
# countData_subs <- countData[,colnames(countData) %in% rownames(colData_filt)]
# } else {
# colData_filt = colData
# countData_subs = countData
# }
# # Ensure that both conditions exist after subsetting
# print("Checking if both groups are still present after filtering:")
# print(table(colData_filt[["condition"]]))

# if (length(unique(colData_filt[[VAR_TO_TEST]])) < 2) {
#     stop("Error: Subsetting removed one of the groups. Adjust filtering to keep both conditions.")
# }
# Ensure only the selected condition and control are included
# Debug: Print the condition name received from Snakemake
print(paste("Filtering for condition:", opt$var_to_test, "and control"))
print("Available conditions before filtering:")
print(table(colData[["condition"]]))

# Ensure only the selected condition and control are included
colData_filt <- colData[colData$condition %in% c(opt$var_to_test, CONTROL_NAME), , drop=FALSE]
countData_subs <- countData[, colnames(countData) %in% rownames(colData_filt), drop=FALSE]

# Debug: Check if filtering removed all samples
if (nrow(colData_filt) == 0) {
    stop("Error: colData_filt is empty! The condition name might be incorrect or missing.")
}

print("Checking if both groups are still present after filtering:")
print(table(colData_filt[["condition"]]))
# Ensure condition column is a factor with only two levels: {condition} and "control"
colData_filt$condition <- relevel(factor(colData_filt$condition), ref = CONTROL_NAME)

# Debug: Print factor levels
print("Final condition levels:")
print(levels(colData_filt$condition))

# Check again if both conditions are present
if (length(unique(colData_filt$condition)) < 2) {
    stop("Error: The condition or control is missing after filtering. DESeq2 requires at least two groups.")
}

# Ensure that both the selected condition and control are still present
if (length(unique(colData_filt[["condition"]])) < 2) {
    stop(paste("Error: The condition", opt$var_to_test, "or control is missing after filtering."))
}


print(paste("The subsetted dataset is :", dim(countData_subs), sep=""))

#filtering:
library(genefilter)
print("####### DATASET PRE-FILTERING ###########")
print("Removing genes whose mean expression is less than 5 FPKM in at least 80% of the samples")
print(paste("Genes before filtering:", dim(countData_subs)[1]))

#Define a function to remove genes whose mean expression is less than 5 in at least 80% of the samples:
fun <- kOverA(round(dim(countData_subs)[2]*50/100),1)

#Apply it to my matrix:
filter1 <- apply( countData_subs, 1, fun)
countData_filt <- countData_subs[filter1,]
print(paste("Genes after filtering:", dim(countData_filt)[1]))


print(paste("The subsetted and filtered dataset dimension is :", dim(countData_filt), sep=""))

##Questo relevel non è sempre fondamentale. Si può saltare se non c'è un vero e proprio campione di controllo (ad. es Pioppo, maiale). E' importante perché il segno del FC dipende da chi si trova al denominatore. Va invece utilizzato ad esempio nel caso dello stress laddove il valore di riferimento DEVE essere "control".
# if ( ! is.null(CONTROL_NAME) )
# colData_filt[["condition"]] <- relevel(factor(colData_filt[["condition"]]), ref = CONTROL_NAME)


# Check if both conditions are present
if (length(unique(colData_filt$condition)) < 2) {
    stop("Error: The condition or control is missing after filtering. DESeq2 requires at least two groups.")
}

###################################
#######  DESeq2 ANALYSIS ###########
###################################

library(DESeq2)
library(sva)
print("loading ggplot2 package")
library(ggplot2)

# Create a DESeqDataSet from count matrix and labels
# (bisogna ricostruire nuovamente l'oggetto DESeq ogni volta che si cambia il design!)
print("#######  DESeq2 ANALYSIS ###########")
# Add a default batch column if missing
if (!"batch" %in% colnames(colData_filt)) {
    print("No batch information found. Assigning a single batch to all samples.")
    colData_filt$batch <- 1  # All samples assigned to a single batch
}

# Ensure the variable has at least two unique values
if (length(unique(colData_filt$condition)) < 2) {
    stop("Error: The variable to test (`condition`) must have at least two unique values (e.g., Control vs. Treatment).")
}

# Debug: Check sample counts per condition before DESeq2
print("Final sample count per condition before DESeq2:")
print(table(colData_filt$condition))

# Debug: Check matrix dimensions before DESeq2
print(paste("Final count matrix dimensions:", dim(countData_filt)[1], "genes x", dim(countData_filt)[2], "samples"))

# Stop if there's only one condition remaining (to prevent DESeq2 errors)
if (length(unique(colData_filt$condition)) < 2) {
    stop("Error: Only one condition remains in the dataset after pre-filtering. DESeq2 requires at least two.")
}

print("Debug: Checking if colData_filt and countData_filt samples match...")

# Check sample names
print("Samples in colData_filt:")
print(rownames(colData_filt))
print("Samples in countData_filt:")
print(colnames(countData_filt))

# Ensure matching sample names
valid_samples <- intersect(rownames(colData_filt), colnames(countData_filt))

print("Matched sample names:")
print(valid_samples)

# Subset colData_filt using the matched samples
colData_filt <- colData_filt[valid_samples, , drop=FALSE]
countData_filt <- countData_filt[, valid_samples, drop=FALSE]

filtered_conditions <- colData_filt$condition

print("Final sample count per condition after filtering:")
print(table(filtered_conditions))

if (length(unique(filtered_conditions)) < 2) {
    stop("Error: Only one condition remains after filtering. DESeq2 requires at least two.")
}


if (length(unique(filtered_conditions)) < 2) {
    stop("Error: Only one condition remains after filtering. DESeq2 requires at least two.")
}

colData_filt$condition <- factor(colData_filt$condition)  # Ensure 'condition' is a factor
print("Debug: Final colData_filt before DESeq2:")
print(colData_filt)  # Print the full metadata table

print("Debug: Checking colData_filt before DESeq2:")
print(colData_filt)  # Print colData_filt

# Ensure 'condition' is a factor
colData_filt$condition <- factor(colData_filt$condition)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = countData_filt, colData = colData_filt, design = ~ condition)

# Run the default analysis for DESeq2 and generate results table
dds <- DESeq(dds)
#Lowe dice che è meglio usare DESeq e non ReplaceOutliers poiché la prima funzione invoca già la seconda. Inoltre sconsiglia caldamente di scendere sotto il 7 per minReplicatesForReplace.

#diagnostic PCA plot
png(paste(PATH,"PCA_diagnostica.png",sep=""))
plotMDS(countData_filt )
dev.off()

# Variance-stabilizing transformation
print("Applying VST...")
vsd <- vst(dds, blind=FALSE)

# Ensure 'condition' exists in colData(dds) before PCA
print("Debug: Full colData(dds) before PCA (ensure 'condition' exists)...")
print(as.data.frame(colData(dds)))  # Print as a data frame for easy reading

print("Columns in colData(dds):")
print(colnames(colData(dds)))

if (!"condition" %in% colnames(colData(dds))) {
    print("Fixing: Adding 'condition' to colData(dds)...")
    colData(dds)$condition <- colData_filt$condition  # Explicitly copy from colData_filt
}

print("Columns in colData(dds) after fix:")
print(colnames(colData(dds)))

if (!"condition" %in% colnames(colData(dds))) {
    stop("Error: The 'condition' column is still missing from colData(dds) after fixing.")
}

# Step 1: Clustering Exploration
print("Performing PCA for exploratory clustering analysis...")
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw()

ggsave(file=paste(PATH, "Exploratory_PCA_Plot_", gsub(" ","_",opt$var_to_test), "_vs_control.png", sep=""), width=10, height=8)

# Optional: Hierarchical Clustering
sample_dist <- dist(t(assay(vsd)))
hc <- hclust(sample_dist)

png(file=paste(PATH, "Sample_Clustering_Dendrogram.png", sep=""))
plot(hc, main="Sample Clustering Dendrogram", xlab="Samples", sub="")
dev.off()

# Step 2: Detect Hidden Batches with SVA (If No Explicit Batch Variable)
print("Detecting hidden batch effects using SVA...")
mod <- model.matrix(~ condition, colData_filt)
mod0 <- model.matrix(~ 1, colData_filt) #modello nullo senza la variabile sperimentale.

# Run SVA
#svobj <- sva(assay(vsd), mod, mod0) #Confronta la varianza spiegata dal modello sperimentale con la varianza residua.
##SVA dovrebbe essere applicato direttamente ai dati raw count, non ai dati trasformati con VST. Sarebbe più corretto:
svobj <- sva(counts(dds, normalized=TRUE), mod, mod0)

# Add surrogate variables to colData
colData_filt$sv1 <- svobj$sv[, 1]  # Add first surrogate variable
if (ncol(svobj$sv) > 1) {
    colData_filt$sv2 <- svobj$sv[, 2]  # Add second surrogate variable (if available)
}

# Print detected surrogate variables
print("Detected surrogate variables (hidden batch effects):")
print(svobj$sv)

# Step 3: Update DESeq2 Object for SVA
print("Updating DESeq2 design to account for surrogate variables...")
dds <- DESeqDataSetFromMatrix(
    countData = countData_filt,
    colData = colData_filt,
    design = ~ sv1 + condition)
##modello aggiornato che assorbe il batch effect

# Run DESeq2 pipeline with updated design
dds <- DESeq(dds)

# Step 4: Differential Expression Analysis
print("Running differential expression analysis...")
res <- results(dds)

# Step 5: Correct for Batch Effects Using ComBat (If Explicit Batch Variable Exists)
if (length(unique(colData_filt$batch)) > 1) {
    print("Correcting for batch effects using ComBat...")
    batch <- colData_filt$batch
    mod <- model.matrix(~ condition, colData_filt)

    # Apply ComBat
    assay_corrected <- ComBat(dat=assay(vsd), batch=batch, mod=mod)
    assay(vsd) <- assay_corrected  # Replace assay with corrected data
    ###ATTENZIONE: Modificare assay(vsd) non ha alcun effetto sull’analisi DESeq2, 
    ##perché vsd viene usato solo per la PCA e non per calcolare i DEGs.
    # Ensure 'condition' exists in colData(dds) before PCA
    print("Debug: Full colData(dds) before batch-corrected PCA (ensure 'condition' exists)...")
    print(as.data.frame(colData(dds)))  # Print as a data frame for easy reading

    print("Columns in colData(dds):")
    print(colnames(colData(dds)))

    if (!"condition" %in% colnames(colData(dds))) {
        print("Fixing: Adding 'condition' to colData(dds)...")
        colData(dds)$condition <- colData_filt$condition  # Explicitly copy from colData_filt
    }

    print("Columns in colData(dds) after fix:")
    print(colnames(colData(dds)))

    if (!"condition" %in% colnames(colData(dds))) {
        stop("Error: The 'condition' column is still missing from colData(dds) after fixing.")
    }

    print("Performing PCA after batch effect correction...")
    pcaData_corrected <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)  # Fix VAR_TO_TEST

    # Save PCA plot after batch correction
    print("Performing PCA after batch effect correction...")
    pcaData_corrected <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
    percentVar_corrected <- round(100 * attr(pcaData_corrected, "percentVar"))

    ggplot(pcaData_corrected, aes(PC1, PC2, color=!!sym(VAR_TO_TEST))) +
      geom_point(size=3) +
      xlab(paste0("PC1: ", percentVar_corrected[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar_corrected[2], "% variance")) +
      coord_fixed() +
      theme_bw()

    ggsave(file=paste(PATH, "Corrected_PCA_Plot_", gsub(" ","_",opt$var_to_test), "_vs_control.png", sep=""), width=10, height=8)
} else {
    print("Only one batch detected. Skipping ComBat correction.")
}

######PARTE NUOVA:
countData_norm <- counts(dds, normalized=TRUE) 
batch <- colData_filt$batch 
mod <- model.matrix(~ condition, colData_filt)
if (length(unique(batch)) < 2) {
    print("Only one batch detected. Skipping ComBat correction.")
    countData_corrected <- round(countData_norm)  # Convert to integers
} else {
    countData_corrected <- ComBat(dat=countData_norm, batch=batch, mod=mod)
    countData_corrected <- round(countData_corrected)  # Ensure integer values
}


 
dds_corrected <- DESeqDataSetFromMatrix(countData = countData_corrected,  colData = colData_filt,  design = ~ condition) 
# Rianalisi con DESeq2 tenendo conto della correzione di COMBAT:
dds_corrected <- DESeq(dds_corrected) 
res <- results(dds_corrected)
##FINE NUOVA PARTE

#add an extra column to mark outliers:
res$outlier = res$baseMean > 0 & is.na(res$pvalue)


# Computing and Exporting only the results which pass an adjusted p value threshold and a particular FC value:
# First of all, we need to subset the matrix:
print("Subsetting dataset")
resSig <- subset(res, padj < PADJ & abs(log2FoldChange)>=FC_THR)
#sort and show by FC value:
resSig <- resSig[order(abs(resSig$log2FoldChange),decreasing=T), ]

##########SAVE WHOLE DATA RESULTS###################
###################################################

# Sort and filt  by adjusted p-value:
resOrdered <- res[order(res$padj), ]

#####################################################
### Compose the output file's name in a parametric manner. ###
#####################################################
# This string gives a short summary of the tested categories and I parse and use it in the output's name:
tested_cat <- resOrdered@elementMetadata@listData$description[5]
tested_cat <- gsub("Wald test p-value: ", "" , tested_cat)

# Exporting NOT filtered results:
write.csv(as.data.frame(resOrdered), file=paste(PATH, paste("gene_expr-", gsub(" ","_",tested_cat), ".csv", sep=""), sep=""))

###################################
######### SIGNIFICANT DATA SAVE  ##############
###################################

CONDITION_STRING = paste("qvalue", PADJ, "log2FC", FC_THR, gsub(" ","_",tested_cat), sep=".")

# Exporting filtered results:
write.csv(as.data.frame(resSig), file=paste(PATH, "DEGs-", CONDITION_STRING, ".csv", sep=""))

print("Written table of DEGs")

##access the value of a specific gene:
##assay(rld)["Plekhg2",]

###################################
######### GRAPHICAL PART #########
###################################

print("######### GRAPHICAL PART #########")


################################################################################Graphic part:
##Plot2: Istogramma di tutti i p.value:
png(file=paste(PATH, "pvalues-histogram" , ".png", sep=""), res=100)
hist(res$pvalue,main="p.values distribution", col="green4", breaks=20, xlab="p-value")
abline(v=0.05,col="red",lwd=3, lty=2)
dev.off()

##Plot2bis: Istogramma di tutti i q.value:
png(file=paste(PATH, "qvalues-histogram" , ".png", sep=""), res=100 )
hist(res$padj,main="q-values distribution", ylab="q-values", col="green4", breaks=20, ylim=c(0, 2500), xlab="q-value")
abline(v=0.05,col="red",lwd=3, lty=2)
dev.off()

# Vulcano plot on DEGs samples:
toptable <- as.data.frame(res)

EnhancedVolcanoDESeq2 <- function(toptable, AdjustedCutoff=PADJ, LabellingCutoff=PADJ, FCCutoff=FC_THR, main="VolcanoPlot")
{
  toptable$Significance <- "NS"
  toptable$Significance[(abs(toptable$log2FoldChange) > FCCutoff)] <- "FC"
  toptable$Significance[(toptable$padj<AdjustedCutoff)] <- "FDR"
  toptable$Significance[(toptable$padj<AdjustedCutoff) & (abs(toptable$log2FoldChange)>FCCutoff)] <- "FC_FDR"
  table(toptable$Significance)
  toptable$Significance <- factor(toptable$Significance, levels=c("NS", "FC", "FDR", "FC_FDR"))

  plot <- ggplot(toptable, aes(x=log2FoldChange, y=-log10(padj))) +
    #scale_y_continuous(breaks = seq(0, 150, by = 25)) +
    #axis.break(axis=1,breakpos=100,pos=NA,bgcol="white",breakcol="black", style="slash",brw=0.02)
    
    #Add points:
    #   Colour based on factors set a few lines up
    #   'alpha' provides gradual shading of colour
    #   Set size of points
    geom_point(aes(color=factor(Significance)), alpha=1/2, size=0.8) +

    #Choose which colours to use; otherwise, ggplot2 choose automatically (order depends on how factors are ordered in toptable$Significance)
    scale_color_manual(values=c(NS="grey30", FC="forestgreen", FDR="royalblue", FC_FDR="red2"), labels=c(NS="NS", FC=paste("LogFC>|", FCCutoff, "|", sep=""), FDR=paste("FDR Q<", AdjustedCutoff, sep=""), FC_FDR=paste("FDR Q<", AdjustedCutoff, " & LogFC>|", FCCutoff, "|", sep=""))) +

    #Set the size of the plotting window
    theme_bw(base_size=24) +

    #Modify various aspects of the plot text and legend
    theme(legend.background=element_rect(),
          plot.title=element_text(angle=0, size=12, face="bold", vjust=1),

          panel.grid.major=element_blank(), #Remove gridlines
          panel.grid.minor=element_blank(), #Remove gridlines

          axis.text.x=element_text(angle=0, size=12, vjust=1),
          axis.text.y=element_text(angle=0, size=12, vjust=1),
          axis.title=element_text(size=12),

          #Legend
          legend.position="top",            #Moves the legend to the top of the plot
          legend.key=element_blank(),       #removes the border
          legend.key.size=unit(0.5, "cm"),  #Sets overall area/size of the legend
          legend.text=element_text(size=8), #Text size
          title=element_text(size=8),       #Title text size
          legend.title=element_blank()) +       #Remove the title

    #Change the size of the icons/symbols in the legend
    guides(colour = guide_legend(override.aes=list(size=2.5))) +

    #Set x- and y-axes labels
    xlab(bquote(~Log[2]~ "fold change")) +
    ylab(bquote(~-Log[10]~adjusted~italic(P))) +

    #Set the axis limits
    #xlim(-6.5, 6.5) +
    #ylim(0, 100) +

    #Set title
    ggtitle(main) +

    #Tidy the text labels for a subset of genes
    geom_text(data=subset(toptable, padj<LabellingCutoff & abs(log2FoldChange)>FCCutoff),
              aes(label=rownames(subset(toptable, padj<LabellingCutoff & abs(log2FoldChange)>FCCutoff))),
              size=2.25,
              #segment.color="black", #This and the next parameter spread out the labels and join them to their points by a line
              #segment.size=0.01,
              check_overlap=TRUE,
              vjust=1.0) +

    #Add a vertical line for fold change cut-offs
    geom_vline(xintercept=c(-FCCutoff, FCCutoff), linetype="longdash", colour="black", size=0.4) +

    #Add a horizontal line for P-value cut-off
    geom_hline(yintercept=-log10(AdjustedCutoff), linetype="longdash", colour="black", size=0.4)

  return(plot)
}

# Create and save the volcano plot (if the matrix contains rows):
if(nrow(toptable)!=0)
{
print("Table of DEGs contains results: I'm plotting a volcano plot")
volc <- EnhancedVolcanoDESeq2(toptable, AdjustedCutoff=PADJ, LabellingCutoff=PADJ, FCCutoff=FC_THR, main=paste("Volcano Plot of differentially expressed genes", gsub(" ","_",tested_cat), sep="\n"))
ggsave(file=paste(PATH, "DEGs_volcanoplot_", CONDITION_STRING, ".png", sep=""), volc)
} else {
print("Table of DEGs does not contain any results: No volcano plot produced")
}

# Enhanced Volcano Plot
library(EnhancedVolcano)

EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = PADJ,
    FCcutoff = FC_THR,
    pointSize = 3.0,
    labSize = 6.0)

ggsave(file=paste(PATH, "Enhanced_Volcano_Plot_", gsub(" ","_",opt$var_to_test), "_vs_control.png", sep=""), width=12, height=10)

# Improved PCA Plot
vsd <- vst(dds, blind=FALSE)

print("Debug: Full colData(dds) before PCA (ensure 'condition' exists)...")
print(as.data.frame(colData(dds)))  # Print as a data frame for easy reading

# Ensure 'condition' is in colData(dds)
if (!"condition" %in% colnames(colData(dds))) {
    print("Fixing: Adding 'condition' to colData(dds)...")
    colData(dds)$condition <- colData_filt$condition  # Explicitly copy from colData_filt
}

print("Columns in colData(dds):")
print(colnames(colData(dds)))

if (!"condition" %in% colnames(colData(dds))) {
    stop("Error: The 'condition' column is still missing from colData(dds) after fixing.")
}




pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  theme_bw()

ggsave(file=paste(PATH, "Improved_PCA_Plot_", gsub(" ","_",opt$var_to_test), "_vs_control.png", sep=""), width=10, height=8)

# Heatmap Generation
library(pheatmap)

# Select top 50 genes by adjusted p-value
top_genes <- head(order(res$padj), 50)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)

pheatmap(mat, 
        annotation_col = colData_filt[, "condition", drop=FALSE],
         show_rownames = FALSE,
         filename = paste(PATH, "Top50_DEGs_Heatmap_", gsub(" ","_",opt$var_to_test), "_vs_control.png", sep=""))


###Sample distances

#png(file=paste(PATH, "samples_dist", ".png", sep=""))
#dist(t(assay(vsd)))
#dev.off()

###################################
#### KEGG and GO Enrichment #######
###################################

library(clusterProfiler)
library(org.Hs.eg.db) # Replace with appropriate OrgDb for your organism
library(ggplot2)

# Convert gene IDs to Entrez IDs
# Split the identifiers
gene_ids <- rownames(resSig)
ensembl_ids <- sapply(strsplit(gene_ids, "\\|"), `[`, 1)

# Map to Entrez IDs
entrez_ids <- mapIds(org.Hs.eg.db, 
                     keys = ensembl_ids,
                     keytype = "ENSEMBL",
                     column = "ENTREZID")

# Remove any NA values
entrez_ids <- entrez_ids[!is.na(entrez_ids)]

# KEGG Pathway Analysis
kegg_result <- enrichKEGG(gene = entrez_ids,
                          organism = 'hsa',
                          pvalueCutoff = 0.05)

# Save KEGG results
write.csv(as.data.frame(kegg_result), file=paste(PATH, "KEGG_enrichment_", gsub(" ","_",opt$var_to_test), "_vs_control.csv", sep=""))

# GO Enrichment Analysis
go_result <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Hs.eg.db,
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

# Save GO results
write.csv(as.data.frame(go_result), file=paste(PATH, "GO_enrichment_", gsub(" ","_",opt$var_to_test), "_vs_control.csv", sep=""))

# Visualizations
# png(file=paste(PATH, "KEGG_dotplot_", gsub(" ","_",opt$var_to_test), "_vs_control.png", sep=""))
# dotplot(kegg_result, showCategory=20)
# dev.off()

# png(file=paste(PATH, "GO_dotplot_", gsub(" ","_",opt$var_to_test), "_vs_control.png", sep=""))
# dotplot(go_result, showCategory=20)
# dev.off()
