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

# Read count matrix and preserve original column names (avoid R prepending 'X' to numeric headers)
print(paste("Reading gene counts from ", GENE_COUNTS))
countData <- as.matrix(read.csv(GENE_COUNTS, row.names = "gene_id", check.names = FALSE))
countData <- countData[complete.cases(countData), ]
print(paste("Input matrix dimension is:", paste(dim(countData), collapse = "x")))

# Read phenodata
print(paste("Reading phenodata from ", PHENO_DATA))
colData <- read.csv(PHENO_DATA, sep = "\t", row.names = 1)
print(paste("Metadata dimension is:", paste(dim(colData), collapse = "x")))

# Confirm matching samples
print("Do phenodata and gene count matrix refer to the same set of samples?")

# Check 1: all samples in phenodata are present in count matrix
print("Check n.1:")
check1 <- all(rownames(colData) %in% colnames(countData))
print(check1)

# Check 2: same order
print("Check n.2:")
check2 <- all(rownames(colData) == colnames(countData))
print(check2)

# Debug if checks fail
if (!check1 || !check2) {
    print("🔍 DEBUG: Listing sample name mismatches")

    pheno_samples <- rownames(colData)
    count_samples <- colnames(countData)

    for (i in seq_along(pheno_samples)) {
        if (i > length(count_samples)) {
            cat(sprintf("⚠️ Extra phenodata sample at position %d: '%s'\n", i, pheno_samples[i]))
            next
        }

        ps <- pheno_samples[i]
        cs <- count_samples[i]
        if (ps != cs) {
            cat(sprintf("❌ Mismatch at position %d:\n", i))
            cat(sprintf("  phenodata: '%s' | ASCII: %s\n", ps, paste(as.integer(charToRaw(ps)), collapse = ' ')))
            cat(sprintf("  countData: '%s' | ASCII: %s\n", cs, paste(as.integer(charToRaw(cs)), collapse = ' ')))
        }
    }

    if (length(count_samples) > length(pheno_samples)) {
        for (j in (length(pheno_samples) + 1):length(count_samples)) {
            cat(sprintf("⚠️ Extra count matrix sample at position %d: '%s'\n", j, count_samples[j]))
        }
    }

    if (!check1) {
        stop("❌ Error: Some samples in phenodata are missing from the count matrix.")
    }
    if (!check2) {
        stop("❌ Error: Sample names are not in the same order. Please reorder columns of countData.")
    }
}

# Subset count matrix to match phenodata order
countData <- countData[, rownames(colData)]

# Final confirmation
print("✅ Final sample order check:")
print(all(rownames(colData) == colnames(countData)))



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


print(paste("The subsetted dataset is:", paste(dim(countData_subs), collapse="x")))

#filtering:
library(genefilter)
print("####### DATASET PRE-FILTERING ###########")
print("Removing genes whose mean expression is less than 5 FPKM in at least 80% of the samples")
print(paste("Genes before filtering:", dim(countData_subs)[1]))

#Define a function to remove genes whose mean expression is less than 5 in at least 80% of the samples:
fun <- kOverA(round(dim(countData_subs)[2]*80/100),5)

#Apply it to my matrix:
filter1 <- apply( countData_subs, 1, fun)
countData_filt <- countData_subs[filter1,]
print(paste("Genes after filtering:", dim(countData_filt)[1]))


print(paste("The subsetted and filtered dataset dimension is :", paste(dim(countData_filt), collapse="x")))

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

# Check if SVA detected any surrogate variables
if (ncol(svobj$sv) > 0) {
    # Add first surrogate variable to colData
    colData_filt$sv1 <- svobj$sv[, 1]
    if (ncol(svobj$sv) > 1) {
        colData_filt$sv2 <- svobj$sv[, 2]  # Add second surrogate variable if available
    }
    
    print("Detected surrogate variables (hidden batch effects):")
    print(svobj$sv)
    
    # Update DESeq2 design to include SVA batch correction
    print("Updating DESeq2 design to account for surrogate variables...")
    dds <- DESeqDataSetFromMatrix(
        countData = countData_filt,
        colData = colData_filt,
        design = ~ sv1 + condition
    )
} else {
    print("No significant surrogate variables detected. Proceeding without SVA correction.")
    
    # Keep the original design without SVA correction
    dds <- DESeqDataSetFromMatrix(
        countData = countData_filt,
        colData = colData_filt,
        design = ~ condition
    )
}

# Run DESeq2 pipeline with the appropriate design
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
#### KEGG, GO, and Reactome Enrichment #######
###################################

library(clusterProfiler)
library(org.Hs.eg.db) # Replace with appropriate OrgDb for your organism
library(ggplot2)
library(goseq)
library(biomaRt)
library(ReactomePA)

# Convert gene IDs to Entrez IDs
gene_ids <- rownames(resSig)
# Subset DEGs
resSig <- subset(res, padj < PADJ & abs(log2FoldChange) >= FC_THR)
resSig <- resSig[order(abs(resSig$log2FoldChange), decreasing = TRUE), ]

# ✅ Stop enrichment early if resSig is empty
if (nrow(resSig) == 0) {
    print("⚠️ No significant DEGs found. Skipping enrichment analysis.")

    # Write empty placeholder CSVs
    write.csv(data.frame(), file=paste0(PATH, "KEGG_enrichment_", gsub(" ", "_", opt$var_to_test), "_vs_control.csv"))
    write.csv(data.frame(), file=paste0(PATH, "GO_enrichment_", gsub(" ", "_", opt$var_to_test), "_vs_control.csv"))
    write.csv(data.frame(), file=paste0(PATH, "GOseq_enrichment_", gsub(" ", "_", opt$var_to_test), "_vs_control.csv"))
    write.csv(data.frame(), file=paste0(PATH, "Reactome_enrichment_", gsub(" ", "_", opt$var_to_test), "_vs_control.csv"))

    # Skip the rest of this script
    quit(save = "no", status = 0)
}

ensembl_ids <- sapply(strsplit(gene_ids, "\\|"), `[`, 1)

# Clean and filter
ensembl_ids <- ensembl_ids[!is.na(ensembl_ids) & ensembl_ids != ""]

if (length(ensembl_ids) == 0) {
    stop("❌ Error: No valid Ensembl IDs found in DEGs.")
}

# Map to Entrez IDs (fail safely if nothing found)
entrez_ids <- tryCatch({
    mapIds(org.Hs.eg.db,
           keys = ensembl_ids,
           keytype = "ENSEMBL",
           column = "ENTREZID",
           multiVals = "first")
}, error = function(e) {
    stop("❌ Error during mapIds(): ", e$message)
})

# Remove any NA values
entrez_ids <- entrez_ids[!is.na(entrez_ids)]

# Final check
if (length(entrez_ids) == 0) {
    stop("❌ Error: No Entrez IDs were mapped. Check gene ID format.")
}

print(paste("✅ Mapped Entrez IDs:", length(entrez_ids)))


# KEGG Pathway Analysis
kegg_result <- enrichKEGG(gene = entrez_ids,
                          organism = 'hsa',
                          pvalueCutoff = 0.05)

# Save KEGG results
write.csv(as.data.frame(kegg_result), file=paste(PATH, "KEGG_enrichment_", gsub(" ", "_", opt$var_to_test), "_vs_control.csv", sep=""))

# GO Enrichment Analysis
go_result <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Hs.eg.db,
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

# Save GO results
write.csv(as.data.frame(go_result), file=paste(PATH, "GO_enrichment_", gsub(" ", "_", opt$var_to_test), "_vs_control.csv", sep=""))

# Define a function to get gene lengths dynamically
get_gene_lengths <- function(entrez_ids, organism) {
    ensembl_datasets <- list(
        "Homo sapiens" = "hsapiens_gene_ensembl",
        "Mus musculus" = "mmusculus_gene_ensembl",
        "Rattus norvegicus" = "rnorvegicus_gene_ensembl",
        "Danio rerio" = "drerio_gene_ensembl",
        "Drosophila melanogaster" = "dmelanogaster_gene_ensembl",
        "Caenorhabditis elegans" = "celegans_gene_ensembl"
    )

    if (!(organism %in% names(ensembl_datasets))) {
        stop(paste("Error: Organism", organism, "is not supported!"))
    }
    
    dataset <- ensembl_datasets[[organism]]
    mirrors <- c("www", "useast", "uswest", "asia")

    for (mirror in mirrors) {
        print(paste("Trying Ensembl mirror:", mirror))
        tryCatch({
            mart <- useEnsembl(biomart = "genes", dataset = dataset, mirror = mirror)
            gene_annotations <- getBM(attributes = c("entrezgene_id", "transcript_length"),
                                      filters = "entrezgene_id",
                                      values = entrez_ids,
                                      mart = mart)

            # Ensure only genes with lengths are kept
            gene_annotations <- gene_annotations[!is.na(gene_annotations$transcript_length), ]

            # Convert to named vector with Entrez IDs
            gene_lengths <- setNames(gene_annotations$transcript_length, gene_annotations$entrezgene_id)
            
            if (length(gene_lengths) == 0) {
                stop("Error: No gene lengths retrieved! Check Ensembl query.")
            }
            
            return(gene_lengths)
        }, error = function(e) {
            print(paste("Failed to connect to", mirror, ":", e$message))
        })
    }
    
    stop("All Ensembl mirrors are unavailable. Try again later.")
}

# Set organism
selected_organism <- "Homo sapiens"

# Ensure Entrez IDs are character
entrez_ids <- as.character(entrez_ids)

# Get gene lengths using corrected function
gene_lengths <- get_gene_lengths(entrez_ids, selected_organism)

# Ensure all gene length IDs are character
names(gene_lengths) <- as.character(names(gene_lengths))

# Ensure valid genes exist before filtering
valid_genes <- intersect(entrez_ids, names(gene_lengths))

if (length(valid_genes) == 0) {
    stop("Error: No valid genes found after intersection! Check ID conversion.")
}

# Assign DE_genes
DE_genes <- as.integer(entrez_ids %in% valid_genes)
names(DE_genes) <- valid_genes

# Remove duplicate transcript lengths
gene_lengths <- gene_lengths[!duplicated(names(gene_lengths))]

# Assign default length to missing genes
default_length <- median(gene_lengths, na.rm = TRUE)
missing_genes <- setdiff(as.character(entrez_ids), as.character(names(gene_lengths)))
for (gene in missing_genes) {
    if (!(gene %in% names(gene_lengths))) {
        gene_lengths[gene] <- default_length
    }
}

# Final matching step
valid_genes <- intersect(names(DE_genes), names(gene_lengths))
DE_genes <- DE_genes[valid_genes]
gene_lengths <- gene_lengths[valid_genes]

# Final validation
if (length(DE_genes) != length(gene_lengths)) {
    stop("❌ Error: DE_genes and gene_lengths must have the same length before GOseq!")
}

print("✅ Successfully matched DE_genes and gene_lengths. Proceeding with GOseq...")

# GOseq Analysis
print("🚀 Running GOseq analysis...")

# Verify genome version before proceeding
available_genomes <- supportedGenomes()
if (!("hg38" %in% available_genomes)) {
    print("⚠️ Warning: hg38 genome not found in GOseq. Trying hg19 instead...")
    genome_version <- "hg19"
} else {
    genome_version <- "hg38"
}
print(paste("✅ Using genome version:", genome_version))

# Run GOseq with enhanced error handling
suppressWarnings({
    print("🚀 Running nullp() for GOseq...")
    pwf <- tryCatch(
        nullp(DE_genes, genome_version, "ensGene", bias.data = gene_lengths),
        error = function(e) {
            print("❌ Error in nullp():")
            print(as.character(e))
            return(NULL)
        }
    )
    
    if (is.null(pwf)) {
        stop("❌ Error: PWF is NULL, possibly due to invalid gene IDs or missing annotations.")
    }
    print("✅ Successfully generated pwf. Proceeding with GO annotation mapping...")
    
    # Ensure GO annotations are correctly mapped before proceeding
    go_map <- tryCatch(
        getgo(names(DE_genes), genome_version, "ensGene"),
        error = function(e) {
            print("❌ Error in getgo():")
            print(as.character(e))
            return(NULL)
        }
    )
    
    if (is.null(go_map)) {
        stop("❌ Error: GO annotation mapping returned NULL. Check genome version and database availability.")
    }
    
    print(paste("✅ Initial GO mappings retrieved:", length(go_map)))
    print("🔍 Sample GO mappings:")
    print(utils::head(go_map, 10))  # Print first 10 mappings for verification
    
    # Validate and clean the GO map before proceeding
    go_map <- go_map[!sapply(go_map, is.null)]
    go_map <- go_map[lengths(go_map) > 0]  # Ensure no empty lists
    
    if (length(go_map) == 0) {
        warning("⚠️ Warning: Filtered GO annotation mapping is empty. Proceeding with an empty GOseq analysis.")
    } else {
        print(paste("✅ Filtered GO mappings available:", length(go_map)))
    }
    
    # Run GOseq only if we have valid mappings
    if (length(go_map) > 0) {
        print("✅ GO annotation mapping cleaned and validated. Running GOseq...")
        goseq_results <- tryCatch(
            goseq(pwf, genome_version, "ensGene", use_genes_without_cat = TRUE),
            error = function(e) {
                print("❌ Error in goseq():")
                print(as.character(e))
                return(NULL)
            }
        )
        
        if (is.null(goseq_results) || nrow(goseq_results) == 0) {
            warning("⚠️ GOseq returned no results. Placeholder file will be created.")
            goseq_results <- data.frame(Term=character(), Ontology=character(), Pvalue=numeric())
        }
    } else {
        goseq_results <- data.frame(Term=character(), Ontology=character(), Pvalue=numeric())
        print("⚠️ GOseq was skipped due to missing GO annotations. Placeholder file will be created.")
    }
})

# Save GOseq results
output_file <- paste(PATH, "GOseq_enrichment_", gsub(" ", "_", opt$var_to_test), "_vs_control.csv", sep="")
write.csv(goseq_results, file=output_file, row.names=FALSE)


print("✅ GOseq analysis completed with placeholder file handling.")

# Run Reactome Pathway Analysis
reactome_result <- enrichPathway(gene = entrez_ids, organism = "human", pvalueCutoff = 0.05)

# Save Reactome results
write.csv(as.data.frame(reactome_result), file=paste(PATH, "Reactome_enrichment_", gsub(" ", "_", opt$var_to_test), "_vs_control.csv", sep=""))

# Generate dotplots for KEGG, GO, GOseq, and Reactome (ensuring placeholders are created if needed)
library(ggplot2)
library(clusterProfiler)

safe_dotplot <- function(enrichment_result, filename, title) {
    png(filename)
    
    if (!is.null(enrichment_result) && nrow(as.data.frame(enrichment_result)) > 0) {
        print(paste("Generating dotplot for:", title))
        print(paste("Saving file to:", filename))
        tryCatch({
            p <- dotplot(enrichment_result, showCategory=20) + ggtitle(title)
            print(p)
            dev.off()
        }, error = function(e) {
            print(paste("Error in generating dotplot for:", title, "-", e$message))
            dev.off()
        })
    } else {
        print(paste("No data for:", title, "Generating placeholder."))
        tryCatch({
            p <- ggplot() + 
                annotate("text", x=1, y=1, label="No Data Available", size=6, color="red") +
                theme_void() + 
                ggtitle(title)
            print(p)
            dev.off()
        }, error = function(e) {
            print(paste("Error in generating placeholder for:", title, "-", e$message))
            dev.off()
        })
    }
}

# Ensure all dotplot files exist even if empty
dotplot_files <- list(
    paste(PATH, "KEGG_dotplot_", gsub(" ", "_", opt$var_to_test), "_vs_control.png", sep=""),
    paste(PATH, "GO_dotplot_", gsub(" ", "_", opt$var_to_test), "_vs_control.png", sep=""),
    paste(PATH, "GOseq_dotplot_", gsub(" ", "_", opt$var_to_test), "_vs_control.png", sep=""),
    paste(PATH, "Reactome_dotplot_", gsub(" ", "_", opt$var_to_test), "_vs_control.png", sep="")
)

# Generate dotplots for KEGG, GO, GOseq, and Reactome
if (exists("kegg_result") && !is.null(kegg_result)) safe_dotplot(kegg_result, dotplot_files[[1]], "KEGG Pathway Enrichment")
if (exists("go_result") && !is.null(go_result)) safe_dotplot(go_result, dotplot_files[[2]], "GO Enrichment")
if (exists("goseq_results") && !is.null(goseq_results)) safe_dotplot(goseq_results, dotplot_files[[3]], "GOseq Enrichment")
if (exists("reactome_result") && !is.null(reactome_result)) safe_dotplot(reactome_result, dotplot_files[[4]], "Reactome Pathway Enrichment")

# Ensure placeholder files exist in case of missing plots
for (file in dotplot_files) {
    if (!file.exists(file)) {
        tryCatch({
            png(file)
            print(paste("Generating placeholder for:", file))
            p <- ggplot() + 
                annotate("text", x=1, y=1, label="No Data Available", size=6, color="red") +
                theme_void() + 
                ggtitle("Placeholder Plot")
            print(p)
            dev.off()
        }, error = function(e) {
            print(paste("Error in creating placeholder for:", file, "-", e$message))
            dev.off()
        })
    }
}

print("✅ Dotplots generated and placeholders ensured.")