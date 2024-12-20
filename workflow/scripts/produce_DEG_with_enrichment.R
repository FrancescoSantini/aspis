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

## costruisco a monte la stringa per il subset che abbia un aspetto come questa:
## SUBS <- quote(condition=="drought" | condition=="watered")
if (! is.null(SUBS) )
{
#Posso così subsettare phenodata e matrice di conta:
colData_filt <- subset(colData, eval(parse(text = SUBS)))
countData_subs <- countData[,colnames(countData) %in% rownames(colData_filt)]
} else {
colData_filt = colData
countData_subs = countData
}

print(paste("The subsetted dataset is :", dim(countData_subs), sep=""))

#filtering:
library(genefilter)
print("####### DATASET PRE-FILTERING ###########")
print("Removing genes whose mean expression is less than 5 FPKM in at least 80% of the samples")
#Define a function to remove genes whose mean expression is less than 5 in at least 80% of the samples:
fun <- kOverA(round(dim(countData_subs)[2]*80/100),5)
#Apply it to my matrix:
filter1 <-  apply( countData_subs, 1, fun)
countData_filt <- countData_subs[filter1,]
print(paste("The subsetted and filtered dataset dimension is :", dim(countData_filt), sep=""))

##Questo relevel non è sempre fondamentale. Si può saltare se non c'è un vero e proprio campione di controllo (ad. es Pioppo, maiale). E' importante perché il segno del FC dipende da chi si trova al denominatore. Va invece utilizzato ad esempio nel caso dello stress laddove il valore di riferimento DEVE essere "control".
if ( ! is.null(CONTROL_NAME) )
colData_filt[[VAR_TO_TEST]] <- relevel(factor(colData_filt[[VAR_TO_TEST]]), ref = CONTROL_NAME)

###################################
#######  DESeq2 ANALYSIS ###########
###################################

library(DESeq2)

# Create a DESeqDataSet from count matrix and labels
# (bisogna ricostruire nuovamente l'oggetto DESeq ogni volta che si cambia il design!)
print("#######  DESeq2 ANALYSIS ###########")
dds <- DESeqDataSetFromMatrix(countData = countData_filt, colData = colData_filt, design = as.formula(paste("~", VAR_TO_TEST)))

# Run the default analysis for DESeq2 and generate results table
dds <- DESeq(dds)
#Lowe dice che è meglio usare DESeq e non ReplaceOutliers poiché la prima funzione invoca già la seconda. Inoltre sconsiglia caldamente di scendere sotto il 7 per minReplicatesForReplace.

#diagnostic plot
png(paste(PATH,"PCA_diagnostica.png",sep=""))
plotMDS(countData_filt )
dev.off()

res <- results(dds)
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
print("loading ggplot2 package")
library(ggplot2)

################################################################################Graphic part:
##Plot2: Istogramma di tutti i p.value:
png(file=paste(PATH, "pvalues-histogram.png" , ".png", sep=""), res=100)
hist(res$pvalue,main="p.values distribution", col="green4", breaks=20, xlab="p-value")
abline(v=0.05,col="red",lwd=3, lty=2)
dev.off()

##Plot2bis: Istogramma di tutti i q.value:
png(file=paste(PATH, "qvalues-histogram.png" , ".png", sep=""), res=100 )
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
ggsave(file=paste(PATH, "DEGs_volcanoplot-", CONDITION_STRING, ".png", sep=""), volc)
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

ggsave(file=paste(PATH, "Enhanced_Volcano_Plot.png", sep=""), width=12, height=10)

# Improved PCA Plot
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c(VAR_TO_TEST), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=!!sym(VAR_TO_TEST))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  theme_bw()

ggsave(file=paste(PATH, "Improved_PCA_Plot.png", sep=""), width=10, height=8)

# Heatmap Generation
library(pheatmap)

# Select top 50 genes by adjusted p-value
top_genes <- head(order(res$padj), 50)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)

pheatmap(mat, 
         annotation_col = colData_filt[, c(VAR_TO_TEST), drop=FALSE],
         show_rownames = FALSE,
         filename = paste(PATH, "Top50_DEGs_Heatmap.png", sep=""))


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
write.csv(as.data.frame(kegg_result), file=paste(PATH, "KEGG_enrichment.csv", sep=""))

# GO Enrichment Analysis
go_result <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Hs.eg.db,
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

# Save GO results
write.csv(as.data.frame(go_result), file=paste(PATH, "GO_enrichment.csv", sep=""))

# Visualizations
png(file=paste(PATH, "KEGG_dotplot.png", sep=""))
dotplot(kegg_result, showCategory=20)
dev.off()

png(file=paste(PATH, "GO_dotplot.png", sep=""))
dotplot(go_result, showCategory=20)
dev.off()
