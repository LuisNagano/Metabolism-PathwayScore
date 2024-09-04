setwd("C:/Users/bruep/OneDrive/Área de Trabalho/Kleiton/Pathways/BROAD_~1/LUSC")

library(Rsubread)
library(limma)
library(edgeR)
options(digits=2)

####################################################################
####################### DEG FILE PREPARATION #######################
####################################################################
# This section prepares the DEG (Differentially Expressed Genes) files
# for multiple datasets, as required for the analysis.

# Loading raw counts dataset for LUSC-TCGA
LUSCgenomicMatrix <- read.delim("C:/Users/bruep/OneDrive/Área de Trabalho/Kleiton/Pathways/BROAD_~1/LUSC/GDACBR~1.0/GDACBR~1.0/LUSC.uncv2.mRNAseq_raw_counts.txt", row.names=1)

# Function to clean gene names by removing the pipe ('|') and any suffixes.
clean_gene_names <- function(gene_names) {
  sapply(gene_names, function(x) strsplit(x, "\\|")[[1]][1])
}

# Apply cleaning function to gene names.
cleaned_gene_names <- clean_gene_names(rownames(LUSCgenomicMatrix))

# Ensure unique gene names by appending numeric suffixes to duplicates.
make_unique <- function(names) {
  unique_names <- names
  dup_count <- table(names)
  for (name in names(dup_count)[dup_count > 1]) {
    indices <- which(names == name)
    unique_names[indices] <- paste0(name, "_", seq_along(indices))
  }
  return(unique_names)
}

unique_cleaned_gene_names <- make_unique(cleaned_gene_names)
rownames(LUSCgenomicMatrix) <- unique_cleaned_gene_names

# Verify the changes to the row names.
head(rownames(LUSCgenomicMatrix))

####################################################################
####################### GROUP SELECTION ############################
####################################################################
# Splitting the dataset into two groups: Normal and Tumor samples.

NormalLUSC <- LUSCgenomicMatrix[,grep(".11", colnames(LUSCgenomicMatrix))]
TumorLUSC <- LUSCgenomicMatrix[,grep(".01", colnames(LUSCgenomicMatrix))]

# Design matrix creation for the analysis (Normal vs Tumor)
NormalDesign <- as.data.frame(colnames(NormalLUSC))
NormalDesign$Pheno <- "Normal"
NormalDesign$Color <- "Purple" # Visualization purposes (Optional)
colnames(NormalDesign) <- c("SampleID", "Pheno", "Color")

TumorDesign <- as.data.frame(colnames(TumorLUSC))
TumorDesign$Pheno <- "Tumor"
TumorDesign$Color <- "Green" # Visualization purposes (Optional)
colnames(TumorDesign) <- c("SampleID", "Pheno", "Color")

# Combine the design for both groups.
Design <- rbind(TumorDesign, NormalDesign)

Samples <- as.matrix(Design$SampleID)
Pheno <- as.matrix(Design$Pheno)
Color <- as.matrix(Design$Color)

LUSCdesign <- model.matrix(~ 0 + factor(c(Pheno)))
colnames(LUSCdesign) <- c("Normal", "Tumor")

####################################################################
####################### DATA CLEANUP AND FILTERING #################
####################################################################
# Removing rows (genes) with zero counts across all samples.

LUSCgenomicMatrix <- round(LUSCgenomicMatrix)
LUSCRSEM.zero <- LUSCgenomicMatrix[!apply(LUSCgenomicMatrix, 1, function(row) all(row == 0)), ]
LUSCRSEM.50zero <- LUSCRSEM.zero[rowSums(LUSCRSEM.zero != 0) > ncol(LUSCRSEM.zero)/2, ]
LUSCRSEM.80zero <- LUSCRSEM.50zero[rowSums(LUSCRSEM.50zero == 0)/ncol(LUSCRSEM.50zero) < 0.8, ]

####################################################################
####################### NORMALIZATION AND SCALING ##################
####################################################################
# Filter the dataset based on expression threshold, removing low-expressed genes.

LUSCRSEMNvT <- LUSCgenomicMatrix[, Samples]
LUSCRSEMNvT <- as.matrix(LUSCRSEMNvT)
isexpr <- rowSums(cpm(LUSCRSEMNvT) > 10) >= 2
LUSCRSEMNvT <- LUSCRSEMNvT[isexpr,]

# Scaling and normalizing counts for further differential expression analysis.
y <- DGEList(counts = LUSCRSEMNvT, lib.size = colSums(LUSCRSEMNvT))
y <- calcNormFactors(y)

# Perform Voom normalization, preparing for linear modeling.
LUSCvoomData <- voom(y, LUSCdesign, plot = TRUE)

####################################################################
##################### DIFFERENTIAL EXPRESSION ######################
####################################################################
# Linear modeling and contrast creation for differential expression analysis.

fit <- lmFit(LUSCvoomData, LUSCdesign)
contrast.matrix <- makeContrasts(Tumor - Normal, levels = LUSCdesign)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extracting the top differentially expressed genes (DEGs).
LUSCTCGATvNDEG <- topTable(fit2, adjust = "fdr", number = nrow(LUSCvoomData))

# Adding a "score" based on -log(p-value) and fold change.
LUSCTCGATvNDEG$Score <- (-log(LUSCTCGATvNDEG$adj.P.Val) * LUSCTCGATvNDEG$logFC)
LUSCTCGATvNDEG$ABScore <- abs(LUSCTCGATvNDEG$Score)

# Saving the DEG results as a CSV file.
write.csv(LUSCTCGATvNDEG, "LUSCTCGATvNDEG.csv")

####################################################################
######################## PATHWAY SCORING ###########################
####################################################################
# This section calculates pathway scores for specific pathways based 
# on the differential expression data.

AllPathwaysKEGGPS <- read.delim("C:/Users/bruep/OneDrive/Área de Trabalho/Kleiton/Pathways/Artigo_suplementares/Supplementary_Dataset_2_txt.txt", header=FALSE)
AllPathwaysKEGGPS <- data.frame(lapply(AllPathwaysKEGGPS, as.character), stringsAsFactors = FALSE)

# Function to calculate the absolute pathway score.
LUSCABSScoreSum <- function(x){
  ABSPathwayScore <- list()
  for(i in 1:ncol(AllPathwaysKEGGPS)) {
    y <- x[rownames(x) %in% as.character(AllPathwaysKEGGPS[,i]), ]
    score <- sum(y[,8]) / sqrt(length(Samples))
    ABSPathwayScore[[colnames(AllPathwaysKEGGPS)[i]]] <- score
    print(score)
    plot(score)
  }
  return(ABSPathwayScore)
}

# Apply absolute pathway score calculation.
ABSPathwayScore <- LUSCABSScoreSum(LUSCTCGATvNDEG)

####################################################################
################## POSITIVE/NEGATIVE PATHWAY SCORES ################
####################################################################
# Function to calculate positive/negative pathway score based on logFC.

LUSCScoreSum <- function(x){
  PathwayScore <- list()
  for(i in 1:ncol(AllPathwaysKEGGPS)) {
    y <- x[rownames(x) %in% as.character(AllPathwaysKEGGPS[,i]), ]
    score <- sum(y[,7]) / sqrt(length(Samples))
    PathwayScore[[colnames(AllPathwaysKEGGPS)[i]]] <- score
    print(score)
    plot(score)
  }
  return(PathwayScore)
}

# Apply pathway score calculation.
PathwayScore <- LUSCScoreSum(LUSCTCGATvNDEG)

# Save the pathway score results.
write.csv(t(as.data.frame(ABSPathwayScore)), file = "ABSPathwayScore.csv", quote = FALSE)
write.csv(t(as.data.frame(PathwayScore)), file = "PathwayScore.csv", quote = FALSE)

####################################################################
############# BOOTSTRAP P-VALUES FOR PATHWAY SCORES ################
####################################################################
# Bootstrap resampling of pathway scores to calculate p-values.

set.seed(1234)

RemainingABSScores <- data.frame(Pathway = character(), ABSScore = numeric())

for (i in 1:ncol(AllPathwaysKEGGPS)) {
  pathway_genes <- as.character(AllPathwaysKEGGPS[, i])
  y <- LUSCTCGATvNDEG[rownames(LUSCTCGATvNDEG) %in% pathway_genes, ]
  ABSPathwayScore <- sum(y$ABScore) / sqrt(length(pathway_genes))
  RemainingABSScores <- rbind(RemainingABSScores, data.frame(Pathway = colnames(AllPathwaysKEGGPS)[i], ABSScore = ABSPathwayScore))
}


# Resample DEG 100000 times to compute empirical p-values.
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval1 <- sum(z >= RemainingABSScores[1, 2]) / 100000
print(pval1)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval2 <- sum(z >= RemainingABSScores[2, 2]) / 100000
print(pval2)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval3 <- sum(z >= RemainingABSScores[3, 2]) / 100000
print(pval3)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval4 <- sum(z >= RemainingABSScores[4, 2]) / 100000
print(pval4)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval5 <- sum(z >= RemainingABSScores[5, 2]) / 100000
print(pval5)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval6 <- sum(z >= RemainingABSScores[6, 2]) / 100000
print(pval6)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval7 <- sum(z >= RemainingABSScores[7, 2]) / 100000
print(pval7)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval8 <- sum(z >= RemainingABSScores[8, 2]) / 100000
print(pval8)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval9 <- sum(z >= RemainingABSScores[9, 2]) / 100000
print(pval9)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval10 <- sum(z >= RemainingABSScores[10, 2]) / 100000
print(pval10)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval11 <- sum(z >= RemainingABSScores[11, 2]) / 100000
print(pval11)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval12 <- sum(z >= RemainingABSScores[12, 2]) / 100000
print(pval12)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval13 <- sum(z >= RemainingABSScores[13, 2]) / 100000
print(pval13)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval14 <- sum(z >= RemainingABSScores[14, 2]) / 100000
print(pval14)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval15 <- sum(z >= RemainingABSScores[15, 2]) / 100000
print(pval15)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval16 <- sum(z >= RemainingABSScores[16, 2]) / 100000
print(pval16)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval17 <- sum(z >= RemainingABSScores[17, 2]) / 100000
print(pval17)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval18 <- sum(z >= RemainingABSScores[18, 2]) / 100000
print(pval18)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval19 <- sum(z >= RemainingABSScores[19, 2]) / 100000
print(pval19)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval20 <- sum(z >= RemainingABSScores[20, 2]) / 100000
print(pval20)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval21 <- sum(z >= RemainingABSScores[21, 2]) / 100000
print(pval21)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval22 <- sum(z >= RemainingABSScores[22, 2]) / 100000
print(pval22)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval23 <- sum(z >= RemainingABSScores[23, 2]) / 100000
print(pval23)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval24 <- sum(z >= RemainingABSScores[24, 2]) / 100000
print(pval24)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval25 <- sum(z >= RemainingABSScores[25, 2]) / 100000
print(pval25)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval26 <- sum(z >= RemainingABSScores[26, 2]) / 100000
print(pval26)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval27 <- sum(z >= RemainingABSScores[27, 2]) / 100000
print(pval27)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval28 <- sum(z >= RemainingABSScores[28, 2]) / 100000
print(pval28)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval29 <- sum(z >= RemainingABSScores[29, 2]) / 100000
print(pval29)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval30 <- sum(z >= RemainingABSScores[30, 2]) / 100000
print(pval30)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval31 <- sum(z >= RemainingABSScores[31, 2]) / 100000
print(pval31)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval32 <- sum(z >= RemainingABSScores[32, 2]) / 100000
print(pval32)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval33 <- sum(z >= RemainingABSScores[33, 2]) / 100000
print(pval33)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval34 <- sum(z >= RemainingABSScores[34, 2]) / 100000
print(pval34)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval35 <- sum(z >= RemainingABSScores[35, 2]) / 100000
print(pval35)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval36 <- sum(z >= RemainingABSScores[36, 2]) / 100000
print(pval36)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval37 <- sum(z >= RemainingABSScores[37, 2]) / 100000
print(pval37)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval38 <- sum(z >= RemainingABSScores[38, 2]) / 100000
print(pval38)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval39 <- sum(z >= RemainingABSScores[39, 2]) / 100000
print(pval39)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval40 <- sum(z >= RemainingABSScores[40, 2]) / 100000
print(pval40)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval41 <- sum(z >= RemainingABSScores[41, 2]) / 100000
print(pval41)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval42 <- sum(z >= RemainingABSScores[42, 2]) / 100000
print(pval42)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval43 <- sum(z >= RemainingABSScores[43, 2]) / 100000
print(pval43)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval44 <- sum(z >= RemainingABSScores[44, 2]) / 100000
print(pval44)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval45 <- sum(z >= RemainingABSScores[45, 2]) / 100000
print(pval45)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval46 <- sum(z >= RemainingABSScores[46, 2]) / 100000
print(pval46)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval47 <- sum(z >= RemainingABSScores[47, 2]) / 100000
print(pval47)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval48 <- sum(z >= RemainingABSScores[48, 2]) / 100000
print(pval48)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval49 <- sum(z >= RemainingABSScores[49, 2]) / 100000
print(pval49)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval50 <- sum(z >= RemainingABSScores[50, 2]) / 100000
print(pval50)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval51 <- sum(z >= RemainingABSScores[51, 2]) / 100000
print(pval51)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval52 <- sum(z >= RemainingABSScores[52, 2]) / 100000
print(pval52)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval53 <- sum(z >= RemainingABSScores[53, 2]) / 100000
print(pval53)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval54 <- sum(z >= RemainingABSScores[54, 2]) / 100000
print(pval54)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval55 <- sum(z >= RemainingABSScores[55, 2]) / 100000
print(pval55)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval56 <- sum(z >= RemainingABSScores[56, 2]) / 100000
print(pval56)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval57 <- sum(z >= RemainingABSScores[57, 2]) / 100000
print(pval57)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval58 <- sum(z >= RemainingABSScores[58, 2]) / 100000
print(pval58)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval59 <- sum(z >= RemainingABSScores[59, 2]) / 100000
print(pval59)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval60 <- sum(z >= RemainingABSScores[60, 2]) / 100000
print(pval60)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval61 <- sum(z >= RemainingABSScores[61, 2]) / 100000
print(pval61)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval62 <- sum(z >= RemainingABSScores[62, 2]) / 100000
print(pval62)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval63 <- sum(z >= RemainingABSScores[63, 2]) / 100000
print(pval63)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval64 <- sum(z >= RemainingABSScores[64, 2]) / 100000
print(pval64)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval65 <- sum(z >= RemainingABSScores[65, 2]) / 100000
print(pval65)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval66 <- sum(z >= RemainingABSScores[66, 2]) / 100000
print(pval66)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval67 <- sum(z >= RemainingABSScores[67, 2]) / 100000
print(pval67)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval68 <- sum(z >= RemainingABSScores[68, 2]) / 100000
print(pval68)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval69 <- sum(z >= RemainingABSScores[69, 2]) / 100000
print(pval69)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval70 <- sum(z >= RemainingABSScores[70, 2]) / 100000
print(pval70)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval71 <- sum(z >= RemainingABSScores[71, 2]) / 100000
print(pval71)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval72 <- sum(z >= RemainingABSScores[72, 2]) / 100000
print(pval72)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval73 <- sum(z >= RemainingABSScores[73, 2]) / 100000
print(pval73)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval74 <- sum(z >= RemainingABSScores[74, 2]) / 100000
print(pval74)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval75 <- sum(z >= RemainingABSScores[75, 2]) / 100000
print(pval75)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval76 <- sum(z >= RemainingABSScores[76, 2]) / 100000
print(pval76)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval77 <- sum(z >= RemainingABSScores[77, 2]) / 100000
print(pval77)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval78 <- sum(z >= RemainingABSScores[78, 2]) / 100000
print(pval78)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval79 <- sum(z >= RemainingABSScores[79, 2]) / 100000
print(pval79)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval80 <- sum(z >= RemainingABSScores[80, 2]) / 100000
print(pval80)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval81 <- sum(z >= RemainingABSScores[81, 2]) / 100000
print(pval81)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval82 <- sum(z >= RemainingABSScores[82, 2]) / 100000
print(pval82)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval83 <- sum(z >= RemainingABSScores[83, 2]) / 100000
print(pval83)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval84 <- sum(z >= RemainingABSScores[84, 2]) / 100000
print(pval84)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval85 <- sum(z >= RemainingABSScores[85, 2]) / 100000
print(pval85)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval86 <- sum(z >= RemainingABSScores[86, 2]) / 100000
print(pval86)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval87 <- sum(z >= RemainingABSScores[87, 2]) / 100000
print(pval87)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval88 <- sum(z >= RemainingABSScores[88, 2]) / 100000
print(pval88)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval89 <- sum(z >= RemainingABSScores[89, 2]) / 100000
print(pval89)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval90 <- sum(z >= RemainingABSScores[90, 2]) / 100000
print(pval90)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval91 <- sum(z >= RemainingABSScores[91, 2]) / 100000
print(pval91)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval92 <- sum(z >= RemainingABSScores[92, 2]) / 100000
print(pval92)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval93 <- sum(z >= RemainingABSScores[93, 2]) / 100000
print(pval93)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval94 <- sum(z >= RemainingABSScores[94, 2]) / 100000
print(pval94)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval95 <- sum(z >= RemainingABSScores[95, 2]) / 100000
print(pval95)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval96 <- sum(z >= RemainingABSScores[96, 2]) / 100000
print(pval96)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval97 <- sum(z >= RemainingABSScores[97, 2]) / 100000
print(pval97)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval98 <- sum(z >= RemainingABSScores[98, 2]) / 100000
print(pval98)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval99 <- sum(z >= RemainingABSScores[99, 2]) / 100000
print(pval99)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval100 <- sum(z >= RemainingABSScores[100, 2]) / 100000
print(pval100)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval101 <- sum(z >= RemainingABSScores[101, 2]) / 100000
print(pval101)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval102 <- sum(z >= RemainingABSScores[102, 2]) / 100000
print(pval102)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval103 <- sum(z >= RemainingABSScores[103, 2]) / 100000
print(pval103)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval104 <- sum(z >= RemainingABSScores[104, 2]) / 100000
print(pval104)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval105 <- sum(z >= RemainingABSScores[105, 2]) / 100000
print(pval105)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval106 <- sum(z >= RemainingABSScores[106, 2]) / 100000
print(pval106)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval107 <- sum(z >= RemainingABSScores[107, 2]) / 100000
print(pval107)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval108 <- sum(z >= RemainingABSScores[108, 2]) / 100000
print(pval108)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval109 <- sum(z >= RemainingABSScores[109, 2]) / 100000
print(pval109)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval110 <- sum(z >= RemainingABSScores[110, 2]) / 100000
print(pval110)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval111 <- sum(z >= RemainingABSScores[111, 2]) / 100000
print(pval111)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval112 <- sum(z >= RemainingABSScores[112, 2]) / 100000
print(pval112)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval113 <- sum(z >= RemainingABSScores[113, 2]) / 100000
print(pval113)

## Resample DE genes 100000 times
z <- vector()
for (i in 1:100000){
  # Sample genes from the pool of DEG
  boot.ABSScores <- sample(LUSCTCGATvNDEG$ABScore, size=length(6), replace = T)
  # Sum ABSScores
  z[i] <- sum(boot.ABSScores) / sqrt(103)
}

pval114 <- sum(z >= RemainingABSScores[114, 2]) / 100000
print(pval114)


Pvals <- as.data.frame(cbind(pval1, pval2, pval3, pval4, pval5, pval6, pval7, pval8, pval9, pval10, pval11, pval12, pval13, pval14, pval15, pval16, pval17, pval18, pval19, pval20, pval21, pval22, pval23, pval24, pval25, pval26, pval27, pval28, pval29, pval30, pval31, pval32, pval33, pval34, pval35, pval36, pval37, pval38, pval39, pval40, pval41, pval42, pval43, pval44, pval45, pval46, pval47, pval48, pval49, pval50, pval51, pval52, pval53, pval54, pval55, pval56, pval57, pval58, pval59, pval60, pval61, pval62, pval63, pval64, pval65, pval66, pval67, pval68, pval69, pval70, pval71, pval72, pval73, pval74, pval75, pval76, pval77, pval78, pval79, pval80, pval81, pval82, pval83, pval84, pval85, pval86, pval87, pval88, pval89, pval90, pval91, pval92, pval93, pval94, pval95, pval96, pval97, pval98, pval99, pval100, pval101, pval102, pval103, pval104, pval105, pval106, pval107, pval108, pval109, pval110, pval111, pval112, pval113, pval114))
###CAN ADD MORE COLUMNS TO PVALS PRIOR TO EXPORT
###Remove Glycolysis/Gluconeogenesis
write.csv(Pvals, "PValues.csv")