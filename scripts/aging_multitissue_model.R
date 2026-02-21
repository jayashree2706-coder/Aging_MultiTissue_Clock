# Multi-Tissue Aging Clock Model
# Author: Jayashree
# Project: Integrative Aging Transcriptomics

###############################################################
# TRANSCRIPTOMIC MODELING OF SYSTEMIC BIOLOGICAL AGING
# Multi-Tissue Elastic Net + Pharmacogene Extension
###############################################################

############################
# 1. LOAD LIBRARIES
############################

library(data.table)
library(dplyr)
library(glmnet)
library(clusterProfiler)
library(org.Hs.eg.db)

###############################################################
# 2. LOAD METADATA
###############################################################

samples <- fread("C:/Users/sjaya/Downloads/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt")
subjects <- fread("C:/Users/sjaya/Downloads/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt")

###############################################################
# 3. SELECT TISSUES (Whole Blood, Liver, Muscle)
###############################################################

selected_tissues <- c("Whole Blood", "Liver", "Muscle - Skeletal")

samples_filtered <- samples %>%
  filter(SMTSD %in% selected_tissues)

# Extract Subject ID from SAMPID
samples_filtered$SUBJID <- sapply(
  strsplit(samples_filtered$SAMPID, "-"),
  function(x) paste(x[1:2], collapse="-")
)

metadata <- samples_filtered %>%
  left_join(subjects, by="SUBJID")

# Convert AGE ranges to numeric
metadata$AGE_NUM <- as.numeric(sub("-.*", "", metadata$AGE))

# Remove missing ages
metadata_clean <- metadata[!is.na(metadata$AGE_NUM), ]

###############################################################
# 4. MATCH EXPRESSION SAMPLES WITH METADATA
###############################################################

header <- fread(
  "C:/Users/sjaya/Downloads/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz",
  skip=2,
  nrows=0
)

header_names <- colnames(header)

expr_sample_ids <- header_names[-c(1,2)]

common_ids <- intersect(expr_sample_ids, metadata_clean$SAMPID)

metadata_ml <- metadata_clean[metadata_clean$SAMPID %in% common_ids, ]

###############################################################
# 5. LOAD EXPRESSION DATA (TOP 300 AGE-ASSOCIATED GENES)
###############################################################

# top_genes must already be created from genome-wide age scan
# Example:
# top_genes <- results_sorted$Gene[1:300]

expr_full_ml <- fread(
  "C:/Users/sjaya/Downloads/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz",
  skip=2,
  select=c("Name", common_ids)
)

expr_subset_ml <- expr_full_ml[Name %in% top_genes]

X <- as.matrix(expr_subset_ml[, -1, with=FALSE])
rownames(X) <- expr_subset_ml$Name
X <- t(X)

y <- metadata_ml$AGE_NUM

# Align order
metadata_ml <- metadata_ml[match(rownames(X), metadata_ml$SAMPID), ]

###############################################################
# 6. TRAIN / TEST SPLIT (80/20)
###############################################################

set.seed(123)

n <- nrow(X)
train_idx <- sample(1:n, size=0.8*n)

X_train <- X[train_idx, ]
y_train <- y[train_idx]

X_test <- X[-train_idx, ]
y_test <- y[-train_idx]

###############################################################
# 7. ELASTIC NET MODEL
###############################################################

cv_enet_train <- cv.glmnet(
  X_train,
  y_train,
  alpha=0.5,
  nfolds=10,
  standardize=TRUE
)

pred_test <- predict(cv_enet_train, X_test, s="lambda.min")

R2_test <- cor(pred_test, y_test)^2
RMSE_test <- sqrt(mean((pred_test - y_test)^2))
MAE_test <- mean(abs(pred_test - y_test))

print(R2_test)
print(RMSE_test)
print(MAE_test)

###############################################################
# 8. SELECTED FEATURES
###############################################################

coef_enet <- coef(cv_enet_train, s="lambda.min")
selected_features <- rownames(coef_enet)[coef_enet[,1] != 0]
selected_features <- selected_features[selected_features != "(Intercept)"]

length(selected_features)

###############################################################
# 9. PHARMACOGENE ANALYSIS
###############################################################

pharmacogene_symbols <- c(
  "CYP3A4","CYP2D6","CYP2C9","CYP2C19",
  "CYP1A2","CYP2B6",
  "ABCB1","ABCC2",
  "SLCO1B1","SLCO1B3",
  "UGT1A1","UGT2B7",
  "VKORC1","DPYD"
)

gene_annotation <- fread(
  "C:/Users/sjaya/Downloads/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz",
  skip=2,
  select=c("Name","Description")
)

pharmacogene_list <- gene_annotation$Name[
  gene_annotation$Description %in% pharmacogene_symbols
]

length(pharmacogene_list)

# Extract pharmacogene expression
expr_pharma <- expr_full_ml[Name %in% pharmacogene_list]

X_pharma <- as.matrix(expr_pharma[, -1, with=FALSE])
rownames(X_pharma) <- expr_pharma$Name
X_pharma <- t(X_pharma)

# Train/test for pharmacogene-only model
set.seed(123)

n_p <- nrow(X_pharma)
train_idx_p <- sample(1:n_p, size=0.8*n_p)

X_train_p <- X_pharma[train_idx_p, ]
y_train_p <- y[train_idx_p]

X_test_p <- X_pharma[-train_idx_p, ]
y_test_p <- y[-train_idx_p]

cv_pharma <- cv.glmnet(X_train_p, y_train_p, alpha=0.5, nfolds=10)

pred_p <- predict(cv_pharma, X_test_p, s="lambda.min")

R2_p <- cor(pred_p, y_test_p)^2
RMSE_p <- sqrt(mean((pred_p - y_test_p)^2))
MAE_p <- mean(abs(pred_p - y_test_p))

print(R2_p)
print(RMSE_p)
print(MAE_p)

###############################################################
# 10. AGE ACCELERATION ANALYSIS
###############################################################

pred_full <- predict(cv_enet_train, X, s="lambda.min")
age_acceleration <- as.numeric(pred_full) - y

metadata_ml$AgeAcceleration <- age_acceleration

print(aggregate(AgeAcceleration ~ SMTSD, data=metadata_ml, mean))

anova_result <- aov(AgeAcceleration ~ SMTSD, data=metadata_ml)
print(summary(anova_result))

###############################################################
# 11. PATHWAY ENRICHMENT ANALYSIS
###############################################################

selected_annotation <- gene_annotation[
  gene_annotation$Name %in% selected_features, ]

selected_symbols <- selected_annotation$Description

gene_entrez <- bitr(selected_symbols,
                    fromType="SYMBOL",
                    toType="ENTREZID",
                    OrgDb=org.Hs.eg.db)

ego <- enrichGO(
  gene=gene_entrez$ENTREZID,
  OrgDb=org.Hs.eg.db,
  ont="BP",
  pAdjustMethod="BH",
  readable=TRUE
)

print(head(ego))

# Plot enrichment results (dotplot) for visualization
dotplot(ego, showCategory=10)

