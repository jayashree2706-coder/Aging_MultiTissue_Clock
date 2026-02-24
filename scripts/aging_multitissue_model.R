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

# =========================================
# 1000-GENE MODEL EXPERIMENT
# =========================================
# =====================================================
# 1️⃣ Create AgeGroup Variable
# =====================================================

metadata_ml$AgeGroup <- factor(
  metadata_ml$AGE,
  levels = c("20-29","30-39","40-49","50-59","60-69","70-79")
)

table(metadata_ml$AgeGroup)
sum(is.na(metadata_ml$AgeGroup))   # MUST be 0


# =====================================================
# 2️⃣ Remove Ensembl Version Numbers
# =====================================================

results$Gene_clean <- sub("\\..*", "", results$Gene)
expr_full_ml$Gene_clean <- sub("\\..*", "", expr_full_ml$Name)


# =====================================================
# 3️⃣ Select Top 1000 Age-Correlated Genes
# =====================================================

top_genes_1000_clean <- results$Gene_clean[1:1000]


# =====================================================
# 4️⃣ Extract 1000-Gene Expression Matrix
# =====================================================

expr_subset_1000 <- expr_full_ml[
  expr_full_ml$Gene_clean %in% top_genes_1000_clean
]

sample_columns <- setdiff(colnames(expr_subset_1000),
                          c("Name","Gene_clean"))

X_1000 <- as.matrix(expr_subset_1000[, ..sample_columns])
rownames(X_1000) <- expr_subset_1000$Gene_clean
X_1000 <- t(X_1000)

dim(X_1000)   # should be 1883 x 1000


# =====================================================
# 5️⃣ Map Ensembl → Gene Symbols
# =====================================================

gene_map <- gene_annotation[, c("Name","Description")]
gene_map$Gene_clean <- sub("\\..*", "", gene_map$Name)

colnames(X_1000) <- sub("\\..*", "", colnames(X_1000))

symbol_lookup <- gene_map$Description[
  match(colnames(X_1000), gene_map$Gene_clean)
]

colnames(X_1000) <- symbol_lookup


# =====================================================
# 6️⃣ Log2 Transform
# =====================================================

X_class <- log2(X_1000 + 1)

all(rownames(X_class) == metadata_ml$SAMPID)


# =====================================================
# 7️⃣ Compute Pathway Score
# =====================================================

pathway_gene_list <- unique(unlist(
  strsplit(ego@result$geneID[1:5], "/")
))

pathway_genes_present <- intersect(colnames(X_class),
                                   pathway_gene_list)

length(pathway_genes_present)

pathway_score <- rowMeans(
  X_class[, pathway_genes_present],
  na.rm = TRUE
)

summary(pathway_score)


# =====================================================
# 8️⃣ Compute Pharmacogene Score
# =====================================================

pharma_score <- rowMeans(log2(X_pharma + 1))
summary(pharma_score)


# =====================================================
# 9️⃣ Train/Test Split
# =====================================================

set.seed(123)

n <- nrow(X_class)
train_idx <- sample(1:n, size = 0.8*n)

X_train_gene <- X_class[train_idx, ]
X_test_gene  <- X_class[-train_idx, ]

y_train <- metadata_ml$AgeGroup[train_idx]
y_test  <- metadata_ml$AgeGroup[-train_idx]

length(y_train)
nrow(X_train_gene)


# =====================================================
# 🔟 Gene-Level Multinomial Classification
# =====================================================

library(glmnet)

cv_gene <- cv.glmnet(
  X_train_gene,
  y_train,
  family = "multinomial",
  alpha = 0.5
)

pred_gene <- predict(
  cv_gene,
  X_test_gene,
  s = "lambda.min",
  type = "class"
)

cat("\n--- 1000-Gene Confusion Matrix ---\n")
print(table(Predicted = pred_gene,
            Actual = y_test))

accuracy_gene_1000 <- mean(pred_gene == y_test)

cat("\n1000-Gene Accuracy:\n")
print(accuracy_gene_1000)


# =====================================================
# 11️⃣ Multi-Layer Classification
# =====================================================

X_multi <- cbind(
  X_class,
  PathwayScore = pathway_score,
  PharmaScore = pharma_score
)

X_train_multi <- X_multi[train_idx, ]
X_test_multi  <- X_multi[-train_idx, ]

cv_multi <- cv.glmnet(
  X_train_multi,
  y_train,
  family = "multinomial",
  alpha = 0.5
)

pred_multi <- predict(
  cv_multi,
  X_test_multi,
  s = "lambda.min",
  type = "class"
)

cat("\n--- 1000-Gene Multi-Layer Confusion Matrix ---\n")
print(table(Predicted = pred_multi,
            Actual = y_test))

accuracy_multi_1000 <- mean(pred_multi == y_test)

cat("\n1000-Gene Multi-Layer Accuracy:\n")
print(accuracy_multi_1000)

library(glmnet)

cv_multi <- cv.glmnet(
  X_train_multi,
  y_train,
  family = "multinomial",
  alpha = 0.5
)

pred_multi <- predict(
  cv_multi,
  X_test_multi,
  s = "lambda.min",
  type = "class"
)

cat("\n--- 1000-Gene + Pharma Confusion Matrix ---\n")
print(table(Predicted = pred_multi,
            Actual = y_test))

accuracy_multi_1000 <- mean(pred_multi == y_test)

cat("\n1000-Gene + Pharma Accuracy:\n")
print(accuracy_multi_1000)

# -----------------------------------------
# Cross-Tissue Split
# Train: Liver + Muscle
# Test: Whole Blood
# -----------------------------------------

train_idx_ct <- metadata_ml$SMTSD %in% c("Liver","Muscle - Skeletal")
test_idx_ct  <- metadata_ml$SMTSD == "Whole Blood"

table(metadata_ml$SMTSD[train_idx_ct])
table(metadata_ml$SMTSD[test_idx_ct])

X_ct_train <- X_log[train_idx_ct, ]
y_ct_train <- metadata_ml$AGE_NUM[train_idx_ct]

X_ct_test  <- X_log[test_idx_ct, ]
y_ct_test  <- metadata_ml$AGE_NUM[test_idx_ct]

library(glmnet)

set.seed(123)

cv_ct <- cv.glmnet(
  X_ct_train,
  y_ct_train,
  alpha = 0.5
)

pred_ct <- predict(
  cv_ct,
  X_ct_test,
  s = "lambda.min"
)

R2_ct   <- cor(pred_ct, y_ct_test)^2
RMSE_ct <- sqrt(mean((pred_ct - y_ct_test)^2))
MAE_ct  <- mean(abs(pred_ct - y_ct_test))

R2_ct
RMSE_ct
MAE_ct

###############################################################
# 12. CROSS-TISSUE EXTERNAL VALIDATION (REGRESSION)
###############################################################

# Use 300-gene log matrix (create if not already done)
X_log <- log2(X + 1)

# Train on Liver + Muscle
train_idx_ct <- metadata_ml$SMTSD %in% c("Liver","Muscle - Skeletal")

# Test on Whole Blood
test_idx_ct  <- metadata_ml$SMTSD == "Whole Blood"

# Split
X_ct_train <- X_log[train_idx_ct, ]
y_ct_train <- metadata_ml$AGE_NUM[train_idx_ct]

X_ct_test  <- X_log[test_idx_ct, ]
y_ct_test  <- metadata_ml$AGE_NUM[test_idx_ct]

# Train model
set.seed(123)

cv_ct <- cv.glmnet(
  X_ct_train,
  y_ct_train,
  alpha = 0.5,
  nfolds = 10
)

# Predict on Blood
pred_ct <- predict(cv_ct, X_ct_test, s = "lambda.min")

# Metrics
R2_ct   <- cor(pred_ct, y_ct_test)^2
RMSE_ct <- sqrt(mean((pred_ct - y_ct_test)^2))
MAE_ct  <- mean(abs(pred_ct - y_ct_test))

cat("\n--- Cross-Tissue Validation ---\n")
print(R2_ct)
print(RMSE_ct)
print(MAE_ct)

###############################################################
# 13. AGE ACCELERATION ANALYSIS
###############################################################

# Use global model trained on all tissues
pred_full <- predict(cv_enet_train, X_log, s = "lambda.min")

# Biological age deviation
metadata_ml$AgeAcceleration <-
  as.numeric(pred_full) - metadata_ml$AGE_NUM

# Summary
cat("\n--- Age Acceleration Summary ---\n")
print(summary(metadata_ml$AgeAcceleration))

# Tissue-level means
cat("\n--- Mean Age Acceleration by Tissue ---\n")
print(aggregate(AgeAcceleration ~ SMTSD,
                data = metadata_ml,
                mean))

# ANOVA
anova_result <- aov(AgeAcceleration ~ SMTSD,
                    data = metadata_ml)

cat("\n--- ANOVA Results ---\n")
print(summary(anova_result))

###############################################################
# 14. MULTI-LAYER LINEAR MODEL
###############################################################

# Transcriptomic predicted age
transcriptomic_score <- as.numeric(pred_full)

# Ensure pharma_score aligns
pharma_score <- rowMeans(log2(X_pharma + 1))

# Match order
pharma_score <- pharma_score[match(metadata_ml$SAMPID,
                                   rownames(X_pharma))]

# Fit linear model
lm_multi <- lm(
  AGE_NUM ~ transcriptomic_score + pharma_score,
  data = metadata_ml
)

cat("\n--- Multi-Layer Linear Model ---\n")
print(summary(lm_multi))
