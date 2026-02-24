library(data.table)

samples <- fread("C:/Users/sjaya/Downloads/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt")
subjects <- fread("C:/Users/sjaya/Downloads/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt")

# filtering 
samples_gtex <- samples[grepl("^GTEX", samples$SAMPID), ]

# extract Extract SUBJID From SAMPID
samples$SUBJID <- sapply(
  strsplit(samples$SAMPID, "-"),
  function(x) paste(x[1:2], collapse = "-")
)

head(samples[, c("SAMPID", "SUBJID")])


# merge samples and subjects
library(dplyr)

metadata <- left_join(samples, subjects, by = "SUBJID")

# filtering 
samples_gtex <- samples[grepl("^GTEX", samples$SAMPID), ]

# extract SUBJID
samples_gtex$SUBJID <- sapply(
  strsplit(samples_gtex$SAMPID, "-"),
  function(x) paste(x[1:2], collapse = "-")
)
# merge again
metadata <- left_join(samples_gtex, subjects, by = "SUBJID")

# age to numeric
metadata$AGE_NUM <- as.numeric(sub("-.*", "", metadata$AGE))

# midpoint
metadata$AGE_NUM <- sapply(metadata$AGE, function(x) {
  parts <- strsplit(x, "-")[[1]]
  mean(as.numeric(parts))
})

# filter only blood
blood_meta <- metadata %>%
  filter(SMTSD == "Whole Blood")

# header only of large gtex file
header <- fread(
  "C:/Users/sjaya/Downloads/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz",
  skip = 2,
  nrows = 0
)

# blood samples present in expressionn file
expr_sample_ids <- colnames(header)[-c(1,2)]

common_blood_ids <- intersect(expr_sample_ids, blood_ids)

length(common_blood_ids)

# load only required columns
expr_blood <- fread(
  "C:/Users/sjaya/Downloads/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz",
  skip = 2,
  select = c("Name", common_blood_ids)
)

# transpose

# Remove gene ID column
X <- as.matrix(expr_blood[, -1, with=FALSE])

# Set gene names as rownames
rownames(X) <- expr_blood$Name

# Transpose
X <- t(X)

# Align meta data
blood_meta_aligned <- blood_meta[
  match(rownames(X), blood_meta$SAMPID),
]
dim(X)

# compute correlation
cor_results <- apply(X, 2, function(gene) {
  cor(gene, y)
})
y <- blood_meta_aligned$AGE_NUM

# rank genes by strength of assocation

# Rank genes by absolute correlation
gene_ranking <- sort(abs(cor_results), decreasing = TRUE)

# Look at top 10
head(gene_ranking, 10)

# build 500 gene model
# Select top 500 genes
top_500_genes <- names(gene_ranking)[1:500]

X_500 <- X[, top_500_genes]

dim(X_500)

# log transform
X_500_log <- log2(X_500 + 1)

# train/test split
set.seed(123)

n <- nrow(X_500_log)
train_idx <- sample(1:n, size = 0.8*n)

X_train <- X_500_log[train_idx, ]
y_train <- y[train_idx]

X_test  <- X_500_log[-train_idx, ]
y_test  <- y[-train_idx]

# elastic net model
library(glmnet)

cv_model_500 <- cv.glmnet(
  X_train,
  y_train,
  alpha = 0.5,
  nfolds = 10
)

pred_500 <- predict(
  cv_model_500,
  X_test,
  s = "lambda.min"
)

R2_500   <- cor(pred_500, y_test)^2
RMSE_500 <- sqrt(mean((pred_500 - y_test)^2))
MAE_500  <- mean(abs(pred_500 - y_test))

R2_500
RMSE_500
MAE_500

# build 1500 gene model
# Select top 1500 genes
top_1500_genes <- names(gene_ranking)[1:1500]

X_1500 <- X[, top_1500_genes]
X_1500_log <- log2(X_1500 + 1)

# Train/test using SAME split
X_train_1500 <- X_1500_log[train_idx, ]
X_test_1500  <- X_1500_log[-train_idx, ]

cv_model_1500 <- cv.glmnet(
  X_train_1500,
  y_train,
  alpha = 0.5,
  nfolds = 10
)

pred_1500 <- predict(
  cv_model_1500,
  X_test_1500,
  s = "lambda.min"
)

R2_1500   <- cor(pred_1500, y_test)^2
RMSE_1500 <- sqrt(mean((pred_1500 - y_test)^2))
MAE_1500  <- mean(abs(pred_1500 - y_test))

R2_1500
RMSE_1500
MAE_1500

# build 300 gene model
# Select top 300 genes
top_300_genes <- names(gene_ranking)[1:300]

X_300 <- X[, top_300_genes]
X_300_log <- log2(X_300 + 1)

# Use SAME split as before
X_train_300 <- X_300_log[train_idx, ]
X_test_300  <- X_300_log[-train_idx, ]

cv_model_300 <- cv.glmnet(
  X_train_300,
  y_train,
  alpha = 0.5,
  nfolds = 10
)

pred_300 <- predict(
  cv_model_300,
  X_test_300,
  s = "lambda.min"
)

R2_300   <- cor(pred_300, y_test)^2
RMSE_300 <- sqrt(mean((pred_300 - y_test)^2))
MAE_300  <- mean(abs(pred_300 - y_test))

R2_300
RMSE_300
MAE_300

# compute age acceleration
# Get predictions for ALL samples using best model (500 genes)
pred_full_500 <- predict(
  cv_model_500,
  X_500_log,
  s = "lambda.min"
)

age_accel_500 <- as.numeric(pred_full_500) - y

summary(age_accel_500)

# visualization
plot(y, pred_full_500,
     xlab = "Actual Age",
     ylab = "Predicted Age",
     main = "500-Gene Blood Aging Model")

abline(0, 1, col = "red", lwd = 2)

# tune alpha properly
alphas <- seq(0, 1, 0.25)
results_alpha <- data.frame()

for (a in alphas) {
  
  set.seed(123)
  
  cv_model <- cv.glmnet(
    X_train,
    y_train,
    alpha = a,
    nfolds = 10
  )
  
  pred <- predict(cv_model, X_test, s = "lambda.min")
  
  R2 <- cor(pred, y_test)^2
  RMSE <- sqrt(mean((pred - y_test)^2))
  MAE <- mean(abs(pred - y_test))
  
  results_alpha <- rbind(
    results_alpha,
    data.frame(alpha = a,
               R2 = R2,
               RMSE = RMSE,
               MAE = MAE)
  )
}

results_alpha

# 5 repeated splits
set.seed(123)

R2_values <- c()

for(i in 1:5){
  
  train_idx <- sample(1:nrow(X_500_log), size = 0.8*nrow(X_500_log))
  
  X_train <- X_500_log[train_idx, ]
  y_train <- y[train_idx]
  
  X_test <- X_500_log[-train_idx, ]
  y_test <- y[-train_idx]
  
  cv_model <- cv.glmnet(
    X_train,
    y_train,
    alpha = 1,
    nfolds = 10
  )
  
  pred <- predict(cv_model, X_test, s = "lambda.min")
  
  R2 <- cor(pred, y_test)^2
  
  R2_values <- c(R2_values, R2)
}

mean(R2_values)
sd(R2_values)

# multi tissue filtering
selected_tissues <- c("Whole Blood", "Liver", "Muscle - Skeletal")

metadata_multi <- metadata %>%
  filter(SMTSD %in% selected_tissues)

table(metadata_multi$SMTSD)

# sample id of multitissue
multi_ids <- metadata_multi$SAMPID
length(multi_ids)

# expression for multi tissue
expr_multi <- fread(
  "C:/Users/sjaya/Downloads/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz",
  skip = 2,
  select = c("Name", multi_ids)
)

# transpose
X_multi <- as.matrix(expr_multi[, -1, with = FALSE])
rownames(X_multi) <- expr_multi$Name
X_multi <- t(X_multi)

dim(X_multi)

# align metadata
metadata_multi_aligned <- metadata_multi[
  match(rownames(X_multi), metadata_multi$SAMPID),
]

all(rownames(X_multi) == metadata_multi_aligned$SAMPID)

# define y for multi tissue
y_multi <- metadata_multi_aligned$AGE_NUM

# correlation
cor_multi <- apply(X_multi, 2, function(gene) {
  cor(gene, y_multi)
})

gene_ranking_multi <- sort(abs(cor_multi), decreasing = TRUE)

# build 500 multi tissue model
top_500_multi <- names(gene_ranking_multi)[1:500]

X_500_multi <- X_multi[, top_500_multi]
X_500_multi_log <- log2(X_500_multi + 1)

# train/ test model
set.seed(123)

n_multi <- nrow(X_500_multi_log)
train_idx_multi <- sample(1:n_multi, size = 0.8*n_multi)

X_train_multi <- X_500_multi_log[train_idx_multi, ]
y_train_multi <- y_multi[train_idx_multi]

X_test_multi  <- X_500_multi_log[-train_idx_multi, ]
y_test_multi  <- y_multi[-train_idx_multi]

cv_multi <- cv.glmnet(
  X_train_multi,
  y_train_multi,
  alpha = 1,
  nfolds = 10
)

pred_multi <- predict(cv_multi, X_test_multi, s = "lambda.min")

R2_multi   <- cor(pred_multi, y_test_multi)^2
RMSE_multi <- sqrt(mean((pred_multi - y_test_multi)^2))
MAE_multi  <- mean(abs(pred_multi - y_test_multi))

R2_multi
RMSE_multi
MAE_multi

# check tissue effect
# Predict age for ALL samples
pred_full_multi <- predict(cv_multi, X_500_multi_log, s="lambda.min")

age_accel_multi <- as.numeric(pred_full_multi) - y_multi

metadata_multi_aligned$AgeAcceleration <- age_accel_multi

aggregate(AgeAcceleration ~ SMTSD,
          data = metadata_multi_aligned,
          mean)

# proper testing
anova_multi <- aov(AgeAcceleration ~ SMTSD,
                   data = metadata_multi_aligned)

summary(anova_multi)

# Cross-Tissue Validation
# Define train and test sets by tissue
train_idx_ct <- metadata_multi_aligned$SMTSD %in% 
  c("Liver", "Muscle - Skeletal")

test_idx_ct <- metadata_multi_aligned$SMTSD == 
  "Whole Blood"

# Prepare matrices
X_train_ct <- X_500_multi_log[train_idx_ct, ]
y_train_ct <- y_multi[train_idx_ct]

X_test_ct  <- X_500_multi_log[test_idx_ct, ]
y_test_ct  <- y_multi[test_idx_ct]

# Train model
cv_ct <- cv.glmnet(
  X_train_ct,
  y_train_ct,
  alpha = 1,
  nfolds = 10
)

pred_ct <- predict(cv_ct, X_test_ct, s = "lambda.min")

R2_ct   <- cor(pred_ct, y_test_ct)^2
RMSE_ct <- sqrt(mean((pred_ct - y_test_ct)^2))
MAE_ct  <- mean(abs(pred_ct - y_test_ct))

R2_ct
RMSE_ct
MAE_ct

# function for cross tissue validation
library(glmnet)

run_loto <- function(train_tissues, test_tissue) {
  
  train_idx <- metadata_multi_aligned$SMTSD %in% train_tissues
  test_idx  <- metadata_multi_aligned$SMTSD == test_tissue
  
  X_train <- X_500_multi_log[train_idx, ]
  y_train <- y_multi[train_idx]
  
  X_test  <- X_500_multi_log[test_idx, ]
  y_test  <- y_multi[test_idx]
  
  cv_model <- cv.glmnet(
    X_train,
    y_train,
    alpha = 1,
    nfolds = 10
  )
  
  pred <- predict(cv_model, X_test, s = "lambda.min")
  
  R2   <- cor(pred, y_test)^2
  RMSE <- sqrt(mean((pred - y_test)^2))
  MAE  <- mean(abs(pred - y_test))
  
  return(c(R2 = R2, RMSE = RMSE, MAE = MAE))
}


# Train: Liver + Muscle → Test: Blood
train_idx_1 <- metadata_multi_aligned$SMTSD %in%
  c("Liver", "Muscle - Skeletal")

test_idx_1 <- metadata_multi_aligned$SMTSD ==
  "Whole Blood"

X_train_1 <- X_500_multi_log[train_idx_1, ]
y_train_1 <- y_multi[train_idx_1]

X_test_1  <- X_500_multi_log[test_idx_1, ]
y_test_1  <- y_multi[test_idx_1]

cv_1 <- cv.glmnet(
  X_train_1,
  y_train_1,
  alpha = 1,
  nfolds = 10
)

pred_1 <- predict(cv_1, X_test_1, s = "lambda.min")

R2_1   <- cor(pred_1, y_test_1)^2
RMSE_1 <- sqrt(mean((pred_1 - y_test_1)^2))
MAE_1  <- mean(abs(pred_1 - y_test_1))

R2_1
RMSE_1
MAE_1

# Train: Blood + Muscle → Test: Liver
train_idx_2 <- metadata_multi_aligned$SMTSD %in%
  c("Whole Blood", "Muscle - Skeletal")

test_idx_2 <- metadata_multi_aligned$SMTSD ==
  "Liver"

X_train_2 <- X_500_multi_log[train_idx_2, ]
y_train_2 <- y_multi[train_idx_2]

X_test_2  <- X_500_multi_log[test_idx_2, ]
y_test_2  <- y_multi[test_idx_2]

cv_2 <- cv.glmnet(
  X_train_2,
  y_train_2,
  alpha = 1,
  nfolds = 10
)

pred_2 <- predict(cv_2, X_test_2, s = "lambda.min")

R2_2   <- cor(pred_2, y_test_2)^2
RMSE_2 <- sqrt(mean((pred_2 - y_test_2)^2))
MAE_2  <- mean(abs(pred_2 - y_test_2))

R2_2
RMSE_2
MAE_2

# Train: Blood + Liver → Test: Muscle
train_idx_3 <- metadata_multi_aligned$SMTSD %in%
  c("Whole Blood", "Liver")

test_idx_3 <- metadata_multi_aligned$SMTSD ==
  "Muscle - Skeletal"

X_train_3 <- X_500_multi_log[train_idx_3, ]
y_train_3 <- y_multi[train_idx_3]

X_test_3  <- X_500_multi_log[test_idx_3, ]
y_test_3  <- y_multi[test_idx_3]

cv_3 <- cv.glmnet(
  X_train_3,
  y_train_3,
  alpha = 1,
  nfolds = 10
)

pred_3 <- predict(cv_3, X_test_3, s = "lambda.min")

R2_3   <- cor(pred_3, y_test_3)^2
RMSE_3 <- sqrt(mean((pred_3 - y_test_3)^2))
MAE_3  <- mean(abs(pred_3 - y_test_3))

R2_3
RMSE_3
MAE_3

# everything together 
loto_results <- data.frame(
  Tissue = c("Blood","Liver","Muscle"),
  R2  = c(R2_1, R2_2, R2_3),
  RMSE = c(RMSE_1, RMSE_2, RMSE_3),
  MAE  = c(MAE_1, MAE_2, MAE_3)
)

loto_results

# solid tissues only
# solid tissues only
selected_tissues <- c(
  "Liver",
  "Muscle - Skeletal",
  "Lung"
)

metadata_solid <- metadata[
  metadata$SMTSD %in% selected_tissues,
]

table(metadata_solid$SMTSD)
# sample id
solid_ids <- metadata_solid$SAMPID
length(solid_ids)

# Load expression
expr_solid <- fread(
  "C:/Users/sjaya/Downloads/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz",
  skip = 2,
  select = c("Name", solid_ids)
)

# transpose
X_solid <- as.matrix(expr_solid[, -1, with = FALSE])
rownames(X_solid) <- expr_solid$Name
X_solid <- t(X_solid)

dim(X_solid)

# align metadata
metadata_solid_aligned <- metadata_solid[
  match(rownames(X_solid), metadata_solid$SAMPID),
]

all(rownames(X_solid) == metadata_solid_aligned$SAMPID)

# define age vector
y_solid <- metadata_solid_aligned$AGE_NUM

# rank gene by age correlation
cor_solid <- apply(X_solid, 2, function(gene) {
  cor(gene, y_solid)
})

gene_ranking_solid <- sort(abs(cor_solid), decreasing = TRUE)

top_500_solid <- names(gene_ranking_solid)[1:500]

X_500_solid <- X_solid[, top_500_solid]
X_500_solid_log <- log2(X_500_solid + 1)

dim(X_500_solid_log)

# multi tissue model for solid
set.seed(123)

n_solid <- nrow(X_500_solid_log)
train_idx_solid <- sample(1:n_solid, size = 0.8*n_solid)

X_train_solid <- X_500_solid_log[train_idx_solid, ]
y_train_solid <- y_solid[train_idx_solid]

X_test_solid  <- X_500_solid_log[-train_idx_solid, ]
y_test_solid  <- y_solid[-train_idx_solid]

cv_solid <- cv.glmnet(
  X_train_solid,
  y_train_solid,
  alpha = 1,
  nfolds = 10
)

pred_solid <- predict(cv_solid, X_test_solid, s = "lambda.min")

R2_solid   <- cor(pred_solid, y_test_solid)^2
RMSE_solid <- sqrt(mean((pred_solid - y_test_solid)^2))
MAE_solid  <- mean(abs(pred_solid - y_test_solid))

R2_solid
RMSE_solid
MAE_solid

# metabolic tissue selection
selected_tissues <- c(
  "Liver",
  "Heart - Left Ventricle",
  "Muscle - Skeletal"
)

metadata_metabolic <- metadata[
  metadata$SMTSD %in% selected_tissues,
]

table(metadata_metabolic$SMTSD)

# load expression
metabolic_ids <- metadata_metabolic$SAMPID

expr_metabolic <- fread(
  "C:/Users/sjaya/Downloads/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz",
  skip = 2,
  select = c("Name", metabolic_ids)
)

X_metabolic <- as.matrix(expr_metabolic[, -1, with = FALSE])
rownames(X_metabolic) <- expr_metabolic$Name
X_metabolic <- t(X_metabolic)

dim(X_metabolic)

metadata_metabolic_aligned <- metadata_metabolic[
  match(rownames(X_metabolic), metadata_metabolic$SAMPID),
]

all(rownames(X_metabolic) == metadata_metabolic_aligned$SAMPID)

y_metabolic <- metadata_metabolic_aligned$AGE_NUM

# rank and select top 500
cor_metabolic <- apply(X_metabolic, 2, function(gene) {
  cor(gene, y_metabolic)
})

gene_ranking_metabolic <- sort(abs(cor_metabolic), decreasing = TRUE)

top_500_metabolic <- names(gene_ranking_metabolic)[1:500]

X_500_metabolic <- X_metabolic[, top_500_metabolic]
X_500_metabolic_log <- log2(X_500_metabolic + 1)

# multi tissue model (metabolic)
set.seed(123)

n_meta <- nrow(X_500_metabolic_log)
train_idx_meta <- sample(1:n_meta, size = 0.8*n_meta)

X_train_meta <- X_500_metabolic_log[train_idx_meta, ]
y_train_meta <- y_metabolic[train_idx_meta]

X_test_meta  <- X_500_metabolic_log[-train_idx_meta, ]
y_test_meta  <- y_metabolic[-train_idx_meta]

cv_meta <- cv.glmnet(
  X_train_meta,
  y_train_meta,
  alpha = 1,
  nfolds = 10
)

pred_meta <- predict(cv_meta, X_test_meta, s="lambda.min")

R2_meta   <- cor(pred_meta, y_test_meta)^2
RMSE_meta <- sqrt(mean((pred_meta - y_test_meta)^2))
MAE_meta  <- mean(abs(pred_meta - y_test_meta))

R2_meta
RMSE_meta
MAE_meta

# Train Liver + Muscle → Test Heart
train_idx_A <- metadata_metabolic_aligned$SMTSD %in%
c("Liver","Muscle - Skeletal")

test_idx_A <- metadata_metabolic_aligned$SMTSD ==
  "Heart - Left Ventricle"

X_train_A <- X_500_metabolic_log[train_idx_A, ]
y_train_A <- y_metabolic[train_idx_A]

X_test_A  <- X_500_metabolic_log[test_idx_A, ]
y_test_A  <- y_metabolic[test_idx_A]

cv_A <- cv.glmnet(X_train_A, y_train_A, alpha=1, nfolds=10)

pred_A <- predict(cv_A, X_test_A, s="lambda.min")

R2_A_meta <- cor(pred_A, y_test_A)^2
RMSE_A_meta <- sqrt(mean((pred_A - y_test_A)^2))
MAE_A_meta <- mean(abs(pred_A - y_test_A))

R2_A_meta
RMSE_A_meta
MAE_A_meta

# Train Liver + Heart → Test Muscle
train_idx_B <- metadata_metabolic_aligned$SMTSD %in%
  c("Liver","Heart - Left Ventricle")

test_idx_B <- metadata_metabolic_aligned$SMTSD ==
  "Muscle - Skeletal"

X_train_B <- X_500_metabolic_log[train_idx_B, ]
y_train_B <- y_metabolic[train_idx_B]

X_test_B  <- X_500_metabolic_log[test_idx_B, ]
y_test_B  <- y_metabolic[test_idx_B]

cv_B <- cv.glmnet(X_train_B, y_train_B, alpha=1, nfolds=10)

pred_B <- predict(cv_B, X_test_B, s="lambda.min")

R2_B_meta <- cor(pred_B, y_test_B)^2
RMSE_B_meta <- sqrt(mean((pred_B - y_test_B)^2))
MAE_B_meta <- mean(abs(pred_B - y_test_B))

R2_B_meta
RMSE_B_meta
MAE_B_meta

# Train Muscle + Heart → Test Liver
train_idx_C <- metadata_metabolic_aligned$SMTSD %in%
  c("Muscle - Skeletal","Heart - Left Ventricle")

test_idx_C <- metadata_metabolic_aligned$SMTSD ==
  "Liver"

X_train_C <- X_500_metabolic_log[train_idx_C, ]
y_train_C <- y_metabolic[train_idx_C]

X_test_C  <- X_500_metabolic_log[test_idx_C, ]
y_test_C  <- y_metabolic[test_idx_C]

cv_C <- cv.glmnet(X_train_C, y_train_C, alpha=1, nfolds=10)

pred_C <- predict(cv_C, X_test_C, s="lambda.min")

R2_C_meta <- cor(pred_C, y_test_C)^2
RMSE_C_meta <- sqrt(mean((pred_C - y_test_C)^2))
MAE_C_meta <- mean(abs(pred_C - y_test_C))

R2_C_meta
RMSE_C_meta
MAE_C_meta

# fibroblast multi tissue
selected_tissues <- c(
  "Cells - Cultured fibroblasts",
  "Liver",
  "Muscle - Skeletal"
)

metadata_cells <- metadata[
  metadata$SMTSD %in% selected_tissues,
]

table(metadata_cells$SMTSD)

# load expression
cell_ids <- metadata_cells$SAMPID

expr_cells <- fread(
  "C:/Users/sjaya/Downloads/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz",
  skip = 2,
  select = c("Name", cell_ids)
)

X_cells <- as.matrix(expr_cells[, -1, with = FALSE])
rownames(X_cells) <- expr_cells$Name
X_cells <- t(X_cells)

dim(X_cells)

metadata_cells_aligned <- metadata_cells[
  match(rownames(X_cells), metadata_cells$SAMPID),
]

all(rownames(X_cells) == metadata_cells_aligned$SAMPID)


cor_cells <- apply(X_cells, 2, function(gene) {
  cor(gene, y_cells)
})

gene_ranking_cells <- sort(abs(cor_cells), decreasing = TRUE)

top_500_cells <- names(gene_ranking_cells)[1:500]

X_500_cells <- X_cells[, top_500_cells]
X_500_cells_log <- log2(X_500_cells + 1)
y_cells <- metadata_cells_aligned$AGE_NUM

# multi tissue model
set.seed(123)

n_cells <- nrow(X_500_cells_log)
train_idx_cells <- sample(1:n_cells, size = 0.8*n_cells)

X_train_cells <- X_500_cells_log[train_idx_cells, ]
y_train_cells <- y_cells[train_idx_cells]

X_test_cells  <- X_500_cells_log[-train_idx_cells, ]
y_test_cells  <- y_cells[-train_idx_cells]

cv_cells <- cv.glmnet(
  X_train_cells,
  y_train_cells,
  alpha = 1,
  nfolds = 10
)

pred_cells <- predict(cv_cells, X_test_cells, s="lambda.min")

R2_cells   <- cor(pred_cells, y_test_cells)^2
RMSE_cells <- sqrt(mean((pred_cells - y_test_cells)^2))
MAE_cells  <- mean(abs(pred_cells - y_test_cells))

R2_cells
RMSE_cells
MAE_cells

# Train Liver + Muscle → Test Fibroblasts
train_idx_A <- metadata_cells_aligned$SMTSD %in%
  c("Liver","Muscle - Skeletal")

test_idx_A <- metadata_cells_aligned$SMTSD ==
  "Cells - Cultured fibroblasts"

X_train_A <- X_500_cells_log[train_idx_A, ]
y_train_A <- y_cells[train_idx_A]

X_test_A  <- X_500_cells_log[test_idx_A, ]
y_test_A  <- y_cells[test_idx_A]

cv_A <- cv.glmnet(X_train_A, y_train_A, alpha=1, nfolds=10)

pred_A <- predict(cv_A, X_test_A, s="lambda.min")

R2_A_cells <- cor(pred_A, y_test_A)^2
R2_A_cells

# Train Fibroblasts + Muscle → Test Liver
train_idx_B <- metadata_cells_aligned$SMTSD %in%
  c("Cells - Cultured fibroblasts","Muscle - Skeletal")

test_idx_B <- metadata_cells_aligned$SMTSD ==
  "Liver"

X_train_B <- X_500_cells_log[train_idx_B, ]
y_train_B <- y_cells[train_idx_B]

X_test_B  <- X_500_cells_log[test_idx_B, ]
y_test_B  <- y_cells[test_idx_B]

cv_B <- cv.glmnet(X_train_B, y_train_B, alpha=1, nfolds=10)

pred_B <- predict(cv_B, X_test_B, s="lambda.min")

R2_B_cells <- cor(pred_B, y_test_B)^2
R2_B_cells

# Train Fibroblasts + Liver → Test Muscle
train_idx_C <- metadata_cells_aligned$SMTSD %in%
  c("Cells - Cultured fibroblasts","Liver")

test_idx_C <- metadata_cells_aligned$SMTSD ==
  "Muscle - Skeletal"

X_train_C <- X_500_cells_log[train_idx_C, ]
y_train_C <- y_cells[train_idx_C]

X_test_C  <- X_500_cells_log[test_idx_C, ]
y_test_C  <- y_cells[test_idx_C]

cv_C <- cv.glmnet(X_train_C, y_train_C, alpha=1, nfolds=10)

pred_C <- predict(cv_C, X_test_C, s="lambda.min")

R2_C_cells <- cor(pred_C, y_test_C)^2
R2_C_cells

# pathway analysis
# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)

# Load annotation file (if not already loaded)

gene_annotation <- fread(
  "C:/Users/sjaya/Downloads/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_tpm.gct.gz",
  skip = 2,
  select = c("Name", "Description")
)

gene_annotation$Gene_clean <- sub("\\..*", "", gene_annotation$Name)

head(gene_annotation)

# Remove version numbers
gene_annotation$Gene_clean <- sub("\\..*", "", gene_annotation$Name)

# solid cluster enrichment
# Clean top genes
top_500_solid_clean <- sub("\\..*", "", top_500_multi)

# Map to symbols
solid_symbols <- gene_annotation$Description[
  match(top_500_solid_clean, gene_annotation$Gene_clean)
]

solid_symbols <- solid_symbols[!is.na(solid_symbols)]

# fibroblast enrichment
top_500_cells_clean <- sub("\\..*", "", top_500_cells)

cell_symbols <- gene_annotation$Description[
  match(top_500_cells_clean, gene_annotation$Gene_clean)
]

cell_symbols <- cell_symbols[!is.na(cell_symbols)]

cell_entrez <- bitr(
  cell_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

ego_cells <- enrichGO(
  gene = cell_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

head(ego_cells)