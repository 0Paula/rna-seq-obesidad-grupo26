setwd("~/omicas")
# 1) Instalación (solo una vez)
install.packages("readr")
install.packages("dplyr")
install.packages("tibble")

install.packages("BiocManager")
BiocManager::install("edgeR")

# 2) Carga de librerías (siempre)
library(readr)
library(dplyr)
library(tibble)
library(edgeR)

counts <- read.delim(
  "matriz_genes_muestras.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

dim(counts)     # debería ser 37 x 5
head(counts)
colnames(counts)

design <- read.delim(
  "design.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

head(design)
colnames(design)

colnames(counts) <- trimws(tolower(colnames(counts)))
design$sample <- trimws(tolower(design$sample))
design$condition <- trimws(tolower(design$condition))

# Si tu design tuviera más filas, aquí filtramos solo las que están en counts
design_sub <- design %>%
  filter(sample %in% colnames(counts)) %>%
  select(sample, condition)

# Reordenar design_sub para que coincida con el orden de counts
design_sub <- design_sub %>%
  arrange(match(sample, colnames(counts)))

# Comprobación: si falla aquí, no sigas
stopifnot(all(design_sub$sample == colnames(counts)))

design_sub

design_sub$condition <- factor(
  design_sub$condition,
  levels = c("normopeso", "obesotipo1")
)

table(design_sub$condition)  # debe dar 3 y 2

y <- DGEList(counts = counts, group = design_sub$condition)
y

keep <- filterByExpr(y, group = design_sub$condition)
y <- y[keep, , keep.lib.sizes = FALSE]
dim(y)

# Normalización
y <- calcNormFactors(y)

# Matriz de diseño del modelo
design_mat <- model.matrix(~ design_sub$condition)

# Estimar dispersión
y <- estimateDisp(y, design_mat)

# Ajuste y test (QL recomendado)
fit <- glmQLFit(y, design_mat)
qlf <- glmQLFTest(fit, coef = 2)   # Obeso1 vs Normopeso


res <- topTags(qlf, n = Inf)$table
res

# Añadir “significativo”
res$Significant <- res$FDR < 0.05
table(res$Significant)

# Guardar
write.table(
  res,
  file = "EdgeR_results_obesotipo1_vs_normopeso.tsv",
  sep = "\t",
  quote = FALSE
)




