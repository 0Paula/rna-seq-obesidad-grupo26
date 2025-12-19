# ------------------------------------------------------------------------------
# 1. Configuración inicial
# ------------------------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

library(DESeq2)

# ------------------------------------------------------------------------------
# 2. Descargamos los datos
# ------------------------------------------------------------------------------

setwd("C:/Users/lidia/Desktop/Bioinformática/OMI - Secuenciación y Ómicas de Próxima Generación/mubio03_act2/TallerGrupal_Ficheros")

# Cargamos la matriz de conteos (genes en filas, muestras en columnas)
counts <- read.delim("matriz_genes_muestras.tsv", row.names = 1)
#counts <- read.delim("matriz_genes_muestras_NORMALIZADA.tsv", row.names = 1)

# Cargamos el diseño experimental (metadatos)
metadata <- read.csv("Design.csv")

# ------------------------------------------------------------------------------
# 3. Preparación de datos para DESeq2
# ------------------------------------------------------------------------------

# Convertimos counts a enteros: DESeq2 requiere números enteros
counts <- round(counts)

# Alineamos los nombres de las muestras de metadata con los de counts
target_samples <- c("AbrahamSimpson", "BartSimpson", "HomerSimpson", "LisaSimpson", "MaggieSimpson")
metadata_filt <- metadata[metadata$Sample %in% target_samples, ]

# Reordenamos las filas de metadata según las columnas de counts
rownames(metadata_filt) <- c("abraham", "homer", "bart", "lisa", "maggie")
metadata_ordered <- metadata_filt[colnames(counts), ]

# Verificación
all(rownames(metadata_ordered) == colnames(counts))

# Aseguramos que la columna Condition sea un factor
metadata_ordered$Condition <- factor(metadata_ordered$Condition)

# Establecemos el nivel de referencia (el denominador de la comparación)
metadata_ordered$Condition <- relevel(metadata_ordered$Condition, ref = "Normopeso")


# ------------------------------------------------------------------------------
# 4. Análisis diferencial
# ------------------------------------------------------------------------------

# Creamos un objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata_ordered,
                              design = ~ Condition)

# Ejecutamos DESeq2
dds <- DESeq(dds)

# Extraemos los resultados para la comparación Obeso1 vs Normopeso
res <- results(dds, contrast=c("Condition", "Sobrepeso/Obeso1", "Normopeso"))

# Ordenamos por el valor p ajustado
res_ordered <- res[order(res$padj), ]

summary(res)

# ------------------------------------------------------------------------------
# 5. Identificación de Genes UP y DOWN regulados
# ------------------------------------------------------------------------------

res_sig <- subset(res_ordered, padj < 0.05)
genes_up <- subset(res_sig, log2FoldChange > 0) # Genes UP (Más expresado en Obesos)
genes_down <- subset(res_sig, log2FoldChange < 0) # Genes DOWN (Menos expresado en Obesos)

cat("Número de genes diferencialmente expresados (padj < 0.05):", nrow(res_sig), "\n")
cat("Genes UP-regulados:", nrow(genes_up), "\n")
cat("Genes DOWN-regulados:", nrow(genes_down), "\n")

# ------------------------------------------------------------------------------
# 6. Preparación de la tabla final
# ------------------------------------------------------------------------------

df_final <- as.data.frame(res_ordered)

# Añadimos columna de clasificación simple
df_final$diff_expressed <- "NO"
df_final$diff_expressed[df_final$log2FoldChange > 0 & df_final$padj < 0.05] <- "UP"
df_final$diff_expressed[df_final$log2FoldChange < 0 & df_final$padj < 0.05] <- "DOWN"

head(df_final)

write.csv(df_final, "Results_DESeq2.csv")
#write.csv(df_final, "Results_DESeq2_norm.csv")

write.table(df_final, 
            file = "Results_DESeq2.tsv", 
            #file = "Results_DESeq2_norm.tsv", 
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = NA)

