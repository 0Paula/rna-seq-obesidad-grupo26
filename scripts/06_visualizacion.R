# HEATMAP CORRECTO PARA DATOS NORMALIZADOS
library(pheatmap)
library(RColorBrewer)

# 1. Cargar datos NORMALIZADOS (ya están procesados)
counts <- read.delim("C:/Users/Usuario/Downloads/matriz_genes_muestras_NORMALIZADA.tsv", row.names=1)

# 2. NO aplicar log2 - los datos ya están en escala adecuada
samples <- c("abraham", "homer", "bart", "lisa", "maggie")
heat_data <- counts[, samples]  # Sin transformación extra

# 3. Anotaciones
annot <- data.frame(
  Grupo = c("Obeso", "Obeso", "Normopeso", "Normopeso", "Normopeso"),
  row.names = colnames(heat_data)
)

ann_colors <- list(
  Grupo = c(Obeso = "#E74C3C", Normopeso = "#3498DB")
)

# 4. Heatmap CORRECTO
dev.new(width=12, height=10)  # Abre ventana

pheatmap(as.matrix(heat_data),
         scale = "row",  # ¡SÍ ESCALAR! Para comparación entre genes
         annotation_col = annot,
         annotation_colors = ann_colors,
         main = "Perfil de expresión génica entre\nindividuos obesos y normopeso (Normalizado)",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 9,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         border_color = NA,
         labels_col = c("Abraham (O)", "Homer (O)", "Bart (N)", "Lisa (N)", "Maggie (N)"))

# VERSIÓN CON PHEATMAP PERO SIN ERRORES
library(pheatmap)
library(RColorBrewer)

counts <- read.delim("C:/Users/Usuario/Downloads/matriz_genes_muestras.tsv", row.names=1)
samples <- c("abraham", "homer", "bart", "lisa", "maggie")
heat_data <- log2(counts[, samples] + 0.1)  # Transformar antes

annot <- data.frame(
  Grupo = c("Obeso", "Obeso", "Normopeso", "Normopeso", "Normopeso"),
  row.names = colnames(heat_data)
)

ann_colors <- list(
  Grupo = c(Obeso = "#E74C3C", Normopeso = "#3498DB")
)

# Heatmap SIN SCALE (para evitar errores)
pheatmap(as.matrix(heat_data),
         scale = "none",  # Cambiado de "row" a "none"
         annotation_col = annot,
         annotation_colors = ann_colors,
         main = "Perfil de expresión génica entre\nindividuos obesos y normopeso (Datos Raw)",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 9,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         border_color = NA,
         labels_col = c("Abraham (O)", "Homer (O)", "Bart (N)", "Lisa (N)", "Maggie (N)"))

# Volcano plot con raw 
library(ggplot2)
library(dplyr)

# Primero abrir device 6
dev.new(width = 10, height = 8)

# O usar windows() específicamente para Windows
# windows(width = 10, height = 8)

deseq_results <- read.delim("C:/Users/Usuario/Downloads/Results_DESeq2.tsv")

deseq_clean <- deseq_results %>%
  rename_with(~ gsub("\\.", "_", .x)) %>%
  mutate(
    log10p = -log10(pvalue),
    significance = case_when(
      padj < 0.05 & log2FoldChange > 0.5 ~ "UP",
      padj < 0.05 & log2FoldChange < -0.5 ~ "DOWN",
      TRUE ~ "NS"
    )
  )

volcano <- ggplot(deseq_clean, aes(x = log2FoldChange, y = log10p)) +
  geom_point(aes(color = significance), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NS" = "gray70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray40", linewidth = 0.8) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Análisis de expresión diferencial entre \nindividuos obesos y normopeso (DESeq2_raw)",
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Significancia"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5,
                              lineheight = 1.1, margin = margin(b = 10)),
    plot.margin = margin(20, 15, 15, 15, "pt")
  )

# Forzar visualización
print(volcano)

# Opción alternativa: usar plot directamente
plot(volcano)

ggsave("volcano_final.png", volcano, width = 10, height = 8, dpi = 300)

# Volcano plot con norm 

  deseq_results <- deseq_norm
  file_used <- "Results_DESeq2_norm.tsv"
  
# Preparar datos
deseq_clean <- deseq_results %>%
  # Asegurar nombres de columnas
  rename_with(~ gsub("\\.", "_", .x)) %>%
  # Calcular -log10(p)
  mutate(
    log10p = -log10(pvalue),
    significance = case_when(
      padj < 0.05 & log2FoldChange > 0.8 ~ "UP",
      padj < 0.05 & log2FoldChange < -0.8 ~ "DOWN",
      TRUE ~ "NS"
     )
   )

# Volcano Plot con título en dos líneas
cat("\n1️ CREANDO VOLCANO PLOT...\n")
volcano <- ggplot(deseq_clean, aes(x = log2FoldChange, y = log10p)) +
  geom_point(aes(color = significance), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NS" = "gray70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40", linewidth = 0.8) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Análisis de expresión diferencial entre \nindividuos obesos y normopeso (DESeq2_norm)",
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Significancia"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5,
                              lineheight = 1.1, margin = margin(b = 10)),
    plot.margin = margin(20, 15, 15, 15, "pt")
  )

ggsave("volcano_final.png", volcano, width = 10, height = 8, dpi = 300)
cat(" volcano_final.png guardado\n")
print(volcano)

# VOLCANO EDGER CON LOS 9 GENES SIGNIFICATIVOS
library(ggplot2)
library(dplyr)

# 1. Cargar EdgeR
edger <- read.delim("C:/Users/Usuario/Downloads/EdgeR_results_obesotipo1_vs_normopeso.tsv")

# 2. Preparar datos con umbral 0.5 (para capturar los 9 genes)
edger_clean <- edger %>%
  rename(
    log2FoldChange = logFC,
    pvalue = PValue,
    padj = FDR
  ) %>%
  mutate(
    log10p = -log10(pvalue),
    # UMBRAL 0.5 para capturar todos los 9 genes significativos
    significance = case_when(
      padj < 0.05 & log2FoldChange > 0.5 ~ "UP",
      padj < 0.05 & log2FoldChange < -0.5 ~ "DOWN",
      TRUE ~ "NS"
    )
  )

# 3. CONTAR - AHORA SÍ 9 GENES
cat(" ANÁLISIS FINAL EDGER:\n")
cat("• Total genes:", nrow(edger_clean), "\n")
cat("• padj < 0.05:", sum(edger_clean$padj < 0.05, na.rm = TRUE), "\n")
cat("• UP (log2FC > 0.5):", sum(edger_clean$significance == "UP"), "\n")
cat("• DOWN (log2FC < -0.5):", sum(edger_clean$significance == "DOWN"), "\n")
cat("• TOTAL significativos:", sum(edger_clean$significance != "NS"), "\n")

# Mostrar los 9 genes
sig_genes <- edger_clean %>%
  filter(padj < 0.05) %>%
  arrange(desc(log2FoldChange))

cat("\n LOS 9 GENES SIGNIFICATIVOS EDGER:\n")
for (i in 1:nrow(sig_genes)) {
  gene_name <- rownames(sig_genes)[i]
  cat(sprintf("%d. %-6s: log2FC = %6.3f, padj = %.2e %s\n",
              i, gene_name, 
              sig_genes$log2FoldChange[i], 
              sig_genes$padj[i],
              ifelse(abs(sig_genes$log2FoldChange[i]) > 0.8, "(>0.8)", "(0.5-0.8)")))
}

# 4. VOLCANO CON ESTILO IDÉNTICO Y LOS 9 GENES
dev.new(width = 10, height = 8)

volcano_edger <- ggplot(edger_clean, aes(x = log2FoldChange, y = log10p)) +
  geom_point(aes(color = significance), alpha = 0.6, size = 2.5) +
  scale_color_manual(
    name = "Significancia",
    values = c("UP" = "#E74C3C",     # Mismo rojo
               "DOWN" = "#3498DB",   # Mismo azul
               "NS" = "#BDC3C7")     # Mismo gris
  ) +
  # Líneas de umbral 0.5 (para capturar los 9 genes)
  geom_vline(xintercept = c(-0.5, 0.5), 
             linetype = "dashed", 
             color = "#2C3E50", 
             alpha = 0.7,
             linewidth = 0.8) +
  # También mostrar línea de 0.8 para referencia
  geom_vline(xintercept = c(-0.8, 0.8), 
             linetype = "dotted", 
             color = "gray60", 
             alpha = 0.5,
             linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed", 
             color = "#2C3E50", 
             alpha = 0.7,
             linewidth = 0.8) +
  # Tema idéntico
  theme_minimal(base_size = 14) +
  labs(
    title = "Análisis de expresión diferencial entre\nindividuos obesos y normopeso (EdgeR)",
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Significancia"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5,
                              lineheight = 1.2, margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 12),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    plot.margin = margin(20, 20, 15, 15, "pt")
  ) +

# Mostrar
print(volcano_edger)

# 5. Guardar
ggsave("volcano_EdgeR_9genes.png", 
       volcano_edger, 
       width = 10, 
       height = 8, 
       dpi = 300,
       bg = "white")

cat("\n Volcano EdgeR con 9 genes guardado: volcano_EdgeR_9genes.png\n")
cat(" Umbral: |log2FC| > 0.5 para capturar todos los genes significativos\n")

# PCA: NORMPESO vs OBESO TIPO 1 (Abraham y Homer)
library(ggplot2)

# 1. Cargar datos
counts_norm <- read.delim("C:/Users/Usuario/Downloads/matriz_genes_muestras_NORMALIZADA.tsv", row.names=1)

# 2. FILTRAR solo las muestras relevantes: Normopeso + Obeso tipo 1
muestras_relevantes <- c("abraham", "homer", "bart", "lisa", "maggie")
counts_filtrado <- counts_norm[, muestras_relevantes]

cat(" MUESTRAS ANALIZADAS EN PCA:\n")
cat("• Obeso tipo 1: Abraham, Homer\n")
cat("• Normopeso: Bart, Lisa, Maggie\n")
cat("• Total:", ncol(counts_filtrado), "muestras\n")

# 3. Calcular PCA
pca_result <- prcomp(t(counts_filtrado), scale. = TRUE)

# 4. Preparar datos
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$Muestra <- rownames(pca_df)
pca_df$Grupo <- ifelse(pca_df$Muestra %in% c("abraham", "homer"), 
                       "Obeso tipo 1", "Normopeso")

# 5. Varianza explicada
var_exp <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

# 6. PCA ESPECÍFICO - CORREGIDO PARA ETIQUETAS VISIBLES
pca_especifico <- ggplot(pca_df, aes(x = PC1, y = PC2, 
                                     color = Grupo, 
                                     shape = Grupo,
                                     label = Muestra)) +
  geom_point(size = 8, alpha = 0.9) +
  
  # ETIQUETAS CON AJUSTES DE POSICIÓN
  geom_text(color = "black", 
            vjust = -0.8,  # Reducir distancia vertical
            hjust = 0.5,   # Centrar horizontalmente
            fontface = "bold", 
            size = 4.5,    # Tamaño moderado
            check_overlap = FALSE,
            position = position_nudge(y = 0.02)) +  # Pequeño ajuste
  
  # Elipses de confianza 95%
  stat_ellipse(aes(fill = Grupo), 
               geom = "polygon", 
               alpha = 0.2, 
               level = 0.95,
               show.legend = FALSE) +
  
  # Escalas
  scale_color_manual(values = c("Obeso tipo 1" = "#E74C3C", 
                                "Normopeso" = "#3498DB")) +
  scale_fill_manual(values = c("Obeso tipo 1" = "#E74C3C", 
                               "Normopeso" = "#3498DB")) +
  scale_shape_manual(values = c("Obeso tipo 1" = 17,  # Triángulo
                                "Normopeso" = 19)) +   # Círculo
  
  # Tema y etiquetas con MÁRGENES AMPLIADOS
  theme_minimal(base_size = 16) +
  labs(
    title = "Análisis de componentes principales de la expresión génica\nentre individuos obesos y normopeso",
    x = paste0("PC1 (", var_exp[1], "% varianza)"),
    y = paste0("PC2 (", var_exp[2], "% varianza)"),
    color = "Condición",
    shape = "Condición"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5,
                              margin = margin(b = 15)),
    plot.subtitle = element_text(size = 14, color = "gray40", hjust = 0.5,
                                 margin = margin(b = 20)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 13),
    legend.box = "vertical",
    legend.spacing.y = unit(0.5, "cm"),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(size = 14),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    # MÁRGENES MÁS GRANDES para evitar cortes
    plot.margin = margin(40, 40, 40, 40, "pt"),
    # Ajustar límites del gráfico
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  # EXPANDIR LÍMITES del gráfico para dar espacio a las etiquetas
  expand_limits(y = c(min(pca_df$PC2) * 1.15, max(pca_df$PC2) * 1.15),
                x = c(min(pca_df$PC1) * 1.15, max(pca_df$PC1) * 1.15))

# 7. Mostrar en ventana grande
dev.new(width = 14, height = 11)  # Ventana más grande
print(pca_especifico)

# 8. Guardar con dimensiones ampliadas
ggsave("pca_obesotipo1_vs_normopeso.png", 
       pca_especifico, 
       width = 14,  # Más ancho
       height = 11, # Más alto
       dpi = 300,
       bg = "white")

cat("\n PCA específico guardado: pca_obesotipo1_vs_normopeso.png\n")
