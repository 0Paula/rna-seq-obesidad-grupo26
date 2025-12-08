# rna-seq- Secuenciación y ómicas de próxima generación
Análisis de expresión diferencial obesidad vs normopeso (RNA-seq)

## Estructura

- `scripts/`: scripts en R para cada etapa del pipeline.
- `data/`: datos de entrada ligeros (matrices de conteo). **No subir FASTQ ni BAM.**
- `results/`: resultados del análisis (tablas, figuras).
- `docs/`: documentación y materiales para el póster. **Revisión bibliográfica**
- `poster/`: contenido del póster.

## Pipeline

1. Control de calidad (FastQC, MultiQC).
2. Pseudoalineamiento / cuantificación (Salmon).
3. Cuantificación por gen y matriz de conteos (tximport en R).
4. Análisis de expresión diferencial (DESeq2 y EdgeR).
5. Visualización (volcano, heatmap, PCA).
6. Interpretación biológica y póster.
