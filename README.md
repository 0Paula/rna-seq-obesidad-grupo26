# rna-seq- Secuenciación y ómicas de próxima generación
Análisis de expresión diferencial obesidad vs normopeso (RNA-seq)

## Estructura

- `scripts/`: scripts en R para cada etapa del pipeline.
- `data/`: datos de entrada ligeros (matrices de conteo). **No subir FASTQ ni BAM.**
- `results/`: resultados del análisis (tablas, figuras).
- `docs/`: documentación y materiales para el póster. **Revisión bibliográfica**
- `poster/`: contenido del póster.

## Pipeline 

1️⃣ Control de calidad (QC)
📂 Archivos de entrada

Carpeta Fastqs/

Archivos .fastq.gz de las 8 muestras (paired-end, R1 y R2).

Se utilizarán especialmente las muestras pertenecientes a:

Grupo Normopeso (Bart, Lisa y Maggie Simpson)

Grupo Obeso tipo 1 (Abraham y Homer Simpson)

🔧 Qué se hace

Evaluación de la calidad de las lecturas:

calidad por base

longitud de lectura

contenido GC

presencia de adaptadores

Recorte de lecturas solo si se detectan problemas de calidad relevantes.

📤 Resultados

Informes de calidad individuales (FastQC).

Informe conjunto (MultiQC).

Figuras de QC para el póster.

🎯 Objetivo

Asegurar que los datos de secuenciación son adecuados para comparar Normopeso vs Obeso tipo 1.

2️⃣ Pseudoalineamiento / Alineamiento
📂 Archivos de entrada

Referencia.fasta

FASTQ de las muestras Normopeso y Obeso tipo 1.

🔧 Qué se hace

Creación del índice a partir de la referencia.

Ejecución del pseudoalineamiento (Salmon) o alineamiento (STAR) para cada muestra.

📤 Resultados

Salmon:

Carpeta por muestra con archivo quant.sf.

STAR (si se usa):

Archivos BAM y métricas de alineamiento (mapping rate, fragment length, GC bias).

🎯 Objetivo

Asignar las lecturas a transcritos/genes para poder cuantificar la expresión en ambos grupos.

3️⃣ Cuantificación por gen y preparación de la matriz
📂 Archivos de entrada

Archivos quant.sf de cada muestra.

Archivo Transcrito_a_Gen.tsv.

🔧 Qué se hace

Conversión de abundancias por transcrito a conteos por gen.

Construcción de una matriz de expresión:

filas → genes

columnas → muestras Normopeso y Obeso tipo 1.

📤 Resultados

Matriz de conteos Genes × Muestras, lista para análisis estadístico.

Matriz evaluada con y sin normalización (comparativa metodológica).

🎯 Objetivo

Obtener la matriz base para analizar diferencias de expresión entre Normopeso y Obeso tipo 1.

4️⃣ Análisis de expresión diferencial con DESeq2
Normopeso vs Obeso tipo 1
📂 Archivos de entrada

Matriz Genes × Muestras.

Design.csv, indicando la condición:

Normopeso

Obeso tipo 1.

🔧 Qué se hace

Definición del contraste experimental (Obeso tipo 1 vs Normopeso).

Ejecución del pipeline DESeq2:

estimación de dispersión

ajuste del modelo

corrección por FDR.

📤 Resultados

Tabla DESeq2 con:

log2FoldChange

p-valor ajustado

Identificación de genes:

sobreexpresados en Obeso tipo 1

infraexpresados respecto a Normopeso.

🎯 Objetivo

Detectar genes diferencialmente expresados asociados al fenotipo de obesidad tipo 1.

5️⃣ Análisis de expresión diferencial con EdgeR
Normopeso vs Obeso tipo 1
📂 Archivos de entrada

Misma matriz Genes × Muestras.

Design.csv.

🔧 Qué se hace

Creación del objeto DGEList.

Estimación de dispersión.

Ajuste del modelo estadístico para el contraste Obeso tipo 1 vs Normopeso.

📤 Resultados

Tabla de resultados EdgeR.

Comparación de genes significativos con DESeq2.

🎯 Objetivo

Validar la robustez de los resultados usando un método estadístico alternativo.

6️⃣ Visualización de resultados
📂 Archivos de entrada

Resultados de DESeq2 y EdgeR.

Listado de genes relacionados con obesidad (carpeta Genes/).

🔧 Qué se hace

Generación de figuras:

Volcano plot (Normopeso vs Obeso tipo 1).

Heatmap de genes de obesidad.

PCA para evaluar la separación entre muestras Normopeso y Obeso tipo 1.

📤 Resultados

Figuras finales en alta resolución para el póster.

🎯 Objetivo

Visualizar patrones de expresión diferencial entre ambos grupos.

7️⃣ Interpretación biológica y póster final
📂 Archivos de entrada

Resultados de expresión diferencial.

Carpeta Genes/.

Bases de datos externas (GeneCards, OMIM, PubMed).

🔧 Qué se hace

Interpretación funcional de los genes diferencialmente expresados.

Relación de los cambios de expresión con el fenotipo Obeso tipo 1 frente a Normopeso.

Redacción y diseño del póster científico.

🎯 Objetivo

Dar sentido biológico a las diferencias de expresión observadas y comunicar los resultados de forma clara.
