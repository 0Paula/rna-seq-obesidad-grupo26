
install.packages("tidyverse")
library(readr)

install.packages("BiocManager")
library(BiocManager)
BiocManager::install("tximport")
library(tximport)

readLines("Transcrito_a_Gen.tsv", n = 1)
print(readLines("Transcrito_a_Gen.tsv", n = 1))

getwd() #Para comprobar de donde estoy cogiendo el tsv
list.files() # para ver los elementos de mi carpeta


tx2gene <- read_tsv("Transcrito_a_Gen.tsv", col_names = FALSE) #coger el .tsv
tx2gene

tx2gene <- as.data.frame(tx2gene) #tximport funciona mejor con dataframe que con tibble, para asegurar que no es un tibble

colnames(tx2gene) <- c("Transcrito", "Gen_ID") #Renombro las columnas que se llamaban X1 y X2
colnames(tx2gene)

#Ahora le decimos dónde está cada quant.sf:

list.files("C:/Users/fernandezperez/Downloads/SECUENCIACION-ACTIVIDAD 2/Output_Parte_2_Lucia/outputs")

samples <- c("abraham", "bart", "homer", "lisa", "maggie")

ruta_base <- "C:/Users/fernandezperez/Downloads/SECUENCIACION-ACTIVIDAD 2/Output_Parte_2_Lucia/outputs"

files <- file.path(
  ruta_base,
  paste0(samples, "_quant.sf")
)


names(files) <- samples
files #PARA COMPROBAR
file.exists(files) #PARA COMPROBAR QUE TODOS DAN TRUE. IMPORTANTE

#Creación de la matriz Genes x Muestras
#Leer los datos por transcrito. Usar tx2gene para agrupar

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene
)

dim(txi$counts) #Le pido que me de la parte llamada counts que está dentro de tx.
#Me va a dar Filas(genes), Columnas (muestras) y Valores (cuántos reads hay para cada gen en cada muestra)

head(txi$counts)

#para guardar la tabla genes x muestras
write.table(
  txi$counts,
  file = "matriz_genes_muestras.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA,
  fileEncoding = "UTF-8"
)

#Repetir tximport con NORMALIZACION. 
#Lee los resultados por transcrito, los agrupa por gen y devuelve conteos normalizados.
txi_norm <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM"
)

dim(txi_norm$counts) #Mirar el resultado normalizado
head(txi_norm$counts)

#Para guardar la matriz NORMALIZADA
write.table(
  txi_norm$counts,
  file = "matriz_genes_muestras_NORMALIZADA.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA,
  fileEncoding = "UTF-8"
)

#Para visualizar los archivos .tsv en R
tabla <- read.delim("matriz_genes_muestras.tsv", sep = "\t", header = TRUE)
View(tabla)

tabla <- read.delim("matriz_genes_muestras_NORMALIZADA.tsv", sep = "\t", header = TRUE)
View(tabla)
