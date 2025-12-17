
setwd("C:/Users/Lucia/Desktop/Actividad 2 (1)/mubio03_act2/TallerGrupal_Ficheros
      /TallerGrupal_Ficheros/Fastqs")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager") 

if (!require("Rsubread"))BiocManager::install("Rsubread")
library(Rsubread) # Library for alignment & mapping

if (!require("Rsamtools")) BiocManager::install("Rsamtools")
library(Rsamtools) # Library for Fragment Lenght

if (!require("GenomicAlignments")) BiocManager::install("GenomicAlignments")
library(GenomicAlignments) # Library for GC bias

# Referencia # 
buildindex(basename="my_genome_index", reference="Referencia.fasta")

# Homer Simpson (Obeso Groupo 1) # 
HomerSimpson_stats <- align(index="my_genome_index", readfile1="HomerSimpson_R1.fastq.gz", 
      readfile2="HomerSimpson_R2.fastq.gz", output_file="HomerSimpson_Aligned.bam")

# Mapping rate
print(HomerSimpson_stats)

mapped_percent_HomerSimpson <- HomerSimpson_stats[2, 1] / 
  HomerSimpson_stats[1, 1] * 100
print(paste("Mapping Rate:", round(mapped_percent_HomerSimpson, 2), "%"))

# Fragment Lenght
HomerSimpson_bam_file <- "HomerSimpson_Aligned.bam"

param <- ScanBamParam(what = c("isize"))
bam_data <- scanBam(HomerSimpson_bam_file, param = param)
frag_lengths <- abs(bam_data[[1]]$isize)

frag_lengths <- frag_lengths[frag_lengths > 50 & frag_lengths < 800]

summary(frag_lengths)
hist(frag_lengths, breaks=100, col="salmon", 
     main="Fragment Length Distribution Homer Simpson"
     , xlab="Size (bp")

# GC bias 
HomerSimpson_read_alignments <- readGAlignments("HomerSimpson_Aligned.bam", 
                                               param=ScanBamParam(what="seq"))
HomerSimpson_read_seqs <- mcols(HomerSimpson_read_alignments)$seq

gc_content <- letterFrequency(HomerSimpson_read_seqs, letters="GC", as.prob=TRUE)

hist(gc_content, breaks=50, col="pink", main="Read-level GC Content Homer Simpson", 
     xlab="GC Proportion (0.0 to 1.0")

# Abraham Simpson (Obeso Grupo 1) #
AbrahamSimpson_stats <- align(index="my_genome_index", 
                              readfile1="AbrahamSimpson_R1.fastq.gz",
                              readfile2="AbrahamSimpson_R2.fastq.gz", 
                              output_file="AbrahamSimpson_Aligned.bam")
      
# Mapping rate
print(AbrahamSimpson_stats)

mapped_percent_AbrahamSimpson <- AbrahamSimpson_stats[2, 1] / 
  AbrahamSimpson_stats[1, 1] * 100
print(paste("Mapping Rate:", round(mapped_percent_AbrahamSimpson, 2), "%"))

# Fragment Lenght
AbrahamSimpson_bam_file <- "AbrahamSimpson_Aligned.bam"

param <- ScanBamParam(what = c("isize"))
bam_data <- scanBam(AbrahamSimpson_bam_file, param = param)
frag_lengths <- abs(bam_data[[1]]$isize)

frag_lengths <- frag_lengths[frag_lengths > 50 & frag_lengths < 800]

summary(frag_lengths)
hist(frag_lengths, breaks=100, col="salmon", 
     main="Fragment Length Distribution Abraham Simpson"
     , xlab="Size (bp")

# GC bias 
AbrahamSimpson_read_alignments <- readGAlignments("AbrahamSimpson_Aligned.bam", 
                                               param=ScanBamParam(what="seq"))
AbrahamSimpson_read_seqs <- mcols(AbrahamSimpson_read_alignments)$seq

gc_content <- letterFrequency(AbrahamSimpson_read_seqs, letters="GC", as.prob=TRUE)

hist(gc_content, breaks=50, col="pink", main="Read-level GC Content Abraham Simpson", 
     xlab="GC Proportion (0.0 to 1.0")

# Bart Simpson (Normopeso) #
BartSimpson_stats <- align(index="my_genome_index", 
      readfile1="BartSimpson_R1.fastq.gz", 
      readfile2="BartSimpson_R2.fastq.gz", output_file="BartSimpson_Aligned.bam")

# Mapping rate
print(BartSimpson_stats)

mapped_percent_BartSimpson <- BartSimpson_stats[2, 1] / BartSimpson_stats[1, 1] * 100
print(paste("Mapping Rate:", round(mapped_percent_BartSimpson, 2), "%"))

# Fragment Lenght
BartSimpson_bam_file <- "BartSimpson_Aligned.bam"

param <- ScanBamParam(what = c("isize"))
bam_data <- scanBam(BartSimpson_bam_file, param = param)
frag_lengths <- abs(bam_data[[1]]$isize)

frag_lengths <- frag_lengths[frag_lengths > 50 & frag_lengths < 800]

summary(frag_lengths)
hist(frag_lengths, breaks=100, col="salmon",
     main="Fragment Length Distribution Bart Simpson"
     , xlab="Size (bp")

# GC bias 
BartSimpson_read_alignments <- readGAlignments("BartSimpson_Aligned.bam", 
                                               param=ScanBamParam(what="seq"))
BartSimpson_read_seqs <- mcols(BartSimpson_read_alignments)$seq

gc_content <- letterFrequency(BartSimpson_read_seqs, letters="GC", as.prob=TRUE)

hist(gc_content, breaks=50, col="pink", main="Read-level GC Content Bart Simpson", 
     xlab="GC Proportion (0.0 to 1.0")

# Lisa Simpson (Normopeso) # 
LisaSimpson_stats <- align(index="my_genome_index", 
      readfile1="LisaSimpson_R1.fastq.gz", 
      readfile2="LisaSimpson_R2.fastq.gz", output_file="LisaSimpson_Aligned.bam")

# Mapping rate
print(LisaSimpson_stats)
mapped_percent_LisaSimpson <- LisaSimpson_stats[2, 1] / LisaSimpson_stats[1, 1] * 100
print(paste("Mapping Rate:", round(mapped_percent_LisaSimpson, 2), "%"))

# Fragment Lenght
LisaSimpson_bam_file <- "LisaSimpson_Aligned.bam"

param <- ScanBamParam(what = c("isize"))
bam_data <- scanBam(LisaSimpson_bam_file, param = param)
frag_lengths <- abs(bam_data[[1]]$isize)

frag_lengths <- frag_lengths[frag_lengths > 50 & frag_lengths < 800]

summary(frag_lengths)
hist(frag_lengths, breaks=100, col="salmon",
     main="Fragment Length Distribution Lisa Simpson"
     , xlab="Size (bp")

# GC bias 
LisaSimpson_read_alignments <- readGAlignments("LisaSimpson_Aligned.bam", 
                                               param=ScanBamParam(what="seq"))
LisaSimpson_read_seqs <- mcols(LisaSimpson_read_alignments)$seq

gc_content <- letterFrequency(LisaSimpson_read_seqs, letters="GC", as.prob=TRUE)

hist(gc_content, breaks=50, col="pink", main="Read-level GC Content Lisa Simpson", 
     xlab="GC Proportion (0.0 to 1.0")

# Magggie Simpson (Normopeso)# 
MaggieSimpson_stats <- align(index="my_genome_index", 
      readfile1="MaggieSimpson_R1.fastq.gz", 
      readfile2="MaggieSimpson_R2.fastq.gz", output_file="MaggieSimpson_Aligned.bam")

# Mapping rate
print(MaggieSimpson_stats)
mapped_percent_MaggieSimpson <- MaggieSimpson_stats[2, 1] / MaggieSimpson_stats[1, 1] * 100
print(paste("Mapping Rate:", round(mapped_percent_MaggieSimpson, 2), "%"))

# Fragment Lenght
MaggieSimpson_bam_file <- "MaggieSimpson_Aligned.bam"

param <- ScanBamParam(what = c("isize"))
bam_data <- scanBam(MaggieSimpson_bam_file, param = param)
frag_lengths <- abs(bam_data[[1]]$isize)

frag_lengths <- frag_lengths[frag_lengths > 50 & frag_lengths < 800]

summary(frag_lengths)
hist(frag_lengths, breaks=100, col="salmon",
     main="Fragment Length Distribution Maggie Simpson"
     , xlab="Size (bp")

# GC bias 
MaggieSimpson_read_alignments <- readGAlignments("MaggieSimpson_Aligned.bam", 
                                               param=ScanBamParam(what="seq"))
MaggieSimpson_read_seqs <- mcols(MaggieSimpson_read_alignments)$seq

gc_content <- letterFrequency(MaggieSimpson_read_seqs, letters="GC", as.prob=TRUE)

hist(gc_content, breaks=50, col="pink", main="Read-level GC Content Maggie Simpson", 
     xlab="GC Proportion (0.0 to 1.0")
