# para realizar el análisis de expresion diferencial existen paqueterias como
# DESeq2, vamos a cargar la libreria

BiocManager::install("DESeq2")

# abrir paqueterias de trabajo
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(readr)

# pasos para realizar experimento de analisis de expresion diferencial
# preparacion de librerias/muestras biologicas > secuenciar lecturas >
# control de calidad > mapeo de genomas > conteo de lecturas asociado a los
# genes > analisis estadistico para identificar genes expresados diferencialmente

wtrawcounts <- read_csv("data/fibrosis_smoc2_rawcounts.csv")
saveRDS(wtrawcounts, file = "data/raw_counts_fibrosis.rds")
