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

# explorando data set de fribrosis

head(raw_counts_fibrosis)
str(raw_counts_fibrosis)
names(raw_counts_fibrosis)

# para determinar si un gen se encuentra diferencialmente expresado  en dos
# o mas grupos se utiliza
# el analisis estadistico de distribucion de conteo "count distribution"

# evaluar distribucion de cuentas

ggplot(raw_counts_fibrosis, aes(x=smoc2_normal1))+
  geom_histogram(stat = "bin", bins = 200)+
  xlab("Raw expression counts")+
  ylab("Number of genes")

# para iniciar el analisis de expresion diferencial se utiliza la funcion de
# DESeq2 DESeqDataSetFromMatrix() para tomar como datos de entrada
# el conteo de las cuentas crudas
DESeqDataSetFromMatrix()

dds <- DESeqDataSetFromMatrix(countData = raw_counts_fibrosis,
                              colData = metadata,
                              design ~ condition)

# se debe crear metadata y hacer un data frame con esa informacion

genotype <- c("wt", "wt", "wt", "wt", "wt", "wt", "wt")
condition <- c("fibrosis", "fibrosis", "normal", "normal",
               "fibrosis", "normal", "fibrosis")

wt_metadata <- data.frame(genotype, condition)

rownames(wt_metadata) <- c("wt_fibrosis1", "wt_fibrosis2", "wt_normal1", "wt_normal2",
                           "wt_fibrosis3", "wt_normal3", "wt_fibrosis4")
