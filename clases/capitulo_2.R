# DESeq2 y EdgeR son paqueterias para evaluar el analisis de expresion diferencial
# ambas herramientas utilizan un modelo negativo binomial para modelar
# las cuentas crudas

# Limma Voom es otro conjunto de herramientas utilizadas en el analisis de
# expresion diferencial, puede ser poco sensible en muestras pequeñas

# con la función vignette() podemos accdeder a la información del flujo
# de trabajo de DESeq2
vignette()
vignette("DESeq2")

# esto se debe hacer con las herramientas de trabajo, para analiza el flujo de trabajo

# los pasos para el analisis con DESeq2 son los siguientes
# contar cuentas asociadas con los genes
# normalizacion              control de calidad
# analisis de agrupacion no supervisado        control de calidad
# modelar cuentas crudas para cada gen      DE
# log2 fold changes                         DE
# evaluar expresion diferencial             DE

### evaluando cuentas crudas

# cargar data set
raw_counts_fibrosis <- readRDS("data/raw_counts_fibrosis.rds")

# visualizacion de data set
View(raw_counts_fibrosis)

# cargar metadata y visualizar
View(wt_metadata)

### IMPORTANTE ###
# para trabajar con DESeq2 el nombre de las muestras (nombre de las filas) en
# el metadata deben coincidir con el nombre de las muestras (nombre columnas)
# en el data set de las cuentas. Así como el orden, en ambos datasets deben
# tener el mismo orden de aparicion


# cambiando la primer columna de raw_counts_fibrosis como nombre de filas
raw_counts_fibrosis <-  column_to_rownames(raw_counts_fibrosis, var = "...1")

# cambiando los nombre de raw_counts_fibrosis para que sean iguales a wt_metadata

rownames(wt_metadata)

colnames(raw_counts_fibrosis) <- c("wt_fibrosis1","wt_fibrosis2","wt_normal1",
                                   "wt_normal2", "wt_fibrosis3","wt_normal3",
                                   "wt_fibrosis4")

# para saber si el orden de ambos data sets coincide se utiliza el siguiente argumento

all(rownames(wt_metadata) == colnames(raw_counts_fibrosis))

# cuando no coincide el orden se utiliza la función match()
match(vector1, vector2)
# vector1 orden deseado
# vector2 vector de valores a ordenar
# output indices de como ordenar el vector2 en el mismo orden que vector1

# esta funcion genera un vector con el orden en que deben acomodarse los datos
# del vector 2 para que coincidan con el vector 1

match(colnames(raw_counts_fibrosis), rownames(wt_metadata))

# con el output de la funcion se ordenan las filas del metadata

idx <- match(colnames(raw_counts_fibrosis), rownames(wt_metadata))

# despues se reordena acorde a los brackets []
wt_metadata <- wt_metadata[idx,]

# para el flujo de trabajo de DESeq2 se necesita crear un objeto DESeq
# para ello tulizamos la funcion DESeqDataSetFromMatrix()
DESeqDataSetFromMatrix()

dds_wt <- DESeqDataSetFromMatrix(countData = raw_counts_fibrosis,
                                 colData = wt_metadata,
                                 design= ~ condition)

# esta funcion genera un objeto Desq2 de la clase
# Ranged summarized experiment


# el primer paso es normalizar las cuentas crudas para el analisis de control
# de calidad
# los principales factores para normalizar son
# library depth, gene length, rna composition
# DESeq2 utiliza el metodo "median of ratios" para normalizacion
# este metodo ajusta el conteo de las cuentas crudas por tamaño de libreria y
#  es resistente a numeros largos de genes diferenciales
# para calcular las cuentas normalizadas se utiliza la funcion sobre
# el objeto dds
estimateSizeFactors()

dds_wt <- estimateSizeFactors(dds_wt)

sizeFactors(dds_wt)
