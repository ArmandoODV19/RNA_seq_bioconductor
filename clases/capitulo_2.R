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
# para calcular las cuentas normalizadas se utiliza la
# funcion estimateSizeFactors() sobre el objeto dds

estimateSizeFactors()

dds_wt <- estimateSizeFactors(dds_wt)

# DESeq2 utiliza sizeFactors() para normalizar las cuentas crudas
# las cuentas crudas de cada muestra se dividen por el tamaño de factor
# asociado especifico para cada muestra para normalizar
# para observar el tamaño de los factores usados para la normalizacion
# se utiliza la funcion sizeFactors()
sizeFactors(dds_wt)

# Una vez que se han calculado el tamaño de los factores y se han añadido al
# objeto DESeq2, se pueden extraer las cuentas normalizadas de el

# Para extraer las cuentas normalizadas se utiliza la funcion counts()

normalized_wt_counts <- counts(dds_wt, normalized = TRUE)
View(normalized_wt_counts)

# si normalized = FALSE, se extraen las cuentas crudas


# con las cuentas normalizadas se puede comparar el numero de
# cuentas entre muestras. Podemos explorar la similitud de las muestras
# con respecto a la expresión génica para evaluar la calidad del
# experimento. para esto se utilizan metodos de visualizacion
# para analisis de agrupamiento no supervisado

# previo a la visualizacion debemos transformar las cuentas normalizadass en log
# paramejorar la visualizacion en el agrupamiento

# Para RNA DESeq2 utiliza variance stabilizing transformation mediante
# la funcion vst(), la cual es una transformacion logritmica que modera
# la varianza a través de la media

vsd_wt <- vst(dds_wt, blind = TRUE)

# el clustering jerarquico con mapa de calor se utiliza para obtener
# la similitud de la expresion genica entre diferentes muestras del data set
# esta tecnica se utiliza para explorar que tan similares son las replicas entre
# cada una
# el heatmap se realiza a partir de un analisis de correlacion entre parejas
# del data set

# se espera que las replicas biologicas se agrupen entre si y las condiciones
# de la muestra aparte. Como la amayoria de los genes no muestran expresion
# diferencial, las muestras tienden a agruparse entre si. Muestras con
# valor de correlacion menor a 0.8 requieren investigacion sobre
# si las muestras tienen valores atipicos o estan contaminadas

# para el heatmap se debe extraer informacion del objeto vst y convertirlo
# en matriz utilizando la funcion assay()

vsd_mat_wt <- assay(vsd_wt)

vsd_cor_wt <- cor(vsd_mat_wt)
View(vsd_cor_wt)

# despues de extraer la informacion se utiliza la funcion pheatmap() para
# realizar el grafico

pheatmap(vsd_cor_wt, annotation = select(wt_metadata, condition))
# el argumento annotation = selecciona los factores en el metadata
# para incluirlos como antocaiones en la barras

### PCA

# PCS es una tecnica que enfatiza la variación de los datos
# para graficar los componentes principales se utiliza la
# funcion plotPCA() con el objeto vsd_wt y el argumento intgroup = "condition"
# el argumento ingroup se utiliza para seleccionar el argumento del
# metadata que se utilizara para colorear el plot


plotPCA(vsd_wt, intgroup = "condition")


