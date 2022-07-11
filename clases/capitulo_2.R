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
