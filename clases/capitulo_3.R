# comenzando análisis de expresion diferencial
# para esto se siguen tres paso
# 1. adecuar la cuentas crudas para cada gen para el modelo binomial negativo
# de DESeq2 y probar y examinar expresion diferencial tomando en cuenta
# los cambios log2 y extraer y visualizar la informacion

# para comenzar el nalisis se necesita el objeto DESeq dds_wt
# para hacer el ajuste del modelo se utiliza la funcion  DESeq()

dds_wt <- DESeq(dds_wt)

# este modelo contiene la informacion necesaria para realizar los
# analisis de expresion diferencial entre grupos especificos de muestras
