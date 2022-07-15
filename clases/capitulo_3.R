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


# explorando que tan adecuado es el modelo
# primero se deben calcular la media y la varianza de cada gen

mean_counts <- apply(raw_counts_fibrosis[,1:3], 1, mean)

variance_counts <- apply(raw_counts_fibrosis[,1:3], 1, var)

# despues se crea un df para graficar

df <- data.frame(mean_counts, variance_counts)

# se utiliza escala base 10 para graficar la media y varianzas

ggplot(df, aes(x = mean_counts, y = variance_counts)) +
  geom_point()+
  scale_y_log10()+
  scale_x_log10()+
  xlab("mean counts per gene") +
  ylab("variance per gene")

# la varianza aumenta con la media, esto se espera en análisis de RNA seq

# para calcular la dispersion de los datos con respecto a la media se
# utiliza la funcion plotDispEsts()

plotDispEsts(dds_wt)

# se espera que la dispersion disminuya al aumentar la media
