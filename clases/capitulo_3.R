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

# 2. evaluacion fold2
# por default DESeq2 genera el test de Wald para parejas para evaluar
# diferencias en la expresion entre dos muestras para la condicion de interes

# el resultado se obtiene con la funcion result() especificando
# un valor de significancia con el argumento alpha =

results(dds_wt, alpha = 0.05)

# estos resultados comparan normla vs fibrosis
# para hacer comparaciones especificas entre condiciones se utilizan
# el argumento contrast

wt_res <- results(dds_wt,
                  contrast = c("condition", "fibrosis",
                               "normal"), alpha = 0.05)
# en contrast se coloca la condicion, el grupo problema y el grupo base o normal
# con este resultado podemos comparar fibrosis vs normal
# los nombres de la condicion de interes y
# los grupos muestra deben de coincidir con el metadata
# de esta forma los resultados muestran el log2 fold change comparando
# fibrosis vs normal

# para explorar mas los datos se pueden graficar

plotMA(wt_res, ylim = c(-8,8))
# este plot grafica la media de la cuentas normalizadas contra
# log2 fold changes, los genes expresados diferencialmente se
# encuentran en otro color

# log2 fold change shrinkage
# para mejorar el log2 fold change se realiza un shrinke (LFC shrinkage)
wt_res <- lfcShrink(dds_wt,
                    contrast = c("condition", "fibrosis",
                                 "normal"),
                    res = wt_res)

plotMA(wt_res, ylim = c(-8,8))
# este LFC suele ser mas precison


### explorando resultados de DESeq2

# para obtener descripciones de la columna de resultado se utiliza
# la funcion mcols()

mcols(wt_res)

head(wt_res, n = 10)
# para determinar los genes diferencialmente expresados
# se utiliza la columna padj que contiene los valores de p ajustados

# summary() muestra el total de genes expresados diferencialmente
# que se ajusta a el valor p
# este summary() es una funcion de DESeq2
summary(wt_res)

# Asimismo, se puede agregar un limite de log2 fold change para
# ser mas precisos en la seleccion de genes
# de esta forma se realiza una combinacion entre alpha y log2
# para seleccionar los genes DE

# en este ejemplo se utilizara un alpha de 0.05 y un log2 de 1.25
# convertido genera 0.32

wt_res <- results(dds_wt,
                  contrast = c("condition", "fibrosis", "normal"),
                  alpha = 0.05,
                  lfcThreshold = 0.32)

# y se realiza de nuevo la funcion lfcShrink() con los nuevos resultados

wt_res <- lfcShrink(dds_wt,
                    contrast = c("condition", "fibrosis", "normal"),
                    res=wt_res,
                    type="normal") # esta linea se agregó porque marcaba error, se debe usar "apeglm" o "ashr"



