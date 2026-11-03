pkg <- c("mctest", "caret", "ggplot2", "sessioninfo", "glmnet",
         "leaps", "MASS", "car", "corrplot", "mlbench", "pls")
for (i in pkg) {
  if(!requireNamespace(i, quietly = T)){
    install.packages(i)
  } 
}
lapply(pkg, character.only = T, library)
set.seed(123)
sink("reporte_sesion.txt") #Creamos reporte_sesion.txt en nuestro directorio
sessioninfo::session_info()
sink()
file.show("reporte_sesion.txt") #Abre el reporte
# 4 -----------------------------------------------------------------------
data(BostonHousing, package = "mlbench")

# Filtrar únicamente medv y los 4 predictores requeridos
d_boston <- BostonHousing[, c("medv", "indus", "rm", "ptratio", "lstat")]

# Convertir medv a numérico si fuera necesario
d_boston$medv <- as.numeric(as.character(d_boston$medv))


# (a) Mejores subconjuntos con 1 a 4 predictores (Búsqueda Exhaustiva)
subconjuntos <- regsubsets(medv ~ ., data = d_boston, nvmax = 4)
resumen_sub <- summary(subconjuntos)

cat("\n--- (a) Mejores subconjuntos por número de predictores ---\n")
print(resumen_sub$outmat)

cat("\nSuma de Cuadrados de Residuos (RSS / SSE):\n")
print(resumen_sub$rss)

#lstat es el predictor que aparece seleccionada en la matriz de selección, luego rm es seleccionada en 3/4 modelos.
#También con la suma de cuadrados de residuos podemos ver que al agregar la segunda variable rm, el sse disminuye de 19742,38 a 15439.31.
#Al agregar una tercera variable que en este caso sería ptratio, baja a 13727.99 pero al agregar una cuarta variable sse solo baja 0.83 unidades que significaría que la variable indus no aporta al modelo


# (b) Modelo elegido según el Cp de Mallows
# El criterio busca minimizar el valor del Cp (aproximándose al número de parámetros p)
idx_cp <- which.min(resumen_sub$cp)

cat("\n--- (b) Modelo elegido según Cp de Mallows ---\n")
cat("Valores de Cp para modelos de 1 a 4 variables:\n")
print(resumen_sub$cp)

cat("\nNúmero óptimo de predictores según Cp:", idx_cp, "\n")
cat("Variables seleccionadas:\n")
print(names(which(resumen_sub$outmat[idx_cp, ] == "*")))

# (c) Modelo elegido según el criterio AIC
modelo_completo <- lm(medv ~ ., data = d_boston)

# stepAIC con k = 2 evalúa mediante el Criterio de Información de Akaike
modelo_aic <- stepAIC(modelo_completo, direction = "both", trace = TRUE)

cat("\n--- (c) Resumen del modelo óptimo por AIC ---\n")
summary(modelo_aic)

# (d) Comparación de los modelos con diagnósticos de residuos
# Reconstruir el modelo elegido por Cp para comparar sus supuestos
vars_cp <- names(which(resumen_sub$outmat[idx_cp, ] == "*"))
formula_cp <- as.formula(paste("medv ~", paste(vars_cp, collapse = " + ")))
modelo_cp <- lm(formula_cp, data = d_boston)

# Gráficos de diagnóstico para el modelo elegido por Cp
par(mfrow = c(2, 2))
plot(modelo_cp, main = "Diagnóstico: Modelo Cp")
par(mfrow = c(1, 1))

# Gráficos de diagnóstico para el modelo elegido por AIC
par(mfrow = c(2, 2))
plot(modelo_aic, main = "Diagnóstico: Modelo AIC")
par(mfrow = c(1, 1))

# 5 ---------------------------------------------------------------------

# 5.A ----------------------------------------------------------------------

# 1. Definir exactamente las 17 variables que indica la diapositiva 7
vars <- c("MPG.city", "Min.Price", "Price", "Max.Price", "EngineSize", "Horsepower",
          "RPM", "Rev.per.mile", "Fuel.tank.capacity", "Passengers", "Length",
          "Wheelbase", "Width", "Turn.circle", "Rear.seat.room", "Luggage.room", "Weight")

# 2. Filtrar el dataset
d <- Cars93[, vars]
d <- d[complete.cases(d), ] 

# 3. Ajustar el modelo inicial usando MPG.city como variable respuesta
modelo_inicial <- lm(MPG.city ~ ., data = d)

# 4. Ejecutar la selección paso a paso
modelo_optimo <- stepAIC(modelo_inicial, direction = "both", trace = TRUE)

# Ver el resultado del mejor subconjunto
summary(modelo_optimo)

# 5. Evaluar si el modelo es adecuado mediante los gráficos de diagnóstico
par(mfrow = c(2, 2))
plot(modelo_optimo)

#Resultados 
#El algoritmo eliminó 6 variables: Turn.circle, Width, Length, Horsepower, Rear.seat.room y Luggage.room.
#El proceso se detuvo en el paso final (AIC = 174.22) porque la mejor opción en la lista era <none> (no hacer nada). Quitar o agregar cualquier otra variable habría empeorado (subido) el AIC.
#En la tabla coefficients, quedan 10 variables y como afectan al rendimiento (mpg.city).
#Los coeficientes negativos como Weight tiene un coeficiente de -0.007, Por cada libra extra de peso, el auto rinde menos. }
#Si es positivo (como Wheelbase con 0.369), aumentar esa característica mejora el rendimiento, distancia entre ejes.
#Y según el Pr(>|t|) y las estrellas podemos ver que tan real es el efecto como el peso, la distancia entre los ejes, la capacidad del tanque, revoluciones por milla, etc
#Y como tenemos un R2 ajustado de 0.7637, los datos explican el 76.37% del modelo y un p-valor mínimo que indica que el modelo es significativo y válido
#Las variables de precio tienen comportamiento extraño porque Price queda positivo y  pero las demás quedan negativas.

# 5.B ---------------------------------------------------------------------
X_cars <- d[, -1]
R_cars <- cor(X_cars)

# Extraer pares con correlación mayor a 0.90

#Aislamos los pares con alta colinealidad > 0.9.
pares_altos <- which(abs(R_cars) > 0.90 & upper.tri(R_cars), arr.ind = TRUE)
pares_nombres <- data.frame(
  Var1 = rownames(R_cars)[pares_altos[, 1]],
  Var2 = colnames(R_cars)[pares_altos[, 2]],
  Correlacion = R_cars[pares_altos]
)

print("Pares con correlación > 0.90:")
print(pares_nombres[order(-abs(pares_nombres$Correlacion)), ])

# 2. Factor de Inflación de la Varianza (VIF) del modelo inicial
#calculamos el VIF, que mide cuanto se multiplica la varianza de cada Bi, debido a la correlacion con el resto de variables
#Valores de VIF superiores a 10, indican ERROR ESTANDAR esta tan inflado que el estadistico t pierde validez
#(las variables paracen no significativas)
vif_cars <- car::vif(modelo_inicial)
print("VIF de los predictores (ordenados):")
print(sort(round(vif_cars, 1), decreasing = TRUE))

# 3. Índices y Número de Condición
# Trabajamos con los autovalores de la matriz de correlación
ev_cars <- eigen(R_cars, symmetric = TRUE)$values
eta_cars <- sqrt(max(ev_cars) / ev_cars)

# El número de condición es el índice máximo
C_cars <- max(eta_cars)

print(paste("Número de condición (C):", round(C_cars, 1)))
print("Índices de condición superiores a 30:")
print(round(eta_cars[eta_cars > 30], 1))



# 5.C ---------------------------------------------------------------------

n <- nrow(d) #Número de registros (autos)
p <- length(coef(modelo_optimo)) # Número de parámetros (incluyendo intercepto)

#Umbral teórico para apalancamiento (mayor a esto refleja valores inusuales)
umbral_hii <- 2 * p / n

#Umbral teórico para influencia (establece el limite, si supera esto se esta sesgando al estimación completa del modelo)
umbral_cook <- 4 / n

#Estandarizamos para evaluar que autos tienen un MPG fuera de lo que el modelo predice para sus caracteristicas.
residuos_est <- rstudent(modelo_optimo)

# Boxplot de residuos para evaluar distribución y atípicos globales
boxplot(residuos_est, main = "Boxplot de Residuos Estudentizados", 
        ylab = "Residuos Estudentizados", col = "lightblue")
abline(h = c(-2, 2), col = "red", lty = 2) # Límites empíricos de atípicos

# 4. Apalancamiento / Leverages (Atípicos en X)
h_ii <- hatvalues(modelo_optimo)

# 5. Distancia de Cook (Influencia global)
dist_cook <- cooks.distance(modelo_optimo)

# Gráficos de diagnóstico con sus umbrales
par(mfrow = c(1, 2))

#Refleja que autos presentan caracteristicas inusuales respecto al resto de la muestra (puntos rojos)
plot(h_ii, type = "h", main = "Apalancamiento (hii)", 
     ylab = "hii", xlab = "Índice de observación")
abline(h = umbral_hii, col = "red", lty = 2)
puntos_leverage <- which(h_ii > umbral_hii)
points(puntos_leverage, h_ii[puntos_leverage], col = "red", pch = 19)


#Combina la informacion del residuo (que tan mal predicho fue el auto) con su apalancamientos (que tan inusuales son sus caracteristicas dada la muestra)
plot(dist_cook, type = "h", main = "Distancia de Cook", 
     ylab = "Cook's D", xlab = "Índice de observación")
abline(h = umbral_cook, col = "red", lty = 2)
puntos_cook <- which(dist_cook > umbral_cook)
points(puntos_cook, dist_cook[puntos_cook], col = "red", pch = 19)

par(mfrow = c(1, 1))

# 6. DFFITS y DFBETAS (Sensibilidad local y marginal)
df_fits <- dffits(modelo_optimo) #Mide cuantas desviaciones estandar cambia la prediccion de y especifica de un auto si se elimina de la base de datos de entrenamiento
df_betas <- dfbetas(modelo_optimo) #Mide cuanto varia y en que dirección el Bi de cada cariable al sacar una observacion especifica

# Identificar las observaciones que superan la distancia de Cook
print("Observaciones con influencia global crítica (Cook > 4/n):")
print(puntos_cook)

# Identificar las observaciones que superan el limite de apalancamiento
print("Observaciones con alto apalancamiento (hii > 2p/n):")
print(puntos_leverage)

# 7. Análisis de Sensibilidad: Reajuste sin las observaciones influyentes
datos_robustos <- d[-puntos_cook, ] #Eliminamos de los datos los autos que superan el umbral critico de cook
modelo_robusto <- lm(formula(modelo_optimo), data = datos_robustos) #Aplicamos el modelo sin esos datos

print("Comparación de coeficientes (Con influyentes vs Sin influyentes):")
comparacion_coef <- data.frame(
  Original = coef(modelo_optimo),
  Robusto = coef(modelo_robusto),
  Cambio_Absoluto = abs(coef(modelo_optimo) - coef(modelo_robusto))
)
print(round(comparacion_coef, 4))
#Veemos como varia cada Bi en el modelo con los autos que producen apalancamiento vs el modelo sin esas observaciones

# 6 -----------------------------------------------------------------------

vars <- c("MPG.city", "Min.Price", "Price", "Max.Price", "EngineSize", "Horsepower",
          "RPM", "Rev.per.mile", "Fuel.tank.capacity", "Passengers", "Length",
          "Wheelbase", "Width", "Turn.circle", "Rear.seat.room", "Luggage.room", "Weight")

d <- Cars93[, vars]
d <- d[complete.cases(d), ]
n <- nrow(d)

modelo_completo <- lm(MPG.city ~ ., data = d)

# 6.A ---------------------------------------------------------------------
# Selección paso a paso con BIC (k = log(n))
# Al fijar k = log(n), stepAIC pasa a optimizar el Criterio Bayesiano (BIC)
modelo_bic <- stepAIC(modelo_completo, direction = "both", k = log(n), trace = FALSE)

cat("\n--- (a) Modelo elegido por BIC ---\n")
summary(modelo_bic)

# 6.B ---------------------------------------------------------------------
# Atípicos, Apalancamiento e Influencia
leverages <- hatvalues(modelo_bic)
# Mide qué tan atípica es una observación en el espacio de las variables independientes X, sin considerar la variable respuesta Y.
cooks_d   <- cooks.distance(modelo_bic)
#Mide la influencia global que tiene cada observación sobre el ajuste del modelo. 
# Evalúa cuánto cambiarían las predicciones si se eliminara esa observación específica del conjunto de datos.
res_stud  <- rstudent(modelo_bic)
# Mide qué tan alejado está el valor observado del predicho en la variable respuest
p         <- length(coef(modelo_bic)) # número de parámetros

par(mfrow = c(2, 2))

# 1. Distancia de Cook (corte empírico 4/n)
plot(cooks_d, type = "h", main = "Distancia de Cook", ylab = "Distancia")
abline(h = 4/n, col = "red", lty = 2)

# 2. Apalancamiento / Leverage (corte 2p/n)
plot(leverages, type = "h", main = "Apalancamiento (Hat values)", ylab = "Leverage")
abline(h = 2 * p / n, col = "red", lty = 2)

# 3. Residuos Estudentizados (banda [-2, 2])
plot(res_stud, main = "Residuos Estudentizados", ylab = "RStudent")
abline(h = c(-2, 2), col = "red", lty = 2)

# 4. Gráfico combinado Residuals vs Leverage
plot(modelo_bic, which = 5)

par(mfrow = c(1, 1))


# 6.C ---------------------------------------------------------------------
# Todas las regresiones posibles con Cp de Mallows
busqueda_cp <- regsubsets(MPG.city ~ ., data = d, nvmax = 16)
sum_cp <- summary(busqueda_cp)

# Identificar el modelo que minimiza el Cp de Mallows
idx_mejor_cp <- which.min(sum_cp$cp)

cat("\n--- (c) Mejor modelo según Cp de Mallows ---\n")
cat("Número óptimo de predictores:", idx_mejor_cp, "\n")
print(sum_cp$outmat[idx_mejor_cp, ])


# 6.D ---------------------------------------------------------------------
# Comparación entre el modelo (c) Cp y el modelo (a) BIC
pred_bic <- names(coef(modelo_bic))[-1]
pred_cp  <- names(which(sum_cp$outmat[idx_mejor_cp, ] == "*"))

cat("\n--- (d) Comparación de Predictores ---\n")
cat("Predictores seleccionados por BIC (a):\n")
print(pred_bic)

cat("\nPredictores seleccionados por Cp de Mallows (c):\n")
print(pred_cp)

# 7. ----------------------------------------------------------------------
# Cargar los datos Hald (4 predictores: x1, x2, x3, x4 y respuesta: y)
Hald <- as.data.frame(Hald)

# 7.A Y 7.B ---------------------------------------------------------------
# Búsqueda Exhaustiva y Criterios de Selección
# regsubsets evalúa todos los modelos posibles combinando los 4 regresores
busqueda_exhaustiva <- regsubsets(y ~ ., data = Hald, nvmax = 4)
resumen <- summary(busqueda_exhaustiva)
cat("\n--- (a) Mejores ecuaciones por número de regresores (1 a 4) ---\n")
# Muestra qué variables se incluyen (marcadas con un asterisco *) en cada tamaño
print(resumen$outmat)
#El mejor modelo con una variable es solo con X4, el mejor modelo con dos variables es con x1 y x2 y así.
cat("\nSSE (Suma de Cuadrados del Error) mínima para cada tamaño de modelo:\n")
#La SSE va disminuyendo obligatoriamente a medida que se agregan variables
print(resumen$rss)
cat("\n--- (b) Mejor modelo según criterios de penalización ---\n")
# Se busca el máximo para R2 ajustado y el mínimo para Cp y BIC
mejor_r2adj <- which.max(resumen$adjr2)
mejor_cp <- which.min(resumen$cp)
mejor_bic <- which.min(resumen$bic)

cat("Mejor R2 ajustado (maximizar):", mejor_r2adj, "variables ->", resumen$adjr2[mejor_r2adj], "\n")
cat("Mejor Cp de Mallows (minimizar):", mejor_cp, "variables ->", resumen$cp[mejor_cp], "\n")
cat("Mejor BIC (minimizar):", mejor_bic, "variables ->", resumen$bic[mejor_bic], "\n")

#Según el R2 Ajustado el modelo debería tener dos variables, aunque según Cp y BIC el modelo debería contar con dos variables.
#En este caso vimos que el mejor modelo con dos variables sería con las variables X1 y X2

# Se definen los límites inferior (vacío) y superior (completo)
modelo_vacio <- lm(y ~ 1, data = Hald)
modelo_completo <- lm(y ~ ., data = Hald)

# 7.C ---------------------------------------------------------------------
# Forward
cat("\n--- (c) Selección Forward ---\n")
mod_forward <- step(modelo_vacio, direction = "forward", 
                    scope = list(lower = modelo_vacio, upper = modelo_completo), 
                    trace = 0) # Cambiar trace = 1 para ver el paso a paso en la bitácora
print(coef(mod_forward))

# 7.D ---------------------------------------------------------------------
# Backward
cat("\n--- (d) Selección Backward ---\n")
mod_backward <- step(modelo_completo, direction = "backward", trace = 0)
print(coef(mod_backward))

# 7.E ---------------------------------------------------------------------
# Stepwise
cat("\n--- (e) Selección Stepwise ---\n")
mod_stepwise <- step(modelo_completo, direction = "both", trace = 0)
print(coef(mod_stepwise))

#Podemos ver que los tres métodos de selección concuerdan con que deberíamos eliminar la variable X3 del modelo, y dejar las otras 3.
#Según los modelos automáticos y AIC deberíamos usar el modelo con 3 variables, X1, X2 y X3.

# 8 -----------------------------------------------------------------------

data("Hald", package = "mctest")
x_hald<-Hald[,2:5]
y_hald<-Hald[,1]
datos_hald <- data.frame(y = y_hald, x_hald)


# 8.A ---------------------------------------------------------------------

#Correlaciones entre las variables predictoras
R_hald<-cor(x_hald);R_hald

# 8.B ---------------------------------------------------------------------

#Ajustar los modelos elegidos en el Ejercicio 7(b)
#Modelo elegido por Cp de Mallows y BIC (2 variables)
modelo_bic <- lm(y ~ X1 + X2, data = datos_hald)
#Modelo elegido por R2 ajustado (3 variables)
modelo_r2adj <- lm(y ~ X1 + X2 + X4, data = datos_hald)

#Calcular los VIF
cat("\n--- (8b) VIF: Modelo de 2 variables (X1, X2) ---\n")
vif_bic <- vif(modelo_bic)
print(vif_bic)

cat("\n--- (8b) VIF: Modelo de 3 variables (X1, X2, X4) ---\n")
vif_r2adj <- vif(modelo_r2adj)
print(vif_r2adj)

#Ambos predictores presentan un VIF de 1.055 aprox, un valor prácticamente igual al mínimo teórico (1.0). 
#Esto demuestra ausencia de multicolinealidad entre X1 y X2, garantizando estimadores estables y de varianza mínima.

#Para el modelo de las 3 variables
#Mientras X1 se mantiene en 1.066, los factores de X2 (18.78) y X4 (18.94) superan el umbral crítico de 10. 
#Esto indica que la varianza de estos coeficientes está inflada casi 19 veces debido a la fuerte colinealidad existente entre X2 y X4.


# 8.C ---------------------------------------------------------------------

#Número de  condición

# Autovalores de la matriz de correlación, si es cercano a 0 la matriz pasa a ser no invertible, aumentando la varianza
#Justificado por la formula de sigma^2 * sum(1/lambda)
ev_hald <- eigen(R_hald, symmetric = TRUE)$values

eta_hald <- sqrt(max(ev_hald) / ev_hald)

# Número de condición, cuantifica la inestabilidad, si es mayor a 30, existe una multicolinealidad severa 
C_hald <- max(eta_hald)


print(signif(ev_hald, 3)) #Autovalores de la matriz

print(round(eta_hald, 2)) #Índices de condición

print(round(C_hald,2)) #Número de condición

# 9 -----------------------------------------------------------------------
X_hald <- scale(as.matrix(x_hald))

# Regla para muestras pequeñas (Diapositiva 52): con n = 13 se exige Leave-One-Out (LOOCV)
n_hald <- nrow(X_hald);n_hald
foldid_hald <- 1:n_hald;foldid_hald


# 9a ----------------------------------------------------------------------


# Rejilla logarítmica amplia para lambda asegurando que el mínimo quede dentro (EjemploR.R)
lam_grid_hald <- exp(seq(log(1e4), log(1e-4), length.out = 100))

# Ajuste del modelo Ridge
fit_ridge_hald <- glmnet(X_hald, y_hald, alpha = 0, standardize = FALSE, lambda = lam_grid_hald)

# Validación cruzada Leave-One-Out (LOOCV)
cv_ridge_hald <- cv.glmnet(X_hald, y_hald, alpha = 0, foldid = foldid_hald,
                           standardize = FALSE, lambda = lam_grid_hald, type.measure = "mse")

cat("EJERCICIO 9(a): RESULTADOS RIDGE (alpha = 0)\n")
cat("lambda.min (Mínimo MSE):", round(cv_ridge_hald$lambda.min, 4), "\n")
cat("lambda.1se (Regla 1 error estándar):", round(cv_ridge_hald$lambda.1se, 4), "\n")

# Gráficos de Ridge: Traza de coeficientes y curva CV
par(mfrow = c(1, 2))
plot(fit_ridge_hald, xvar = "lambda", label = TRUE, main = "Ridge: Trazas de Coeficientes")
plot(cv_ridge_hald, main = "Ridge: CV (LOOCV)")
par(mfrow = c(1, 1))

# Grados de libertad efectivos df(lambda) = sum(d_j^2 / (d_j^2 + lambda)) (Diapositiva 28)
ev_hald_n <- eigen(crossprod(X_hald) / n_hald, symmetric = TRUE)$values
df_ridge_calc <- function(l) sum(ev_hald_n / (ev_hald_n + l))
cat("Grados de libertad efectivos df(lambda.min):", round(df_ridge_calc(cv_ridge_hald$lambda.min), 2), "\n")
cat("Grados de libertad efectivos df(lambda.1se):", round(df_ridge_calc(cv_ridge_hald$lambda.1se), 2), "\n")



# 9b ----------------------------------------------------------------------


# Ajuste del modelo Lasso
fit_lasso_hald <- glmnet(X_hald, y_hald, alpha = 1, standardize = FALSE)

# Validación cruzada Leave-One-Out (LOOCV)
cv_lasso_hald <- cv.glmnet(X_hald, y_hald, alpha = 1, foldid = foldid_hald,
                           standardize = FALSE, type.measure = "mse")


cat("Ejercicio 9b: Resultados Lasso (alpha = 1)\n")

cat("lambda.min (Mínimo MSE):", round(cv_lasso_hald$lambda.min, 4), "\n")

cat("lambda.1se (Regla 1 error estándar):", round(cv_lasso_hald$lambda.1se, 4), "\n")

# Gráficos de Lasso: Traza de coeficientes y curva CV
par(mfrow = c(1, 2))
plot(fit_lasso_hald, xvar = "lambda", label = TRUE, main = "Lasso: Trazas de Coeficientes")
plot(cv_lasso_hald, main = "Lasso: CV (LOOCV)")
par(mfrow = c(1, 1))

# Verificación de Esparcidad (coeficientes anulados exactamente a cero en lambda.1se)
coef_lasso_1se <- as.matrix(coef(cv_lasso_hald, s = "lambda.1se"))
cat("Predictores seleccionados (coef != 0):", paste(names(which(coef_lasso_1se[-1, 1] != 0)), collapse = ", "), "\n")
cat("Predictores eliminados (coef == 0):", paste(names(which(coef_lasso_1se[-1, 1] == 0)), collapse = ", "), "\n")



# 9c ----------------------------------------------------------------------
# Ajuste de Elastic Net combinando L1 y L2 
fit_enet_hald <- glmnet(X_hald, y_hald, alpha = 0.5, standardize = FALSE)

# Validación cruzada Leave-One-Out (LOOCV)
cv_enet_hald <- cv.glmnet(X_hald, y_hald, alpha = 0.5, foldid = foldid_hald,
                          standardize = FALSE, type.measure = "mse")


cat("9c: Resultados elastic net (alpha = 0.5)\n")
cat("lambda.min (Mínimo MSE):", round(cv_enet_hald$lambda.min, 4), "\n")
cat("lambda.1se (Regla 1 error estándar):", round(cv_enet_hald$lambda.1se, 4), "\n")

# Gráficos elastic net
par(mfrow = c(1, 2))
plot(fit_enet_hald, xvar = "lambda", label = TRUE, main = "Elastic Net (alpha=0.5): Trazas")
plot(cv_enet_hald, main = "Elastic Net: CV (LOOCV)")
par(mfrow = c(1, 1))



# 9d ----------------------------------------------------------------------

#MCO sobre predictores estandarizados 
mod_ols_hald_std <- lm(y_hald ~ ., data = data.frame(y_hald = y_hald, X_hald))

# Función auxiliar para extraer coeficientes de cv.glmnet
bt <- function(cv_obj, s_corte) as.numeric(coef(cv_obj, s = s_corte))

# Tabla comparativa bajo el criterio parsimonioso lambda.1se
tabla_hald_comp <- cbind(
  MCO_OLS   = coef(mod_ols_hald_std),
  Ridge_1se = bt(cv_ridge_hald, "lambda.1se"),
  Lasso_1se = bt(cv_lasso_hald, "lambda.1se"),
  ENet_1se  = bt(cv_enet_hald, "lambda.1se")
)
rownames(tabla_hald_comp) <- c("(Intercepto)", "X1_Aluminato", "X2_SilicatoTri", 
                               "X3_Ferrito", "X4_SilicatoDi")


cat("9d: tabla comparativa Hald (lambda.1se)\n")
print(round(tabla_hald_comp, 3))

cat("\nPredictores activos retenidos (coef != 0):\n")
print(colSums(tabla_hald_comp[-1, ] != 0))



# extension de clase ------------------------------------------------------


# Predictores del modelo seleccionado por BIC en Fase 3
vars_bic_cars <- c("Price", "Rev.per.mile", "Fuel.tank.capacity", "Weight", "Min.Price", "Max.Price", "EngineSize", "RPM", "Wheelbase")

# Función de Error Cuadrático Medio de Predicción (RMSE)
rmse <- function(actual, pred) sqrt(mean((actual - pred)^2))

B <- 100  # 100 particiones independientes 80/20
mat_rmse <- matrix(NA, nrow = B, ncol = 5,
                   dimnames = list(NULL, c("OLS_Completo", "BIC_Fase3", "Ridge", "Lasso", "ElasticNet")))

for (b in 1:B) {
  # Partición 80% entrenamiento y 20% prueba
  tr <- sample(nrow(d), round(0.8 * nrow(d)))
  te <- setdiff(1:nrow(d), tr)
  
  # Reestandarización correcta usando estrictamente la media y sd del conjunto de entrenamiento
  mu_tr <- colMeans(d[tr, -1])
  sd_tr <- apply(d[tr, -1], 2, sd)
  
  X_tr <- scale(d[tr, -1], center = mu_tr, scale = sd_tr)
  X_te <- scale(d[te, -1], center = mu_tr, scale = sd_tr)
  
  y_tr <- d$MPG.city[tr]
  y_te <- d$MPG.city[te]
  
  # 1. MCO Completo
  m_ols <- lm(y_tr ~ X_tr)
  p_ols <- cbind(1, X_te) %*% coef(m_ols)
  mat_rmse[b, "OLS_Completo"] <- rmse(y_te, p_ols)
  
  # 2. Modelo BIC de Fase 3
  m_bic <- lm(MPG.city ~ ., data = d[tr, c("MPG.city", vars_bic_cars)])
  p_bic <- predict(m_bic, newdata = d[te, c("MPG.city", vars_bic_cars)])
  mat_rmse[b, "BIC_Fase3"] <- rmse(y_te, p_bic)
  
  # 3. Ridge (alpha = 0)
  cv_r <- cv.glmnet(X_tr, y_tr, alpha = 0, standardize = FALSE)
  p_r  <- predict(cv_r, newx = X_te, s = "lambda.1se")
  mat_rmse[b, "Ridge"] <- rmse(y_te, p_r)
  
  # 4. Lasso (alpha = 1)
  cv_l <- cv.glmnet(X_tr, y_tr, alpha = 1, standardize = FALSE)
  p_l  <- predict(cv_l, newx = X_te, s = "lambda.1se")
  mat_rmse[b, "Lasso"] <- rmse(y_te, p_l)
  
  # 5. Elastic Net (alpha = 0.5)
  cv_e <- cv.glmnet(X_tr, y_tr, alpha = 0.5, standardize = FALSE)
  p_e  <- predict(cv_e, newx = X_te, s = "lambda.1se")
  mat_rmse[b, "ElasticNet"] <- rmse(y_te, p_e)
}

# Tabla de RMSE promedio fuera de muestra
promedios_rmse <- colMeans(mat_rmse)
cat("\nRMSE Promedio Fuera de Muestra (100 particiones 80/20):\n")
print(round(promedios_rmse, 3))

# Gráfico Boxplot comparativo de RMSE fuera de muestra 
# Ajustar los márgenes de la figura: c(abajo, izquierda, arriba, derecha)
# El margen izquierdo (segundo valor) se aumenta a 6 para dar espacio suficiente.
par(mar = c(8, 6, 4, 2))

# Gráfico Boxplot con etiquetas visibles
boxplot(mat_rmse, 
        col = c("gray90", "#D69E2E", "#2B6CB0", "#319795", "#805AD5"),
        main = "Comparación de Error Fuera de Muestra",
        ylab = "", # Dejamos ylab vacío para colocarlo con title() con más espacio
        xlab = "Método de Regularización",
        las = 2, 
        cex.axis = 0.8)

# Agregar la etiqueta del eje Y ajustando su distancia con mgp o line
title(ylab = "RMSE (Millas por galón)", line = 4)

# Restaurar márgenes por defecto después del gráfico
par(mar = c(5, 4, 4, 2) + 0.1)
cat("\nModelo con menor error predictivo fuera de muestra:", names(which.min(promedios_rmse)), "\n")

# Desafio -----------------------------------------------------------------

# 10 ----------------------------------------------------------------------

mod_pcr_hald <- pcr(y ~ X1 + X2 + X3 + X4, data = datos_hald, scale = T, validation = "LOO")
summary(mod_pcr_hald)
validationplot(mod_pcr_hald, val.type = "MSEP", main = "PCR Hald: Error de predicción",
               legendpos = "topright")
#estimacion betas
ncomp_opt <- 2
coef_pcr_2comp <- coef(mod_pcr_hald, ncomp = ncomp_opt)
cat("\n--- Coeficientes estimados con", ncomp_opt, "Componentes Principales ---\n")
print(round(coef_pcr_2comp, 4))

#comparación formal entre MCO vs PCR
mod_mco_hald <- lm(y ~ X1 + X2 + X3 + X4, data = datos_hald)

comp_mco_pcr <- data.frame(
  MCO_Original = coef(mod_mco_hald)[-1],
  PCR_2Comp    = as.vector(coef_pcr_2comp)
)
rownames(comp_mco_pcr) <- c("X1", "X2", "X3", "X4")

cat("\n--- Comparación de Coeficientes (MCO vs PCR) ---\n")
print(round(comp_mco_pcr, 6))

# 11 ----------------------------------------------------------------------
url_deforest <- "https://raw.githubusercontent.com/GGavilan07/Modelos-Estadisticos-Activos/main/datos/deforest.csv"
deforest <- read.csv(url(url_deforest, method = "libcurl"))
colSums(is.na(deforest))

#transformacion de la variable respuesta
deforest$sqrtfcover <- sqrt(deforest$fcover)
vars_elimnar <- c("id", "ename", "code", "fcover")
datos_deforest <- deforest[,!colnames(deforest) %in% vars_elimnar]
train_prop <- 0.8
trainset <- sort(sample(1:nrow(datos_deforest), ceiling(nrow(datos_deforest) * train_prop)))

#conjuntos de entrenamiento 80% y 20%
train.set <- datos_deforest[trainset,]
test.set <- datos_deforest[-trainset,]

#pls
def_pls <- plsr(sqrtfcover ~ .,data = train.set, scale = TRUE, validation = "CV")

#seleccion de numero optimo de componentes mediante MSEP
msep_pls <- MSEP(def_pls)
ncomp_opt_pls <- which.min(msep_pls$val[1, 1, -1])

#predicción y cálculo del RMSE en test.set
pred_pls <- predict(def_pls, newdata = test.set, ncomp = ncomp_opt_pls)
rmse_pls <- sqrt(mean((test.set$sqrtfcover - pred_pls)^2))

cat("--- PLS (pls::plsr) ---\n")
cat("Componentes óptimos según CV:", ncomp_opt_pls, "\n")
cat("RMSE en conjunto de prueba:", round(rmse_pls, 4), "\n\n")

#ajuste del modelo PCR con la función pcr()
def_pcr <- pcr(sqrtfcover ~ ., data = train.set, scale = TRUE, validation = "CV")

#selección del número óptimo de componentes por MSEP
msep_pcr <- MSEP(def_pcr)
ncomp_opt_pcr <- which.min(msep_pcr$val[1, 1, -1])

#predicción y cálculo del RMSE en test.set
pred_pcr <- predict(def_pcr, newdata = test.set, ncomp = ncomp_opt_pcr)
rmse_pcr <- sqrt(mean((test.set$sqrtfcover - pred_pcr)^2))

cat("--- PCR (pls::pcr) ---\n")
cat("Componentes óptimos según CV:", ncomp_opt_pcr, "\n")
cat("RMSE en conjunto de prueba:", round(rmse_pcr, 4), "\n\n")

#tabla comparativa de rendimiento en el conjunto de prueba
tabla_comparativa_d3 <- data.frame(
  Metodo = c("PLS", "PCR"),
  Componentes = c(ncomp_opt_pls, ncomp_opt_pcr),
  RMSE_Test = c(rmse_pls, rmse_pcr)
)

cat("--- TABLA DE EVALUACIÓN EN TEST SET ---\n")
print(tabla_comparativa_d3)

#gráficos de error de validación cruzada (MSEP)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 4, 1.5))
validationplot(def_pls, val.type = "MSEP", main = "PLS: MSEP vs componentes")
validationplot(def_pcr, val.type = "MSEP", main = "PCR: MSEP vs componentes")
par(mfrow = c(1, 1))
# 12 ----------------------------------------------------------
url_github <- "https://raw.githubusercontent.com/GGavilan07/Modelos-Estadisticos-Activos/main/datos/glcamd.rds"
datos_glc_amd <- readRDS(url(url_github, method = "libcurl"))
X <- datos_glc_amd$x
Y <- as.factor(datos_glc_amd$y)
#lasso
cv.lasso <- cv.glmnet(X, Y, family = "binomial", alpha = 1)
lambda_opt_lasso <- cv.lasso$lambda.min
coef_lasso <- coef(cv.lasso, s = "lambda.min")
total_variables <- ncol(X)
coef_lasso_activos <- sum(coef_lasso !=0) -1
nulos_lasso <- total_variables - coef_lasso_activos
esparcidad_lasso <- nulos_lasso / total_variables
cat("Lambda óptimo:", round(lambda_opt_lasso, 6), "\n")
cat("Variables seleccionadas:", coef_lasso_activos, "de", total_variables, "\n")
cat("Grado de esparcidad:", esparcidad_lasso * 100, "%\n")
#enet
cv.enet_glauc <- cv.glmnet(X,Y, family = "binomial", alpha = 0.5) #No se específica alpha, usamos el estándar
lambda_opt_enet <- cv.enet_glauc$lambda.min
coef_enet <- coef(cv.enet_glauc, s = "lambda.min")
coef_enet_activos <- sum(coef_enet !=0) -1
nulos_enet <- total_variables - coef_enet_activos
esparcidad_enet <- nulos_enet / total_variables
cat("Lambda óptimo:", round(lambda_opt_enet, 6), "\n")
cat("Variables seleccionadas:", coef_enet_activos, "de", total_variables, "\n")
cat("Grado de esparcidad:", esparcidad_enet * 100, "%\n")
#Los coeficientes no nulos son el subconjunto de marcadores genéticos que tienen el mayor poder predictivo,
#para diferenciar ambas enfermedades.
#En este caso "AMD" es la clase 0, y "GLC" es la clase 1
#AMD es la degeneración macular asociada con la edad
#GLC grupo de afecciones oculares que dañan progresivamente el nervio óptico.
#Los betas > 0 Son los factores de riesgo para el glaucoma, la presencia de ese marcador genético aumenta
#la posibilidad de que el paciente tenga glaucoma en vez de AMD.
#Los betas < 0 son los factores que están en contra de la presencia del glaucoma, lo que sugiere que están
#relacionados con la presencia de AMD.
#Como en este caso aplicamos regresión logística al utilizar la familia binomial, podemos obtener la razón
#de posibilidades aplicando la función exponencial a nuestros betas
#Por ejemplo:
betas_lasso_act <- coef_lasso[coef_lasso[,1] !=0,]
tabla_razon_psb_lasso <- data.frame(
  Betas_lasso = betas_lasso_act,
  Razon_lasso = exp(betas_lasso_act)
)
tabla_razon_psb_lasso
#Si Razon_lasso = 1, la presencia del gen no altera las posibilidades de tener la enfermedad
#Si Razon_lasso > 1, significa que presentar ese gen multiplica las posibilidades de padecer glaucoma con respecto a los que no lo tienen
#Si Razon_lasso < 1, La presencia del marcador genético reduce las posibilidade de padecer glaucoma (actuando como factor protector), y se vincula   con la presencia de AMD.
#graficos
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 5, 1.5))
#error de cv lasso
plot(cv.lasso, main = "")
title("Lasso (alpha = 1): error CV", line = 3, cex.main = 1)
#error de cv elastic net
plot(cv.enet_glauc, main = "")
title("Elastic Net (alpha = 0.5): error CV", line = 3, cex.main = 1)
par(mfrow = c(1, 1))