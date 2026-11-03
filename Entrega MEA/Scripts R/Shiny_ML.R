# ==============================================================================
# Modelos lineales: regularización y diagnóstico de supuestos
# Modelos estadísticos activos | Análisis de bases Hald y Cars93
# Código autónomo basado en Avance_Modelos_Lineales.R (mismos modelos, semilla y parámetros)
# Framework: shiny + bslib (bootstrap 5) + plotly
# ==============================================================================

# ---------------------------------------------------------------------------- #
# 0. Dependencias y paquetes                                                    #
# ---------------------------------------------------------------------------- #

pkg <- c(
  "shiny", "bslib", "tidyr", "ggplot2", "plotly",
  "glmnet", "MASS", "mctest", "car", "leaps", "pls"
)

for (i in pkg) {
  if (!requireNamespace(i, quietly = TRUE)) {
    install.packages(i)
  }
}

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(tidyr)
  library(ggplot2)
  library(plotly)
  library(glmnet)
  library(MASS)
  library(mctest)
  library(car)
  library(leaps)
  library(pls)
})

# ---------------------------------------------------------------------------- #
# 1. Procedimientos estadísticos sobre Cars93 y Hald (secciones 5 a 10 del script)
# ---------------------------------------------------------------------------- #

set.seed(123)

# Paleta institucional (estilo UCM)
ucm <- list(
  azul       = "#0A2A5E",
  azul_medio = "#1F4E9C",
  azul_claro = "#DCE6F5",
  blanco     = "white",
  gris       = "#F4F6FA",
  gris_texto = "#5B6678",
  texto      = "#1E2A3A",
  rojo       = "#B3261E",
  verde      = "#2E7D32",
  turquesa   = "#00796B",
  violeta    = "#6A4C93",
  dorado     = "#D69E2E",
  dorado_osc = "#B7791F",
  dorado_cla = "#FBF3DC"
)

# Formateadores numéricos
f3 <- function(x) formatC(x, format = "f", digits = 3)
f4 <- function(x) formatC(x, format = "f", digits = 4)
rmse_calc <- function(act, prd) sqrt(mean((act - prd)^2))
extract_coef <- function(cv_obj, s_corte) {
  coefs <- coef(cv_obj, s = s_corte)
  z <- as.numeric(coefs)
  names(z) <- rownames(coefs)
  z
}
## ---- Base Cars93 (secciones 5 y 6) ----
vars_cars <- c(
  "MPG.city", "Min.Price", "Price", "Max.Price", "EngineSize", "Horsepower",
  "RPM", "Rev.per.mile", "Fuel.tank.capacity", "Passengers", "Length",
  "Wheelbase", "Width", "Turn.circle", "Rear.seat.room", "Luggage.room", "Weight"
)

d_cars <- MASS::Cars93[, vars_cars]
d_cars <- d_cars[complete.cases(d_cars), ]
n_cars <- nrow(d_cars)

# Modelo inicial completo
mod_inicial_cars <- lm(MPG.city ~ ., data = d_cars)

# Modelo según AIC (stepAIC both con k = 2)
mod_aic_cars <- stepAIC(mod_inicial_cars, direction = "both", trace = FALSE)

# Modelo según BIC (stepAIC both con k = log(n))
mod_bic_cars <- stepAIC(mod_inicial_cars, direction = "both", k = log(n_cars), trace = FALSE)
vars_bic_step_cars <- names(coef(mod_bic_cars))[-1]

# Subconjunto BIC
vars_bic_cars <- c("Price", "Rev.per.mile", "Fuel.tank.capacity", "Weight", "Min.Price", "Max.Price", "EngineSize", "RPM", "Wheelbase")
if (!setequal(vars_bic_cars, vars_bic_step_cars)) warning("vars_bic_cars no coincide con las variables del modelo BIC de stepAIC")

# Multicolinealidad
X_cars <- d_cars[, -1]
R_cars <- cor(X_cars)
pares_altos_idx <- which(abs(R_cars) > 0.90 & upper.tri(R_cars), arr.ind = TRUE)
pares_altos_cars <- data.frame(
  "Variable 1" = rownames(R_cars)[pares_altos_idx[, 1]],
  "Variable 2" = colnames(R_cars)[pares_altos_idx[, 2]],
  "Correlación" = R_cars[pares_altos_idx],
  check.names = FALSE, stringsAsFactors = FALSE
)
pares_altos_cars <- pares_altos_cars[order(-pares_altos_cars$`Correlación`), ]
pares_altos_cars$`Correlación` <- f4(pares_altos_cars$`Correlación`)

vif_inicial_cars <- car::vif(mod_inicial_cars)
vif_max_cars     <- max(vif_inicial_cars)
vif_max_nom_cars <- names(which.max(vif_inicial_cars))

ev_cars  <- eigen(R_cars, symmetric = TRUE)$values
eta_cars <- sqrt(max(ev_cars) / ev_cars)
C_cars   <- max(eta_cars)

# Influencia y atípicos
resumen_influencia <- function(m, n) {
  p  <- length(coef(m))
  ck <- cooks.distance(m)
  h  <- hatvalues(m)
  rs <- rstudent(m)
  list(
    p = p, umbral_hii = 2 * p / n, umbral_cook = 4 / n,
    cook = ck, hii = h, rstud = rs,
    idx_cook = which(ck > 4 / n), idx_hii = which(h > 2 * p / n),
    idx_rstud = which(abs(rs) > 2)
  )
}
infl_aic_cars <- resumen_influencia(mod_aic_cars, n_cars)
infl_bic_cars <- resumen_influencia(mod_bic_cars, n_cars)

# Análisis de sensibilidad
datos_robustos_cars <- if (length(infl_aic_cars$idx_cook) > 0) d_cars[-infl_aic_cars$idx_cook, ] else d_cars
mod_robusto_cars    <- lm(formula(mod_aic_cars), data = datos_robustos_cars)
comp_sensibilidad_cars <- data.frame(
  "Predictor"       = names(coef(mod_aic_cars)),
  "Modelo original" = f4(coef(mod_aic_cars)),
  "Sin influyentes" = f4(coef(mod_robusto_cars)),
  "Cambio absoluto" = f4(abs(coef(mod_aic_cars) - coef(mod_robusto_cars))),
  check.names = FALSE, stringsAsFactors = FALSE
)

# Búsqueda exhaustiva con Cp de Mallows
busq_cp_cars <- regsubsets(MPG.city ~ ., data = d_cars, nvmax = 16)
sum_cp_cars  <- summary(busq_cp_cars)
idx_cp_cars  <- which.min(sum_cp_cars$cp)
vars_cp_cars <- names(which(sum_cp_cars$outmat[idx_cp_cars, ] == "*"))
mod_cp_cars  <- lm(as.formula(paste("MPG.city ~", paste(vars_cp_cars, collapse = " + "))), data = d_cars)

comp_bic_cp_cars <- data.frame(
  "Variable" = names(d_cars)[-1],
  setNames(list(ifelse(names(d_cars)[-1] %in% vars_bic_step_cars, "Sí", "No")),
           paste0("En modelo BIC (", length(vars_bic_step_cars), " variables)")),
  setNames(list(ifelse(names(d_cars)[-1] %in% vars_cp_cars, "Sí", "No")),
           paste0("En modelo Cp (", length(vars_cp_cars), " variables)")),
  check.names = FALSE, stringsAsFactors = FALSE
)

# Tabla de coeficientes
tabla_coeficientes <- function(m) {
  tb <- as.data.frame(summary(m)$coefficients)
  data.frame(
    "Predictor"       = rownames(tb),
    "Estimación"      = f4(tb[[1]]),
    "Error estándar"  = f4(tb[[2]]),
    "Valor t"         = f3(tb[[3]]),
    "Pr(>|t|)"        = sprintf("%.4f", tb[[4]]),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}
tabla_coef_bic <- tabla_coeficientes(mod_bic_cars)

# Resumen comparativo de modelos de Cars93
resumen_modelo <- function(m, nombre, origen) {
  data.frame(
    "Modelo" = nombre, "Origen" = origen,
    "Predictores" = length(coef(m)) - 1,
    "R² ajustado" = sprintf("%.2f%%", 100 * summary(m)$adj.r.squared),
    "RSE (mpg)" = f3(summary(m)$sigma),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}
tabla_modelos_cars <- rbind(
  resumen_modelo(mod_inicial_cars, "Modelo inicial", "todas las covariables"),
  resumen_modelo(mod_aic_cars, "Modelo AIC", "stepAIC con k = 2"),
  resumen_modelo(mod_bic_cars, "Modelo BIC", "stepAIC con k = log(n)"),
  resumen_modelo(mod_cp_cars, "Modelo Cp de Mallows", "búsqueda exhaustiva")
)

# Comparación predictiva fuera de muestra
B_rep <- 100
mat_rmse_cars <- matrix(
  NA, nrow = B_rep, ncol = 5,
  dimnames = list(NULL, c("OLS_Completo", "BIC_Fase3", "Ridge", "Lasso", "ElasticNet"))
)

set.seed(123)
for (b in 1:B_rep) {
  tr <- sample(nrow(d_cars), round(0.8 * nrow(d_cars)))
  te <- setdiff(1:nrow(d_cars), tr)
  
  mu_tr <- colMeans(d_cars[tr, -1])
  sd_tr <- apply(d_cars[tr, -1], 2, sd)
  
  X_tr <- scale(d_cars[tr, -1], center = mu_tr, scale = sd_tr)
  X_te <- scale(d_cars[te, -1], center = mu_tr, scale = sd_tr)
  
  y_tr <- d_cars$MPG.city[tr]
  y_te <- d_cars$MPG.city[te]
  
  m_ols <- lm(y_tr ~ X_tr)
  mat_rmse_cars[b, "OLS_Completo"] <- rmse_calc(y_te, cbind(1, X_te) %*% coef(m_ols))
  
  m_bic <- lm(MPG.city ~ ., data = d_cars[tr, c("MPG.city", vars_bic_cars)])
  mat_rmse_cars[b, "BIC_Fase3"] <- rmse_calc(y_te, predict(m_bic, newdata = d_cars[te, c("MPG.city", vars_bic_cars)]))
  
  cv_r <- cv.glmnet(X_tr, y_tr, alpha = 0, standardize = FALSE)
  mat_rmse_cars[b, "Ridge"] <- rmse_calc(y_te, predict(cv_r, newx = X_te, s = "lambda.1se"))
  
  cv_l <- cv.glmnet(X_tr, y_tr, alpha = 1, standardize = FALSE)
  mat_rmse_cars[b, "Lasso"] <- rmse_calc(y_te, predict(cv_l, newx = X_te, s = "lambda.1se"))
  
  cv_e <- cv.glmnet(X_tr, y_tr, alpha = 0.5, standardize = FALSE)
  mat_rmse_cars[b, "ElasticNet"] <- rmse_calc(y_te, predict(cv_e, newx = X_te, s = "lambda.1se"))
}

prom_rmse_cars <- colMeans(mat_rmse_cars)
ganador_cars   <- names(which.min(prom_rmse_cars))
etiquetas_rmse <- c(
  OLS_Completo = "MCO completo", BIC_Fase3 = sprintf("Modelo BIC (%d variables)", length(vars_bic_cars)),
  Ridge = "Ridge", Lasso = "Lasso", ElasticNet = "Elastic Net"
)

# Trazas y CV globales de Cars93
set.seed(123)
X_cars_std <- scale(as.matrix(d_cars[, -1]))
cv_ridge_cars_all <- cv.glmnet(X_cars_std, d_cars$MPG.city, alpha = 0, standardize = FALSE)
cv_lasso_cars_all <- cv.glmnet(X_cars_std, d_cars$MPG.city, alpha = 1, standardize = FALSE)
cv_enet_cars_all  <- cv.glmnet(X_cars_std, d_cars$MPG.city, alpha = 0.5, standardize = FALSE)

## ---- Base Hald (secciones 7 a 10) ----
data("Hald", package = "mctest")
Hald <- as.data.frame(Hald)
x_hald <- Hald[, 2:5]
y_hald <- Hald[, 1]
datos_hald <- data.frame(y = y_hald, x_hald)
n_hald <- nrow(datos_hald)

# Búsqueda exhaustiva
busq_hald    <- regsubsets(y ~ ., data = datos_hald, nvmax = 4)
resumen_hald <- summary(busq_hald)
mejor_r2adj_hald <- which.max(resumen_hald$adjr2)
mejor_cp_hald    <- which.min(resumen_hald$cp)
mejor_bic_hald   <- which.min(resumen_hald$bic)

vars_por_tam_hald <- apply(resumen_hald$outmat, 1, function(r) paste(names(r)[r == "*"], collapse = ", "))
tabla_subconjuntos_hald <- data.frame(
  "Predictores"        = 1:4,
  "Variables elegidas" = unname(vars_por_tam_hald),
  "SSE"                = f4(resumen_hald$rss),
  "R² ajustado"        = f4(resumen_hald$adjr2),
  "Cp de Mallows"      = f4(resumen_hald$cp),
  "BIC"                = f4(resumen_hald$bic),
  check.names = FALSE, stringsAsFactors = FALSE
)

# Selección automática
modelo_vacio_hald    <- lm(y ~ 1, data = datos_hald)
mod_completo_hald    <- lm(y ~ ., data = datos_hald)
mod_forward_hald  <- step(modelo_vacio_hald, direction = "forward",
                          scope = list(lower = modelo_vacio_hald, upper = mod_completo_hald), trace = 0)
mod_backward_hald <- step(mod_completo_hald, direction = "backward", trace = 0)
mod_stepwise_hald <- step(mod_completo_hald, direction = "both", trace = 0)
tabla_auto_hald <- data.frame(
  "Método" = c("Forward", "Backward", "Stepwise"),
  "Variables retenidas" = c(
    paste(names(coef(mod_forward_hald))[-1], collapse = ", "),
    paste(names(coef(mod_backward_hald))[-1], collapse = ", "),
    paste(names(coef(mod_stepwise_hald))[-1], collapse = ", ")
  ),
  check.names = FALSE, stringsAsFactors = FALSE
)

# Modelos y VIF
mod_hald_2v <- lm(y ~ X1 + X2, data = datos_hald)
mod_hald_3v <- lm(y ~ X1 + X2 + X4, data = datos_hald)
vif_hald_2v <- car::vif(mod_hald_2v)
vif_hald_3v <- car::vif(mod_hald_3v)

# Correlaciones y número de condición
R_hald   <- cor(x_hald)
ev_hald  <- eigen(R_hald, symmetric = TRUE)$values
eta_hald <- sqrt(max(ev_hald) / ev_hald)
C_hald   <- max(eta_hald)
suma_masa_hald <- rowSums(x_hald)

# Regularización con LOOCV
X_hald_std    <- scale(as.matrix(x_hald))
foldid_hald   <- 1:n_hald
lam_grid_hald <- exp(seq(log(1e4), log(1e-4), length.out = 100))

fit_ridge_hald <- glmnet(X_hald_std, y_hald, alpha = 0, standardize = FALSE, lambda = lam_grid_hald)
cv_ridge_hald  <- cv.glmnet(X_hald_std, y_hald, alpha = 0, foldid = foldid_hald,
                            standardize = FALSE, lambda = lam_grid_hald, type.measure = "mse",
                            grouped = FALSE)

fit_lasso_hald <- glmnet(X_hald_std, y_hald, alpha = 1, standardize = FALSE)
cv_lasso_hald  <- cv.glmnet(X_hald_std, y_hald, alpha = 1, foldid = foldid_hald,
                            standardize = FALSE, type.measure = "mse", grouped = FALSE)

fit_enet_hald  <- glmnet(X_hald_std, y_hald, alpha = 0.5, standardize = FALSE)
cv_enet_hald   <- cv.glmnet(X_hald_std, y_hald, alpha = 0.5, foldid = foldid_hald,
                            standardize = FALSE, type.measure = "mse", grouped = FALSE)

# Grados de libertad efectivos de Ridge
ev_hald_n     <- eigen(crossprod(X_hald_std) / n_hald, symmetric = TRUE)$values
df_ridge_calc <- function(l) sum(ev_hald_n / (ev_hald_n + l))

# MCO estandarizado y tabla comparativa
mod_ols_hald_std <- lm(y_hald ~ ., data = data.frame(y_hald = y_hald, X_hald_std))
tabla_hald_num <- cbind(
  MCO   = coef(mod_ols_hald_std),
  Ridge = extract_coef(cv_ridge_hald, "lambda.1se"),
  Lasso = extract_coef(cv_lasso_hald, "lambda.1se"),
  ENet  = extract_coef(cv_enet_hald, "lambda.1se")
)
activos_hald <- colSums(tabla_hald_num[-1, ] != 0)
lasso_cero_hald <- names(which(tabla_hald_num[-1, "Lasso"] == 0))
tabla_hald_comp <- data.frame(
  "Predictor estandarizado" = c("(Intercepto)", "X1 (aluminato)", "X2 (silicato tricálcico)", "X3 (ferrito)", "X4 (silicato dicálcico)", "Predictores activos"),
  "MCO"             = c(f3(tabla_hald_num[, "MCO"]),   activos_hald[["MCO"]]),
  "Ridge (1se)"     = c(f3(tabla_hald_num[, "Ridge"]), activos_hald[["Ridge"]]),
  "Lasso (1se)"     = c(f3(tabla_hald_num[, "Lasso"]), activos_hald[["Lasso"]]),
  "Elastic Net (1se)" = c(f3(tabla_hald_num[, "ENet"]), activos_hald[["ENet"]]),
  check.names = FALSE, stringsAsFactors = FALSE
)

# Parámetros de regularización
tabla_lambda_hald <- data.frame(
  "Método"      = c("Ridge (alpha = 0)", "Lasso (alpha = 1)", "Elastic Net (alpha = 0.5)"),
  "lambda.min"  = f4(c(cv_ridge_hald$lambda.min, cv_lasso_hald$lambda.min, cv_enet_hald$lambda.min)),
  "lambda.1se"  = f4(c(cv_ridge_hald$lambda.1se, cv_lasso_hald$lambda.1se, cv_enet_hald$lambda.1se)),
  "gl efectivos (min)" = c(sprintf("%.2f", df_ridge_calc(cv_ridge_hald$lambda.min)), "no aplica", "no aplica"),
  "gl efectivos (1se)" = c(sprintf("%.2f", df_ridge_calc(cv_ridge_hald$lambda.1se)), "no aplica", "no aplica"),
  check.names = FALSE, stringsAsFactors = FALSE
)

# PCR
mod_pcr_hald   <- pcr(y ~ X1 + X2 + X3 + X4, data = datos_hald, scale = TRUE, validation = "LOO")
ncomp_opt      <- 2
coef_pcr_2comp <- coef(mod_pcr_hald, ncomp = ncomp_opt)
comp_mco_pcr <- data.frame(
  "Predictor"    = c("X1", "X2", "X3", "X4"),
  "MCO original" = f4(coef(mod_completo_hald)[-1]),
  "PCR con 2 componentes" = f4(as.vector(coef_pcr_2comp)),
  check.names = FALSE, stringsAsFactors = FALSE
)

# Diagnósticos homogéneos: VIF y C se calculan sobre los predictores realmente
# incluidos/activos; Cook y RMSE corresponden al ajuste auxiliar de esos predictores.
calc_diag_metrics <- function(m, rmse_override = NULL) {
  X <- model.matrix(m)[, -1, drop = FALSE]
  p <- ncol(X)
  if (p > 1) {
    vif_vals <- tryCatch(car::vif(m), error = function(e) NA_real_)
    vif_val <- if (all(is.na(vif_vals))) NA_real_ else max(vif_vals, na.rm = TRUE)
    R <- suppressWarnings(cor(X))
    ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    c_val <- if (min(ev) <= 1e-10) Inf else sqrt(max(ev) / min(ev))
  } else {
    vif_val <- if (p == 1) 1 else NA_real_
    c_val <- if (p == 1) 1 else NA_real_
  }
  c(
    VIF = vif_val,
    C = c_val,
    Cook = max(cooks.distance(m), na.rm = TRUE),
    RMSE = if (is.null(rmse_override)) sqrt(mean(resid(m)^2)) else rmse_override
  )
}

fmt_big <- function(x) formatC(x, format = "f", digits = 2, big.mark = " ")

# Tabla Cars93
m_cars_list <- list(mod_inicial_cars, mod_aic_cars, mod_bic_cars, mod_cp_cars)
nombres_m_cars <- c("Modelo inicial completo", "Modelo AIC", "Modelo BIC", "Modelo Cp de Mallows")
ind_cars_mat <- t(sapply(m_cars_list, calc_diag_metrics))
tabla_ind_cars <- data.frame(
  "Modelo" = nombres_m_cars,
  "VIF Máximo" = fmt_big(ind_cars_mat[, "VIF"]),
  "Número de Condición (C)" = f3(ind_cars_mat[, "C"]),
  "Cook Máximo" = f4(ind_cars_mat[, "Cook"]),
  "RMSE en muestra" = f3(ind_cars_mat[, "RMSE"]),
  check.names = FALSE, stringsAsFactors = FALSE
)

# Referencias para el simulador
mu_x_hald <- colMeans(x_hald)
sd_x_hald <- apply(x_hald, 2, sd)

# ---------------------------------------------------------------------------- #
# 2. Tema visual y helpers de interfaz                                          #
# ---------------------------------------------------------------------------- #

tema_app <- bs_theme(
  version      = 5,
  primary      = ucm$azul,
  secondary    = ucm$blanco,
  success      = ucm$verde,
  info         = ucm$azul_medio,
  warning      = ucm$blanco,
  danger       = ucm$rojo,
  bg           = "#FFFFFF",
  fg           = ucm$texto,
  base_font    = font_google("Open Sans"),
  heading_font = font_google("Montserrat"),
  code_font    = font_google("Fira Code"),
  font_scale   = 0.95
)

css_ucm <- paste0("
  body { background: ", ucm$gris, "; }
  .navbar { background: ", ucm$azul, " !important; border-bottom: 4px solid ", ucm$blanco, "; }
  .navbar .navbar-brand, .navbar .nav-link, .navbar .navbar-text { color: #FFFFFF !important; }
  .navbar .nav-link { font-weight: 600; letter-spacing: 0.01em; }
  .navbar .nav-link:hover { color: ", ucm$blanco, " !important; }
  .navbar .nav-link.active { color: ", ucm$blanco, " !important; border-bottom: 2px solid ", ucm$blanco, "; }
  .navbar-toggler { border-color: rgba(255,255,255,0.6); background: rgba(255,255,255,0.9); }
  .ucm-sello { background: ", ucm$blanco, "; color: ", ucm$azul, "; font-family: Montserrat, sans-serif;
               font-weight: 800; font-size: 1.05rem; padding: 4px 10px; border-radius: 4px; letter-spacing: 0.04em; }
  .ucm-titulo { display: block; color: #FFFFFF; font-family: Montserrat, sans-serif; font-weight: 700; font-size: 1.05rem; line-height: 1.1; }
  .ucm-subtitulo { display: block; color: ", ucm$azul_claro, "; font-size: 0.74rem; line-height: 1.2; }
  .cabecera-clara { background: #FFFFFF; color: ", ucm$azul, "; border-bottom: 1px solid ", ucm$azul_claro, "; border-left: 5px solid ", ucm$blanco, "; }
  .cabecera-azul { background: ", ucm$azul, "; color: #FFFFFF; border-left: 5px solid ", ucm$blanco, "; }
  .card { border-radius: 8px; overflow: hidden; }
  .texto-ucm { color: ", ucm$azul, "; }
  .tabla-ucm { font-size: 0.85rem; }
  .tabla-ucm thead th { background: ", ucm$azul, "; color: #FFFFFF; font-weight: 600; white-space: nowrap; border-color: ", ucm$azul, "; }
  .tabla-ucm thead th.col-dest { background: ", ucm$blanco, "; color: ", ucm$azul, "; }
  .tabla-ucm tbody td.col-dest { background: ", ucm$dorado_cla, "; font-weight: 600; }
  .tabla-ucm tr.fila-dest > td { background: ", ucm$dorado_cla, " !important; font-weight: 600; }
  .tabla-ucm tr.fila-dest > td:first-child { border-left: 4px solid ", ucm$blanco, "; }
  .badge-ucm { background: ", ucm$azul, "; color: #FFFFFF; font-weight: 600; font-size: 0.72rem; }
  .caja-destacada { background: ", ucm$azul_claro, "; border-left: 4px solid ", ucm$azul, "; border-radius: 6px; }
  .caja-ecuacion { background: #FFFFFF; border: 1px solid ", ucm$azul_claro, "; border-radius: 6px;
                   font-family: 'Fira Code', monospace; font-size: 0.82rem; line-height: 1.7; overflow-wrap: anywhere; word-break: normal; white-space: normal; color: ", ucm$azul, "; }
  .nav-pills .nav-link { color: ", ucm$azul, "; font-weight: 600; }
  .nav-pills .nav-link.active { background: ", ucm$azul, "; color: #FFFFFF; border-bottom: 3px solid ", ucm$blanco, "; }
  .form-check-input:checked { background-color: ", ucm$azul, "; border-color: ", ucm$azul, "; }
  .irs--shiny .irs-bar { background: ", ucm$azul_medio, "; border-color: ", ucm$azul_medio, "; }
  .irs--shiny .irs-single, .irs--shiny .irs-from, .irs--shiny .irs-to { background: ", ucm$azul, "; }
  .btn-outline-secondary { color: ", ucm$azul, "; border-color: ", ucm$azul, "; }
  .btn-outline-secondary:hover { background: ", ucm$azul, "; color: #FFFFFF; }
  .pie-ucm { font-size: 0.80rem; color: ", ucm$azul_claro, "; }
  .cols-2 { column-count: 2; column-gap: 1.5rem; }
  .cols-2 .form-check { break-inside: avoid; }
  @media (max-width: 768px) { .cols-2 { column-count: 1; } }
  .sim-inactivo { opacity: 0.35; pointer-events: none; }
  .sim-inactivo .control-label::after { content: ' (no interviene en este modelo)'; font-style: italic; color: ", ucm$gris_texto, "; }
")

niveles_tarjeta <- list(
  success = list(borde = ucm$verde,      fondo = "#EEF7EE"),
  warning = list(borde = ucm$dorado_osc, fondo = ucm$dorado_cla),
  danger  = list(borde = ucm$rojo,       fondo = "#FCEEEC"),
  info    = list(borde = ucm$azul_medio, fondo = ucm$azul_claro),
  primary = list(borde = ucm$azul,       fondo = "#EEF1F7")
)

tarjeta_metrica <- function(titulo, valor_texto, estado_umbral, descripcion, nivel = "danger") {
  est <- if (nivel %in% names(niveles_tarjeta)) niveles_tarjeta[[nivel]] else niveles_tarjeta[["primary"]]
  div(
    class = "card h-100 shadow-sm border-0",
    style = paste0("border-left: 5px solid ", est$borde, " !important; background: ", est$fondo, ";"),
    div(
      class = "card-body py-3 px-3",
      p(class = "fw-semibold mb-1 text-uppercase",
        style = paste0("font-size:0.75rem; color:", est$borde, "; letter-spacing:0.04em;"),
        titulo),
      h4(class = "fw-bold mb-0", style = paste0("color:", est$borde, "; font-size:1.45rem;"), valor_texto),
      span(class = "badge mt-1 mb-2",
           style = paste0("background:", est$borde, "; color:white; font-size:0.75rem;"),
           estado_umbral),
      p(class = "mb-0", style = paste0("font-size:0.80rem; line-height:1.25; color:", ucm$gris_texto, ";"), descripcion)
    )
  )
}

# Umbrales de interpretación
nivel_vif <- function(v) {
  if (is.na(v)) return(list(txt = "No aplica", nivel = "info"))
  if (v < 5) list(txt = "Bajo (umbral < 5)", nivel = "success")
  else if (v <= 10) list(txt = "Moderado (umbral 5 a 10)", nivel = "warning")
  else list(txt = "Alto (umbral > 10)", nivel = "danger")
}
nivel_cond <- function(cn) {
  if (is.na(cn)) return(list(txt = "No aplica", nivel = "info"))
  if (cn < 15) list(txt = "Bajo (umbral < 15)", nivel = "success")
  else if (cn <= 30) list(txt = "Moderado (umbral 15 a 30)", nivel = "warning")
  else list(txt = "Alto (umbral > 30)", nivel = "danger")
}

# Tabla HTML uniforme
tabla_ucm <- function(df, fila_resaltada = integer(0), insignia = rep("", nrow(df)), col_resaltadas = integer(0)) {
  encabezado <- tags$tr(lapply(seq_along(df), function(j) {
    tags$th(names(df)[j], class = if (j %in% col_resaltadas) "col-dest")
  }))
  filas <- lapply(seq_len(nrow(df)), function(i) {
    celdas <- lapply(seq_along(df), function(j) {
      contenido <- as.character(df[i, j])
      if (j == 1 && nzchar(insignia[i])) {
        contenido <- tagList(contenido, span(class = "badge badge-ucm ms-2", insignia[i]))
      }
      tags$td(contenido, class = if (j %in% col_resaltadas) "col-dest")
    })
    tags$tr(class = if (i %in% fila_resaltada) "fila-dest", celdas)
  })
  div(
    class = "table-responsive",
    tags$table(class = "table table-hover table-bordered mb-0 tabla-ucm",
               tags$thead(encabezado), tags$tbody(filas))
  )
}

panel_ucm <- function(titulo, ..., clase = "mb-4", cabecera = "clara") {
  card(
    class = paste("border-0 shadow-sm", clase),
    card_header(class = paste("fw-bold", if (cabecera == "azul") "cabecera-azul" else "cabecera-clara"), titulo),
    card_body(...)
  )
}

selector_base <- function(id, titulo, descripcion) {
  card(
    class = "mb-4 border-0 shadow-sm",
    card_body(
      class = "py-3",
      div(
        class = "d-flex flex-wrap align-items-center justify-content-between gap-3",
        div(
          h6(class = "fw-bold mb-0 texto-ucm", titulo),
          p(class = "mb-0", style = paste0("font-size:0.83rem; color:", ucm$gris_texto, ";"), descripcion)
        ),
        radioButtons(
          inputId = id, label = NULL,
          choices = c("Datos Cars93" = "cars", "Datos Hald" = "hald"),
          selected = "cars", inline = TRUE
        )
      )
    )
  )
}

def_cars <- list()
rango_cars <- list()
paso_cars  <- list()
banda_cars <- range(d_cars$MPG.city)
slider_cars <- function(id, etiqueta, variable, paso) {
  x <- d_cars[[variable]]
  minimo <- round(floor(round(min(x) / paso, 6)) * paso, 4)
  maximo <- round(ceiling(round(max(x) / paso, 6)) * paso, 4)
  valor  <- min(max(round(round(mean(x) / paso, 6)) * paso, minimo), maximo)
  def_cars[[id]] <<- valor
  rango_cars[[id]] <<- c(minimo, maximo)
  paso_cars[[id]] <<- paso
  sliderInput(id, etiqueta, min = minimo, max = maximo, value = valor, step = paso)
}
def_hald <- list()
slider_hald <- function(id, etiqueta, variable, paso) {
  x <- x_hald[[variable]]
  minimo <- round(floor(round(min(x) / paso, 6)) * paso, 4)
  maximo <- round(ceiling(round(max(x) / paso, 6)) * paso, 4)
  valor  <- min(max(round(round(mean(x) / paso, 6)) * paso, minimo), maximo)
  def_hald[[id]] <<- valor
  sliderInput(id, etiqueta, min = minimo, max = maximo, value = valor, step = paso)
}

tema_graf <- function(base = 10) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(), text = element_text(color = ucm$texto))
}
a_plotly <- function(p) {
  suppressWarnings(suppressMessages(ggplotly(p))) |>
    layout(paper_bgcolor = "transparent", plot_bgcolor = "transparent") |>
    config(displayModeBar = FALSE)
}

ecuacion_lm <- function(m, resp = "MPG.city") {
  b <- coef(m)
  signo <- ifelse(b[-1] < 0, "-", "+")
  mag <- sapply(abs(b[-1]), f4)
  paste0(resp, " = ", f4(b[1]),
         paste0(" ", signo, " ", mag, "*", names(b[-1]), collapse = ""))
}


# Ecuación uniforme para los modelos de Cars93, expresada en escala original.
coeficientes_cars <- function(clave) {
  if (clave %in% c("cars_bic", "cars_aic", "cars_ols", "cars_cp")) {
    m <- switch(clave, cars_bic = mod_bic_cars, cars_aic = mod_aic_cars,
                cars_ols = mod_inicial_cars, cars_cp = mod_cp_cars)
    return(coef(m))
  }
  cv_obj <- switch(clave, cars_ridge = cv_ridge_cars_all,
                   cars_lasso = cv_lasso_cars_all, cars_enet = cv_enet_cars_all)
  bb <- as.matrix(coef(cv_obj, s = "lambda.1se"))[, 1]
  names(bb) <- rownames(as.matrix(coef(cv_obj, s = "lambda.1se")))
  pendientes <- bb[-1] / esc_cars_std
  intercepto <- bb[1] - sum(bb[-1] * cen_cars_std / esc_cars_std)
  c("(Intercept)" = intercepto, pendientes)
}

ecuacion_cars_modelo <- function(clave) {
  b <- coeficientes_cars(clave)
  intercepto <- unname(b[1])
  bp <- b[-1]
  keep <- abs(bp) > 1e-8
  bp <- bp[keep]
  paste0("MPG.city = ", f4(intercepto),
         if (length(bp)) paste0(" ", ifelse(bp < 0, "-", "+"), " ",
           vapply(abs(bp), f4, character(1)), " · ", names(bp), collapse = "") else "")
}

etq_ganador <- unname(etiquetas_rmse[ganador_cars])
txt_otros_cars <- paste0(
  unname(etiquetas_rmse[setdiff(names(prom_rmse_cars), ganador_cars)]), " (",
  sprintf("%.3f", prom_rmse_cars[setdiff(names(prom_rmse_cars), ganador_cars)]), " mpg)",
  collapse = ", "
)
mod_ecuacion_cars <- mod_bic_cars
tit_ecuacion_cars <- "Ecuación del modelo BIC"

idx_par_hald <- which(upper.tri(R_hald), arr.ind = TRUE)
pares_hald <- data.frame(
  v1 = rownames(R_hald)[idx_par_hald[, 1]],
  v2 = colnames(R_hald)[idx_par_hald[, 2]],
  r  = R_hald[idx_par_hald],
  stringsAsFactors = FALSE
)
pares_hald <- pares_hald[order(-abs(pares_hald$r)), ]
idx_muy_alta_hald <- which(abs(pares_hald$r) > 0.90)
pares_hald_tab <- data.frame(
  "Variable 1" = pares_hald$v1,
  "Variable 2" = pares_hald$v2,
  "Correlación" = f4(pares_hald$r),
  "Lectura" = ifelse(abs(pares_hald$r) > 0.90, "muy alta: supera el umbral de 0.90",
                     ifelse(abs(pares_hald$r) > 0.80, "alta: cercana al umbral",
                            ifelse(abs(pares_hald$r) > 0.50, "moderada", "baja"))),
  check.names = FALSE, stringsAsFactors = FALSE
)

vif_c <- nivel_vif(vif_max_cars)
cond_c <- nivel_cond(C_cars)
vif2_h <- nivel_vif(max(vif_hald_2v))
vif3_h <- nivel_vif(max(vif_hald_3v))
cond_h <- nivel_cond(C_hald)

lasso_txt_hald <- if (length(lasso_cero_hald) > 0) {
  sprintf("%s = 0.000", paste(lasso_cero_hald, collapse = ", "))
} else {
  "Sin variables anuladas"
}
txt_hald_x3 <- sprintf(
  "En mínimos cuadrados con predictores estandarizados, X3 presenta una pendiente de %+.3f (p = %.2f en el modelo con predictores originales). Bajo Lasso con lambda.1se = %.3f, %s. Elastic Net con lambda.1se = %.3f retiene %d variables y asigna a X3 un coeficiente de %+.3f y a X4 de %+.3f.",
  tabla_hald_num["X3", "MCO"], summary(mod_completo_hald)$coefficients["X3", "Pr(>|t|)"],
  cv_lasso_hald$lambda.1se,
  if (length(lasso_cero_hald) > 0) {
    sprintf("%s queda anulado exactamente en cero y el submodelo resultante es {%s}",
            paste(lasso_cero_hald, collapse = ", "),
            paste(setdiff(rownames(tabla_hald_num)[-1], lasso_cero_hald), collapse = ", "))
  } else "ningún coeficiente se anula",
  cv_enet_hald$lambda.1se, activos_hald[["ENet"]],
  tabla_hald_num["X3", "ENet"], tabla_hald_num["X4", "ENet"]
)
txt_hald_just <- sprintf(
  "En la búsqueda exhaustiva, el modelo de %d variables {%s} minimiza Cp (%.2f) y BIC (%.2f), con VIF = %.3f. El modelo de %d variables {%s} maximiza el R² ajustado (%.4f), pero sus VIF para X2 y X4 son %.2f y %.2f. Los métodos regularizados penalizan esa inestabilidad sin forzar el descarte previo de componentes químicos.",
  mejor_cp_hald, vars_por_tam_hald[mejor_cp_hald], resumen_hald$cp[mejor_cp_hald], resumen_hald$bic[mejor_bic_hald],
  max(vif_hald_2v),
  mejor_r2adj_hald, vars_por_tam_hald[mejor_r2adj_hald], resumen_hald$adjr2[mejor_r2adj_hald],
  vif_hald_3v[["X2"]], vif_hald_3v[["X4"]]
)

sub_lasso_hald <- setdiff(rownames(tabla_hald_num)[-1], lasso_cero_hald)
txt_res_cars <- sprintf(
  "Para anticipar el rendimiento urbano de un auto que no se usó al ajustar, el método con menor error fue %s, con un error típico de %.2f millas por galón en %d particiones de prueba.",
  etq_ganador, min(prom_rmse_cars), B_rep
)
txt_res_cars_medido <- sprintf(
  "cuántas millas por galón rinde un auto en ciudad según %d características técnicas y de precio, sobre %d autos del mercado estadounidense.",
  ncol(d_cars) - 1, n_cars
)
txt_res_cars_alerta <- sprintf(
  "el precio mínimo, el promedio y el máximo miden casi lo mismo (VIF máximo de %s), por lo que el efecto individual de cada uno no es confiable aunque la predicción conjunta sí lo sea.",
  format(round(vif_max_cars, 0), big.mark = " ")
)
txt_res_cars_limite <- sprintf(
  "el modelo sirve para predecir, no para atribuir causas: sus valores p son optimistas porque las variables se eligieron con los mismos datos. Además, %d autos son influyentes y condicionan los coeficientes.",
  length(infl_bic_cars$idx_cook)
)
txt_res_hald <- sprintf(
  "Según el criterio experto, se adopta Lasso como modelo de 3 variables {%s}. Aunque la búsqueda por Cp/BIC favorece el submodelo de 2 variables, la decisión final prioriza conservar la tercera variable químicamente relevante; la regularización reduce la magnitud de los coeficientes.",
  paste(sub_lasso_hald, collapse = ", ")
)
txt_res_hald_medido <- sprintf(
  "cuánto calor libera el cemento a los 180 días según la proporción de %d compuestos químicos, con solo %d ensayos.",
  ncol(x_hald), n_hald
)
txt_res_hald_alerta <- sprintf(
  "los cuatro compuestos suman casi %.1f%% de la mezcla, así que subir uno obliga a bajar otros. Por eso X2 y X4 (correlación de %.3f) no se pueden separar con mínimos cuadrados.",
  mean(suma_masa_hald), R_hald["X2", "X4"]
)
txt_res_hald_limite <- sprintf(
  "con solo %d ensayos las conclusiones son orientativas y deben confirmarse con más experimentos antes de decidir sobre una formulación.",
  n_hald
)

medias_cars  <- colMeans(d_cars[, -1])
cen_cars_std <- attr(X_cars_std, "scaled:center")
esc_cars_std <- attr(X_cars_std, "scaled:scale")

clave_elegido_cars <- switch(ganador_cars,
                             OLS_Completo = "cars_ols", BIC_Fase3 = "cars_bic", Ridge = "cars_ridge",
                             Lasso = "cars_lasso", ElasticNet = "cars_enet"
)
etq_sim_cars <- c(
  cars_bic   = sprintf("Modelo BIC (%d variables)", length(coef(mod_bic_cars)) - 1),
  cars_aic   = sprintf("Modelo AIC (%d variables)", length(coef(mod_aic_cars)) - 1),
  cars_ols   = sprintf("Modelo inicial o MCO completo (%d variables)", length(coef(mod_inicial_cars)) - 1),
  cars_cp    = sprintf("Modelo Cp de Mallows (%d variables)", length(coef(mod_cp_cars)) - 1),
  cars_ridge = "Ridge (lambda.1se)",
  cars_lasso = "Lasso (lambda.1se)",
  cars_enet  = "Elastic Net (lambda.1se)"
)
grupos_sim_cars <- list(
  "Selección de variables" = c("cars_bic", "cars_aic", "cars_ols", "cars_cp"),
  "Regularización" = c("cars_ridge", "cars_lasso", "cars_enet")
)

mod_hald_1v <- lm(reformulate(strsplit(vars_por_tam_hald[[1]], ", ")[[1]], response = "y"), data = datos_hald)
clave_elegido_hald <- "hald_lasso"
vars_modelo_hald <- function(m) paste(names(coef(m))[-1], collapse = ", ")
etq_sim_hald <- c(
  hald_lasso = "Lasso (lambda.1se)",
  hald_enet  = "Elastic Net (lambda.1se)",
  hald_ridge = "Ridge (lambda.1se)",
  hald_1v    = sprintf("Modelo de 1 variable (%s)", vars_modelo_hald(mod_hald_1v)),
  hald_2v    = sprintf("Modelo de 2 variables (%s)", vars_modelo_hald(mod_hald_2v)),
  hald_3v    = sprintf("Modelo de 3 variables (%s)", vars_modelo_hald(mod_hald_3v)),
  hald_fwd   = sprintf("Forward (%s)", vars_modelo_hald(mod_forward_hald)),
  hald_bwd   = sprintf("Backward (%s)", vars_modelo_hald(mod_backward_hald)),
  hald_step  = sprintf("Stepwise (%s)", vars_modelo_hald(mod_stepwise_hald)),
  hald_ols   = "MCO completo (4 variables)",
  hald_pcr   = sprintf("Componentes principales (%d componentes)", ncomp_opt)
)
grupos_sim_hald <- list(
  "Regularización" = c("hald_lasso", "hald_enet", "hald_ridge"),
  "Búsqueda exhaustiva" = c("hald_1v", "hald_2v", "hald_3v"),
  "Selección automática" = c("hald_fwd", "hald_bwd", "hald_step"),
  "Modelo completo y componentes principales" = c("hald_ols", "hald_pcr")
)

fmt_clp <- function(x) paste0("$", formatC(round(x), format = "d", big.mark = "."))

nombres_cars <- paste(MASS::Cars93$Manufacturer, MASS::Cars93$Model)[as.integer(rownames(d_cars))]
opciones_auto_cars <- c("Elegir un auto" = "",
                        setNames(as.character(seq_len(nrow(d_cars))), sprintf("%s (%g mpg)", nombres_cars, d_cars$MPG.city)))
opciones_ensayo_hald <- c("Elegir un ensayo" = "",
                          setNames(as.character(seq_len(n_hald)), sprintf("Ensayo %d (%g cal/g)", seq_len(n_hald), y_hald)))

opciones_sim <- function(etq, grupos, elegido) {
  lapply(grupos, function(claves) {
    v <- claves
    names(v) <- ifelse(claves == elegido, paste0(etq[claves], " (modelo elegido)"), unname(etq[claves]))
    v
  })
}

pred_sim_cars <- function(clave, nd) {
  nd_full <- nd
  for (v in setdiff(names(medias_cars), names(nd))) nd_full[[v]] <- medias_cars[[v]]
  nd_full <- nd_full[, names(medias_cars), drop = FALSE]
  reg <- list(cars_ridge = cv_ridge_cars_all, cars_lasso = cv_lasso_cars_all, cars_enet = cv_enet_cars_all)
  if (clave %in% names(reg)) {
    x <- matrix((as.numeric(nd_full[1, ]) - cen_cars_std) / esc_cars_std, nrow = 1,
                dimnames = list(NULL, names(medias_cars)))
    return(as.numeric(predict(reg[[clave]], newx = x, s = "lambda.1se")))
  }
  m <- switch(clave, cars_bic = mod_bic_cars,
              cars_aic = mod_aic_cars, cars_ols = mod_inicial_cars, cars_cp = mod_cp_cars)
  as.numeric(predict(m, newdata = nd_full))
}

pred_sim_hald <- function(clave, nd_h) {
  nd_h <- as.data.frame(nd_h)
  if (clave %in% c("hald_lasso", "hald_enet", "hald_ridge")) {
    cv_obj <- switch(clave, hald_lasso = cv_lasso_hald, hald_enet = cv_enet_hald, hald_ridge = cv_ridge_hald)
    x_new_std <- (as.numeric(nd_h[1, ]) - mu_x_hald) / sd_x_hald
    b <- extract_coef(cv_obj, "lambda.1se")
    return(b[1] + sum(b[-1] * x_new_std))
  }
  if (clave == "hald_pcr") return(as.numeric(predict(mod_pcr_hald, newdata = nd_h, ncomp = ncomp_opt)))
  m <- switch(clave, hald_1v = mod_hald_1v, hald_2v = mod_hald_2v, hald_3v = mod_hald_3v,
              hald_fwd = mod_forward_hald, hald_bwd = mod_backward_hald, hald_step = mod_stepwise_hald,
              hald_ols = mod_completo_hald)
  as.numeric(predict(m, newdata = nd_h))
}

act_reg <- function(cv_obj, nombres) nombres[extract_coef(cv_obj, "lambda.1se")[-1] != 0]
ids_sim_cars <- c(
  Price = "sim_price", Min.Price = "sim_minprice", Max.Price = "sim_maxprice", EngineSize = "sim_engine",
  Horsepower = "sim_hp", RPM = "sim_rpm", Rev.per.mile = "sim_rev", Fuel.tank.capacity = "sim_tank",
  Passengers = "sim_pass", Wheelbase = "sim_wheel", Weight = "sim_weight",
  Length = "sim_length", Width = "sim_width", Turn.circle = "sim_turn",
  Rear.seat.room = "sim_rear", Luggage.room = "sim_luggage"
)
ids_sim_hald <- c(X1 = "sim_x1", X2 = "sim_x2", X3 = "sim_x3", X4 = "sim_x4")
vars_activas_cars <- list(
  cars_bic   = names(coef(mod_bic_cars))[-1],
  cars_aic   = names(coef(mod_aic_cars))[-1],
  cars_ols   = names(d_cars)[-1],
  cars_cp    = names(coef(mod_cp_cars))[-1],
  cars_ridge = act_reg(cv_ridge_cars_all, names(d_cars)[-1]),
  cars_lasso = act_reg(cv_lasso_cars_all, names(d_cars)[-1]),
  cars_enet  = act_reg(cv_enet_cars_all, names(d_cars)[-1])
)
vars_activas_hald <- list(
  hald_lasso = act_reg(cv_lasso_hald, names(x_hald)),
  hald_enet  = act_reg(cv_enet_hald, names(x_hald)),
  hald_ridge = act_reg(cv_ridge_hald, names(x_hald)),
  hald_1v    = names(coef(mod_hald_1v))[-1],
  hald_2v    = names(coef(mod_hald_2v))[-1],
  hald_3v    = names(coef(mod_hald_3v))[-1],
  hald_fwd   = names(coef(mod_forward_hald))[-1],
  hald_bwd   = names(coef(mod_backward_hald))[-1],
  hald_step  = names(coef(mod_stepwise_hald))[-1],
  hald_ols   = names(x_hald),
  hald_pcr   = names(x_hald)
)

coef_lineales_hald <- function(clave) {
  ceros <- data.frame(X1 = 0, X2 = 0, X3 = 0, X4 = 0)
  a <- pred_sim_hald(clave, ceros)
  b <- vapply(names(ceros), function(v) {
    nd <- ceros
    nd[[v]] <- 1
    pred_sim_hald(clave, nd) - a
  }, numeric(1))
  b[abs(b) < 1e-8] <- 0
  list(a = a, b = b)
}
ecuacion_hald_modelo <- function(clave) {
  cf <- coef_lineales_hald(clave)
  nz <- cf$b != 0
  paste0("y = ", f4(cf$a),
         paste0(" ", ifelse(cf$b[nz] < 0, "-", "+"), " ", sapply(abs(cf$b[nz]), f4), " · ", names(cf$b)[nz], collapse = ""))
}
opciones_ec_hald <- do.call(c, unname(opciones_sim(etq_sim_hald, grupos_sim_hald, clave_elegido_hald)))

# Tabla Hald: se incluyen todos los modelos disponibles en el selector.
modelos_diag_hald <- list(
  hald_1v = mod_hald_1v,
  hald_2v = mod_hald_2v,
  hald_3v = mod_hald_3v,
  hald_fwd = mod_forward_hald,
  hald_bwd = mod_backward_hald,
  hald_step = mod_stepwise_hald,
  hald_ols = mod_completo_hald
)
for (cl in c("hald_ridge", "hald_lasso", "hald_enet")) {
  cv_obj <- switch(cl, hald_ridge = cv_ridge_hald, hald_lasso = cv_lasso_hald, hald_enet = cv_enet_hald)
  active <- names(extract_coef(cv_obj, "lambda.1se")[-1])
  active <- active[extract_coef(cv_obj, "lambda.1se")[-1] != 0]
  if (!length(active)) active <- character(0)
  modelos_diag_hald[[cl]] <- if (length(active)) {
    lm(reformulate(active, response = "y"), data = datos_hald)
  } else lm(y ~ 1, data = datos_hald)
}
# PCR: diagnóstico sobre los componentes retenidos.
scores_pcr_diag <- as.data.frame(scores(mod_pcr_hald)[, seq_len(ncomp_opt), drop = FALSE])
names(scores_pcr_diag) <- paste0("PC", seq_len(ncomp_opt))
modelos_diag_hald$hald_pcr <- lm(y_hald ~ ., data = data.frame(y_hald = y_hald, scores_pcr_diag))

nombres_diag_hald <- c(
  hald_1v = sprintf("Modelo de 1 variable (%s)", paste(names(coef(mod_hald_1v))[-1], collapse = ", ")),
  hald_2v = "Modelo de 2 variables (X1, X2)",
  hald_3v = "Modelo de 3 variables (X1, X2, X4)",
  hald_fwd = "Forward",
  hald_bwd = "Backward",
  hald_step = "Stepwise",
  hald_ols = "MCO completo (4 variables)",
  hald_ridge = "Ridge (lambda.1se)",
  hald_lasso = "Lasso (lambda.1se) — modelo elegido",
  hald_enet = "Elastic Net (lambda.1se)",
  hald_pcr = sprintf("Componentes principales (%d componentes)", ncomp_opt)
)
rmse_hald_cv <- c(
  hald_ridge = sqrt(cv_ridge_hald$cvm[cv_ridge_hald$index["1se", ]]),
  hald_lasso = sqrt(cv_lasso_hald$cvm[cv_lasso_hald$index["1se", ]]),
  hald_enet = sqrt(cv_enet_hald$cvm[cv_enet_hald$index["1se", ]])
)
ind_hald_mat <- t(vapply(names(modelos_diag_hald), function(cl) {
  calc_diag_metrics(modelos_diag_hald[[cl]],
                    if (cl %in% names(rmse_hald_cv)) unname(rmse_hald_cv[[cl]]) else NULL)
}, numeric(4)))
tabla_ind_hald <- data.frame(
  "Modelo" = unname(nombres_diag_hald[rownames(ind_hald_mat)]),
  "VIF Máximo" = f3(ind_hald_mat[, "VIF"]),
  "Número de Condición (C)" = f3(ind_hald_mat[, "C"]),
  "Cook Máximo" = f4(ind_hald_mat[, "Cook"]),
  "RMSE en muestra / LOOCV" = f3(ind_hald_mat[, "RMSE"]),
  check.names = FALSE, stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------- #
# 3. Interfaz de usuario (UI)                                                   #
# ---------------------------------------------------------------------------- #

ui <- page_navbar(
  title = div(
    class = "d-flex align-items-center gap-2",
    span(class = "ucm-sello", "UCM"),
    div(
      span(class = "ucm-titulo", "Modelos lineales"),
      span(class = "ucm-subtitulo", "Universidad Católica del Maule")
    )
  ),
  window_title = "Modelos lineales | UCM",
  theme = tema_app,
  fillable = FALSE,
  header = tags$head(
    tags$style(HTML(css_ucm)),
    tags$script(HTML("
      (function() {
        function registrar() {
          Shiny.addCustomMessageHandler('sim_activos', function(m) {
            for (var i = 0; i < m.ids.length; i++) {
              $('#' + m.ids[i]).closest('.shiny-input-container').toggleClass('sim-inactivo', !m.activo[i]);
            }
          });
        }
        if (window.Shiny && Shiny.addCustomMessageHandler) { registrar(); }
        else { document.addEventListener('DOMContentLoaded', registrar); }
      })();
    "))
  ),
  
  # ========================================================================== #
  # Pestaña 1: contextualización de la base de datos                           #
  # ========================================================================== #
  nav_panel(
    title = "Contextualización de la base de datos",
    div(
      class = "container-fluid py-3",
      
      selector_base(
        "tab1_dataset",
        "Selección de base de datos: elija el conjunto de estudio",
        "Alterne entre el conjunto Cars93 y el conjunto Hald para revisar su contexto y observaciones."
      ),
      
      # Vista Cars93
      conditionalPanel(
        condition = "input.tab1_dataset == 'cars'",
        
        div(
          class = "row g-3 mb-4",
          div(class = "col-lg-6",
              panel_ucm(
                "Contexto analítico: rendimiento de combustible en automóviles", clase = "h-100",
                p(style = "font-size:0.87rem;",
                  sprintf("La base Cars93 contiene registros sobre especificaciones técnicas de vehículos del mercado estadounidense. Tras depurar casos incompletos, se cuenta con %d observaciones y %d covariables numéricas.",
                          n_cars, ncol(d_cars) - 1)),
                p(style = "font-size:0.85rem; line-height:1.6;",
                  strong("Variable de respuesta: "),
                  sprintf("rendimiento urbano en millas por galón (MPG.city). La distribución presenta una media de %.2f mpg con valores entre %g y %g mpg.",
                          mean(d_cars$MPG.city), min(d_cars$MPG.city), max(d_cars$MPG.city))),
                p(style = "font-size:0.85rem; line-height:1.6;",
                  strong("Problema de multicolinealidad: "),
                  "existe una relación de combinación lineal exacta en precios, donde Price ≈ (Min.Price + Max.Price) / 2. ",
                  "Asimismo, variables dimensionales y de motor como EngineSize y Weight presentan dependencias mecánicas que impiden la estimación estable por mínimos cuadrados.")
              )
          ),
          div(class = "col-lg-6",
              panel_ucm(
                "Correlaciones lineales altas: pares que superan el umbral 0.90", clase = "h-100",
                p(class = "mb-2", style = paste0("font-size:0.82rem; color:", ucm$gris_texto, ";"),
                  "Pares de covariables con valor absoluto de correlación bivariada superior a 0.90:"),
                uiOutput("tab1_tabla_pares_cars")
              )
          )
        ),
        
        accordion(
          id = "acc_resumen_cars", open = FALSE, class = "mb-3",
          accordion_panel(
            "Resumen para la decisión: recomendación y límites en Cars93",
            div(
              class = "p-3 rounded-2 mb-3 caja-destacada",
              h6(class = "fw-bold texto-ucm", "Recomendación"),
              p(style = "font-size:0.88rem; line-height:1.6; margin-bottom:0;", txt_res_cars)
            ),
            tags$ul(
              class = "ps-3 mb-0",
              style = "font-size:0.85rem; line-height:1.65;",
              tags$li(strong("Qué se midió: "), txt_res_cars_medido),
              tags$li(strong("Qué advertimos: "), txt_res_cars_alerta),
              tags$li(strong("Límites: "), txt_res_cars_limite)
            )
          )
        ),
        accordion(
          id = "acc_datos_cars", open = FALSE,
          accordion_panel(
            "Muestra de datos: primeros registros de Cars93",
            div(style = "overflow-x: auto;", uiOutput("tab1_tabla_datos_cars"))
          )
        )
      ),
      
      # Vista Hald
      conditionalPanel(
        condition = "input.tab1_dataset == 'hald'",
        
        div(
          class = "row g-3 mb-4",
          div(class = "col-lg-6",
              panel_ucm(
                "Contexto analítico: calor de hidratación en cemento Portland", clase = "h-100",
                p(style = "font-size:0.87rem;",
                  sprintf("La base de datos Hald contiene %d observaciones experimentales donde se mide el calor acumulado liberado a 180 días de fraguado (y, en calorías por gramo) ", n_hald),
                  "a partir del porcentaje en peso de 4 compuestos químicos activos (X1 a X4)."),
                p(style = "font-size:0.85rem; line-height:1.6;",
                  strong("Restricción física de masa cerrada: "),
                  sprintf("la suma de los cuatro compuestos representa casi la totalidad de la muestra (X1 + X2 + X3 + X4 ≈ %.1f%%).",
                          mean(suma_masa_hald))),
                p(style = "font-size:0.85rem; line-height:1.6;",
                  "Esta restricción implica que X4 quede casi determinado por los restantes componentes, generando correlaciones bivariadas ",
                  sprintf("(Cor(X2, X4) = %.3f y Cor(X1, X3) = %.3f) y distorsionando las pendientes en mínimos cuadrados.",
                          R_hald["X2", "X4"], R_hald["X1", "X3"]))
              )
          ),
          div(class = "col-lg-6",
              panel_ucm(
                "Definición de componentes químicos: significado de las variables", clase = "h-100",
                tags$ul(
                  class = "ps-3 mb-0",
                  style = "font-size:0.85rem; line-height:1.65;",
                  tags$li(strong("X1 (aluminato tricálcico): "), "aporte de calor rápido en el fraguado inicial."),
                  tags$li(strong("X2 (silicato tricálcico): "), "aporte de resistencia temprana con calor intermedio."),
                  tags$li(strong("X3 (ferrito tetracálcico): "), "bajo desprendimiento calórico."),
                  tags$li(strong("X4 (silicato dicálcico): "), "fraguado lento con muy bajo calor por gramo.")
                )
              )
          )
        ),
        
        div(
          class = "mb-4",
          panel_ucm(
            "Correlaciones lineales", clase = "mb-0",
            uiOutput("tab1_tabla_pares_hald")
          )
        ),
        
        accordion(
          id = "acc_resumen_hald", open = FALSE, class = "mb-3",
          accordion_panel(
            "Resumen para la decisión: recomendación y límites en Hald",
            div(
              class = "p-3 rounded-2 mb-3 caja-destacada",
              h6(class = "fw-bold texto-ucm", "Recomendación"),
              p(style = "font-size:0.88rem; line-height:1.6; margin-bottom:0;", txt_res_hald)
            ),
            tags$ul(
              class = "ps-3 mb-0",
              style = "font-size:0.85rem; line-height:1.65;",
              tags$li(strong("Qué se midió: "), txt_res_hald_medido),
              tags$li(strong("Qué advertimos: "), txt_res_hald_alerta),
              tags$li(strong("Límites: "), txt_res_hald_limite)
            )
          )
        ),
        accordion(
          id = "acc_datos_hald", open = FALSE,
          accordion_panel(
            "Registros completos: observaciones y balance de masa en Hald",
            div(style = "overflow-x: auto;", uiOutput("tab1_tabla_datos_hald"))
          )
        )
      )
    )
  ),
  
  # ========================================================================== #
  # Pestaña 2: métricas y selección de modelos                                 #
  # ========================================================================== #
  nav_panel(
    title = "Métricas y selección de modelos",
    div(
      class = "container-fluid py-3",
      
      selector_base(
        "tab2_dataset",
        "Selección de base de datos: evaluación de modelos",
        "Consulte el modelo elegido, las métricas reactivas por modelo y las comparaciones técnicas."
      ),
      
      # VISTA CARS93
      conditionalPanel(
        condition = "input.tab2_dataset == 'cars'",
        
        panel_ucm(
          paste0(etq_ganador, ", modelo con menor error fuera de muestra"), cabecera = "azul",
          div(
            class = "p-3 rounded-2 mb-3 caja-destacada",
            uiOutput("tab2_titulo_ecuacion_cars"),
            uiOutput("tab2_ecuacion_cars")
          ),
          h6(class = "fw-bold texto-ucm", "Justificación de la elección según las métricas"),
          p(style = "font-size:0.87rem; line-height:1.6;",
            sprintf("En la evaluación de %d particiones independientes 80/20, ", B_rep),
            sprintf("%s alcanza el menor error cuadrático medio fuera de muestra, con ", etq_ganador),
            strong(sprintf("%.3f mpg", min(prom_rmse_cars))),
            sprintf(", frente a %s. ", txt_otros_cars),
            sprintf("El modelo BIC conserva %d predictores (%s) y el modelo AIC conserva %d.",
                    length(vars_bic_step_cars), paste(vars_bic_step_cars, collapse = ", "),
                    length(coef(mod_aic_cars)) - 1))
        ),
        
        panel_ucm(
          "Selección del modelo para evaluación de indicadores", cabecera = "clara", clase = "mb-3",
          p(class = "mb-2", style = paste0("font-size:0.83rem; color:", ucm$gris_texto, ";"),
            "Seleccione el modelo para recalcular dinámicamente los 4 cuadros de indicadores de diagnóstico:"),
          radioButtons("tab2_cars_modelo_sel", label = NULL, width = "100%",
                       choices = c(
                         "Modelo BIC (9 variables - Elegido)" = "cars_bic",
                         "Modelo AIC (10 variables)" = "cars_aic",
                         "Modelo inicial MCO completo (16 variables)" = "cars_ols",
                         "Modelo Cp de Mallows (10 variables)" = "cars_cp",
                         "Ridge (lambda.1se)" = "cars_ridge",
                         "Lasso (lambda.1se)" = "cars_lasso",
                         "Elastic Net (lambda.1se)" = "cars_enet"
                       ),
                       selected = "cars_bic", inline = TRUE)
        ),
        
        uiOutput("tab2_kpi_cards_cars"),
        
        accordion(
          id = "acc_ind_cars", open = FALSE, class = "mb-4",
          accordion_panel(
            "Tabla comparativa: indicadores de diagnóstico por modelo (Cars93)",
            p(class = "mb-2", style = paste0("font-size:0.82rem; color:", ucm$gris_texto, ";"),
              "Comparativa de VIF, número de Condición (C), distancia de Cook máxima y error RMSE para cada especificación:"),
            uiOutput("tab2_tabla_indicadores_cars")
          )
        ),
        
        panel_ucm(
          "Comparación de modelos de Cars93 según criterio de selección",
          p(class = "mb-2", style = paste0("font-size:0.82rem; color:", ucm$gris_texto, ";"),
            "R² ajustado y error residual (RSE) de cada modelo ajustado sobre la muestra completa:"),
          uiOutput("tab2_tabla_modelos_cars")
        ),
        panel_ucm(
          "Comparación de métodos: error cuadrático medio fuera de muestra", clase = "mb-4",
          p(class = "mb-2", style = paste0("font-size:0.82rem; color:", ucm$gris_texto, ";"),
            sprintf("RMSE promedio evaluado sobre %d particiones independientes, 80%% entrenamiento y 20%% prueba:", B_rep)),
          uiOutput("tab2_tabla_rmse_cars_html")
        ),
        accordion(
          id = "acc_tecnico_cars", open = FALSE,
          accordion_panel(
            "Detalle técnico para consultas: umbrales, coeficientes y comparación con Cp",
            panel_ucm(
              "Reglas de decisión y umbrales de interpretación",
              div(
                class = "row g-3",
                div(class = "col-md-4",
                    h6(class = "fw-bold texto-ucm", "Factor de inflación de varianza (VIF)"),
                    p(style = "font-size:0.83rem;",
                      "Mide el incremento de varianza de un estimador debido a la colinealidad con el resto de covariables:"),
                    tags$ul(style = "font-size:0.80rem; padding-left:1.1rem;",
                            tags$li(strong("Menor a 5: "), "bajo (estimadores estables)."),
                            tags$li(strong("Entre 5 y 10: "), "moderado (precaución en la inferencia)."),
                            tags$li(strong("Mayor a 10: "), "alto (error estándar inflado, pruebas t no válidas)."))
                ),
                div(class = "col-md-4",
                    h6(class = "fw-bold texto-ucm", "Número de condición (C)"),
                    p(style = "font-size:0.83rem;",
                      "Cuantifica la cercanía de la matriz de correlación a ser singular o no invertible:"),
                    tags$ul(style = "font-size:0.80rem; padding-left:1.1rem;",
                            tags$li(strong("Menor a 15: "), "bajo (matriz bien condicionada)."),
                            tags$li(strong("Entre 15 y 30: "), "moderado (colinealidad emergente)."),
                            tags$li(strong("Mayor a 30: "), "alto (colinealidad severa e inestabilidad)."))
                ),
                div(class = "col-md-4",
                    h6(class = "fw-bold texto-ucm", "Medidas de influencia (Cook y apalancamiento)"),
                    p(style = "font-size:0.83rem;",
                      "Permiten detectar casos que alteran los parámetros de regresión:"),
                    tags$ul(style = "font-size:0.80rem; padding-left:1.1rem;",
                            tags$li(strong("Cook mayor a 4/n: "), "alto (observación influyente sobre el ajuste global)."),
                            tags$li(strong("Apalancamiento hii mayor a 2p/n: "), "alto (valor atípico en las covariables X)."),
                            tags$li(strong("|Rstudent| mayor a 2: "), "alto (residuo atípico en la variable Y)."))
                )
              )
            ),
            panel_ucm(
              paste0("Tabla de coeficientes: modelo BIC de ", length(vars_bic_step_cars), " predictores"),
              clase = "mb-4",
              div(style = "overflow-x: auto;", uiOutput("tab2_tabla_coef_bic"))
            ),
            panel_ucm(
              "Comparación de variables seleccionadas: modelo BIC frente a Cp de Mallows",
              p(class = "mb-2", style = paste0("font-size:0.82rem; color:", ucm$gris_texto, ";"),
                sprintf("BIC penaliza con k = log(n) = %.2f y conserva %d variables. Cp de Mallows penaliza con menor intensidad y conserva %d:",
                        log(n_cars), length(vars_bic_step_cars), length(vars_cp_cars))),
              uiOutput("tab2_tabla_comp_bic_cp")
            )
          )
        )
      ),
      
      # VISTA HALD
      conditionalPanel(
        condition = "input.tab2_dataset == 'hald'",
        
        panel_ucm(
          "Resultados clave y justificación: regularización frente a selección de subconjuntos", cabecera = "azul",
          div(
            class = "p-3 rounded-2 mb-3 caja-destacada",
            uiOutput("tab2_titulo_ecuacion_hald"),
            uiOutput("tab2_ecuacion_hald")
          ),
          div(
            class = "p-3 rounded-2 mb-3 caja-destacada",
            h6(class = "fw-bold texto-ucm", "Comportamiento de X3 bajo regularización"),
            p(style = "font-size:0.86rem; margin-bottom:0;", txt_hald_x3)
          ),
          h6(class = "fw-bold texto-ucm", "Justificación de la elección según las métricas"),
          p(style = "font-size:0.87rem; line-height:1.6;",
            "La búsqueda por Cp/BIC favorece el submodelo de 2 variables desde el punto de vista estadístico; sin embargo, según el criterio experto se conserva el modelo Lasso de 3 variables. ",
            txt_hald_just)
        ),

        panel_ucm(
          "Selección del modelo para evaluación de indicadores", cabecera = "clara", clase = "mb-3",
          p(class = "mb-2", style = paste0("font-size:0.83rem; color:", ucm$gris_texto, ";"),
            "Seleccione el modelo para actualizar su ecuación y recalcular dinámicamente los 4 cuadros de indicadores:"),
          radioButtons("tab2_hald_modelo_ec", label = NULL, width = "100%",
                       choices = opciones_ec_hald, selected = clave_elegido_hald, inline = TRUE)
        ),

        uiOutput("tab2_kpi_cards_hald"),
        
        accordion(
          id = "acc_ind_hald", open = FALSE, class = "mb-4",
          accordion_panel(
            "Tabla comparativa: Indicadores de diagnóstico por modelo (Hald)",
            p(class = "mb-2", style = paste0("font-size:0.82rem; color:", ucm$gris_texto, ";"),
              "Comparativa de VIF, Número de Condición (C), Distancia de Cook máxima y Error RMSE para cada submodelo:"),
            uiOutput("tab2_tabla_indicadores_hald")
          )
        ),
        
        panel_ucm(
          "Comparación de coeficientes estandarizados: contraste entre MCO y las penalizaciones Ridge, Lasso y Elastic Net",
          p(class = "mb-2", style = paste0("font-size:0.82rem; color:", ucm$gris_texto, ";"),
            "Estimaciones bajo la regla de 1 error estándar (lambda.1se):"),
          uiOutput("tab2_tabla_comp_hald_html")
        ),
        accordion(
          id = "acc_tecnico_hald", open = FALSE,
          accordion_panel(
            "Detalle técnico para consultas: criterios, búsqueda exhaustiva, selección automática y PCR",
            panel_ucm(
              "Reglas de decisión y umbrales de interpretación: criterios en Hald",
              div(
                class = "row g-3",
                div(class = "col-md-4",
                    h6(class = "fw-bold texto-ucm", "Criterios de penalización en subconjuntos"),
                    p(style = "font-size:0.83rem;",
                      "R² ajustado busca maximizar capacidad explicativa penalizando grados de libertad. Cp de Mallows busca un valor cercano a p. BIC impone mayor parsimonia al penalizar con log(n).")
                ),
                div(class = "col-md-4",
                    h6(class = "fw-bold texto-ucm", "Umbrales de inflación de varianza"),
                    p(style = "font-size:0.83rem;",
                      sprintf("El modelo de 2 variables omite X4. Al incorporarla, el VIF máximo sube a %.2f (%s), de modo que la colinealidad no desaparece sino que queda condicionada a las variables elegidas.",
                              max(vif_hald_3v), tolower(sub(" \\(.*", "", vif3_h$txt))))
                ),
                div(class = "col-md-4",
                    h6(class = "fw-bold texto-ucm", "Selección de lambda (min frente a 1se)"),
                    p(style = "font-size:0.83rem;",
                      "lambda.min minimiza el error muestral de validación cruzada. lambda.1se applies la regla de 1 error estándar para seleccionar el modelo más parsimonioso sin pérdida significativa de precisión.")
                )
              )
            ),
            panel_ucm(
              "Búsqueda exhaustiva: evaluación de 1 a 4 predictores en Hald",
              uiOutput("tab2_tabla_subconjuntos_hald_html")
            ),
            panel_ucm(
              "Selección automática: forward, backward y stepwise",
              p(class = "mb-2", style = paste0("font-size:0.82rem; color:", ucm$gris_texto, ";"),
                "Variables retenidas por cada método según el criterio AIC de la función step:"),
              uiOutput("tab2_tabla_auto_hald")
            ),
            panel_ucm(
              "Regresión por componentes principales: comparación entre MCO y PCR con 2 componentes", clase = "mb-0",
              uiOutput("tab2_tabla_pcr_hald")
            )
          )
        )
      )
    )
  ),
  
  # ========================================================================== #
  # Pestaña 3: gráficos y supuestos                                            #
  # ========================================================================== #
  nav_panel(
    title = "Gráficos y supuestos",
    div(
      class = "container-fluid py-3",
      
      selector_base(
        "tab3_dataset",
        "Selección de base de datos: visualización gráfica",
        "Alterne entre el conjunto Cars93 y Hald para acceder a los módulos interactivos."
      ),
      
      navset_card_pill(
        id = "tab3_submodulo",
        
        # SUBPESTAÑA A: SIMULADOR DE ESCENARIOS
        nav_panel(
          title = "Simulador de escenarios con covariables",
          div(
            class = "py-3",
            div(
              class = "row g-3",
              
              conditionalPanel(
                condition = "input.tab3_dataset == 'cars'",
                div(
                  class = "col-lg-4",
                  panel_ucm(
                    "Covariables de los modelos de Cars93", clase = "h-100 mb-0",
                    selectInput("sim_auto", "Cargar un auto real de Cars93:", choices = opciones_auto_cars, width = "100%"),
                    uiOutput("sim_auto_real"),
                    slider_cars("sim_minprice", "Precio mínimo (Min.Price, miles USD):", "Min.Price", 0.1),
                    slider_cars("sim_maxprice", "Precio máximo (Max.Price, miles USD):", "Max.Price", 0.1),
                    div(id = "sim_price", class = "shiny-input-container mb-3",
                        tags$label(class = "control-label", "Precio promedio (Price, miles USD):"),
                        div(class = "fw-semibold", textOutput("sim_price_valor", inline = TRUE),
                            span(style = paste0("font-size:0.78rem; font-weight:400; color:", ucm$gris_texto, ";"),
                                 " (promedio entre el mínimo y el máximo)"))),
                    slider_cars("sim_engine", "Cilindrada (EngineSize, litros):", "EngineSize", 0.1),
                    slider_cars("sim_hp", "Potencia (Horsepower, hp):", "Horsepower", 5),
                    slider_cars("sim_rpm", "Revoluciones máximas (RPM):", "RPM", 50),
                    slider_cars("sim_rev", "Revoluciones por milla (Rev.per.mile):", "Rev.per.mile", 25),
                    slider_cars("sim_tank", "Capacidad del estanque (Fuel.tank.capacity, galones):", "Fuel.tank.capacity", 0.5),
                    slider_cars("sim_pass", "Pasajeros (Passengers):", "Passengers", 1),
                    slider_cars("sim_wheel", "Distancia entre ejes (Wheelbase, pulgadas):", "Wheelbase", 1),
                    slider_cars("sim_weight", "Peso del vehículo (Weight, libras):", "Weight", 50),
                    slider_cars("sim_length", "Largo (Length, pulgadas):", "Length", 1),
                    slider_cars("sim_width", "Ancho (Width, pulgadas):", "Width", 1),
                    slider_cars("sim_turn", "Diámetro de giro (Turn.circle, pies):", "Turn.circle", 1),
                    slider_cars("sim_rear", "Espacio del asiento trasero (Rear.seat.room, pulgadas):", "Rear.seat.room", 0.5),
                    slider_cars("sim_luggage", "Espacio de equipaje (Luggage.room, pies cúbicos):", "Luggage.room", 1),
                    div(
                      class = "d-grid mt-2",
                      actionButton("btn_reset_cars", "Restaurar valores promedio", class = "btn btn-outline-secondary btn-sm")
                    )
                  )
                )
              ),
              
              conditionalPanel(
                condition = "input.tab3_dataset == 'hald'",
                div(
                  class = "col-lg-4",
                  panel_ucm(
                    "Composición química del cemento", clase = "h-100 mb-0",
                    selectInput("sim_ensayo", "Cargar un ensayo real de Hald:", choices = opciones_ensayo_hald, width = "100%"),
                    uiOutput("sim_ensayo_real"),
                    slider_hald("sim_x1", "X1 - aluminato tricálcico (%):", "X1", 0.5),
                    slider_hald("sim_x2", "X2 - silicato tricálcico (%):", "X2", 0.5),
                    slider_hald("sim_x3", "X3 - ferrito tetracálcico (%):", "X3", 0.5),
                    slider_hald("sim_x4", "X4 - silicato dicálcico (%):", "X4", 0.5),
                    div(
                      class = "d-grid mt-2",
                      actionButton("btn_reset_hald", "Restaurar mezcla promedio", class = "btn btn-outline-secondary btn-sm")
                    )
                  )
                )
              ),
              
              div(
                class = "col-lg-8",
                panel_ucm(
                  "Selección del modelo: elija con cuál calcular el escenario", clase = "mb-3",
                  p(class = "mb-2", style = paste0("font-size:0.83rem; color:", ucm$gris_texto, ";"),
                    "Puede probar cualquiera de los modelos que aparecen en las tablas de métricas. El modelo elegido por el análisis viene marcado:"),
                  uiOutput("ui_selector_modelo_sim"),
                  uiOutput("sim_nota_variables")
                ),
                div(
                  class = "row g-3 mb-3",
                  div(class = "col-sm-6", uiOutput("sim_kpi_1")),
                  div(class = "col-sm-6", uiOutput("sim_kpi_2"))
                ),
                panel_ucm(
                  "Comparación de predicciones ante el escenario configurado", clase = "mb-3",
                  uiOutput("sim_tabla_comparativa_html")
                ),
                panel_ucm(
                  "Posición del valor simulado sobre la distribución observada", clase = "mb-3",
                  plotlyOutput("sim_grafico_distribucion", height = "320px")
                ),
                panel_ucm(
                  "Qué variable pesa más en la predicción", clase = "mb-3",
                  p(class = "mb-2", style = paste0("font-size:0.83rem; color:", ucm$gris_texto, ";"),
                    "Muestra cuánto cambia la predicción al pasar cada variable de un valor bajo a uno alto de la muestra (percentiles 10 y 90), manteniendo el resto del escenario. Las barras hacia la derecha la aumentan y las de la izquierda la reducen."),
                  plotlyOutput("sim_grafico_tornado", height = "330px"),
                  uiOutput("sim_nota_tornado")
                ),
                conditionalPanel(
                  condition = "input.tab3_dataset == 'cars'",
                  panel_ucm(
                    "Conversión a unidades métricas: de millas por galón a kilómetros por litro", clase = "mb-3",
                    uiOutput("sim_seccion_km")
                  ),
                  panel_ucm(
                    "Costo de combustible en Chile: estimación con bencina de 95 octanos", clase = "mb-0",
                    div(
                      class = "p-3 rounded-2 mb-3 caja-destacada",
                      h6(class = "fw-bold texto-ucm", "Datos de Estados Unidos, lectura para Chile"),
                      p(class = "mb-0", style = "font-size:0.86rem; line-height:1.6;",
                        "Los datos de Cars93 provienen del mercado estadounidense de 1993, pero varias de estas marcas y modelos (o sus versiones actuales) también se comercializan en Chile. Por eso el costo se estima con un valor aproximado de la bencina de 95 octanos. Es una referencia: el rendimiento real depende del año, la versión y la forma de manejar.")
                    ),
                    div(
                      class = "row g-3 mb-2",
                      div(class = "col-sm-6",
                          numericInput("costo_precio", "Precio del litro de bencina de 95 octanos (CLP):",
                                       value = 1600, min = 500, max = 5000, step = 10, width = "100%")),
                      div(class = "col-sm-6",
                          sliderInput("costo_km", "Kilómetros recorridos al año:",
                                      min = 5000, max = 40000, value = 15000, step = 1000, width = "100%"))
                    ),
                    uiOutput("sim_costo_tarjetas"),
                    plotlyOutput("sim_grafico_costo", height = "300px")
                  )
                )
              )
            )
          )
        ),
        
        # SUBPESTAÑA B: DIAGNÓSTICO DE SUPUESTOS
        nav_panel(
          title = "Diagnóstico de supuestos estadísticos",
          div(
            class = "py-3",
            
            card(
              class = "mb-3 border-0 shadow-sm",
              card_body(
                class = "py-3",
                h6(class = "fw-bold texto-ucm mb-2", "Seleccione el modelo"),
                p(class = "mb-3", style = paste0("font-size:0.83rem; color:", ucm$gris_texto, ";"),
                  "Elija el modelo que desea auditar:"),
                uiOutput("ui_selector_modelo_supuestos")
              )
            ),
            
            div(
              class = "row g-3 mb-3",
              div(class = "col-lg-6",
                  panel_ucm("1. Linealidad: residuos frente a valores ajustados", clase = "h-100 mb-0",
                            plotlyOutput("plot_sup_linealidad", height = "280px"))),
              div(class = "col-lg-6",
                  panel_ucm("2. Normalidad: gráfico cuantil-cuantil", clase = "h-100 mb-0",
                            plotlyOutput("plot_sup_normalidad", height = "280px")))
            ),
            div(
              class = "row g-3",
              div(class = "col-lg-6",
                  panel_ucm("3. Homocedasticidad: escala y localización", clase = "h-100 mb-0",
                            plotlyOutput("plot_sup_homocedasticidad", height = "280px"))),
              div(class = "col-lg-6",
                  panel_ucm("4. Apalancamiento: residuos estandarizados frente a leverage", clase = "h-100 mb-0",
                            plotlyOutput("plot_sup_leverage", height = "280px")))
            )
          )
        ),
        
        # SUBPESTAÑA C: INFLUYENTES Y ATÍPICOS
        nav_panel(
          title = "Observaciones influyentes y apalancamiento",
          div(
            class = "py-3",
            
            card(
              class = "mb-3 border-0 shadow-sm",
              card_body(
                class = "py-3",
                h6(class = "fw-bold texto-ucm mb-2", "Seleccione el modelo"),
                p(class = "mb-3", style = paste0("font-size:0.83rem; color:", ucm$gris_texto, ";"),
                  "Los tres gráficos y el análisis de sensibilidad corresponden al modelo que elija aquí:"),
                uiOutput("ui_selector_modelo_influencia")
              )
            ),
            div(
              class = "row g-3 mb-3",
              div(class = "col-md-4",
                  panel_ucm("Distancia de Cook: corte crítico en 4/n", clase = "h-100 mb-0",
                            plotlyOutput("plot_diag_cook", height = "260px"))),
              div(class = "col-md-4",
                  panel_ucm("Apalancamiento hii: corte crítico en 2p/n", clase = "h-100 mb-0",
                            plotlyOutput("plot_diag_hii", height = "260px"))),
              div(class = "col-md-4",
                  panel_ucm("Residuos estudentizados: banda de atípicos [-2, 2]", clase = "h-100 mb-0",
                            plotlyOutput("plot_diag_rstud", height = "260px")))
            ),
            
            panel_ucm(
              "Análisis de sensibilidad: reajuste de coeficientes sin observaciones influyentes", clase = "mb-0",
              uiOutput("tab3_texto_sensibilidad"),
              uiOutput("tab3_tabla_sensibilidad")
            )
          )
        ),
        
        # SUBPESTAÑA D: REGULARIZACIÓN (TRAZAS Y CV)
        nav_panel(
          title = "Trazas y curvas de regularización",
          div(
            class = "py-3",
            
            card(
              class = "mb-3 border-0 shadow-sm",
              card_body(
                class = "py-2",
                div(
                  class = "d-flex flex-wrap align-items-center gap-3",
                  span(class = "fw-bold texto-ucm", "Seleccione la penalización a graficar"),
                  radioButtons(
                    inputId  = "tab3_reg_metodo",
                    label    = NULL,
                    choices  = c("Lasso (alpha = 1)" = "lasso",
                                 "Ridge (alpha = 0)" = "ridge",
                                 "Elastic Net (alpha = 0.5)" = "enet"),
                    selected = "lasso",
                    inline   = TRUE
                  )
                ),
                conditionalPanel(
                  condition = "input.tab3_dataset == 'cars'",
                  p(class = "mb-1", style = paste0("font-size:0.80rem; color:", ucm$gris_texto, ";"),
                    "En Cars93 las trazas aplican el procedimiento del script sobre la muestra completa estandarizada, con validación cruzada de 10 pliegues.")
                ),
                conditionalPanel(
                  condition = "input.tab3_dataset == 'hald'",
                  p(class = "mb-1", style = paste0("font-size:0.80rem; color:", ucm$gris_texto, ";"),
                    "En Hald la validación cruzada es de tipo leave-one-out (LOOCV), con 13 pliegues de una observación.")
                )
              )
            ),
            
            div(
              class = "row g-3 mb-3",
              div(class = "col-lg-6",
                  panel_ucm("Trayectoria de coeficientes: evolución frente a log(lambda)", clase = "h-100 mb-0",
                            plotlyOutput("plot_reg_trazas", height = "320px"))),
              div(class = "col-lg-6",
                  panel_ucm("Validación cruzada: error cuadrático medio frente a log(lambda)", clase = "h-100 mb-0",
                            plotlyOutput("plot_reg_cv", height = "320px")))
            ),
            
            conditionalPanel(
              condition = "input.tab3_dataset == 'hald'",
              panel_ucm(
                "Parámetros de regularización: lambda y grados de libertad efectivos en Hald", clase = "mb-0",
                uiOutput("tab3_tabla_lambda_hald")
              )
            )
          )
        ),
        
        # SUBPESTAÑA E: ERROR FUERA DE MUESTRA
        nav_panel(
          title = "Distribución de error fuera de muestra",
          div(
            class = "py-3",
            
            conditionalPanel(
              condition = "input.tab3_dataset == 'cars'",
              panel_ucm(
                sprintf("Error al predecir autos nuevos: cuánto se equivoca cada método en %d pruebas", B_rep), clase = "mb-0",
                p(class = "mb-2", style = paste0("font-size:0.83rem; color:", ucm$gris_texto, ";"),
                  "Cada caja resume las pruebas en las que el método predijo autos que no había visto. Mientras más abajo y más angosta esté la caja, más preciso y estable es el método:"),
                plotlyOutput("plot_reg_box_rmse_cars", height = "340px")
              )
            ),
            
            conditionalPanel(
              condition = "input.tab3_dataset == 'hald'",
              panel_ucm(
                "Error al predecir ensayos nuevos: cuánto se equivoca cada método", clase = "mb-0",
                p(class = "mb-2", style = paste0("font-size:0.83rem; color:", ucm$gris_texto, ";"),
                  "Cada barra muestra cuánto se equivoca el método al predecir un ensayo que se dejó fuera del ajuste. Se comparan dos niveles de penalización: el de menor error y el más simple que rinde casi igual. Mientras más baja la barra, mejor la predicción:"),
                plotlyOutput("plot_reg_cv_hald_barras", height = "340px")
              )
            )
          )
        )
      )
    )
  ),
  
  nav_spacer(),
  nav_item(
    tags$div(
      class = "py-2 px-3 pie-ucm",
      span("Universidad Católica del Maule | Modelos estadísticos activos")
    )
  )
)

# ---------------------------------------------------------------------------- #
# 4. Servidor (Server)                                                          #
# ---------------------------------------------------------------------------- #

server <- function(input, output, session) {
  
  # Sincronizar selección de base de datos entre pestañas
  observeEvent(input$tab1_dataset, {
    updateRadioButtons(session, "tab2_dataset", selected = input$tab1_dataset)
    updateRadioButtons(session, "tab3_dataset", selected = input$tab1_dataset)
  })
  
  observeEvent(input$tab2_dataset, {
    updateRadioButtons(session, "tab1_dataset", selected = input$tab2_dataset)
    updateRadioButtons(session, "tab3_dataset", selected = input$tab2_dataset)
  })
  
  observeEvent(input$tab3_dataset, {
    updateRadioButtons(session, "tab1_dataset", selected = input$tab3_dataset)
    updateRadioButtons(session, "tab2_dataset", selected = input$tab3_dataset)
  })
  
  # ========================================================================== #
  # PESTAÑA 1: TABLAS                                                          #
  # ========================================================================== #
  
  output$tab1_tabla_pares_cars <- renderUI({
    tabla_ucm(pares_altos_cars)
  })
  
  output$tab2_titulo_ecuacion_cars <- renderUI({
    req(input$tab2_cars_modelo_sel)
    h6(class = "fw-bold texto-ucm",
       if (input$tab2_cars_modelo_sel == clave_elegido_cars) "Ecuación del modelo elegido" else "Ecuación del modelo seleccionado")
  })

  output$tab2_ecuacion_cars <- renderUI({
    req(input$tab2_cars_modelo_sel)
    cl <- input$tab2_cars_modelo_sel
    req(cl %in% names(etq_sim_cars))
    div(class = "py-2 px-3 text-center caja-ecuacion",
        tags$div(style = "font-size:0.84rem; font-weight:600; margin-bottom:0.4rem;", etq_sim_cars[[cl]]),
        tags$div(style = "overflow-wrap:anywhere; white-space:normal;", ecuacion_cars_modelo(cl)))
  })

  output$tab2_titulo_ecuacion_hald <- renderUI({
    req(input$tab2_hald_modelo_ec)
    h6(class = "fw-bold texto-ucm",
       if (input$tab2_hald_modelo_ec == clave_elegido_hald) "Ecuación del modelo elegido" else "Ecuación del modelo alternativo")
  })
  
  output$tab2_ecuacion_hald <- renderUI({
    req(input$tab2_hald_modelo_ec)
    cl <- input$tab2_hald_modelo_ec
    req(cl %in% names(etq_sim_hald))
    act   <- vars_activas_hald[[cl]]
    fuera <- setdiff(names(x_hald), act)
    txt_vars <- if (cl == "hald_pcr") {
      sprintf("Este modelo resume las 4 variables en %d componentes, por lo que todas intervienen.", ncomp_opt)
    } else if (length(fuera) == 0) {
      "Este modelo conserva las 4 variables."
    } else {
      sprintf("Este modelo conserva %s y deja fuera %s.", paste(act, collapse = ", "), paste(fuera, collapse = ", "))
    }
    estilo <- paste0("font-size:0.82rem; color:", ucm$gris_texto, ";")
    tagList(
      p(class = "mb-2", style = "font-size:0.84rem;", strong(etq_sim_hald[[cl]])),
      div(class = "py-2 px-3 text-center caja-ecuacion",
          tags$div(style = "font-size:0.84rem; font-weight:600; margin-bottom:0.4rem;", etq_sim_hald[[cl]]),
          tags$div(style = "overflow-wrap:anywhere; white-space:normal;", ecuacion_hald_modelo(cl))),
      p(class = "mt-2 mb-1", style = estilo,
        "En la ecuación, y es el calor liberado a los 180 días (cal/g) y X1, X2, X3 y X4 son los porcentajes en peso del aluminato tricálcico, el silicato tricálcico, el ferrito tetracálcico y el silicato dicálcico, respectivamente."),
      p(class = "mb-0", style = estilo, txt_vars,
        " Cada coeficiente indica cuántas calorías por gramo cambia el calor liberado por cada punto porcentual adicional del compuesto, manteniendo los demás constantes.")
    )
  })
  
  output$tab1_tabla_pares_hald <- renderUI({
    tabla_ucm(pares_hald_tab, fila_resaltada = idx_muy_alta_hald)
  })
  
  output$tab1_tabla_datos_cars <- renderUI({
    df <- head(d_cars, 8)
    df[] <- lapply(df, function(col) format(col, trim = TRUE))
    tabla_ucm(df)
  })
  
  output$tab1_tabla_datos_hald <- renderUI({
    df <- datos_hald
    df$`Suma de masa` <- rowSums(datos_hald[, 2:5])
    df[] <- lapply(df, function(col) format(col, trim = TRUE))
    tabla_ucm(df)
  })
  
  # ========================================================================== #
  # PESTAÑA 2: MÉTRICAS DINÁMICAS Y TABLAS                                     #
  # ========================================================================== #
  
  # Tarjetas dinámicas de indicadores para Cars93
  output$tab2_kpi_cards_cars <- renderUI({
    req(input$tab2_cars_modelo_sel)
    cl <- input$tab2_cars_modelo_sel
    act_vars <- vars_activas_cars[[cl]]
    is_reg <- cl %in% c("cars_ridge", "cars_lasso", "cars_enet")

    if (!is_reg) {
      m <- switch(cl, cars_bic = mod_bic_cars, cars_aic = mod_aic_cars,
                  cars_ols = mod_inicial_cars, cars_cp = mod_cp_cars)
      X <- model.matrix(m)[, -1, drop = FALSE]
      rmse_val <- sqrt(mean(resid(m)^2))
      ck <- cooks.distance(m)
      max_ck <- max(ck, na.rm = TRUE)
      umb_ck <- 4 / n_cars
      n_over <- sum(ck > umb_ck, na.rm = TRUE)
    } else {
      # Para regularización, los diagnósticos se calculan sobre los predictores
      # con coeficiente no nulo; Cook se obtiene de un ajuste auxiliar explícito.
      if (!length(act_vars)) act_vars <- names(d_cars)[-1][1]
      m <- lm(reformulate(act_vars, response = "MPG.city"), data = d_cars)
      X <- as.matrix(d_cars[, act_vars, drop = FALSE])
      rmse_map <- c(cars_ridge = "Ridge", cars_lasso = "Lasso", cars_enet = "ElasticNet")
      rmse_val <- unname(prom_rmse_cars[rmse_map[[cl]]])
      ck <- cooks.distance(m)
      max_ck <- max(ck, na.rm = TRUE)
      umb_ck <- 4 / n_cars
      n_over <- sum(ck > umb_ck, na.rm = TRUE)
    }

    p_cnt <- ncol(X)
    if (p_cnt > 1) {
      vifs <- tryCatch(car::vif(m), error = function(e) NA_real_)
      v_max <- if (all(is.na(vifs))) NA_real_ else max(vifs, na.rm = TRUE)
      v_nom <- if (all(is.na(vifs))) "No calculable" else names(which.max(vifs))
      R <- suppressWarnings(cor(X))
      ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
      C_val <- if (min(ev) <= 1e-10) Inf else sqrt(max(ev) / min(ev))
    } else {
      v_max <- 1
      v_nom <- if (p_cnt == 1) colnames(X)[1] else "Sin predictores activos"
      C_val <- 1
    }
    v_lvl <- nivel_vif(v_max)
    c_lvl <- nivel_cond(C_val)

    card_vif <- tarjeta_metrica(if (is_reg) "VIF variables activas" else "VIF máximo",
      if (is.finite(v_max)) sprintf("VIF = %.2f", v_max) else "No calculable",
      v_lvl$txt, sprintf("Variable con mayor VIF: %s", v_nom), v_lvl$nivel)
    card_cond <- tarjeta_metrica("Número de condición",
      if (is.finite(C_val)) sprintf("C = %.1f", C_val) else "C = ∞",
      c_lvl$txt, sprintf("Calculado con %d predictores activos", p_cnt), c_lvl$nivel)
    card_cook <- tarjeta_metrica("Distancia de Cook máx.", sprintf("Cook = %.3f", max_ck),
      sprintf("Corte 4/n = %.3f", umb_ck), sprintf("%d observaciones sobre el umbral", n_over),
      if (n_over > 0) "warning" else "success")
    card_rmse <- tarjeta_metrica(if (is_reg) "Error fuera de muestra" else "RMSE en muestra",
      sprintf("RMSE = %.3f mpg", rmse_val),
      if (cl == clave_elegido_cars) "Modelo elegido" else "Modelo seleccionado",
      if (is_reg) sprintf("Promedio de %d particiones de prueba", B_rep) else "Error residual del modelo ajustado",
      if (cl == clave_elegido_cars) "success" else "primary")

    div(class = "row g-3 mb-4",
      div(class = "col-sm-6 col-xl-3", card_vif),
      div(class = "col-sm-6 col-xl-3", card_cond),
      div(class = "col-sm-6 col-xl-3", card_cook),
      div(class = "col-sm-6 col-xl-3", card_rmse))
  })

  # Tarjetas dinámicas de indicadores para Hald  # Tarjetas dinámicas de indicadores para Hald
  output$tab2_kpi_cards_hald <- renderUI({
    req(input$tab2_hald_modelo_ec)
    cl <- input$tab2_hald_modelo_ec
    
    if (cl %in% c("hald_1v", "hald_2v", "hald_3v", "hald_fwd", "hald_bwd", "hald_step", "hald_ols")) {
      m <- switch(cl,
                  hald_1v   = mod_hald_1v,
                  hald_2v   = mod_hald_2v,
                  hald_3v   = mod_hald_3v,
                  hald_fwd  = mod_forward_hald,
                  hald_bwd  = mod_backward_hald,
                  hald_step = mod_stepwise_hald,
                  hald_ols  = mod_completo_hald)
      
      p_cnt <- length(coef(m)) - 1
      if (p_cnt > 1) {
        vifs  <- car::vif(m)
        v_max <- max(vifs)
        v_nom <- names(which.max(vifs))
        v_lvl <- nivel_vif(v_max)
        X <- model.matrix(m)[, -1, drop = FALSE]
        R <- cor(X)
        ev <- eigen(R, symmetric = TRUE)$values
        C_val <- sqrt(max(ev) / min(ev))
        c_lvl <- nivel_cond(C_val)
      } else {
        v_max <- 1.0; v_nom <- names(coef(m))[-1]; v_lvl <- list(txt = "Bajo (umbral < 5)", nivel = "success")
        C_val <- 1.0; c_lvl <- list(txt = "Bajo (umbral < 15)", nivel = "success")
      }
      
      ck <- cooks.distance(m)
      max_ck <- max(ck)
      umb_ck <- 4 / n_hald
      n_over <- sum(ck > umb_ck)
      rmse_val <- sqrt(mean(resid(m)^2))
      
      card_vif  <- tarjeta_metrica("VIF máximo", sprintf("VIF = %.3f", v_max), v_lvl$txt, sprintf("Variable con mayor colinealidad: %s", v_nom), v_lvl$nivel)
      card_cond <- tarjeta_metrica("Número de condición", sprintf("C = %.2f", C_val), c_lvl$txt, sprintf("Matriz de %d compuestos", p_cnt), c_lvl$nivel)
      card_cook <- tarjeta_metrica("Distancia de Cook máx.", sprintf("Cook = %.3f", max_ck), sprintf("Corte 4/n = %.3f", umb_ck), sprintf("%d ensayos superan el umbral", n_over), if(n_over > 0) "warning" else "success")
      card_rmse <- tarjeta_metrica("Error RMSE en muestra", sprintf("RMSE = %.3f cal/g", rmse_val), if(cl == clave_elegido_hald) "Modelo elegido (criterio experto)" else "Evaluado", "Ajuste muestral", if(cl == clave_elegido_hald) "success" else "primary")
      
    } else if (cl == "hald_pcr") {
      pred_pcr <- as.numeric(predict(mod_pcr_hald, newdata = datos_hald, ncomp = ncomp_opt))
      rmse_pcr <- sqrt(mean((y_hald - pred_pcr)^2))
      scores_pcr <- as.data.frame(scores(mod_pcr_hald)[, seq_len(ncomp_opt), drop = FALSE])
      names(scores_pcr) <- paste0("PC", seq_len(ncomp_opt))
      diag_pcr <- lm(y_hald ~ ., data = data.frame(y_hald = y_hald, scores_pcr))
      cook_pcr <- cooks.distance(diag_pcr)
      over_pcr <- sum(cook_pcr > 4 / n_hald)
      card_vif  <- tarjeta_metrica("VIF máximo", sprintf("VIF = %.3f", max(car::vif(diag_pcr))), "Componentes principales", "VIF calculado sobre los componentes retenidos", "success")
      card_cond <- tarjeta_metrica("Número de condición", sprintf("C = %.2f", calc_diag_metrics(diag_pcr)[["C"]]), "Componentes principales", "Número de condición de los componentes retenidos", "success")
      card_cook <- tarjeta_metrica("Distancia de Cook máx.", sprintf("Cook = %.3f", max(cook_pcr)), sprintf("Corte 4/n = %.3f", 4 / n_hald), sprintf("%d ensayos superan el umbral", over_pcr), if (over_pcr > 0) "warning" else "success")
      card_rmse <- tarjeta_metrica("Error RMSE en muestra", sprintf("RMSE = %.3f cal/g", rmse_pcr), "Modelo seleccionado", "Error calculado con las predicciones PCR", "primary")
      
    } else {
      cv_obj <- switch(cl, hald_ridge = cv_ridge_hald, hald_lasso = cv_lasso_hald, hald_enet = cv_enet_hald)
      b <- extract_coef(cv_obj, "lambda.1se")[-1]
      act_vars <- names(b)[b != 0]
      
      m_aux <- modelos_diag_hald[[cl]]
      diag_act <- calc_diag_metrics(m_aux)
      v_max <- unname(diag_act[["VIF"]])
      C_val <- unname(diag_act[["C"]])
      v_lvl <- nivel_vif(v_max)
      c_lvl <- nivel_cond(C_val)
      
      rmse_cv <- sqrt(cv_obj$cvm[cv_obj$index["1se", ]])
      
      card_vif  <- tarjeta_metrica("VIF variables activas", sprintf("VIF = %.3f", v_max), "Penalizado", sprintf("%d compuestos activos; λ = %.3f", length(act_vars), cv_obj$lambda.1se), v_lvl$nivel)
      card_cond <- tarjeta_metrica("Número de condición", sprintf("C = %.2f", C_val), c_lvl$txt, "Evaluado sobre compuestos activos", c_lvl$nivel)
      cook_aux <- cooks.distance(m_aux)
      max_cook_aux <- max(cook_aux, na.rm = TRUE)
      n_cook_aux <- sum(cook_aux > 4 / n_hald, na.rm = TRUE)
      card_cook <- tarjeta_metrica("Distancia de Cook máx.", sprintf("Cook = %.3f", max_cook_aux),
        sprintf("Corte 4/n = %.3f", 4 / n_hald),
        sprintf("%d ensayos superan el umbral (ajuste auxiliar)", n_cook_aux),
        if (n_cook_aux > 0) "warning" else "success")
      card_rmse <- tarjeta_metrica("Error LOOCV fuera de muestra", sprintf("RMSE = %.3f cal/g", rmse_cv),
        if(cl == clave_elegido_hald) "Modelo elegido (criterio experto)" else "Evaluado",
        "Validación cruzada Leave-One-Out", if(cl == clave_elegido_hald) "success" else "primary")
    }
    
    div(
      class = "row g-3 mb-4",
      div(class = "col-sm-6 col-xl-3", card_vif),
      div(class = "col-sm-6 col-xl-3", card_cond),
      div(class = "col-sm-6 col-xl-3", card_cook),
      div(class = "col-sm-6 col-xl-3", card_rmse)
    )
  })
  
  output$tab2_tabla_modelos_cars <- renderUI({
    tabla_ucm(tabla_modelos_cars, fila_resaltada = which(tabla_modelos_cars$Modelo == "Modelo BIC" & ganador_cars == "BIC_Fase3"),
              insignia = ifelse(tabla_modelos_cars$Modelo == "Modelo BIC" & ganador_cars == "BIC_Fase3",
                                "Menor error fuera de muestra", ""))
  })
  
  output$tab2_tabla_coef_bic <- renderUI({
    tabla_ucm(tabla_coef_bic)
  })
  
  output$tab2_tabla_comp_bic_cp <- renderUI({
    tabla_ucm(comp_bic_cp_cars)
  })
  
  output$tab2_tabla_rmse_cars_html <- renderUI({
    proms <- prom_rmse_cars
    df <- data.frame(
      "Método / modelo"       = unname(etiquetas_rmse[names(proms)]),
      "RMSE fuera de muestra" = sprintf("%.4f mpg", proms),
      "Diferencia frente a MCO" = sprintf("%+.4f mpg", proms - proms[["OLS_Completo"]]),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    es_ganador <- names(proms) == ganador_cars
    tabla_ucm(df, fila_resaltada = which(es_ganador), insignia = ifelse(es_ganador, "Menor error", ""))
  })
  
  output$tab2_tabla_subconjuntos_hald_html <- renderUI({
    etiq <- rep("", 4)
    for (i in 1:4) {
      partes <- c(
        if (i == mejor_cp_hald && i == mejor_bic_hald) "Mínimo Cp y BIC"
        else if (i == mejor_cp_hald) "Mínimo Cp"
        else if (i == mejor_bic_hald) "Mínimo BIC",
        if (i == mejor_r2adj_hald) "Máximo R² ajustado"
      )
      etiq[i] <- paste(partes, collapse = " | ")
    }
    tabla_ucm(tabla_subconjuntos_hald, fila_resaltada = which(nzchar(etiq)), insignia = etiq)
  })
  
  output$tab2_tabla_auto_hald <- renderUI({
    tabla_ucm(tabla_auto_hald)
  })
  
  output$tab2_tabla_comp_hald_html <- renderUI({
    tabla_ucm(tabla_hald_comp, col_resaltadas = c(4, 5), fila_resaltada = nrow(tabla_hald_comp))
  })
  
  output$tab2_tabla_pcr_hald <- renderUI({
    tabla_ucm(comp_mco_pcr)
  })
  
  output$tab2_tabla_indicadores_cars <- renderUI({
    tabla_ucm(tabla_ind_cars, fila_resaltada = 3, insignia = c("", "", "Modelo seleccionado", ""))
  })
  
  output$tab2_tabla_indicadores_hald <- renderUI({
    tabla_ucm(tabla_ind_hald,
              fila_resaltada = which(grepl("Lasso", tabla_ind_hald$Modelo, ignore.case = TRUE)),
              insignia = ifelse(grepl("Lasso", tabla_ind_hald$Modelo, ignore.case = TRUE),
                                "Modelo elegido (criterio experto)", ""))
  })
  
  # ========================================================================== #
  # PESTAÑA 3: SIMULADOR                                                       #
  # ========================================================================== #
  
  observeEvent(input$btn_reset_cars, {
    for (id in names(def_cars)) updateSliderInput(session, id, value = def_cars[[id]])
    updateSelectInput(session, "sim_auto", selected = "")
  })
  
  observeEvent(input$btn_reset_hald, {
    for (id in names(def_hald)) updateSliderInput(session, id, value = def_hald[[id]])
    updateSelectInput(session, "sim_ensayo", selected = "")
  })
  
  observeEvent(input$sim_minprice, {
    if (input$sim_minprice > input$sim_maxprice)
      updateSliderInput(session, "sim_maxprice", value = input$sim_minprice)
  }, ignoreInit = TRUE)
  
  observeEvent(input$sim_maxprice, {
    if (input$sim_maxprice < input$sim_minprice)
      updateSliderInput(session, "sim_minprice", value = input$sim_maxprice)
  }, ignoreInit = TRUE)
  
  observeEvent(input$sim_auto, {
    req(input$sim_auto != "")
    fila <- d_cars[as.integer(input$sim_auto), ]
    for (v in names(ids_sim_cars)) {
      id <- ids_sim_cars[[v]]
      if (is.null(def_cars[[id]])) next
      valor <- round(fila[[v]] / paso_cars[[id]]) * paso_cars[[id]]
      valor <- min(max(valor, rango_cars[[id]][1]), rango_cars[[id]][2])
      updateSliderInput(session, id, value = round(valor, 4))
    }
  })
  
  output$sim_auto_real <- renderUI({
    req(input$sim_auto != "")
    i <- as.integer(input$sim_auto)
    p(class = "mb-3", style = paste0("font-size:0.8rem; color:", ucm$gris_texto, ";"),
      sprintf("Rendimiento real de %s: %g mpg (%.2f km/L). Los deslizadores se redondean a su paso.",
              nombres_cars[i], d_cars$MPG.city[i], d_cars$MPG.city[i] * 1.609344 / 3.785411784))
  })
  
  observeEvent(input$sim_ensayo, {
    req(input$sim_ensayo != "")
    fila <- x_hald[as.integer(input$sim_ensayo), ]
    for (v in names(ids_sim_hald)) updateSliderInput(session, ids_sim_hald[[v]], value = fila[[v]])
  })
  
  output$sim_ensayo_real <- renderUI({
    req(input$sim_ensayo != "")
    i <- as.integer(input$sim_ensayo)
    p(class = "mb-3", style = paste0("font-size:0.8rem; color:", ucm$gris_texto, ";"),
      sprintf("Calor real del ensayo %d: %g cal/g.", i, y_hald[i]))
  })
  
  output$sim_price_valor <- renderText({
    req(input$sim_minprice, input$sim_maxprice)
    sprintf("%.2f", (input$sim_minprice + input$sim_maxprice) / 2)
  })
  
  nd_cars <- reactive({
    data.frame(
      Min.Price          = input$sim_minprice,
      Price              = (input$sim_minprice + input$sim_maxprice) / 2,
      Max.Price          = input$sim_maxprice,
      EngineSize         = input$sim_engine,
      Horsepower         = input$sim_hp,
      RPM                = input$sim_rpm,
      Rev.per.mile       = input$sim_rev,
      Fuel.tank.capacity = input$sim_tank,
      Passengers         = input$sim_pass,
      Wheelbase          = input$sim_wheel,
      Weight             = input$sim_weight,
      Length             = input$sim_length,
      Width              = input$sim_width,
      Turn.circle        = input$sim_turn,
      Rear.seat.room     = input$sim_rear,
      Luggage.room       = input$sim_luggage
    )
  })
  
  nd_hald <- reactive({
    data.frame(X1 = input$sim_x1, X2 = input$sim_x2, X3 = input$sim_x3, X4 = input$sim_x4)
  })
  
  output$ui_selector_modelo_sim <- renderUI({
    req(input$tab3_dataset)
    if (input$tab3_dataset == "cars") {
      radioButtons("sim_modelo", label = NULL, width = "100%",
                   choices = do.call(c, unname(opciones_sim(etq_sim_cars, grupos_sim_cars, clave_elegido_cars))),
                   selected = clave_elegido_cars)
    } else {
      radioButtons("sim_modelo", label = NULL, width = "100%",
                   choices = do.call(c, unname(opciones_sim(etq_sim_hald, grupos_sim_hald, clave_elegido_hald))),
                   selected = clave_elegido_hald)
    }
  })
  
  modelo_sim <- reactive({
    req(input$tab3_dataset, input$sim_modelo)
    claves <- if (input$tab3_dataset == "cars") names(etq_sim_cars) else names(etq_sim_hald)
    req(input$sim_modelo %in% claves)
    input$sim_modelo
  })
  
  observe({
    req(input$tab3_dataset)
    cl <- modelo_sim()
    if (input$tab3_dataset == "cars") {
      act <- vars_activas_cars[[cl]]
      ids <- ids_sim_cars
    } else {
      act <- vars_activas_hald[[cl]]
      ids <- ids_sim_hald
    }
    session$sendCustomMessage("sim_activos", list(ids = unname(ids), activo = unname(names(ids) %in% act)))
  })
  
  output$sim_nota_variables <- renderUI({
    req(input$tab3_dataset)
    cl <- modelo_sim()
    estilo <- paste0("font-size:0.82rem; color:", ucm$gris_texto, ";")
    if (input$tab3_dataset == "cars") {
      act    <- vars_activas_cars[[cl]]
      n_ctrl <- sum(names(ids_sim_cars) %in% act)
      txt <- if (n_ctrl == length(ids_sim_cars)) {
        sprintf("Este modelo utiliza las %d covariables, por lo que todos los controles intervienen.", n_ctrl)
      } else {
        sprintf("Este modelo utiliza %d de las %d covariables. Los controles atenuados no influyen en la predicción.",
                n_ctrl, length(ids_sim_cars))
      }
    } else if (cl == "hald_pcr") {
      txt <- sprintf("Este modelo resume las 4 variables en %d componentes, por lo que todos los controles intervienen.", ncomp_opt)
    } else {
      act <- vars_activas_hald[[cl]]
      txt <- sprintf("Este modelo utiliza %d de las 4 variables (%s). Los controles atenuados no influyen en la predicción.",
                     length(act), paste(act, collapse = ", "))
    }
    p(class = "mt-2 mb-0", style = estilo, txt)
  })
  
  pred_cars_sel <- reactive({
    req(input$sim_minprice)
    pred_sim_cars(modelo_sim(), nd_cars())
  })
  
  pred_hald_sel <- reactive({
    req(input$sim_x1)
    pred_sim_hald(modelo_sim(), nd_hald())
  })
  
  output$sim_grafico_tornado <- renderPlotly({
    req(input$tab3_dataset)
    cl <- modelo_sim()
    if (input$tab3_dataset == "cars") {
      req(input$sim_minprice)
      act <- vars_activas_cars[[cl]]
      nd  <- nd_cars()
      etq <- c(EngineSize = "Cilindrada", Horsepower = "Potencia", RPM = "Revoluciones máximas",
               Rev.per.mile = "Revoluciones por milla", Fuel.tank.capacity = "Capacidad del estanque",
               Passengers = "Pasajeros", Wheelbase = "Distancia entre ejes", Weight = "Peso",
               Length = "Largo", Width = "Ancho", Turn.circle = "Diámetro de giro",
               Rear.seat.room = "Espacio del asiento trasero", Luggage.room = "Espacio de equipaje")
      grupos <- c(list("Precio" = c("Min.Price", "Max.Price")), setNames(as.list(names(etq)), unname(etq)))
      ef <- vapply(names(grupos), function(g) {
        vv <- grupos[[g]]
        if (!any(c(vv, if (g == "Precio") "Price") %in% act)) return(NA_real_)
        lo <- nd
        hi <- nd
        for (v in vv) {
          q <- quantile(d_cars[[v]], c(0.1, 0.9))
          lo[[v]] <- q[[1]]
          hi[[v]] <- q[[2]]
        }
        if (g == "Precio") {
          lo$Price <- (lo$Min.Price + lo$Max.Price) / 2
          hi$Price <- (hi$Min.Price + hi$Max.Price) / 2
        }
        pred_sim_cars(cl, hi) - pred_sim_cars(cl, lo)
      }, numeric(1))
      eje <- "Cambio en la eficiencia (mpg)"
    } else {
      req(input$sim_x1)
      act <- vars_activas_hald[[cl]]
      nd  <- nd_hald()
      etq <- c(X1 = "X1 (aluminato tricálcico)", X2 = "X2 (silicato tricálcico)",
               X3 = "X3 (ferrito tetracálcico)", X4 = "X4 (silicato dicálcico)")
      grupos <- setNames(as.list(names(etq)), unname(etq))
      ef <- vapply(names(grupos), function(g) {
        v <- grupos[[g]]
        if (!v %in% act) return(NA_real_)
        lo <- nd
        hi <- nd
        q <- quantile(x_hald[[v]], c(0.1, 0.9))
        lo[[v]] <- q[[1]]
        hi[[v]] <- q[[2]]
        pred_sim_hald(cl, hi) - pred_sim_hald(cl, lo)
      }, numeric(1))
      eje <- "Cambio en el calor liberado (cal/g)"
    }
    df <- data.frame(var = names(ef), efecto = as.numeric(ef), stringsAsFactors = FALSE)
    df <- df[!is.na(df$efecto), ]
    req(nrow(df) > 0)
    df <- head(df[order(-abs(df$efecto)), ], 8)
    df$var   <- factor(df$var, levels = rev(df$var))
    df$Sentido <- ifelse(df$efecto >= 0, "Aumenta", "Reduce")
    p <- ggplot(df, aes(x = var, y = efecto, fill = Sentido)) +
      geom_col(width = 0.7) +
      geom_hline(yintercept = 0, color = ucm$gris_texto, linewidth = 0.4) +
      coord_flip() +
      scale_fill_manual(values = c("Aumenta" = ucm$azul_medio, "Reduce" = ucm$rojo)) +
      labs(x = NULL, y = eje) +
      tema_graf(11)
    a_plotly(p)
  })
  
  output$sim_nota_tornado <- renderUI({
    req(input$tab3_dataset == "hald")
    p(class = "mt-2 mb-0", style = paste0("font-size:0.8rem; color:", ucm$gris_texto, ";"),
      "Cada compuesto se mueve por separado, sin respetar que la mezcla suma cerca de 97 %. Úselo solo como orientación.")
  })
  
  costo_datos <- reactive({
    req(input$tab3_dataset == "cars", input$sim_minprice, input$costo_precio, input$costo_km)
    req(input$costo_precio > 0)
    pred <- pred_cars_sel()
    req(pred > 0)
    conv <- 1.609344 / 3.785411784
    km_l <- pred * conv
    km_l_prom <- mean(d_cars$MPG.city) * conv
    list(precio = input$costo_precio, km = input$costo_km, km_l = km_l, km_l_prom = km_l_prom,
         c100 = 100 * input$costo_precio / km_l,
         anual = input$costo_km * input$costo_precio / km_l,
         anual_prom = input$costo_km * input$costo_precio / km_l_prom)
  })
  
  output$sim_costo_tarjetas <- renderUI({
    cd <- costo_datos()
    dif <- cd$anual - cd$anual_prom
    div(
      class = "row g-3 mb-3",
      div(class = "col-sm-6",
          tarjeta_metrica("Costo cada 100 km",
                          fmt_clp(cd$c100),
                          "Referencia Chile",
                          sprintf("Con un rendimiento de %.2f km/L y bencina de 95 octanos a %s por litro.", cd$km_l, fmt_clp(cd$precio)),
                          nivel = "primary")),
      div(class = "col-sm-6",
          tarjeta_metrica("Costo anual de combustible",
                          fmt_clp(cd$anual),
                          sprintf("%s km al año", format(cd$km, big.mark = ".")),
                          sprintf("%s al año %s que el auto promedio de Cars93 (%s).",
                                  fmt_clp(abs(dif)), if (dif > 0) "más" else "menos", fmt_clp(cd$anual_prom)),
                          nivel = if (dif > 0) "warning" else "success"))
    )
  })
  
  output$sim_grafico_costo <- renderPlotly({
    cd <- costo_datos()
    km <- seq(0, cd$km, length.out = 40)
    df <- rbind(
      data.frame(km = km, costo = km * cd$precio / cd$km_l, Auto = "Escenario simulado"),
      data.frame(km = km, costo = km * cd$precio / cd$km_l_prom, Auto = "Promedio de Cars93")
    )
    p <- ggplot(df, aes(x = km, y = costo, color = Auto, linetype = Auto)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("Escenario simulado" = ucm$azul, "Promedio de Cars93" = ucm$gris_texto)) +
      scale_x_continuous(labels = function(x) formatC(x, format = "d", big.mark = ".")) +
      scale_y_continuous(labels = function(x) paste0("$", formatC(round(x), format = "d", big.mark = "."))) +
      labs(x = "Kilómetros recorridos", y = "Costo acumulado de combustible (CLP)") +
      tema_graf(11)
    a_plotly(p)
  })
  
  output$sim_seccion_km <- renderUI({
    req(input$tab3_dataset == "cars", input$sim_minprice)
    pred <- pred_cars_sel()
    km_milla <- 1.609344
    l_galon  <- 3.785411784
    km_l <- pred * km_milla / l_galon
    estilo <- paste0("font-size:0.84rem; color:", ucm$gris_texto, "; line-height:1.6;")
    div(
      class = "row g-3",
      div(class = "col-sm-6",
          tarjeta_metrica("Rendimiento en km/litro",
                          if (pred > 0) sprintf("%.2f km/L", km_l) else "No válido",
                          "Equivalente métrico",
                          if (pred > 0) sprintf("Equivale a %.1f litros cada 100 km, con el mismo escenario y modelo de la eficiencia proyectada.", 100 / km_l)
                          else "La eficiencia proyectada no es válida para este escenario.",
                          nivel = "primary")),
      div(class = "col-sm-6",
          p(style = estilo, sprintf("Las millas por galón son la unidad original de Cars93. Una milla equivale a %.3f km y un galón estadounidense a %.3f litros, por lo que se multiplica por %.3f y se divide por %.3f (es decir, por %.4f).",
                                    km_milla, l_galon, km_milla, l_galon, km_milla / l_galon)),
          p(class = "mb-0", style = estilo, strong("En este escenario: "),
            sprintf("%.2f mpg × %.3f ÷ %.3f = %.2f km/L.", pred, km_milla, l_galon, km_l)))
    )
  })
  
  output$sim_kpi_1 <- renderUI({
    req(input$tab3_dataset)
    cl <- modelo_sim()
    if (input$tab3_dataset == "cars") {
      pred  <- round(pred_cars_sel(), 2)
      fuera <- pred < banda_cars[1] || pred > banda_cars[2]
      nivel <- if (fuera) "danger" else if (pred >= 24) "success" else if (pred >= 19) "info" else "warning"
      tarjeta_metrica("Eficiencia proyectada",
                      if (pred > 0) paste0(pred, " mpg") else "No válido",
                      if (fuera) "Fuera de lo observado" else if (cl == clave_elegido_cars) "Modelo elegido" else "Modelo alternativo",
                      paste0(etq_sim_cars[[cl]], ". ",
                             if (fuera) sprintf("El escenario queda fuera del rango observado (%g a %g mpg), por lo que la predicción no es confiable.", banda_cars[1], banda_cars[2])
                             else "Rendimiento en ciclo urbano según las covariables ingresadas."),
                      nivel = nivel)
    } else {
      pred_h <- round(pred_hald_sel(), 2)
      tarjeta_metrica("Calor proyectado",
                      paste0(pred_h, " cal/g"),
                      if (cl == clave_elegido_hald) "Modelo elegido" else "Modelo alternativo",
                      paste0(etq_sim_hald[[cl]], ". Calor acumulado liberado a 180 días de fraguado."),
                      nivel = "primary")
    }
  })
  
  output$sim_kpi_2 <- renderUI({
    req(input$tab3_dataset)
    if (input$tab3_dataset == "cars") {
      pred <- pred_cars_sel()
      pct  <- round(mean(d_cars$MPG.city < pred) * 100)
      tarjeta_metrica("Posición en la muestra",
                      paste0("Supera al ", pct, "%"),
                      "Rango percentil",
                      "Porcentaje de vehículos del conjunto Cars93 superados en eficiencia",
                      nivel = if (pct >= 50) "success" else "info")
    } else {
      req(input$sim_x1)
      act     <- vars_activas_hald[[modelo_sim()]]
      fuera_h <- setdiff(names(x_hald), act)
      suma_masa <- input$sim_x1 + input$sim_x2 + input$sim_x3 + input$sim_x4
      en_rango  <- suma_masa >= min(suma_masa_hald) && suma_masa <= max(suma_masa_hald)
      desc <- if (length(fuera_h) == 0) {
        "Suma porcentual de los 4 componentes químicos"
      } else {
        sprintf("Suma de los 4 compuestos; %s no influye en la predicción de este modelo, pero sigue formando parte de la mezcla",
                paste(fuera_h, collapse = ", "))
      }
      tarjeta_metrica("Balance de masa cerrada",
                      paste0(round(suma_masa, 1), "%"),
                      if (suma_masa > 100) "Supera el 100 %: no es posible"
                      else if (en_rango) sprintf("Dentro del rango observado (%g%% a %g%%)", min(suma_masa_hald), max(suma_masa_hald))
                      else sprintf("Fuera del rango observado (%g%% a %g%%)", min(suma_masa_hald), max(suma_masa_hald)),
                      desc,
                      nivel = if (en_rango) "success" else "danger")
    }
  })
  
  output$sim_tabla_comparativa_html <- renderUI({
    req(input$tab3_dataset)
    cl <- modelo_sim()
    if (input$tab3_dataset == "cars") {
      req(input$sim_minprice)
      nd <- nd_cars()
      claves <- names(etq_sim_cars)
      preds  <- vapply(claves, function(k) pred_sim_cars(k, nd), numeric(1))
      n_act  <- function(cv_obj) sum(extract_coef(cv_obj, "lambda.1se")[-1] != 0)
      modelos_lm <- list(cars_bic = mod_bic_cars, cars_aic = mod_aic_cars, cars_ols = mod_inicial_cars, cars_cp = mod_cp_cars)
      clave_rmse <- c(cars_bic = "BIC_Fase3", cars_ols = "OLS_Completo", cars_ridge = "Ridge", cars_lasso = "Lasso", cars_enet = "ElasticNet")
      
      vars_txt <- c(
        cars_bic   = sprintf("%d predictores (%s)", length(coef(mod_bic_cars)) - 1, paste(names(coef(mod_bic_cars))[-1], collapse = ", ")),
        cars_aic   = sprintf("%d predictores (%s)", length(coef(mod_aic_cars)) - 1, paste(names(coef(mod_aic_cars))[-1], collapse = ", ")),
        cars_ols   = sprintf("%d predictores (todas las covariables)", length(coef(mod_inicial_cars)) - 1),
        cars_cp    = sprintf("%d predictores (%s)", length(coef(mod_cp_cars)) - 1, paste(names(coef(mod_cp_cars))[-1], collapse = ", ")),
        cars_ridge = sprintf("%d de %d covariables, con coeficientes reducidos", n_act(cv_ridge_cars_all), ncol(d_cars) - 1),
        cars_lasso = sprintf("%d de %d covariables activas", n_act(cv_lasso_cars_all), ncol(d_cars) - 1),
        cars_enet  = sprintf("%d de %d covariables activas", n_act(cv_enet_cars_all), ncol(d_cars) - 1)
      )
      
      df <- data.frame(
        "Modelo evaluado" = unname(etq_sim_cars[claves]),
        "Variables consideradas" = unname(vars_txt[claves]),
        "Predicción" = sprintf("%.2f mpg", preds),
        "R² ajustado" = vapply(claves, function(k) {
          if (k %in% names(modelos_lm)) sprintf("%.2f%%", 100 * summary(modelos_lm[[k]])$adj.r.squared) else "no aplica"
        }, character(1)),
        "Error fuera de muestra" = vapply(claves, function(k) {
          if (k %in% names(clave_rmse)) sprintf("%.3f mpg", prom_rmse_cars[[clave_rmse[[k]]]]) else "no evaluado"
        }, character(1)),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      i_elegido <- match(clave_elegido_cars, claves)
      i_sel     <- match(cl, claves)
      insig <- rep("", length(claves))
      insig[i_sel] <- "Modelo seleccionado"
      insig[i_elegido] <- "Modelo elegido"
      tabla_ucm(df, fila_resaltada = unique(c(i_elegido, i_sel)), insignia = insig)
    } else {
      req(input$sim_x1)
      nd_h <- nd_hald()
      claves <- names(etq_sim_hald)
      preds  <- vapply(claves, function(k) pred_sim_hald(k, nd_h), numeric(1))
      
      vars_txt <- c(
        hald_lasso = if (length(lasso_cero_hald) > 0) sprintf("%s (anula %s)", paste(setdiff(names(x_hald), lasso_cero_hald), collapse = ", "), paste(lasso_cero_hald, collapse = ", "))
        else "X1, X2, X3, X4",
        hald_enet  = sprintf("%d de 4 componentes", activos_hald[[ "ENet" ]]),
        hald_ridge = sprintf("%d de 4 componentes", activos_hald[[ "Ridge" ]]),
        hald_1v    = vars_modelo_hald(mod_hald_1v),
        hald_2v    = sprintf("%s (Cp = %.2f, BIC = %.2f)", vars_modelo_hald(mod_hald_2v), resumen_hald$cp[2], resumen_hald$bic[2]),
        hald_3v    = sprintf("%s (R² ajustado = %.4f)", vars_modelo_hald(mod_hald_3v), resumen_hald$adjr2[3]),
        hald_fwd   = vars_modelo_hald(mod_forward_hald),
        hald_bwd   = vars_modelo_hald(mod_backward_hald),
        hald_step  = vars_modelo_hald(mod_stepwise_hald),
        hald_ols   = "X1, X2, X3, X4",
        hald_pcr   = sprintf("%d componentes que resumen X1, X2, X3 y X4", ncomp_opt)
      )
      interp <- c(
        hald_lasso = "Penalización L1: aplica esparcidad y reduce el número de componentes activos.",
        hald_enet  = "Penalización L1 + L2: distribuye el peso entre componentes correlacionados.",
        hald_ridge = "Penalización L2: contrae todos los coeficientes sin anular ninguno.",
        hald_1v    = sprintf("Modelo mínimo: explica menos que los demás (SSE = %.1f).", resumen_hald$rss[1]),
        hald_2v    = sprintf("VIF = %.3f (bajo). Evita la colinealidad, pero omite X4.", max(vif_hald_2v)),
        hald_3v    = sprintf("VIF máximo = %.2f (alto). Estimación inestable por la correlación entre X2 y X4.", max(vif_hald_3v)),
        hald_fwd   = sprintf("Selección automática según AIC: retiene %d variables.", length(coef(mod_forward_hald)) - 1),
        hald_bwd   = sprintf("Selección automática según AIC: retiene %d variables.", length(coef(mod_backward_hald)) - 1),
        hald_step  = sprintf("Selección automática según AIC: retiene %d variables.", length(coef(mod_stepwise_hald)) - 1),
        hald_ols   = sprintf("Matriz casi singular (C = %.2f); varianza inflada por la masa cerrada.", C_hald),
        hald_pcr   = "Resume las 4 variables en componentes independientes, lo que evita la colinealidad."
      )
      
      df <- data.frame(
        "Modelo evaluado" = unname(etq_sim_hald[claves]),
        "Variables utilizadas" = unname(vars_txt[claves]),
        "Calor predicho (y)" = sprintf("%.2f cal/g", preds),
        "Interpretación práctica" = unname(interp[claves]),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      i_elegido <- match(clave_elegido_hald, claves)
      i_sel     <- match(cl, claves)
      insig <- rep("", length(claves))
      insig[i_sel] <- "Modelo seleccionado"
      insig[i_elegido] <- "Modelo elegido"
      tabla_ucm(df, fila_resaltada = unique(c(i_elegido, i_sel)), insignia = insig)
    }
  })
  
  output$sim_grafico_distribucion <- renderPlotly({
    req(input$tab3_dataset)
    if (input$tab3_dataset == "cars") {
      req(input$sim_minprice)
      pred <- pred_cars_sel()
      p <- ggplot(d_cars, aes(x = MPG.city)) +
        geom_histogram(aes(y = after_stat(density)), bins = 16, fill = ucm$azul_claro, color = "#FFFFFF", alpha = 0.95) +
        geom_density(color = ucm$azul, linewidth = 1) +
        geom_vline(xintercept = pred, color = ucm$rojo, linewidth = 1.2, linetype = "dashed") +
        annotate("label", x = pred, y = Inf, label = paste0("Escenario: ", round(pred, 1), " mpg"),
                 vjust = 1.5, fill = "#FCEEEC", color = ucm$rojo, size = 3.3) +
        labs(x = "Rendimiento urbano (millas por galón)", y = "Densidad") +
        tema_graf(11)
      a_plotly(p)
    } else {
      req(input$sim_x1)
      pred_h <- pred_hald_sel()
      p <- ggplot(datos_hald, aes(x = y)) +
        geom_density(fill = ucm$azul_claro, color = ucm$azul, alpha = 0.7, linewidth = 1) +
        geom_point(aes(y = 0), size = 3, color = ucm$azul_medio, alpha = 0.8) +
        geom_vline(xintercept = pred_h, color = ucm$rojo, linewidth = 1.2, linetype = "dashed") +
        annotate("label", x = pred_h, y = Inf, label = paste0("Calor: ", round(pred_h, 1), " cal/g"),
                 vjust = 1.5, fill = "#FCEEEC", color = ucm$rojo, size = 3.3) +
        labs(x = "Calor liberado a 180 días (cal/g)", y = "Densidad") +
        tema_graf(11)
      a_plotly(p)
    }
  })
  
  # ========================================================================== #
  # PESTAÑA 3: DIAGNÓSTICO DE SUPUESTOS                                        #
  # ========================================================================== #
  
  opciones_diag <- function(ds) {
    if (ds == "cars") {
      c(setNames("cars_bic", sprintf("Cars93: modelo BIC (%d predictores)", length(coef(mod_bic_cars)) - 1)),
        setNames("cars_aic", sprintf("Cars93: modelo según AIC (%d predictores)", length(coef(mod_aic_cars)) - 1)),
        setNames("cars_ols", sprintf("Cars93: modelo inicial completo (%d predictores)", length(coef(mod_inicial_cars)) - 1)))
    } else {
      c("Hald: Lasso de 3 variables (modelo elegido por criterio experto)" = "hald_lasso",
        "Hald: modelo de 2 variables X1 + X2 (mínimo Cp y BIC)" = "hald_2v",
        "Hald: modelo de 3 variables X1 + X2 + X4 (máximo R² ajustado)" = "hald_3v",
        "Hald: modelo inicial completo (4 predictores)" = "hald_ols")
    }
  }
  
  selector_diag <- function(id, ds) {
    ops    <- opciones_diag(ds)
    actual <- isolate(input$tab3_sup_modelo)
    sel    <- if (!is.null(actual) && actual %in% ops) actual else if (ds == "cars") "cars_aic" else "hald_lasso"
    radioButtons(id, label = NULL, choices = ops, selected = sel, width = "100%")
  }
  
  output$ui_selector_modelo_supuestos <- renderUI({
    req(input$tab3_dataset)
    selector_diag("tab3_sup_modelo", input$tab3_dataset)
  })
  
  output$ui_selector_modelo_influencia <- renderUI({
    req(input$tab3_dataset)
    selector_diag("tab3_inf_modelo", input$tab3_dataset)
  })
  
  observeEvent(input$tab3_inf_modelo, {
    if (!identical(input$tab3_sup_modelo, input$tab3_inf_modelo))
      updateRadioButtons(session, "tab3_sup_modelo", selected = input$tab3_inf_modelo)
  }, ignoreInit = TRUE)
  
  observeEvent(input$tab3_sup_modelo, {
    if (!is.null(input$tab3_inf_modelo) && !identical(input$tab3_inf_modelo, input$tab3_sup_modelo))
      updateRadioButtons(session, "tab3_inf_modelo", selected = input$tab3_sup_modelo)
  }, ignoreInit = TRUE)
  
  modelo_auditado <- reactive({
    req(input$tab3_sup_modelo)
    switch(input$tab3_sup_modelo,
           "cars_bic" = mod_bic_cars,
           "cars_aic" = mod_aic_cars,
           "cars_ols" = mod_inicial_cars,
           "hald_2v"  = mod_hald_2v,
           "hald_3v"  = mod_hald_3v,
           "hald_lasso" = modelos_diag_hald[["hald_lasso"]],
           "hald_ols" = mod_completo_hald)
  })
  
  output$plot_sup_linealidad <- renderPlotly({
    m <- modelo_auditado()
    df <- data.frame(Ajustados = fitted(m), Residuos = resid(m))
    p <- ggplot(df, aes(x = Ajustados, y = Residuos)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = ucm$gris_texto) +
      geom_point(color = ucm$azul, size = 2.2, alpha = 0.8) +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = ucm$rojo, linewidth = 0.9) +
      labs(x = "Valores ajustados", y = "Residuos ordinarios") +
      tema_graf()
    a_plotly(p)
  })
  
  output$plot_sup_normalidad <- renderPlotly({
    m <- modelo_auditado()
    r_std <- rstandard(m)
    probs <- ppoints(length(r_std))
    df <- data.frame(Teoricos = qnorm(probs), Estandarizados = sort(r_std))
    p <- ggplot(df, aes(x = Teoricos, y = Estandarizados)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = ucm$rojo, linewidth = 0.8) +
      geom_point(color = ucm$azul_medio, size = 2.2, alpha = 0.85) +
      labs(x = "Cuantiles teóricos normales", y = "Residuos estandarizados") +
      tema_graf()
    a_plotly(p)
  })
  
  output$plot_sup_homocedasticidad <- renderPlotly({
    m <- modelo_auditado()
    r_std <- rstandard(m)
    df <- data.frame(Ajustados = fitted(m), RaizResid = sqrt(abs(r_std)))
    p <- ggplot(df, aes(x = Ajustados, y = RaizResid)) +
      geom_point(color = ucm$turquesa, size = 2.2, alpha = 0.8) +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = ucm$rojo, linewidth = 0.9) +
      labs(x = "Valores ajustados", y = "sqrt(|residuos estandarizados|)") +
      tema_graf()
    a_plotly(p)
  })
  
  output$plot_sup_leverage <- renderPlotly({
    m <- modelo_auditado()
    df <- data.frame(Leverage = hatvalues(m), Residuos = rstandard(m), Cook = cooks.distance(m))
    p <- ggplot(df, aes(x = Leverage, y = Residuos)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = ucm$gris_texto) +
      geom_point(aes(size = Cook), color = ucm$dorado_osc, alpha = 0.75) +
      labs(x = "Apalancamiento (hii)", y = "Residuos estandarizados") +
      tema_graf() + theme(legend.position = "none")
    a_plotly(p)
  })
  
  output$plot_diag_cook <- renderPlotly({
    m <- modelo_auditado()
    inf <- resumen_influencia(m, length(resid(m)))
    df <- data.frame(Obs = seq_along(inf$cook), Cook = as.numeric(inf$cook), Alerta = as.numeric(inf$cook) > inf$umbral_cook)
    p <- ggplot(df, aes(x = Obs, y = Cook, color = Alerta)) +
      geom_segment(aes(xend = Obs, yend = 0), linewidth = 0.5) +
      geom_point(size = 2) +
      geom_hline(yintercept = inf$umbral_cook, color = ucm$rojo, linetype = "dashed", linewidth = 0.7) +
      scale_color_manual(values = c("FALSE" = ucm$azul, "TRUE" = ucm$rojo)) +
      labs(x = "Índice de observación", y = "Distancia de Cook") +
      tema_graf() + theme(legend.position = "none")
    a_plotly(p)
  })
  
  output$plot_diag_hii <- renderPlotly({
    m <- modelo_auditado()
    inf <- resumen_influencia(m, length(resid(m)))
    df <- data.frame(Obs = seq_along(inf$hii), Hii = as.numeric(inf$hii), Alerta = as.numeric(inf$hii) > inf$umbral_hii)
    p <- ggplot(df, aes(x = Obs, y = Hii, color = Alerta)) +
      geom_segment(aes(xend = Obs, yend = 0), linewidth = 0.5) +
      geom_point(size = 2) +
      geom_hline(yintercept = inf$umbral_hii, color = ucm$rojo, linetype = "dashed", linewidth = 0.7) +
      scale_color_manual(values = c("FALSE" = ucm$azul_medio, "TRUE" = ucm$rojo)) +
      labs(x = "Índice de observación", y = "Apalancamiento (hii)") +
      tema_graf() + theme(legend.position = "none")
    a_plotly(p)
  })
  
  output$plot_diag_rstud <- renderPlotly({
    m <- modelo_auditado()
    rs <- rstudent(m)
    df <- data.frame(Obs = seq_along(rs), Rstud = as.numeric(rs), Alerta = abs(as.numeric(rs)) > 2)
    p <- ggplot(df, aes(x = Obs, y = Rstud, color = Alerta)) +
      geom_segment(aes(xend = Obs, yend = 0), linewidth = 0.5) +
      geom_point(size = 2) +
      geom_hline(yintercept = c(-2, 2), color = ucm$rojo, linetype = "dashed", linewidth = 0.7) +
      scale_color_manual(values = c("FALSE" = ucm$turquesa, "TRUE" = ucm$rojo)) +
      labs(x = "Índice de observación", y = "Residuo estudentizado") +
      tema_graf() + theme(legend.position = "none")
    a_plotly(p)
  })
  
  etq_sup <- c(
    cars_bic = "modelo BIC",
    cars_aic = "modelo AIC",
    cars_ols = "modelo inicial completo",
    hald_2v  = "modelo de 2 variables (X1, X2)",
    hald_3v  = "modelo de 3 variables (X1, X2, X4)",
    hald_ols = "modelo completo de 4 variables"
  )
  
  sens_sel <- reactive({
    req(input$tab3_sup_modelo)
    clave <- input$tab3_sup_modelo
    m     <- modelo_auditado()
    datos <- if (startsWith(clave, "cars")) d_cars else datos_hald
    inf   <- resumen_influencia(m, nrow(datos))
    res   <- list(clave = clave, idx = inf$idx_cook, umbral = inf$umbral_cook, df = NULL)
    if (length(inf$idx_cook) > 0) {
      m2 <- lm(formula(m), data = datos[-inf$idx_cook, ])
      b1 <- coef(m)
      b2 <- coef(m2)
      cambio <- abs(b1 - b2)
      res$df <- data.frame(
        "Predictor"       = names(b1),
        "Modelo original" = f4(b1),
        "Sin influyentes" = f4(b2),
        "Cambio absoluto" = f4(cambio),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      res$mayor_var    <- names(which.max(cambio[-1]))
      res$mayor_cambio <- max(cambio[-1])
    }
    res
  })
  
  output$tab3_texto_sensibilidad <- renderUI({
    sv <- sens_sel()
    estilo <- paste0("font-size:0.82rem; color:", ucm$gris_texto, ";")
    if (length(sv$idx) == 0) {
      p(class = "mb-0", style = estilo,
        sprintf("En el %s ningún caso supera la distancia de Cook crítica (4/n = %.4f), por lo que no hay observaciones que retirar y los coeficientes no cambian.",
                etq_sup[[sv$clave]], sv$umbral))
    } else {
      p(class = "mb-2", style = estilo,
        sprintf("Se reajusta el %s sin los %d casos que superan la distancia de Cook crítica (4/n = %.4f), observaciones %s:",
                etq_sup[[sv$clave]], length(sv$idx), sv$umbral, paste(sv$idx, collapse = ", ")))
    }
  })
  
  output$tab3_tabla_sensibilidad <- renderUI({
    sv <- sens_sel()
    req(!is.null(sv$df))
    tagList(
      tabla_ucm(sv$df),
      p(class = "mt-2 mb-0", style = paste0("font-size:0.82rem; color:", ucm$gris_texto, ";"),
        sprintf("El coeficiente que más se mueve es %s (cambia %s). Mientras mayor es el cambio, menos estable es la conclusión sobre esa variable frente a esos casos.",
                sv$mayor_var, f4(sv$mayor_cambio)))
    )
  })
  
  # ========================================================================== #
  # PESTAÑA 3: TRAZAS Y CURVAS DE REGULARIZACIÓN                               #
  # ========================================================================== #
  
  output$plot_reg_trazas <- renderPlotly({
    req(input$tab3_dataset, input$tab3_reg_metodo)
    ds  <- input$tab3_dataset
    met <- input$tab3_reg_metodo
    
    if (ds == "hald") {
      fit_obj <- switch(met, "lasso" = fit_lasso_hald, "ridge" = fit_ridge_hald, "enet" = fit_enet_hald)
      cv_obj  <- switch(met, "lasso" = cv_lasso_hald, "ridge" = cv_ridge_hald, "enet" = cv_enet_hald)
      vnames  <- c("X1 - aluminato", "X2 - silicato tricálcico", "X3 - ferrito", "X4 - silicato dicálcico")
    } else {
      cv_obj  <- switch(met, "lasso" = cv_lasso_cars_all, "ridge" = cv_ridge_cars_all, "enet" = cv_enet_cars_all)
      fit_obj <- cv_obj$glmnet.fit
      vnames  <- colnames(d_cars)[-1]
    }
    
    coef_m  <- as.matrix(fit_obj$beta)
    lambdas <- fit_obj$lambda
    
    df_trazas <- do.call(rbind, lapply(seq_len(nrow(coef_m)), function(i) {
      data.frame(
        log_lambda  = log(lambdas),
        Coeficiente = coef_m[i, ],
        Variable    = vnames[i],
        stringsAsFactors = FALSE
      )
    }))
    
    lmin <- log(cv_obj$lambda.min)
    l1se <- log(cv_obj$lambda.1se)
    
    p <- ggplot(df_trazas, aes(x = log_lambda, y = Coeficiente, color = Variable, group = Variable)) +
      geom_line(linewidth = 0.85, alpha = 0.85) +
      geom_vline(xintercept = lmin, color = ucm$verde, linetype = "dashed", linewidth = 0.7) +
      geom_vline(xintercept = l1se, color = ucm$violeta, linetype = "dotted", linewidth = 0.7) +
      labs(x = "log(lambda): nivel de penalización", y = "Coeficiente estandarizado") +
      tema_graf() + theme(legend.position = "right")
    
    a_plotly(p)
  })
  
  output$plot_reg_cv <- renderPlotly({
    req(input$tab3_dataset, input$tab3_reg_metodo)
    ds  <- input$tab3_dataset
    met <- input$tab3_reg_metodo
    
    cv_obj <- if (ds == "hald") {
      switch(met, "lasso" = cv_lasso_hald, "ridge" = cv_ridge_hald, "enet" = cv_enet_hald)
    } else {
      switch(met, "lasso" = cv_lasso_cars_all, "ridge" = cv_ridge_cars_all, "enet" = cv_enet_cars_all)
    }
    
    df_cv <- data.frame(
      log_lambda = log(cv_obj$lambda),
      MSE        = cv_obj$cvm,
      Upper      = cv_obj$cvup,
      Lower      = cv_obj$cvlo
    )
    
    lmin <- log(cv_obj$lambda.min)
    l1se <- log(cv_obj$lambda.1se)
    
    p <- ggplot(df_cv, aes(x = log_lambda)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = ucm$azul_claro, alpha = 0.6) +
      geom_line(aes(y = MSE), color = ucm$azul, linewidth = 0.9) +
      geom_vline(xintercept = lmin, color = ucm$verde, linetype = "dashed", linewidth = 0.8) +
      geom_vline(xintercept = l1se, color = ucm$violeta, linetype = "dotted", linewidth = 0.8) +
      labs(x = "log(lambda)", y = "Error cuadrático medio (MSE)") +
      tema_graf()
    
    a_plotly(p)
  })
  
  output$tab3_tabla_lambda_hald <- renderUI({
    tabla_ucm(tabla_lambda_hald)
  })
  
  output$plot_reg_box_rmse_cars <- renderPlotly({
    df_long <- tidyr::pivot_longer(
      as.data.frame(mat_rmse_cars),
      cols = everything(),
      names_to = "Clave",
      values_to = "RMSE"
    )
    df_long$Metodo <- factor(unname(etiquetas_rmse[df_long$Clave]), levels = unname(etiquetas_rmse))
    
    colores_box <- c(ucm$gris_texto, ucm$dorado, ucm$azul_medio, ucm$turquesa, ucm$violeta)
    names(colores_box) <- unname(etiquetas_rmse)
    
    p <- ggplot(df_long, aes(x = Metodo, y = RMSE, fill = Metodo)) +
      geom_boxplot(alpha = 0.85, outlier.size = 1.3) +
      scale_fill_manual(values = colores_box) +
      labs(x = "Método evaluado", y = "Error típico de predicción (mpg)") +
      tema_graf() +
      theme(axis.text.x = element_text(angle = 15, hjust = 1), legend.position = "none")
    
    a_plotly(p)
  })
  
  output$plot_reg_cv_hald_barras <- renderPlotly({
    cv_lista <- list("Ridge" = cv_ridge_hald, "Lasso" = cv_lasso_hald, "Elastic Net" = cv_enet_hald)
    df <- do.call(rbind, lapply(names(cv_lista), function(nm) {
      cv_obj <- cv_lista[[nm]]
      data.frame(
        Metodo = nm,
        Criterio = c("Menor error", "Modelo más simple"),
        RMSE = sqrt(c(cv_obj$cvm[cv_obj$index["min", ]], cv_obj$cvm[cv_obj$index["1se", ]])),
        stringsAsFactors = FALSE
      )
    }))
    df$Metodo <- factor(df$Metodo, levels = names(cv_lista))
    
    p <- ggplot(df, aes(x = Metodo, y = RMSE, fill = Criterio)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.9) +
      scale_fill_manual(values = c("Menor error" = ucm$azul, "Modelo más simple" = ucm$dorado)) +
      labs(x = "Método evaluado", y = "Error típico de predicción (cal/g)", fill = "Nivel de penalización") +
      tema_graf()
    
    a_plotly(p)
  })
}

# ---------------------------------------------------------------------------- #
# 5. Ejecución de la aplicación                                                 #
# ---------------------------------------------------------------------------- #

shinyApp(ui = ui, server = server)