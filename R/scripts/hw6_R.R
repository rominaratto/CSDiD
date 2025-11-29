#################################
# Question 1 — TWFE & Event-Study
#################################

# install.packages(c("sandwich", "lmtest"))

library(readr)   # para read_csv
library(sandwich)
library(lmtest)
library(ggplot2)
library(dplyr)
library(tidyr)

## a) TWFE

### Load the dataset available. Note that asmrs is the outcome variable, pcinc, asmrh and cases are controls. The dataset already includes treatment and post-treatment variables.
url <- "https://raw.githubusercontent.com/LOST-STATS/LOST-STATS.github.io/master/Model_Estimation/Data/Event_Study_DiD/bacon_example.csv"
df <- read_csv(url)
print(head(df))
print(dim(df))

### Estimate a Two-Way Fixed Effects (TWFE) regression with unit and time fixed effects.
twfe_lm <- lm(asmrs ~ post + pcinc + asmrh + cases + factor(stfips) + factor(year), data = df)
vcov_cl <- vcovCL(twfe_lm, cluster = ~ stfips)
coeftest(twfe_lm, vcov = vcov_cl)

## b) Cleaning for Event-Study

### Create the relative time variable (event time = time − treatment time).
df <- as.data.frame(df)
df_post <- df[df$post == 1, ]     # Filtrar solo filas donde post == 1
treat_time_by_state <- aggregate( # Sacar el primer año de tratamiento por estado
  year ~ stfips,                  # fórmula (year por stfips)
  data = df_post,     
  FUN = min)          
treat_vec <- setNames(treat_time_by_state$year, treat_time_by_state$stfips) # Crear vector con nombres (stfips → año)
df$treat_time <- treat_vec[as.character(df$stfips)]                         # Mapear treat_time a df
df$event_time <- df$year - df$treat_time                                    # Crear event_time
head(df[, c("stfips", "year", "post", "treat_time", "event_time")], 15)     # Verificar

### Provide a frequency table or descriptive summary of event time.
freq_full <- table(df$event_time, useNA = "ifany")
freq_full

event_nonmissing <- df$event_time[!is.na(df$event_time)] # Remover NA 
summary(event_nonmissing)

### Based on the distribution, choose reasonable upper and lower bounds for grouping extreme event times.
valid <- df$event_time[!is.na(df$event_time)]
ggplot(data.frame(valid), aes(x = valid)) +
  geom_histogram(binwidth = 1, color = "black", fill = "lightgray") +
  scale_x_continuous(breaks = seq(floor(min(valid)),
                                  ceiling(max(valid)),
                                  by = 2)) +
  labs(title = "Distribución de event_time",
       x = "event_time",
       y = "Frecuencia") +
  theme_minimal()

lower <- -10
upper <- 20
df$event_time_group <- df$event_time
df$event_time_group[df$event_time <= lower] <- lower
df$event_time_group[df$event_time >= upper] <- upper

### Question: Why do we usually group very distant event times together?
#### We group very distant event times because there are too few observations in those periods, 
#### which makes the coefficients noisy, unstable, hard to interpret, and irrelevant for policy. 
#### Binning the extremes improves precision, readability, and the reliability of the event-study.

### Create relative-time dummy variables.
df$event_time_group_int <- as.integer(df$event_time_group) # Asegurar que event_time_group es numérico
event_dummies <- model.matrix(                             # Crear dummies (model.matrix produce una matriz con columnas "event_…")
  ~ factor(event_time_group_int),   # formula
  data = df)
colnames(event_dummies) <- sub("factor\\(event_time_group_int\\)", "event", colnames(event_dummies)) # Renombrar columnas para que se vean como "event_-10", "event_0", etc.

## c) Event-Study Estimation

### Estimate an Event-Study model using the relative-time dummies.

### 1) Filtrar solo observaciones con event_time_group definido
df_es <- subset(df, !is.na(event_time_group_int))
### 2) Factor + periodo base = -1
df_es$event_time_group_int <- factor(df_es$event_time_group_int)
df_es$event_time_group_int <- relevel(df_es$event_time_group_int, ref = "-1")
### 3) Fórmula del Event-Study (NO necesitamos construir dummy_cols)
# Esto genera internamente las dummies para cada periodo relativo,
# tomando -1 como periodo base (igual que en tu código de Python).
es_formula <- asmrs ~ event_time_group_int + pcinc + asmrh + cases +
  factor(stfips) + factor(year)
cat("Fórmula usada:\n")
print(es_formula)
### 4) Estimar el modelo OLS con FE vía dummies
es_lm <- lm(es_formula, data = df_es)
### 5) Errores estándar clusterizados por stfips
vcov_es_cl <- vcovCL(es_lm, cluster = ~ stfips)
### 6) Resumen (coeficientes + SE clusterizados)
es_results <- coeftest(es_lm, vcov = vcov_es_cl)
es_results

### Store coefficient estimates and standard errors.
coefs <- coef(es_lm)
ses   <- sqrt(diag(vcov_es_cl))
event_idx  <- grep("^event_time_group_int", names(coefs))
event_coef <- coefs[event_idx]
event_se   <- ses[event_idx]
event_names <- names(event_coef)
event_time <- as.numeric(sub("event_time_group_int", "", event_names))
es_df <- data.frame(
  event_time = event_time,
  beta       = as.numeric(event_coef),
  se         = as.numeric(event_se))
es_df <- rbind(
  es_df,
  data.frame(event_time = -1, beta = 0, se = 0))
es_df <- es_df[order(es_df$event_time), ]
row.names(es_df) <- NULL
z <- 1.96
es_df$ci_low  <- es_df$beta - z * es_df$se
es_df$ci_high <- es_df$beta + z * es_df$se
es_df

### Plot the event-study coefficients with confidence intervals.
ggplot(es_df, aes(x = event_time, y = beta)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +   # IC 95% (errorbars)
  geom_line() +                                                      # Línea que une los coeficientes
  geom_point() +                                                     # Puntos
  geom_hline(yintercept = 0, linetype = "dashed") +                  # Línea horizontal en 0 (sin efecto)
  geom_vline(xintercept = -1, linetype = "dotted") +                 # Línea vertical en el periodo base (t = -1)
  labs(
    x = "Event time (años relativos al tratamiento)",
    y = "Efecto en asmrs (vs t = -1)",
    title = "Event-study: coeficientes e IC 95%"
  ) +
  theme_minimal()


#################################
# Question 2 — CSDiD
#################################

## a) Estimation

### Using CSDiD, estimate ATT(g,t). Present results in a clean table. (You may choose which controls to include.)

# G = cohorte de tratamiento:
# G = 0 para nunca tratados, treat_time para tratados
df$G <- ifelse(is.na(df$treat_time), 0, df$treat_time)
df$G <- as.integer(df$G)

cat("Valores únicos de G (cohortes):\n")
print(sort(unique(df$G)))

# Lista para ir guardando los resultados ATT(g,t)
rows <- list()
row_id <- 1

# Cohortes tratadas (excluimos G = 0)
cohortes_tratadas <- sort(unique(df$G[df$G != 0]))

for (g in cohortes_tratadas) {
  s <- g - 1   # año base (último pretratamiento)
  
  # Tratados en el año base s
  treated_base <- subset(df, G == g & year == s)
  if (nrow(treated_base) == 0) {
    next  # si no hay datos en (g,s), no se puede hacer DID 2x2
  }
  
  # Años en que aparece la cohorte g
  years_g <- sort(unique(df$year[df$G == g]))
  
  for (t in years_g) {
    if (t <= s) next  # solo años posteriores al base
    
    # ---- Tratados ----
    treated_t <- subset(df, G == g & year == t)
    if (nrow(treated_t) == 0) next
    
    # ---- Controles en t: nunca tratados + not-yet-treated en t ----
    control_t <- subset(
      df,
      year == t & (is.na(treat_time) | treat_time > t)
    )
    
    # ---- Controles en s: nunca tratados + tratados después de g ----
    control_base <- subset(
      df,
      year == s & (is.na(treat_time) | treat_time > g)
    )
    
    if (nrow(control_t) == 0 || nrow(control_base) == 0) next
    
    # Promedios de outcome
    mu_g_s <- mean(treated_base$asmrs, na.rm = TRUE)
    mu_g_t <- mean(treated_t$asmrs,   na.rm = TRUE)
    
    mu_c_s <- mean(control_base$asmrs, na.rm = TRUE)
    mu_c_t <- mean(control_t$asmrs,    na.rm = TRUE)
    
    # DID 2x2 estilo Callaway & Sant'Anna:
    # ATT(g,t) = (Y_g,t - Y_g,s) - (Y_c,t - Y_c,s)
    att_gt <- (mu_g_t - mu_g_s) - (mu_c_t - mu_c_s)
    
    rows[[row_id]] <- data.frame(
      g           = g,
      t           = t,
      ATT_gt      = att_gt,
      n_treated_s = nrow(treated_base),
      n_treated_t = nrow(treated_t),
      n_control_s = nrow(control_base),
      n_control_t = nrow(control_t)
    )
    row_id <- row_id + 1
  }
}

# Unir todo en una tabla
att_gt_table <- do.call(rbind, rows)

# Ordenar y redondear ATT
att_gt_table <- att_gt_table[order(att_gt_table$g, att_gt_table$t), ]
att_gt_table$ATT_gt <- round(att_gt_table$ATT_gt, 3)
row.names(att_gt_table) <- NULL

# Ver primeras filas (tabla limpia)
head(att_gt_table, 20)

## b) Aggregations
### Compute and present:
### Aggregation by group
att_by_group <- att_gt_table |>
  group_by(g) |>
  summarise(
    ATT_g = mean(ATT_gt, na.rm = TRUE),
    n_gt  = n(),
    .groups = "drop"
  ) |>
  mutate(ATT_g = round(ATT_g, 3)) |>
  arrange(g)
cat("Aggregation by group (cohort):\n")
print(att_by_group)

### Aggregation by period
att_by_period <- att_gt_table |>
  group_by(t) |>
  summarise(
    ATT_t = mean(ATT_gt, na.rm = TRUE),
    n_gt  = n(),          # cuántos (g,t) entran en el promedio
    .groups = "drop"
  ) |>
  mutate(ATT_t = round(ATT_t, 3)) |>
  arrange(t)
cat("Aggregation by period (calendar time):\n")
print(att_by_period)

### Aggregation by event-time
att_gt_table <- att_gt_table |>
  mutate(event_time = t - g)
att_by_event_time <- att_gt_table |>
  group_by(event_time) |>
  summarise(
    ATT_event = mean(ATT_gt, na.rm = TRUE),
    n_gt      = n(),
    .groups   = "drop"
  ) |>
  mutate(ATT_event = round(ATT_event, 3)) |>
  arrange(event_time)
cat("Aggregation by event-time (t - g):\n")
print(att_by_event_time)

## c) Explanation of Aggregations
### Explain:
### Meaning of aggregation by group
### Meaning of aggregation by period
### Meaning of aggregation by event-time
### Question: Which aggregation is most comparable to the TWFE Event-Study coefficients?

#### 1. Aggregation by Group (ATT(g))
#### Averages (ATT(g,t)) over all post-treatment periods for each cohort (g). Meaning: Measures how large the treatment effect was on average for each cohort across time.

#### 2. Aggregation by Period (ATT(t))
#### Averages (ATT(g,t)) across all cohorts already treated in calendar year (t). Meaning: Measures the average effect in a given year, pooling all treated cohorts.

#### 3. Aggregation by Event-Time (ATT(k))
#### Averages (ATT(g, g+k)) across cohorts by time relative to treatment ((k = t - g)). Meaning: Shows dynamic effects: how the treatment impact evolves (k) years before/after treatment.

#### Which aggregation matches TWFE Event-Study coefficients?
#### Event-time aggregation: TWFE event-study coefficients also estimate effects by time relative to treatment ((k)).

## d) Compare CSDiD and Event-Study Results
### Create a table comparing:(i) Event-time aggregated ATT from CSDiD & (ii) Event-study coefficients from the TWFE model
### Produce a combined coefficient plot.
### Provide a brief comparison.

# 1. TWFE event-study
twfe_event <- es_df |>
  select(
    event_time,
    beta_TWFE = beta,
    se_TWFE   = se)

# 2. CSDiD event-time aggregation
csdid_event <- att_by_event_time |>
  select(
    event_time,
    ATT_CSDiD      = ATT_event,
    n_cells_CSDiD  = n_gt)

# 3. MERGE: tabla conjunta TWFE + CSDiD por event_time
comparison <- twfe_event |>
  inner_join(csdid_event, by = "event_time") |>
  arrange(event_time) |>
  mutate(
    beta_TWFE  = round(beta_TWFE, 3),
    se_TWFE    = round(se_TWFE, 3),
    ATT_CSDiD  = round(ATT_CSDiD, 3))

cat("\n==============================\n")
cat("Comparación TWFE vs CSDiD\n")
cat("==============================\n")
print(comparison)

z <- 1.96                                                           # Parámetro: IC 95% TWFE
comparison <- comparison |>                                         # Calcular IC para TWFE
  mutate(
    ci_low_TWFE  = beta_TWFE - z * se_TWFE,
    ci_high_TWFE = beta_TWFE + z * se_TWFE)

ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Línea horizontal en 0
  
  # --- 1) TWFE event-study (puntos + barras de error) ---
  geom_errorbar(
    data = comparison,
    aes(x = event_time, ymin = ci_low_TWFE, ymax = ci_high_TWFE),
    width = 0.2,
    color = "#1f77b4"
  ) +
  geom_line(
    data = comparison,
    aes(x = event_time, y = beta_TWFE, color = "TWFE event-study")
  ) +
  geom_point(
    data = comparison,
    aes(x = event_time, y = beta_TWFE, color = "TWFE event-study")
  ) +
  
  # --- 2) CSDiD agregado (línea discontinua + cuadrados) ---
  geom_line(
    data = comparison,
    aes(x = event_time, y = ATT_CSDiD, color = "CSDiD (ATT event-time)"),
    linetype = "dashed"
  ) +
  geom_point(
    data = comparison,
    aes(x = event_time, y = ATT_CSDiD, color = "CSDiD (ATT event-time)"),
    shape = 15, # cuadradito
    size = 2.8
  ) +
  
  labs(
    x = "Event time (t - g)",
    y = "Efecto estimado (vs período base)",
    title = "Comparación: TWFE Event-Study vs CSDiD por event-time",
    color = "Modelo"
  ) +
  theme_minimal()

### Brief Comparison: TWFE exaggerates negative effects due to bad comparisons; CSDiD provides cleaner, more credible dynamic effects by using valid control groups at each g,t.
