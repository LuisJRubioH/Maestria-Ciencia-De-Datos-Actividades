# ============================================================
# ACTIVIDAD 3 — APLICACIÓN SHINY
# Modelación del Rendimiento Académico en Educación Secundaria
# Portuguesa mediante Regresión Lineal Simple
#
# Autor : Luis Javier Rubio Hernández
# Curso : Métodos y Simulación Estadística
# Fuente: Cortez & Silva (2008) — UCI Repository
# ============================================================
#
# Basado en el template del curso (Introduccion_a_Shiny_2.R),
# extendido con:
#   - Dos datasets reales (MAT y POR)
#   - Dos predictores seleccionables (G2 y failures)
#   - Exclusión interactiva de observaciones
#   - 6 pestañas: Datos · Correlaciones · Modelo ·
#                 Diagnóstico · Atípicos · Predicción
#
# Instalación de paquetes (solo primera vez):
# install.packages(c("shiny","ggplot2","DT","dplyr","corrplot",
#                    "lmtest","nortest","gridExtra"))
# ============================================================


# ------------------------------------------------------------
# 1. LIBRERÍAS
# ------------------------------------------------------------
library(shiny)
library(ggplot2)
library(DT)
library(dplyr)
library(corrplot)
library(lmtest)    # bptest, dwtest
library(nortest)   # lillie.test
library(gridExtra) # grid.arrange para paneles de diagnóstico


# ------------------------------------------------------------
# 2. CARGA Y PREPARACIÓN DE DATOS
# ------------------------------------------------------------
# Los CSV usan ";" como separador (Cortez & Silva, 2008)

mat_raw <- read.table("data/studentmat.csv", sep = ";",
                      header = TRUE, stringsAsFactors = FALSE)
por_raw <- read.table("data/studentpor.csv", sep = ";",
                      header = TRUE, stringsAsFactors = FALSE)

# Forzar conversión numérica de G1 y G2 (el CSV las trae entre comillas)
for (col in c("G1", "G2")) {
  mat_raw[[col]] <- suppressWarnings(as.numeric(as.character(mat_raw[[col]])))
  por_raw[[col]] <- suppressWarnings(as.numeric(as.character(por_raw[[col]])))
}

# Validación defensiva: aborta si la conversión introdujo NAs
if (anyNA(mat_raw$G1) || anyNA(mat_raw$G2))
  stop("Conversión de G1/G2 en studentmat.csv introdujo valores NA.")
if (anyNA(por_raw$G1) || anyNA(por_raw$G2))
  stop("Conversión de G1/G2 en studentpor.csv introdujo valores NA.")

# Excluir G3 = 0 (abandono tardío, no rendimiento ordinario)
# Justificación: Cortez & Silva (2008); ver Sección 1.3 del informe
mat_base <- mat_raw[mat_raw$G3 > 0, ]
por_base <- por_raw[por_raw$G3 > 0, ]

# Lista de datasets disponibles en la aplicación
lista_datos <- list(
  "Matemáticas (MAT)" = mat_base,
  "Portugués (POR)"   = por_base
)

# Variables cuantitativas candidatas (para la matriz de correlación)
vars_cuant <- c("age", "Medu", "Fedu", "traveltime", "studytime",
                "failures", "famrel", "freetime", "goout",
                "Dalc", "Walc", "health", "absences", "G1", "G2", "G3")

# Predictores disponibles en la interfaz (los dos mejores no-coincidentes)
# Rangos de correlación con G3 tras depurar G3 = 0 (verificado empíricamente):
#   G2:        r = 0.966 (MAT),  r = 0.934 (POR)
#   failures:  r = -0.294 (MAT), r = -0.388 (POR)
predictores_disponibles <- c(
  "G2       (r ≈ 0.93–0.97 con G3)"   = "G2",
  "failures (r ≈ -0.29 a -0.39 con G3)" = "failures"
)


# ------------------------------------------------------------
# 3. FUNCIONES AUXILIARES
# ------------------------------------------------------------

# Formatear p-valor con asteriscos de significancia
fmt_p <- function(p) {
  if (is.na(p)) return("—")
  sig <- dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
  paste0(formatC(p, format = "f", digits = 4), " ", sig)
}

# Umbrales de influencia (Belsley, Kuh y Welsch, 1980; p = 1 predictor)
#   Cook:    D_i > 4/n
#   Leverage:  h_ii > 2(p+1)/n  =  4/n  (con p = 1)
#   DFFITS:  |DFFITS_i| > 2 * sqrt((p+1)/n)  =  2 * sqrt(2/n)  (con p = 1)
umbral_cook  <- function(n) 4 / n
umbral_lev   <- function(n) 4 / n
umbral_dff   <- function(n) 2 * sqrt(2 / n)


# ------------------------------------------------------------
# 3.1 SUAVIZADO ADAPTATIVO
# ------------------------------------------------------------
# LOESS requiere variabilidad local; con predictores fuertemente discretos
# (p.ej. failures: 4 valores únicos, 80% concentrados en 0) las matrices
# locales se vuelven casi singulares y aparecen warnings tipo:
#   "reciprocal condition number 0", "pseudoinverse used", "near singularities".
# Estrategia: si el predictor tiene < 6 valores únicos, sustituir LOESS por
# una línea recta (geom_smooth method = "lm"). Para diagnósticos de residuales
# vs. ajustados, se aplica el mismo criterio de seguridad.
es_predictor_singular <- function(x, min_unicos = 6) {
  length(unique(x[!is.na(x)])) < min_unicos
}

# Capa ggplot adaptativa para overlay sobre dispersión
capa_suavizado <- function(x_vector, ...) {
  if (es_predictor_singular(x_vector)) {
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, ...)
  } else {
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, ...)
  }
}


# ------------------------------------------------------------
# 4. INTERFAZ DE USUARIO (UI)
# ------------------------------------------------------------
ui <- fluidPage(

  # ── Título ──────────────────────────────────────────────────
  titlePanel(
    div(
      h2("Regresión Lineal Simple — Rendimiento Académico"),
      h5(style = "color:#666;",
         "Cortez & Silva (2008) · Métodos y Simulación Estadística")
    )
  ),

  # ── Layout principal ─────────────────────────────────────────
  sidebarLayout(

    # ── Panel de control ───────────────────────────────────────
    sidebarPanel(
      width = 3,

      # Dataset (variable respuesta)
      selectInput(
        inputId  = "dataset",
        label    = tags$b("Asignatura (variable respuesta G3):"),
        choices  = names(lista_datos),
        selected = "Matemáticas (MAT)"
      ),

      # Predictor
      selectInput(
        inputId  = "predictor",
        label    = tags$b("Variable predictora (X):"),
        choices  = predictores_disponibles,
        selected = "G2"
      ),

      # Aviso cuando se selecciona failures (variable fuertemente discreta)
      conditionalPanel(
        condition = "input.predictor == 'failures'",
        div(
          style = "background:#fff3cd; padding:8px; border-radius:4px;
                   font-size:0.85em; color:#856404; margin-top:-5px;",
          icon("info-circle"),
          " 'failures' toma solo valores {0,1,2,3} y ~80% son 0. ",
          "El modelo es ilustrativo; interprétese con cautela."
        )
      ),

      hr(),

      # Nivel de confianza
      sliderInput(
        inputId = "nivel_conf",
        label   = tags$b("Nivel de confianza (1 - α):"),
        min = 0.80, max = 0.99, value = 0.95, step = 0.01
      ),

      hr(),

      # Exclusión de observaciones
      tags$b("Excluir observaciones del modelo:"),
      helpText("Seleccione índices para excluirlos en tiempo real."),
      uiOutput("ui_excluir"),

      hr(),

      # Resumen rápido del modelo
      tags$b("Modelo actual:"),
      verbatimTextOutput("modelo_texto"),

      hr(),
      helpText(
        icon("info-circle"),
        "Fuente: Cortez & Silva (2008).",
        tags$br(),
        "Notas: *** p<0.001, ** p<0.01, * p<0.05"
      )
    ),

    # ── Panel principal ────────────────────────────────────────
    mainPanel(
      width = 9,

      tabsetPanel(
        id = "pestanas",

        # ── PESTAÑA 1: DATOS ──────────────────────────────────
        tabPanel(
          title = "📋 Datos",
          br(),
          fluidRow(
            column(12,
              h4("Vista del conjunto de datos (observaciones activas)"),
              p("Las observaciones excluidas en el panel de control no aparecen en esta tabla."),
              DTOutput("tabla_datos")
            )
          )
        ),

        # ── PESTAÑA 2: CORRELACIONES ──────────────────────────
        tabPanel(
          title = "🔗 Correlaciones",
          br(),
          fluidRow(
            column(6,
              h4("Correlaciones de Pearson con G3"),
              p("Variables cuantitativas ordenadas por |r| descendente."),
              DTOutput("tabla_correlaciones")
            ),
            column(6,
              h4("Mapa de calor — variables cuantitativas"),
              plotOutput("plot_corrplot", height = "500px")
            )
          ),
          br(),
          fluidRow(
            column(12,
              h4("Diagrama de dispersión: G3 ~ predictor seleccionado"),
              plotOutput("plot_scatter", height = "400px")
            )
          )
        ),

        # ── PESTAÑA 3: MODELO ─────────────────────────────────
        tabPanel(
          title = "📈 Modelo",
          br(),
          fluidRow(
            column(6,
              h4("Ecuación ajustada"),
              verbatimTextOutput("ecuacion_modelo"),
              br(),
              h4("Tabla de coeficientes (prueba t)"),
              DTOutput("tabla_coeficientes"),
              br(),
              h4("Tabla ANOVA (prueba F)"),
              DTOutput("tabla_anova")
            ),
            column(6,
              h4("Métricas de bondad de ajuste"),
              DTOutput("tabla_r2"),
              br(),
              h4("Intervalos de confianza para parámetros (95%)"),
              DTOutput("tabla_ic_params"),
              br(),
              h4("Recta de regresión ajustada"),
              plotOutput("plot_regresion", height = "350px")
            )
          )
        ),

        # ── PESTAÑA 4: DIAGNÓSTICO ────────────────────────────
        tabPanel(
          title = "🔍 Diagnóstico",
          br(),
          fluidRow(
            column(8,
              h4("Gráficos de diagnóstico"),
              plotOutput("plot_diagnostico", height = "550px")
            ),
            column(4,
              h4("Pruebas formales de supuestos"),
              br(),
              DTOutput("tabla_supuestos"),
              br(),
              div(
                style = "background:#fff3cd; padding:10px; border-radius:5px;",
                tags$b("Nota metodológica:"),
                tags$p("Para n > 300, Shapiro-Wilk y Lilliefors tienen alta potencia
                        y pueden rechazar H₀ ante desviaciones triviales. El TCL
                        garantiza la validez asintótica de la inferencia sobre β₁
                        (Casella y Berger, 2002, Teo. 10.1.6); sin embargo, ",
                       tags$b("la cobertura nominal del 95% de los IP individuales
                        sí asume normalidad de ε"),
                       " y puede subcubrir bajo desviaciones marcadas (verificación
                        bootstrap en Anexo A10 del informe).")
              )
            )
          )
        ),

        # ── PESTAÑA 5: ATÍPICOS E INFLUYENTES ────────────────
        tabPanel(
          title = "⚠️ Atípicos e influyentes",
          br(),
          fluidRow(
            column(12,
              h4("Umbrales de detección"),
              uiOutput("ui_umbrales"),
              br()
            )
          ),
          fluidRow(
            column(7,
              h4("Distancia de Cook por observación"),
              plotOutput("plot_cook", height = "300px"),
              br(),
              h4("Leverage vs Residual estandarizado"),
              plotOutput("plot_bubble", height = "350px")
            ),
            column(5,
              h4("Observaciones influyentes detectadas"),
              DTOutput("tabla_influyentes"),
              br(),
              h4("Resumen por criterio"),
              DTOutput("tabla_flags")
            )
          )
        ),

        # ── PESTAÑA 6: PREDICCIÓN ─────────────────────────────
        tabPanel(
          title = "🎯 Predicción",
          br(),
          fluidRow(
            column(4,
              h4("Predicción puntual"),
              numericInput(
                inputId = "x_nuevo",
                label   = "Valor del predictor (X₀):",
                value   = 10, min = 0, max = 20, step = 1
              ),
              br(),
              DTOutput("tabla_prediccion_puntual"),
              br(),
              div(
                style = "background:#d4edda; padding:10px; border-radius:5px;",
                h5("Interpretación:"),
                uiOutput("texto_prediccion")
              )
            ),
            column(8,
              h4("Recta ajustada con bandas IC y IP"),
              plotOutput("plot_prediccion", height = "420px"),
              br(),
              h4("Tabla de predicciones en rango del predictor"),
              DTOutput("tabla_pred_rango")
            )
          )
        )

      ) # fin tabsetPanel
    )   # fin mainPanel
  )     # fin sidebarLayout
)       # fin fluidPage


# ------------------------------------------------------------
# 5. SERVIDOR
# ------------------------------------------------------------
server <- function(input, output, session) {

  # ── 5.1 DATOS REACTIVOS ──────────────────────────────────────
  datos_base <- reactive({
    lista_datos[[input$dataset]]
  })

  # Reset de índices excluidos cuando cambia el dataset
  # (evita arrastrar selecciones del dataset anterior con n distinto)
  observeEvent(input$dataset, {
    updateSelectInput(session, "excluir", selected = character(0))
  }, ignoreInit = TRUE)

  # UI de exclusión: se regenera si cambia el dataset
  output$ui_excluir <- renderUI({
    n <- nrow(datos_base())
    selectInput(
      inputId  = "excluir",
      label    = NULL,
      choices  = seq_len(n),
      selected = NULL,
      multiple = TRUE,
      selectize = TRUE
    )
  })

  # Datos filtrados (sin las observaciones excluidas)
  datos <- reactive({
    df <- datos_base()
    excl <- as.integer(input$excluir)
    # Filtrar índices válidos (defensivo ante índices obsoletos)
    excl <- excl[!is.na(excl) & excl >= 1 & excl <= nrow(df)]
    if (length(excl) > 0) df <- df[-excl, , drop = FALSE]
    df
  })

  # Ajuste dinámico del rango de x_nuevo según el predictor seleccionado
  # (G2 ∈ [5,19] vs failures ∈ [0,3] — evita extrapolaciones absurdas)
  observe({
    req(input$predictor)
    df <- datos_base()
    xv <- input$predictor
    xr <- range(df[[xv]], na.rm = TRUE)
    # Valor por defecto: mediana redondeada, dentro del rango observado
    val_def <- round(median(df[[xv]], na.rm = TRUE))
    updateNumericInput(
      session, "x_nuevo",
      value = val_def,
      min   = xr[1],
      max   = xr[2],
      step  = 1
    )
  })

  # ── 5.2 MODELO REACTIVO ──────────────────────────────────────
  modelo <- reactive({
    req(input$predictor)
    df  <- datos()
    frm <- as.formula(paste("G3 ~", input$predictor))
    lm(frm, data = df)
  })

  # ── 5.3 TEXTO RESUMEN EN SIDEBAR ─────────────────────────────
  output$modelo_texto <- renderPrint({
    m   <- modelo()
    b0  <- round(coef(m)[1], 3)
    b1  <- round(coef(m)[2], 3)
    r2  <- round(summary(m)$r.squared, 3)
    sig <- round(summary(m)$sigma, 3)
    cat(paste0(
      "G3 = ", b0,
      ifelse(b1 >= 0, " + ", " - "),
      abs(b1), " * ", input$predictor, "\n",
      "R² = ", r2, "   σ̂ = ", sig, "\n",
      "n  = ", nrow(datos())
    ))
  })

  # ── 5.4 PESTAÑA DATOS ────────────────────────────────────────
  output$tabla_datos <- renderDT({
    datatable(
      datos(),
      options   = list(pageLength = 10, scrollX = TRUE),
      rownames  = TRUE,
      class     = "stripe hover compact"
    )
  })

  # ── 5.5 PESTAÑA CORRELACIONES ─────────────────────────────────

  # Tabla de correlaciones con G3
  output$tabla_correlaciones <- renderDT({
    df <- datos()
    vars_disp <- intersect(vars_cuant, names(df))
    vars_pred <- setdiff(vars_disp, "G3")

    tab <- data.frame(
      Variable  = vars_pred,
      `r Pearson` = sapply(vars_pred, function(v)
        round(cor(df[[v]], df$G3, use = "complete.obs", method = "pearson"), 4)),
      `ρ Spearman` = sapply(vars_pred, function(v)
        round(cor(df[[v]], df$G3, use = "complete.obs", method = "spearman"), 4)),
      check.names = FALSE
    )
    tab <- tab[order(abs(tab$`r Pearson`), decreasing = TRUE), ]

    datatable(
      tab,
      rownames  = FALSE,
      options   = list(pageLength = 16, dom = "t"),
      class     = "stripe hover compact"
    ) |>
      formatStyle(
        "r Pearson",
        background = styleInterval(
          c(-0.5, -0.3, 0.3, 0.5),
          c("#f8d7da", "#ffeeba", "white", "#d4edda", "#c3e6cb")
        )
      )
  })

  # Mapa de calor de correlaciones
  output$plot_corrplot <- renderPlot({
    df      <- datos()
    vars_ok <- intersect(vars_cuant, names(df))
    mat_cor <- cor(df[, vars_ok], use = "complete.obs")

    corrplot::corrplot(
      mat_cor,
      method      = "color",
      type        = "lower",
      tl.col      = "#2C3E50",
      tl.srt      = 45,
      tl.cex      = 0.75,
      addCoef.col = "#333333",
      number.cex  = 0.55,
      col         = colorRampPalette(c("#C0392B","white","#2980B9"))(200),
      title       = paste("Correlaciones —", input$dataset),
      mar         = c(0, 0, 2, 0)
    )
  })

  # Dispersograma G3 ~ predictor
  output$plot_scatter <- renderPlot({
    df  <- datos()
    xv  <- input$predictor
    m   <- modelo()
    b0  <- round(coef(m)[1], 3)
    b1  <- round(coef(m)[2], 3)
    r2  <- round(summary(m)$r.squared, 3)

    # Banda IC al 95% sobre la recta (no LOESS) — válido para cualquier
    # predictor, evita las casi-singularidades con failures.
    suppressWarnings(
      ggplot(df, aes(x = .data[[xv]], y = G3)) +
        geom_jitter(width = 0.15, height = 0.15,
                    alpha = 0.45, size = 1.8, color = "#2980B9") +
        geom_smooth(method = "lm", formula = y ~ x,
                    se = TRUE, color = "#E74C3C",
                    fill = "#F1948A", alpha = 0.25, linewidth = 1.1) +
        annotate("text",
                 x = min(df[[xv]], na.rm = TRUE),
                 y = max(df$G3, na.rm = TRUE),
                 label = paste0("Ĝ3 = ", b0,
                                ifelse(b1 >= 0," + "," - "),
                                abs(b1), " · ", xv,
                                "\nR² = ", r2),
                 hjust = 0, vjust = 1, size = 4,
                 color = "#922B21", fontface = "bold") +
        labs(title = paste(input$dataset, "· G3 ~", xv),
             x = xv, y = "Calificación final (G3)") +
        theme_minimal(base_size = 13) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
    )
  })

  # ── 5.6 PESTAÑA MODELO ───────────────────────────────────────

  output$ecuacion_modelo <- renderPrint({
    m  <- modelo()
    b0 <- round(coef(m)[1], 4)
    b1 <- round(coef(m)[2], 4)
    cat(paste0(
      "Ĝ3 = ", b0,
      ifelse(b1 >= 0, " + ", " - "), abs(b1),
      " · ", input$predictor, "\n\n",
      "Interpretación de β̂₁:\n",
      "Por cada punto adicional en ", input$predictor, ",\n",
      "G3 cambia en promedio ", b1, " puntos."
    ))
  })

  output$tabla_coeficientes <- renderDT({
    m   <- modelo()
    sm  <- summary(m)
    co  <- coef(sm)
    tab <- data.frame(
      Parámetro    = rownames(co),
      Estimación   = round(co[, "Estimate"],   4),
      `Error Est.` = round(co[, "Std. Error"], 4),
      `t`          = round(co[, "t value"],    4),
      `p-valor`    = sapply(co[, "Pr(>|t|)"], fmt_p),
      check.names  = FALSE
    )
    # Solo aplicar formatStyle si hay matches de significancia,
    # evitando error de longitud 0 en styleEqual. Se usa unique() por
    # robustez ante duplicados (si ambos coeficientes comparten p-valor).
    matches_sig <- unique(grep("\\*", tab$`p-valor`, value = TRUE))
    dt <- datatable(tab, rownames = FALSE,
                    options = list(dom = "t", ordering = FALSE),
                    class = "stripe compact")
    if (length(matches_sig) > 0) {
      dt <- dt |> formatStyle(
        "p-valor",
        color = styleEqual(matches_sig,
                           rep("darkgreen", length(matches_sig)))
      )
    }
    dt
  })

  output$tabla_anova <- renderDT({
    m    <- modelo()
    av   <- anova(m)
    n    <- nrow(datos())
    SSR  <- round(av[1, "Sum Sq"],  4)
    SSE  <- round(av[2, "Sum Sq"],  4)
    SST  <- round(SSR + SSE,        4)
    MSR  <- round(av[1, "Mean Sq"], 4)
    MSE  <- round(av[2, "Mean Sq"], 4)
    Fv   <- round(av[1, "F value"], 4)
    pF   <- fmt_p(av[1, "Pr(>F)"])

    tab <- data.frame(
      `Fuente`  = c("Regresión", "Error", "Total"),
      SC        = c(SSR, SSE, SST),
      gl        = c(1, n - 2, n - 1),
      CM        = c(MSR, MSE, NA),
      F         = c(Fv,  NA,  NA),
      `p-valor` = c(pF,  "—", "—"),
      check.names = FALSE
    )
    datatable(tab, rownames = FALSE,
              options = list(dom = "t", ordering = FALSE),
              class = "stripe compact")
  })

  output$tabla_r2 <- renderDT({
    sm  <- summary(modelo())
    tab <- data.frame(
      Métrica = c("R²", "R² ajustado",
                  "σ̂ (Error est. residual)", "n"),
      Valor   = c(round(sm$r.squared,     4),
                  round(sm$adj.r.squared, 4),
                  round(sm$sigma,         4),
                  nrow(datos()))
    )
    datatable(tab, rownames = FALSE,
              options = list(dom = "t", ordering = FALSE),
              class = "stripe compact")
  })

  output$tabla_ic_params <- renderDT({
    m   <- modelo()
    nc  <- input$nivel_conf
    ic  <- confint(m, level = nc)
    co  <- coef(m)
    pct <- c((1 - nc)/2, 1 - (1 - nc)/2) * 100

    tab <- data.frame(
      Parámetro = rownames(ic),
      Estimación = round(co, 4),
      `LC inf.`  = round(ic[, 1], 4),
      `LC sup.`  = round(ic[, 2], 4),
      check.names = FALSE
    )
    names(tab)[3:4] <- paste0(c("LC inf.", "LC sup."),
                              " (", round(pct, 1), "%)")
    datatable(tab, rownames = FALSE,
              options = list(dom = "t", ordering = FALSE),
              class = "stripe compact")
  })

  output$plot_regresion <- renderPlot({
    df  <- datos()
    m   <- modelo()
    xv  <- input$predictor
    nc  <- input$nivel_conf
    xr  <- range(df[[xv]], na.rm = TRUE)
    xs  <- data.frame(x = seq(xr[1], xr[2], length.out = 200))
    names(xs) <- xv
    ic  <- as.data.frame(predict(m, newdata = xs,
                                 interval = "confidence", level = nc))
    bnd <- cbind(xs, ic)

    ggplot() +
      geom_ribbon(data = bnd,
                  aes(x = .data[[xv]], ymin = lwr, ymax = upr),
                  fill = "#2980B9", alpha = 0.25) +
      geom_jitter(data = df,
                  aes(x = .data[[xv]], y = G3),
                  width = 0.12, height = 0.12,
                  alpha = 0.40, size = 1.5, color = "#1A5276") +
      geom_line(data = bnd,
                aes(x = .data[[xv]], y = fit),
                color = "#2980B9", linewidth = 1.2) +
      labs(title = paste0("G3 ~ ", xv, "  [IC ", round(nc*100), "%]"),
           x = xv, y = "G3") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
  })

  # ── 5.7 PESTAÑA DIAGNÓSTICO ──────────────────────────────────

  output$plot_diagnostico <- renderPlot({
    m   <- modelo()
    df  <- datos()
    xv  <- input$predictor
    res <- residuals(m)
    fit <- fitted(m)
    rst <- rstandard(m)

    df_d <- data.frame(
      idx         = seq_along(res),
      resid       = res,
      fitted      = fit,
      rst         = rst,
      sqrt_abs_r  = sqrt(abs(rst))
    )

    # Si los valores ajustados tienen pocos valores únicos (predictor muy
    # discreto, p.ej. failures), LOESS no es estable. Usar lm en su lugar.
    metodo_suave <- if (es_predictor_singular(df_d$fitted)) "lm" else "loess"

    p1 <- suppressWarnings(
      ggplot(df_d, aes(x = fitted, y = resid)) +
        geom_point(alpha = 0.45, size = 1.6, color = "#2980B9") +
        geom_hline(yintercept = 0, linetype = "dashed",
                   color = "#C0392B", linewidth = 0.8) +
        geom_smooth(method = metodo_suave, formula = y ~ x,
                    color = "#C0392B", se = FALSE, linewidth = 0.9) +
        labs(title = "(S1) Residuales vs Ajustados",
             x = "Valores ajustados", y = "Residuales") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
    )

    p2 <- ggplot(df_d, aes(sample = rst)) +
      stat_qq(alpha = 0.45, size = 1.6, color = "#2980B9") +
      stat_qq_line(color = "#C0392B", linewidth = 0.9, linetype = "dashed") +
      labs(title = "(S2) Q-Q Normal",
           x = "Cuantiles teóricos", y = "Residuales estand.") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))

    p3 <- suppressWarnings(
      ggplot(df_d, aes(x = fitted, y = sqrt_abs_r)) +
        geom_point(alpha = 0.45, size = 1.6, color = "#2980B9") +
        geom_smooth(method = metodo_suave, formula = y ~ x,
                    color = "#C0392B", se = FALSE, linewidth = 0.9) +
        labs(title = "(S3) Scale-Location",
             x = "Valores ajustados",
             y = expression(sqrt("|Res. estand.|"))) +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
    )

    p4 <- ggplot(df_d, aes(x = idx, y = resid)) +
      geom_line(alpha = 0.50, color = "#2980B9", linewidth = 0.6) +
      geom_point(alpha = 0.35, size = 1.2, color = "#2980B9") +
      geom_hline(yintercept = 0, linetype = "dashed",
                 color = "#C0392B", linewidth = 0.8) +
      labs(title = "(S4) Residuales vs Índice",
           x = "Índice de observación", y = "Residuales") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))

    gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
  })

  output$tabla_supuestos <- renderDT({
    m   <- modelo()
    res <- residuals(m)
    n   <- length(res)

    # (S1) Linealidad — RESET de Ramsey (1969)
    # Para predictores con muy pocos valores únicos (failures), el RESET
    # puede no estar bien definido (falta variación en fitted^2, fitted^3).
    # Se omite con mensaje informativo en ese caso.
    rs  <- if (es_predictor_singular(fitted(m))) {
      NULL
    } else {
      tryCatch(lmtest::resettest(m, power = 2:3, type = "fitted"),
               error = function(e) NULL)
    }
    sw  <- tryCatch(shapiro.test(res),  error = function(e) NULL)
    lf  <- tryCatch(nortest::lillie.test(res), error = function(e) NULL)
    bp  <- tryCatch(lmtest::bptest(m),  error = function(e) NULL)
    dw  <- tryCatch(lmtest::dwtest(m, alternative = "two.sided"),
                    error = function(e) NULL)

    tab <- data.frame(
      Supuesto  = c("(S1) Linealidad",
                    "(S2) Normalidad", "(S2) Normalidad",
                    "(S3) Homocedasticidad", "(S4) Independencia"),
      Prueba    = c("RESET (Ramsey)",
                    "Shapiro-Wilk", "Lilliefors",
                    "Breusch-Pagan", "Durbin-Watson"),
      Estadístico = c(
        if (!is.null(rs)) round(rs$statistic, 4)  else "—",
        if (!is.null(sw)) round(sw$statistic, 4)  else NA,
        if (!is.null(lf)) round(lf$statistic, 4)  else NA,
        if (!is.null(bp)) round(bp$statistic, 4)  else NA,
        if (!is.null(dw)) round(dw$statistic, 4)  else NA
      ),
      `p-valor` = c(
        if (!is.null(rs)) fmt_p(rs$p.value) else "n/a",
        if (!is.null(sw)) fmt_p(sw$p.value) else "—",
        if (!is.null(lf)) fmt_p(lf$p.value) else "—",
        if (!is.null(bp)) fmt_p(bp$p.value) else "—",
        if (!is.null(dw)) fmt_p(dw$p.value) else "—"
      ),
      `Decisión (α=0.05)` = c(
        if (!is.null(rs))
          ifelse(rs$p.value < 0.05, "Rechaza H₀", "No rechaza H₀")
        else "n/a (predictor con pocos niveles)",
        if (!is.null(sw))
          ifelse(sw$p.value < 0.05, "Rechaza H₀", "No rechaza H₀") else "—",
        if (!is.null(lf))
          ifelse(lf$p.value < 0.05, "Rechaza H₀", "No rechaza H₀") else "—",
        if (!is.null(bp))
          ifelse(bp$p.value < 0.05, "Rechaza H₀", "No rechaza H₀") else "—",
        if (!is.null(dw))
          ifelse(dw$p.value < 0.05, "Rechaza H₀", "No rechaza H₀") else "—"
      ),
      check.names = FALSE
    )

    datatable(tab, rownames = FALSE,
              options = list(dom = "t", ordering = FALSE),
              class = "stripe compact")
  })

  # ── 5.8 PESTAÑA ATÍPICOS E INFLUYENTES ───────────────────────

  output$ui_umbrales <- renderUI({
    n  <- nrow(datos())
    div(
      style = "background:#f8f9fa; padding:10px; border-radius:5px;",
      fluidRow(
        column(3, tags$b("Cook D_i >"), br(), round(umbral_cook(n), 4)),
        column(3, tags$b("Leverage h_ii >"), br(), round(umbral_lev(n), 4)),
        column(3, tags$b("|r_i| >"), br(), "2.0000"),
        column(3, tags$b("|DFFITS| >"), br(), round(umbral_dff(n), 4))
      )
    )
  })

  # Indicadores de influencia reactivos
  infl_data <- reactive({
    m   <- modelo()
    df  <- datos()
    n   <- nrow(df)

    cook  <- cooks.distance(m)
    lev   <- hatvalues(m)
    rstd  <- rstandard(m)
    dffit <- dffits(m)

    flag_c <- cook  > umbral_cook(n)
    flag_l <- lev   > umbral_lev(n)
    flag_r <- abs(rstd)  > 2
    flag_d <- abs(dffit) > umbral_dff(n)
    flag_a <- flag_c | flag_l | flag_r | flag_d

    list(
      cook = cook, lev = lev, rstd = rstd, dffit = dffit,
      flag_c = flag_c, flag_l = flag_l,
      flag_r = flag_r, flag_d = flag_d, flag_a = flag_a,
      n = n
    )
  })

  output$plot_cook <- renderPlot({
    id  <- infl_data()
    n   <- id$n
    df_c <- data.frame(
      idx  = seq_along(id$cook),
      cook = id$cook,
      flag = id$flag_a
    )
    ggplot(df_c, aes(x = idx, y = cook, fill = flag)) +
      geom_col(width = 0.7, alpha = 0.80) +
      geom_hline(yintercept = umbral_cook(n),
                 color = "#C0392B", linetype = "dashed", linewidth = 0.9) +
      geom_text(data = df_c[df_c$flag, ],
                aes(label = idx), vjust = -0.4,
                size = 2.8, color = "#C0392B") +
      scale_fill_manual(values = c("FALSE" = "#85C1E9",
                                   "TRUE"  = "#C0392B"),
                        guide = "none") +
      labs(title = "Distancia de Cook por observación",
           x = "Índice", y = expression(D[i])) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
  })

  output$plot_bubble <- renderPlot({
    id  <- infl_data()
    n   <- id$n
    df_b <- data.frame(
      lev  = id$lev,
      rstd = id$rstd,
      cook = id$cook,
      flag = id$flag_a,
      idx  = seq_along(id$lev)
    )
    ggplot(df_b, aes(x = lev, y = rstd,
                     size = cook, color = flag)) +
      geom_point(alpha = 0.65) +
      geom_hline(yintercept = c(-2, 2),
                 color = "#E67E22", linetype = "dashed") +
      geom_vline(xintercept = umbral_lev(n),
                 color = "#8E44AD", linetype = "dashed") +
      geom_text(data = df_b[df_b$flag, ],
                aes(label = idx), size = 2.8,
                vjust = -0.8, color = "#C0392B") +
      scale_color_manual(values = c("FALSE" = "#5DADE2",
                                    "TRUE"  = "#C0392B"),
                         guide = "none") +
      scale_size_continuous(name = expression(D[i]),
                            range = c(1, 8)) +
      labs(title = "Influencia: leverage vs residual estandarizado",
           x = expression(h[ii]~"(leverage)"),
           y = expression(r[i]~"(residual estandarizado)")) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
  })

  output$tabla_influyentes <- renderDT({
    id  <- infl_data()
    df  <- datos()
    xv  <- input$predictor
    idx <- which(id$flag_a)

    if (length(idx) == 0) {
      return(datatable(
        data.frame(Resultado = "No se detectaron observaciones influyentes."),
        rownames = FALSE, options = list(dom = "t")
      ))
    }

    tab <- data.frame(
      Obs           = idx,
      G2_o_pred     = df[[xv]][idx],
      G3            = df$G3[idx],
      `r_i`         = round(id$rstd[idx],  3),
      `h_ii`        = round(id$lev[idx],   4),
      `D_i (Cook)`  = round(id$cook[idx],  5),
      `DFFITS`      = round(id$dffit[idx], 3),
      `Cook ✓`      = ifelse(id$flag_c[idx], "✓", ""),
      `Lev ✓`       = ifelse(id$flag_l[idx], "✓", ""),
      `|r|>2 ✓`     = ifelse(id$flag_r[idx], "✓", ""),
      `DFFITS ✓`    = ifelse(id$flag_d[idx], "✓", ""),
      check.names   = FALSE
    )
    names(tab)[2] <- xv

    datatable(tab, rownames = FALSE,
              options = list(pageLength = 8, scrollX = TRUE,
                             dom = "tip"),
              class = "stripe compact")
  })

  output$tabla_flags <- renderDT({
    id <- infl_data()
    n  <- id$n
    tab <- data.frame(
      Criterio    = c("|r_i| > 2", "Leverage", "Cook", "DFFITS", "≥ 1 criterio"),
      `N° obs.`   = c(sum(id$flag_r), sum(id$flag_l),
                      sum(id$flag_c), sum(id$flag_d), sum(id$flag_a)),
      `% total`   = round(c(sum(id$flag_r), sum(id$flag_l),
                             sum(id$flag_c), sum(id$flag_d),
                             sum(id$flag_a)) / n * 100, 2),
      check.names = FALSE
    )
    datatable(tab, rownames = FALSE,
              options = list(dom = "t", ordering = FALSE),
              class = "stripe compact")
  })

  # ── 5.9 PESTAÑA PREDICCIÓN ───────────────────────────────────

  output$tabla_prediccion_puntual <- renderDT({
    m   <- modelo()
    xv  <- input$predictor
    x0  <- input$x_nuevo
    nc  <- input$nivel_conf
    nd  <- setNames(data.frame(x0), xv)

    ic  <- predict(m, newdata = nd, interval = "confidence", level = nc)
    ip  <- predict(m, newdata = nd, interval = "prediction", level = nc)

    tab <- data.frame(
      Tipo        = c("Ĝ3 (predicción puntual)",
                      paste0("IC media (", round(nc*100), "%)"),
                      paste0("IP individual (", round(nc*100), "%)")),
      Inferior    = c("—",
                      round(ic[,"lwr"], 3),
                      round(ip[,"lwr"], 3)),
      Centro      = c(round(ic[,"fit"], 3), "—", "—"),
      Superior    = c("—",
                      round(ic[,"upr"], 3),
                      round(ip[,"upr"], 3)),
      Amplitud    = c("—",
                      round(ic[,"upr"] - ic[,"lwr"], 3),
                      round(ip[,"upr"] - ip[,"lwr"], 3)),
      check.names = FALSE
    )
    datatable(tab, rownames = FALSE,
              options = list(dom = "t", ordering = FALSE),
              class = "stripe compact")
  })

  output$texto_prediccion <- renderUI({
    m   <- modelo()
    df  <- datos()
    xv  <- input$predictor
    x0  <- input$x_nuevo
    nc  <- input$nivel_conf
    nd  <- setNames(data.frame(x0), xv)
    ic  <- predict(m, newdata = nd, interval = "confidence", level = nc)
    ip  <- predict(m, newdata = nd, interval = "prediction", level = nc)
    xr  <- range(df[[xv]], na.rm = TRUE)
    extrapolando <- (x0 < xr[1]) || (x0 > xr[2])

    fuera_cota <- (xv == "G2" || xv == "G1") &&
                  (ip[,"upr"] > 20 || ip[,"lwr"] < 0)

    tagList(
      tags$p(
        paste0(
          "Para un estudiante con ", xv, " = ", x0,
          ", se predice una calificación final de ",
          round(ic[,"fit"], 2), " puntos. ",
          "Con un ", round(nc*100), "% de confianza, ",
          tags$b("la calificación media"),
          " de todos los estudiantes con este valor se ubica entre ",
          round(ic[,"lwr"], 2), " y ", round(ic[,"upr"], 2),
          "; y ", tags$b("la calificación de un estudiante específico"),
          " entre ", round(ip[,"lwr"], 2), " y ", round(ip[,"upr"], 2), "."
        )
      ),
      if (extrapolando) {
        tags$p(
          style = "color:#922B21; font-size:0.9em; margin-top:8px;",
          icon("exclamation-triangle"),
          " ", tags$b("Advertencia:"),
          " x₀ = ", x0, " está fuera del rango muestral [",
          xr[1], ", ", xr[2], "]. La predicción es una ",
          tags$b("extrapolación"), " y debe interpretarse con cautela."
        )
      },
      if (fuera_cota) {
        tags$p(
          style = "color:#7E5109; font-size:0.9em; margin-top:8px;",
          icon("info-circle"),
          " El IP rebasa la cota física G3 ∈ [0, 20]; artefacto de la
            regresión lineal sobre variable acotada (alternativa: regresión
            beta de Ferrari y Cribari-Neto, 2004)."
        )
      }
    )
  })

  output$plot_prediccion <- renderPlot({
    m   <- modelo()
    df  <- datos()
    xv  <- input$predictor
    x0  <- input$x_nuevo
    nc  <- input$nivel_conf
    xr  <- range(df[[xv]], na.rm = TRUE)
    xs  <- setNames(data.frame(
      seq(xr[1], xr[2], length.out = 300)), xv)

    ic  <- as.data.frame(predict(m, newdata = xs,
                                 interval = "confidence", level = nc))
    ip  <- as.data.frame(predict(m, newdata = xs,
                                 interval = "prediction", level = nc))
    bnd <- cbind(xs, fit = ic$fit,
                 ic_lwr = ic$lwr, ic_upr = ic$upr,
                 ip_lwr = ip$lwr, ip_upr = ip$upr)

    # Predicción en x0
    nd0 <- setNames(data.frame(x0), xv)
    pred0 <- predict(m, newdata = nd0)

    ggplot() +
      geom_ribbon(data = bnd,
                  aes(x = .data[[xv]], ymin = ip_lwr, ymax = ip_upr),
                  fill = "#AED6F1", alpha = 0.30) +
      geom_ribbon(data = bnd,
                  aes(x = .data[[xv]], ymin = ic_lwr, ymax = ic_upr),
                  fill = "#2980B9", alpha = 0.30) +
      geom_jitter(data = df,
                  aes(x = .data[[xv]], y = G3),
                  width = 0.12, height = 0.12,
                  alpha = 0.35, size = 1.4, color = "#1A5276") +
      geom_line(data = bnd,
                aes(x = .data[[xv]], y = fit),
                color = "#2980B9", linewidth = 1.2) +
      # Punto de predicción
      geom_vline(xintercept = x0,
                 color = "#E74C3C", linetype = "dashed", linewidth = 0.9) +
      geom_point(aes(x = x0, y = pred0),
                 color = "#E74C3C", size = 4, shape = 18) +
      annotate("text", x = x0 + 0.3, y = pred0 + 0.5,
               label = paste0("Ĝ3 = ", round(pred0, 2)),
               color = "#C0392B", size = 4, fontface = "bold") +
      labs(
        title = paste0(input$dataset, " — Recta ajustada con IC y IP (",
                       round(nc*100), "%)"),
        subtitle = paste0("Banda oscura: IC para E(G3|X=x₀)  ·  ",
                          "Banda clara: IP para obs. individual"),
        x = xv, y = "Calificación final (G3)"
      ) +
      theme_minimal(base_size = 12) +
      theme(plot.title    = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(color = "#555", hjust = 0.5))
  })

  output$tabla_pred_rango <- renderDT({
    m  <- modelo()
    xv <- input$predictor
    df <- datos()
    nc <- input$nivel_conf
    xr <- range(df[[xv]], na.rm = TRUE)
    x0s <- setNames(data.frame(
      seq(ceiling(xr[1]), floor(xr[2]), by = 1)), xv)

    ic  <- as.data.frame(predict(m, newdata = x0s,
                                 interval = "confidence", level = nc))
    ip  <- as.data.frame(predict(m, newdata = x0s,
                                 interval = "prediction", level = nc))

    tab <- data.frame(
      X0             = x0s[[1]],
      `Ĝ3`           = round(ic$fit, 3),
      `IC inf.`      = round(ic$lwr, 3),
      `IC sup.`      = round(ic$upr, 3),
      `IP inf.`      = round(ip$lwr, 3),
      `IP sup.`      = round(ip$upr, 3),
      check.names    = FALSE
    )
    names(tab)[1] <- xv

    datatable(tab, rownames = FALSE,
              options = list(pageLength = 10, dom = "tip"),
              class   = "stripe hover compact") |>
      formatRound(columns = 2:6, digits = 3)
  })

} # fin server


# ------------------------------------------------------------
# 6. EJECUTAR LA APLICACIÓN
# ------------------------------------------------------------
shinyApp(ui = ui, server = server)
