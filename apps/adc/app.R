# ============================================================================
# Antibody-drug conjugate: T-DM1 from mouse efficacy to patients
# (Singh & Shah, AAPS J 2017)
#
# Plasma PK of total trastuzumab, T-DM1 and released DM1; tumour disposition
# of the ADC (vascular and surface exchange, HER2 binding, internalisation,
# lysosomal release of DM1, tubulin binding); intracellular DM1 drives
# tumour-cell killing. Mouse tumour models give the killing constants; human
# PK is scaled from monkeys.
# ============================================================================

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(scales)

shared <- function(f) source(if (file.exists(file.path("R", f))) file.path("R", f) else file.path("..", "..", "shared", f))
for (f in c("theme.R", "ode_engine.R", "app_helpers.R")) shared(f)
source(file.path("R", "tdm1.R"))

M <- load_model("adc_tdm1_singh2017.cpp")
TOL <- list(rtol = 1e-5, atol = 1e-9)

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = tags$div(tags$span("QSP", style = "font-weight: bold; font-size: 24px;"),
                                   tags$span(" Antibody-drug conjugate", style = "font-size: 18px;")),
                  titleWidth = 350),
  dashboardSidebar(
    width = 320, custom_css,
    sidebarMenu(
      id = "tabs",
      menuItem("Mouse efficacy", tabName = "mouse", icon = icon("flask")),
      menuItem("Patients: exposure", tabName = "human", icon = icon("user")),
      menuItem("Patients: tumour course", tabName = "course", icon = icon("chart-line")),
      menuItem("Model setup", tabName = "setup", icon = icon("sliders")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    tags$div(
      style = "padding: 18px; padding-top: 8px;",
      conditionalPanel(
        "input.tabs == 'mouse'",
        tags$h4("Mouse tumour model"),
        selectInput("mmodel", NULL, choices = MOUSE_MODELS$model, selected = "KPL-4", width = "100%"),
        uiOutput("mmodel_note"),
        tags$h4("T-DM1 regimen (IV)"),
        fluidRow(column(6, numericInput("mdose", "Dose (mg/kg)", 3, min = 0, max = 60, step = 0.5)),
                 column(6, numericInput("mn", "Doses", 1, min = 1, max = 12, step = 1))),
        fluidRow(column(6, numericInput("mii", "Every (days)", 21, min = 1, max = 42, step = 1)),
                 column(6, numericInput("mtv0", "Start (mm3)", 200, min = 20, max = 2000, step = 20))),
        numericInput("mdays", "Follow (days)", 63, min = 14, max = 180, step = 7, width = "100%")
      ),
      conditionalPanel(
        "input.tabs != 'mouse'",
        tags$h4("Patient regimen"),
        selectInput("reg", NULL, choices = names(REGIMENS), selected = names(REGIMENS)[1], width = "100%"),
        numericInput("cycles", "Cycles (3 weeks each)", 6, min = 1, max = 17, step = 1, width = "100%"),
        radioButtons("her2", "Tumour HER2 (IHC)", c("3+" = "3+", "2+" = "2+", "1+" = "1+"), inline = TRUE),
        numericInput("bw", "Body weight (kg)", 70, min = 35, max = 150, step = 5, width = "100%")
      )
    )
  ),
  dashboardBody(
    custom_css,
    tabItems(
      tabItem(
        tabName = "mouse",
        fluidRow(valueBoxOutput("m_tv", width = 4), valueBoxOutput("m_tgi", width = 4), valueBoxOutput("m_dm1", width = 4)),
        fluidRow(
          box(title = tags$div(tags$strong("Tumour volume, 60 virtual mice"),
                               tags$span(textOutput("engine_note", inline = TRUE),
                                         style = "color: #6b7078; font-size: 13px; font-weight: normal;")),
              status = "primary", solidHeader = TRUE, width = 7, plotOutput("m_tv_plot", height = "400px"),
              muted("Median and 50% / 90% ranges from the killing constant's inter-animal variability (Table I); grey: untreated.")),
          box(title = "T-DM1 and DM1 in the tumour", status = "primary", solidHeader = TRUE, width = 5,
              plotOutput("m_tum_plot", height = "400px"),
              muted("Intracellular DM1 catabolites (free + tubulin-bound) drive killing: rate = KKILL x DM1."))
        )
      ),
      tabItem(
        tabName = "human",
        fluidRow(valueBoxOutput("h_cmax", width = 3), valueBoxOutput("h_trough", width = 3),
                 valueBoxOutput("h_dm1", width = 3), valueBoxOutput("h_tum", width = 3)),
        fluidRow(
          box(title = "Plasma: total trastuzumab, T-DM1 and DM1 catabolites", status = "primary", solidHeader = TRUE, width = 6,
              plotOutput("h_pk_plot", height = "400px"),
              muted("Human PK scaled from the monkey fit (exponent 1 for the antibody, 0.75 for DM1 clearances). The paper notes that DM1 is under-predicted in patients.")),
          box(title = "Tumour: intracellular DM1 by regimen", status = "primary", solidHeader = TRUE, width = 6,
              plotOutput("h_reg_plot", height = "400px"),
              muted("All four regimens of the paper's Fig. 6 at the chosen HER2 level. Fractionated dosing keeps tumour DM1 higher between doses."))
        )
      ),
      tabItem(
        tabName = "course",
        fluidRow(
          box(title = "Tumour course in patients - illustrative", status = "warning", solidHeader = TRUE, width = 4,
              tags$p("The paper translates to patients with clinically reported growth rates (Table III), but it
                      does not give the conversion from the linear-phase doubling time to the growth
                      equation's rate, and the published progression-free survival could not be reproduced
                      from the text. This tab therefore uses exponential growth with a doubling time you
                      choose, and the paper's clinical killing constants (HER2 3+: 1.8e-4, 1+: 3.8e-5 per day
                      per nM, 60% variability)."),
              numericInput("dt", "Tumour doubling time (days)", 120, min = 10, max = 1000, step = 10),
              numericInput("len", "Lesion length at start (mm)", 19, min = 5, max = 100, step = 1),
              fluidRow(column(6, numericInput("npat", "Virtual patients", 500, min = 20, max = 5000, step = 100)),
                       column(6, numericInput("months", "Follow (months)", 18, min = 3, max = 36, step = 3))),
              muted("Treatment continues on the chosen regimen for the whole follow-up, as in the trials. For speed, every patient's tumour is driven by the intracellular DM1 profile of a reference tumour (tumour size changes its own exposure only through the small surface-exchange term). Tumour volume = 1/2 x L x B^2 with B = 16/19 L (the median lengths of Table III); lengths vary with the CV of Table III (101%). Progression: diameter 20% above its smallest value.")),
          box(title = "Tumour diameter and progression-free fraction", status = "primary", solidHeader = TRUE, width = 8,
              plotOutput("course_plot", height = "440px"))
        )
      ),
      tabItem(
        tabName = "setup",
        fluidRow(
          box(title = "Parameters", status = "primary", solidHeader = TRUE, width = 7, DTOutput("par_table")),
          box(title = "Mouse tumour models (Table I)", status = "primary", solidHeader = TRUE, width = 5,
              DTOutput("model_table"))
        )
      ),
      about_tab(
        "T-DM1: a multiscale PK-PD model for translating ADC efficacy",
        "The Shah lab's strategy for taking an antibody-drug conjugate from mouse tumour models to
         patients, validated with trastuzumab emtansine. Plasma PK of three analytes is joined to a
         mechanistic tumour model in which the ADC crosses tumour vessels and the tumour surface, binds
         HER2, is internalised and degraded, and releases its DM1 payload inside the cell. The
         intracellular DM1 concentration drives tumour-cell killing; the killing constants come from
         eleven mouse models, and the human PK from allometric scaling of monkey data.",
        list("Plasma: total trastuzumab and conjugated T-DM1 (two compartments each; T-DM1 also loses payload by nonspecific deconjugation), DM1 catabolites (two compartments), average DAR",
             "Tumour: Krogh-cylinder vascular exchange and spherical surface exchange for ADC and DM1; HER2 binding, internalisation, lysosomal degradation releasing DAR x DM1; tubulin binding and diffusion of DM1",
             "Growth and killing: exponential growth in mice; dTV/dt = (growth - KKILL x intracellular DM1) x TV",
             "Mass-balance corrections of printed typos (Eqs. 3, 4, 8, 9, 13, 16) are listed in the model file; DAR is reset to 3.5 with every dose",
             "Parameters: Table I (cellular, tumour, systemic PK in mouse and human, TGI in 11 mouse models), Table III (clinical tumour size)"),
        tags$span("Singh AP, Shah DK. Application of a PK-PD modeling and simulation-based strategy for
                   clinical translation of antibody-drug conjugates: a case study with trastuzumab emtansine
                   (T-DM1). AAPS J 2017;19:1054-1070. Singh AP, Maass KF, Betts AM, et al. Evolution of
                   antibody-drug conjugate tumor disposition model to predict preclinical tumor
                   pharmacokinetics of trastuzumab-emtansine (T-DM1). AAPS J 2016;18:861-875."),
        M$engine$name,
        extra = tagList(
          tags$h4("Further reading"),
          tags$ul(
            tags$li("Shah DK, Haddish-Berhane N, Betts A. Bench to bedside translation of antibody drug conjugates using a multiscale mechanistic PK/PD model: a case study with brentuximab-vedotin. J Pharmacokinet Pharmacodyn 2012;39:643-659."),
            tags$li("Betts AM, Haddish-Berhane N, Tolsma J, et al. Preclinical to clinical translation of antibody-drug conjugates using PK/PD modeling: a retrospective analysis of inotuzumab ozogamicin. AAPS J 2016;18:1101-1116."),
            tags$li("Singh AP, Seigel GM, Guo L, et al. Evolution of the systems PK-PD model for antibody-drug conjugates to characterize tumor heterogeneity and in vivo bystander effect. J Pharmacol Exp Ther 2020;374:184-199."),
            tags$li("Chang HP, Shah DK. A translational physiologically-based pharmacokinetic model for MMAE-based antibody-drug conjugates. J Pharmacokinet Pharmacodyn 2025;52:27.")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {

  output$engine_note <- renderText(sprintf("  (%s)", M$engine$name))
  mm <- reactive(MOUSE_MODELS[MOUSE_MODELS$model == input$mmodel, ])

  output$mmodel_note <- renderUI({
    x <- mm()
    muted(sprintf("%s, HER2 %s. Doubling time %.1f days; killing constant %.2e per day per nM (IIV %.0f%%). Studied: %s.",
                  x$type, x$her2, x$dt, x$kkill, 100 * x$iiv, x$studied))
  })

  # --- mouse -----------------------------------------------------------------------

  mouse_sim <- reactive({
    x <- mm()
    shiny::req(input$mdose >= 0, input$mn >= 1, input$mii >= 1, input$mtv0 > 0, input$mdays > 0)
    set.seed(11)
    n <- 60
    base <- c(PK_SETS$mouse, list(GLIN = 0, DTEXP = x$dt, AG = HER2_AG[[x$her2]], TV0 = input$mtv0 * 1e-6))
    P <- as.data.frame(base)[rep(1, n + 1), ]
    P$KKILL <- c(0, draw_ln(n, x$kkill, x$iiv))
    times <- c(0, seq(0.01, input$mdays, by = 1))       # just after, not at, dose times
    dt <- seq(0, by = input$mii, length.out = input$mn)
    ev <- tdm1_doses(dt, input$mdose, base$BW)
    ev_ctrl <- ev; ev_ctrl$ID <- 1; ev_ctrl$amt[ev_ctrl$cmt != "DAR"] <- 0
    ev_trt <- do.call(rbind, lapply(2:(n + 1), function(i) transform(ev, ID = i)))
    r <- withProgress(message = "Simulating 60 mice", value = 0.4,
                      M$engine$solve(M$model, P, rbind(ev_ctrl, ev_trt), times, rtol = 1e-4, atol = 1e-8, nonneg = TRUE))
    list(r = r, t = times, n = n)
  }) |> debounce(600)

  output$m_tv <- renderValueBox({
    x <- mouse_sim(); k <- length(x$t)
    stat_box(sprintf("%.0f", stats::median(x$r$TV_MM3[k, -1])), " mm3", "Median tumour at the end", sprintf("Untreated: %.0f mm3", x$r$TV_MM3[k, 1]), "circle")
  })
  output$m_tgi <- renderValueBox({
    x <- mouse_sim(); k <- length(x$t)
    v0 <- x$r$TV_MM3[1, 1]
    tgi <- 100 * (1 - (stats::median(x$r$TV_MM3[k, -1]) - v0) / (x$r$TV_MM3[k, 1] - v0))
    stat_box(sprintf("%.0f", tgi), "%", "Tumour growth inhibition", "Median, relative to untreated growth", "percent")
  })
  output$m_dm1 <- renderValueBox({
    x <- mouse_sim()
    stat_box(sprintf("%.0f", max(x$r$DM1_TUMOUR[, 2])), " nM", "Peak intracellular DM1", "Tumour, first animal", "arrow-up")
  })

  output$m_tv_plot <- renderPlot({
    x <- mouse_sim()
    b <- summarise_bands(x$t, x$r$TV_MM3[, -1, drop = FALSE])
    band_plot(b, "Day", "Tumour volume (mm3)") +
      geom_line(data = data.frame(time = x$t, v = x$r$TV_MM3[, 1]), aes(time, v), colour = PAL$ink_3, linewidth = 0.9, linetype = "22")
  })

  output$m_tum_plot <- renderPlot({
    x <- mouse_sim()
    lv <- c("T-DM1 in plasma (nM)", "Intracellular DM1 (nM)", "HER2 occupied (%)")
    d <- rbind(data.frame(time = x$t, value = x$r$ADC_UGML[, 2] * 1e3 / 148.5, series = lv[1]),
               data.frame(time = x$t, value = x$r$DM1_TUMOUR[, 2], series = lv[2]),
               data.frame(time = x$t, value = 100 * x$r$RO_HER2[, 2], series = lv[3]))
    d$value <- ifelse(d$value > 1e-3, d$value, NA)
    series_plot(d, "Day", "Log scale", lv) + scale_y_log10(labels = plain_number)
  })

  # --- patients: exposure -------------------------------------------------------------

  human_P <- reactive({
    P <- as.data.frame(PK_SETS$human); P$BW <- input$bw
    P$AG <- HER2_AG[[input$her2]]; P$KKILL <- 0
    P
  })

  human_sim <- reactive({
    shiny::req(input$cycles >= 1, input$bw > 0)
    reg <- REGIMENS[[input$reg]](input$cycles)
    reg <- reg[reg$time < 21 * input$cycles, ]
    times <- seq(0.01, 21 * input$cycles + 21, by = 0.25)
    r <- withProgress(message = "Simulating patient PK", value = 0.4,
                      M$engine$solve(M$model, human_P(), regimen_events(reg, input$bw), times, rtol = TOL$rtol, atol = TOL$atol, nonneg = TRUE))
    list(r = r, t = times, reg = reg)
  }) |> debounce(500)

  output$h_cmax <- renderValueBox({
    x <- human_sim(); stat_box(sprintf("%.0f", max(x$r$ADC_UGML[x$t <= 21, 1])), " ug/mL", "T-DM1 peak", "First cycle", "arrow-up")
  })
  output$h_trough <- renderValueBox({
    x <- human_sim(); k <- which.min(abs(x$t - 20.75))
    stat_box(sprintf("%.2f", x$r$ADC_UGML[k, 1]), " ug/mL", "T-DM1 on day 21", "Before the next cycle", "arrow-down")
  })
  output$h_dm1 <- renderValueBox({
    x <- human_sim(); stat_box(sprintf("%.1f", max(x$r$DM1_NGML[x$t <= 21, 1])), " ng/mL", "DM1 peak in plasma", "First cycle", "vial")
  })
  output$h_tum <- renderValueBox({
    x <- human_sim(); w <- x$t > 21 & x$t <= 21 * input$cycles
    if (!any(w)) w <- x$t > 0
    stat_box(sprintf("%.0f", mean(x$r$DM1_TUMOUR[w, 1])), " nM", "Mean intracellular DM1", "Tumour, after the first cycle", "bullseye")
  })

  output$h_pk_plot <- renderPlot({
    x <- human_sim()
    lv <- c("Total trastuzumab (ug/mL)", "T-DM1 (ug/mL)", "DM1 catabolites (ng/mL)")
    d <- rbind(data.frame(time = x$t, value = x$r$TT_UGML[, 1], series = lv[1]),
               data.frame(time = x$t, value = x$r$ADC_UGML[, 1], series = lv[2]),
               data.frame(time = x$t, value = x$r$DM1_NGML[, 1], series = lv[3]))
    d$value <- ifelse(d$value > 1e-3, d$value, NA)
    series_plot(d, "Day", "Plasma (log scale)", lv) + scale_y_log10(labels = plain_number)
  })

  output$h_reg_plot <- renderPlot({
    P <- human_P(); cyc <- max(2, min(input$cycles, 6))
    times <- seq(0.01, 21 * cyc, by = 0.5)
    d <- do.call(rbind, lapply(names(REGIMENS), function(nm) {
      reg <- REGIMENS[[nm]](cyc); reg <- reg[reg$time < 21 * cyc, ]
      r <- M$engine$solve(M$model, P, regimen_events(reg, input$bw), times, rtol = TOL$rtol, atol = TOL$atol, nonneg = TRUE)
      data.frame(time = times, value = r$DM1_TUMOUR[, 1], series = nm)
    }))
    series_plot(d, "Day", "Intracellular DM1 in tumour (nM)", names(REGIMENS))
  })

  # --- patients: illustrative tumour course --------------------------------------------

  # One reference tumour gives the intracellular DM1 profile; each patient's
  # tumour then follows dTV/dt = (ln2/DT - KKILL x DM1(t)) x TV exactly. The
  # full model feeds tumour size back into exposure only through the surface
  # exchange term (6D/R^2), which is small beside vascular exchange for
  # tumours of clinical size, so this is fast and close; the mouse tab keeps
  # the full per-animal model.
  course <- reactive({
    shiny::req(input$dt > 0, input$len > 0, input$npat >= 10, input$months > 0)
    set.seed(7)
    n <- input$npat
    k <- KKILL_CLIN[[input$her2]]
    if (is.na(k)) k <- sqrt(KKILL_CLIN[["3+"]] * KKILL_CLIN[["1+"]])
    L <- draw_ln(n, input$len, 1.01); B <- L * 16 / 19
    tv0 <- 0.5 * L * B^2                                  # mm3
    kk <- draw_ln(n, k, 0.6)
    days <- input$months * 30.44
    reg <- REGIMENS[[input$reg]](ceiling(days / 21)); reg <- reg[reg$time < days, ]
    P <- human_P(); P$GLIN <- 0; P$DTEXP <- 1e12
    P$TV0 <- 0.5 * input$len * (input$len * 16 / 19)^2 * 1e-6
    times <- c(0, seq(0.01, days, by = 1))
    r <- withProgress(message = "Simulating tumour exposure", value = 0.4,
                      M$engine$solve(M$model, P, regimen_events(reg, input$bw), times, rtol = 1e-4, atol = 1e-8, nonneg = TRUE))
    cum <- c(0, cumsum(diff(times) * (utils::head(r$DM1_TUMOUR[, 1], -1) + utils::tail(r$DM1_TUMOUR[, 1], -1)) / 2))
    g <- log(2) / input$dt
    tv_ctrl <- outer(exp(g * times), tv0)
    tv_trt <- exp(outer(g * times, rep(1, n)) - outer(cum, kk)) * rep(tv0, each = length(times))
    diam <- function(tv) 2 * (3 * tv / (4 * pi))^(1 / 3)
    prog <- function(D) apply(D, 2, function(d) { i <- which(d >= 1.2 * cummin(d) & seq_along(d) > 1)[1]; if (is.na(i)) Inf else times[i] })
    list(t = times, n = n, pfs_ctrl = prog(diam(tv_ctrl)), pfs_trt = prog(diam(tv_trt)))
  }) |> debounce(800)

  output$course_plot <- renderPlot({
    x <- course()
    surv <- function(e) vapply(x$t, function(tt) mean(e > tt), 0)
    d <- rbind(data.frame(time = x$t / 30.44, value = 100 * surv(x$pfs_ctrl), series = "Untreated"),
               data.frame(time = x$t / 30.44, value = 100 * surv(x$pfs_trt), series = paste0("T-DM1, HER2 ", input$her2)))
    series_plot(d, "Month", "Progression-free (%)", unique(d$series)) + scale_y_continuous(limits = c(0, 100))
  })

  # --- setup -------------------------------------------------------------------------

  output$par_table <- renderDT({
    p <- M$model$param
    small_table(data.frame(Parameter = names(p), Default = signif(unname(p), 4)), page = 60)
  })
  output$model_table <- renderDT({
    x <- MOUSE_MODELS[, c("model", "her2", "dt", "kkill", "iiv")]
    names(x) <- c("Model", "HER2", "Doubling (d)", "KKILL", "IIV")
    small_table(x)
  })
}

shinyApp(ui, server)
