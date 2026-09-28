# ============================================================================
# PROTAC targeted protein degradation: the kcat model (Pharmaceutics 2023)
#
# A PROTAC brings a target protein and an E3 ligase together; the ternary
# complex gets the target ubiquitinated and destroyed, and the PROTAC is
# released to do it again. Binary complexes win at high concentrations, so
# degradation falls off again: the hook effect. The app shows the hook, how
# binding affinities and cooperativity set Dmax, DC50 and DCmax, and how
# degradation plus occupancy-driven inhibition combine into the downstream
# response.
# ============================================================================

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(scales)

shared <- function(f) source(if (file.exists(file.path("R", f))) file.path("R", f) else file.path("..", "..", "shared", f))
for (f in c("theme.R", "ode_engine.R", "app_helpers.R")) shared(f)

M <- load_model("protac_degrader_kcat2023.cpp")

# Zorba et al. (PNAS 2018) BTK degrader series, Table S1 of the paper
# (surface plasmon resonance, equilibrium analysis)
COMPOUNDS <- data.frame(
  id = LETTERS[1:9],
  KDP = c(1535, 489, 1150, 71, 79, 80, 74, 61, 138),
  KDE = c(15700, 5300, 8800, 2500, 2700, 3600, 3200, 3000, 3100),
  COOP = c(0.89, 0.47, 2.50, 0.86, 0.83, 1.05, 1.21, 1.02, 1.34)
)
# Table S2: E3 ligase (cereblon) and BTK levels, BTK half-life
CELLS <- data.frame(
  name = c("Ramos cells", "THP-1 cells", "Rat splenocytes"),
  E0 = c(203, 120, 120), P0 = c(1231, 1311, 427), THALFP = c(16, 16, 70)
)
KCAT_FIT <- 4.6   # 1/h, fitted to Cpd. A in Ramos cells

#' Steady-state degradation and its hook-curve descriptors (Appendix A)
kcat_summary <- function(KDP, KDE, COOP, E0, THALFP, KCAT) {
  k <- COOP * E0 * KCAT * THALFP / log(2)
  s <- k + COOP * E0 + KDE + KDP + 4 * sqrt(KDE * KDP)
  list(
    dss = function(C) k / (k + COOP * E0 + KDP * KDE / C + KDP + KDE + C),
    dmax = k / (k + COOP * E0 + KDE + KDP + 2 * sqrt(KDE * KDP)),
    dcmax = sqrt(KDE * KDP),
    dc50 = 0.5 * (s - sqrt(s^2 - 4 * KDE * KDP))
  )
}

fmt_nm <- function(x) if (x >= 1000) sprintf("%.2f uM", x / 1000) else if (x >= 10) sprintf("%.0f nM", x) else sprintf("%.2g nM", x)

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = tags$div(tags$span("QSP", style = "font-weight: bold; font-size: 24px;"),
                                   tags$span(" PROTAC degraders", style = "font-size: 18px;")),
                  titleWidth = 350),
  dashboardSidebar(
    width = 320, custom_css,
    sidebarMenu(
      menuItem("In vitro: the hook effect", tabName = "vitro", icon = icon("flask")),
      menuItem("In vivo projection", tabName = "vivo", icon = icon("pills")),
      menuItem("Compound series", tabName = "series", icon = icon("table")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    tags$div(
      style = "padding: 18px; padding-top: 8px;",
      tags$h4("PROTAC"),
      selectInput("cpd", NULL, width = "100%", selected = "D",
                  choices = stats::setNames(COMPOUNDS$id, paste("BTK degrader, Cpd.", COMPOUNDS$id))),
      fluidRow(
        column(4, numericInput("kdp", "KD,P (nM)", 71, min = 0.01, step = 10)),
        column(4, numericInput("kde", "KD,E (nM)", 2500, min = 0.01, step = 100)),
        column(4, numericInput("coop", "Coop. alpha", 0.86, min = 0.01, step = 0.1))
      ),
      muted("Binding to the target (KD,P) and to the E3 ligase (KD,E), and how much the ternary complex is stabilised (alpha > 1) or destabilised (alpha < 1)."),
      tags$h4("Cell system"),
      selectInput("cell", NULL, choices = CELLS$name, width = "100%"),
      fluidRow(
        column(4, numericInput("e0", "E3 ligase (nM)", 203, min = 0.1, step = 10)),
        column(4, numericInput("thalf", "Target t1/2 (h)", 16, min = 0.5, step = 1)),
        column(4, numericInput("kcat", "kcat (1/h)", KCAT_FIT, min = 0.01, step = 0.5))
      ),
      tags$h4("Downstream response"),
      checkboxInput("inhib", "Binding also inhibits the target", TRUE),
      fluidRow(
        column(4, numericInput("pdmin", "PDmin", 0, min = 0, max = 0.99, step = 0.05)),
        column(4, numericInput("p50", "P50", 0.5, min = 0.01, max = 0.99, step = 0.05)),
        column(4, numericInput("npd", "n", 1, min = 0.1, max = 5, step = 0.5))
      ),
      muted("P50: fraction of functional target left when the response is half-way. Needs P50 <= 0.5^(1/n); P50 = 0.5 with n = 1 means the response is proportional to the functional target.")
    )
  ),
  dashboardBody(
    custom_css,
    tabItems(
      tabItem(
        tabName = "vitro",
        fluidRow(valueBoxOutput("dmax", width = 3), valueBoxOutput("dc50", width = 3),
                 valueBoxOutput("dcmax", width = 3), valueBoxOutput("te", width = 3)),
        fluidRow(
          box(title = tags$div(tags$strong("Degradation against concentration"),
                               tags$span(textOutput("engine_note", inline = TRUE),
                                         style = "color: #6b7078; font-size: 13px; font-weight: normal;")),
              status = "primary", solidHeader = TRUE, width = 7, plotOutput("hook_plot", height = "360px"),
              muted("Lines at 6 h and 24 h are simulated from the differential equation; the steady-state curve is the closed form. Degradation peaks at DCmax = sqrt(KD,P x KD,E) and falls at higher concentrations as binary complexes take over.")),
          box(title = "Time course at fixed concentrations", status = "primary", solidHeader = TRUE, width = 5,
              plotOutput("time_plot", height = "360px"),
              muted("Target protein relative to baseline. The approach to the new steady state is faster when degradation is deeper."))
        ),
        fluidRow(
          box(title = "Degradation, inhibition and the downstream response at 24 h", status = "primary",
              solidHeader = TRUE, width = 12, plotOutput("tm_plot", height = "320px"),
              muted("Total target modulation TM = D + I - D x I. When binding also inhibits the target, the inhibition rising at high concentration fills in the loss of degradation, so the hook disappears from the downstream response."))
        )
      ),
      tabItem(
        tabName = "vivo",
        fluidRow(
          box(title = "Dosing and illustrative pharmacokinetics", status = "warning", solidHeader = TRUE, width = 12,
              fluidRow(
                column(2, numericInput("dose", "Oral dose (mg)", 200, min = 1, step = 25)),
                column(2, numericInput("tau", "Every (h)", 24, min = 4, step = 4)),
                column(2, numericInput("ndose", "Number of doses", 14, min = 1, max = 60, step = 1)),
                column(2, numericInput("cl", "CL/F (L/h)", 12, min = 0.1, step = 1)),
                column(2, numericInput("vd", "V/F (L)", 150, min = 1, step = 10)),
                column(1, numericInput("fu", "fu", 0.05, min = 0.001, max = 1, step = 0.01)),
                column(1, numericInput("ka", "ka (1/h)", 0.8, min = 0.05, step = 0.1))
              ),
              muted("The paper's framework is in vitro. The one-compartment PK here is an illustration, not a fitted model: it only turns a dose into an unbound concentration so you can see how turnover and the hook play out over time. Molecular weight 950 g/mol."))
        ),
        fluidRow(valueBoxOutput("cu_range", width = 3), valueBoxOutput("deg_max", width = 3),
                 valueBoxOutput("deg_end", width = 3), valueBoxOutput("pd_end", width = 3)),
        fluidRow(
          box(title = "Unbound PROTAC concentration", status = "primary", solidHeader = TRUE, width = 6,
              plotOutput("pk_plot", height = "340px"),
              muted("Dashed: DC50 and DCmax at steady state. Above DCmax more drug gives less degradation.")),
          box(title = "Target protein and downstream response", status = "primary", solidHeader = TRUE, width = 6,
              plotOutput("pd_plot", height = "340px"),
              muted("Target protein P/P0 lags the concentration by the protein's turnover; the response adds occupancy-driven inhibition when that is switched on."))
        )
      ),
      tabItem(
        tabName = "series",
        box(title = "The BTK degrader series in the selected cell system", status = "primary", solidHeader = TRUE,
            width = 12, DTOutput("series_table"),
            muted("Steady-state Dmax, DC50 and DCmax predicted from the binding constants alone, with kcat and the cell parameters above. With kcat fitted to Cpd. A, the paper predicted the Dmax of the other eight compounds in Ramos cells with Pearson rho = 0.97, and anticipated the values in THP-1 cells, where E3 ligase is lower."))
      ),
      about_tab(
        "The kcat model of PROTAC-mediated degradation",
        "A mechanistic pharmacodynamic framework for proteolysis-targeting chimeras. Ternary
         complex formation follows a rapid-equilibrium binding model with cooperativity; the
         ternary complex drives catalytic degradation on top of the target's own turnover.
         Degradation and occupancy-driven inhibition are combined into total target
         modulation, which drives the downstream response. The same equations give closed
         forms for Dmax, DC50 and DCmax, the parameters of the empirical hook model.",
        list("Target engagement TE = TC / P = alpha E0 / (alpha E0 + KD,P + KD,E + KD,P KD,E / C + C), assuming the ternary complex is small relative to E0",
             "Target turnover: d(P/P0)/dt = kdeg,P (1 - P/P0) - kcat TE P/P0, with kdeg,P = ln 2 / t1/2,P",
             "Degradation D = 1 - P/P0; inhibition I = C / (KD,P + C); total modulation TM = D + I - D I",
             "Downstream response PD = PDmin + (1 - PDmin) (1 - TM)^n (1 - P50^n) / (P50^n + (1 - TM)^n (1 - 2 P50^n))",
             "Defaults: BTK degraders of Zorba et al. (PNAS 2018) with cereblon levels and BTK half-life from the paper's Tables S1-S2, kcat = 4.6 1/h",
             "In vivo tab only: an illustrative one-compartment oral PK model (not from the paper) supplies the unbound concentration"),
        tags$span("A Mechanistic Pharmacodynamic Modeling Framework for the Assessment and Optimization
                   of Proteolysis Targeting Chimeras (PROTACs). Pharmaceutics 2023;15:195. Equations (3),
                   (6)-(8), (19), (20) and Appendix A; parameter values from its Supplementary Tables S1-S2."),
        M$engine$name
      )
    )
  )
)

server <- function(input, output, session) {

  observeEvent(input$cpd, {
    r <- COMPOUNDS[COMPOUNDS$id == input$cpd, ]
    updateNumericInput(session, "kdp", value = r$KDP)
    updateNumericInput(session, "kde", value = r$KDE)
    updateNumericInput(session, "coop", value = r$COOP)
  }, ignoreInit = TRUE)
  observeEvent(input$cell, {
    r <- CELLS[CELLS$name == input$cell, ]
    updateNumericInput(session, "e0", value = r$E0)
    updateNumericInput(session, "thalf", value = r$THALFP)
  }, ignoreInit = TRUE)

  pars <- reactive({
    shiny::req(input$kdp > 0, input$kde > 0, input$coop > 0, input$e0 > 0, input$thalf > 0, input$kcat > 0,
               input$pdmin >= 0, input$pdmin < 1, input$p50 > 0, input$npd > 0)
    shiny::validate(shiny::need(input$p50^input$npd <= 0.5 + 1e-9, "P50 must be at most 0.5^(1/n)."))
    data.frame(KDP = input$kdp, KDE = input$kde, COOP = input$coop, E0 = input$e0, THALFP = input$thalf,
               KCAT = input$kcat, INHIB = as.numeric(input$inhib), PDMIN = input$pdmin, P50 = input$p50,
               NPD = input$npd)
  }) |> debounce(400)

  summ <- reactive({ p <- pars(); kcat_summary(p$KDP, p$KDE, p$COOP, p$E0, p$THALFP, p$KCAT) })

  output$engine_note <- renderText(sprintf("  (%s)", M$engine$name))

  # ---- In vitro -----------------------------------------------------------
  CONC <- 10^seq(-1, 5, length.out = 49)
  vitro <- reactive({
    p <- pars()
    P <- cbind(p[rep(1, length(CONC)), ], CFIX = CONC)
    withProgress(message = "Simulating degradation", value = 0.3, {
      M$engine$solve(M$model, P, times = c(0, 6, 24), rtol = 1e-6, atol = 1e-9)
    })
  })

  output$dmax <- renderValueBox({
    s <- summ(); stat_box(sprintf("%.0f", 100 * s$dmax), " %", "Dmax", "Deepest steady-state degradation", "arrow-down")
  })
  output$dc50 <- renderValueBox({
    s <- summ(); stat_box(fmt_nm(s$dc50), "", "DC50", "Half of Dmax, rising side", "crosshairs")
  })
  output$dcmax <- renderValueBox({
    s <- summ(); stat_box(fmt_nm(s$dcmax), "", "DCmax", "Most degradation; hook above", "mountain")
  })
  output$te <- renderValueBox({
    p <- pars(); s <- summ(); C <- s$dcmax
    te <- p$COOP * p$E0 / (p$COOP * p$E0 + p$KDP + p$KDE + p$KDP * p$KDE / C + C)
    stat_box(sprintf("%.2g", 100 * te), " %", "Target in ternary complex", "At DCmax", "link")
  })

  output$hook_plot <- renderPlot({
    r <- vitro(); s <- summ()
    d <- rbind(data.frame(conc = CONC, value = 100 * r$DEG[2, ], series = "6 h"),
               data.frame(conc = CONC, value = 100 * r$DEG[3, ], series = "24 h"),
               data.frame(conc = CONC, value = 100 * s$dss(CONC), series = "Steady state"))
    lv <- c("6 h", "24 h", "Steady state")
    d$series <- factor(d$series, levels = lv)
    ggplot(d, aes(conc, value, colour = series)) +
      geom_vline(xintercept = c(s$dc50, s$dcmax), colour = PAL$ink_3, linetype = "22") +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = stats::setNames(SERIES[1:3], lv)) +
      scale_x_log10(labels = label_number(big.mark = ",", drop0trailing = TRUE)) +
      coord_cartesian(ylim = c(0, 100)) +
      labs(x = "Unbound PROTAC (nM)", y = "Degradation (% of target lost)",
           subtitle = "Dashed lines: DC50 and DCmax") +
      theme_sim(12)
  })

  output$time_plot <- renderPlot({
    p <- pars(); s <- summ()
    cs <- signif(c(s$dc50, s$dcmax, 30 * s$dcmax), 2)
    P <- cbind(p[rep(1, 3), ], CFIX = cs)
    tt <- seq(0, max(72, 4 * p$THALFP), length.out = 121)
    r <- M$engine$solve(M$model, P, times = tt, rtol = 1e-6, atol = 1e-9)
    lv <- sprintf("%s (%s)", c("DC50", "DCmax", "30 x DCmax"), vapply(cs, fmt_nm, ""))
    d <- do.call(rbind, lapply(1:3, function(j) data.frame(time = tt, value = 100 * (1 - r$DEG[, j]), series = lv[j])))
    series_plot(d, "Hours", "Target protein (% of baseline)", lv) + coord_cartesian(ylim = c(0, 100))
  })

  output$tm_plot <- renderPlot({
    r <- vitro()
    lv <- c("Degradation D", "Inhibition I", "Total modulation TM", "Response lost, 1 - PD")
    d <- rbind(data.frame(conc = CONC, value = 100 * r$DEG[3, ], series = lv[1]),
               data.frame(conc = CONC, value = 100 * r$INH[3, ], series = lv[2]),
               data.frame(conc = CONC, value = 100 * r$TM[3, ], series = lv[3]),
               data.frame(conc = CONC, value = 100 * (1 - r$PDR[3, ]), series = lv[4]))
    d$series <- factor(d$series, levels = lv)
    ggplot(d, aes(conc, value, colour = series)) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = stats::setNames(SERIES[1:4], lv)) +
      scale_x_log10(labels = label_number(big.mark = ",", drop0trailing = TRUE)) +
      coord_cartesian(ylim = c(0, 100)) +
      labs(x = "Unbound PROTAC (nM)", y = "% at 24 h") + theme_sim(12)
  })

  # ---- In vivo -----------------------------------------------------------
  vivo <- reactive({
    p <- pars()
    shiny::req(input$dose > 0, input$tau > 0, input$ndose >= 1, input$cl > 0, input$vd > 0,
               input$fu > 0, input$fu <= 1, input$ka > 0)
    P <- cbind(p, CFIX = -1, CL = input$cl, VD = input$vd, FU = input$fu, KA = input$ka)
    ev <- data.frame(time = 0, cmt = "GUT", amt = input$dose, ii = input$tau, addl = input$ndose - 1)
    end <- input$tau * input$ndose + max(72, 3 * p$THALFP)
    tt <- seq(0, end, length.out = 601)
    withProgress(message = "Simulating dosing", value = 0.3, {
      r <- M$engine$solve(M$model, P, ev, tt, rtol = 1e-6, atol = 1e-9)
    })
    list(t = tt, r = r, last = input$tau * (input$ndose - 1), tau = input$tau)
  }) |> debounce(400)

  output$cu_range <- renderValueBox({
    v <- vivo(); w <- v$t >= v$last & v$t <= v$last + v$tau
    rng <- range(v$r$CUNB[w, 1])
    stat_box(sprintf("%.3g-%.3g", rng[1], rng[2]), " nM", "Unbound, last interval", "Trough to peak", "vial")
  })
  output$deg_max <- renderValueBox({
    v <- vivo()
    stat_box(sprintf("%.0f", 100 * max(v$r$DEG)), " %", "Deepest degradation", sprintf("At %.0f h", v$t[which.max(v$r$DEG)]), "arrow-down")
  })
  output$deg_end <- renderValueBox({
    v <- vivo(); i <- which.min(abs(v$t - (v$last + v$tau)))
    stat_box(sprintf("%.0f", 100 * v$r$DEG[i, 1]), " %", "Degradation at trough", "End of the last interval", "hourglass-end")
  })
  output$pd_end <- renderValueBox({
    v <- vivo(); i <- which.min(abs(v$t - (v$last + v$tau)))
    stat_box(sprintf("%.0f", 100 * (1 - v$r$PDR[i, 1])), " %", "Response suppressed", "End of the last interval", "wave-square")
  })

  output$pk_plot <- renderPlot({
    v <- vivo(); s <- summ()
    d <- data.frame(time = v$t, value = pmax(v$r$CUNB[, 1], 1e-3))
    ggplot(d, aes(time, value)) +
      geom_hline(yintercept = c(s$dc50, s$dcmax), colour = PAL$ink_3, linetype = "22") +
      geom_line(colour = PAL$blue_ink, linewidth = 1) +
      scale_y_log10(labels = label_number(big.mark = ",", drop0trailing = TRUE)) +
      labs(x = "Hours", y = "Unbound PROTAC (nM)") + theme_sim(12)
  })

  output$pd_plot <- renderPlot({
    v <- vivo()
    lv <- c("Target protein P/P0", "Downstream response PD")
    d <- rbind(data.frame(time = v$t, value = 100 * (1 - v$r$DEG[, 1]), series = lv[1]),
               data.frame(time = v$t, value = 100 * v$r$PDR[, 1], series = lv[2]))
    series_plot(d, "Hours", "% of baseline", lv) + coord_cartesian(ylim = c(0, 105))
  })

  # ---- Compound series ----------------------------------------------------
  output$series_table <- renderDT({
    p <- pars()
    rows <- lapply(seq_len(nrow(COMPOUNDS)), function(i) {
      c0 <- COMPOUNDS[i, ]
      s <- kcat_summary(c0$KDP, c0$KDE, c0$COOP, p$E0, p$THALFP, p$KCAT)
      data.frame(Compound = paste("Cpd.", c0$id), `KD,P (nM)` = c0$KDP, `KD,E (nM)` = c0$KDE,
                 alpha = c0$COOP, `Dmax (%)` = round(100 * s$dmax), `DC50 (nM)` = signif(s$dc50, 3),
                 `DCmax (nM)` = signif(s$dcmax, 3), check.names = FALSE)
    })
    small_table(do.call(rbind, rows))
  })
}

shinyApp(ui, server)
