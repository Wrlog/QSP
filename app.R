# ============================================================================
# QSP Simulator: monoclonal antibody with target-mediated drug disposition
#
# Two-compartment antibody PK with IV or subcutaneous dosing, a membrane
# target with its own turnover, quasi-steady-state binding and
# internalisation, receptor occupancy and a biomarker driven by free target.
# Used to ask the usual first-in-human and dose-selection questions: which
# dose and interval keep occupancy above target at trough, and what does
# that do to the biomarker?
#
# The structural model lives in models/tmdd.cpp (mrgsolve).
# ============================================================================

# Two engines, same answers. Run locally with mrgsolve installed, the app
# simulates through models/tmdd.cpp. In the browser (webR) nothing can be
# compiled, so it integrates the same ODEs with the Dormand-Prince solver in
# R/tmdd_engine.R. tests/test_engine_vs_mrgsolve.R holds the two together
# and the deploy fails if they drift apart.
library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(scales)

for (f in c("theme", "tmdd_engine", "population")) {
  source(file.path("R", paste0(f, ".R")))
}

ENGINE <- list(name = "Dormand-Prince (base R)", simulate = tmdd_simulate)
if (file.exists(file.path("reference", "mrgsolve_engine.R"))) {
  source(file.path("reference", "mrgsolve_engine.R"))
  if (mrgsolve_ready()) {
    ok <- tryCatch({ tmdd_mrgsolve_model(); TRUE }, error = function(e) FALSE)
    if (ok) {
      ENGINE <- list(name = paste("mrgsolve", utils::packageVersion("mrgsolve")),
                     simulate = tmdd_simulate_mrgsolve)
    }
  }
}

muted <- function(...) tags$p(..., style = "color: #6b7078; font-size: 12px;")
plain_number <- function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)

DOSE_GRID <- c(0.01, 0.03, 0.1, 0.3, 1, 3, 10)

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = tags$div(
      tags$span("QSP", style = "font-weight: bold; font-size: 24px;"),
      tags$span(" mAb TMDD Simulator", style = "font-size: 18px; opacity: 0.9;")
    ),
    titleWidth = 350
  ),

  dashboardSidebar(
    width = 320,
    custom_css,
    sidebarMenu(
      menuItem("Simulation", tabName = "dashboard", icon = icon("chart-line")),
      menuItem("Dose selection", tabName = "doses", icon = icon("bullseye")),
      menuItem("TMDD mechanism", tabName = "mechanism", icon = icon("diagram-project")),
      menuItem("Model setup", tabName = "setup", icon = icon("sliders")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),

    tags$div(
      class = "sidebar-inputs",
      style = "padding: 18px; padding-top: 8px;",

      tags$h4("Patient population", style = "font-size: 15px;"),
      sliderInput("wt_range", "Body weight (kg)", min = 30, max = 150,
                  value = c(50, 100), step = 1, width = "100%"),

      tags$h4("Dosing regimen", style = "font-size: 15px; margin-top: 10px;"),
      radioButtons("route", NULL, choices = c("Subcutaneous" = "sc", "Intravenous" = "iv"),
                   selected = "sc", inline = TRUE),
      fluidRow(
        column(7, numericInput("dose", "Dose", value = 0.3, min = 0, step = 0.1)),
        column(5, selectInput("dose_basis", "Unit", choices = c("mg/kg" = "mgkg", "mg" = "flat")))
      ),
      numericInput("load", "Loading dose (same unit, 0 = none)", value = 0, min = 0, step = 0.1, width = "100%"),
      fluidRow(
        column(6, numericInput("interval", "Interval (days)", value = 14, min = 1, step = 7)),
        column(6, numericInput("n_doses", "Doses", value = 6, min = 1, max = 52, step = 1))
      ),
      numericInput("duration", "Simulation duration (days)", value = 84, min = 7, step = 7, width = "100%"),

      tags$h4("Target", style = "font-size: 15px; margin-top: 10px;"),
      sliderInput("target_ro", "Receptor occupancy target at trough (%)",
                  min = 50, max = 99, value = 90, step = 1, width = "100%")
    )
  ),

  dashboardBody(
    custom_css,
    tabItems(

      # --- Simulation --------------------------------------------------------
      tabItem(
        tabName = "dashboard",
        fluidRow(
          valueBoxOutput("ro_box", width = 3),
          valueBoxOutput("pta_box", width = 3),
          valueBoxOutput("ctrough_box", width = 3),
          valueBoxOutput("bio_box", width = 3)
        ),
        fluidRow(
          box(
            title = tags$div(
              tags$strong("Free drug concentration"),
              tags$span(textOutput("plot_subtitle", inline = TRUE),
                        style = "color: #6b7078; font-size: 14px; font-weight: normal;")
            ),
            status = "primary", solidHeader = TRUE, width = 12, collapsible = TRUE,
            checkboxInput("log_conc", "Log concentration axis", value = TRUE),
            plotOutput("conc_plot", height = "400px")
          )
        ),
        fluidRow(
          box(
            title = "Receptor occupancy", status = "primary", solidHeader = TRUE,
            width = 6, collapsible = TRUE,
            plotOutput("ro_plot", height = "340px")
          ),
          box(
            title = "Biomarker", status = "primary", solidHeader = TRUE,
            width = 6, collapsible = TRUE,
            plotOutput("bio_plot", height = "340px")
          )
        ),
        fluidRow(
          box(
            title = "Response over the final dosing interval", status = "success",
            solidHeader = TRUE, width = 12, collapsible = TRUE,
            DT::dataTableOutput("summary_table"),
            muted("Occupancy is complex / total target. Trough values are at the end of the final interval; minimum occupancy is the lowest point anywhere in it. Medians with the 5th-95th percentile across virtual subjects.")
          )
        )
      ),

      # --- Dose selection ----------------------------------------------------
      tabItem(
        tabName = "doses",
        fluidRow(
          box(
            title = "Minimum occupancy over the final interval, by dose",
            status = "primary", solidHeader = TRUE, width = 6,
            plotOutput("dr_ro_plot", height = "380px")
          ),
          box(
            title = "Share of subjects above the occupancy target",
            status = "primary", solidHeader = TRUE, width = 6,
            plotOutput("dr_pta_plot", height = "380px")
          )
        ),
        fluidRow(
          box(
            title = "Dose-ranging summary", status = "success", solidHeader = TRUE,
            width = 12,
            DT::dataTableOutput("dr_table"),
            muted(textOutput("dr_note", inline = TRUE))
          )
        )
      ),

      # --- Mechanism -----------------------------------------------------------
      tabItem(
        tabName = "mechanism",
        fluidRow(
          box(
            title = "Dose-normalised free drug after a single dose", status = "primary",
            solidHeader = TRUE, width = 6,
            plotOutput("tmdd_plot", height = "420px"),
            muted("If kinetics were linear these curves would sit on top of each other. At low doses the target binds and internalises a large share of the drug, so concentrations collapse once the drug falls below the target capacity. That bend is the TMDD signature.")
          ),
          box(
            title = "Target under the current regimen, typical subject", status = "primary",
            solidHeader = TRUE, width = 6,
            plotOutput("target_plot", height = "420px"),
            muted("Total target accumulates or falls depending on whether the complex is internalised faster (KINT) or slower than free target is degraded (KDEG). Free target is what drives the biomarker.")
          )
        )
      ),

      # --- Setup ---------------------------------------------------------------
      tabItem(
        tabName = "setup",
        fluidRow(
          box(
            title = "Antibody PK (reference weight)", status = "primary",
            solidHeader = TRUE, width = 4,
            fluidRow(
              column(6, numericInput("CL", "CL (L/day)", value = TYPICAL$CL, min = 0.001, step = 0.02)),
              column(6, numericInput("V1", "V1 (L)", value = TYPICAL$V1, min = 0.1, step = 0.1))
            ),
            fluidRow(
              column(6, numericInput("Q", "Q (L/day)", value = TYPICAL$Q, min = 0, step = 0.05)),
              column(6, numericInput("V2", "V2 (L)", value = TYPICAL$V2, min = 0.1, step = 0.1))
            ),
            fluidRow(
              column(6, numericInput("KA", "ka SC (1/day)", value = TYPICAL$KA, min = 0.01, step = 0.05)),
              column(6, numericInput("F", "F SC", value = TYPICAL$F, min = 0.01, max = 1, step = 0.05))
            ),
            fluidRow(
              column(6, numericInput("WT_EXP_CL", "Weight exp. on CL, Q", value = TYPICAL$WT_EXP_CL, step = 0.05)),
              column(6, numericInput("WT_EXP_V", "Weight exp. on V1, V2", value = TYPICAL$WT_EXP_V, step = 0.05))
            ),
            numericInput("MW", "Molecular weight (kDa)", value = TYPICAL$MW / 1000, min = 1, step = 5)
          ),
          box(
            title = "Target and binding", status = "primary", solidHeader = TRUE, width = 4,
            fluidRow(
              column(6, numericInput("R0", "Baseline target R0 (nM)", value = TYPICAL$R0, min = 0.001, step = 0.5)),
              column(6, numericInput("KSS", "Kss (nM)", value = TYPICAL$KSS, min = 0.001, step = 0.1))
            ),
            fluidRow(
              column(6, numericInput("KDEG", "kdeg, free (1/day)", value = TYPICAL$KDEG, min = 0.001, step = 0.05)),
              column(6, numericInput("KINT", "kint, complex (1/day)", value = TYPICAL$KINT, min = 0, step = 0.1))
            ),
            muted("Kss = (koff + kint) / kon. It is the free drug concentration that gives 50% occupancy at quasi-steady state, so 90% occupancy needs free drug at 9 x Kss."),
            tags$hr(),
            tags$h5(tags$strong("Biomarker")),
            fluidRow(
              column(4, numericInput("KOUT", "kout (1/day)", value = TYPICAL$KOUT, min = 0.001, step = 0.05)),
              column(4, numericInput("IMAX", "Imax", value = TYPICAL$IMAX, min = 0, max = 1, step = 0.05)),
              column(4, numericInput("GAMMA", "Gamma", value = TYPICAL$GAMMA, min = 0.1, max = 5, step = 0.1))
            )
          ),
          box(
            title = "Variability and simulation", status = "info", solidHeader = TRUE, width = 4,
            fluidRow(
              column(6, numericInput("cv_CL", "CV on CL (%)", value = 30, min = 0, max = 150, step = 5)),
              column(6, numericInput("cv_V1", "CV on V1 (%)", value = 20, min = 0, max = 150, step = 5))
            ),
            fluidRow(
              column(6, numericInput("cv_R0", "CV on R0 (%)", value = 30, min = 0, max = 150, step = 5)),
              column(6, numericInput("cv_KA", "CV on ka (%)", value = 20, min = 0, max = 150, step = 5))
            ),
            fluidRow(
              column(6, numericInput("n_subjects", "Virtual subjects", value = 200, min = 1, max = 2000, step = 50)),
              column(6, numericInput("seed", "Random seed", value = 123, min = 1, step = 1))
            ),
            numericInput("delta", "Output step (days)", value = 0.5, min = 0.05, max = 7, step = 0.25),
            muted("The dose-selection tab uses up to 300 subjects per dose to stay responsive in the browser.")
          )
        )
      ),

      # --- About ---------------------------------------------------------------
      tabItem(
        tabName = "about",
        box(
          title = "About this application", status = "primary", solidHeader = TRUE, width = 12,
          tags$div(
            style = "padding: 20px;",
            tags$h3("QSP mAb TMDD Simulator"),
            tags$p("A quantitative systems pharmacology model of a monoclonal
                    antibody binding a membrane target. It links dose to
                    exposure, exposure to target engagement and target
                    engagement to a downstream biomarker, which is the chain of
                    reasoning behind most antibody dose selection."),
            tags$h4("Model structure"),
            tags$ul(
              tags$li("Two-compartment antibody PK with linear clearance; IV bolus or first-order SC absorption with bioavailability"),
              tags$li("Target synthesised at kdeg x R0, degraded at kdeg when free; the drug-target complex is internalised at kint"),
              tags$li("Quasi-steady-state binding (Gibiansky et al. 2008): free drug from total drug and total target"),
              tags$li("Biomarker as an indirect response to free target, with maximum suppression Imax"),
              tags$li("Allometric weight effects on CL, Q, V1 and V2; log-normal variability on CL, V1, R0 and ka")
            ),
            tags$h4("Engine"),
            tags$p("This session is using: ", tags$strong(ENGINE$name), ". ",
                   "The model is defined in models/tmdd.cpp for mrgsolve. The
                    browser build integrates the same equations with an
                    adaptive Dormand-Prince 5(4) solver in base R, vectorised
                    across subjects; the test suite checks the two against each
                    other before every deployment."),
            tags$h4("Units"),
            tags$p("Time in days. Drug and target in nM internally; free drug is shown in ug/mL using the molecular weight."),
            tags$hr(),
            tags$p(
              tags$strong("For research and teaching only. "),
              "The default parameters are round numbers typical of an IgG1,
               not estimates for any real antibody or target. Nothing here is
               validated for clinical use.",
              style = "color: #d1453b;"
            )
          )
        )
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  typ <- reactive({
    ids <- c("CL", "V1", "Q", "V2", "KA", "F", "R0", "KSS", "KDEG", "KINT",
             "KOUT", "IMAX", "GAMMA", "WT_EXP_CL", "WT_EXP_V", "MW")
    vals <- lapply(ids, function(i) input[[i]])
    shiny::req(all(vapply(vals, function(v) !is.null(v) && is.finite(v), TRUE)))
    validate(
      need(input$CL > 0 && input$V1 > 0 && input$V2 > 0, "Clearance and volumes must be positive."),
      need(input$R0 > 0 && input$KSS > 0 && input$KDEG > 0, "R0, Kss and kdeg must be positive."),
      need(input$F > 0 && input$F <= 1, "Bioavailability must be between 0 and 1."),
      need(input$IMAX >= 0 && input$IMAX <= 1, "Imax must be between 0 and 1.")
    )
    t <- TYPICAL
    for (i in ids) t[[i]] <- input[[i]]
    t$MW <- input$MW * 1000
    t
  }) |> debounce(300)

  cv <- reactive(list(CL = input$cv_CL, V1 = input$cv_V1, R0 = input$cv_R0, KA = input$cv_KA))

  population <- reactive({
    n <- input$n_subjects
    shiny::req(n, input$seed)
    validate(need(n >= 1 && n <= 2000, "Use between 1 and 2000 virtual subjects."))
    tmdd_individual(n, input$wt_range, typ(), cv(), input$seed)
  })

  regimen <- reactive({
    shiny::req(input$dose, input$interval, input$n_doses, input$duration)
    validate(
      need(input$dose > 0, "Dose must be positive."),
      need(input$interval > 0, "Dosing interval must be positive."),
      need(input$duration >= input$interval, "Simulate at least one dosing interval.")
    )
    list(route = input$route, dose = input$dose, basis = input$dose_basis,
         load = if (isTRUE(input$load > 0)) input$load else NA,
         interval = input$interval, n_doses = as.integer(input$n_doses))
  }) |> debounce(300)

  # Per-subject amounts in mg.
  with_amounts <- function(reg, P) {
    scale <- if (reg$basis == "mgkg") P$WT else rep(1, nrow(P))
    reg$amt <- reg$dose * scale
    reg$load <- if (is.na(reg$load)) NA else reg$load * scale
    reg
  }

  times <- reactive({
    shiny::req(input$delta, input$duration)
    validate(need(input$delta > 0, "Output step must be positive."))
    seq(0, input$duration, by = min(input$delta, input$duration / 20))
  })

  sim <- reactive({
    P <- population()
    reg <- with_amounts(regimen(), P)
    res <- withProgress(message = "Simulating", value = 0.5,
                        ENGINE$simulate(P, reg, times()))
    list(res = res, fi = final_interval(res, reg), reg = reg)
  })

  target <- reactive(input$target_ro / 100)

  # --- Value boxes ------------------------------------------------------------

  stat_box <- function(value, unit, label, sub, icon_name, color = "blue") {
    valueBox(
      value = tags$div(
        tags$span(value, style = "font-size: 34px; font-weight: bold;"),
        tags$span(unit, style = "font-size: 18px;")
      ),
      subtitle = tags$div(tags$strong(label), tags$br(),
                          tags$span(sub, style = "font-size: 11px;")),
      icon = icon(icon_name), color = color, width = NULL
    )
  }

  output$ro_box <- renderValueBox({
    fi <- sim()$fi
    stat_box(sprintf("%.1f", 100 * stats::median(fi$trough_ro)), "%", "Trough occupancy",
             "Median, end of final interval", "bullseye")
  })

  output$pta_box <- renderValueBox({
    fi <- sim()$fi
    v <- 100 * mean(fi$min_ro >= target())
    stat_box(sprintf("%.0f", v), "%", sprintf("Subjects at >= %d%% occupancy", input$target_ro),
             "Throughout the final interval", "check-circle",
             color = if (v >= 90) "green" else if (v >= 70) "yellow" else "red")
  })

  output$ctrough_box <- renderValueBox({
    fi <- sim()$fi
    stat_box(formatC(stats::median(fi$trough_conc), digits = 3, format = "fg", flag = "#"),
             " ug/mL", "Trough free drug", "Median, end of final interval", "tint")
  })

  output$bio_box <- renderValueBox({
    fi <- sim()$fi
    stat_box(sprintf("%+.0f", 100 * (stats::median(fi$bio_end) - 1)), "%", "Biomarker change",
             "Median vs baseline, end of simulation", "arrow-down")
  })

  output$plot_subtitle <- renderText({
    sprintf(" (n = %s virtual subjects, %s)", format(input$n_subjects, big.mark = ","), ENGINE$name)
  })

  # --- Plots -----------------------------------------------------------------------

  band_plot <- function(d, y_lab, title = NULL, subtitle = NULL) {
    ggplot(d, aes(x = time)) +
      geom_ribbon(aes(ymin = q05, ymax = q95), fill = PAL$blue, alpha = 0.14, na.rm = TRUE) +
      geom_ribbon(aes(ymin = q25, ymax = q75), fill = PAL$blue, alpha = 0.28, na.rm = TRUE) +
      geom_line(aes(y = med), colour = PAL$blue_ink, linewidth = 1.1, na.rm = TRUE) +
      labs(x = "Time (days)", y = y_lab, title = title, subtitle = subtitle) +
      scale_x_continuous(breaks = scales::breaks_width(if (max(d$time) > 120) 28 else 7),
                         expand = expansion(c(0, 0.02))) +
      theme_sim(12)
  }

  output$conc_plot <- renderPlot({
    s <- sim()
    d <- summarise_profiles(s$res$time, s$res$conc)
    log_y <- isTRUE(input$log_conc)
    if (log_y) {
      floor_c <- max(d$q95, na.rm = TRUE) * 1e-5
      d[, -1] <- lapply(d[, -1], function(v) ifelse(v > floor_c, v, NA))
    }
    # Free drug giving the occupancy target at quasi-steady state.
    c_target <- typ()$KSS * target() / (1 - target()) * typ()$MW / 1e6
    p <- band_plot(d, "Free drug (ug/mL)",
                   subtitle = "Median with 50% and 90% prediction intervals") +
      geom_hline(yintercept = c_target, colour = PAL$amber, linetype = "22", linewidth = 0.8) +
      annotate("text", x = max(d$time), y = c_target, vjust = -0.5, hjust = 1, size = 3.6,
               colour = PAL$amber_ink,
               label = sprintf("Free drug for %d%% occupancy  %.3g ug/mL", input$target_ro, c_target))
    if (log_y) p + scale_y_log10(labels = plain_number) else p + scale_y_continuous(limits = c(0, NA))
  })

  output$ro_plot <- renderPlot({
    s <- sim()
    d <- summarise_profiles(s$res$time, 100 * s$res$ro)
    band_plot(d, "Receptor occupancy (%)") +
      geom_hline(yintercept = input$target_ro, colour = PAL$amber, linetype = "22", linewidth = 0.8) +
      annotate("text", x = max(d$time), y = input$target_ro, vjust = 1.6, hjust = 1, size = 3.6,
               colour = PAL$amber_ink, label = sprintf("Target %d%%", input$target_ro)) +
      scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))
  })

  output$bio_plot <- renderPlot({
    s <- sim()
    d <- summarise_profiles(s$res$time, 100 * s$res$bio_rel)
    band_plot(d, "Biomarker (% of baseline)") +
      geom_hline(yintercept = 100, colour = PAL$ink_3, linetype = "22") +
      geom_hline(yintercept = 100 * (1 - typ()$IMAX), colour = PAL$ink_3, linetype = "dotted") +
      annotate("text", x = max(d$time), y = 100 * (1 - typ()$IMAX), vjust = -0.6, hjust = 1,
               size = 3.4, colour = PAL$ink_3, label = "Maximum suppression (Imax)") +
      scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.05)))
  })

  output$summary_table <- DT::renderDataTable({
    fi <- sim()$fi
    row <- function(x, fmt) sprintf(paste0(fmt, " (", fmt, " - ", fmt, ")"), stats::median(x),
                                    stats::quantile(x, 0.05), stats::quantile(x, 0.95))
    df <- data.frame(
      Metric = c("Trough occupancy (%)", "Minimum occupancy in final interval (%)",
                 "Trough free drug (ug/mL)", "Free target at end (% of baseline)",
                 "Biomarker at end (% of baseline)",
                 sprintf("Subjects at >= %d%% occupancy throughout the final interval", input$target_ro)),
      Value = c(row(100 * fi$trough_ro, "%.1f"), row(100 * fi$min_ro, "%.1f"),
                row(fi$trough_conc, "%.3g"), row(100 * fi$rfree_end, "%.1f"),
                row(100 * fi$bio_end, "%.1f"),
                sprintf("%.1f%%", 100 * mean(fi$min_ro >= target())))
    )
    DT::datatable(df, options = list(dom = "t", ordering = FALSE), rownames = FALSE) |>
      DT::formatStyle("Metric", fontWeight = "bold", color = "#16181d")
  })

  # --- Dose selection -----------------------------------------------------------------

  dr <- reactive({
    P <- population()
    if (nrow(P) > 300) P <- P[seq_len(300), ]
    reg <- regimen()
    reg$load <- NA
    withProgress(message = "Running the dose sweep", value = 0.5,
                 dose_ranging(DOSE_GRID, P, reg, times(), target(), ENGINE$simulate))
  })

  output$dr_ro_plot <- renderPlot({
    d <- dr()
    ggplot(d, aes(dose, 100 * ro_med)) +
      geom_hline(yintercept = input$target_ro, colour = PAL$amber, linetype = "22", linewidth = 0.8) +
      geom_linerange(aes(ymin = 100 * ro_lo, ymax = 100 * ro_hi), colour = PAL$blue, linewidth = 2.2, alpha = 0.35) +
      geom_line(colour = PAL$blue_ink, linewidth = 1) +
      geom_point(colour = PAL$blue_ink, size = 3) +
      annotate("text", x = min(d$dose), y = input$target_ro, vjust = -0.6, hjust = 0,
               size = 3.6, colour = PAL$amber_ink, label = sprintf("Target %d%%", input$target_ro)) +
      scale_x_log10(breaks = DOSE_GRID, labels = plain_number) +
      scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
      labs(x = "Dose (mg/kg, log scale)", y = "Minimum occupancy, final interval (%)",
           subtitle = "Median with the 5th-95th percentile") +
      theme_sim(12)
  })

  output$dr_pta_plot <- renderPlot({
    d <- dr()
    ggplot(d, aes(factor(plain_number(dose), levels = plain_number(DOSE_GRID)), 100 * pta)) +
      geom_col(fill = PAL$blue, width = 0.65) +
      geom_text(aes(label = sprintf("%.0f%%", 100 * pta)), vjust = -0.5, colour = PAL$ink_2, size = 3.6) +
      geom_hline(yintercept = 90, colour = PAL$ink_3, linetype = "22") +
      scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 20), expand = expansion(c(0, 0))) +
      labs(x = "Dose (mg/kg)", y = sprintf("Subjects at >= %d%% occupancy (%%)", input$target_ro),
           subtitle = "Dashed line: 90% of subjects") +
      theme_sim(12) + theme(panel.grid.major.x = element_blank())
  })

  output$dr_table <- DT::renderDataTable({
    d <- dr()
    df <- data.frame(
      `Dose (mg/kg)` = plain_number(d$dose),
      `Minimum occupancy, median (5th-95th) %` = sprintf("%.1f (%.1f - %.1f)", 100 * d$ro_med, 100 * d$ro_lo, 100 * d$ro_hi),
      `Subjects at target` = sprintf("%.0f%%", 100 * d$pta),
      `Trough free drug (ug/mL)` = sprintf("%.3g", d$ctrough),
      `Biomarker at end (% of baseline)` = sprintf("%.0f", 100 * d$bio_med),
      check.names = FALSE
    )
    DT::datatable(df, options = list(dom = "t", ordering = FALSE), rownames = FALSE) |>
      DT::formatStyle(1, fontWeight = "bold", color = "#16181d")
  })

  output$dr_note <- renderText({
    reg <- regimen()
    lowest <- dr()$dose[dr()$pta >= 0.9]
    sprintf("%s every %g days, %d doses, no loading dose. %s",
            if (reg$route == "sc") "Subcutaneous" else "Intravenous", reg$interval, reg$n_doses,
            if (length(lowest)) sprintf("Lowest dose on the grid with at least 90%% of subjects at target: %s mg/kg.", plain_number(min(lowest)))
            else "No dose on the grid gets 90% of subjects to target.")
  })

  # --- Mechanism ---------------------------------------------------------------------

  output$tmdd_plot <- renderPlot({
    t1 <- typ()
    doses <- DOSE_GRID[DOSE_GRID >= 0.03]
    P1 <- tmdd_individual(1, c(t1$WT_REF, t1$WT_REF), t1)
    Pn <- P1[rep(1, length(doses)), ]
    Pn$ID <- seq_along(doses)
    tt <- seq(0, 56, by = 0.25)
    reg <- list(route = "iv", amt = doses * t1$WT_REF, load = NA, interval = 1e3, n_doses = 1)
    r <- ENGINE$simulate(Pn, reg, tt)
    d <- do.call(rbind, lapply(seq_along(doses), function(k) {
      data.frame(time = tt, y = r$conc[, k] / doses[k], dose = paste(plain_number(doses[k]), "mg/kg"))
    }))
    lv <- paste(plain_number(doses), "mg/kg")
    d$dose <- factor(d$dose, levels = lv)
    d$y <- ifelse(d$y > max(d$y) * 1e-5, d$y, NA)
    ggplot(d, aes(time, y, colour = dose)) +
      geom_line(linewidth = 1, na.rm = TRUE) +
      scale_colour_manual(values = stats::setNames(SERIES[seq_along(lv)], lv)) +
      scale_y_log10(labels = plain_number) +
      scale_x_continuous(breaks = seq(0, 56, 7)) +
      labs(x = "Time (days)", y = "Free drug per mg/kg (ug/mL per mg/kg, log)",
           subtitle = "Typical 70 kg subject, single IV dose") +
      theme_sim(12)
  })

  output$target_plot <- renderPlot({
    t1 <- typ()
    P1 <- tmdd_individual(1, c(t1$WT_REF, t1$WT_REF), t1)
    reg <- with_amounts(regimen(), P1)
    r <- ENGINE$simulate(P1, reg, times())
    lv <- c("Total target", "Free target")
    d <- rbind(
      data.frame(time = r$time, y = 100 * r$rtot[, 1] / P1$R0, what = lv[1]),
      data.frame(time = r$time, y = 100 * r$rfree_rel[, 1], what = lv[2])
    )
    d$what <- factor(d$what, levels = lv)
    ggplot(d, aes(time, y, colour = what)) +
      geom_hline(yintercept = 100, colour = PAL$ink_3, linetype = "22") +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = stats::setNames(SERIES[1:2], lv)) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.05))) +
      labs(x = "Time (days)", y = "% of baseline target",
           subtitle = "Typical 70 kg subject, current regimen") +
      theme_sim(12)
  })
}

shinyApp(ui, server)
