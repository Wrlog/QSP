# ============================================================================
# mRNA vaccine: from lipid nanoparticle uptake to antibodies (Dasti 2025)
#
# A multiscale QSP model of mRNA vaccines: LNP uptake and antigen expression
# by innate cells at the injection site, dendritic-cell maturation and
# migration to the lymph node, helper T-cell activation, B-cell activation in
# 17 affinity classes, germinal centres, memory and plasma cells, and the
# antibody response in blood.
#
# Derivative of the COSBI "Multiscale QSP model for mRNA vaccines" (Fondazione
# The Microsoft Research - University of Trento Centre for Computational and
# Systems Biology). COSBI-SSLA licence, non-commercial use only; see LICENSE
# in this folder. Modified 2026-09-28 (translated from MATLAB to R/mrgsolve).
# ============================================================================

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(scales)

shared <- function(f) source(if (file.exists(file.path("R", f))) file.path("R", f) else file.path("..", "..", "shared", f))
for (f in c("theme.R", "ode_engine.R", "app_helpers.R")) shared(f)
source(file.path("R", "mrna_rhs.R"))
source(file.path("R", "vaccines.R"))

M <- load_model("mrna_vaccine_dasti2025.cpp")
RHS <- mrna_rhs(M$model)
PRESETS <- readRDS(file.path("data", "presets.rds"))
SAHIN <- utils::read.csv(file.path("data", "sahin2020_bnt162b2.csv"))

num <- label_number(scale_cut = cut_short_scale(), accuracy = 0.1, drop0trailing = TRUE)

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = tags$div(tags$span("QSP", style = "font-weight: bold; font-size: 24px;"),
                                   tags$span(" mRNA vaccine", style = "font-size: 18px;")),
                  titleWidth = 350),
  dashboardSidebar(
    width = 320, custom_css,
    sidebarMenu(
      menuItem("Antibody response", tabName = "main", icon = icon("shield-alt")),
      menuItem("Cells and affinity", tabName = "cells", icon = icon("project-diagram")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    tags$div(
      style = "padding: 18px; padding-top: 8px;",
      tags$h4("Vaccine"),
      selectInput("vaccine", NULL, width = "100%",
                  choices = stats::setNames(names(VACCINES), vapply(VACCINES, `[[`, "", "label"))),
      fluidRow(
        column(6, numericInput("dose", "Dose (ug mRNA)", 30, min = 1, max = 250, step = 5)),
        column(6, numericInput("interval", "Second dose on day", 21, min = 7, max = 120, step = 7))
      ),
      checkboxInput("boost", "Add a booster", FALSE),
      conditionalPanel("input.boost", numericInput("boost_day", "Booster on day", 180, min = 60, max = 360, step = 30)),
      muted("The standard schedules open instantly (computed with mrgsolve). Other settings are simulated in your browser: 220 equations, which takes up to a minute."),
      actionButton("go", "Simulate", icon = icon("play"), width = "100%")
    )
  ),
  dashboardBody(
    custom_css,
    tabItems(
      tabItem(
        tabName = "main",
        fluidRow(valueBoxOutput("peak", width = 3), valueBoxOutput("prime", width = 3),
                 valueBoxOutput("m6", width = 3), valueBoxOutput("aff", width = 3)),
        fluidRow(
          box(title = tags$div(tags$strong("Antibodies against the encoded antigen (IgG in blood)"),
                               tags$span(textOutput("engine_note", inline = TRUE),
                                         style = "color: #6b7078; font-size: 13px; font-weight: normal;")),
              status = "primary", solidHeader = TRUE, width = 12, plotOutput("igg_plot", height = "400px"),
              uiOutput("igg_note"))
        ),
        fluidRow(
          box(title = "Dose-response: peak and 6-month antibody, BNT162b2 fit", status = "primary", solidHeader = TRUE,
              width = 12, plotOutput("dose_plot", height = "260px"),
              muted("From the standard two-dose schedule (day 1 and day 22) at 1, 10, 20 and 30 ug. The first dose primes; the second, given while germinal centres are active, drives most of the antibody."))
        )
      ),
      tabItem(
        tabName = "cells",
        fluidRow(
          box(title = "The cellular cascade", status = "primary", solidHeader = TRUE, width = 7,
              plotOutput("cells_plot", height = "420px"),
              muted("Antigen-presenting cells reach the lymph node within days of each dose; helper T cells and germinal-centre B cells follow; short-lived plasma cells give the early antibody peak and long-lived plasma cells the slow tail.")),
          box(title = "Affinity maturation", status = "primary", solidHeader = TRUE, width = 5,
              plotOutput("aff_plot", height = "420px"),
              muted("Share of antibody in each of the 17 affinity classes (association constant doubling per class). As free antigen falls, only high-affinity B cells keep their receptors occupied, so the antibody pool shifts towards high affinity, most after the second dose."))
        )
      ),
      about_tab(
        "Multiscale QSP model of mRNA vaccines",
        "A mechanistic model of the immune response to lipid-nanoparticle mRNA vaccines, from
         the injection site to the lymph node and blood. Innate-cell parameters were fitted to
         the non-human primate data of Liang et al.; the adaptive parameters to antibody data for
         BNT162b2 (general population and adults over 60) and mRNA-1273. It was used to explore
         dose, dosing interval and population differences.",
        list("Injection site: LNP-mRNA taken up by neutrophils, monocytes, myeloid and plasmacytoid dendritic cells, which are recruited from blood, express antigen and migrate",
             "Dendritic cells mature reversibly through low, medium and high antigen presentation, in tissue and in the lymph node",
             "Helper T cells: naive, activated, memory and functional; activation depends on presented MHC-II-antigen complexes",
             "B cells in 17 affinity classes: receptor occupancy by free antigen (a binding equilibrium solved at every step) and T-cell help drive activation, germinal centres, memory, short- and long-lived plasma cells",
             "Plasma cells migrate to blood and secrete antibody (IgG, 150 kDa)",
             "220 differential equations"),
        tags$span("Dasti A, et al. A multiscale quantitative systems pharmacology model for mRNA vaccines.
                   CPT Pharmacometrics Syst Pharmacol 2025. Code and fitted parameters:
                   github.com/cosbi-research/QSPmRNAVaccines (tissue layer), translated to an mrgsolve model
                   by tools/build_mrna_vaccine.py. Antibody data: Sahin U, et al. Nature 2020;586:594-599,
                   as distributed with the COSBI code."),
        M$engine$name,
        extra = tags$div(
          tags$h4("Licence"),
          tags$p("This app and its model are a derivative of the ", tags$em("Multiscale QSP model for mRNA vaccines"),
                 ", COSBI, Fondazione The Microsoft Research - University of Trento Centre for Computational and Systems
                 Biology, and are distributed under the COSBI Licence Terms for non-commercial use only (see LICENSE in the
                 app folder). Changed on 2026-09-28: translated from MATLAB to an mrgsolve model file and R; free lymph-node
                 antigen by Newton iteration instead of fzero. The fitted parameters are unchanged. The software comes as is,
                 with no warranties.")
        )
      )
    )
  )
)

server <- function(input, output, session) {

  observeEvent(input$vaccine, {
    v <- VACCINES[[input$vaccine]]
    updateNumericInput(session, "dose", value = v$dose)
    updateNumericInput(session, "interval", value = v$interval)
  }, ignoreInit = TRUE)

  settings <- reactive({
    shiny::req(input$dose > 0, input$interval >= 7)
    list(key = input$vaccine, dose = input$dose, interval = input$interval,
         booster = if (isTRUE(input$boost) && isTruthy(input$boost_day) && input$boost_day > input$interval) input$boost_day else NA)
  })

  # Standard schedules come from the precomputed file; anything else is
  # simulated when the button is pressed.
  result <- reactiveVal(NULL)
  observe({
    s <- settings()
    id <- preset_id(s$key, s$dose, s$interval, s$booster)
    if (!is.null(PRESETS[[id]])) result(list(s = s, r = PRESETS[[id]], source = "precomputed with mrgsolve"))
  })
  observeEvent(input$go, {
    s <- settings()
    su <- vaccine_setup(s$key, s$dose, s$interval, s$booster)
    withProgress(message = "Simulating 220 equations", detail = "up to a minute in the browser", value = 0.2, {
      r <- M$engine$solve(M$model, su$P, su$ev, OUT_TIMES, rtol = 1e-4, atol = 1e-4, nonneg = TRUE, rhs = RHS)
    })
    result(list(s = s, r = vaccine_outputs(r), source = M$engine$name))
  })
  stale <- reactive({
    x <- result(); s <- settings()
    is.null(x) || !identical(preset_id(s$key, s$dose, s$interval, s$booster),
                             preset_id(x$s$key, x$s$dose, x$s$interval, x$s$booster))
  })
  current <- reactive({
    x <- result()
    shiny::validate(shiny::need(!is.null(x), "Press Simulate to run this schedule."))
    x
  })

  output$engine_note <- renderText({
    x <- current()
    sprintf("  (%s%s)", x$source, if (stale()) "; settings changed - press Simulate" else "")
  })

  output$peak <- renderValueBox({
    x <- current(); r <- x$r
    stat_box(num(max(r$IGG) / 1000), " ug/mL", "Peak IgG", sprintf("Day %.0f after the first dose", r$time[which.max(r$IGG)]), "arrow-up")
  })
  output$prime <- renderValueBox({
    x <- current(); r <- x$r
    pre <- r$IGG[which.min(abs(r$time - x$s$interval))]
    stat_box(sprintf("%.0f", max(r$IGG) / max(pre, 1e-9)), " x", "Boost from the second dose",
             sprintf("Peak vs level on day %g", x$s$interval), "layer-group")
  })
  output$m6 <- renderValueBox({
    x <- current(); r <- x$r; i <- which.min(abs(r$time - 180))
    stat_box(num(r$IGG[i] / 1000), " ug/mL", "IgG at 6 months", sprintf("%.0f%% of the peak", 100 * r$IGG[i] / max(r$IGG)), "hourglass-half")
  })
  output$aff <- renderValueBox({
    x <- current(); r <- x$r
    i1 <- which(r$IGG > 0.01 * max(r$IGG))[1]
    gain <- 2^(r$AFFINITY[length(r$time)] - r$AFFINITY[i1])
    stat_box(sprintf("%.1f", gain), " x", "Affinity maturation", "Mean antibody affinity, 1 year vs first response", "bullseye")
  })

  output$igg_plot <- renderPlot({
    x <- current(); r <- x$r
    d <- data.frame(time = r$time, value = pmax(r$IGG, 1))
    days <- c(0, x$s$interval, if (!is.na(x$s$booster)) x$s$booster)
    g <- ggplot(d, aes(time, value)) +
      geom_vline(xintercept = days, colour = PAL$ink_3, linetype = "22") +
      geom_line(colour = PAL$blue_ink, linewidth = 1.1) +
      scale_y_log10(labels = label_number(scale_cut = cut_short_scale(), drop0trailing = TRUE)) +
      scale_x_continuous(breaks = seq(0, 360, 30)) +
      labs(x = "Days after the first dose", y = "IgG (ng/mL)", subtitle = "Dashed lines: doses") + theme_sim(12)
    obs <- SAHIN[SAHIN$dose_ug == x$s$dose, ]
    if (x$s$key == "bnt" && x$s$interval == 21 && nrow(obs)) {
      g <- g + geom_errorbar(data = obs, aes(x = day, ymin = lower95_ng_per_mL, ymax = upper95_ng_per_mL),
                             inherit.aes = FALSE, width = 4, colour = SERIES[2]) +
        geom_point(data = obs, aes(x = day, y = geomean_ng_per_mL), inherit.aes = FALSE, colour = SERIES[2], size = 2.6)
    }
    g
  })
  output$igg_note <- renderUI({
    x <- current()
    if (x$s$key == "bnt" && x$s$interval == 21 && any(SAHIN$dose_ug == x$s$dose)) {
      muted("Orange: geometric mean and 95% CI of RBD-binding IgG in the BNT162b2 phase 1/2 trial (Sahin et al. 2020), converted to ng/mL with the factor used in the source code. The model was fitted to the 30 ug arm; it reproduces the pre-boost and later levels, but its rise after the second dose peaks about a week later than the data.")
    } else {
      muted("Trial data are overlaid for BNT162b2 at 1, 10, 20 or 30 ug with the second dose on day 22.")
    }
  })

  output$dose_plot <- renderPlot({
    ids <- c(1, 10, 20, 30)
    d <- do.call(rbind, lapply(ids, function(dz) {
      r <- PRESETS[[preset_id("bnt", dz, 21)]]
      data.frame(dose = dz, what = c("Peak", "6 months"), value = c(max(r$IGG), r$IGG[which.min(abs(r$time - 180))]))
    }))
    d$what <- factor(d$what, levels = c("Peak", "6 months"))
    ggplot(d, aes(dose, value / 1000, colour = what)) +
      geom_line(linewidth = 1) + geom_point(size = 2.6) +
      scale_colour_manual(values = stats::setNames(SERIES[1:2], levels(d$what)), name = NULL) +
      scale_x_continuous(breaks = ids) +
      labs(x = "Dose (ug mRNA per injection)", y = "IgG (ug/mL)") + theme_sim(12)
  })

  output$cells_plot <- renderPlot({
    x <- current(); r <- x$r
    lv <- c("Antigen-presenting cells, lymph node", "Functional helper T cells", "Germinal-centre B cells",
            "Short-lived plasma cells", "Long-lived plasma cells", "Memory B cells")
    vv <- list(r$APC_LN, r$THELP, r$GCBTOT, r$SPTOT, r$LPTOT, r$MBTOT)
    d <- do.call(rbind, lapply(seq_along(lv), function(j) data.frame(time = r$time, value = pmax(vv[[j]], 1e-2), series = lv[j])))
    series_plot(d, "Days after the first dose", "Cells (log scale)", lv) +
      scale_y_log10(labels = label_number(scale_cut = cut_short_scale(), drop0trailing = TRUE)) +
      coord_cartesian(ylim = c(1, NA)) + guides(colour = guide_legend(ncol = 2))
  })

  output$aff_plot <- renderPlot({
    x <- current(); r <- x$r
    when <- c(x$s$interval, min(x$s$interval + 21, 364), 364)
    lv <- sprintf("Day %g", when)
    d <- do.call(rbind, lapply(seq_along(when), function(k) {
      a <- r$AB[which.min(abs(r$time - when[k])), ]
      data.frame(cls = 1:17, value = 100 * a / sum(a), series = lv[k])
    }))
    d <- d[d$cls >= 10, ]   # classes 1-9 hold well under 1% of antibody
    d$series <- factor(d$series, levels = lv)
    ggplot(d, aes(cls, value, fill = series)) +
      geom_col(position = position_dodge(width = 0.85), width = 0.8) +
      scale_fill_manual(values = stats::setNames(SERIES[seq_along(lv)], lv), name = NULL) +
      scale_x_continuous(breaks = 10:17, labels = c("10", 11:16, "17\nhighest")) +
      labs(x = "Affinity class", y = "% of antibody") + theme_sim(12)
  })
}

shinyApp(ui, server)
