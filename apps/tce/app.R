# ============================================================================
# CD20xCD3 T-cell engager: mosunetuzumab and cytokine release (Hosseini 2020)
#
# A bispecific antibody tethers T cells to B cells. T cells activate, kill
# B cells, proliferate and release cytokines. Because activation needs both
# drug and B cells, the first dose - when B cells are plentiful - gives the
# largest IL-6 surge. Step-up dosing depletes B cells at a low dose first,
# the strategy the model supported for the Phase I trial of mosunetuzumab.
# ============================================================================

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(scales)

shared <- function(f) source(if (file.exists(file.path("R", f))) file.path("R", f) else file.path("..", "..", "shared", f))
for (f in c("theme.R", "ode_engine.R", "app_helpers.R")) shared(f)

M <- load_model("tce_mosunetuzumab_hosseini2020.cpp")

# 'Tumor compartment - DLBCL' variant of the authors' project
TUMOUR <- data.frame(tumor_on = 1, Vtumor = 100, KTrptumor = 100, KBptumor = 1600, Kptumor = 0.05,
                     kBtumorprolif = 0.025, IL6_tiss_contribution = 0.0004)

PRESETS <- list(
  stepup = list(d1 = 1, d8 = 2, d15 = 60, c2 = 60, maint = 30),
  full = list(d1 = 60, d8 = 0, d15 = 0, c2 = 60, maint = 30),
  low = list(d1 = 0.4, d8 = 1, d15 = 2.8, c2 = 2.8, maint = 2.8)
)

#' Doses: C1D1, C1D8, C1D15, C2D1 (day 21), then every 21 days
regimen <- function(d1, d8, d15, c2, maint, cycles, bw, id = NA) {
  days <- c(0, 7, 14, 21, if (cycles > 2) 21 * (2:(cycles - 1)))
  mg <- c(d1, d8, d15, c2, rep(maint, max(0, cycles - 2)))
  keep <- mg > 0
  days <- days[keep]; mg <- mg[keep]
  rbind(data.frame(time = days, cmt = "TDBc_ugperkg", amt = mg * 1000 / bw, ID = id),
        data.frame(time = days, cmt = "injection_effect", amt = 1, ID = id),
        data.frame(time = days, cmt = "drug_effect", amt = 1, ID = id))
}

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = tags$div(tags$span("QSP", style = "font-weight: bold; font-size: 24px;"),
                                   tags$span(" T-cell engager", style = "font-size: 18px;")),
                  titleWidth = 350),
  dashboardSidebar(
    width = 320, custom_css,
    sidebarMenu(
      menuItem("Dosing and cytokines", tabName = "main", icon = icon("syringe")),
      menuItem("Cells by tissue", tabName = "cells", icon = icon("th")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    tags$div(
      style = "padding: 18px; padding-top: 8px;",
      tags$h4("Mosunetuzumab regimen (mg, IV)"),
      selectInput("preset", NULL, width = "100%", selected = "stepup",
                  choices = c("Step-up: 1 / 2 / 60 mg, then 60, then 30 q3w" = "stepup",
                              "No step-up: 60 mg from day 1" = "full",
                              "Low fixed dosing" = "low", "Custom" = "custom")),
      fluidRow(
        column(4, numericInput("d1", "Day 1", 1, min = 0, step = 0.5)),
        column(4, numericInput("d8", "Day 8", 2, min = 0, step = 0.5)),
        column(4, numericInput("d15", "Day 15", 60, min = 0, step = 5))
      ),
      fluidRow(
        column(4, numericInput("c2", "Day 22", 60, min = 0, step = 5)),
        column(4, numericInput("maint", "Then q3w", 30, min = 0, step = 5)),
        column(4, numericInput("cycles", "Cycles", 3, min = 1, max = 8, step = 1))
      ),
      checkboxInput("compare", "Compare with 60 mg from day 1", TRUE),
      tags$h4("Patient"),
      fluidRow(
        column(6, numericInput("bw", "Weight (kg)", 70, min = 30, max = 150, step = 5)),
        column(6, numericInput("bcell", "Blood B cells (/uL)", 500, min = 1, max = 5000, step = 50))
      ),
      numericInput("tcell", "Blood CD8+ T cells (/uL)", 500, min = 50, max = 3000, step = 50, width = "50%"),
      muted("Many lymphoma patients have few blood B cells after rituximab; B cells in spleen, nodes and marrow still start at their normal levels here."),
      checkboxInput("tumour", "Add a lymphoma mass (100 mL, DLBCL setting)", FALSE)
    )
  ),
  dashboardBody(
    custom_css,
    tabItems(
      tabItem(
        tabName = "main",
        fluidRow(valueBoxOutput("il6", width = 3), valueBoxOutput("il6_cmp", width = 3),
                 valueBoxOutput("bdep", width = 3), valueBoxOutput("tnadir", width = 3)),
        fluidRow(
          box(title = tags$div(tags$strong("Peak IL-6 after each dose"),
                               tags$span(textOutput("engine_note", inline = TRUE),
                                         style = "color: #6b7078; font-size: 13px; font-weight: normal;")),
              status = "primary", solidHeader = TRUE, width = 7, plotOutput("il6_plot", height = "360px"),
              muted("Highest blood IL-6 between one dose and the next. IL-6 is released when activated T cells engage B cells in the presence of drug, so it tracks the product of drug exposure and B cells remaining. The IL-6 production rate was calibrated in monkeys; compare regimens and doses rather than read absolute values.")),
          box(title = "Mosunetuzumab in plasma", status = "primary", solidHeader = TRUE, width = 5,
              plotOutput("pk_plot", height = "360px"),
              muted("Two-compartment linear PK (human fit: CL 5.4 mL/kg/day, Vc 36.8 mL/kg)."))
        ),
        fluidRow(
          box(title = "B cells and T cells in blood", status = "primary", solidHeader = TRUE, width = 12,
              plotOutput("blood_plot", height = "300px"),
              muted("B cells are killed within days. CD8+ T cells leave the blood after each injection (margination) and return; activated T cells expand after the first doses."))
        )
      ),
      tabItem(
        tabName = "cells",
        fluidRow(
          box(title = "B cells by compartment, relative to baseline", status = "primary", solidHeader = TRUE,
              width = 12, plotOutput("tissue_plot", height = "420px"),
              muted("Activation in tissues is scaled down (fTact = 0.25) and killing needs more activated T cells per B cell (fKmTB_kill = 10), so B cells in spleen, lymph nodes and bone marrow - and a tumour, if added - are cleared more slowly than in blood."))
        )
      ),
      about_tab(
        "CD20xCD3 T-cell engager: mosunetuzumab",
        "A translational QSP model of a CD20/CD3 T-cell-dependent bispecific antibody, calibrated
         on B-cell depletion, T-cell kinetics and IL-6 in cynomolgus monkeys and translated to
         patients. It was used to choose the step-up (fractionated) dosing of the first cycle in
         the Phase I trial of mosunetuzumab in non-Hodgkin lymphoma, to limit cytokine release
         syndrome while keeping B-cell killing.",
        list("Compartments: blood, spleen, lymph nodes, bone marrow (CD19+CD20+ and CD19+CD20- B cells) and an optional tumour",
             "CD8+ T cells: resting, activated and post-activated; activation depends on drug concentration and on the local B:T ratio (Hill functions)",
             "Activated T cells kill B cells (depends on drug and on the activated-T:B ratio), proliferate, and die by activation-induced cell death",
             "Trafficking between blood and tissues with transient margination after each injection",
             "IL-6 produced by activated T cells engaging B cells, cleared with a 20-minute half-life",
             "Two-compartment linear PK of mosunetuzumab (human)"),
        tags$span("Hosseini I, Gadkar K, Stefanich E, Li CC, Sun LL, Chu YW, Ramanujan S. Mitigating the risk
                   of cytokine release syndrome in a Phase I trial of CD20/CD3 bispecific antibody
                   mosunetuzumab in NHL: impact of translational system modeling. npj Syst Biol Appl
                   2020;6:28. Generated from the authors' SimBiology project (Supplementary Software) by
                   tools/build_tce.py, with the variants their scripts activate and human physiology and
                   PK; parameters checked against Supplementary Table 2."),
        M$engine$name
      )
    )
  )
)

server <- function(input, output, session) {

  observeEvent(input$preset, {
    p <- PRESETS[[input$preset]]
    if (is.null(p)) return()
    for (k in names(p)) updateNumericInput(session, k, value = p[[k]])
  }, ignoreInit = TRUE)

  setup <- reactive({
    shiny::req(input$bw > 0, input$bcell > 0, input$tcell > 0, input$cycles >= 1,
               input$d1 >= 0, input$d8 >= 0, input$d15 >= 0, input$c2 >= 0, input$maint >= 0)
    cycles <- min(round(input$cycles), 8)
    P <- data.frame(Bpbo_perml = input$bcell * 1000, Trpbo_perml = input$tcell * 1000,
                    Trpbref_perml = input$tcell * 1000)
    if (input$tumour) P <- cbind(P, TUMOUR)
    compare <- input$compare && !(input$d1 == 60 && input$d8 == 0 && input$d15 == 0)
    ev <- regimen(input$d1, input$d8, input$d15, input$c2, input$maint, cycles, input$bw, id = 1)
    if (compare) {
      P <- P[c(1, 1), , drop = FALSE]
      ev <- rbind(ev, regimen(60, 0, 0, 60, input$maint, cycles, input$bw, id = 2))
    }
    shiny::validate(shiny::need(any(ev$cmt == "TDBc_ugperkg" & ev$ID == 1), "Enter at least one dose."))
    list(P = P, ev = ev, end = 21 * cycles, compare = compare, bw = input$bw)
  }) |> debounce(600)

  sim <- reactive({
    s <- setup()
    # IL-6 has a 20-minute half-life: sample densely in the day after each dose
    dt <- unique(s$ev$time)
    tt <- sort(unique(round(c(seq(0, s$end, by = 0.1), outer(dt, c(seq(0.002, 0.2, by = 0.004), seq(0.25, 1, by = 0.05)), "+")), 4)))
    tt <- tt[tt <= s$end]
    withProgress(message = "Simulating T cells, B cells and IL-6", value = 0.3, {
      r <- M$engine$solve(M$model, s$P, s$ev, tt, rtol = 1e-5, atol = 1e-3)
    })
    list(s = s, t = tt, r = r)
  })

  labels <- function(x) c("Selected regimen", "60 mg from day 1")[seq_len(if (x$s$compare) 2 else 1)]
  output$engine_note <- renderText(sprintf("  (%s)", M$engine$name))

  output$il6 <- renderValueBox({
    x <- sim(); w <- x$t <= 21
    stat_box(sprintf("%.0f", max(x$r$IL6combo[w, 1])), " pg/mL", "IL-6 peak, cycle 1",
             sprintf("On day %.1f", x$t[w][which.max(x$r$IL6combo[w, 1])]), "fire")
  })
  output$il6_cmp <- renderValueBox({
    x <- sim(); w <- x$t <= 21
    if (!x$s$compare) return(stat_box("-", "", "vs 60 mg from day 1", "Comparison off", "balance-scale"))
    ratio <- max(x$r$IL6combo[w, 1]) / max(x$r$IL6combo[w, 2])
    stat_box(sprintf("%.0f", 100 * ratio), " %", "Of the no-step-up peak", "Cycle 1 IL-6 maximum", "balance-scale")
  })
  output$bdep <- renderValueBox({
    x <- sim(); i <- which.min(abs(x$t - 21))
    stat_box(sprintf("%.1f", 100 * (1 - x$r$Bpb_perml[i, 1] / x$r$Bpb_perml[1, 1])), " %", "Blood B cells depleted",
             "Day 21", "arrow-down")
  })
  output$tnadir <- renderValueBox({
    x <- sim(); w <- x$t <= 21
    stat_box(sprintf("%.0f", min(x$r$totTpb_perml[w, 1]) / 1000), " /uL", "Blood CD8+ T-cell nadir",
             "Margination after dosing, cycle 1", "user-minus")
  })

  output$il6_plot <- renderPlot({
    x <- sim(); lv <- labels(x)
    d <- do.call(rbind, lapply(seq_along(lv), function(j) {
      e <- x$s$ev[x$s$ev$cmt == "TDBc_ugperkg" & x$s$ev$ID == j, ]
      ends <- c(e$time[-1], x$s$end + 1e-9)
      pk <- vapply(seq_len(nrow(e)), function(k) max(x$r$IL6combo[x$t >= e$time[k] & x$t < ends[k], j]), 0)
      data.frame(dose = sprintf("Day %s
%s mg", plain_number(e$time + 1), plain_number(signif(e$amt * x$s$bw / 1000, 3))),
                 day = e$time, value = pk, series = lv[j])
    }))
    d$series <- factor(d$series, levels = lv)
    ggplot(d, aes(day, value, fill = series)) +
      geom_col(position = position_dodge(width = 5.4, preserve = "single"), width = 5) +
      scale_fill_manual(values = stats::setNames(SERIES[seq_along(lv)], lv), name = NULL) +
      scale_x_continuous(breaks = unique(d$day), labels = function(b) paste("Day", b + 1)) +
      labs(x = NULL, y = "Peak IL-6 (pg/mL)") + theme_sim(12)
  })

  output$pk_plot <- renderPlot({
    x <- sim(); lv <- labels(x)
    d <- do.call(rbind, lapply(seq_along(lv), function(j) data.frame(time = x$t, value = pmax(x$r$TDBc_ugperml[, j], 1e-4), series = lv[j])))
    series_plot(d, "Day", "Mosunetuzumab (ug/mL)", lv) +
      scale_y_log10(labels = label_number(drop0trailing = TRUE)) + scale_x_continuous(breaks = seq(0, x$s$end, 7))
  })

  output$blood_plot <- renderPlot({
    x <- sim()
    lv <- c("B cells", "CD8+ T cells, total", "CD8+ T cells, activated")
    d <- rbind(data.frame(time = x$t, value = pmax(x$r$Bpb_perml[, 1] / 1000, 1e-3), series = lv[1]),
               data.frame(time = x$t, value = x$r$totTpb_perml[, 1] / 1000, series = lv[2]),
               data.frame(time = x$t, value = pmax(x$r$Tafraction_pb[, 1] * x$r$totTpb_perml[, 1] / 1000, 1e-3), series = lv[3]))
    series_plot(d, "Day", "Cells per uL (log scale)", lv, subtitle = "Selected regimen") +
      scale_y_log10(labels = label_number(drop0trailing = TRUE)) + scale_x_continuous(breaks = seq(0, x$s$end, 7))
  })

  output$tissue_plot <- renderPlot({
    x <- sim()
    comps <- c("Blood" = "Bpb_perml", "Spleen" = "Btiss_perml", "Lymph nodes" = "Btiss2_perml",
               "Bone marrow (CD20+)" = "B1920tiss3_perml", if (x$s$P$tumor_on[1] %in% 1) c("Tumour" = "Btumor_perml"))
    comps <- comps[comps %in% names(x$r)]
    d <- do.call(rbind, lapply(names(comps), function(nm) {
      v <- x$r[[comps[[nm]]]][, 1]
      data.frame(time = x$t, value = pmax(100 * v / v[1], 1e-3), series = nm)
    }))
    series_plot(d, "Day", "B cells (% of baseline, log scale)", names(comps), subtitle = "Selected regimen") +
      scale_y_log10(labels = label_number(drop0trailing = TRUE)) + scale_x_continuous(breaks = seq(0, x$s$end, 7))
  })
}

shinyApp(ui, server)
