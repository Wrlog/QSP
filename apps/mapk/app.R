# ============================================================================
# MAPK pathway inhibitors in BRAF(V600E) colorectal cancer (Kirouac 2017)
#
# EGFR -> RAS -> BRAF/CRAF -> MEK -> ERK signalling with its feedback loops,
# a PI3K/AKT arm, a proliferation signal and tumour growth, driven by the
# PK of cetuximab, vemurafenib, cobimetinib and the ERK inhibitor
# GDC-0994. Virtual patients from the published population give a
# distribution of tumour responses for any combination - the question the
# model was built to answer.
# ============================================================================

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(scales)

shared <- function(f) source(if (file.exists(file.path("R", f))) file.path("R", f) else file.path("..", "..", "shared", f))
for (f in c("theme.R", "ode_engine.R", "app_helpers.R")) shared(f)
source(file.path("R", "mapk.R"))

M <- load_model("mapk_kirouac2017.cpp")

ARMS <- list(
  "Control" = character(0), "Cetuximab" = "CETUX", "Vemurafenib" = "VEMU",
  "Cobimetinib" = "COBI", "GDC-0994" = "ERKI", "Cetuximab + vemurafenib" = c("CETUX", "VEMU"),
  "Vemurafenib + cobimetinib" = c("VEMU", "COBI"), "Cobimetinib + GDC-0994" = c("COBI", "ERKI"),
  "Cetuximab + vemurafenib + cobimetinib" = c("CETUX", "VEMU", "COBI"),
  "Cetuximab + vemurafenib + GDC-0994" = c("CETUX", "VEMU", "ERKI")
)

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = tags$div(tags$span("QSP", style = "font-weight: bold; font-size: 24px;"),
                                   tags$span(" MAPK in BRAF CRC", style = "font-size: 18px;")),
                  titleWidth = 350),
  dashboardSidebar(
    width = 320, custom_css,
    sidebarMenu(
      menuItem("Virtual trial", tabName = "trial", icon = icon("users")),
      menuItem("Signalling and PK", tabName = "signal", icon = icon("diagram-project")),
      menuItem("Compare regimens", tabName = "compare", icon = icon("scale-balanced")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    tags$div(
      style = "padding: 18px; padding-top: 8px;",
      tags$h4("Regimen"),
      checkboxGroupInput("drugs", NULL, choices = stats::setNames(names(MAPK_DRUGS), vapply(MAPK_DRUGS, `[[`, "", "label")),
                         selected = c("CETUX", "VEMU")),
      fluidRow(
        column(6, numericInput("d_CETUX", "Cetuximab mg IV weekly", 450, min = 0, step = 50)),
        column(6, numericInput("d_VEMU", "Vemurafenib mg PO BID", 960, min = 0, step = 120))
      ),
      fluidRow(
        column(6, numericInput("d_COBI", "Cobimetinib mg PO QD", 60, min = 0, step = 20)),
        column(6, numericInput("d_ERKI", "GDC-0994 mg PO QD", 400, min = 0, step = 50))
      ),
      muted("Cobimetinib and GDC-0994 are given 21 days on, 7 off, as in the trials the model was calibrated to."),
      tags$h4("Virtual patients"),
      fluidRow(
        column(6, numericInput("n", "Patients", 12, min = 4, max = 60, step = 4)),
        column(6, numericInput("seed", "Seed", 1, min = 1, step = 1))
      )
    )
  ),
  dashboardBody(
    custom_css,
    tabItems(
      tabItem(
        tabName = "trial",
        fluidRow(valueBoxOutput("orr", width = 3), valueBoxOutput("sd", width = 3),
                 valueBoxOutput("pd", width = 3), valueBoxOutput("med", width = 3)),
        fluidRow(
          box(title = tags$div(tags$strong("Tumour size change at day 56, by patient"),
                               tags$span(textOutput("engine_note", inline = TRUE),
                                         style = "color: #6b7078; font-size: 13px; font-weight: normal;")),
              status = "primary", solidHeader = TRUE, width = 7,
              plotOutput("waterfall", height = "400px"),
              muted("RECIST-like categories from tumour size relative to baseline: complete response < 10%, partial response < 70%, progressive disease > 120%.")),
          box(title = "Tumour size over time", status = "primary", solidHeader = TRUE, width = 5,
              plotOutput("tumour", height = "400px"))
        )
      ),
      tabItem(
        tabName = "signal",
        fluidRow(
          box(title = "Pathway activity, median patient", status = "primary", solidHeader = TRUE, width = 7,
              plotOutput("pathway", height = "420px"),
              muted("Activities are relative to the untreated steady state. Watch pERK recover under a single agent as the feedback loops (DUSP, Sprouty) relax and RTK signalling rebounds; that adaptive resistance is why combinations work better.")),
          box(title = "Drug concentrations, median patient", status = "primary", solidHeader = TRUE, width = 5,
              plotOutput("pk", height = "420px"))
        )
      ),
      tabItem(
        tabName = "compare",
        fluidRow(
          box(title = "Compare regimens in the same virtual patients", status = "primary", solidHeader = TRUE, width = 12,
              fluidRow(
                column(8, checkboxGroupInput("arms", NULL, choices = names(ARMS), inline = TRUE,
                                             selected = c("Control", "Vemurafenib", "Cetuximab + vemurafenib",
                                                          "Cetuximab + vemurafenib + cobimetinib"))),
                column(4, actionButton("run_cmp", "Run comparison", icon = icon("play")),
                       muted("Each regimen is simulated in the same patients at the standard doses; this takes a while in the browser."))
              ),
              plotOutput("cmp_plot", height = "420px"),
              DT::dataTableOutput("cmp_table"))
        )
      ),
      about_tab(
        "MAPK pathway inhibition in BRAF(V600E) colorectal cancer",
        "BRAF inhibitors that work well in BRAF-mutant melanoma do little in BRAF-mutant
         colorectal cancer, because inhibiting the pathway relieves feedback and EGFR
         re-activates it. This QSP model captures that biology quantitatively and was
         used to predict clinical responses to an ERK inhibitor before the trial read out.",
        list("Signalling (algebraic, quasi-steady state): RTK1 (EGFR), alternative RTKs, RAS, BRAF, CRAF, MEK, ERK, PI3K, AKT, S6",
             "Four negative feedback loops (DUSP-, Sprouty- and MYC-type, and an AKT loop) with their own time constants",
             "A proliferation signal driven by S6 and logistic tumour growth",
             "PK: cetuximab (IV), vemurafenib, cobimetinib (two-compartment), GDC-0994 and an AKT inhibitor (oral)",
             "Virtual population of 1000 parameter sets with prevalence weights; population-PK variability for each drug"),
        tags$span("Kirouac DC, Schaefer G, Chan J, et al. Clinical responses to ERK inhibition in
                   BRAF(V600E)-mutant colorectal cancer predicted using a computational model.
                   npj Syst Biol Appl 2017;3:14 (CC BY 4.0). Converted from the published SimBiology
                   SBML and Supplementary Tables by tools/build_mapk.py."),
        M$engine$name
      )
    )
  )
)

server <- function(input, output, session) {

  patients <- reactive({
    shiny::req(input$n >= 5, input$n <= 60, input$seed)
    mapk_patients(input$n, M$model, input$seed)
  })

  doses <- reactive(list(CETUX = input$d_CETUX, VEMU = input$d_VEMU, COBI = input$d_COBI, ERKI = input$d_ERKI))

  regimen <- reactive(list(drugs = input$drugs, doses = doses())) |> debounce(600)

  simulate <- function(P, drugs, doses, times = seq(0, 56, by = 1)) {
    M$engine$solve(M$model, P, mapk_events(drugs, doses), times, rtol = 1e-4, atol = 1e-7)
  }

  sim <- reactive({
    rg <- regimen()
    withProgress(message = "Simulating virtual patients", value = 0.3, {
      r <- simulate(patients(), rg$drugs, rg$doses)
    })
    cells <- r$states[nrow(r$states), "CELLS", ]
    list(r = r, cells = cells, cat = recist(cells))
  })

  output$engine_note <- renderText(sprintf("  (%s)", M$engine$name))

  share <- function(cat, lv) 100 * mean(cat %in% lv)
  output$orr <- renderValueBox({
    s <- sim(); v <- share(s$cat, c("Complete response", "Partial response"))
    stat_box(sprintf("%.0f", v), "%", "Objective response rate", "CR + PR at day 56", "bullseye")
  })
  output$sd <- renderValueBox({
    s <- sim(); stat_box(sprintf("%.0f", share(s$cat, "Stable disease")), "%", "Stable disease", "Day 56", "equals")
  })
  output$pd <- renderValueBox({
    s <- sim(); stat_box(sprintf("%.0f", share(s$cat, "Progressive disease")), "%", "Progressive disease", "Day 56", "arrow-up")
  })
  output$med <- renderValueBox({
    s <- sim(); stat_box(sprintf("%+.0f", 100 * (stats::median(s$cells) - 1)), "%", "Median tumour change", "Day 56 vs baseline", "chart-column")
  })

  output$waterfall <- renderPlot({
    s <- sim()
    d <- data.frame(change = 100 * (s$cells - 1), cat = s$cat)
    d <- d[order(-d$change), ]
    d$rank <- seq_len(nrow(d))
    ggplot(d, aes(rank, change, fill = cat)) +
      geom_col(width = 0.8) +
      geom_hline(yintercept = c(-30, 20), colour = PAL$ink_3, linetype = "22") +
      scale_fill_manual(values = RECIST_COLS, drop = FALSE) +
      labs(x = "Virtual patient (sorted)", y = "Tumour size change at day 56 (%)", fill = NULL) +
      theme_sim(12) + theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
  })

  output$tumour <- renderPlot({
    s <- sim()
    band_plot(summarise_bands(s$r$time, 100 * s$r$states[, "CELLS", ]), "Day", "Tumour size (% of baseline)",
              "Median with 50% and 90% ranges across patients") +
      geom_hline(yintercept = c(70, 120), colour = PAL$ink_3, linetype = "22")
  })

  median_patient <- reactive({
    s <- sim(); k <- which.min(abs(s$cells - stats::median(s$cells)))
    list(s = s, k = k)
  })

  output$pathway <- renderPlot({
    x <- median_patient(); r <- x$s$r
    nodes <- c("RTK1", "RAS", "BRAF", "MEK", "ERK", "AKT", "S6")
    d <- do.call(rbind, lapply(nodes, function(nd) {
      v <- r[[nd]][, x$k]; data.frame(time = r$time, value = 100 * v / v[1], series = nd)
    }))
    series_plot(d, "Day", "% of untreated activity", nodes)
  })

  output$pk <- renderPlot({
    x <- median_patient(); r <- x$s$r
    lab <- c(RTK1i_C = "Cetuximab", RAFi_C = "Vemurafenib", MEKi_C = "Cobimetinib", ERKi_C = "GDC-0994")
    d <- do.call(rbind, lapply(names(lab), function(v) data.frame(time = r$time, value = r[[v]][, x$k], series = lab[[v]])))
    d <- d[d$series %in% lab[paste0(c(CETUX = "RTK1i_C", VEMU = "RAFi_C", COBI = "MEKi_C", ERKI = "ERKi_C")[input$drugs])], ]
    if (!nrow(d)) d <- data.frame(time = r$time, value = 0, series = "No drug")
    series_plot(d, "Day", "Concentration (mg/L)", unique(d$series)) + scale_y_continuous(trans = "sqrt")
  })

  cmp <- eventReactive(input$run_cmp, {
    arms <- input$arms
    shiny::req(length(arms) > 0)
    P <- patients()
    std <- list(CETUX = 450, VEMU = 960, COBI = 60, ERKI = 400)
    withProgress(message = "Comparing regimens", value = 0, {
      do.call(rbind, lapply(seq_along(arms), function(i) {
        incProgress(1 / length(arms), detail = arms[i])
        r <- simulate(P, ARMS[[arms[i]]], std, times = c(0, 56))
        cat <- recist(r$states[2, "CELLS", ])
        data.frame(arm = arms[i], ORR = share(cat, c("Complete response", "Partial response")),
                   SD = share(cat, "Stable disease"), PD = share(cat, "Progressive disease"))
      }))
    })
  })

  output$cmp_plot <- renderPlot({
    d <- cmp()
    d$arm <- factor(d$arm, levels = d$arm[order(d$ORR)])
    ggplot(d, aes(ORR, arm)) +
      geom_col(fill = PAL$blue, width = 0.65) +
      geom_text(aes(label = sprintf("%.0f%%", ORR)), hjust = -0.2, colour = PAL$ink_2, size = 3.8) +
      scale_x_continuous(limits = c(0, max(100, max(d$ORR) + 10)), expand = expansion(c(0, 0))) +
      labs(x = "Objective response rate at day 56 (%)", y = NULL,
           subtitle = sprintf("%d virtual patients per regimen, same patients in every arm", nrow(patients()))) +
      theme_sim(12) + theme(panel.grid.major.y = element_blank())
  })

  output$cmp_table <- DT::renderDataTable({
    d <- cmp()
    small_table(data.frame(Regimen = d$arm, `ORR %` = sprintf("%.0f", d$ORR), `SD %` = sprintf("%.0f", d$SD),
                           `PD %` = sprintf("%.0f", d$PD), check.names = FALSE))
  })
}

shinyApp(ui, server)
