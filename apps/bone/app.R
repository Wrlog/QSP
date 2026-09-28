# ============================================================================
# Bone and mineral homeostasis: denosumab and teriparatide
# (Peterson & Riggs 2010/2012; OpenBoneMin)
#
# Calcium and phosphate balance, PTH, calcitriol, and bone remodelling by
# osteoblasts and osteoclasts under RANK-RANKL-OPG and TGF-beta control.
# Denosumab binds RANKL and switches off osteoclasts; teriparatide is
# intermittent PTH, which stimulates osteoblasts. The outputs are the
# markers and endpoints used in osteoporosis trials: CTx, bone-specific
# ALP, serum calcium, PTH and lumbar-spine BMD.
#
# This app and models/bone_peterson_riggs.cpp are licensed under GPL-3 (see
# LICENSE in this directory), because the model is OpenBoneMin (GPL-3).
# ============================================================================

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(scales)

shared <- function(f) source(if (file.exists(file.path("R", f))) file.path("R", f) else file.path("..", "..", "shared", f))
for (f in c("theme.R", "ode_engine.R", "app_helpers.R")) shared(f)

M <- load_model("bone_peterson_riggs.cpp")

TERI_PMOL_PER_UG <- 1e6 / 4117.8        # teriparatide: 4117.8 g/mol
TERI_MAX_DAYS <- if (isTRUE(M$engine$mrgsolve)) 730 else 56

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = tags$div(tags$span("QSP", style = "font-weight: bold; font-size: 24px;"),
                                   tags$span(" Bone & mineral", style = "font-size: 18px;")),
                  titleWidth = 350),
  dashboardSidebar(
    width = 320, custom_css,
    sidebarMenu(
      menuItem("Treatment response", tabName = "main", icon = icon("bone")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    tags$div(
      style = "padding: 18px; padding-top: 8px;",
      tags$h4("Treatment"),
      radioButtons("drug", NULL, choices = c("Denosumab (anti-RANKL)" = "den",
                                             "Teriparatide (PTH 1-34)" = "teri",
                                             "Untreated" = "none"), selected = "den"),
      conditionalPanel(
        "input.drug == 'den'",
        fluidRow(
          column(6, numericInput("den_dose", "Dose (mg SC)", value = 60, min = 1, max = 210, step = 10)),
          column(6, numericInput("den_ii", "Every (months)", value = 6, min = 1, max = 12, step = 1))
        ),
        numericInput("den_n", "Number of doses", value = 4, min = 1, max = 10, step = 1, width = "100%"),
        numericInput("den_years", "Follow-up (years)", value = 3, min = 1, max = 6, step = 1, width = "100%"),
        muted("Stopping denosumab before the end of follow-up shows the rebound in bone resorption and loss of BMD seen clinically.")
      ),
      conditionalPanel(
        "input.drug == 'teri'",
        fluidRow(
          column(6, numericInput("teri_dose", "Dose (ug SC)", value = 20, min = 5, max = 80, step = 5)),
          column(6, numericInput("teri_days", "Daily for (days)", value = 28, min = 1, max = TERI_MAX_DAYS, step = 7))
        ),
        muted(if (TERI_MAX_DAYS < 100) "Daily PTH pulses need many solver steps; in the browser the course is limited to 8 weeks. Run the app locally with mrgsolve for courses up to 2 years." else "Courses up to 2 years.")
      ),
      conditionalPanel(
        "input.drug == 'none'",
        numericInput("none_years", "Follow-up (years)", value = 2, min = 1, max = 6, step = 1, width = "100%")
      )
    )
  ),
  dashboardBody(
    custom_css,
    tabItems(
      tabItem(
        tabName = "main",
        fluidRow(valueBoxOutput("b1", width = 3), valueBoxOutput("b2", width = 3),
                 valueBoxOutput("b3", width = 3), valueBoxOutput("b4", width = 3)),
        fluidRow(
          box(title = tags$div(tags$strong("Lumbar-spine bone mineral density"),
                               tags$span(textOutput("engine_note", inline = TRUE),
                                         style = "color: #6b7078; font-size: 13px; font-weight: normal;")),
              status = "primary", solidHeader = TRUE, width = 6, plotOutput("bmd_plot", height = "340px")),
          box(title = "Bone turnover markers", status = "primary", solidHeader = TRUE, width = 6,
              plotOutput("marker_plot", height = "340px"),
              muted("Osteoclasts track serum CTx (resorption); osteoblasts track bone-specific ALP / P1NP (formation)."))
        ),
        fluidRow(
          box(title = "Serum calcium and PTH", status = "primary", solidHeader = TRUE, width = 6,
              plotOutput("ca_plot", height = "320px")),
          box(title = "Drug concentration", status = "primary", solidHeader = TRUE, width = 6,
              plotOutput("pk_plot", height = "320px"))
        )
      ),
      about_tab(
        "Bone and mineral homeostasis",
        "A multiscale systems pharmacology model that joins whole-body calcium and
         phosphate homeostasis to bone remodelling at the cellular level, so that a drug
         acting on one signalling protein (RANKL, the PTH receptor) can be followed
         through bone cells and mineral fluxes to bone mineral density years later.",
        list("Calcium: gut absorption (calcitriol-dependent), plasma, exchangeable and non-exchangeable bone pools, renal excretion (PTH-dependent reabsorption)",
             "Phosphate: gut, extracellular, intracellular and bone exchange; renal excretion",
             "PTH gland pool and capacity, PTH secretion regulated by calcium and calcitriol; 1-alpha-hydroxylase and calcitriol",
             "Bone cells: responding and active osteoblasts (fast and slow), osteoclasts; RANK, RANKL, OPG binding; TGF-beta latent and active pools; RUNX2/CREB/BCL2 signalling",
             "Lumbar-spine and femoral-neck BMD driven by osteoblast and osteoclast activity",
             "Denosumab: two-compartment PK with SC absorption and nonlinear elimination; binds RANKL",
             "Teriparatide: SC absorption into the PTH pool"),
        tags$span("Peterson MC, Riggs MM. A physiologically based mathematical model of integrated
                   calcium homeostasis and bone remodeling. Bone 2010;46:49-63. Peterson MC,
                   Riggs MM. Predicting nonlinear changes in bone mineral density over time using
                   a multiscale systems pharmacology model. CPT Pharmacometrics Syst Pharmacol
                   2012;1:e14. Model code: OpenBoneMin (Metrum Research Group), GPL-3."),
        M$engine$name
      )
    )
  )
)

server <- function(input, output, session) {

  scenario <- reactive({
    d <- input$drug
    if (d == "den") {
      shiny::req(input$den_dose > 0, input$den_ii > 0, input$den_n >= 1, input$den_years >= 1)
      ii <- input$den_ii * 24 * 365.25 / 12
      list(drug = d, ev = data.frame(time = 0, cmt = "DENSC", amt = input$den_dose, ii = ii, addl = input$den_n - 1),
           end = input$den_years * 24 * 365.25, step = 24 * 7, stop_t = (input$den_n - 1) * ii + ii,
           units = "months", tscale = 24 * 365.25 / 12)
    } else if (d == "teri") {
      shiny::req(input$teri_dose > 0, input$teri_days >= 1)
      days <- min(input$teri_days, TERI_MAX_DAYS)
      list(drug = d, ev = data.frame(time = 0, cmt = "TERISC", amt = input$teri_dose * TERI_PMOL_PER_UG, ii = 24, addl = days - 1),
           end = (days + 7) * 24, step = 3, stop_t = days * 24, units = "days", tscale = 24)
    } else {
      shiny::req(input$none_years >= 1)
      list(drug = d, ev = NULL, end = input$none_years * 24 * 365.25, step = 24 * 14, stop_t = NA,
           units = "months", tscale = 24 * 365.25 / 12)
    }
  }) |> debounce(500)

  sim <- reactive({
    s <- scenario()
    times <- seq(0, s$end, by = s$step)
    r <- withProgress(message = "Simulating bone and mineral homeostasis", value = 0.4,
                      M$engine$solve(M$model, data.frame(row.names = 1), s$ev, times, rtol = 1e-4, atol = 1e-8))
    bmd <- if (s$drug == "den") r$BMDlsDENchange[, 1] else 100 * (r$states[, "BMDls", 1] - 1)
    list(s = s, t = times / s$tscale, r = r, bmd = bmd)
  })

  output$engine_note <- renderText(sprintf("  (%s)", M$engine$name))

  at_stop <- function(x, v) {
    if (is.na(x$s$stop_t)) return(utils::tail(v, 1))
    v[which.min(abs(x$t * x$s$tscale - x$s$stop_t))]
  }

  output$b1 <- renderValueBox({
    x <- sim()
    stat_box(sprintf("%+.1f", at_stop(x, x$bmd)), "%", "BMD at end of treatment", "Lumbar spine, vs baseline", "chart-line")
  })
  output$b2 <- renderValueBox({
    x <- sim()
    stat_box(sprintf("%+.1f", utils::tail(x$bmd, 1)), "%", "BMD at end of follow-up", "Lumbar spine, vs baseline", "flag-checkered")
  })
  output$b3 <- renderValueBox({
    x <- sim(); v <- x$r$OCchange[, 1] - 100
    k <- which.max(abs(v))
    stat_box(sprintf("%+.0f", v[k]), "%", "Largest CTx change", "Osteoclasts (resorption)", "arrow-down")
  })
  output$b4 <- renderValueBox({
    x <- sim(); v <- x$r$OBchange[, 1] - 100
    k <- which.max(abs(v))
    stat_box(sprintf("%+.0f", v[k]), "%", "Largest BSAP change", "Osteoblasts (formation)", "arrow-up")
  })

  xlab <- function(x) if (x$s$units == "days") "Day" else "Month"
  stop_line <- function(x) {
    if (is.na(x$s$stop_t)) NULL else geom_vline(xintercept = x$s$stop_t / x$s$tscale, colour = PAL$ink_3, linetype = "22")
  }

  output$bmd_plot <- renderPlot({
    x <- sim()
    ggplot(data.frame(time = x$t, v = x$bmd), aes(time, v)) +
      geom_hline(yintercept = 0, colour = PAL$ink_3) + stop_line(x) +
      geom_line(colour = PAL$blue_ink, linewidth = 1.2) +
      labs(x = xlab(x), y = "% change from baseline",
           subtitle = if (is.na(x$s$stop_t)) NULL else "Dashed line: end of the dosing period") +
      theme_sim(12)
  })

  output$marker_plot <- renderPlot({
    x <- sim()
    d <- rbind(data.frame(time = x$t, value = x$r$OCchange[, 1] - 100, series = "CTx (osteoclasts)"),
               data.frame(time = x$t, value = x$r$OBchange[, 1] - 100, series = "BSAP (osteoblasts)"))
    series_plot(d, xlab(x), "% change from baseline", c("CTx (osteoclasts)", "BSAP (osteoblasts)")) +
      geom_hline(yintercept = 0, colour = PAL$ink_3) + stop_line(x)
  })

  output$ca_plot <- renderPlot({
    x <- sim()
    d <- rbind(data.frame(time = x$t, value = x$r$CAchange[, 1] - 100, series = "Serum calcium"),
               data.frame(time = x$t, value = 100 * x$r$PTHpm[, 1] / x$r$PTHpm[1, 1] - 100, series = "PTH"))
    series_plot(d, xlab(x), "% change from baseline", c("Serum calcium", "PTH")) +
      geom_hline(yintercept = 0, colour = PAL$ink_3)
  })

  output$pk_plot <- renderPlot({
    x <- sim()
    if (x$s$drug == "den") {
      d <- data.frame(time = x$t, v = x$r$DENCP[, 1])
      ylab <- "Denosumab (ug/mL)"
    } else if (x$s$drug == "teri") {
      d <- data.frame(time = x$t, v = x$r$PTHpm[, 1])
      ylab <- "Plasma PTH incl. teriparatide (pM)"
    } else {
      d <- data.frame(time = x$t, v = 0)
      ylab <- "No drug"
    }
    ggplot(d, aes(time, v)) + geom_line(colour = PAL$blue_ink, linewidth = 1) +
      labs(x = xlab(x), y = ylab) + theme_sim(12)
  })
}

shinyApp(ui, server)
