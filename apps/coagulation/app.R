# ============================================================================
# Warfarin and the coagulation network (Wajima et al. 2009)
#
# A 54-species model of the humoral coagulation cascade with warfarin PK,
# the vitamin K cycle, and the in-vitro prothrombin time test. Daily
# warfarin inhibits vitamin K epoxide reduction, the vitamin K-dependent
# factors (II, VII, IX, X, protein C) fall at their own turnover rates, and
# the clotting time of a plasma sample - reported as INR - rises.
# ============================================================================

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(scales)

shared <- function(f) source(if (file.exists(file.path("R", f))) file.path("R", f) else file.path("..", "..", "shared", f))
for (f in c("theme.R", "ode_engine.R", "app_helpers.R")) shared(f)
source(file.path("R", "coagulation.R"))

M <- load_model("coagulation_wajima2009.cpp")

FACTORS <- c("Factor II (prothrombin)" = "FII", "Factor VII" = "VII", "Factor IX" = "IX",
             "Factor X" = "X", "Protein C" = "PC")

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = tags$div(tags$span("QSP", style = "font-weight: bold; font-size: 24px;"),
                                   tags$span(" Coagulation", style = "font-size: 18px;")),
                  titleWidth = 350),
  dashboardSidebar(
    width = 320, custom_css,
    sidebarMenu(
      menuItem("Warfarin and INR", tabName = "main", icon = icon("droplet")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    tags$div(
      style = "padding: 18px; padding-top: 8px;",
      tags$h4("Warfarin"),
      sliderInput("warf", "Daily oral dose (mg)", min = 0, max = 15, value = 5, step = 0.5, width = "100%"),
      sliderInput("days", "Days of treatment", min = 1, max = 28, value = 14, step = 1, width = "100%"),
      tags$h4("Patient"),
      sliderInput("cl_mult", "Warfarin clearance (x typical)", min = 0.3, max = 2, value = 1, step = 0.1, width = "100%"),
      sliderInput("ic50_mult", "Warfarin IC50 (x typical)", min = 0.3, max = 2, value = 1, step = 0.1, width = "100%"),
      muted("Illustrative levers for slow metabolisers (lower clearance) and sensitive patients (lower IC50); the published model has one typical patient."),
      tags$h4("Vitamin K reversal"),
      fluidRow(
        column(6, numericInput("vitk", "IV dose (mg)", value = 0, min = 0, max = 20, step = 1)),
        column(6, numericInput("vitk_day", "On day", value = 10, min = 0, max = 60, step = 1))
      ),
      tags$h4("Simulation"),
      muted("Each INR point is a separate stiff simulation of the clotting test; in the browser a run takes about a minute."),
      fluidRow(
        column(6, numericInput("sim_days", "Length (days)", value = 21, min = 3, max = 60, step = 1)),
        column(6, numericInput("every", "INR every (days)", value = 3, min = 1, max = 7, step = 1))
      )
    )
  ),
  dashboardBody(
    custom_css,
    tabItems(
      tabItem(
        tabName = "main",
        fluidRow(
          valueBoxOutput("inr_end", width = 3), valueBoxOutput("inr_max", width = 3),
          valueBoxOutput("t_2", width = 3), valueBoxOutput("in_range", width = 3)
        ),
        fluidRow(
          box(title = tags$div(tags$strong("INR"), tags$span(textOutput("engine_note", inline = TRUE),
                               style = "color: #6b7078; font-size: 13px; font-weight: normal;")),
              status = "primary", solidHeader = TRUE, width = 7,
              plotOutput("inr_plot", height = "380px"),
              muted("Each point is a simulated prothrombin-time test: the plasma at that time is diluted 1:3, tissue factor is added, and the clotting time is when the fibrin integral reaches 1500 nM*s (Wajima 2009). INR = PT / PT at baseline.")),
          box(title = "Warfarin in plasma", status = "primary", solidHeader = TRUE, width = 5,
              plotOutput("warf_plot", height = "380px"))
        ),
        fluidRow(
          box(title = "Vitamin K-dependent factors", status = "primary", solidHeader = TRUE, width = 7,
              plotOutput("factor_plot", height = "360px"),
              muted("Each factor falls at its own degradation rate: factor VII (half-life ~6 h) first, prothrombin (~70 h) last. That lag is why INR rises for days after the first dose.")),
          box(title = "Clotting test results", status = "info", solidHeader = TRUE, width = 5,
              DT::dataTableOutput("inr_table"))
        )
      ),
      about_tab(
        "Coagulation network with warfarin",
        "A quantitative systems pharmacology model of human blood coagulation. It links
         warfarin pharmacokinetics to the vitamin K cycle, the synthesis of the vitamin
         K-dependent clotting factors, the full extrinsic, intrinsic and common pathways,
         fibrin formation and fibrinolysis, and finally to the clinical read-out, the
         prothrombin time / INR.",
        list("54 species and 116 reactions: clotting factors, their activated forms and complexes, inhibitors (antithrombin, TFPI, protein C/S), fibrin(ogen), plasmin, D-dimer",
             "Vitamin K cycle (vitamin K, hydroquinone, epoxide) driving synthesis of factors II, VII, IX, X, protein C and protein S",
             "Warfarin: one-compartment oral PK; inhibition of vitamin K epoxide reductase (Imax 1, IC50 0.34 mg/L)",
             "In-vitro prothrombin time test simulated on the plasma composition at each sampling time"),
        tags$span("Wajima T, Isbister GK, Duffull SB. A comprehensive model for the humoral
                   coagulation network in humans. Clin Pharmacol Ther 2009;86:290-298.
                   Converted from BioModels BIOMD0000000340 (in vivo) and BIOMD0000000339
                   (PT test) by tools/build_coagulation.py."),
        M$engine$name
      )
    )
  )
)

server <- function(input, output, session) {

  inputs <- reactive({
    shiny::req(input$warf >= 0, input$days >= 1, input$sim_days >= 3, input$every >= 1)
    list(warf = input$warf, days = input$days, cl = input$cl_mult, ic50 = input$ic50_mult,
         vitk = input$vitk, vitk_day = input$vitk_day, sim_days = input$sim_days, every = input$every)
  }) |> debounce(500)

  sim <- reactive({
    x <- inputs()
    P <- data.frame(warfarin_ke = M$model$param[["warfarin_ke"]] * x$cl,
                    IC50 = M$model$param[["IC50"]] * x$ic50)
    ev <- coag_events(x$warf, x$days, x$vitk, x$vitk_day)
    times <- seq(0, x$sim_days * 24, by = 6)
    withProgress(message = "Simulating the coagulation network", value = 0.2, {
      r <- M$engine$solve(M$model, P, ev, times, rtol = 1e-5, atol = 1e-8, nonneg = TRUE)
      days <- seq(0, x$sim_days, by = x$every)
      idx <- match(days * 24, times)
      incProgress(0.3, message = "Running prothrombin-time tests")
      if (isTRUE(M$engine$mrgsolve)) {
        pt <- vapply(idx, function(i) coag_pt_mrgsolve(M$model, r$states[i, , 1]), 0)
      } else {
        prep <- mrg_prepare(M$model, data.frame(INVITRO = 1), nonneg = TRUE)
        pt <- vapply(idx, function(i) coag_pt(prep, r$states[i, , 1]), 0)
      }
      list(r = r, times = times, days = days, pt = pt, inr = pt / pt[1], x = x)
    })
  })

  output$engine_note <- renderText(sprintf("  (%s)", M$engine$name))

  output$inr_end <- renderValueBox({
    s <- sim(); d <- min(s$x$days, max(s$days))
    v <- stats::approx(s$days, s$inr, xout = d)$y
    stat_box(sprintf("%.2f", v), "", sprintf("INR on day %d", d), "End of warfarin course", "vial")
  })
  output$inr_max <- renderValueBox({
    s <- sim()
    stat_box(sprintf("%.2f", max(s$inr, na.rm = TRUE)), "", "Peak INR",
             sprintf("Day %s", plain_number(s$days[which.max(s$inr)])), "arrow-up",
             color = if (max(s$inr, na.rm = TRUE) > 4) "red" else "blue")
  })
  output$t_2 <- renderValueBox({
    s <- sim(); k <- which(s$inr >= 2)[1]
    stat_box(if (is.na(k)) "-" else plain_number(s$days[k]), if (is.na(k)) "" else " days",
             "Time to INR >= 2", "First sampled day in range", "clock")
  })
  output$in_range <- renderValueBox({
    s <- sim()
    stat_box(sprintf("%d / %d", sum(s$inr >= 2 & s$inr <= 3, na.rm = TRUE), length(s$inr)), "",
             "Tests with INR 2-3", "Typical therapeutic range", "bullseye")
  })

  output$inr_plot <- renderPlot({
    s <- sim()
    d <- data.frame(day = s$days, inr = s$inr)
    ggplot(d, aes(day, inr)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2, ymax = 3, fill = PAL$sunken) +
      annotate("text", x = max(d$day), y = 3, label = "Therapeutic range 2-3", hjust = 1, vjust = -0.5,
               colour = PAL$ink_3, size = 3.5) +
      geom_vline(xintercept = s$x$days, colour = PAL$ink_3, linetype = "22") +
      geom_line(colour = PAL$blue_ink, linewidth = 1.1) +
      geom_point(colour = PAL$blue_ink, size = 2.6) +
      scale_y_continuous(limits = c(0.8, max(3.5, max(d$inr, na.rm = TRUE) * 1.08))) +
      labs(x = "Day", y = "INR", subtitle = sprintf("Dashed line: last warfarin dose (day %d)", s$x$days)) +
      theme_sim(12)
  })

  output$warf_plot <- renderPlot({
    s <- sim()
    d <- data.frame(time = s$times / 24, value = s$r$states[, "C_warf", 1])
    ggplot(d, aes(time, value)) + geom_line(colour = PAL$blue_ink, linewidth = 1) +
      labs(x = "Day", y = "Warfarin (mg/L)") + theme_sim(12)
  })

  output$factor_plot <- renderPlot({
    s <- sim()
    d <- do.call(rbind, lapply(names(FACTORS), function(nm) {
      v <- s$r$states[, FACTORS[[nm]], 1]
      data.frame(time = s$times / 24, value = 100 * v / v[1], series = nm)
    }))
    series_plot(d, "Day", "% of baseline", names(FACTORS)) +
      scale_y_continuous(limits = c(0, NA))
  })

  output$inr_table <- DT::renderDataTable({
    s <- sim()
    idx <- match(s$days * 24, s$times)
    st <- s$r$states[idx, , 1]
    small_table(data.frame(Day = s$days, `PT (s)` = sprintf("%.1f", s$pt), INR = sprintf("%.2f", s$inr),
                           `Factor II %` = sprintf("%.0f", 100 * st[, "FII"] / st[1, "FII"]),
                           `Factor VII %` = sprintf("%.0f", 100 * st[, "VII"] / st[1, "VII"]),
                           check.names = FALSE), page = 12)
  })
}

shinyApp(ui, server)
