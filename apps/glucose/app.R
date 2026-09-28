# ============================================================================
# Glucose, insulin and the incretins: the 4GI model (Bosch et al. 2022)
#
# Glucose, insulin, GLP-1, glucagon and GIP with the feedbacks between them:
# glucose drives insulin, GLP-1 and GIP amplify glucose-dependent insulin
# secretion, GLP-1 slows gastric emptying and suppresses glucagon, and
# glucagon drives hepatic glucose output. Liraglutide, a GLP-1 receptor
# agonist, acts on the same receptors as endogenous GLP-1.
# ============================================================================

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(scales)

shared <- function(f) source(if (file.exists(file.path("R", f))) file.path("R", f) else file.path("..", "..", "shared", f))
for (f in c("theme.R", "ode_engine.R", "app_helpers.R")) shared(f)

M <- load_model("glucose_4gi_bosch2022.cpp")

LIRA_PMOL_PER_MG <- 1e9 / 3751.2
MMOL_PER_G_GLUCOSE <- 1000 / 180.16

HORMONES <- c(Glucose = "GLC", Insulin = "INS", "GLP-1" = "GLP1", Glucagon = "GLG", GIP = "GIP")
UNITS <- c(GLC = "mM", INS = "pM", GLP1 = "pM", GLG = "pM", GIP = "pM")

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = tags$div(tags$span("QSP", style = "font-weight: bold; font-size: 24px;"),
                                   tags$span(" Glucose & incretins", style = "font-size: 18px;")),
                  titleWidth = 350),
  dashboardSidebar(
    width = 320, custom_css,
    sidebarMenu(
      menuItem("Meal response", tabName = "main", icon = icon("utensils")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    tags$div(
      style = "padding: 18px; padding-top: 8px;",
      tags$h4("Population"),
      radioButtons("pat", NULL, choices = c("Type 2 diabetes" = "0", "Healthy" = "1"), selected = "0", inline = TRUE),
      fluidRow(
        column(6, numericInput("bsl_glc", "Fasting glucose (mM)", 7.8, min = 3.5, max = 15, step = 0.1)),
        column(6, numericInput("bsl_ins", "Fasting insulin (pM)", 52, min = 10, max = 300, step = 5))
      ),
      tags$h4("Test meal"),
      radioButtons("meal", NULL, choices = c("75 g oral glucose (OGTT)" = "ogtt", "Three meals a day" = "meals"),
                   selected = "meals"),
      conditionalPanel(
        "input.meal == 'meals'",
        fluidRow(
          column(4, numericInput("g_b", "Breakfast g", 50, min = 0, max = 150, step = 5)),
          column(4, numericInput("g_l", "Lunch g", 70, min = 0, max = 150, step = 5)),
          column(4, numericInput("g_d", "Dinner g", 80, min = 0, max = 150, step = 5))
        ),
        muted("Carbohydrate per meal. Glucose appearance uses the meal bioavailabilities the model estimated for mixed meals (Jauslin study).")
      ),
      tags$h4("Liraglutide"),
      radioButtons("lira", NULL, choices = c("None" = "0", "0.6 mg" = "0.6", "1.2 mg" = "1.2", "1.8 mg" = "1.8"),
                   selected = "1.8", inline = TRUE),
      numericInput("lira_days", "Daily SC injections before the test day", 13, min = 0, max = 28, step = 1, width = "100%"),
      numericInput("bw", "Body weight (kg)", 90, min = 40, max = 200, step = 5, width = "100%")
    )
  ),
  dashboardBody(
    custom_css,
    tabItems(
      tabItem(
        tabName = "main",
        fluidRow(valueBoxOutput("mean_glc", width = 3), valueBoxOutput("peak_glc", width = 3),
                 valueBoxOutput("auc_glc", width = 3), valueBoxOutput("peak_ins", width = 3)),
        fluidRow(
          box(title = tags$div(tags$strong("Plasma glucose on the test day"),
                               tags$span(textOutput("engine_note", inline = TRUE),
                                         style = "color: #6b7078; font-size: 13px; font-weight: normal;")),
              status = "primary", solidHeader = TRUE, width = 7, plotOutput("glc_plot", height = "360px")),
          box(title = "Liraglutide", status = "primary", solidHeader = TRUE, width = 5,
              plotOutput("lira_plot", height = "360px"),
              muted("Watson 2010 PK: SC absorption half-life ~4.5 h, elimination half-life ~13 h; steady state after 3-4 days."))
        ),
        fluidRow(
          box(title = "Hormones on the test day, with and without liraglutide", status = "primary",
              solidHeader = TRUE, width = 12, plotOutput("hormones", height = "330px"),
              muted("Blue: with the selected liraglutide regimen. Orange: the same meals without drug. Liraglutide raises glucose-dependent insulin secretion, lowers glucagon and slows gastric emptying, which flattens the endogenous GLP-1 and GIP responses."))
        )
      ),
      about_tab(
        "4GI: glucose, insulin, GLP-1, glucagon and GIP",
        "An integrated quantitative systems pharmacology model of glucose regulation,
         built to support development of a glucagon/GLP-1 dual agonist. It was fitted to
         published intravenous and oral glucose tests, meal tests and hormone infusions
         in healthy volunteers and people with type 2 diabetes, and validated externally.",
        list("Glucose: gut transit and absorption, central and peripheral pools, insulin-dependent and insulin-independent disposal, glucagon-driven hepatic output",
             "Insulin: glucose-dependent secretion amplified by GLP-1 and GIP (incretin effect), effect compartment for insulin action",
             "GLP-1: food-stimulated secretion, saturable elimination (DPP-4); slows glucose absorption, stimulates insulin, suppresses glucagon",
             "Glucagon: glucose feedback (switched off below baseline in T2DM), suppressed by GLP-1, stimulated by GIP",
             "GIP: food-stimulated secretion, two-compartment kinetics; insulinotropic only in healthy volunteers",
             "Liraglutide: one-compartment SC PK (Watson 2010); free fraction 0.5%, acts on the GLP-1 receptor with its in-vitro potency"),
        tags$span("Bosch R, Petrone M, Arends R, Vicini P, Sijbrands EJG, Hoefman S, Snelder N. A novel
                   integrated QSP model of in vivo human glucose regulation to support the development
                   of a glucagon/GLP-1 dual agonist. CPT Pharmacometrics Syst Pharmacol 2022;11:302-317.
                   Transcribed from the published NONMEM control stream (Supplementary Material 2) with
                   the final parameter estimates."),
        M$engine$name
      )
    )
  )
)

server <- function(input, output, session) {

  observeEvent(input$pat, {
    if (input$pat == "1") {
      updateNumericInput(session, "bsl_glc", value = 4.65); updateNumericInput(session, "bsl_ins", value = 49.1)
    } else {
      updateNumericInput(session, "bsl_glc", value = 7.8); updateNumericInput(session, "bsl_ins", value = 52)
    }
  }, ignoreInit = TRUE)

  setup <- reactive({
    shiny::req(input$bsl_glc > 0, input$bsl_ins > 0, input$bw > 0, input$lira_days >= 0)
    healthy <- input$pat == "1"
    P <- data.frame(PAT = as.numeric(healthy), BSLglc = input$bsl_glc, BSLins = input$bsl_ins, BW = input$bw,
                    BSLglg = if (healthy) 8.85 else 28.2)
    day <- input$lira_days                     # test day index (0-based)
    t0 <- day * 24
    meals <- if (input$meal == "ogtt") {
      data.frame(time = t0 + 8, cmt = "Dglc", amt = 75 * MMOL_PER_G_GLUCOSE * 0.776)
    } else {
      data.frame(time = t0 + c(8, 13, 19), cmt = "Dglc",
                 amt = c(input$g_b, input$g_l, input$g_d) * MMOL_PER_G_GLUCOSE * c(0.446, 0.294, 0.287))
    }
    meals <- meals[meals$amt > 0, , drop = FALSE]
    dose <- as.numeric(input$lira)
    lira <- if (dose > 0) data.frame(time = seq(0, t0, by = 24), cmt = "Ddrug", amt = dose * LIRA_PMOL_PER_MG) else NULL
    list(P = P, meals = meals, lira = lira, t0 = t0, dose = dose)
  }) |> debounce(500)

  sim <- reactive({
    s <- setup()
    times <- sort(unique(c(seq(0, s$t0, by = 2), s$t0 + seq(0, 24, by = 0.1))))
    withProgress(message = "Simulating glucose regulation", value = 0.3, {
      with_drug <- M$engine$solve(M$model, s$P, rbind(s$lira, s$meals), times, rtol = 1e-5, atol = 1e-8)
      no_drug <- if (is.null(s$lira)) with_drug else M$engine$solve(M$model, s$P, s$meals, times, rtol = 1e-5, atol = 1e-8)
    })
    w <- times >= s$t0
    list(s = s, times = times, w = w, a = with_drug, b = no_drug)
  })

  output$engine_note <- renderText(sprintf("  (%s)", M$engine$name))

  day_stats <- function(x, r) {
    g <- r$GLC[x$w, 1]; tt <- x$times[x$w]
    list(mean = mean(g), peak = max(g), auc = sum(diff(tt) * (utils::head(g, -1) + utils::tail(g, -1)) / 2),
         ins = max(r$INS[x$w, 1]))
  }
  delta <- function(a, b, fmt) if (identical(a, b)) "" else sprintf(fmt, a - b)

  output$mean_glc <- renderValueBox({
    x <- sim(); a <- day_stats(x, x$a); b <- day_stats(x, x$b)
    stat_box(sprintf("%.1f", a$mean), " mM", "Mean glucose, test day",
             if (x$s$dose > 0) sprintf("%+.1f mM vs no drug", a$mean - b$mean) else "No drug", "droplet")
  })
  output$peak_glc <- renderValueBox({
    x <- sim(); a <- day_stats(x, x$a); b <- day_stats(x, x$b)
    stat_box(sprintf("%.1f", a$peak), " mM", "Peak glucose",
             if (x$s$dose > 0) sprintf("%+.1f mM vs no drug", a$peak - b$peak) else "No drug", "arrow-up")
  })
  output$auc_glc <- renderValueBox({
    x <- sim(); a <- day_stats(x, x$a); b <- day_stats(x, x$b)
    stat_box(sprintf("%.0f", a$auc), " mM*h", "Glucose AUC, 24 h",
             if (x$s$dose > 0) sprintf("%+.0f%% vs no drug", 100 * (a$auc / b$auc - 1)) else "No drug", "chart-area")
  })
  output$peak_ins <- renderValueBox({
    x <- sim(); a <- day_stats(x, x$a)
    stat_box(sprintf("%.0f", a$ins), " pM", "Peak insulin", "Test day", "syringe")
  })

  test_day <- function(x, r, var) data.frame(time = x$times[x$w] - x$s$t0, value = r[[var]][x$w, 1])

  output$glc_plot <- renderPlot({
    x <- sim()
    lv <- if (x$s$dose > 0) c(sprintf("Liraglutide %s mg", plain_number(x$s$dose)), "No drug") else "No drug"
    d <- rbind(transform(test_day(x, x$a, "GLC"), series = lv[1]),
               if (x$s$dose > 0) transform(test_day(x, x$b, "GLC"), series = lv[2]))
    series_plot(d, "Hour of the test day", "Glucose (mM)", lv) +
      geom_hline(yintercept = c(3.9, 10), colour = PAL$ink_3, linetype = "22") +
      scale_x_continuous(breaks = seq(0, 24, 4)) +
      labs(subtitle = "Dashed lines: 3.9 and 10 mM, the usual target range")
  })

  output$lira_plot <- renderPlot({
    x <- sim()
    d <- data.frame(time = x$times / 24, value = x$a$LIRA[, 1] / 1000)
    ggplot(d, aes(time, value)) + geom_line(colour = PAL$blue_ink, linewidth = 1) +
      labs(x = "Day", y = "Liraglutide, total (nM)") + theme_sim(12)
  })

  output$hormones <- renderPlot({
    x <- sim()
    d <- do.call(rbind, lapply(names(HORMONES)[-1], function(nm) {
      v <- HORMONES[[nm]]
      rbind(transform(test_day(x, x$a, v), series = "With drug", panel = sprintf("%s (%s)", nm, UNITS[[v]])),
            if (x$s$dose > 0) transform(test_day(x, x$b, v), series = "No drug", panel = sprintf("%s (%s)", nm, UNITS[[v]])))
    }))
    lv <- c("With drug", "No drug")[c(TRUE, x$s$dose > 0)]
    series_plot(d, "Hour of the test day", NULL, lv) +
      facet_wrap(~panel, nrow = 1, scales = "free_y") +
      scale_x_continuous(breaks = seq(0, 24, 8)) +
      theme(strip.text = element_text(colour = PAL$ink, face = "bold", hjust = 0))
  })
}

shinyApp(ui, server)
