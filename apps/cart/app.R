# ============================================================================
# CAR-T cellular kinetics: tisagenlecleucel (Stein et al. 2019)
#
# CAR-T cells are a living drug. After a single infusion they expand
# several-thousand-fold, peak at about day 9-10, contract as activated
# effector cells die, and a small memory-like fraction persists for years.
# The population model from the pediatric/young-adult ALL trials (ELIANA,
# ENSIGN) with its between-patient variability, covariates on Cmax and the
# effect of tocilizumab and corticosteroids on expansion.
# ============================================================================

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(scales)

shared <- function(f) source(if (file.exists(file.path("R", f))) file.path("R", f) else file.path("..", "..", "shared", f))
for (f in c("theme.R", "ode_engine.R", "app_helpers.R")) shared(f)

M <- load_model("cart_tisagenlecleucel_stein2019.cpp")

# Table 1: population estimates, log-normal between-patient SDs and
# covariate effects on log Cmax
TV <- c(FOLDX = 3900, TMAX = 9.3, CMAX = 24000, ALPHA = 0.16, FB = 0.0079, BETA = 0.0032)
OMEGA <- c(FOLDX = 2.4, TMAX = 0.38, CMAX = 0.65, ALPHA = 0.91, FB = 0.8, BETA = 0.86)
COV <- c(female = 0.25, asian = 0.13, other = 0.33, downs = 0.25, hsct = 0.29, nofluda = -0.63,
         b2205j = -0.11, te = 0.22, dose = 0.093, toci = 0.44, ster = -0.36)
MEDIAN_TE <- 19      # transduction efficiency, %
MEDIAN_DOSE <- 3.1   # 10^6 cells/kg

log_axis <- scale_y_log10(labels = label_number(scale_cut = cut_short_scale(), drop0trailing = TRUE))

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = tags$div(tags$span("QSP", style = "font-weight: bold; font-size: 24px;"),
                                   tags$span(" CAR-T cell kinetics", style = "font-size: 18px;")),
                  titleWidth = 350),
  dashboardSidebar(
    width = 320, custom_css,
    sidebarMenu(
      menuItem("Virtual patients", tabName = "main", icon = icon("users")),
      menuItem("Covariates", tabName = "covs", icon = icon("sliders-h")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    tags$div(
      style = "padding: 18px; padding-top: 8px;",
      tags$h4("Product and patient"),
      fluidRow(
        column(6, numericInput("dose", "Dose (10^6 cells/kg)", MEDIAN_DOSE, min = 0.1, max = 10, step = 0.5)),
        column(6, numericInput("te", "Transduction eff. (%)", MEDIAN_TE, min = 1, max = 80, step = 1))
      ),
      fluidRow(
        column(6, radioButtons("sex", NULL, c("Male" = "0", "Female" = "1"), inline = TRUE)),
        column(6, selectInput("race", NULL, c("White" = "white", "Asian" = "asian", "Other/unknown" = "other")))
      ),
      checkboxInput("fluda", "Fludarabine lymphodepletion", TRUE),
      checkboxInput("hsct", "Previous stem cell transplant", FALSE),
      checkboxInput("downs", "Down syndrome", FALSE),
      tags$h4("Cytokine release syndrome treatment"),
      fluidRow(
        column(6, checkboxInput("toci", "Tocilizumab", FALSE)),
        column(6, numericInput("toci_day", "from day", 5.7, min = 0, max = 30, step = 0.5))
      ),
      fluidRow(
        column(6, checkboxInput("ster", "Corticosteroids", FALSE)),
        column(6, numericInput("ster_day", "from day", 7.5, min = 0, max = 30, step = 0.5))
      ),
      tags$h4("Simulation"),
      fluidRow(
        column(6, numericInput("n", "Virtual patients", 100, min = 10, max = 300, step = 10)),
        column(6, selectInput("horizon", "Follow-up", c("28 days" = "28", "3 months" = "91", "1 year" = "365", "2 years" = "730"),
                              selected = "365"))
      ),
      checkboxInput("iiv", "Between-patient variability", TRUE),
      numericInput("seed", "Random seed", 2019, min = 1, step = 1, width = "50%")
    )
  ),
  dashboardBody(
    custom_css,
    tabItems(
      tabItem(
        tabName = "main",
        fluidRow(valueBoxOutput("cmax", width = 3), valueBoxOutput("tmax", width = 3),
                 valueBoxOutput("auc", width = 3), valueBoxOutput("persist", width = 3)),
        fluidRow(
          box(title = tags$div(tags$strong("Transgene in peripheral blood"),
                               tags$span(textOutput("engine_note", inline = TRUE),
                                         style = "color: #6b7078; font-size: 13px; font-weight: normal;")),
              status = "primary", solidHeader = TRUE, width = 7, plotOutput("pop_plot", height = "370px"),
              muted("Median with 50% and 90% ranges across the virtual patients. Expansion, contraction and persistence vary more between patients than almost any small-molecule PK: the fold expansion alone has a log-SD of 2.4.")),
          box(title = "The typical patient: two populations of cells", status = "primary", solidHeader = TRUE, width = 5,
              plotOutput("phase_plot", height = "370px"),
              muted("Both populations expand together at the same rate. After Tmax, 99% contract with a half-life of about 4 days; the rest persist with a half-life of about 7 months."))
        ),
        fluidRow(
          box(title = "Tocilizumab and corticosteroids", status = "primary", solidHeader = TRUE, width = 12,
              plotOutput("comed_plot", height = "300px"),
              muted("The typical patient with tocilizumab or steroids started during expansion. The rate factors were estimated at 1.2 (tocilizumab) and 1.0 (steroids), so neither slowed expansion. The Cmax covariates are separate: patients who needed tocilizumab had higher peaks (effect +0.44 on log Cmax), most likely because high expansion causes the cytokine release syndrome that tocilizumab treats, not the other way round."))
        )
      ),
      tabItem(
        tabName = "covs",
        fluidRow(
          box(title = "What moves Cmax: covariate effects against between-patient variability", status = "primary",
              solidHeader = TRUE, width = 12, plotOutput("cov_plot", height = "430px"),
              muted("Bars: the estimated effect of each covariate on Cmax (Table 1), as a fold change for a patient who has it. Continuous covariates are shown at the ends of the range studied (Table 2). Shaded band: the middle 90% of Cmax from between-patient variability alone (log-SD 0.65). Every covariate effect sits inside the unexplained variability, and most were estimated imprecisely (relative standard errors of 59-250%): across the 27-fold range of doses studied, dose per kg moves Cmax by less than 25%."))
        )
      ),
      about_tab(
        "Tisagenlecleucel cellular kinetics",
        "A population model of CAR-T cell kinetics fitted to transgene levels (qPCR, copies per
         ug of genomic DNA) from 90 children and young adults with relapsed or refractory B-cell
         acute lymphoblastic leukemia treated with tisagenlecleucel. It describes the three phases
         that set cell therapies apart from classical drugs: expansion, contraction and persistence.",
        list("Expansion: exponential growth at rate rho = ln(fold expansion) / Tmax until Tmax",
             "Contraction: most cells (1 - FB) decline at rate alpha after Tmax",
             "Persistence: a fraction FB declines at the much slower rate beta (memory-like cells)",
             "Tocilizumab or corticosteroids started before Tmax multiply the expansion rate by Ftoci or Fster",
             "Between-patient variability: log-normal on fold expansion, Tmax, Cmax, alpha, FB and beta",
             "Covariates on Cmax: sex, race, Down syndrome, previous transplant, fludarabine, study, transduction efficiency, dose per kg, tocilizumab and steroid use"),
        tags$span("Stein AM, Grupp SA, Levine JE, Laetsch TW, Pulsipher MA, Boyer MW, et al. Tisagenlecleucel
                   model-based cellular kinetic analysis of chimeric antigen receptor-T cells. CPT
                   Pharmacometrics Syst Pharmacol 2019;8:285-295. Parameters from Table 1; comedication
                   model from the MLXTRAN code in the Supplementary Material, written here as differential
                   equations for the log levels."),
        M$engine$name
      )
    )
  )
)

server <- function(input, output, session) {

  setup <- reactive({
    shiny::req(input$dose > 0, input$te > 0, input$n >= 10, input$toci_day >= 0, input$ster_day >= 0)
    n <- min(round(input$n), 300)
    set.seed(input$seed)
    eta <- if (input$iiv) sapply(names(TV), function(k) stats::rnorm(n, 0, OMEGA[[k]])) else matrix(0, n, 6, dimnames = list(NULL, names(TV)))
    if (n == 1) eta <- matrix(eta, 1, dimnames = list(NULL, names(TV)))
    toci <- input$toci; ster <- input$ster
    cov_cmax <- exp(COV[["female"]] * (input$sex == "1") + COV[["asian"]] * (input$race == "asian") +
                    COV[["other"]] * (input$race == "other") + COV[["downs"]] * input$downs +
                    COV[["hsct"]] * input$hsct + COV[["nofluda"]] * !input$fluda +
                    COV[["toci"]] * toci + COV[["ster"]] * ster) *
                (input$te / MEDIAN_TE)^COV[["te"]] * (input$dose / MEDIAN_DOSE)^COV[["dose"]]
    P <- as.data.frame(lapply(names(TV), function(k) TV[[k]] * exp(eta[, k])))
    names(P) <- names(TV)
    P$CMAX <- P$CMAX * cov_cmax
    P$FB <- pmin(P$FB, 0.95)
    P$TTOCI <- if (toci) input$toci_day else 99999
    P$TSTER <- if (ster) input$ster_day else 99999
    list(P = P, horizon = as.numeric(input$horizon), cov_cmax = cov_cmax)
  }) |> debounce(500)

  sim_times <- function(horizon) sort(unique(c(seq(0, min(28, horizon), by = 0.25), seq(28, horizon, length.out = 150))))

  sim <- reactive({
    s <- setup(); tt <- sim_times(s$horizon)
    withProgress(message = "Simulating CAR-T kinetics", value = 0.3, {
      r <- M$engine$solve(M$model, s$P, times = tt, rtol = 1e-6, atol = 1e-8)
    })
    list(s = s, t = tt, r = r)
  })

  output$engine_note <- renderText(sprintf("  (%s)", M$engine$name))

  qfmt <- function(x) sprintf("%s-%s", label_number(scale_cut = cut_short_scale(), accuracy = 1)(stats::quantile(x, 0.05)),
                              label_number(scale_cut = cut_short_scale(), accuracy = 1)(stats::quantile(x, 0.95)))
  per_patient <- reactive({
    x <- sim(); w <- x$t <= 28
    cmax <- apply(x$r$CART, 2, max)
    tmax <- x$t[apply(x$r$CART, 2, which.max)]
    tw <- x$t[w]; y <- x$r$CART[w, , drop = FALSE]
    auc <- colSums(diff(tw) * (y[-1, , drop = FALSE] + y[-nrow(y), , drop = FALSE]) / 2)
    last <- x$r$CART[length(x$t), ]
    list(cmax = cmax, tmax = tmax, auc = auc, last = last, h = x$s$horizon)
  })

  output$cmax <- renderValueBox({
    p <- per_patient()
    stat_box(label_number(scale_cut = cut_short_scale(), accuracy = 0.1)(stats::median(p$cmax)), " copies/ug",
             "Cmax, median", if (length(p$cmax) > 1) paste("90% of patients:", qfmt(p$cmax)) else "", "arrow-up")
  })
  output$tmax <- renderValueBox({
    p <- per_patient()
    stat_box(sprintf("%.1f", stats::median(p$tmax)), " days", "Tmax, median",
             sprintf("90%% of patients: %.0f-%.0f days", stats::quantile(p$tmax, 0.05), stats::quantile(p$tmax, 0.95)), "clock")
  })
  output$auc <- renderValueBox({
    p <- per_patient()
    stat_box(label_number(scale_cut = cut_short_scale(), accuracy = 0.1)(stats::median(p$auc)), " copies/ug*d",
             "AUC day 0-28, median", "The exposure metric linked to response", "chart-area")
  })
  output$persist <- renderValueBox({
    p <- per_patient()
    stat_box(label_number(accuracy = 1, big.mark = ",")(stats::median(p$last)), " copies/ug",
             sprintf("Transgene at day %d", p$h), sprintf("%.0f%% of patients above 100 copies/ug", 100 * mean(p$last > 100)), "hourglass-half")
  })

  output$pop_plot <- renderPlot({
    x <- sim()
    d <- summarise_bands(x$t, x$r$CART)
    band_plot(d, "Days after infusion", "Transgene (copies/ug DNA)") + log_axis
  })

  typical <- function(toci_day = 99999, ster_day = 99999, cmax_factor = 1, horizon = 91) {
    P <- as.data.frame(as.list(TV)); P$CMAX <- P$CMAX * cmax_factor
    P$TTOCI <- toci_day; P$TSTER <- ster_day
    tt <- sim_times(horizon)
    list(t = tt, r = M$engine$solve(M$model, P, times = tt, rtol = 1e-6, atol = 1e-8))
  }

  output$phase_plot <- renderPlot({
    s <- setup(); y <- typical(horizon = s$horizon)
    lv <- c("Total", "Contracting (effector)", "Persisting (memory-like)")
    d <- rbind(data.frame(time = y$t, value = y$r$CART[, 1], series = lv[1]),
               data.frame(time = y$t, value = y$r$EFF[, 1], series = lv[2]),
               data.frame(time = y$t, value = y$r$PER[, 1], series = lv[3]))
    d <- d[d$value > 0.1, ]
    series_plot(d, "Days after infusion", "Transgene (copies/ug DNA)", lv) + log_axis
  })

  output$comed_plot <- renderPlot({
    td <- input$toci_day; sd <- input$ster_day
    shiny::req(td >= 0, sd >= 0)
    lv <- c("Neither", sprintf("Tocilizumab from day %s", plain_number(td)), sprintf("Steroids from day %s", plain_number(sd)),
            "Tocilizumab: rate effect only")
    runs <- list(typical(horizon = 28), typical(toci_day = td, cmax_factor = exp(COV[["toci"]]), horizon = 28),
                 typical(ster_day = sd, cmax_factor = exp(COV[["ster"]]), horizon = 28), typical(toci_day = td, horizon = 28))
    d <- do.call(rbind, lapply(seq_along(runs), function(j) data.frame(time = runs[[j]]$t, value = runs[[j]]$r$CART[, 1], series = lv[j])))
    series_plot(d, "Days after infusion", "Transgene (copies/ug DNA)", lv) + log_axis +
      labs(subtitle = "Typical patient; the first three include the Cmax covariate for patients who received the drug")
  })

  # ---- Covariates ---------------------------------------------------------
  output$cov_plot <- renderPlot({
    d <- data.frame(
      label = c("Female (vs male)", "Asian (vs white)", "Other/unknown race (vs white)", "Down syndrome",
                "Previous stem cell transplant", "No fludarabine lymphodepletion", "ENSIGN (vs ELIANA)",
                "Transduction efficiency 2.3% (vs 19%)", "Transduction efficiency 56% (vs 19%)",
                "Dose 0.2 x 10^6/kg (vs 3.1)", "Dose 5.4 x 10^6/kg (vs 3.1)", "Received tocilizumab",
                "Received corticosteroids"),
      fold = c(exp(COV[c("female", "asian", "other", "downs", "hsct", "nofluda", "b2205j")]),
               (c(2.3, 56) / MEDIAN_TE)^COV[["te"]], (c(0.2, 5.4) / MEDIAN_DOSE)^COV[["dose"]],
               exp(COV[c("toci", "ster")]))
    )
    d$label <- factor(d$label, levels = rev(d$label))
    iiv <- exp(c(-1, 1) * stats::qnorm(0.95) * OMEGA[["CMAX"]])
    ggplot(d, aes(fold, label)) +
      annotate("rect", xmin = iiv[1], xmax = iiv[2], ymin = -Inf, ymax = Inf, fill = PAL$blue, alpha = 0.12) +
      geom_vline(xintercept = 1, colour = PAL$ink_3) +
      geom_segment(aes(x = 1, xend = fold, yend = label), colour = PAL$blue_ink, linewidth = 2.5) +
      scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.25, 4)) +
      labs(x = "Fold change in Cmax (log scale)", y = NULL) + theme_sim(12)
  })
}

shinyApp(ui, server)
