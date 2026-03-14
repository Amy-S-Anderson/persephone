# =============================================================================
# Persephone Shiny App
# Agent-Based Model for Bioarchaeology
# Based on Anderson & DeWitte, "Known Unknowns and the Osteological Paradox"
# =============================================================================

library(shiny)
remotes::install_github('Amy-S-Anderson/persephone')
library(persephone)
library(ggplot2)
library(dplyr)
library(survival)

# =============================================================================
# Mortality Regimes: named list for UI dropdown
# Data comes from persephone package
# =============================================================================

mortality_regimes <- list(
  "CDW Level 3 (high mortality)" = CoaleDemenyWestF3,
  "CDW Level 5" = CoaleDemenyWestF5,
  "CDW Level 11" = CoaleDemenyWestF11,
  "CDW Level 15" = CoaleDemenyWestF15,
  "CDW Level 17" = CoaleDemenyWestF17,
  "CDW Level 21 (low mortality)" = CoaleDemenyWestF21
)

taphonomy_regimes <- list(
  "no_decay" = no_decay,
  "weak_decay" = weak_decay,
  "moderate_decay" = moderate_decay,
  "strong_decay" = strong_decay
)

# =============================================================================
# Helper: age intervals matching the paper
# =============================================================================

assign_age_interval <- function(age) {
  factor(
    case_when(
      age < 2 ~ "0-1", age < 6 ~ "2-5", age < 10 ~ "6-9",
      age < 15 ~ "10-14", age < 20 ~ "15-19", age < 30 ~ "20-29",
      age < 40 ~ "30-39", age < 50 ~ "40-49", age < 60 ~ "50-59",
      TRUE ~ "60+"
    ),
    levels = c("0-1", "2-5", "6-9", "10-14", "15-19",
               "20-29", "30-39", "40-49", "50-59", "60+")
  )
}

# =============================================================================
# Color palette
# =============================================================================

col_no_lesion <- "#5B8DB8"
col_lesion    <- "#CC6677"
col_dark      <- "#2c3e50"

# =============================================================================
# UI
# =============================================================================

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; }
    h4.param-header {
      color: #495057; border-bottom: 1px solid #dee2e6;
      padding-bottom: 5px; margin-top: 18px; margin-bottom: 10px; font-size: 15px;
    }
    .help-text { font-size: 11px; color: #888; margin-top: -5px; margin-bottom: 8px; }
    .run-btn { margin-top: 15px; width: 100%; font-size: 16px; padding: 10px; }
    .interp { padding: 12px; border-radius: 5px; margin-top: 10px; font-size: 13px; line-height: 1.5; }
    .interp-warn   { background: #f8d7da; border: 1px solid #f5c6cb; }
    .interp-ok     { background: #d4edda; border: 1px solid #c3e6cb; }
    .summary-text  { font-size: 14px; }
    .tab-content   { padding-top: 5px; }
  "))),

  titlePanel(
    div(
      "Persephone",
      tags$small(style = "color: #6c757d; display: block; font-size: 14px; margin-top: 2px;",
                 "An agent-based model for bioarchaeology (Anderson & DeWitte)")
    ),
    windowTitle = "Persephone ABM"
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Population", class = "param-header"),
      numericInput("cohort_size", "Cohort size", 1000, min = 50, max = 10000, step = 50),

      h4("Mortality", class = "param-header"),
      selectInput("mortality_regime", "Mortality regime",
        choices = names(mortality_regimes), selected = "CDW Level 5"),
      selectInput("risk_type", "Lesion-mortality interaction",
        choices = c(
          "Proportional" = "proportional",
          "Decreasing with age" = "time_decreasing",
          "Increasing with age" = "time_increasing"
        )),
      sliderInput("rmr", "Relative mortality risk (RMR)", 1, 5, 1, 0.1),
      div(class = "help-text", "1 = no effect on mortality. 2 = doubles mortality risk."),

      h4("Skeletal Lesions", class = "param-header"),
      sliderInput("lesion_rate", "Annual formation probability", 0.01, 0.20, 0.05, 0.01),
      fluidRow(
        column(6, numericInput("window_opens", "Window opens", 0, min = 0, max = 30)),
        column(6, numericInput("window_closes", "Window closes", 9, min = 0, max = 60))
      ),
      div(class = "help-text", "Age range in which new skeletal lesions can form."),

      h4("Deposition Filters", class = "param-header"),
      sliderInput(
        "remove_children_below",
        "Remove Children below age...",
        min = 0,
        max = 10,
        value = 0,
        step = 1
      ),
      # sliderInput(
      #   "remove_adults_older",
      #   "Remove Adults older than age...",
      #   min = 50,
      #   max = 100,
      #   value = 100,
      #   step = 10
      # ),

      h4("Taphonomic Filter", class = "param-header"),
      selectInput(
        "taphonomic_filter",
        "Taphonomic strength",
        choices = c(
          "No Preservation Loss" = "no_decay",
          "Weak Preservation Loss" = "weak_decay",
          "Moderate Preservation Loss" = "moderate_decay",
          "Strong Preservation Loss" = "strong_decay"
        ),
        selected = "No Preservation Loss"
      ),

      h4("Age Estimation Error", class = "param-header"),
      checkboxInput(
        "age_error_exists",
        "Age Error Exists",
        value = FALSE
      ),

      h4("Simulation", class = "param-header"),
      numericInput("seed", "Random seed (blank = random)", NA, min = 1),
      actionButton("run", "Run Simulation", class = "btn-primary run-btn")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        tabPanel("Set-up",
          br(),
          fluidRow(
            column(12, plotOutput("plot_mort_siler_setup", height = "500px"))
          ),
          fluidRow(
            column(12, plotOutput("plot_death_pmf", height = "400px"))
          ),
          fluidRow(
            column(12, plotOutput("plot_taph_siler_setup", height = "500px"))
          ),
          conditionalPanel(
            condition = "input.age_error_exists == true",
            fluidRow(
              column(12, plotOutput("plot_age_error", height = "500px"))
            )
          )
        ),

        tabPanel("Observed Sample",
          br(),

          # fluidRow(
          #   column(12, plotOutput("plot_alive_observed", height = "400px"))
          # ),

          # fluidRow(
          #   column(12, plotOutput("plot_lesion_pct_observed", height = "400px"))
          # ),

          hr(),

          fluidRow(
            column(12, plotOutput("plot_cem_age_observed", height = "400px"))
          ),

          fluidRow(
            column(12, plotOutput("plot_cem_lesion_observed", height = "400px"))
          ),

          hr(),
          div(class = "summary-text", verbatimTextOutput("cem_summary_observed")),

          hr(),

          fluidRow(
            column(8, plotOutput("plot_km_observed", height = "460px")),
            column(4,
              wellPanel(
                h4("Age Filter", class = "param-header"),
                sliderInput("min_age_observed", "Minimum estimated age to include", 0, 25, 0, 1),
                div(class = "help-text",
                  "Exclude individuals estimated to die within the lesion formation window."
                )
              ),
              wellPanel(
                h4("Log-Rank Test", class = "param-header"),
                verbatimTextOutput("logrank_text_observed")
              ),
              uiOutput("interpretation_observed")
            )
          ),
          fluidRow(
            column(12, plotOutput("plot_true_vs_observed_age", height = "450px"))
          )
        ),

        tabPanel("True Population",
          br(),

          fluidRow(
            column(12, plotOutput("plot_alive", height = "400px"))
          ),

          fluidRow(
            column(12, plotOutput("plot_lesion_pct", height = "400px"))
          ),

          hr(),

          fluidRow(
            column(12, plotOutput("plot_cem_age", height = "400px"))
          ),

          fluidRow(
            column(12, plotOutput("plot_cem_lesion", height = "400px"))
          ),

          hr(),
          div(class = "summary-text", verbatimTextOutput("cem_summary"))
        ),

        # =====================================================================
        # Tab: Survival Analysis
        # =====================================================================
        tabPanel("Survival Analysis",
          br(),
          fluidRow(
            column(8, plotOutput("plot_km", height = "460px")),
            column(4,
              wellPanel(
                h4("Age Filter", class = "param-header"),
                sliderInput("min_age", "Minimum age to include", 0, 25, 0, 1),
                div(class = "help-text",
                  "Exclude individuals who died within the lesion formation window.",
                  "This is the key insight of the paper.")
              ),
              wellPanel(
                h4("Log-Rank Test", class = "param-header"),
                verbatimTextOutput("logrank_text")
              ),
              uiOutput("interpretation")
            )
          )
        ),

        tabPanel("Data",
          br(),
          fluidRow(
            column(6,
              h4("Cemetery (Individual Outcomes)"),
              downloadButton("dl_cemetery", "Download CSV"),
              div(style = "max-height: 500px; overflow-y: auto; margin-top: 10px;",
                  tableOutput("table_cemetery"))
            ),
            column(6,
              h4("Survivors by Age"),
              downloadButton("dl_survivors", "Download CSV"),
              div(style = "max-height: 500px; overflow-y: auto; margin-top: 10px;",
                  tableOutput("table_survivors"))
            )
          )
        )
      )
    )
  )
)

# =============================================================================
# Server
# =============================================================================

server <- function(input, output, session) {

  # -----------------------------------
  # Run simulation on button click
  # -----------------------------------
  sim <- eventReactive(input$run, {
    validate(
      need(input$cohort_size >= 50, "Cohort size must be at least 50."),
      need(input$window_closes >= input$window_opens,
           "Formation window closes must be >= window opens.")
    )

    if (!is.na(input$seed)) set.seed(input$seed)

    withProgress(message = "Simulating cohort...", value = 0.3, {
      result <- Simulate_Cemetery(
        cohort_size = input$cohort_size,
        lesion_formation_rate = input$lesion_rate,
        formation_window_opens = input$window_opens,
        formation_window_closes = input$window_closes,
        mortality_risk_type = input$risk_type,
        relative_mortality_risk = input$rmr,
        mortality_regime = mortality_regimes[[input$mortality_regime]],
        deposition_param = input$remove_children_below, # missing upper cutoff
        taphonomy_regime = taphonomy_regimes[[input$taphonomic_filter]],
        loss_strength = input$taphonomic_filter,
        age_noise = input$age_error_exists
      )
      setProgress(1, message = "Done.")
    })
    result
  })

  # Update survival analysis age filter range when simulation runs
  observeEvent(sim(), {
    new_max <- min(input$window_closes + 15, 40)
    updateSliderInput(session, "min_age", max = new_max)
  })

  # Convenience reactive: cemetery data with age intervals
  cemetery <- reactive({
    req(sim())
    sim()$individual_outcomes %>%
      mutate(
        Age_Interval = assign_age_interval(age),
        Lesion_Status = factor(
          ifelse(lesion == 1, "Present", "Absent"),
          levels = c("Absent", "Present")
        )
      )
  })


  cemetery_observed <- reactive({
    req(sim())
    sim()$individual_outcomes %>%
      mutate(
        Age_Interval = assign_age_interval(estimated_age),
        Lesion_Status = factor(
          ifelse(lesion == 1, "Present", "Absent"),
          levels = c("Absent", "Present")
        )
      )
  })

  filtered_observed <- reactive({
    req(sim())
    d <- sim()$individual_outcomes %>% filter(estimated_age >= input$min_age_observed)
    d$dead <- 1
    d
  })

  # =========================================================================
  # Total Population Tab
  # =========================================================================

  output$plot_cem_age <- renderPlot({
    validate(need(sim(), "Click 'Run Simulation' to begin."))
    ggplot(cemetery(), aes(x = Age_Interval, fill = Lesion_Status)) +
      geom_bar(position = "stack", width = 0.8) +
      scale_fill_manual(values = c("Absent" = col_no_lesion, "Present" = col_lesion)) +
      labs(title = "Age-at-Death Distribution",
           x = "Age at Death", y = "Count", fill = "lesion") +
      theme_bw(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      )
  })

  output$plot_cem_lesion <- renderPlot({
    validate(need(sim(), ""))
    prev <- cemetery() %>%
      group_by(Age_Interval) %>%
      summarise(Pct = sum(lesion) / n() * 100, n = n(), .groups = "drop")

    ggplot(prev, aes(x = Age_Interval, y = Pct)) +
      geom_col(fill = col_lesion, alpha = 0.85, width = 0.8) +
      geom_text(aes(label = paste0("n=", n)), vjust = -0.5, size = 3.2, color = "#555") +
      labs(title = "Lesion Prevalence by Age Group",
           x = "Age at Death", y = "% with Lesions") +
      theme_bw(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      coord_cartesian(ylim = c(0, max(prev$Pct, na.rm = TRUE) * 1.15))
  })

  output$cem_summary <- renderPrint({
    validate(need(sim(), ""))
    d <- sim()$individual_outcomes
    with_les <- sum(d$lesion == 1)
    cat(sprintf("Total individuals: %d\n", nrow(d)))
    cat(sprintf("With lesions: %d (%.1f%%)\n", with_les, with_les / nrow(d) * 100))
    cat(sprintf("Mean Age at death:  %.1f years (all)\n", mean(d$age)))
    if (with_les > 0) {
      cat(sprintf("                    %.1f years (with lesion)\n", mean(d$age[d$lesion == 1])))
    }
    cat(sprintf("                    %.1f years (without lesion)\n", mean(d$age[d$lesion == 0])))
  })

  output$plot_alive <- renderPlot({
    validate(need(sim(), "Click 'Run Simulation' to begin."))
    ggplot(sim()$survivors, aes(x = Age, y = Alive)) +
      geom_area(fill = col_no_lesion, alpha = 0.3) +
      geom_line(color = col_dark, linewidth = 1) +
      labs(title = "Surviving Population Over Time",
           x = "Age (years)", y = "Number Alive") +
      theme_bw(base_size = 13) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  })

  output$plot_lesion_pct <- renderPlot({
    validate(need(sim(), ""))
    surv <- sim()$survivors %>% filter(!is.na(Lesion_perc))
    ggplot(surv, aes(x = Age, y = Lesion_perc)) +
      geom_area(fill = col_lesion, alpha = 0.2) +
      geom_line(color = col_lesion, linewidth = 1) +
      labs(title = "Lesion Prevalence Among Living",
           x = "Age (years)", y = "% with Lesions") +
      theme_bw(base_size = 13) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      coord_cartesian(ylim = c(0, NA))
  })

  filtered <- reactive({
    req(sim())
    d <- sim()$individual_outcomes %>% filter(age >= input$min_age)
    d$dead <- 1
    d
  })

  output$plot_km <- renderPlot({
    validate(need(sim(), "Click 'Run Simulation' to begin."))
    d <- filtered()
    validate(need(nrow(d) > 10, "Too few individuals remain after filtering."))
    validate(need(length(unique(d$lesion)) == 2,
                  "All remaining individuals have the same lesion status."))

    fit <- survfit(Surv(age, dead) ~ lesion, data = d)
    s <- summary(fit)

    surv_df <- data.frame(
      time = s$time,
      survival = s$surv,
      group = ifelse(grepl("=0", as.character(s$strata)), "No Lesion", "Has Lesion")
    )

    # Add starting points so step curves begin at survival = 1
    starts <- surv_df %>%
      group_by(group) %>%
      summarise(time = min(time) - 0.5, .groups = "drop") %>%
      mutate(survival = 1.0, time = pmax(time, 0))
    surv_df <- bind_rows(starts, surv_df) %>% arrange(group, time)

    title <- "Kaplan-Meier Survival Curves"
    if (input$min_age > 0) title <- paste0(title, "  (ages \u2265 ", input$min_age, ")")

    ggplot(surv_df, aes(x = time, y = survival, color = group)) +
      geom_step(linewidth = 1.1) +
      scale_color_manual(values = c("No Lesion" = col_dark, "Has Lesion" = col_lesion)) +
      labs(title = title, x = "Age (years)", y = "Survival Probability", color = "") +
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 12)
      )
  })

  output$logrank_text <- renderPrint({
    validate(need(sim(), ""))
    d <- filtered()
    validate(need(length(unique(d$lesion)) == 2, "Need both groups for comparison."))

    test <- survdiff(Surv(age, dead) ~ lesion, data = d)
    p <- 1 - pchisq(test$chisq, df = length(test$n) - 1)

    cat(sprintf("Chi-squared: %.2f\n", test$chisq))
    cat(sprintf("p-value:     %.4f\n", p))
    cat(sprintf("\n%s at alpha = 0.05",
                ifelse(p < 0.05, "SIGNIFICANT", "Not significant")))
  })

  output$interpretation <- renderUI({
    req(sim())
    d <- filtered()
    if (length(unique(d$lesion)) < 2) return(NULL)

    test <- tryCatch(
      survdiff(Surv(age, dead) ~ lesion, data = d),
      error = function(e) NULL
    )
    if (is.null(test)) return(NULL)

    p <- 1 - pchisq(test$chisq, df = length(test$n) - 1)
    sig <- p < 0.05
    rmr <- input$rmr

    if (rmr == 1 && sig) {
      cls <- "interp interp-warn"
      msg <- paste(
        "FALSE POSITIVE: The test finds a significant difference,",
        "but lesions have no mortality effect in this simulation (RMR = 1).",
        "This is driven by the age structure of lesion formation &mdash;",
        "the osteological paradox in action.",
        "Try raising the minimum age above the formation window."
      )
    } else if (rmr == 1 && !sig) {
      cls <- "interp interp-ok"
      msg <- paste(
        "Correct inference: no significant difference detected,",
        "and lesions truly have no effect on mortality in this simulation."
      )
    } else if (rmr > 1 && sig) {
      cls <- "interp interp-ok"
      msg <- paste(
        "Significant difference detected &mdash; and lesions truly do increase",
        "mortality in this simulation (RMR =", rmr, ").",
        "Examine the survival curves to check whether the detected direction",
        "of effect matches the true direction."
      )
    } else {
      cls <- "interp interp-warn"
      msg <- paste(
        "FALSE NEGATIVE: No significant difference detected, despite lesions",
        "increasing mortality risk in this simulation (RMR =", rmr, ").",
        "The true effect may be masked by age-structured confounding.",
        "Try adjusting the minimum age filter."
      )
    }

    div(class = cls, HTML(msg))
  })

  # =========================================================================
  # Observed Sample Tab
  # =========================================================================

  # plot_alive_observed
  
  # plot_lesion_pct_observed

  # plot_cem_age_observed

  output$plot_cem_age_observed <- renderPlot({
    validate(need(sim(), "Click 'Run Simulation' to begin."))
    ggplot(cemetery_observed(), aes(x = Age_Interval, fill = Lesion_Status)) +
      geom_bar(position = "stack", width = 0.8) +
      scale_fill_manual(values = c("Absent" = col_no_lesion, "Present" = col_lesion)) +
      labs(title = "Age-at-Death Distribution",
           x = "Age at Death", y = "Count", fill = "lesion") +
      theme_bw(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      )
  })

  # plot_cem_lesion_observed

  output$plot_cem_lesion_observed <- renderPlot({
    validate(need(sim(), ""))
    prev <- cemetery_observed() %>%
      group_by(Age_Interval) %>%
      summarise(Pct = sum(lesion) / n() * 100, n = n(), .groups = "drop")

    ggplot(prev, aes(x = Age_Interval, y = Pct)) +
      geom_col(fill = col_lesion, alpha = 0.85, width = 0.8) +
      geom_text(aes(label = paste0("n=", n)), vjust = -0.5, size = 3.2, color = "#555") +
      labs(title = "Lesion Prevalence by Age Group",
           x = "Age at Death", y = "% with Lesions") +
      theme_bw(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      coord_cartesian(ylim = c(0, max(prev$Pct, na.rm = TRUE) * 1.15))
  })

  # cem_summary_observed

  output$cem_summary_observed <- renderPrint({
    validate(need(sim(), ""))
    d <- sim()$individual_outcomes
    d <- filter(d, in_sample)
    with_les <- sum(d$lesion == 1)
    cat(sprintf("Total individuals: %d\n", nrow(d)))
    cat(sprintf("With lesions: %d (%.1f%%)\n", with_les, with_les / nrow(d) * 100))
    cat(sprintf("Mean Age at death:  %.1f years (all)\n", mean(d$age)))
    if (with_les > 0) {
      cat(sprintf("                    %.1f years (with lesion)\n", mean(d$age[d$lesion == 1])))
    }
    cat(sprintf("                    %.1f years (without lesion)\n", mean(d$age[d$lesion == 0])))
  })


  # plot_km_observed

  output$plot_km_observed <- renderPlot({
    validate(need(sim(), "Click 'Run Simulation' to begin."))
    d <- filtered()
    d <- dplyr::filter(d, in_sample)
    validate(need(nrow(d) > 10, "Too few individuals remain after filtering."))
    validate(need(length(unique(d$lesion)) == 2,
                  "All remaining individuals have the same lesion status."))

    fit <- survfit(Surv(age, dead) ~ lesion, data = d)
    s <- summary(fit)

    surv_df <- data.frame(
      time = s$time,
      survival = s$surv,
      group = ifelse(grepl("=0", as.character(s$strata)), "No Lesion", "Has Lesion")
    )

    # Add starting points so step curves begin at survival = 1
    starts <- surv_df %>%
      group_by(group) %>%
      summarise(time = min(time) - 0.5, .groups = "drop") %>%
      mutate(survival = 1.0, time = pmax(time, 0))
    surv_df <- bind_rows(starts, surv_df) %>% arrange(group, time)

    title <- "Kaplan-Meier Survival Curves"
    if (input$min_age > 0) title <- paste0(title, "  (ages \u2265 ", input$min_age, ")")

    ggplot(surv_df, aes(x = time, y = survival, color = group)) +
      geom_step(linewidth = 1.1) +
      scale_color_manual(values = c("No Lesion" = col_dark, "Has Lesion" = col_lesion)) +
      labs(title = title, x = "Age (years)", y = "Survival Probability", color = "") +
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 12)
      )
  })



  # =========================================================================
  # Data Tab
  # =========================================================================

  output$table_cemetery <- renderTable({
    validate(need(sim(), "Run a simulation to see data."))
    sim()$individual_outcomes
  }, digits = 0)

  output$table_survivors <- renderTable({
    validate(need(sim(), ""))
    sim()$survivors
  }, digits = 1)

  output$dl_cemetery <- downloadHandler(
    filename = function() {
      paste0("persephone_cemetery_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(sim())
      write.csv(sim()$individual_outcomes, file, row.names = FALSE)
    }
  )

  output$dl_survivors <- downloadHandler(
    filename = function() {
      paste0("persephone_survivors_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(sim())
      write.csv(sim()$survivors, file, row.names = FALSE)
    }
  )

  # =========================================================================
  # Set-up Tab
  # =========================================================================

  output$plot_mort_siler_setup <- renderPlot({
  
    regime <- mortality_regimes[[input$mortality_regime]]

    age_max <- 100

    mort_siler_df <- data.frame(
      age = 0:age_max
    ) %>%
      mutate(
        Juvenile = regime$a1 * exp(-regime$b1 * age),
        Background = regime$a2,
        Senescent = regime$a3 * exp(regime$b3 * age),
        Total = Juvenile + Background + Senescent
      )

    mort_siler_df <- mort_siler_df %>%
      mutate(
        Lesion_Modified = case_when(
          input$risk_type == "proportional" ~
            Total * input$rmr,

          input$risk_type == "time_decreasing" ~
            Total * input$rmr / ((age / 10) + input$rmr),

          input$risk_type == "time_increasing" ~
            Total * ((age / 10) + input$rmr) / input$rmr
        )
      )

    ggplot(mort_siler_df, aes(x = age)) +
      geom_line(aes(y = Total, color = "Baseline hazard"), linewidth = 1.2) +
      geom_line(aes(y = Lesion_Modified, color = "Lesion-modified hazard"), linewidth = 1.2) +
      scale_color_manual(
        values = c(
          "Baseline hazard" = col_dark,
          "Lesion-modified hazard" = col_lesion
        )
      ) +
      labs(
        title = paste0("Mortality Hazards: ", input$mortality_regime),
        x = "Age (years)",
        y = "Annual mortality hazard",
        color = NULL
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      )
  })

  output$plot_death_pmf <- renderPlot({
    regime <- mortality_regimes[[input$mortality_regime]]

    d_base <- make_discrete_death_distribution(
      prob_fun = function(age) {
        baseline_death_prob(age, regime)
      },
      age_max = 100
    ) %>%
      mutate(Group = "Baseline")

    d_lesion <- make_discrete_death_distribution(
      prob_fun = function(age) {
        lesion_death_prob(
          age = age,
          regime = regime,
          risk_type = input$risk_type,
          rmr = input$rmr
        )
      },
      age_max = 100
    ) %>%
      mutate(Group = "Lesion-modified")

    d_plot <- bind_rows(d_base, d_lesion)

    ggplot(d_plot, aes(x = age, y = dx, fill = Group)) +
      geom_col(position = "identity", alpha = 0.45) +
      scale_fill_manual(values = c(
        "Baseline" = col_dark,
        "Lesion-modified" = col_lesion
      )) +
      labs(
        title = paste0("Discrete Age-at-Death Distribution: ", input$mortality_regime),
        x = "Age (years)",
        y = "Probability of death at age x",
        fill = NULL
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      )
  })


  output$plot_age_error <- renderPlot({

    # hard coded for now
    sd_per_year <- 0.1375
    sd_at_20 <- 2
    bias_start <- 60
    bias_at_90 <- 25

    cfact <- data.frame(
      age = 0:100
    )

    cfact$error_sd <- case_when(
      cfact$age <= 20 ~ cfact$age / 20 * sd_at_20,
      cfact$age > 20 ~ sd_at_20 + (cfact$age - 20) * sd_per_year
    )

    bias_slope <- -(bias_at_90 / (90 - bias_start))

    cfact$error_bias <- case_when(
      cfact$age <= bias_start ~ 0,
      cfact$age > bias_start ~ (cfact$age - bias_start) * bias_slope
    )

    cfact$q25 <- qnorm(0.25, cfact$error_bias, cfact$error_sd)
    cfact$q75 <- qnorm(0.75, cfact$error_bias, cfact$error_sd)
    cfact$q05 <- qnorm(0.05, cfact$error_bias, cfact$error_sd)
    cfact$q95 <- qnorm(0.95, cfact$error_bias, cfact$error_sd)

    plot(NULL, ylim = c(min(cfact$q05), max(cfact$q95)), xlim = c(0, 100),
    xlab = "true age at death", ylab = "error in age estimation")
    polygon(c(cfact$age, rev(cfact$age)), c(cfact$q05, rev(cfact$q95)), border = NA, col = "orange")
    polygon(c(cfact$age, rev(cfact$age)), c(cfact$q25, rev(cfact$q75)), border = NA, col = "red")
    abline(h = 0, lty = 2)
    points(cfact$age, cfact$error_bias, type = "l")
    abline(h = seq(-50, 50, by = 10), col = gray(0.3, 0.3))
    abline(v = seq(0, 100, by = 10), col = gray(0.3, 0.3))

  })

  output$plot_taph_siler_setup <- renderPlot({

    age_max <- 100

    taph_siler_df <- bind_rows(
      lapply(names(taphonomy_regimes), function(reg_name) {
        regime <- taphonomy_regimes[[reg_name]]

        data.frame(age = 0:age_max) %>%
          mutate(
            Juvenile = regime$a1 * exp(-regime$b1 * age),
            Background = regime$a2,
            Senescent = regime$a3 * exp(regime$b3 * age),
            Total = Juvenile + Background + Senescent,
            Regime = reg_name
          )
      })
    )

    ggplot(taph_siler_df, aes(x = age, y = Total, color = Regime)) +
      geom_line(linewidth = 1.2) +
      labs(
        title = "Taphonomic Siler Hazards by Regime",
        x = "Age (years)",
        y = "Annual hazard of loss",
        color = NULL
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      )
  })

  output$plot_true_vs_observed_age <- renderPlot({

    validate(need(sim(), "Click 'Run Simulation' to begin."))

    d <- sim()$individual_outcomes

    validate(need("age" %in% names(d), "Column `age` missing"))
    validate(need("estimated_age" %in% names(d), "Column `estimated_age` missing"))
    validate(need("in_sample" %in% names(d), "Column `in_sample` missing"))

    true_df <- d %>%
      count(age) %>%
      rename(plot_age = age, count = n) %>%
      mutate(Source = "True age distribution")

    observed_df <- d %>%
      filter(in_sample) %>%
      count(estimated_age) %>%
      rename(plot_age = estimated_age, count = n) %>%
      mutate(Source = "Observed sample (estimated age)")

    plot_df <- bind_rows(true_df, observed_df)

    ggplot(plot_df, aes(x = plot_age, y = count, fill = Source)) +
      geom_col(position = "identity", alpha = 0.45, width = 0.9) +
      scale_fill_manual(values = c(
        "True age distribution" = col_dark,
        "Observed sample (estimated age)" = col_lesion
      )) +
      labs(
        title = "True vs Observed Age-at-Death Distributions",
        x = "Age at death",
        y = "Count",
        fill = NULL
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      )

  })

} # end server() function

# =============================================================================
# Launch
# =============================================================================

shinyApp(ui = ui, server = server)
