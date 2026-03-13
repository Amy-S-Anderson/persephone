# =============================================================================
# Persephone Shiny App
# Agent-Based Model for Bioarchaeology
# Based on Anderson & DeWitte, "Known Unknowns and the Osteological Paradox"
# =============================================================================

library(shiny)
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
      sliderInput(
        "remove_adults_older",
        "Remove Adults older than age...",
        min = 50,
        max = 100,
        value = 100,
        step = 10
      ),

      h4("Taphonomic Filter", class = "param-header"),
      selectInput(
        "taphonomic_filter",
        "Taphonomic strength",
        choices = c("Weak", "Medium", "Strong"),
        selected = "Medium"
      ),

      h4("Age Estimation Bias", class = "param-header"),
      checkboxInput(
        "age_bias_exists",
        "Age Bias Exists",
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
            column(12, plotOutput("plot_siler_setup", height = "500px"))
          ),
          fluidRow(
            column(12, plotOutput("plot_death_pmf", height = "400px"))
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

        # =====================================================================
        # Tab: Data
        # =====================================================================
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
        mortality_regime = mortality_regimes[[input$mortality_regime]]
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
        Age_Interval = assign_age_interval(Age),
        Lesion_Status = factor(
          ifelse(Lesion == 1, "Present", "Absent"),
          levels = c("Absent", "Present")
        )
      )
  })

  # =========================================================================
  # Cemetery Tab
  # =========================================================================

  output$plot_cem_age <- renderPlot({
    validate(need(sim(), "Click 'Run Simulation' to begin."))
    ggplot(cemetery(), aes(x = Age_Interval, fill = Lesion_Status)) +
      geom_bar(position = "stack", width = 0.8) +
      scale_fill_manual(values = c("Absent" = col_no_lesion, "Present" = col_lesion)) +
      labs(title = "Age-at-Death Distribution",
           x = "Age at Death", y = "Count", fill = "Lesion") +
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
      summarise(Pct = sum(Lesion) / n() * 100, n = n(), .groups = "drop")

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
    with_les <- sum(d$Lesion == 1)
    cat(sprintf("Total individuals: %d\n", nrow(d)))
    cat(sprintf("With lesions: %d (%.1f%%)\n", with_les, with_les / nrow(d) * 100))
    cat(sprintf("Mean age at death:  %.1f years (all)\n", mean(d$Age)))
    if (with_les > 0) {
      cat(sprintf("                    %.1f years (with lesion)\n", mean(d$Age[d$Lesion == 1])))
    }
    cat(sprintf("                    %.1f years (without lesion)\n", mean(d$Age[d$Lesion == 0])))
  })

  # =========================================================================
  # Living Population Tab
  # =========================================================================

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

  # =========================================================================
  # Survival Analysis Tab
  # =========================================================================

  filtered <- reactive({
    req(sim())
    d <- sim()$individual_outcomes %>% filter(Age >= input$min_age)
    d$Dead <- 1
    d
  })

  output$plot_km <- renderPlot({
    validate(need(sim(), "Click 'Run Simulation' to begin."))
    d <- filtered()
    validate(need(nrow(d) > 10, "Too few individuals remain after filtering."))
    validate(need(length(unique(d$Lesion)) == 2,
                  "All remaining individuals have the same lesion status."))

    fit <- survfit(Surv(Age, Dead) ~ Lesion, data = d)
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
    validate(need(length(unique(d$Lesion)) == 2, "Need both groups for comparison."))

    test <- survdiff(Surv(Age, Dead) ~ Lesion, data = d)
    p <- 1 - pchisq(test$chisq, df = length(test$n) - 1)

    cat(sprintf("Chi-squared: %.2f\n", test$chisq))
    cat(sprintf("p-value:     %.4f\n", p))
    cat(sprintf("\n%s at alpha = 0.05",
                ifelse(p < 0.05, "SIGNIFICANT", "Not significant")))
  })

  output$interpretation <- renderUI({
    req(sim())
    d <- filtered()
    if (length(unique(d$Lesion)) < 2) return(NULL)

    test <- tryCatch(
      survdiff(Surv(Age, Dead) ~ Lesion, data = d),
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

  output$plot_siler_setup <- renderPlot({
  regime <- mortality_regimes[[input$mortality_regime]]

  age_max <- 100

  siler_df <- data.frame(
    Age = 0:age_max
  ) %>%
    mutate(
      Juvenile = regime$a1 * exp(-regime$b1 * Age),
      Background = regime$a2,
      Senescent = regime$a3 * exp(regime$b3 * Age),
      Total = Juvenile + Background + Senescent
    )

  siler_df <- siler_df %>%
    mutate(
      Lesion_Modified = case_when(
        input$risk_type == "proportional" ~
          Total * input$rmr,

        input$risk_type == "time_decreasing" ~
          Total * input$rmr / ((Age / 10) + input$rmr),

        input$risk_type == "time_increasing" ~
          Total * ((Age / 10) + input$rmr) / input$rmr
      )
    )

  ggplot(siler_df, aes(x = Age)) +
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

  ggplot(d_plot, aes(x = Age, y = dx, fill = Group)) +
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

}

# =============================================================================
# Launch
# =============================================================================

shinyApp(ui = ui, server = server)
