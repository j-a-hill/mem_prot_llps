library(shiny)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)
source("R/data.R")

PLLPS_COLORS <- c("High" = "#e74c3c", "Medium" = "#f39c12", "Low" = "#3498db")

# ── UI ───────────────────────────────────────────────────────────────────────

ui <- fluidPage(
  tags$head(
    tags$title("LLPS Protein Data Explorer"),
    tags$style(HTML("
      body { font-family: 'Segoe UI', sans-serif; }
      .metric-box { background: #fff; border: 1px solid #dee2e6; border-radius: 6px;
                    padding: 14px 18px; text-align: center; margin-bottom: 12px; }
      .metric-label { font-size: 12px; color: #6c757d; font-weight: 600; text-transform: uppercase; }
      .metric-val   { font-size: 28px; font-weight: 700; color: #212529; }
    "))
  ),
  div(style = "max-width: 1400px; margin: auto;",
    h2("LLPS Protein Data Explorer", style = "margin-top: 1rem;"),
    p("Explore protein phase-separation (pLLPS) data. Use the full dataset or upload your own XLSX/CSV file."),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h5("Data Source"),
        radioButtons("data_source", NULL,
          choices = c("Full dataset" = "default", "Upload file" = "upload"),
          selected = "default"
        ),
        conditionalPanel(
          "input.data_source === 'upload'",
          fileInput("file_upload", "XLSX or CSV", accept = c(".xlsx", ".csv"))
        ),
        uiOutput("data_status"),
        hr(),
        h5("Filters"),
        p(style = "font-size: 12px; color: #6c757d; margin-bottom: 8px;",
          "All plots and the table update together when filters change."),
        uiOutput("filter_controls")
      ),
      mainPanel(
        width = 9,
        uiOutput("metrics_row"),
        hr(),
        h4("Protein Table"),
        textInput("search_text", "Search by name / entry ID",
          placeholder = "e.g. PCLO, Q9Y6V0, kinase..."),
        DTOutput("data_table"),
        hr(),
        h4("Score & Feature Distributions"),
        fluidRow(
          column(4, plotlyOutput("plot_pllps_dist", height = "280px")),
          column(4, plotlyOutput("plot_length_dist", height = "280px")),
          column(4, plotlyOutput("plot_tmd_dist", height = "280px"))
        ),
        hr(),
        h4("p(LLPS) by Category"),
        p(style = "color: #6c757d; font-size: 13px;",
          "Box plots (median, IQR, 1.5x IQR whiskers) sorted by median p(LLPS) descending."),
        fluidRow(
          column(6, plotlyOutput("plot_pllps_by_location", height = "320px")),
          column(6, plotlyOutput("plot_pllps_by_function", height = "320px"))
        ),
        plotlyOutput("plot_pllps_by_tmd", height = "300px"),
        hr(),
        h4("Scatter"),
        fluidRow(
          column(3, selectInput("scatter_x", "X axis",
            choices = c("p(LLPS)", "Length", "TMD_count", "n(DPR=> 25)"))),
          column(3, selectInput("scatter_y", "Y axis",
            choices = c("Length", "p(LLPS)", "TMD_count", "n(DPR=> 25)")))
        ),
        plotlyOutput("plot_scatter", height = "350px"),
        hr(),
        h4("Category Overview"),
        fluidRow(
          column(6, plotlyOutput("plot_locations", height = "300px")),
          column(6, plotlyOutput("plot_functions", height = "300px"))
        ),
        hr(),
        h4("Download Filtered Data"),
        uiOutput("export_info"),
        downloadButton("download_csv", "Download as CSV", class = "btn btn-primary")
      )
    )
  )
)

# ── Server ───────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  raw_df <- reactive({
    if (input$data_source == "default") {
      load_default()
    } else {
      req(input$file_upload)
      load_data(input$file_upload$datapath)
    }
  })

  filtered_df <- reactive({
    df <- raw_df()

    if ("p(LLPS)" %in% names(df) && !is.null(input$filter_pllps)) {
      lo <- input$filter_pllps[1]; hi <- input$filter_pllps[2]
      df <- df[!is.na(df[["p(LLPS)"]]) & df[["p(LLPS)"]] >= lo & df[["p(LLPS)"]] <= hi, ]
    }

    if ("Length" %in% names(df) && !is.null(input$filter_length)) {
      lo <- input$filter_length[1]; hi <- input$filter_length[2]
      df <- df[!is.na(df$Length) & df$Length >= lo & df$Length <= hi, ]
    }

    if ("TMD_count" %in% names(df) && !is.null(input$filter_tmd)) {
      lo <- input$filter_tmd[1]; hi <- input$filter_tmd[2]
      df <- df[!is.na(df$TMD_count) & df$TMD_count >= lo & df$TMD_count <= hi, ]
    }

    if ("pLLPS_class" %in% names(df) && !is.null(input$filter_class) && length(input$filter_class) > 0) {
      df <- df[as.character(df$pLLPS_class) %in% input$filter_class, ]
    }

    if ("location_categories" %in% names(df) && length(input$filter_locations) > 0) {
      sel <- input$filter_locations
      df <- df[sapply(df$location_categories, function(x) any(x %in% sel)), ]
    }

    if ("function_categories" %in% names(df) && length(input$filter_functions) > 0) {
      sel <- input$filter_functions
      df <- df[sapply(df$function_categories, function(x) any(x %in% sel)), ]
    }

    q <- trimws(input$search_text)
    if (nzchar(q)) {
      q <- tolower(q)
      search_cols <- intersect(c("Entry", "Entry name", "Protein names"), names(df))
      mask <- Reduce(`|`, lapply(search_cols, function(col) {
        grepl(q, tolower(as.character(df[[col]])), fixed = TRUE)
      }))
      df <- df[mask, ]
    }

    df
  })

  # ── Sidebar ──────────────────────────────────────────────────────────────

  output$data_status <- renderUI({
    df <- raw_df()
    label <- if (input$data_source == "default") "Full dataset" else "Uploaded data"
    div(class = "alert alert-success mt-2",
        sprintf("%s: %d proteins", label, nrow(df)))
  })

  output$filter_controls <- renderUI({
    df <- raw_df()
    controls <- list()

    if ("p(LLPS)" %in% names(df)) {
      lo <- floor(min(df[["p(LLPS)"]], na.rm = TRUE) * 100) / 100
      hi <- ceiling(max(df[["p(LLPS)"]], na.rm = TRUE) * 100) / 100
      controls <- c(controls, list(
        sliderInput("filter_pllps", "p(LLPS) range", min = lo, max = hi,
                    value = c(lo, hi), step = 0.01)
      ))
    }

    if ("Length" %in% names(df)) {
      lo <- min(df$Length, na.rm = TRUE); hi <- max(df$Length, na.rm = TRUE)
      controls <- c(controls, list(
        sliderInput("filter_length", "Protein length", min = lo, max = hi,
                    value = c(lo, hi), step = 10)
      ))
    }

    if ("TMD_count" %in% names(df)) {
      lo <- min(df$TMD_count, na.rm = TRUE); hi <- max(df$TMD_count, na.rm = TRUE)
      controls <- c(controls, list(
        sliderInput("filter_tmd", "Transmembrane domains", min = lo, max = hi,
                    value = c(lo, hi), step = 1)
      ))
    }

    if ("pLLPS_class" %in% names(df)) {
      controls <- c(controls, list(
        checkboxGroupInput("filter_class", "pLLPS class",
          choices = c("High", "Medium", "Low"),
          selected = c("High", "Medium", "Low"))
      ))
    }

    if ("location_categories" %in% names(df)) {
      all_locs <- sort(unique(unlist(df$location_categories)))
      all_locs <- all_locs[nzchar(all_locs)]
      if (length(all_locs) > 0) {
        controls <- c(controls, list(
          selectizeInput("filter_locations", "Subcellular location",
            choices = all_locs, multiple = TRUE,
            options = list(placeholder = "All locations"))
        ))
      }
    }

    if ("function_categories" %in% names(df)) {
      all_funcs <- sort(unique(unlist(df$function_categories)))
      all_funcs <- all_funcs[nzchar(all_funcs)]
      if (length(all_funcs) > 0) {
        controls <- c(controls, list(
          selectizeInput("filter_functions", "Functional category",
            choices = all_funcs, multiple = TRUE,
            options = list(placeholder = "All functions"))
        ))
      }
    }

    do.call(tagList, controls)
  })

  # ── Metrics ───────────────────────────────────────────────────────────────

  output$metrics_row <- renderUI({
    df <- filtered_df()
    total <- nrow(df)
    high  <- if ("pLLPS_class" %in% names(df)) sum(as.character(df$pLLPS_class) == "High") else "N/A"
    avg_p <- if ("p(LLPS)" %in% names(df)) sprintf("%.3f", mean(df[["p(LLPS)"]], na.rm = TRUE)) else "N/A"
    avg_l <- if ("Length"   %in% names(df)) sprintf("%.0f",  mean(df$Length,      na.rm = TRUE)) else "N/A"

    box <- function(label, value) {
      column(3, div(class = "metric-box",
        div(class = "metric-label", label),
        div(class = "metric-val", as.character(value))
      ))
    }

    fluidRow(
      box("Proteins shown", total),
      box("High pLLPS", high),
      box("Mean p(LLPS)", avg_p),
      box("Mean length (aa)", avg_l)
    )
  })

  # ── Table ─────────────────────────────────────────────────────────────────

  output$data_table <- renderDT({
    df <- filtered_df()
    display_cols <- intersect(
      c("Entry", "Entry name", "Protein names", "p(LLPS)", "pLLPS_class",
        "Length", "TMD_count", "Organism"),
      names(df)
    )
    datatable(df[display_cols], rownames = FALSE,
              options = list(pageLength = 15, scrollX = TRUE))
  })

  # ── Distribution plots ────────────────────────────────────────────────────

  output$plot_pllps_dist <- renderPlotly({
    df <- filtered_df()
    req("p(LLPS)" %in% names(df))
    plot_ly(df, x = df[["p(LLPS)"]], color = ~pLLPS_class, type = "histogram",
            colors = PLLPS_COLORS, nbinsx = 30, alpha = 0.9) %>%
      layout(barmode = "stack",
             title = list(text = "Distribution of p(LLPS) scores", font = list(size = 13)),
             xaxis = list(title = "p(LLPS) score"),
             yaxis = list(title = "Count"),
             legend = list(title = list(text = "pLLPS class")))
  })

  output$plot_length_dist <- renderPlotly({
    df <- filtered_df()
    req("Length" %in% names(df))
    plot_ly(df, x = ~Length, color = ~pLLPS_class, type = "histogram",
            colors = PLLPS_COLORS, nbinsx = 30, alpha = 0.9) %>%
      layout(barmode = "stack",
             title = list(text = "Distribution of protein lengths", font = list(size = 13)),
             xaxis = list(title = "Protein length (aa)"),
             yaxis = list(title = "Count"))
  })

  output$plot_tmd_dist <- renderPlotly({
    df <- filtered_df()
    req("TMD_count" %in% names(df))
    tmd_data <- df %>%
      group_by(TMD_count, pLLPS_class) %>%
      summarise(Count = n(), .groups = "drop")
    plot_ly(tmd_data, x = ~TMD_count, y = ~Count, color = ~pLLPS_class,
            type = "bar", colors = PLLPS_COLORS) %>%
      layout(barmode = "stack",
             title = list(text = "Transmembrane domain count", font = list(size = 13)),
             xaxis = list(title = "Number of TM domains"),
             yaxis = list(title = "Protein count"))
  })

  # ── Box plots by category ─────────────────────────────────────────────────

  output$plot_pllps_by_location <- renderPlotly({
    df <- filtered_df()
    req("location_categories" %in% names(df), "p(LLPS)" %in% names(df))

    loc_df <- df %>%
      select(location_categories, pllps = `p(LLPS)`) %>%
      unnest(cols = location_categories) %>%
      filter(!is.na(location_categories), nzchar(location_categories))

    req(nrow(loc_df) >= 3)

    top15 <- loc_df %>% count(location_categories) %>%
      slice_max(n, n = 15) %>% pull(location_categories)

    loc_df <- loc_df %>% filter(location_categories %in% top15)

    order_locs <- loc_df %>%
      group_by(location_categories) %>%
      summarise(med = median(pllps, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(med)) %>%
      pull(location_categories)

    plot_ly(loc_df, x = ~location_categories, y = ~pllps, type = "box",
            color = ~location_categories, showlegend = FALSE) %>%
      layout(
        title = list(text = "p(LLPS) by subcellular location (top 15)", font = list(size = 13)),
        xaxis = list(title = "", categoryorder = "array",
                     categoryarray = order_locs, tickangle = -40),
        yaxis = list(title = "p(LLPS)", range = c(0, 1))
      )
  })

  output$plot_pllps_by_function <- renderPlotly({
    df <- filtered_df()
    req("function_categories" %in% names(df), "p(LLPS)" %in% names(df))

    func_df <- df %>%
      select(function_categories, pllps = `p(LLPS)`) %>%
      unnest(cols = function_categories) %>%
      filter(!is.na(function_categories), nzchar(function_categories))

    req(nrow(func_df) >= 3)

    order_funcs <- func_df %>%
      group_by(function_categories) %>%
      summarise(med = median(pllps, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(med)) %>%
      pull(function_categories)

    plot_ly(func_df, x = ~function_categories, y = ~pllps, type = "box",
            color = ~function_categories, showlegend = FALSE) %>%
      layout(
        title = list(text = "p(LLPS) by functional category", font = list(size = 13)),
        xaxis = list(title = "", categoryorder = "array",
                     categoryarray = order_funcs, tickangle = -40),
        yaxis = list(title = "p(LLPS)", range = c(0, 1))
      )
  })

  output$plot_pllps_by_tmd <- renderPlotly({
    df <- filtered_df()
    req("TMD_count" %in% names(df), "p(LLPS)" %in% names(df))

    plot_df <- df[df$TMD_count <= 20, ]
    req(nrow(plot_df) >= 3)

    plot_ly(plot_df, x = ~as.factor(TMD_count), y = plot_df[["p(LLPS)"]],
            type = "box", color = ~as.factor(TMD_count),
            showlegend = FALSE, colors = "Viridis") %>%
      layout(
        title = list(text = sprintf("p(LLPS) by TM domain count (n=%d)", nrow(plot_df)),
                     font = list(size = 13)),
        xaxis = list(title = "Number of TM domains"),
        yaxis = list(title = "p(LLPS)", range = c(0, 1))
      )
  })

  # ── Scatter ───────────────────────────────────────────────────────────────

  output$plot_scatter <- renderPlotly({
    df <- filtered_df()
    x_col <- input$scatter_x; y_col <- input$scatter_y
    req(x_col %in% names(df), y_col %in% names(df))

    x_data <- df[[x_col]]; y_data <- df[[y_col]]

    hover_text <- paste0(
      if ("Entry"         %in% names(df)) paste0("Entry: ", df$Entry, "<br>") else "",
      if ("Protein names" %in% names(df)) paste0(df[["Protein names"]], "<br>") else "",
      x_col, ": ", round(x_data, 3), "<br>",
      y_col, ": ", round(y_data, 3)
    )

    plot_ly(
      x = x_data, y = y_data, type = "scatter", mode = "markers",
      color = if ("pLLPS_class" %in% names(df)) df$pLLPS_class else NULL,
      colors = PLLPS_COLORS,
      text = hover_text, hoverinfo = "text",
      marker = list(opacity = 0.5, size = 5)
    ) %>%
      layout(
        title = paste(y_col, "vs", x_col),
        xaxis = list(title = x_col),
        yaxis = list(title = y_col)
      )
  })

  # ── Category count bars ────────────────────────────────────────────────────

  output$plot_locations <- renderPlotly({
    df <- filtered_df()
    req("location_categories" %in% names(df))

    counts <- df %>%
      unnest(cols = location_categories) %>%
      filter(!is.na(location_categories), nzchar(location_categories)) %>%
      count(location_categories, name = "Count") %>%
      arrange(desc(Count))

    req(nrow(counts) > 0)

    plot_ly(counts, x = ~reorder(location_categories, -Count), y = ~Count,
            type = "bar", color = ~Count, colors = "Blues", showlegend = FALSE) %>%
      layout(title = list(text = "Proteins per subcellular location", font = list(size = 13)),
             xaxis = list(title = "", tickangle = -40),
             yaxis = list(title = "Count"))
  })

  output$plot_functions <- renderPlotly({
    df <- filtered_df()
    req("function_categories" %in% names(df))

    counts <- df %>%
      unnest(cols = function_categories) %>%
      filter(!is.na(function_categories), nzchar(function_categories)) %>%
      count(function_categories, name = "Count") %>%
      arrange(desc(Count))

    req(nrow(counts) > 0)

    plot_ly(counts, x = ~reorder(function_categories, -Count), y = ~Count,
            type = "bar", color = ~Count, colors = "Greens", showlegend = FALSE) %>%
      layout(title = list(text = "Proteins per functional category", font = list(size = 13)),
             xaxis = list(title = "", tickangle = -40),
             yaxis = list(title = "Count"))
  })

  # ── Export ────────────────────────────────────────────────────────────────

  output$export_info <- renderUI({
    df <- filtered_df()
    p(sprintf("%d proteins · %d columns (after filters)", nrow(df), ncol(df)))
  })

  output$download_csv <- downloadHandler(
    filename = "llps_proteins_filtered.csv",
    content = function(file) {
      df <- filtered_df()
      export_cols <- setdiff(names(df), c("location_categories", "function_categories"))
      write.csv(df[export_cols], file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
