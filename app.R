#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(quicR)
library(tidyverse)
library(readxl)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Getcha QuIC Data Heah"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fileInput(
              "file",
              "Choose File",
              accept = ".xlsx"
            ),
            textInput(
              "threshold",
              "Choose a threshold for calculating RAF."
            ),
            downloadButton(
              "downlad_excel",
              "Download Data to Excel"
            )
        ),

        # Show a plot of the generated distribution
        mainPanel(
           # tableOutput("df_")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  print(1)
  output$df_ <- renderTable({
    file <- input$file
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    
    # Get the real-time data.
    quicR::get_real(file$datapath, ordered = FALSE)
    
    # Ask the user which real-time data set they want to use.
    df_id <- 1
    
    # Ask the user for the threshold.
    threshold <- as.numeric(input$threshold)
    
    # Select the real-time data set that the user signified.
    df <- df[[df_id]]
    
    
    
    # Get the metadata of each sample -----------------------------------------
    
    
    print(2)
    # Export the tables in the first sheet of the file.
    dic <- quicR::organize_tables(file$datapath)
    
    column_names <- c("Time")
    for (i in t(dic[["Sample IDs"]])) {
      for (j in i) {
        if (!is.na(j)) {
          column_names <- cbind(column_names, j)
        }
      }
    }
    
    # Apply the column names.
    colnames(df) <- column_names
    
    # Determine if there is a dilutions table.
    dilution_bool <- "Dilutions" %in% names(dic)
    
    
    
    # Calculate the normalized real-time data ---------------------------------
    
    
    print(3)
    # Calculate the normalized real-time data.
    df_norm <- quicR::normalize_RFU(df)
    
    # Add dilution factors if applicable.
    if (dilution_bool) {
      dilutions <- c()
      for (i in t(dic[["Dilutions"]])) {
        for (j in i) {
          if (!is.na(j)) {
            dilutions <- rbind(dilutions, j)
          }
        }
      }
    }
    
    
    
    # Calculate the relevant metrics ------------------------------------------
    
    
    
    # Define the number of hours that the rxn ran for.
    hours <- as.numeric(colnames(df_norm)[ncol(df_norm)])
    
    # Initialized the dataframe with the calculated metrics.
    df_analyzed <- data.frame(`Sample_ID` = df_norm$`Sample ID`) %>%
      
      mutate(
        # Add dilutions if applicable.
        Dilutions = if (dilution_bool) -log10(as.numeric(dilutions)),
        # Maxpoint Ratio
        MPR = quicR::calculate_MPR(df_norm, start_col=3, data_is_norm=TRUE),
        # Max Slope
        MS  = quicR::calculate_MS(df_norm, start_col=3),
        # Time to Threshold
        TtT = quicR::calculate_TtT(df_norm, threshold=threshold, start_col=3, run_time=hours)
      ) %>%
      
      mutate(
        # Rate of Amyloid Formation
        RAF = ifelse(TtT == run_time, 0, 1 / (3600 * TtT)),
        # Crossed threshold?
        crossed = TtT != run_time
      ) %>%
      
      # Order the data frame based on Sample_ID.
      arrange(Sample_ID)
    
    
    
    # Summarize the data ------------------------------------------------------
    
    
    
    # Create a summary data frame.
    summary <- (
      if (dilution_bool) {
        summary <- df_analyzed %>%
          group_by(Sample_ID, Dilutions)
      } else {
        summary <- df_analyzed %>%
          group_by(Sample_ID)
      }
    ) %>%
      
      summarise(
        reps      = n(),
        mean_TtT  = mean(TtT),
        sd_TtT    = sd(TtT),
        mean_RAF  = mean(RAF),
        sd_RAF    = sd(RAF),
        mean_MPR  = mean(MPR),
        sd_MPR    = sd(MPR),
        mean_MS   = mean(MS),
        sd_MS     = sd(MS),
        thres_pos = sum(crossed) / n() > 0.5
      )
    
    
    
    # Run the statistical analysis against the negative control ---------------
    
    
    
    metrics <- c("MPR", "MS")
    for (metric in metrics) {
      
      # Create a dataframe of the individual comparisons.
      comps <- LSD.test( # Perform the post-hoc multiple comparisons test.
        # Create the statistical model using ANOVA.
        aov(as.formula(paste0(metric, " ~ ", "Sample_ID")),
            data = df_analyzed),
        "Sample_ID",  p.adj = "holm", group = F
      )[["comparison"]]
      
      # Initialize columns which will hold unique IDs for each sample compared.
      comps <- comps %>%
        
        cbind(
          rownames(comps) %>%
            strsplit(" - ") %>%
            as.data.frame() %>%
            t() %>%
            as.data.frame()
        ) %>%
        
        select(-difference) %>%
        
        # Remove all comparisons that are not against "N".
        subset(V1 == "N" | V2 == "N") %>%
        
        rename(
          "{metric}_pvalue" := pvalue,
          "{metric}_significance" := signif.
        ) %>%
        
        mutate(
          V1 = replace(V1, V1 == "N", NA),
          V2 = replace(V2, V2 == "N", NA)
        ) %>%
        
        unite(
          Sample_ID,
          c("V1", "V2"),
          sep = "",
          na.rm = T
        ) %>%
        
        rbind(c(NA, NA, "N"))
      
      summary <- left_join(summary, comps)
    }
    
    summary <- summary %>%
      mutate(Positive = thres_pos & MPR_pvalue <= 0.05 & MS_pvalue <= 0.05)
    
    
    
    # Save the data to an Excel workbook --------------------------------------
    
    
    
    # Initialize the workbook for Excel.
    wb <- createWorkbook()
    
    # Add the sheets.
    addWorksheet(wb, "Total")
    addWorksheet(wb, "Summary")
    
    # Write the "summary" df to the "Summary" sheet.
    writeData(wb, "Total", df_analyzed)
    writeData(wb, "Summary", summary)
    
    # Save the Excel file.
    # saveWorkbook(wb, paste0(file$datapath, "/summary.xlsx"), overwrite = TRUE)
    
  })
  
}

# Run the application 
shinyApp(ui, server)
