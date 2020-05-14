library(shiny) #install.packages('shiny')
library(shinydashboard) #install.packages('shinydashboard')

source("global.R")

## ui

ui <- dashboardPage(
  dashboardHeader(title = "Fly Blood Explorer"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Clusters", tabName = "cluster"),
      menuItem("Markers", tabName = "markers")
    )
  ),
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "cluster",
              fluidPage(
                
                # Application title
                titlePanel("Drosophila blood single-cell RNA-seq explorer: Cluster"),
                
                # Sidebar with text inputs for number of bins 
                sidebarLayout(
                  sidebarPanel(
                    selectInput(inputId = "stage", label ="Stage", choices = c("Adult_circulating"), selected = "Adult_circulating"),
                    numericInput(inputId="numInputs", label="Number of target genes", 4),
                    textInput(inputId = "gene1", label = "Target gene1", value="Hml"),
                    textInput(inputId = "gene2", label = "Target gene2", value="Pxn"),
                    textInput(inputId = "gene3", label = "Target gene3", value="nSyb"),
                    textInput(inputId = "gene4", label = "Target gene4", value="SPARC"),
                    selectInput(inputId = "reduction", label ="Reduction method", choices = c("umap","tsne"), selected = "UMAP"),
                    downloadButton("downloadData1", "Download raw data"),
                    width=3
                  ),
                  
                  # Show a plot of the scRNAseq clustering landscape
                  mainPanel(
                    tabsetPanel(type = "tabs",
                                tabPanel("Clustering plot", plotOutput(outputId = "clustering",width=800, height=700)),
                                tabPanel("Expression level_violin plot", plotOutput(outputId = "expression",width=800, height=700)),
                                tabPanel("Expression level_raw data", tableOutput("table1")),
                                tabPanel("Expression level_summary", tableOutput("table2"))),
                    width=9
                  )
                )
              )
      ),
      
      # Second tab content
      tabItem(tabName = "markers",
              fluidPage(
                
                # Application title
                titlePanel("Drosophila blood single-cell RNA-seq explorer: Cluster markers"),
                
                # Sidebar with text inputs for number of bins 
                sidebarLayout(
                  sidebarPanel(
                    selectInput(inputId = "stage", label ="Stage", choices = c("Adult_circulating"), selected = "Adult_circulating"),
                    selectInput(inputId = "cluster", label = "Target cluster", choices = unique(adult_markers$cluster)),
                    sliderInput(inputId = "ngenes", label = "Number of marker genes", min=1, max=200, value=10),
                    downloadButton("downloadData2", "Download csv"),
                    width=3
                  ),
                  
                  # Show a plot of the scRNAseq clustering landscape
                  mainPanel(
                    tabsetPanel(type = "tabs",
                                tabPanel("Cluster markers", tableOutput("table"))),
                    width=9
                  )
                )
              )
      )
    )
  )
)



# Define server logic required to draw a histogram
server <- function(input, output) {
   
  datasetInput1 <- reactive(adult_circulating_exp %>%
                             dplyr::filter(gene %in% c(input$gene1, input$gene2, input$gene3, input$gene4)) %>%
                             tidyr::spread(key="gene", value="expression") %>%
                             dplyr::arrange(cluster))
  
   output$clustering <- renderPlot({
   
     FeaturePlot(adult_circulating, features = c(input$gene1, input$gene2, input$gene3, input$gene4), 
                 reduction=input$reduction, ncol=2)
     
   })
   
   output$expression <- renderPlot({
      
      VlnPlot(adult_circulating, features = c(input$gene1, input$gene2, input$gene3, input$gene4), ncol=2)
      
   })
   
   output$table1 <- renderTable({
      
     datasetInput1()
      
   })
   
   output$table2 <- renderTable({
      
      adult_circulating_exp %>%
         dplyr::filter(gene %in% c(input$gene1, input$gene2, input$gene3, input$gene4)) %>%
         dplyr::arrange(cluster) %>%
         dplyr::group_by(cluster, gene) %>%
         dplyr::summarise(expression_mean = mean(expression), expression_sd = sd(expression)) %>% 
         dplyr::ungroup() %>%
          dplyr::arrange(gene, cluster)
      
      
   })
   
   
   output$downloadData1 <- downloadHandler(
     filename = function() {
       paste(input$stage, "_", input$gene1, "_", input$gene2, "_", input$gene3, "_", input$gene4, ".csv", sep = "")
     },
     content = function(file) {
       write.csv(datasetInput1(), file, row.names = FALSE)
     }
   )
   
   datasetInput2 <- reactive(dplyr::arrange(dplyr::filter(adult_markers, cluster == input$cluster),
                                           P_adjust, abs(Fraction_target-Fraction_others), Mean_log_FC))
   
   output$table <- renderTable(digits=7, {
     
     datasetInput2()[1:input$ngenes,]
     
   })
   
   output$downloadData2 <- downloadHandler(
     filename = function() {
       paste(input$stage, "_", input$cluster, ".csv", sep = "")
     },
     content = function(file) {
       write.csv(datasetInput2(), file, row.names = FALSE)
     }
   )
   
}

# Run the application 
shinyApp(ui = ui, server = server)

