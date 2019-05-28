library(shiny)

ui <- fluidPage(
  selectInput("pdb", label="choose PDB", choices=unique(Table_DDG_residues$PDB)),
  titlePanel("Binding energy analysis"),
  fluidRow(
    column(2,
      tableOutput("DDG_binding"),
      tableOutput(outputId = "PTM_count")
    ),
    column(2,
           tableOutput(outputId="ChainRange"),
           tableOutput(outputId="ChainName")
    ),
    column(8,
      plotOutput(outputId="DDG_res_all", height = "100%")
    )
  ),
  tabsetPanel(type = "tabs",
              tabPanel("System", plotOutput(outputId="DDG_res")),
              tabPanel("GDP", plotOutput(outputId="DDG_res_GDP"))
  )
  # plotOutput(outputId="DDG_res"),
  # titlePanel("GDP binding"),
  # plotOutput(outputId="DDG_res_GDP")
  #plotOutput(outputId ="DDG_res_all_PTM", height ="100%")
)

server <- function(input, output) {
  output$DDG_res_all <- renderPlot({
    ResidueDDGplot_all(input$pdb, "all")
  }, height = 400)
  output$DDG_res_all_PTM <- renderPlot({
    ResidueDDGplot_all(input$pdb, "PTM")
  }, height =400)
  output$DDG_res <- renderPlot({
    ResidueDDGplot(input$pdb, FALSE)
  })
  output$PTM_count <- renderTable({
        x = table(Table_DDG_residues$PTM[Table_DDG_residues$PDB == input$pdb])
        x / (length(unique(Table_DDG_residues$COMP[Table_DDG_residues == input$pdb]))-1)
  })
  output$DDG_binding <- renderTable({
    subset(Table_DDG_complex, subset = PDB == input$pdb, select=c("COMP", "DDG"))
  })
  output$ChainRange <- renderTable({
    subset(Table_chain_range, PDB == input$pdb, select=c("CHAIN", "BEG", "END"))
  })
  output$ChainName <- renderTable({
    cbind(c("A","B","C","D","E","F"),t(subset(Table_Complex_Build, PDB==input$pdb)[,1:6]))
  })
  output$DDG_res_GDP <- renderPlot({
    if(sum(Table_DDG_complex$GDP[Table_DDG_complex$PDB==input$pdb]) > 0){
    ResidueDDGplot(input$pdb,TRUE)
    }
    else{
      ggplot() + theme_bw() + labs(title = "There is no GDP in this system")
    }
  })
}
shinyApp(ui=ui, server=server)
