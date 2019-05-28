library(shiny)

ui <- fluidPage(
  selectInput("unp", label="Choose Uniprot ID", choices=c("all",unique(Table_UNP_EMBL$UNP))),
  titlePanel("Conservation analysis"),
  plotOutput(outputId="boxPlot_Res_Cons"),
  # plotOutput(outputId = "smoothPlot_ConsBck"),
  # plotOutput(outputId = "boxPlot_PTM_ConsGen"),
  # plotOutput(outputId = "boxPlot_PTM_ConsAdj"),
  # plotOutput(outputId = "boxPlot_PTM_ConsRel_Adj"),
  # plotOutput(outputId = "boxPlot_PTM_ConsRel_Sha"),
  plotOutput(outputId = "Boxplot_PTM_ConsDSSP")
  
)

server <- function(input, output) {
  output$boxPlot_Res_Cons <- renderPlot({
    if (input$unp == "all"){
      df = Table_conservation
    }
    else{
    df = subset(Table_conservation, UNP == input$unp)
    }
    boxPlot_Res_Cons(df)
  })
  output$boxPlot_PTM_ConsGen <- renderPlot({
    if (input$unp == "all"){
      df = Table_conservation
    }
    else{
    df = subset(Table_conservation, UNP == input$unp)
    }
    boxPlot_PTM_ConsGen(df)
  })
  output$boxPlot_PTM_ConsAdj <- renderPlot({
    if (input$unp == "all"){
      df = Table_conservation
    }
    else{
      df = subset(Table_conservation, UNP == input$unp)
    }
    boxPlot_PTM_ConsAdj(df)
  })
  output$boxPlot_PTM_ConsRel_Adj <- renderPlot({
    if (input$unp == "all"){
      df = Table_conservation
    }
    else{
      df = subset(Table_conservation, UNP == input$unp)
    }
    boxPlot_PTM_ConsRel(df)
  })
  output$boxPlot_PTM_ConsRel_Sha <- renderPlot({
    if (input$unp == "all"){
      df = Table_conservation
    }
    else{
      df = subset(Table_conservation, UNP == input$unp)
    }
    boxPlot_PTM_ConsRel_Bck(df)
  })
  output$Boxplot_PTM_ConsDSSP<- renderPlot({
    if (input$unp == "all"){
      Boxplot_PTM_ConsDSSP(Table_conservation)
    }
    else{
      df = subset(Table_conservation, UNP == input$unp)
      Boxplot_PTM_ConsDSSP(df)
    }
  })
  
}
shinyApp(ui=ui, server=server)
