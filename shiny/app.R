# Load packages ----
library(shiny)
library(DT)
library(ggplot2)
library(dplyr)

##
# Load data ----
data(res)

# User interface ----
ui <- fluidPage(
    title = "volcano, DGE and FPKMs",
    fluidPage(DT::dataTableOutput("tbl"),
              textOutput("text_out"),
              shinydashboard::box( plotOutput("volcano_out"),
                   downloadButton("download_volcano","Download"))
    )
)

plot_volcano <- function(x) {
    res <- res %>% dplyr::mutate(sig = case_when( padj > 0.05 ~  1, padj < 0.05 ~ 0))
    res$sig <- as.factor(res$sig)
    ggplot(res,aes(log2FoldChange,-log10(padj),colour = sig, label = symbol)) + geom_point() +
        ylab(expression(paste("-", log[10], "(padj)", sep=""))) + scale_colour_manual(values = c("red","black")) +
        xlab(expression(paste( log[2], "(fold change)", sep=""))) + theme_bw() + theme(legend.position = "none") +
        geom_label(data = subset(res[x,]))
}

server <- function(input, output) {
    output$tbl <- DT::renderDataTable(res, selection = list(target = 'row' , selected = 1))
    output$volcano_out <- renderPlot(plot_volcano(input$tbl_rows_selected))
    output$text_out <- renderText(input$tbl_rows_selected)
    output$download_volcano <- downloadHandler(
        filename = function(){
            paste("volcano",".pdf",sep="")},
        content = function(file) {
            ggsave(file,plot = plot_volcano(input$tbl_rows_selected),device = "pdf")
        }
    )
}

# Run app ----
shinyApp(ui, server)
