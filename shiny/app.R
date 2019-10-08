# Load packages ----
library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(shinydashboard)
library(hypeR)
library(clusterProfiler)
library(gridExtra)
library(tibble)


# User interface ----
ui <- dashboardPage(
      dashboardHeader(title = "volcano plots"),
      
      dashboardSidebar(
        sidebarMenu(id="tabs",
                    menuItem(
                      "DGE",
                      tabName = "dge",
                      icon = icon("bar-chart-o")),
                    conditionalPanel(
                      "input.tabs == 'dge'",
                      radioButtons("res1", "Comparisons", 
                                   list("HSC CS v Adult",
                                        "IL7R CS v Adult",
                                        "proB CS v Adult",
                                        "HSC W12-13 v Adult",
                                        "IL7R W12-13 v Adult",
                                        "ProB W12-13 CS v Adult"),
                                   select = "HSC CS v Adult")
                    ),
                    menuItem("PCA", tabName = "pca_tab", icon = icon("barcode"))
        )
       ),
      dashboardBody(tabItems(tabItem("dge",
          fluidRow(
              box( title = "DGE", DT::dataTableOutput("tbl"), width = 12),
              box( title = "Volcano", plotOutput("volcano_out"),
                   downloadButton("download_volcano","Download")),
              box( title = "Halmark GSEA", plotOutput('hallmark_out'),
                   downloadButton("download_hallmark","Download")),
              box( title = "Laurenti GSEA", plotOutput('laurenti_out'),
                   downloadButton("download_laurenti","Download")),
              box( title = "GO analysis", plotOutput('ego_out'), width = 12,
                  downloadButton("download_ego","Download"))
         )),
         tabItem("pca_tab",
                 fluidRow(
                   box (title = "PCA", plotOutput("pca_plot")), 
                   downloadButton("download_pca","Download")
                 )
                 )
      )
   )
)

plot_volcano <- function(res,row_sel) {
    res <- res %>% rownames_to_column(var = "symbol")
    res <- res %>% dplyr::mutate(sig = case_when( padj > 0.05 ~  1, padj < 0.05 ~ 0))
    res$sig <- as.factor(res$sig)
    ggplot(res,aes(log2FoldChange,-log10(padj),colour = sig, label = symbol)) + geom_point() +
        ylab(expression(paste("-", log[10], "(padj)", sep=""))) + scale_colour_manual(values = c("red","black")) +
        xlab(expression(paste( log[2], "(fold change)", sep=""))) + theme_bw() + theme(legend.position = "none")  +
        geom_label(data = subset(res[row_sel,]))
}

pca_plot <- readRDS('data/pca_plot.rds')

server <- function(input, output) {
    rv <- reactiveValues(results=character(),
                         res=data.frame(),
                         hallmark=character(),
                         laurenti=character(),
                         ego=character()
                         )
    observeEvent(input$res1,{
      if(input$res1=="HSC CS v Adult"){
      rv$results<-'data/hsc_cs_v_adult.rds'
      rv$hallmark <- readRDS('data/hsc_cs_v_adult_hallmark.rds')
      rv$laurenti <- readRDS('data/hsc_cs_v_adult_laurenti.rds')
      rv$ego <- readRDS('data/hsc_cs_v_adult_ego.rds')
      }
      else if(input$res1=="IL7R CS v Adult"){
        rv$results<-'data/il7r_cs_v_adult.rds'
        rv$hallmark <- readRDS('data/il7r_cs_v_adult_hallmark.rds')
        rv$laurenti <- readRDS('data/il7r_cs_v_adult_laurenti.rds')
        rv$ego <- readRDS('data/il7r_cs_v_adult_ego.rds')
        }
        else if(input$res1=="proB CS v Adult"){
          rv$results<-'data/proB_cs_v_adult.rds'
          rv$hallmark <- readRDS('data/proB_cs_v_adult_hallmark.rds')
          rv$laurenti <- readRDS('data/proB_cs_v_adult_laurenti.rds')
          rv$ego <- readRDS('data/proB_cs_v_adult_ego.rds')
          }
        else if(input$res1=="HSC W12-13 v Adult"){
          rv$results<-'data/hsc_w12_13_v_adult.rds'
          rv$hallmark <- readRDS('data/hsc_w12_13_v_adult_hallmark.rds')
          rv$laurenti <- readRDS('data/hsc_w12_13_v_adult_laurenti.rds')
          rv$ego <- readRDS('data/hsc_w12_13_v_adult_ego.rds')
          }
        else if(input$res1=="IL7R W12-13 v Adult"){
          rv$results<-'data/il7r_w12_13_v_adult.rds'
          rv$hallmark <- readRDS('data/il7r_w12_13_v_adult_hallmark.rds')
          rv$laurenti <- readRDS('data/il7r_w12_13_v_adult_laurenti.rds')
          rv$ego <- readRDS('data/il7r_w12_13_v_adult_ego.rds')
          }
        else if(input$res1=="ProB W12-13 CS v Adult"){
          rv$results<-'data/proB_w12_13_v_adult.rds'
          rv$hallmark <- readRDS('data/proB_w12_13_v_adult_hallmark.rds')
          rv$laurenti <- readRDS('data/proB_w12_13_v_adult_laurenti.rds')
          rv$ego <- readRDS('data/proB_w12_13_v_adult_ego.rds')
          }
        rv$res <- readRDS(rv$results)
     })
    output$tbl <- DT::renderDataTable(rv$res, selection = list(target = 'row' , selected = 1))
    output$volcano_out <- renderPlot(plot_volcano(rv$res,input$tbl_rows_selected))
    output$hallmark_out <- renderPlot(hyp_dots(rv$hallmark, show_plots=FALSE, return_plots=TRUE))
    output$laurenti_out <- renderPlot(hyp_dots(rv$laurenti))
    output$ego_out <- renderPlot(barplot(rv$ego, showCategory=20))
    
    output$pca_plot <- renderPlot(grid.arrange(pca_plot))
    
    output$download_volcano <- downloadHandler(
        filename = function(){
            paste("volcano",".pdf",sep="")},
        content = function(file) {
            ggsave(file,plot = plot_volcano(rv$res,input$tbl_rows_selected),device = "pdf")
        }
    )
  
    output$download_hallmark <- downloadHandler(
      filename = function(){
        paste("hallmark",".pdf",sep="")},
      content = function(file) {
        ggsave(file,plot = hyp_dots(rv$hallmark, show_plots=FALSE, return_plots=TRUE),device = "pdf")
      }
    )
    
    output$download_laurenti <- downloadHandler(
      filename = function(){
        paste("laurenti",".pdf",sep="")},
      content = function(file) {
        ggsave(file,plot = hyp_dots(rv$laurenti),device = "pdf")
      }
    )
    
    output$download_ego <- downloadHandler(
      filename = function(){
        paste("go_analysis",".pdf",sep="")},
      content = function(file) {
        ggsave(file,plot = barplot(rv$ego, showCategory=20),device = "pdf")
      }
    )
    output$download_pca <- downloadHandler(
      filename = function(){
        paste("pca",".pdf",sep="")},
      content = function(file) {
        ggsave(file,plot = grid.arrange(pca_plot),device = "pdf")
      }
    )
}

# Run app ----
shinyApp(ui, server)
