#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    
    #CSS Styles
    tags$style(HTML(
      paste0(".shiny-input-container {background-color: #f5f5f5; border: 1px solid #e3e3e3;",
             "padding-left: 10px; padding-right: 10px; padding-top: 10px; padding-bottom: 10px; border-radius: 3px;}")
    )),
  
    
    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins
    column(
        width = 5,
          
        sliderInput("features",
                        "Number of Features:",
                        min = 2,
                        max = 1000,
                        value = 30),
        
        actionButton(inputId = "fit_model",
                     label = "Fit Model"),
        
        hr(),
        tableOutput("model_record")
        
        ),

        # Show a plot of the generated distribution
        column(
          width = 7,
          tabsetPanel(
            tabPanel("Overview",
                     h2("Title"),
                     h3("Feature Selection Methods"),
                     h4("Limma"),
                     h5("Explanation of Limma"),
                     h4("Variance"),
                     h5("Explanation of Variance"),
                     h4("Random"),
                     h5("Explanation of Random"),
                     h3("Regularised Logistic Regression"),
                     h4("Explanation of Model")
            ),
            tabPanel("Model",
                     h3("Model Accuracy"),
                     plotOutput("ROC Curve")
                     
            )
          )
          
          
          
          #textOutput("AUC_score")
        )
))