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
    titlePanel("Comparison of Feature Selection Methods"),

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
        h3("Effect of the Number of Features on Classifier"),
        tableOutput("model_record")
        
        ),

        # Show a plot of the generated distribution
        column(
          width = 7,
          tabsetPanel(
            tabPanel("Overview",
                     h3("Feature Selection Methods"),
                     h4("Limma"),
                     ("The limma method involves performing multiple t tests moderated with Empirical Bayes using the LIMMA package. These t tests investigate whether all gene-wise contrasts are zero. A linear model was fitted using the RNA-sequencing data and a design matrix with known entries of 0 or 1 corresponding to the patient outcome. However, since numerous tests were executed, the probability of observing a significant result due to chance is high. Hence, using the eBayes function, the test statistics were moderated with Empirical Bayes which estimates the prior from all features.  This method assumes that the variance comes from a gamma distribution. Each gene’s calculated variance is used to estimate the mean and variance of this distribution. The gamma distribution represents the empirical prior distribution for the variances. Using this moderated variance, the degrees of freedom increased to reflect the additional information from the prior resulting in increased statistical power to detect differential expression. "),
                     h4("Ordered Variance"),
                     ("As the variance of each gene vary greatly, the variance of each gene was calculated. The top n genes with the highest variance were selected. Note that n stands for the number of feature selected by the user."),
                     h4("Random"),
                     ("As a baseline to compare the other feature selection methods, in this method, genes were chosen randomly from the set of all genes."),
                     h3("LASSO Logistic Regression Model"),
                     ("Logistic regression models the class probabilities. Any input produces a number between 0 and 1 which represents the probability that the patient will reject their kidney transplant. A penalty term was added to the model to shrink the regression coefficient of the features that contribute most to the error to zero. This is called LASSO regularisation.")
            ),
            tabPanel("Model",
                     h3("Effect of Feature Selection Methods on Model Performance"),
                     plotOutput("ROC Curve")
                     
            )
          )
        )
))