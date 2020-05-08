#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(caret)
library(glmnet)
library(readr)
library(ggplot2)
library(dplyr)
library(limma)
library(plotROC)
library(pROC)
library(ggpubr)
library(BiocManager)
options(repos = BiocManager::repositories())



# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    Model <- reactiveValues(
      Record = matrix(ncol = 3, dimnames = list("N",c("Feature Selection Method","Number of Features", "AUC Score")))
    )  
  
    Limma <- observeEvent(input[["fit_model"]],{
      
      #Load processed datasets
      gse <- read_csv("GSE120396_expression_matrix.txt", col_names = T)
      gse <- data.frame(gse[,-1], row.names = gse$X1)
      
      rejection_status <- read_csv("rejection_status.txt", col_names = T)
      rejection_status <- as.character(rejection_status$response)
      
      
      #Processing
      processed_gse <- as.matrix(gse)
      processed_gse <- as.data.frame(processed_gse) %>% tibble::rownames_to_column("Gene")
      
      n <- input$features
      
      #Feature Selection using Limma
      
      groupname <- factor(rejection_status)
      design <- model.matrix(~ groupname + 0)
        
      fit <- lmFit(gse, design)
      cont.matrix <- makeContrasts(groupnameYes-groupnameNo, levels=design)
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit2)
        
      tT <- topTable(fit2, number = n) #select the top 200 genes from the fit
        
      #Create list of chosen features
      genes <- round(tT[1:n,], 2)
      genes <- tibble::rownames_to_column(genes, "Gene")
      features <- genes$Gene
        
      #Select those features from the dataset
      processed_gse <- processed_gse %>% filter(Gene %in% features)
      processed_gse <- data.frame(processed_gse[,-1], row.names = processed_gse[,1]) 
      processed_gse <- as.matrix(processed_gse)
      
      #Create the dataframe that will be loaded into the model
      X <- t(processed_gse)
      y <- rejection_status
      
      X_df <- as.data.frame(X)
      y_df <- as.data.frame(y)
      
      dataset <- cbind(X_df, response = y_df$y)
      
      set.seed(71)
      gse_i <- createDataPartition(dataset$response, p = 0.8, list = FALSE)
      
      x <- dataset[gse_i,] %>% select(-response)
      x <- as.matrix(x)
      y <- dataset[gse_i,] %>% select(response) %>% mutate(response = ifelse(response == "Yes", 1, 0))
      y <- as.matrix(y)
      
      cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial", type.measure='auc', nfolds = 5)
      
      # Fit the final model on the training data
      gse_model <- glmnet(x, y, alpha = 1, family = "binomial",
                          lambda = cv.lasso$lambda.min)
      
      test_x = as.matrix(dataset[-gse_i,] %>% select(-response))
      test_y = dataset[-gse_i,]%>% select(response)
      
      
      probabilities <- predict(gse_model, newx = test_x, s=cv.lasso$lambda.min, type = "response")
      probabilities <- as.data.frame(probabilities) %>% rename(prob = 1)
      pred <- predict(gse_model, newx = test_x, type = "class") 
      
      
      vector1 <- cbind(probabilities, test_y) %>% rename(obs = response) %>% mutate(obs = ifelse(obs == "Yes", 1, 0)) 
      vector1 <- cbind(vector1, pred) 
      vector1 <- vector1 %>% rename(pred = s0)
      
      gse_ggROC <- ggplot(vector1, aes(d=obs, m=prob)) + geom_roc()
      AUC = round(calc_auc(gse_ggROC)$AUC, 2)
      
      Model$Record <- rbind(Model$Record, c("Limma",input$features,AUC))
      
      return(vector1)
    })
        
      
    Random <- observeEvent(input[["fit_model"]],{
      
      #Load processed datasets
      gse <- read_csv("GSE120396_expression_matrix.txt", col_names = T)
      gse <- data.frame(gse[,-1], row.names = gse$X1)
      
      rejection_status <- read_csv("rejection_status.txt", col_names = T)
      rejection_status <- as.character(rejection_status$response)
      
      
      #Processing
      processed_gse <- as.matrix(gse)
      processed_gse <- as.data.frame(processed_gse) %>% tibble::rownames_to_column("Gene")
      
      n <- input$features
      
      #Create list of chosen features
      features <- processed_gse$Gene
      n <- input$features
      features <- sample(features, n, replace = F)
        
      #Select those features from the dataset
      processed_gse <- processed_gse %>% filter(Gene %in% features)
      processed_gse <- data.frame(processed_gse[,-1], row.names = processed_gse[,1]) 
      processed_gse <- as.matrix(processed_gse)
      
      
      #Create the dataframe that will be loaded into the model
      X <- t(processed_gse)
      y <- rejection_status
      
      X_df <- as.data.frame(X)
      y_df <- as.data.frame(y)
      
      dataset <- cbind(X_df, response = y_df$y)
      
      set.seed(71)
      gse_i <- createDataPartition(dataset$response, p = 0.8, list = FALSE)
      
      x <- dataset[gse_i,] %>% select(-response)
      x <- as.matrix(x)
      y <- dataset[gse_i,] %>% select(response) %>% mutate(response = ifelse(response == "Yes", 1, 0))
      y <- as.matrix(y)
      
      cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial", type.measure='auc', nfolds = 5)
      
      # Fit the final model on the training data
      gse_model <- glmnet(x, y, alpha = 1, family = "binomial",
                          lambda = cv.lasso$lambda.min)
      
      test_x = as.matrix(dataset[-gse_i,] %>% select(-response))
      test_y = dataset[-gse_i,]%>% select(response)
      
      
      probabilities <- predict(gse_model, newx = test_x, s=cv.lasso$lambda.min, type = "response")
      probabilities <- as.data.frame(probabilities) %>% rename(prob = 1)
      pred <- predict(gse_model, newx = test_x, type = "class") 
      
      
      vector2 <- cbind(probabilities, test_y) %>% rename(obs = response) %>% mutate(obs = ifelse(obs == "Yes", 1, 0)) 
      vector2 <- cbind(vector2, pred) 
      vector2 <- vector2 %>% rename(pred = s0)
      
      gse_ggROC <- ggplot(vector2, aes(d=obs, m=prob)) + geom_roc()
      AUC = round(calc_auc(gse_ggROC)$AUC, 2)
    
      Model$Record <- rbind(Model$Record, c("Random",input$features,AUC))
      return(vector2)
    })
  
    Variance <- observeEvent(input[["fit_model"]],{
                            
      #Load processed datasets
      gse <- read_csv("GSE120396_expression_matrix.txt", col_names = T)
      gse <- data.frame(gse[,-1], row.names = gse$X1)
      
      rejection_status <- read_csv("rejection_status.txt", col_names = T)
      rejection_status <- as.character(rejection_status$response)
      
      
      n <- input$features
      
      #Feature Selection using Variance
      
      #Select top variance
      varValue = apply(gse, 1, var, na.rm=TRUE)
      cutoffvalue = sort(varValue, decreasing = TRUE)[n]
      varid = which(varValue >= cutoffvalue)
      
      #Select those features from the dataset
      processed_gse <- gse[varid,]
      processed_gse <- as.matrix(processed_gse)
      processed_gse
      
      #Create the dataframe that will be loaded into the model
      X <- t(processed_gse)
      y <- rejection_status
      
      X_df <- as.data.frame(X)
      y_df <- as.data.frame(y)
      
      dataset <- cbind(X_df, response = y_df$y)
      
      set.seed(71)
      gse_i <- createDataPartition(dataset$response, p = 0.8, list = FALSE)
      dataset
      
      x <- dataset[gse_i,] %>% select(-response)
      x <- as.matrix(x)
      y <- dataset[gse_i,] %>% select(response) %>% mutate(response = ifelse(response == "Yes", 1, 0))
      y <- as.matrix(y)
      
      cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial", type.measure='auc', nfolds = 5)
      
      # Fit the final model on the training data
      gse_model <- glmnet(x, y, alpha = 1, family = "binomial",
                          lambda = cv.lasso$lambda.min)
      
      test_x = as.matrix(dataset[-gse_i,] %>% select(-response))
      test_y = dataset[-gse_i,]%>% select(response)
      
      
      probabilities <- predict(gse_model, newx = test_x, s=cv.lasso$lambda.min, type = "response")
      probabilities <- as.data.frame(probabilities) %>% rename(prob = 1)
      pred <- predict(gse_model, newx = test_x, type = "class") 
      
      
      vector3 <- cbind(probabilities, test_y) %>% rename(obs = response) %>% mutate(obs = ifelse(obs == "Yes", 1, 0)) 
      vector3 <- cbind(vector3, pred) 
      vector3 <- vector3 %>% rename(pred = s0)
      
      gse_ggROC <- ggplot(vector3, aes(d=obs, m=prob)) + geom_roc()
      AUC = round(calc_auc(gse_ggROC)$AUC, 2)
      
      Model$Record <- rbind(Model$Record, c("Variance",input$features,AUC))
        
      Model$Record <- na.omit(Model$Record)
      return(vector3)
        
        
      })

    feature1 <- reactive({
      
      #Load processed datasets
      gse <- read_csv("GSE120396_expression_matrix.txt", col_names = T)
      gse <- data.frame(gse[,-1], row.names = gse$X1)
      
      rejection_status <- read_csv("rejection_status.txt", col_names = T)
      rejection_status <- as.character(rejection_status$response)
      
      
      #Processing
      processed_gse <- as.matrix(gse)
      processed_gse <- as.data.frame(processed_gse) %>% tibble::rownames_to_column("Gene")
      
      n <- input$features
      
      #Feature Selection using Limma
      
      groupname <- factor(rejection_status)
      design <- model.matrix(~ groupname + 0)
      
      fit <- lmFit(gse, design)
      cont.matrix <- makeContrasts(groupnameYes-groupnameNo, levels=design)
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit2)
      
      tT <- topTable(fit2, number = n) #select the top 200 genes from the fit
      
      #Create list of chosen features
      genes <- round(tT[1:n,], 2)
      genes <- tibble::rownames_to_column(genes, "Gene")
      features <- genes$Gene
      
      #Select those features from the dataset
      processed_gse <- processed_gse %>% filter(Gene %in% features)
      processed_gse <- data.frame(processed_gse[,-1], row.names = processed_gse[,1]) 
      processed_gse <- as.matrix(processed_gse)
      
      #Create the dataframe that will be loaded into the model
      X <- t(processed_gse)
      y <- rejection_status
      
      X_df <- as.data.frame(X)
      y_df <- as.data.frame(y)
      
      dataset <- cbind(X_df, response = y_df$y)
      
      set.seed(71)
      gse_i <- createDataPartition(dataset$response, p = 0.8, list = FALSE)
      
      x <- dataset[gse_i,] %>% select(-response)
      x <- as.matrix(x)
      y <- dataset[gse_i,] %>% select(response) %>% mutate(response = ifelse(response == "Yes", 1, 0))
      y <- as.matrix(y)
      
      cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial", type.measure='auc', nfolds = 5)
      
      # Fit the final model on the training data
      gse_model <- glmnet(x, y, alpha = 1, family = "binomial",
                          lambda = cv.lasso$lambda.min)
      
      test_x = as.matrix(dataset[-gse_i,] %>% select(-response))
      test_y = dataset[-gse_i,]%>% select(response)
      
      
      probabilities <- predict(gse_model, newx = test_x, s=cv.lasso$lambda.min, type = "response")
      probabilities <- as.data.frame(probabilities) %>% rename(prob = 1)
      pred <- predict(gse_model, newx = test_x, type = "class") 
      
      
      vector1 <- cbind(probabilities, test_y) %>% rename(obs = response) %>% mutate(obs = ifelse(obs == "Yes", 1, 0)) 
      vector1 <- cbind(vector1, pred) 
      vector1 <- vector1 %>% rename(pred = s0)
      
      return(vector1)
    })
    
    
    feature2 <- reactive({
      
      #Load processed datasets
      gse <- read_csv("GSE120396_expression_matrix.txt", col_names = T)
      gse <- data.frame(gse[,-1], row.names = gse$X1)
      
      rejection_status <- read_csv("rejection_status.txt", col_names = T)
      rejection_status <- as.character(rejection_status$response)
      
      
      #Processing
      processed_gse <- as.matrix(gse)
      processed_gse <- as.data.frame(processed_gse) %>% tibble::rownames_to_column("Gene")
      
      n <- input$features
      
      #Create list of chosen features
      features <- processed_gse$Gene
      n <- input$features
      features <- sample(features, n, replace = F)
      
      #Select those features from the dataset
      processed_gse <- processed_gse %>% filter(Gene %in% features)
      processed_gse <- data.frame(processed_gse[,-1], row.names = processed_gse[,1]) 
      processed_gse <- as.matrix(processed_gse)
      
      
      #Create the dataframe that will be loaded into the model
      X <- t(processed_gse)
      y <- rejection_status
      
      X_df <- as.data.frame(X)
      y_df <- as.data.frame(y)
      
      dataset <- cbind(X_df, response = y_df$y)
      
      set.seed(71)
      gse_i <- createDataPartition(dataset$response, p = 0.8, list = FALSE)
      
      x <- dataset[gse_i,] %>% select(-response)
      x <- as.matrix(x)
      y <- dataset[gse_i,] %>% select(response) %>% mutate(response = ifelse(response == "Yes", 1, 0))
      y <- as.matrix(y)
      
      cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial", type.measure='auc', nfolds = 5)
      
      # Fit the final model on the training data
      gse_model <- glmnet(x, y, alpha = 1, family = "binomial",
                          lambda = cv.lasso$lambda.min)
      
      test_x = as.matrix(dataset[-gse_i,] %>% select(-response))
      test_y = dataset[-gse_i,]%>% select(response)
      
      
      probabilities <- predict(gse_model, newx = test_x, s=cv.lasso$lambda.min, type = "response")
      probabilities <- as.data.frame(probabilities) %>% rename(prob = 1)
      pred <- predict(gse_model, newx = test_x, type = "class") 
      
      
      vector2 <- cbind(probabilities, test_y) %>% rename(obs = response) %>% mutate(obs = ifelse(obs == "Yes", 1, 0)) 
      vector2 <- cbind(vector2, pred) 
      vector2 <- vector2 %>% rename(pred = s0)
      
      
      return(vector2)
    })
    
    feature3 <- reactive({
      
      #Load processed datasets
      gse <- read_csv("GSE120396_expression_matrix.txt", col_names = T)
      gse <- data.frame(gse[,-1], row.names = gse$X1)
      
      rejection_status <- read_csv("rejection_status.txt", col_names = T)
      rejection_status <- as.character(rejection_status$response)
      
      
      n <- input$features
      
      #Feature Selection using Variance
      
      #Select top variance
      varValue = apply(gse, 1, var, na.rm=TRUE)
      cutoffvalue = sort(varValue, decreasing = TRUE)[n]
      varid = which(varValue >= cutoffvalue)
      
      #Select those features from the dataset
      processed_gse <- gse[varid,]
      processed_gse <- as.matrix(processed_gse)
      
      #Create the dataframe that will be loaded into the model
      X <- t(processed_gse)
      y <- rejection_status
      
      X_df <- as.data.frame(X)
      y_df <- as.data.frame(y)
      
      dataset <- cbind(X_df, response = y_df$y)
      
      set.seed(71)
      gse_i <- createDataPartition(dataset$response, p = 0.8, list = FALSE)
      
      x <- dataset[gse_i,] %>% select(-response)
      x <- as.matrix(x)
      y <- dataset[gse_i,] %>% select(response) %>% mutate(response = ifelse(response == "Yes", 1, 0))
      y <- as.matrix(y)
      
      cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial", type.measure='auc', nfolds = 5)
      
      # Fit the final model on the training data
      gse_model <- glmnet(x, y, alpha = 1, family = "binomial",
                          lambda = cv.lasso$lambda.min)
      
      test_x = as.matrix(dataset[-gse_i,] %>% select(-response))
      test_y = dataset[-gse_i,]%>% select(response)
      
      
      probabilities <- predict(gse_model, newx = test_x, s=cv.lasso$lambda.min, type = "response")
      probabilities <- as.data.frame(probabilities) %>% rename(prob = 1)
      pred <- predict(gse_model, newx = test_x, type = "class") 
      
      
      vector3 <- cbind(probabilities, test_y) %>% rename(obs = response) %>% mutate(obs = ifelse(obs == "Yes", 1, 0)) 
      vector3 <- cbind(vector3, pred) 
      vector3 <- vector3 %>% rename(pred = s0)
      
      
      return(vector3)
      
      
    })
  
  
      
      
      
    output$`ROC Curve` <- renderPlot({
      req(input$fit_model)
      
      vector1 <- feature1()
      vector2 <- feature2()
      vector3 <- feature3()

      

      gse_ggROC1 <- ggplot(vector1, aes(d=obs, m=prob)) + geom_roc(colour = "darksalmon") + style_roc(theme = theme_gray()) + ggtitle("Feature Selection Using Limma")
      gse_ggROC_styled1 <- gse_ggROC1 +  annotate("text", x = .75, y = .25, 
                                                  label = paste("AUC =", round(calc_auc(gse_ggROC1)$AUC, 2)))
      
      gse_ggROC2 <- ggplot(vector2, aes(d=obs, m=prob)) + geom_roc(colour = "dodgerblue4") + style_roc(theme = theme_gray()) + ggtitle("Feature Selection Using Random Selection")
      gse_ggROC_styled2 <- gse_ggROC2 +  annotate("text", x = .75, y = .25, 
                                                  label = paste("AUC =", round(calc_auc(gse_ggROC2)$AUC, 2)))
      
      gse_ggROC3 <- ggplot(vector3, aes(d=obs, m=prob)) + geom_roc(colour = "darkslategray") + style_roc(theme = theme_gray()) + ggtitle("Feature Selection Using Ordered Variance")
      gse_ggROC_styled3 <- gse_ggROC3 +  annotate("text", x = .75, y = .25, 
                                                  label = paste("AUC =", round(calc_auc(gse_ggROC3)$AUC, 2)))
      
      library(ggpubr)
      figure <- ggarrange(gse_ggROC_styled1, gse_ggROC_styled2, gse_ggROC_styled3, ncol = 2, nrow = 2)
      return(figure)
    }, height = 800, width = 800)
    
    output$model_record <- 
      renderTable({
        req(input$fit_model)
        mat <- Model$Record
        mat = as.data.frame(mat) %>% tidyr::pivot_wider(names_from = `Feature Selection Method`, values_from = `AUC Score`)
        return(mat)
      })

})
