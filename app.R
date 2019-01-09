library(seqinr)
library(dplyr)
library(ggplot2)
library(DT)
library(shiny)
library(biogram)
library(reshape2)

source("generate-prop-df.R")

ui <- fluidPage(
   
   titlePanel("Amino acid properties"),
   fluidRow(
     column(width = 5, 
            dataTableOutput("prop-table"),
            h3("Selected properties"),
            verbatimTextOutput("prop-vector")
     ),
     column(width = 7, 
            plotOutput("prop-hm", width = "950px", height = "950px")
     )
   )
   
)

server <- function(input, output) {
  
   output[["prop-table"]] <- renderDataTable({
     datatable(prop_df, style = "bootstrap", 
               rownames = FALSE, filter = "top", 
               extensions = "Buttons",
               options = list(pageLength = 10, dom = "Brtip",
                              buttons = c("copy", "csv", "excel", "print")))
   })
   
   
   cor_df <- reactive({
     validate(
       need(length(input[["prop-table_rows_selected"]]) > 2, "Select at least three properties from the table on the right")
     )
     
     selected_props <- prop_df[input[["prop-table_rows_selected"]], "ID"]
     cors_grid <- expand.grid(selected_props, selected_props) %>% 
       split(f = 1L:nrow(.)) %>% 
       lapply(function(i) as.character(unlist(i)))
     
     safe_cor_test <- function(x, y) {
       res <- try(cor.test(x, y, method = "pearson")[["p.value"]], silent = TRUE)
       
       if(class(res) == "try-error") {
         NA
       } else {
         res
       }
     }
     
     all_cors <- lapply(cors_grid, function(i)
       data.frame(feature1 = i[1],
                  feature2 = i[2], 
                  pearson = cor(aaprop[i[1], ], aaprop[i[2], ], method = "pearson"),
                  pval = safe_cor_test(aaprop[i[1], ], aaprop[i[2], ]),
                  stringsAsFactors = FALSE)) %>% 
       do.call(rbind, .) %>% 
       arrange(feature1, feature2)
     
     par_order <- all_cors %>% 
       select(-pval) %>% 
       dcast(feature1 ~ feature2) %>%
       select(-feature1) %>% 
       dist(method = "euclidean") %>%
       hclust(method = "ward.D2") %>% 
       getElement("order")
     
     all_cors %>% 
       right_join(combn(x = unique(all_cors[, 1])[par_order], m = 2, simplify = TRUE) %>%
                    t %>%
                    data.frame(),
                  by = c("feature1" = "X1", "feature2" = "X2")) %>%
       mutate(feature1 = factor(feature1, levels = unique(all_cors[, 1])[par_order]),
              feature2 = factor(feature2, levels = rev(unique(all_cors[, 1])[par_order])),
              pval = p.adjust(pval, method = "fdr"))
   })
   
   output[["prop-hm"]] <- renderPlot({
     ggplot(cor_df(), aes(x = feature1, y = feature2, fill = pearson, color = pval < 0.05)) +
       geom_tile(size = 2) +
       #geom_point() +
       scale_fill_gradient2() +
       scale_color_manual(values = c(NA, "black", NA)) +
       theme_bw(base_size = 14) + 
       theme(axis.text.x = element_text(angle = 90, hjust = 1),
             legend.position = "bottom")
   })
   
   output[["prop-vector"]] <- renderPrint({
     dput(prop_df[input[["prop-table_rows_selected"]], "ID"])
   })
   
}

shinyApp(ui = ui, server = server)

