library(tidyverse)

# Home directory
dir_home <- here::here() %>% dirname()  %>% dirname()  %>% dirname() 

#Data location of metabolomics
dir_metab <- fs::path(dir_home,
                      "2_Cleaned Data",
                      "Metabolomics") 

# Data location of outcomes and covariates
dir_tl_data <- fs::path(dir_home,
                     "2_Cleaned Data") 

# Analysis data
dir_data <- fs::path(here::here() %>% dirname(),
                     "0_data")

#Results for Mummichog
dir_results <- fs::path(dir_data %>% dirname(),
                        "2_results")

#figures
dir_figures <- fs::path(dir_data %>% dirname(),
                        "3_figures")


