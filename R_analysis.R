library(TCGAbiolinks)
library(tidyverse)
library(plotly)
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(dplyr)

# create theme for plots
mynamestheme <- theme(
  plot.title = element_text(family = "Helvetica", face = "bold", size = (15)),
  legend.title = element_text(colour = "steelblue", face = "bold.italic", family = "Helvetica"),
  legend.text = element_text(face = "italic", colour = "steelblue4", family = "Helvetica"),
  axis.title = element_text(family = "Helvetica", size = (10), colour = "steelblue4"),
  axis.text = element_text(family = "Courier", colour = "cornflowerblue", size = (10))
)

# read clinic data
clinic_1 = read_delim(file="clinical.csv", delim = "\t")


# downloading clinical data
clin_query <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Clinical",
                       file.type = "xml")

 GDCdownload(clin_query)
# get the clinic dataframe
clinic_2 <- GDCprepare_clinic(clin_query, clinical.info = "patient")



colnames(clinic_2)[1] <- "case_submitter_id"
clinic <- merge(x = clinic_2, y = clinic_1[ , c("case_submitter_id", "ajcc_pathologic_stage")], by = "case_submitter_id", all.x=TRUE)

clinic <- clinic %>%
  filter(gender == "FEMALE")


# extract empty value in age
clinic_age_no_NAs <- clinic[!(clinic$age_at_initial_pathologic_diagnosis == "" | is.na(clinic$age_at_initial_pathologic_diagnosis)), ]

# create age_category variable
clinic_age_no_NAs$age_category = ifelse(clinic_age_no_NAs$age_at_initial_pathologic_diagnosis <50, "young", "old")

# create graph to see the distribution of stages in each category
distribution_graph = ggplot(clinic_age_no_NAs, aes(age_category, after_stat(count))) + geom_bar(aes(fill = ajcc_pathologic_stage), position = "dodge") + scale_fill_viridis_d() + labs(title="Distribution of pathologic stages in each age category", x = "age category", y= "count") + mynamestheme

ggplotly(distribution_graph)


# survival analysis
sum_exp_dataframe_age<- as.data.frame(clinic_age_no_NAs)

# replace the NA balue in days_to_death to days_to_last_follow_up
sum_exp_dataframe_age$days_to_death = ifelse(is.na(sum_exp_dataframe_age$days_to_death), 
                                             sum_exp_dataframe_age$days_to_last_followup, sum_exp_dataframe_age$days_to_death)

# make the days_to_death column numeric 
sum_exp_dataframe_age$days_to_death = as.numeric(sum_exp_dataframe_age$days_to_death)

# create the death_event column here
sum_exp_dataframe_age$death_event = ifelse(sum_exp_dataframe_age$vital_status == "Alive", 0, 1)

# examine the relationship between age, lymph_node_count based on vital_status
lymph_node_no_NAs <- sum_exp_dataframe_age[!(sum_exp_dataframe_age$lymph_node_examined_count == "" | is.na(sum_exp_dataframe_age$lymph_node_examined_count)), ]

lymph_node_no_NAs %>% 
  plot_ly(x = ~age_at_initial_pathologic_diagnosis, y = ~lymph_node_examined_count, 
          type = 'scatter',
          mode = 'markers',
          color = ~vital_status,
          sizes = c(5, 70),
          marker = list(sizemode='diameter', opacity=0.5), hoverinfo = 'text',
          text = ~paste( 
                         paste("Age: ",age_at_initial_pathologic_diagnosis , sep="") ,
                         paste("Lymph Node Count: ", lymph_node_examined_count, sep=""),
                         paste("Vital Status: ", vital_status,sep=""),
                         sep = "<br>")) %>%
  layout(title = "Ages, Lymph Node Count for Female Breast Cancer Patients",
         xaxis = list(title = "age at diagnosis"), 
         yaxis = list(title = "lymph node count"),
         hovermode = "compare")
  




# Replace N/A value in person_neoplasm_cancer_status to "unknown"
sum_exp_dataframe_age <- sum_exp_dataframe_age %>% mutate_at(c('person_neoplasm_cancer_status'), ~na_if(., ''))
sum_exp_dataframe_age$person_neoplasm_cancer_status <- as.character(sum_exp_dataframe_age$person_neoplasm_cancer_status)
sum_exp_dataframe_age$person_neoplasm_cancer_status <- sum_exp_dataframe_age$person_neoplasm_cancer_status %>% replace_na('UNKNOWN')

sum_exp_dataframe_age$person_neoplasm_cancer_status<-gsub(" ", "_", sum_exp_dataframe_age$person_neoplasm_cancer_status)

# examine different status of tumor at different age 
with_tumor <- sum_exp_dataframe_age %>%
  filter(person_neoplasm_cancer_status == "WITH_TUMOR")

tumor_free <- sum_exp_dataframe_age %>%
  filter(person_neoplasm_cancer_status == "TUMOR_FREE")

# We initialize a 'survival' object first, which contains the data we need.
surv_object_with_tumor <- Surv(time = with_tumor$days_to_death, 
                               event = with_tumor$death_event)

surv_object_tumor_free <- Surv(time = tumor_free$days_to_death, 
                               event = tumor_free$death_event)

# We then create a fit object to do the gender analysis 
age_fit_with_tumor <- surv_fit( surv_object_with_tumor ~ with_tumor$age_category, data = with_tumor )

age_fit_tumor_free <- surv_fit( surv_object_tumor_free ~ tumor_free$age_category, data = tumor_free )

#the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot_with_tumor = ggsurvplot(age_fit_with_tumor, 
                                 legend = "right",
                                 legend.title = c("Age"),
                                 legend.labs = c("old (N=61)", "young (N=38)"),
                                 title="Kaplan-Meier Curve for Breast Cancer Survival with tumor (N=99)",
                                 ggtheme = theme_bw()) 

plotly::ggplotly(survplot_with_tumor[[1]])

survplot_tumor_free = ggsurvplot(age_fit_tumor_free, 
                                 legend = "right",
                                 legend.title = c("Age"),
                                 legend.labs = c("old (N=690)", "young (N=237)"),
                                 title="Kaplan-Meier Curve for Breast Cancer Survival without tumor (N=928)",
                                 ggtheme = theme_bw()) 

plotly::ggplotly(survplot_tumor_free[[1]])
