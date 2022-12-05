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

# GDCdownload(clin_query)
# get the clinic dataframe
clinic_2 <- GDCprepare_clinic(clin_query, clinical.info = "patient")



colnames(clinic_2)[1] <- "case_submitter_id"
clinic <- merge(x = clinic_2, y = clinic_1[ , c("case_submitter_id", "ajcc_pathologic_stage")], by = "case_submitter_id", all.x=TRUE)

clinic <- clinic %>%
  filter(gender == "FEMALE") %>%
  filter(ajcc_pathologic_stage !="\'--" & ajcc_pathologic_stage !="Stage X")


# extract empty value in age
clinic_age_no_NAs <- clinic[!(clinic$age_at_initial_pathologic_diagnosis == "" | is.na(clinic$age_at_initial_pathologic_diagnosis)), ]

# create age_category variable
clinic_age_no_NAs$age_category = ifelse(clinic_age_no_NAs$age_at_initial_pathologic_diagnosis <62, "young", "old")

# create graph to see the distribution of stages in each category
distribution_graph = ggplot(clinic_age_no_NAs, aes(age_category, after_stat(count))) + geom_bar(aes(fill = ajcc_pathologic_stage), position = "dodge") + scale_fill_viridis_d() + labs(title="Distribution of pathologic stages in each age category", x = "age category", y= "count") + mynamestheme

ggplotly(distribution_graph)

# create graph to see the distribution of vital status in each stage
distribution_graph_2 = ggplot(clinic_age_no_NAs, aes(ajcc_pathologic_stage, after_stat(count), fill = vital_status)) + geom_bar(position="fill", stat="count") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5,hjust = 0.5)) + labs(title="Distribution of vital status in each pathologic stage", x = "pathologic stage", y= "percentage") + mynamestheme

ggplotly(distribution_graph_2)


# survival analysis
sum_exp_dataframe_age<- as.data.frame(clinic_age_no_NAs)

# replace the NA balue in days_to_death to days_to_last_follow_up
sum_exp_dataframe_age$days_to_death = ifelse(is.na(sum_exp_dataframe_age$days_to_death), 
                                             sum_exp_dataframe_age$days_to_last_followup, sum_exp_dataframe_age$days_to_death)

# make the days_to_death column numeric 
sum_exp_dataframe_age$days_to_death = as.numeric(sum_exp_dataframe_age$days_to_death)

# create the death_event column here
sum_exp_dataframe_age$death_event = ifelse(sum_exp_dataframe_age$vital_status == "Alive", 0, 1)

# examine the relationship between age, pathologic_stage based on vital_status
sum_exp_dataframe_age %>%
  plot_ly(x = ~age_at_initial_pathologic_diagnosis, y = ~ajcc_pathologic_stage, 
          type = 'scatter',
          mode = 'markers',
          color = ~vital_status,
          sizes = c(5, 70),
          marker = list(sizemode='diameter', opacity=0.5), hoverinfo = 'text',
          text = ~paste( 
            paste("Age: ",age_at_initial_pathologic_diagnosis , sep="") ,
            paste("Pathologic Stage: ", ajcc_pathologic_stage, sep=""),
            paste("Vital Status: ", vital_status,sep=""),
            sep = "<br>")) %>%
  layout(title = "Ages, Pathologic Stage, and Vital Status for Female Breast Cancer Patients",
         xaxis = list(title = "age at diagnosis"), 
         yaxis = list(title = "Pathologic Stage"),
         hovermode = "compare") 


# extract data in stage I 
data_stageI <- sum_exp_dataframe_age %>%
  filter( ajcc_pathologic_stage == "Stage I" | ajcc_pathologic_stage == "Stage IA" | ajcc_pathologic_stage == "Stage IB")

# We initialize a 'survival' object first, which contains the data we need.
 surv_object_data_stageI <- Surv(time = data_stageI$days_to_death, 
                                event = data_stageI$death_event)

# We then create a fit object to do the gender analysis 
age_fit_data_stageI <- surv_fit( surv_object_data_stageI ~ data_stageI$age_category, data = data_stageI )

#the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot_data_stageI = ggsurvplot(age_fit_data_stageI, 
                                  legend = "right",
                                  legend.title = c("Age"),
                                  legend.labs = c("old (N=186)", "young (N=210)"),
                                  title="Survival Analysis for Female Patients at Stage I (N=396)",
                                  ggtheme = theme_bw()) 
 
 plotly::ggplotly(survplot_data_stageI[[1]])

 
 # extract data in stage II 
 data_stageII <- sum_exp_dataframe_age %>%
   filter( ajcc_pathologic_stage == "Stage II" | ajcc_pathologic_stage == "Stage IIA" | ajcc_pathologic_stage == "Stage IIB")
 
 # We initialize a 'survival' object first, which contains the data we need.
 surv_object_data_stageII <- Surv(time = data_stageII$days_to_death, 
                                 event = data_stageII$death_event)
 
 # We then create a fit object to do the gender analysis 
 age_fit_data_stageII <- surv_fit( surv_object_data_stageII ~ data_stageII$age_category, data = data_stageII )
 
 #the ggtheme and legend arguments are for formatting. 
 # Feel free to play around with the margins and legend placement
 survplot_data_stageII = ggsurvplot(age_fit_data_stageII, 
                                   legend = "right",
                                   legend.title = c("Age"),
                                   legend.labs = c("old (N=544)", "young (N=758)"),
                                   title="Survival Analysis for Female Patients at Stage II (N=1302)",
                                   ggtheme = theme_bw()) 
 
 plotly::ggplotly(survplot_data_stageII[[1]])
 
 # extract data in stage III 
 data_stageIII <- sum_exp_dataframe_age %>%
   filter( ajcc_pathologic_stage == "Stage III" | ajcc_pathologic_stage == "Stage IIIA" | ajcc_pathologic_stage == "Stage IIIB" | ajcc_pathologic_stage == "Stage IIIC")
 
 # We initialize a 'survival' object first, which contains the data we need.
 surv_object_data_stageIII <- Surv(time = data_stageIII$days_to_death, 
                                  event = data_stageIII$death_event)
 
 # We then create a fit object to do the gender analysis 
 age_fit_data_stageIII <- surv_fit( surv_object_data_stageIII ~ data_stageIII$age_category, data = data_stageIII )
 
 #the ggtheme and legend arguments are for formatting. 
 # Feel free to play around with the margins and legend placement
 survplot_data_stageIII = ggsurvplot(age_fit_data_stageIII, 
                                    legend = "right",
                                    legend.title = c("Age"),
                                    legend.labs = c("old (N=216)", "young (N=314)"),
                                    title="Survival Analysis for Female Patients at Stage III (N=530)",
                                    ggtheme = theme_bw()) 
 
 plotly::ggplotly(survplot_data_stageIII[[1]])
 
 # extract data in stage IV 
 data_stageIV <- sum_exp_dataframe_age %>%
   filter( ajcc_pathologic_stage == "Stage IV")
 
 # We initialize a 'survival' object first, which contains the data we need.
 surv_object_data_stageIV <- Surv(time = data_stageIV$days_to_death, 
                                   event = data_stageIV$death_event)
 
 # We then create a fit object to do the gender analysis 
 age_fit_data_stageIV <- surv_fit( surv_object_data_stageIV ~ data_stageIV$age_category, data = data_stageIV )
 
 #the ggtheme and legend arguments are for formatting. 
 # Feel free to play around with the margins and legend placement
 survplot_data_stageIV = ggsurvplot(age_fit_data_stageIV, 
                                     legend = "right",
                                     legend.title = c("Age"),
                                     legend.labs = c("old (N=22)", "young (N=18)"),
                                     title="Survival Analysis for Female Patients at Stage IV (N=40)",
                                     ggtheme = theme_bw()) 
 
 plotly::ggplotly(survplot_data_stageIV[[1]])
 
 