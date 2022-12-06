How to access to the data:

There are two clinic data that need to be accessed:

1. the clinical.csv which can be directly accessed through read.csv function

# read clinic data
clinic_1 = read_delim(file="GDCdata/clinical.csv", delim = "\t")

2. Get the TCGA data by accession code “BRCA.”

# downloading clinical data
clin_query <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Clinical",
                       file.type = "xml")

# GDCdownload(clin_query)
# get the clinic dataframe
clinic_2 <- GDCprepare_clinic(clin_query, clinical.info = "patient")

3. Merge the pathologic stage column from clinic_2 to clinic_1 by the same id "case_submitter_id"

colnames(clinic_2)[1] <- "case_submitter_id"
clinic <- merge(x = clinic_2, y = clinic_1[ , c("case_submitter_id", "ajcc_pathologic_stage")], by = "case_submitter_id", all.x=TRUE)

