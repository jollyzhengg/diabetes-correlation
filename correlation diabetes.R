#Correlation Diabetes Data----
#Author: Jolly Zheng
#Date Created: January 12, 2021
#Data Last Modified: December 10, 2021
#libraries----
library(reshape2)
library(openxlsx)
library(tidyverse)
#library(magrittr)
library(ggplot2)
library(dplyr)
library(pheatmap)
#working directory----
setwd("/Users/jollyzheng/Desktop/r/diabetes")
load("Diabetes_Data2.rda")
#p-vals----
new_clinicals <- Clinicals[c("Body mass index (kg/m2) [CLIN]", "Waist circumference (cm) [CLIN]",
                                 "Hip circumference (cm) [CLIN]", "Waist-to-hip ratio [CLIN]",                 
                                 "Mean systolic blood pressure (mmHg) [CLIN]", "Mean AF score [CLIN]",                      
                                 "HbA1c (%) [CLIN]","Triglycerides (mmol/L) [CLIN]",
                                 "Corrected calcium (mmol/L) [CLIN]","Calcium (mmol/L) [CLIN]",
                                 "Urea (mmol/L) [CLIN]", "Creatinine (umol/L) [CLIN]",
                                 "Alkaline phosphotase (U/L) [CLIN]", "Sodium (mmol/L) [CLIN]",
                                 "Potasium (mmol/L) [CLIN]", "WBC (10*3 /uL) [CLIN]",
                                 "Neutrophils (10*3 /uL) [CLIN]","Basophils (10*3 /uL) [CLIN]"),] 
total <- rbind(new_clinicals,Metabolites, Proteins, Diabetes_diagnosis)
total <- t(total)

pvals<- NULL

for (i in c(1:18)){
  for (j in c(19:2003)){
    pvals2 <- cor.test(total[,i],total[,j],use="pairwise.complete.obs")$p.value
    pvals <- c(pvals,pvals2)
  }
}

variable_one <- rep(c(colnames(total)[1:18]),each=1985)
variable_two <- rep(c(colnames(total)[19:2003]),18)
pval_data <- data.frame(variable_one,variable_two, pvals, padj= p.adjust(pvals, method = 'fdr'))
sig_pval_data <- filter(pval_data,padj <= 0.05)
sig_pval_data <- sig_pval_data[-c(3)]

hb1ac_sig_pval_data <- subset(sig_pval_data,variable_one == "HbA1c (%) [CLIN]") 
hb1ac <- hb1ac_sig_pval_data[,2]
hb1ac <- gsub(pattern = "[SOMA]",replacement = "", hb1ac, fixed = TRUE)
hb1ac <- gsub(pattern = "[HD4]",replacement = "", hb1ac, fixed = TRUE )
hb1ac <- sapply(strsplit(hb1ac," :",fixed = TRUE),function(i){i[1]})
hb1ac_sig_pval_data <- data.frame(hb1ac_sig_pval_data[,1],hb1ac, hb1ac_sig_pval_data[,3])
write.csv(hb1ac_sig_pval_data, file = "pitracerfile/hb1ac.csv")

w2h_sig_pval_data <- subset(sig_pval_data,variable_one == "Waist-to-hip ratio [CLIN]")
w2h <- w2h_sig_pval_data[,2]
w2h <- gsub(pattern = "[SOMA]",replacement = "", w2h, fixed = TRUE)
w2h <- gsub(pattern = "[HD4]",replacement = "", w2h, fixed = TRUE)
w2h <- sapply(strsplit(w2h," :",fixed = TRUE),function(i){i[1]})
w2h_sig_pval_data <- data.frame(w2h_sig_pval_data[,1],w2h, w2h_sig_pval_data[,3])
write.csv(w2h_sig_pval_data, file = "pitracerfile/waist_to_hip_ratio.csv")

pressure_sig_pval_data <- subset(sig_pval_data,variable_one == "Mean systolic blood pressure (mmHg) [CLIN]")
pressure <- pressure_sig_pval_data[,2]
pressure <- gsub(pattern = "[SOMA]",replacement = "", pressure, fixed = TRUE)
pressure <- gsub(pattern = "[HD4]",replacement = "", pressure, fixed = TRUE)
pressure <- sapply(strsplit(pressure," :",fixed = TRUE),function(i){i[1]})
pressure_sig_pval_data <- data.frame(pressure_sig_pval_data[,1],pressure, pressure_sig_pval_data[,3])
write.csv(pressure_sig_pval_data, file = "pitracerfile/blood_pressure.csv")

urea_sig_pval_data <- subset(sig_pval_data,variable_one == "Urea (mmol/L) [CLIN]")
urea <- urea_sig_pval_data[,2]
urea <- gsub(pattern = "[SOMA]",replacement = "", urea, fixed = TRUE)
urea <- gsub(pattern = "[HD4]",replacement = "", urea, fixed = TRUE)
urea <- sapply(strsplit(urea," :",fixed = TRUE),function(i){i[1]})
urea_sig_pval_data <- data.frame(urea_sig_pval_data[,1],urea, urea_sig_pval_data[,3])
write.csv(urea_sig_pval_data, file = "pitracerfile/urea.csv")

sod_sig_pval_data <- subset(sig_pval_data,variable_one == "Sodium (mmol/L) [CLIN]")
sod <- sod_sig_pval_data[,2]
sod <- gsub(pattern = "[SOMA]",replacement = "", sod, fixed = TRUE)
sod <- gsub(pattern = "[HD4]",replacement = "", sod, fixed = TRUE)
sod <- gsub(pattern = "*",replacement = "", sod, fixed = TRUE)
sod <- sapply(strsplit(sod," :",fixed = TRUE),function(i){i[1]})
sod_sig_pval_data <- data.frame(sod_sig_pval_data[,1],sod, sod_sig_pval_data[,3])
write.csv(sod_sig_pval_data, file = "pitracerfile/sodium.csv")

#total_and_pvals_prep for diab
total_v2 <- rbind(Diabetes_diagnosis, new_clinicals, Metabolites, Proteins)
total_v2 <- t(total_v2)

pvalsdiab <- NULL

for (j in c(19:2003)){
  pvals2diab <- cor.test(total_v2[,1],total_v2[,j],use="pairwise.complete.obs")$p.value
  pvalsdiab <- c(pvalsdiab,pvals2diab)
}

variable_one_v2 <- rep(c(colnames(total_v2)[1]),each=1985)
variable_two_v2 <- c(colnames(total_v2)[19:2003])
diab_sig_pval_data <- data.frame(variable_one_v2,variable_two_v2, pvalsdiab, padj= p.adjust(pvalsdiab, method = 'fdr'))
diab_sig_pval_data <- filter(diab_sig_pval_data,padj <= 0.05)
diab <- diab_sig_pval_data[,2]
diab <- gsub(pattern = "[SOMA]",replacement = "", diab, fixed = TRUE)
#diab<- gsub(pattern = "[HD4]",replacement = "", diab, fixed = TRUE)
#diab <- gsub(pattern = "*",replacement = "", diab, fixed = TRUE)
diab <- sapply(strsplit(diab," :",fixed = TRUE),function(i){i[1]})
diab_sig_pval_data <- data.frame(diab_sig_pval_data[,1],diab
                                 , diab_sig_pval_data[,3])

#parttwo----
load("HMDBs_table.rda")
colnames(diab_sig_pval_data)[3]= "diabpadj"
left_join(diab_sig_pval_data, HMDB_table, by = c('diab' = 'name')) %>%
  transmute(diab = coalesce(HMDB, diab),
            diabpadj) -> diab_sig_pval_data
write.csv(diab_sig_pval_data, file = "pitracerfile/diabetes.csv") 

pressure_sig_pval_data <- subset(sig_pval_data,variable_one == "Mean systolic blood pressure (mmHg) [CLIN]")
pressure <- pressure_sig_pval_data[,2]
# pressure <- gsub(pattern = "[SOMA]",replacement = "", pressure, fixed = TRUE)
# pressure <- gsub(pattern = "[HD4]",replacement = "", pressure, fixed = TRUE)
pressure <- sapply(strsplit(pressure," :",fixed = TRUE),function(i){i[1]})
pressure_sig_pval_data <- data.frame(pressure_sig_pval_data[,1],pressure, pressure_sig_pval_data[,3])
colnames(pressure_sig_pval_data)[3]= "pressurepadj"
left_join(pressure_sig_pval_data, HMDB_table, by = c('pressure' = 'name')) %>%
  transmute(pressure = coalesce(HMDB, pressure),
            pressurepadj) -> pressure_sig_pval_data
write.csv(pressure_sig_pval_data, file = "pitracerfile/blood_pressure.csv")

#tracer analysis----
#files
genes_blood_pressure <- read.csv("tracer_results/genes_blood_pressure.csv")
genes_diabetes<- read.csv("tracer_results/genes_diabetes.csv")
genes_hb1ac <- read.csv("tracer_results/genes_hb1ac.csv")
genes_sodium <- read.csv("tracer_results/genes_sodium.csv")
genes_urea <- read.csv("tracer_results/genes_urea.csv")
genes_waist_to_hip <- read.csv("tracer_results/genes_waist_to_hip.csv")

metabolites_blood_pressure <- read.csv("tracer_results/metabolites_blood_pressure.csv")
metabolites_diabetes<- read.csv("tracer_results/metabolites_diabetes.csv")
metabolites_hb1ac <- read.csv("tracer_results/metabolites_hb1ac.csv")
metabolites_sodium <- read.csv("tracer_results/metabolites_sodium.csv")
metabolites_urea <- read.csv("tracer_results/metabolites_urea.csv")
metabolites_waist_to_hip <- read.csv("tracer_results/metabolites_waist_to_hip.csv")

#prep_dataframe
prep_dataframe <- function(x){
  x[,8] <- x[,3]/x[,5]
  x <- subset(x, x[,7] <= 0.05)
#  x <- subset(x,x[,8] >= median(x[,8]))
}

genes_blood_pressure <- prep_dataframe(genes_blood_pressure)
genes_diabetes <- prep_dataframe(genes_diabetes)
genes_hb1ac <- prep_dataframe(genes_hb1ac)
genes_sodium <- prep_dataframe(genes_sodium)
genes_urea <- prep_dataframe(genes_urea)
genes_waist_to_hip <- prep_dataframe(genes_waist_to_hip)

metabolites_blood_pressure <- prep_dataframe(metabolites_blood_pressure)
metabolites_diabetes <- prep_dataframe(metabolites_diabetes)
metabolites_hb1ac <- prep_dataframe(metabolites_hb1ac)
metabolites_sodium <- prep_dataframe(metabolites_sodium)
metabolites_urea <- prep_dataframe(metabolites_urea)
metabolites_waist_to_hip <- prep_dataframe(metabolites_waist_to_hip)

#vs_dataframe
vs_dataframe_creation <- function(x,y){
  subset(merge(x,y, by.x = colnames(x)[1], colnames(y)[1]), 
         merge(x,y, by.x = colnames(x)[1], colnames(y)[1])[,8] >= median(merge(x,y, by.x = colnames(x)[1], colnames(y)[1])[,8]))
}

#genes
genes_diabetes_vs_blood_pressure <- vs_dataframe_creation(genes_diabetes,genes_blood_pressure)
save(genes_diabetes_vs_blood_pressure,file= "pathway_analysis/genes_diabetes_vs_blood_pressure.rda")
genes_diabetes_vs_hb1ac <- vs_dataframe_creation(genes_diabetes,genes_hb1ac)
save(genes_diabetes_vs_hb1ac,file= "pathway_analysis/genes_diabetes_vs_hb1ac.rda")
genes_diabetes_vs_sodium <- vs_dataframe_creation(genes_diabetes,genes_sodium)
save(genes_diabetes_vs_sodium,file= "pathway_analysis/genes_diabetes_vs_sodium.rda")
genes_diabetes_vs_urea <- vs_dataframe_creation(genes_diabetes,genes_urea)
save(genes_diabetes_vs_urea,file= "pathway_analysis/genes_diabetes_vs_urea.rda")
genes_diabetes_vs_waist_to_hip <- vs_dataframe_creation(genes_diabetes,genes_waist_to_hip)
save(genes_diabetes_vs_waist_to_hip,file= "pathway_analysis/genes_diabetes_vs_waist_to_hip.rda")

#metabolites
metabolites_diabetes_vs_blood_pressure <- vs_dataframe_creation(metabolites_diabetes, metabolites_blood_pressure)
save(metabolites_diabetes_vs_blood_pressure,file= "pathway_analysis/metabolites_diabetes_vs_blood_pressure.rda")
metabolites_diabetes_vs_hb1ac <- vs_dataframe_creation(metabolites_diabetes, metabolites_hb1ac)
save(metabolites_diabetes_vs_hb1ac,file= "pathway_analysis/metabolites_diabetes_vs_hb1ac.rda")
metabolites_diabetes_vs_sodium <- vs_dataframe_creation(metabolites_diabetes, metabolites_sodium)
save(metabolites_diabetes_vs_sodium,file= "pathway_analysis/metabolites_diabetes_vs_sodium.rda")
metabolites_diabetes_vs_urea <- vs_dataframe_creation(metabolites_diabetes, metabolites_urea)
save(metabolites_diabetes_vs_urea,file= "pathway_analysis/metabolites_diabetes_vs_urea.rda")
metabolites_diabetes_vs_waist_to_hip <- vs_dataframe_creation(metabolites_diabetes, metabolites_waist_to_hip)
save(metabolites_diabetes_vs_waist_to_hip,file= "pathway_analysis/metabolites_diabetes_vs_waist_to_hip.rda")

#rest of variables----
#hmdb setup
hmdb_setup <- function(x,y,z,aa,bb){
  x <- subset(sig_pval_data, variable_one == y)
  z <- x[,2]
  z <- sapply(strsplit(z," :",fixed = TRUE),function(i){i[1]})
  x <- data.frame(x[,1], z, x[,3])
  colnames(x)[3]= aa
  left_join(x, HMDB_table, by = c(z = 'name')) %>%
    transmute(z = coalesce(HMDB, z),
              colnames(x)[3]) -> x
  z <- gsub(pattern = "[SOMA]",replacement = "", z, fixed = TRUE)
  z <- gsub(pattern = "[HD4]",replacement = "", z, fixed = TRUE)
  write.csv(x, file = bb)
}

hmdb_setup(bmi_sig_pval_data, "Body mass index (kg/m2) [CLIN]",
           bmi, "bmiadj", "pitracerfile/bmi.csv")
hmdb_setup(wcircum_sig_pval_data, "Waist circumference (cm) [CLIN]",
           wcircum, "wcircum", "pitracerfile/wcircum.csv")
hmdb_setup(hcircum_sig_pval_data, "Hip circumference (cm) [CLIN]",
           hcircum, "hcircumadj", "pitracerfile/hcircum.csv")
hmdb_setup(meanaf_sig_pval_data, "Mean AF score [CLIN]",
           meanaf, "meanafadj", "pitracerfile/meanaf.csv")
hmdb_setup(triglyc_sig_pval_data, "Triglycerides (mmol/L) [CLIN]",
           triglyc, "triglycadj", "pitracerfile/triglyc.csv")
hmdb_setup(corrcalcium_sig_pval_data, "Corrected calcium (mmol/L) [CLIN]",
           corrcalcium, "corrcalciumadj", "pitracerfile/corrcalcium.csv")
hmdb_setup(calcium_sig_pval_data, "Calcium (mmol/L) [CLIN]",
           calcium, "calciumadj", "pitracerfile/calcium.csv")
hmdb_setup(creatinine_sig_pval_data, "Creatinine (umol/L) [CLIN]",
           creatinine, "creatinineadj", "pitracerfile/creatinine.csv")
hmdb_setup(alkaline_sig_pval_data, "Alkaline phosphotase (U/L) [CLIN]",
           alkaline, "alkalineadj", "pitracerfile/alkaline.csv")
hmdb_setup(potassium_sig_pval_data, "Potasium (mmol/L) [CLIN]",
           potassium, "potassiumadj", "pitracerfile/potassium.csv")
hmdb_setup(wbc_sig_pval_data, "WBC (10*3 /uL) [CLIN]",
           wbc, "wbcadj", "pitracerfile/wbc.csv")
hmdb_setup(neutrophil_sig_pval_data, "Neutrophils (10*3 /uL) [CLIN]",
           neutrophil, "neutrophiladj", "pitracerfile/neutrophil.csv")
hmdb_setup(basophil_sig_pval_data, "Basophils (10*3 /uL) [CLIN]",
           basophil, "basophiladj", "pitracerfile/basophil.csv")

prep_dataframe <- function(x){
  x[,8] <- x[,3]/x[,5]
  x <- subset(x, x[,7] <= 0.05)
  #  x <- subset(x,x[,8] >= median(x[,8]))
}
vs_dataframe_creation <- function(x,y){
  subset(merge(x,y, by.x = colnames(x)[1], colnames(y)[1]), 
         merge(x,y, by.x = colnames(x)[1], colnames(y)[1])[,8] >= median(merge(x,y, by.x = colnames(x)[1], colnames(y)[1])[,8]))
}

read_prep_dataframe <- function(x,y){
  x <- read.csv(y)
  x[,8] <- x[,3]/x[,5]
  x <- subset(x, x[,7] <= 0.05)
  #  x <- subset(x,x[,8] >= median(x[,8]))
}

#metabolites...
m_vs_dataframe_creation <- function(x,y,z){
  y <- subset(merge(metabolites_diabetes, x, by.x = colnames(metabolites_diabetes)[1], colnames(x)[1]), 
         merge(metabolites_diabetes, x, by.x = colnames(metabolites_diabetes)[1], colnames(x)[1])[,8]
         >= median(merge(metabolites_diabetes, x, by.x = colnames(metabolites_diabetes)[1], colnames(x)[1])[,8]))
  save(y, file= z)
}

metabolites_bmi <- read_prep_dataframe(metabolites_bmi,"tracer_results/metabolites_bmi.csv")
m_vs_dataframe_creation(metabolites_bmi, metabolites_diabetes_vs_bmi, "pathway_analysis/metabolites_diabetes_vs_bmi.rda")

metabolites_wcircum <- read_prep_dataframe(metabolites_wcircum,"tracer_results/metabolites_wcircum.csv")
m_vs_dataframe_creation(metabolites_wcircum, metabolites_diabetes_vs_wcircum, "pathway_analysis/metabolites_diabetes_vs_wcircum.rda")

metabolites_hcircum <- read_prep_dataframe(metabolites_hcircum,"tracer_results/metabolites_hcircum.csv")
m_vs_dataframe_creation(metabolites_hcircum, metabolites_diabetes_vs_hcircum, "pathway_analysis/metabolites_diabetes_vs_hcircum.rda")

metabolites_meanaf <- read_prep_dataframe(metabolites_meanaf,"tracer_results/metabolites_meanaf.csv")
m_vs_dataframe_creation(metabolites_meanaf, metabolites_diabetes_vs_meanaf, "pathway_analysis/metabolites_diabetes_vs_meanaf.rda")

metabolites_triglyc <- read_prep_dataframe(metabolites_triglyc,"tracer_results/metabolites_triglyc.csv")
m_vs_dataframe_creation(metabolites_triglyc, metabolites_diabetes_vs_triglyc, "pathway_analysis/metabolites_diabetes_vs_triglyc.rda")

metabolites_corrcalcium <- read_prep_dataframe(metabolites_corrcalcium,"tracer_results/metabolites_corrcalcium.csv")
m_vs_dataframe_creation(metabolites_corrcalcium, metabolites_diabetes_vs_corrcalcium, "pathway_analysis/metabolites_diabetes_vs_corrcalcium.rda")

metabolites_calcium <- read_prep_dataframe(metabolites_calcium,"tracer_results/metabolites_calcium.csv")
m_vs_dataframe_creation(metabolites_calcium, metabolites_diabetes_vs_calcium, "pathway_analysis/metabolites_diabetes_vs_calcium.rda")

metabolites_creatinine <- read_prep_dataframe(metabolites_creatinine,"tracer_results/metabolites_creatinine.csv")
m_vs_dataframe_creation(metabolites_creatinine, metabolites_diabetes_vs_creatinine, "pathway_analysis/metabolites_diabetes_vs_creatinine.rda")

metabolites_alkaline <- read_prep_dataframe(metabolites_alkaline,"tracer_results/metabolites_alkaline.csv")
m_vs_dataframe_creation(metabolites_alkaline, metabolites_diabetes_vs_alkaline, "pathway_analysis/metabolites_diabetes_vs_alkaline.rda")

metabolites_potassium <- read_prep_dataframe(metabolites_potassium,"tracer_results/metabolites_potassium.csv")
m_vs_dataframe_creation(metabolites_potassium, metabolites_diabetes_vs_potassium, "pathway_analysis/metabolites_diabetes_vs_potassium.rda")

metabolites_wbc <- read_prep_dataframe(metabolites_wbc,"tracer_results/metabolites_wbc.csv")
m_vs_dataframe_creation(metabolites_wbc, metabolites_diabetes_vs_wbc, "pathway_analysis/metabolites_diabetes_vs_wbc.rda")

metabolites_neutrophil <- read_prep_dataframe(metabolites_neutrophil,"tracer_results/metabolites_neutrophil.csv")
m_vs_dataframe_creation(metabolites_neutrophil, metabolites_diabetes_vs_neutrophil, "pathway_analysis/metabolites_diabetes_vs_neutrophil.rda")

metabolites_basophil <- read_prep_dataframe(metabolites_basophil,"tracer_results/metabolites_basophil.csv")
m_vs_dataframe_creation(metabolites_basophil,metabolites_diabetes_vs_basophil, "pathway_analysis/metabolites_diabetes_vs_basophil.rda")

#genes....

g_vs_dataframe_creation <- function(x,y,z){
  y <- subset(merge(genes_diabetes, x, by.x = colnames(genes_diabetes)[1], colnames(x)[1]), 
         merge(genes_diabetes, x, by.x = colnames(genes_diabetes)[1], colnames(x)[1])[,8]
         >= median(merge(genes_diabetes, x, by.x = colnames(genes_diabetes)[1], colnames(x)[1])[,8]))
  save(y, file= z)
}

genes_bmi <- read_prep_dataframe(genes_bmi,"tracer_results/genes_bmi.csv")
g_vs_dataframe_creation(genes_bmi, genes_diabetes_vs_bmi, "pathway_analysis/genes_diabetes_vs_bmi.rda")

genes_wcircum <- read_prep_dataframe(genes_wcircum,"tracer_results/genes_wcircum.csv")
g_vs_dataframe_creation(genes_wcircum, genes_diabetes_vs_wcircum, "pathway_analysis/genes_diabetes_vs_wcircum.rda")

genes_hcircum <- read_prep_dataframe(genes_hcircum,"tracer_results/genes_hcircum.csv")
g_vs_dataframe_creation(genes_hcircum, genes_diabetes_vs_hcircum, "pathway_analysis/genes_diabetes_vs_hcircum.rda")

genes_meanaf <- read_prep_dataframe(genes_meanaf,"tracer_results/genes_meanaf.csv")
g_vs_dataframe_creation(genes_meanaf, genes_diabetes_vs_meanaf, "pathway_analysis/genes_diabetes_vs_meanaf.rda")

genes_triglyc <- read_prep_dataframe(genes_triglyc,"tracer_results/genes_triglyc.csv")
g_vs_dataframe_creation(genes_triglyc, genes_diabetes_vs_triglyc, "pathway_analysis/genes_diabetes_vs_triglyc.rda")

genes_corrcalcium <- read_prep_dataframe(genes_corrcalcium,"tracer_results/genes_corrcalcium.csv")
g_vs_dataframe_creation(genes_corrcalcium, genes_diabetes_vs_corrcalcium, "pathway_analysis/genes_diabetes_vs_corrcalcium.rda")

genes_calcium <- read_prep_dataframe(genes_calcium,"tracer_results/genes_calcium.csv")
g_vs_dataframe_creation(genes_calcium, genes_diabetes_vs_calcium, "pathway_analysis/genes_diabetes_vs_calcium.rda")

genes_creatinine <- read_prep_dataframe(genes_creatinine,"tracer_results/genes_creatinine.csv")
g_vs_dataframe_creation(genes_creatinine, genes_diabetes_vs_creatinine, "pathway_analysis/genes_diabetes_vs_creatinine.rda")

genes_alkaline <- read_prep_dataframe(genes_alkaline,"tracer_results/genes_alkaline.csv")
g_vs_dataframe_creation(genes_alkaline, genes_diabetes_vs_alkaline, "pathway_analysis/genes_diabetes_vs_alkaline.rda")

genes_potassium <- read_prep_dataframe(genes_potassium,"tracer_results/genes_potassium.csv")
g_vs_dataframe_creation(genes_potassium, genes_diabetes_vs_potassium, "pathway_analysis/genes_diabetes_vs_potassium.rda")

genes_wbc <- read_prep_dataframe(genes_wbc,"tracer_results/genes_wbc.csv")
g_vs_dataframe_creation(genes_wbc, genes_diabetes_vs_wbc, "pathway_analysis/genes_diabetes_vs_wbc.rda")

genes_neutrophil <- read_prep_dataframe(genes_neutrophil,"tracer_results/genes_neutrophil.csv")
g_vs_dataframe_creation(genes_neutrophil, genes_diabetes_vs_neutrophil, "pathway_analysis/genes_diabetes_vs_neutrophil.rda")

genes_basophil <- read_prep_dataframe(genes_basophil,"tracer_results/genes_basophil.csv")
g_vs_dataframe_creation(genes_basophil, genes_diabetes_vs_basophil, "pathway_analysis/genes_diabetes_vs_basophil.rda")


