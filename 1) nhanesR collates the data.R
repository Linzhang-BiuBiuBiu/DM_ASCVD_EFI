library(nhanesR)
library(reshape2)
library(dplyr)
library(plyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##Acquire the diagnosis information
ASCVD <- diag_ASCVD()
DM <- diag_DM()

CKD <- diag_CKD()
COPD <- diag_COPD()
Hype <- diag_Hypertension()
Hyli <- diag_Hyperlipidemia()

##Merge and rebuild the diagnosis information
all_disease <- left_join(ASCVD,DM)
all_disease <- left_join(all_disease,CKD)
all_disease <- left_join(all_disease,COPD)
all_disease <- left_join(all_disease,CKD)
all_disease <- left_join(all_disease,COPD)
all_disease <- left_join(all_disease,Hype)
all_disease <- left_join(all_disease,Hyli)

#Acquire the age, sex, and ethnic details. 
d_demo <- db_demo(all_disease,
                  ageyr = "age",
                  sex="sex",
                  eth1="eth",
                  wtmec2yr = T,
                  wtmec4yr = T,
                  Year = T)

#Acquire the information of blood routine.
d_cbc <- db_cbc(d_demo,
                wbc_1000cells.ul = "WBC",
                Lymphocyte_percent="LymP",
                Monocyte_percent="MonP",
                Segmented_neutrophils_percent="SegneP",
                Eosinophils_percent="EoP",
                Basophils_percent="BaP",
                lymphocyte_number_1000cells.ul="lym",
                Monocyte_number_1000cells.ul="Mon",
                Segmented_neutrophils_number_1000cells.ul="SeNeP",
                Eosinophils_number_1000cells.ul="Eo",
                Basophils_number_1000cells.ul="Ba",
                Red_blood_cell_count_MillionCells.uL="RBC",
                hemoglobin_g.dl="Hg",
                hematocrit="Hem",
                Red_cell_distribution_width="RDW",
                Platelet_count_1000cells.uL="Plt",
                Mean_platelet_volume_fL="MPV",
                Mean_cell_hemoglobin_pg="MCH",
                Mean_cell_hemoglobin_concentration_g.dL="MCHC",
                Mean_cell_volume_fL="MCV")

#Add the biochemical test files
nhs.tsv <- nhs_tsv("lab18|l40_b|l40_c|biopro_d|biopro_e|biopro_f|biopro_g|biopro_h|biopro_i|biopro_j|p_biopro|ssbnp_a")

##Add the  biochemical test information. 
nhs_bio <- nhs_read(nhs.tsv,"LBXSAL:ALB","LBXSATSI:ALT","LBXSASSI:AST",
                    "LBDSAPSI:ALP","LBXSBU:BUN","LBXSCA:CA","LBXSCH:TC","LBXSC3SI:HCO3",
                    "LBXSGTSI:GGT","LBXSGL:GLU","LBXSIR:Fe","LBDSLDSI:LDH",
                    "LBDSPH:P","LBDSTB:TBIL","LBXSTP:TP","LBXSTR:TG","LBXSUA:UA",
                    "LBDSCR:SCR","LBXSNASI:Na","LBXSCLSI:CL","LBXSGB:glb","LBXFSH:FSH",
                    "LBDLHSI:LH","SSBNP:ntProBNP")

dis_lab_data <- full_join(d_cbc,nhs_bio)

#Add CRE and uACR information.
svrvival_data <- db_urine.alb.cr(dis_lab_data,
                                 creatinine_urine_mg.dl="CRE",
                                 uACR_mg.g="uACR")

#Add drug use information
svrvival_data <- drug_anti.Diabetic(svrvival_data)

colnames(svrvival_data) <- gsub("take_drug","Anti_Diabetic",colnames(svrvival_data))

svrvival_data <- drug_anti.Hypertensive(svrvival_data)
colnames(svrvival_data) <- gsub("take_drug","Anti_Hyperten",colnames(svrvival_data))

svrvival_data <- drug_anti.Hyperlipidemic(svrvival_data)
colnames(svrvival_data) <- gsub("take_drug","Anti_Hyperlip",colnames(svrvival_data))

#Calculate the weight details
n <- length(unique(svrvival_data$Year))
svrvival_data$nhs_wt <- ifelse(svrvival_data$Year %in% c("1999-2000","2001-2002"),
                               2/n*svrvival_data$wtmec4yr,1/n*svrvival_data$wtmec2yr)

#Add the follow-up time details
svrvival_data <- db_mort(svrvival_data)

#Save the merged data
save(svrvival_data,file="svrvival_data.RData")
