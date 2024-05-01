library(reshape2)
library(dplyr)
library(plyr)
library(nhanesR)
library(mice)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##Load the primary data derived from NHANES
load("svrvival_data.RData")

#Process the primary data and calculate the EFI
svrvival_data1 <- svrvival_data
svrvival_data1$MCHC <- svrvival_data1$MCHC/10
svrvival_data1 <- mutate(svrvival_data1,EFI=RDW/MCHC)

#Filter the DM individuals
DM_svrvival_data <- svrvival_data1 |> filter(DM %in%  c("DM","no"))

#Delete the data with zero weight
DM_svrvival_data <- DM_svrvival_data[(DM_svrvival_data$nhs_wt!=0),]

#Delete the individuals with NA value in DM, MCHC, RDW, ASCVD
DM_svrvival_data <- DM_svrvival_data[!is.na(DM_svrvival_data$DM),]
DM_svrvival_data <- DM_svrvival_data[!is.na(DM_svrvival_data$MCHC),]
DM_svrvival_data <- DM_svrvival_data[!is.na(DM_svrvival_data$RDW),]
DM_svrvival_data <- DM_svrvival_data[!is.na(DM_svrvival_data$ASCVD),]

#Delete the zero data in EFI, RDW, MCHC
DM_svrvival_data <- DM_svrvival_data |> filter(EFI != 0)
DM_svrvival_data <- DM_svrvival_data |> filter(RDW != 0)
DM_svrvival_data <- DM_svrvival_data |> filter(MCHC != 0)

#Built the function to delete the NA in lab test
na_delet_data <- function(data,number){
    for(i in colnames(data)){
        na_count_tem <- sum(is.na(data[,i]))
        if(na_count_tem >= number){
            data <- data[,which(colnames(data)!= i)]
        }
    }
    return(data)
}

#Remove clinical features that are missing more than 20%
DM_svrvival_data <- na_delet_data(DM_svrvival_data,nrow(DM_svrvival_data)*0.2)

#Recode the factor
DM_svrvival_data$sex <- Recode(DM_svrvival_data$sex,
                           "Male::1", 
                           "Female::0",
                           to.numeric = F)

# Recode(DM_svrvival_data$CKD)
DM_svrvival_data$CKD <- Recode(DM_svrvival_data$CKD,
	"no::0", 
	"yes::1", 
	"NA::NA",
	to.numeric = F)

# Recode(DM_svrvival_data$COPD)
DM_svrvival_data$COPD <- Recode(DM_svrvival_data$COPD,
	"no::0", 
	"yes::1", 
	"NA::NA",
	to.numeric = F)

# Recode(DM_svrvival_data$Hypertension)
DM_svrvival_data$Hypertension <- Recode(DM_svrvival_data$Hypertension,
	"no::0", 
	"yes::1", 
	"NA::NA",
	to.numeric = F)

# Recode(DM_svrvival_data$Hyperlipidemia)
DM_svrvival_data$Hyperlipidemia <- Recode(DM_svrvival_data$Hyperlipidemia,
	"yes::1", 
	"no::0", 
	"NA::NA",
	to.numeric = F)

# Recode(DM_svrvival_data$Anti_Diabetic)
DM_svrvival_data$Anti_Diabetic <- Recode(DM_svrvival_data$Anti_Diabetic,
	"other::2", 
	"no::0", 
	"yes::1", 
	"NA::NA",
	to.numeric = F)

# Recode(DM_svrvival_data$Anti_Hyperten)     
DM_svrvival_data$Anti_Hyperten <- Recode(DM_svrvival_data$Anti_Hyperten,
	"yes::1", 
	"no::0", 
	"other::2", 
	"NA::NA",
	to.numeric = F)

# Recode(DM_svrvival_data$Anti_Hyperlip)
DM_svrvival_data$Anti_Hyperlip <- Recode(DM_svrvival_data$Anti_Hyperlip,
	"other::2", 
	"no::0", 
	"yes::1", 
	"NA::NA",
	to.numeric = F)

#mice for supplementation of missing values
mice_tem_data1 <- mice(DM_svrvival_data[,-5],m=5,seed = 3,method = "norm")
DM_mice_data <- complete(mice_tem_data1,action = 1)

#Calculate the quartiles of EFI and delete the NA in quartiles of EFI.
DM_mice_data$EFIQ <- quant(DM_mice_data$EFI, n = 4,Q = TRUE,round=2)
DM_mice_data <- DM_mice_data[!is.na(DM_mice_data$EFIQ),]

#Count the NA value
ckkk <- na_count(DM_mice_data,2) |> as.data.frame()

DM_mice_data[is.na(DM_mice_data)] <- 0

DM_mice_data$sex <- Recode(DM_mice_data$sex,
	"Male::1", 
	"Female::0",
	to.numeric = T)

save(DM_mice_data,file="DM_mice_data.RData")
