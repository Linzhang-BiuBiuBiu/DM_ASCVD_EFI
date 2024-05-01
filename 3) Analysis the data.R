library(nhanesR)
library(reshape2)
library(dplyr)
library(survey)
library(openxlsx)
library(plyr)
library(dplyr)
library(rms)

#Load the cleaned data
load("DM_mice_data.RData")

#Summary the Continuity and subtype characteristics
cla_va <- c("CKD", "COPD", "Hypertension", "Hyperlipidemia", "sex", "eth")

gre_va <- c("EFI","RDW","MCHC","age","WBC", "GLU","LymP", "SegneP", "EoP",  
            "lym", "Mon", "SeNeP", "Plt", "MPV",  "ALT", "AST", "BUN", 
            "CA", "TC", "HCO3", "GGT", "TP", "TG",  "Na","CL")


#Calculate the quartile of EFI
DM_mice_data$EFIQ <- quant(DM_mice_data$EFI, n = 4,Q = TRUE,round=2)

#Analysis of differences between different groups
nhs <- svy_design(DM_mice_data)
differ_variable <- svy_tableone(design = nhs,
                                c_meanPMse =T,
                                cv=gre_va,
                                gv=c("EFIQ","DM",cla_va),
                                by = "ASCVD",
                                xlsx = "ASCVD_difference.xlsx")

differ_with_EFIQ <- svy_tableone(design = nhs,
                               c_meanPMse=T,
                               cv=gre_va,
                               gv=c("ASCVD","DM",cla_va),
                               by = "EFIQ",
                               xlsx = "EFIQ_difference.xlsx")

differ_with_EFIQ_DM <- svy_tableone(design = nhs,
                                 c_meanPMse=T,
                                 cv=gre_va,
                                 gv=c("ASCVD",cla_va),
                                 by = "DM",
                                 xlsx = "DM_difference.xlsx")

DM_mice_data_ASCVD <- DM_mice_data |> filter(DM=="DM")
DM_mice_data_ASCVD_svy <- DM_mice_data_ASCVD |> svy_design()
svy_tableone(design = DM_mice_data_ASCVD_svy,
                                 c_meanPMse=T,
                                 cv=gre_va,
                                 gv=c(cla_va,"EFIQ"),
                                 by = "ASCVD",
                                 xlsx = "DM_ASCVD_difference.xlsx")


#Recode the variables
DM_mice_data$DM <- Recode(DM_mice_data$DM,
	"no::0", 
	"DM::1",
	to.numeric = T)

# Recode(DM_mice_data$ASCVD)
DM_mice_data$ASCVD <- Recode(DM_mice_data$ASCVD,
	"no::0", 
	"yes::1",
    	to.numeric = T)

#Single factor logistic regression for DM
nhs <- svy_design(DM_mice_data)
single_logi <- svy_uv.logit(design = nhs,y="DM",
                            x=c(gre_va,"EFIQ",cla_va))

#The adjusted variable in Model4
single_logi1 <- c( "age", "sex", "eth","WBC", "LymP", "SegneP", "EoP",  
                   "lym", "Mon", "SeNeP", "Plt", "MPV",  "ALT", "AST", "BUN", 
                   "CA", "TC", "HCO3", "GGT", "TP", "TG",  "Na", "GLU",
                   "CL")

#The variable for various group
variab <- c("RDW","EFI","MCHC","EFIQ")

# The adjustment of the three models in DM
adj_log_result <- list()
for(target in variab){
    m1 <- c("age","sex")
    m2 <- unique(c(m1, "eth", "CKD", "COPD", "Hypertension", 
                   "Hyperlipidemia"))
    
    m3 <- drop_vector(unique(c(m2,single_logi1)),variab)
    
    log_adi_m1 <- formula(paste0("DM~",paste(c(target,m1),collapse = "+")))
    log_adi_m2 <- formula(paste0("DM~",paste(c(target,m2),collapse = "+")))
    log_adi_m3 <- formula(paste0("DM~",paste(c(target,m3),collapse = "+")))
    log_adj_list <- list(log_adi_m1,log_adi_m2,log_adi_m3)
    
    adj_log_tem_result <- lapply(log_adj_list, 
                                 function(x) 
                                     svyglm(x, design=nhs,
                                            family =quasibinomial() ) 
                                 |> reg_table())
    
    adj_log_result[[target]] <- as.data.frame(lapply(adj_log_tem_result, function(x) x[1:5,c(1,7,5)])) 
    
}

#Save the result
adj_log_result[["single_logi"]] <- single_logi
write.xlsx(adj_log_result,file = "logistic_DM_result.xlsx")

#Single factor logistic regression for ASCVD
single_logi <- svy_uv.logit(design = nhs,y="ASCVD",
                            x=c(gre_va,"EFIQ",cla_va))

# The adjustment of the three models in ASCVD
adj_log_result <- list()
for(target in variab){
    m1 <- c("age","sex")
    m2 <- unique(c(m1, "eth", "CKD", "COPD", "Hypertension", 
                   "Hyperlipidemia"))
    
    m3 <- drop_vector(unique(c(m2,single_logi1)),variab)
    
    log_adi_m1 <- formula(paste0("ASCVD~",paste(c(target,m1),collapse = "+")))
    log_adi_m2 <- formula(paste0("ASCVD~",paste(c(target,m2),collapse = "+")))
    log_adi_m3 <- formula(paste0("ASCVD~",paste(c(target,m3),collapse = "+")))
    log_adj_list <- list(log_adi_m1,log_adi_m2,log_adi_m3)
    
    adj_log_tem_result <- lapply(log_adj_list, 
                                 function(x) 
                                     svyglm(x, design=nhs,
                                            family =quasibinomial() ) 
                                 |> reg_table())
    
    adj_log_result[[target]] <- as.data.frame(lapply(adj_log_tem_result, function(x) x[1:5,c(1,7,5)])) 
    
}
#Save the result
adj_log_result[["single_logi"]] <- single_logi
write.xlsx(adj_log_result,file = "logistic_ASCVD_result.xlsx")
write.xlsx(single_logi,file = "single_logi.XLSX")


#Single factor logistic regression for ASCVD in DM patients
nhs <- DM_mice_data |> filter(DM==1) |> svy_design()
single_logi <- svy_uv.logit(design = nhs,y="ASCVD",
                            x=c(gre_va,"EFIQ",cla_va))

# The adjustment of the three models for ASCVD in DM patients
adj_log_result <- list()
for(target in variab){
    m1 <- c("age","sex")
    m2 <- unique(c(m1, "eth", "CKD", "COPD", "Hypertension", 
                   "Hyperlipidemia"))
    
    m3 <- drop_vector(unique(c(m2,single_logi1)),variab)
    
    log_adi_m1 <- formula(paste0("ASCVD~",paste(c(target,m1),collapse = "+")))
    log_adi_m2 <- formula(paste0("ASCVD~",paste(c(target,m2),collapse = "+")))
    log_adi_m3 <- formula(paste0("ASCVD~",paste(c(target,m3),collapse = "+")))
    log_adj_list <- list(log_adi_m1,log_adi_m2,log_adi_m3)
    
    adj_log_tem_result <- lapply(log_adj_list, 
                                 function(x) 
                                     svyglm(x, design=nhs,
                                            family =quasibinomial() ) 
                                 |> reg_table())
    
    adj_log_result[[target]] <- as.data.frame(lapply(adj_log_tem_result, function(x) x[1:5,c(1,7,5)])) 
    
}
#save the result
adj_log_result[["single_logi"]] <- single_logi
write.xlsx(adj_log_result,file = "logistic_DM_ASCVD_result.xlsx")


#Readd the follow-up information
follow_data <- db_mort(DM_mice_data[,c(-58:-60)]) |> filter(!is.na(mortstat)) |> filter(ASCVD == "1")
follow_data$EFIQ <- quant(follow_data$EFI, n = 4,Q = TRUE,round=2)

#Analysis the difference between various follow-up group
nhs <- svy_design(follow_data)
differ_variable <- svy_tableone(design = nhs,
                                c_meanPMse =T,
                                cv=c("permth_int",gre_va),
                                gv=c(cla_va),
                                by = "mortstat",
                                xlsx = "mortstat_difference.xlsx")

differ_variable <- svy_tableone(design = nhs,
                                c_meanPMse =T,
                                cv=c("permth_int",gre_va),
                                gv=c("mortstat",cla_va),
                                by = "EFIQ",
                                xlsx = "mortstat_EFIQ_difference.xlsx")


# Recode the Survival situation
follow_data$mortstat <- Recode(follow_data$mortstat,
                                  "Assumed deceased::1", 
                                  "Assumed alive::0",
                                  to.numeric = T)

#Single cox regression
nhs <- svy_design(follow_data)

uni_cox <- svy_uv.cox(nhs,time = "permth_int",
                      status = "mortstat",
                      x=c(gre_va,"EFIQ",drop_vector(cla_va,"mortstat")))

# The adjustment of the three models for ASCVD in DM patients
variab <- c("RDW","EFI","EFIQ","MCHC")
adj_cox_result <- list()

for(target in variab){
    m1 <- c("age","sex")
    m2 <- unique(c(m1, "eth", "CKD", "COPD", "Hypertension", 
                   "Hyperlipidemia"))
    
    m3 <- drop_vector(unique(c(m2,single_logi1)),variab)
    
    cox_adi_m1 <- formula(paste0("Surv(permth_int,mortstat)~",paste(c(target,m1),collapse = "+")))
    cox_adi_m2 <- formula(paste0("Surv(permth_int,mortstat)~",paste(c(target,m2),collapse = "+")))
    cox_adi_m3 <- formula(paste0("Surv(permth_int,mortstat)~",paste(c(target,m3),collapse = "+")))
    
    cox_adi_list <- list(cox_adi_m1,cox_adi_m2,cox_adi_m3)
    
    cox_tem_result <- lapply(cox_adi_list, 
                             function(x) 
                                 svycoxph(x, design=nhs) 
                             |> reg_table())
    
    names_index <- paste0(target,"_")
    
    adj_cox_result[[names_index]] <- as.data.frame(lapply(cox_tem_result, function(x) x[1:5,c(1,8,6)]))
    
}
#Save the result
adj_cox_result[["uni_cox"]] <- uni_cox
write.xlsx(adj_cox_result,file = "adj_cox_result.xlsx")

#RCS for cut-off value 
Hyp_survi <-follow_data

serch_index <- "EFI"
adju_index <-  "age"

#Optimaze the number of knot
S <- Surv(Hyp_survi$permth_int,Hyp_survi$mortstat==1)
for (knot in 3:10) {
    cph_formula <- as.formula(paste0("S~rcs(",serch_index,",",knot,")+",adju_index))
    
    fit <- cph(cph_formula,data=Hyp_survi,x= TRUE, y= TRUE, surv = TRUE)
    tmp <- extractAIC(fit)
    if(knot==3){AIC=tmp[2];nk=3}
    if(tmp[2]<AIC){AIC=tmp[2];nk=knot}
}

#With the optimal knot for RCS
pacman::p_load(rms,survminer,ggplot2,ggsci)
cph_formula <- as.formula(paste0("S~rcs(",serch_index,",",nk,")+",adju_index))

fit.RAR <- cph(cph_formula, x=TRUE, y=TRUE,data=Hyp_survi)
cox.zph(fit.RAR, "rank")    
ggcoxzph(cox.zph(fit.RAR, "rank"))

#Print the plot
pdf_name <- paste0("RCS_",serch_index,"-",adju_index,".pdf")
pdf(pdf_name,width = 5,height = 5)
ggcoxzph(cox.zph(fit.RAR, "rank")) 
dev.off()

#Explore the cut-off value
anova(fit.RAR)                      
p <-round(anova(fit.RAR)[,3],3)
Pre_HR.RAR <-rms::Predict(fit.RAR,EFI,fun=exp,type="predictions",ref.zero=T,conf.int = 0.95,digits=2)
ggplot(Pre_HR.RAR)


# Print the processment
rcs_all <-  ggplot()+
    geom_line(data=Pre_HR.RAR,
              aes(EFI,yhat),
              linetype="solid",
              size=1,
              alpha=0.7,
              colour="purple")+
    scale_color_nejm()+
    geom_ribbon(data=Pre_HR.RAR,
                aes(EFI, ymin=lower,ymax=upper,fill="purple"),alpha=0.1)+
    theme_classic()+
    scale_fill_nejm()+
    geom_hline(yintercept=1,linetype=2,size=0.75)
labs(title ="风险随LnAl变化曲RCS",
     x="EFI", 
     y="HR (95%CI)")

pdf("rcs_EFI.pdf",width = 6,height = 4)
rcs_all
dev.off()

#With the cut-off value to discriminate the hign or low group
Hyp_survi$EFI_4.02 <- ifelse(Hyp_survi$EFI>4.02,"High","Low")
Hyp_survi$RDW_13.61 <- ifelse(Hyp_survi$RDW>13.61,"High","Low")
Hyp_survi$MCHC_3.37 <- ifelse(Hyp_survi$MCHC>3.37,"High","Low")


#KM curve with high or low group
nhs <- svy_design(Hyp_survi)
svy_population(nhs)

pdf(file = "survival3.pdf",width = 9,height = 4.9)
EFI <- svykm(Surv(permth_int,mortstat)~EFI_4.02,design = nhs)
svy_kmplot(EFI,ci = F)

RDW <- svykm(Surv(permth_int,mortstat)~RDW_13.61,design = nhs)
svy_kmplot(RDW,ci = F)

MCHC <- svykm(Surv(permth_int,mortstat)~MCHC_3.37,design = nhs)
svy_kmplot(MCHC,ci = F)

EFIQ <- svykm(Surv(permth_int,mortstat)~EFIQ,design = nhs)
svy_kmplot(EFIQ,ci = F)

dev.off()