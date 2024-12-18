# # Logistic regression

# # get ORs
# ORs = exp(log.out$coeff)

# # get p-values
# summary(log.out)

# # get 95% CIs
# exp(confint(log.out))

# read in and organize the data
#AF_data = read.csv("/athena/masonlab/scratch/users/nai2008/Anesthesiology/SafavyniaLab/AF/AF_data_to_read_in.csv")
AF_data = read.csv("/Users/nikolayivanov/Desktop/Anesthesiology_copyFromCluster/SafavyniaLab/AF/AF_data_to_read_in_with_sedation_info.csv")
AF_data$Sex = factor(AF_data$Sex, levels=c("F","M"))
AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg = AF_data$Weight_adjusted_total_AF_dose_g_per_kg * 10
AF_data$Weight_adjusted_total_AF_dose_10mg_per_kg = AF_data$Weight_adjusted_total_AF_dose_g_per_kg * 100

AF_data=AF_data[-which(AF_data$Acute_Infarction_brain=="N/A" | AF_data$Acute_Hemorrhage_brain=="N/A"),]

AF_data$Acute_Hemorrhage_brain[which(AF_data$Acute_Hemorrhage_brain != 0)] = 1
AF_data$Acute_Hemorrhage_brain = as.integer(AF_data$Acute_Hemorrhage_brain)

AF_data$Acute_Infarction_brain = as.integer(AF_data$Acute_Infarction_brain)

AF_data$CPB_yes_or_no = AF_data$Total_CPB_Time_mins
izero = which(AF_data$CPB_yes_or_no==0)
inonzero = which(AF_data$CPB_yes_or_no != 0)
AF_data$CPB_yes_or_no[izero] = 0
AF_data$CPB_yes_or_no[inonzero] = 1

AF_data$Race = NA
AF_data$Race[which(AF_data$White !=0)] = 'White'
AF_data$Race[which(AF_data$Black !=0)] = 'Black'
AF_data$Race[which(AF_data$Asian !=0)] = 'Asian'
AF_data$Race[is.na(AF_data$Race)] = 'Other'
AF_data$Race = factor(AF_data$Race, levels = c("White", "Black", "Asian", "Other"))

length(which(AF_data$Race =='Other')) # Unknown race, Other race, or no race data
# 80

#write.csv(AF_data$Race, file='/Users/nikolayivanov/Desktop/Anesthesiology_copyFromCluster/SafavyniaLab/AF/race.csv', row.names = FALSE)

# add markers of cortical hyperexcitability
AF_data$Markers_of_cortical_hyperexcitability = rep(0, times = nrow(AF_data))
AF_data$Generalized_Marker_Cortical_Hyperexcitability = rep(0, times = nrow(AF_data))
AF_data$Lateralized_Marker_Cortical_Hyperexcitability = rep(0, times = nrow(AF_data))

generalized_CH = c('GPDs','GSWs')
lateralized_CH = c('MfSWs', 'LPDs', 'BIPDs', 'LRDA', 'BIRDA', 'LSWs', 'BISWs')
other_CH = c('GRDA', 'Triphasics')

AF_data_CH = AF_data[,30:40]
temp = rowSums(AF_data_CH)
AF_data$Markers_of_cortical_hyperexcitability[which(temp != 0)] = 1

AF_data_CH_gen = AF_data_CH[,colnames(AF_data_CH) %in% generalized_CH]
temp = rowSums(AF_data_CH_gen)
AF_data$Generalized_Marker_Cortical_Hyperexcitability[which(temp != 0)] = 1

AF_data_CH_lat = AF_data_CH[,colnames(AF_data_CH) %in% lateralized_CH]
temp = rowSums(AF_data_CH_lat)
AF_data$Lateralized_Marker_Cortical_Hyperexcitability[which(temp != 0)] = 1

# 224 samples; 26 seizure events; 102 have any markers of cortical hyperexcitability

# cohort: 224: + seizure 26; - seizure 198

# of the pts who had an electrographic seizure, how much had witnessed convuslions
table(AF_data$Witnessed_Convulsion[which(AF_data$Electrographic_Seizure == 1)])
 # 0  1 
 # 8 18

table(AF_data$Witnessed_Convulsion[which(AF_data$Electrographic_Seizure == 0)])
#   0   1 
# 153  45

# of the pts who had an electrographic seizure, how much had a witnessed event suspicious for seizure
table(AF_data$Witnessed_Event_Suspicious_for_Seizure[which(AF_data$Electrographic_Seizure == 1)])
 # 0  1 
 # 4 22 

table(AF_data$Witnessed_Event_Suspicious_for_Seizure[which(AF_data$Electrographic_Seizure == 0)])
#  0  1  2 
# 92 94 12 

#sedation?
AF_data$any_benzo_coded = AF_data$Benzodiazepine_coded_yes_no + AF_data$Midazolam_yes_no

length(which(AF_data$Dexmedetomidine_yes_no != 0))
(length(which(AF_data$Dexmedetomidine_yes_no != 0))/nrow(AF_data))*100

length(which(AF_data$Propofol_yes_no != 0))
(length(which(AF_data$Propofol_yes_no != 0))/nrow(AF_data))*100

length(which(AF_data$any_benzo_coded != 0))
(length(which(AF_data$any_benzo_coded != 0))/nrow(AF_data))*100

library(bda)
mediation.test(AF_data$Acute_Infarction_brain, AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg, AF_data$Electrographic_Seizure)
#              Sobel     Aroian    Goodman
# z.value 2.02090409 1.97863364 2.06600484
# p.value 0.04328969 0.04785727 0.03882802

#######################################################################################################
### Demographics table
#######################################################################################################

indexACA = which(AF_data$Antifibrinolytic == 1)
indexNoACA = which(AF_data$Antifibrinolytic == 0)

nrow(AF_data)
length(indexACA)
length(indexNoACA)

################## Age

# All
median(AF_data$Age)
quantile(AF_data$Age)

# + ACA
median(AF_data$Age[indexACA])
quantile(AF_data$Age[indexACA])

# - ACA
median(AF_data$Age[indexNoACA])
quantile(AF_data$Age[indexNoACA])

shapiro.test(AF_data$Age) # p-value = 3.375e-07
shapiro.test(AF_data$Age[indexACA]) # p-value = 2.842e-06
shapiro.test(AF_data$Age[indexNoACA]) # p-value = 0.009158

wilcox.test(AF_data$Age[indexACA], AF_data$Age[indexNoACA], alternative = "two.sided") # p-value = 0.03885
#t.test(AF_data$Age[indexACA], AF_data$Age[indexNoACA], alternative = "two.sided")$p.value

################## Male sex

# All
length(which(AF_data$Sex == 'M'))
(length(which(AF_data$Sex == 'M'))/nrow(AF_data))*100

# + ACA
length(which(AF_data$Sex[indexACA] == 'M'))
(length(which(AF_data$Sex[indexACA] == 'M'))/length(indexACA))*100

# - ACA
length(which(AF_data$Sex[indexNoACA] == 'M'))
(length(which(AF_data$Sex[indexNoACA] == 'M'))/length(indexNoACA))*100

fisher.test(table(AF_data$Sex, AF_data$Antifibrinolytic))$p.value

################## White race

# All
length(which(AF_data$White_race == 1))
(length(which(AF_data$White_race == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$White_race[indexACA] == 1))
(length(which(AF_data$White_race[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$White_race[indexNoACA] == 1))
(length(which(AF_data$White_race[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$White_race, AF_data$Antifibrinolytic))$p.value

################## Asian race

# All
length(which(AF_data$Asian_race == 1))
(length(which(AF_data$Asian_race == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$Asian_race[indexACA] == 1))
(length(which(AF_data$Asian_race[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$Asian_race[indexNoACA] == 1))
(length(which(AF_data$Asian_race[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$Asian_race, AF_data$Antifibrinolytic))$p.value

################## Black race

# All
length(which(AF_data$Black_race == 1))
(length(which(AF_data$Black_race == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$Black_race[indexACA] == 1))
(length(which(AF_data$Black_race[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$Black_race[indexNoACA] == 1))
(length(which(AF_data$Black_race[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$Black_race, AF_data$Antifibrinolytic))$p.value

################## Weight

# All
median(AF_data$Weight_Kg)
quantile(AF_data$Weight_Kg)

# + ACA
median(AF_data$Weight_Kg[indexACA])
quantile(AF_data$Weight_Kg[indexACA])

# - ACA
median(AF_data$Weight_Kg[indexNoACA])
quantile(AF_data$Weight_Kg[indexNoACA])

shapiro.test(AF_data$Weight_Kg) # p-value = 8.672e-09
shapiro.test(AF_data$Weight_Kg[indexACA]) # p-value = 0.0004118
shapiro.test(AF_data$Weight_Kg[indexNoACA]) # p-value = 3.839e-08

wilcox.test(AF_data$Weight_Kg[indexACA], AF_data$Weight_Kg[indexNoACA], alternative = "two.sided") # p-value = 0.9629
#t.test(AF_data$Weight_Kg[indexACA], AF_data$Weight_Kg[indexNoACA], alternative = "two.sided")$p.value

################## CBP (yes/no)

# All
length(which(AF_data$CPB_yes_or_no == 1))
(length(which(AF_data$CPB_yes_or_no == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$CPB_yes_or_no[indexACA] == 1))
(length(which(AF_data$CPB_yes_or_no[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$CPB_yes_or_no[indexNoACA] == 1))
(length(which(AF_data$CPB_yes_or_no[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$CPB_yes_or_no, AF_data$Antifibrinolytic))$p.value

################## CPB time (min)

# All
median(AF_data$Total_CPB_Time_mins)
quantile(AF_data$Total_CPB_Time_mins)

# + ACA
median(AF_data$Total_CPB_Time_mins[indexACA])
quantile(AF_data$Total_CPB_Time_mins[indexACA])

# - ACA
median(AF_data$Total_CPB_Time_mins[indexNoACA])
quantile(AF_data$Total_CPB_Time_mins[indexNoACA])

shapiro.test(AF_data$Total_CPB_Time_mins) # p-value < 2.2e-16
shapiro.test(AF_data$Total_CPB_Time_mins[indexACA]) # p-value = 2.666e-07
shapiro.test(AF_data$Total_CPB_Time_mins[indexNoACA]) # p-value < 2.2e-16

wilcox.test(AF_data$Total_CPB_Time_mins[indexACA], AF_data$Total_CPB_Time_mins[indexNoACA], alternative = "two.sided") # p-value < 2.2e-16
#t.test(AF_data$Total_CPB_Time_mins[indexACA], AF_data$Total_CPB_Time_mins[indexNoACA], alternative = "two.sided")$p.value

################## Hx of stroke prior to admission

# All
length(which(AF_data$Stroke_Prior_to_Admit == 1))
(length(which(AF_data$Stroke_Prior_to_Admit == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$Stroke_Prior_to_Admit[indexACA] == 1))
(length(which(AF_data$Stroke_Prior_to_Admit[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$Stroke_Prior_to_Admit[indexNoACA] == 1))
(length(which(AF_data$Stroke_Prior_to_Admit[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$Stroke_Prior_to_Admit, AF_data$Antifibrinolytic))$p.value

################## Hx of seizure/epilepsy prior to admission

# All
length(which(AF_data$Seizure_or_Epilepsy_Hx == 1))
(length(which(AF_data$Seizure_or_Epilepsy_Hx == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$Seizure_or_Epilepsy_Hx[indexACA] == 1))
(length(which(AF_data$Seizure_or_Epilepsy_Hx[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$Seizure_or_Epilepsy_Hx[indexNoACA] == 1))
(length(which(AF_data$Seizure_or_Epilepsy_Hx[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$Seizure_or_Epilepsy_Hx, AF_data$Antifibrinolytic))$p.value

################## Brain infarct during index admission

# All
length(which(AF_data$Acute_Infarction_brain == 1))
(length(which(AF_data$Acute_Infarction_brain == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$Acute_Infarction_brain[indexACA] == 1))
(length(which(AF_data$Acute_Infarction_brain[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$Acute_Infarction_brain[indexNoACA] == 1))
(length(which(AF_data$Acute_Infarction_brain[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$Acute_Infarction_brain, AF_data$Antifibrinolytic))$p.value

################## Brain hemorrhage during index admission

# All
length(which(AF_data$Acute_Hemorrhage_brain == 1))
(length(which(AF_data$Acute_Hemorrhage_brain == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$Acute_Hemorrhage_brain[indexACA] == 1))
(length(which(AF_data$Acute_Hemorrhage_brain[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$Acute_Hemorrhage_brain[indexNoACA] == 1))
(length(which(AF_data$Acute_Hemorrhage_brain[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$Acute_Hemorrhage_brain, AF_data$Antifibrinolytic))$p.value

################## Electrographic seizure

# All
length(which(AF_data$Electrographic_Seizure == 1))
(length(which(AF_data$Electrographic_Seizure == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$Electrographic_Seizure[indexACA] == 1))
(length(which(AF_data$Electrographic_Seizure[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$Electrographic_Seizure[indexNoACA] == 1))
(length(which(AF_data$Electrographic_Seizure[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$Electrographic_Seizure, AF_data$Antifibrinolytic))$p.value

################## (Any) cortical hyperexcitability markers

# All
length(which(AF_data$Markers_of_cortical_hyperexcitability == 1))
(length(which(AF_data$Markers_of_cortical_hyperexcitability == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$Markers_of_cortical_hyperexcitability[indexACA] == 1))
(length(which(AF_data$Markers_of_cortical_hyperexcitability[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$Markers_of_cortical_hyperexcitability[indexNoACA] == 1))
(length(which(AF_data$Markers_of_cortical_hyperexcitability[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$Markers_of_cortical_hyperexcitability, AF_data$Antifibrinolytic))$p.value

################## Generalized cortical hyperexcitability markers

# All
length(which(AF_data$Generalized_Marker_Cortical_Hyperexcitability == 1))
(length(which(AF_data$Generalized_Marker_Cortical_Hyperexcitability == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$Generalized_Marker_Cortical_Hyperexcitability[indexACA] == 1))
(length(which(AF_data$Generalized_Marker_Cortical_Hyperexcitability[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$Generalized_Marker_Cortical_Hyperexcitability[indexNoACA] == 1))
(length(which(AF_data$Generalized_Marker_Cortical_Hyperexcitability[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$Generalized_Marker_Cortical_Hyperexcitability, AF_data$Antifibrinolytic))$p.value

################## Lateralized cortical hyperexcitability markers

# All
length(which(AF_data$Lateralized_Marker_Cortical_Hyperexcitability == 1))
(length(which(AF_data$Lateralized_Marker_Cortical_Hyperexcitability == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$Lateralized_Marker_Cortical_Hyperexcitability[indexACA] == 1))
(length(which(AF_data$Lateralized_Marker_Cortical_Hyperexcitability[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$Lateralized_Marker_Cortical_Hyperexcitability[indexNoACA] == 1))
(length(which(AF_data$Lateralized_Marker_Cortical_Hyperexcitability[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$Lateralized_Marker_Cortical_Hyperexcitability, AF_data$Antifibrinolytic))$p.value

################## 'Other' cortical hyperexcitability markers
AF_data$Other_Marker_Cortical_Hyperexcitability = rep(0, times=nrow(AF_data))
AF_data$Other_Marker_Cortical_Hyperexcitability[which(AF_data$GRDA !=0 | AF_data$Triphasics !=0)] = 1

# All
length(which(AF_data$Other_Marker_Cortical_Hyperexcitability == 1))
(length(which(AF_data$Other_Marker_Cortical_Hyperexcitability == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$Other_Marker_Cortical_Hyperexcitability[indexACA] == 1))
(length(which(AF_data$Other_Marker_Cortical_Hyperexcitability[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$Other_Marker_Cortical_Hyperexcitability[indexNoACA] == 1))
(length(which(AF_data$Other_Marker_Cortical_Hyperexcitability[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$Other_Marker_Cortical_Hyperexcitability, AF_data$Antifibrinolytic))$p.value

################## GPDs (generalized CH marker)

# All
length(which(AF_data$GPDs == 1))
(length(which(AF_data$GPDs == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$GPDs[indexACA] == 1))
(length(which(AF_data$GPDs[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$GPDs[indexNoACA] == 1))
(length(which(AF_data$GPDs[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$GPDs, AF_data$Antifibrinolytic))$p.value

################## GSWs (generalized CH marker)

# All
length(which(AF_data$GSWs == 1))
(length(which(AF_data$GSWs == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$GSWs[indexACA] == 1))
(length(which(AF_data$GSWs[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$GSWs[indexNoACA] == 1))
(length(which(AF_data$GSWs[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$GSWs, AF_data$Antifibrinolytic))$p.value

################## MfSWs (lateralized CH marker)

# All
length(which(AF_data$MfSWs == 1))
(length(which(AF_data$MfSWs == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$MfSWs[indexACA] == 1))
(length(which(AF_data$MfSWs[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$MfSWs[indexNoACA] == 1))
(length(which(AF_data$MfSWs[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$MfSWs, AF_data$Antifibrinolytic))$p.value

################## LPDs (lateralized CH marker)

# All
length(which(AF_data$LPDs == 1))
(length(which(AF_data$LPDs == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$LPDs[indexACA] == 1))
(length(which(AF_data$LPDs[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$LPDs[indexNoACA] == 1))
(length(which(AF_data$LPDs[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$LPDs, AF_data$Antifibrinolytic))$p.value

################## BIPDs (lateralized CH marker)

# All
length(which(AF_data$BIPDs == 1))
(length(which(AF_data$BIPDs == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$BIPDs[indexACA] == 1))
(length(which(AF_data$BIPDs[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$BIPDs[indexNoACA] == 1))
(length(which(AF_data$BIPDs[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$BIPDs, AF_data$Antifibrinolytic))$p.value

################## LRDA (lateralized CH marker)

# All
length(which(AF_data$LRDA == 1))
(length(which(AF_data$LRDA == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$LRDA[indexACA] == 1))
(length(which(AF_data$LRDA[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$LRDA[indexNoACA] == 1))
(length(which(AF_data$LRDA[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$LRDA, AF_data$Antifibrinolytic))$p.value

################## BIRDA (lateralized CH marker)

# All
length(which(AF_data$BIRDA == 1))
(length(which(AF_data$BIRDA == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$BIRDA[indexACA] == 1))
(length(which(AF_data$BIRDA[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$BIRDA[indexNoACA] == 1))
(length(which(AF_data$BIRDA[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$BIRDA, AF_data$Antifibrinolytic))$p.value

################## LSWs (lateralized CH marker)

# All
length(which(AF_data$LSWs == 1))
(length(which(AF_data$LSWs == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$LSWs[indexACA] == 1))
(length(which(AF_data$LSWs[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$LSWs[indexNoACA] == 1))
(length(which(AF_data$LSWs[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$LSWs, AF_data$Antifibrinolytic))$p.value

################## BISWs (lateralized CH marker)

# All
length(which(AF_data$BISWs == 1))
(length(which(AF_data$BISWs == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$BISWs[indexACA] == 1))
(length(which(AF_data$BISWs[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$BISWs[indexNoACA] == 1))
(length(which(AF_data$BISWs[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$BISWs, AF_data$Antifibrinolytic))$p.value

################## GRDA ('other' CH marker)

# All
length(which(AF_data$GRDA == 1))
(length(which(AF_data$GRDA == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$GRDA[indexACA] == 1))
(length(which(AF_data$GRDA[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$GRDA[indexNoACA] == 1))
(length(which(AF_data$GRDA[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$GRDA, AF_data$Antifibrinolytic))$p.value

################## Triphasics ('other' CH marker)

# All
length(which(AF_data$Triphasics == 1))
(length(which(AF_data$Triphasics == 1))/nrow(AF_data))*100

# + ACA
length(which(AF_data$Triphasics[indexACA] == 1))
(length(which(AF_data$Triphasics[indexACA] == 1))/length(indexACA))*100

# - ACA
length(which(AF_data$Triphasics[indexNoACA] == 1))
(length(which(AF_data$Triphasics[indexNoACA] == 1))/length(indexNoACA))*100

fisher.test(table(AF_data$Triphasics, AF_data$Antifibrinolytic))$p.value

################## Total ACA dose (grams)

# All
#N/A

# + ACA
median(AF_data$Total_AF_Dose_grams[indexACA])
quantile(AF_data$Total_AF_Dose_grams[indexACA])

# - ACA
#N/A

#######################################################################################################
### Univariate analyses (outcome: Electrographic_Seizure)
#######################################################################################################

###### AF: yes or no
log.out = glm(AF_data$Electrographic_Seizure ~ AF_data$Antifibrinolytic, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))


###### AF dose [g/kg]
log.out = glm(AF_data$Electrographic_Seizure ~ AF_data$Weight_adjusted_total_AF_dose_g_per_kg, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### AF dose [100mg/kg]
log.out = glm(AF_data$Electrographic_Seizure ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### AF dose [10mg/kg]
log.out = glm(AF_data$Electrographic_Seizure ~ AF_data$Weight_adjusted_total_AF_dose_10mg_per_kg, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### CPB_yes_or_no
log.out = glm(AF_data$Electrographic_Seizure ~ AF_data$CPB_yes_or_no, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Total_CPB_Time_mins
log.out = glm(AF_data$Electrographic_Seizure ~ AF_data$Total_CPB_Time_mins, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Sex
log.out = glm(AF_data$Electrographic_Seizure ~ AF_data$Sex, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Age
log.out = glm(AF_data$Electrographic_Seizure ~ AF_data$Age, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Stroke_Prior_to_Admit
log.out = glm(AF_data$Electrographic_Seizure ~ AF_data$Stroke_Prior_to_Admit, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Acute_Hemorrhage_brain
log.out = glm(AF_data$Electrographic_Seizure ~ AF_data$Acute_Hemorrhage_brain, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Acute_Infarction_brain
log.out = glm(AF_data$Electrographic_Seizure ~ AF_data$Acute_Infarction_brain, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Seizure_or_Epilepsy_Hx
log.out = glm(AF_data$Electrographic_Seizure ~ AF_data$Seizure_or_Epilepsy_Hx, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Race
sz = AF_data$Electrographic_Seizure[-which(AF_data$Race == 'Other')]
rr = AF_data$Race[-which(AF_data$Race == 'Other')]
rr = factor(rr, levels=c('White', 'Black', 'Asian'))

log.out = glm(sz ~ rr, family = binomial (link=logit))

# get stats
exp(log.out$coeff)
summary(log.out)$coeff
exp(confint(log.out))

#######################################################################################################
### Univariate analyses (outcome: Markers_of_cortical_hyperexcitability)
#######################################################################################################

###### AF: yes or no
log.out = glm(AF_data$Markers_of_cortical_hyperexcitability ~ AF_data$Antifibrinolytic, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))


###### AF dose [g/kg]
log.out = glm(AF_data$Markers_of_cortical_hyperexcitability ~ AF_data$Weight_adjusted_total_AF_dose_g_per_kg, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### AF dose [100mg/kg]
log.out = glm(AF_data$Markers_of_cortical_hyperexcitability ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### AF dose [10mg/kg]
log.out = glm(AF_data$Markers_of_cortical_hyperexcitability ~ AF_data$Weight_adjusted_total_AF_dose_10mg_per_kg, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### CPB_yes_or_no
log.out = glm(AF_data$Markers_of_cortical_hyperexcitability ~ AF_data$CPB_yes_or_no, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Total_CPB_Time_mins
log.out = glm(AF_data$Markers_of_cortical_hyperexcitability ~ AF_data$Total_CPB_Time_mins, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Sex
log.out = glm(AF_data$Markers_of_cortical_hyperexcitability ~ AF_data$Sex, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Age
log.out = glm(AF_data$Markers_of_cortical_hyperexcitability ~ AF_data$Age, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Stroke_Prior_to_Admit
log.out = glm(AF_data$Markers_of_cortical_hyperexcitability ~ AF_data$Stroke_Prior_to_Admit, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Acute_Hemorrhage_brain
log.out = glm(AF_data$Markers_of_cortical_hyperexcitability ~ AF_data$Acute_Hemorrhage_brain, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Acute_Infarction_brain
log.out = glm(AF_data$Markers_of_cortical_hyperexcitability ~ AF_data$Acute_Infarction_brain, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Seizure_or_Epilepsy_Hx
log.out = glm(AF_data$Markers_of_cortical_hyperexcitability ~ AF_data$Seizure_or_Epilepsy_Hx, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Race
ch = AF_data$Markers_of_cortical_hyperexcitability[-which(AF_data$Race == 'Other')]
rr = AF_data$Race[-which(AF_data$Race == 'Other')]
rr = factor(rr, levels=c('White', 'Black', 'Asian'))

log.out = glm(ch ~ rr, family = binomial (link=logit))

# get stats
exp(log.out$coeff)
summary(log.out)$coeff
exp(confint(log.out))

#######################################################################################################
### Univariate analyses (outcome: Generalized_Marker_Cortical_Hyperexcitability)
#######################################################################################################

###### AF: yes or no
log.out = glm(AF_data$Generalized_Marker_Cortical_Hyperexcitability ~ AF_data$Antifibrinolytic, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))


###### AF dose [g/kg]
log.out = glm(AF_data$Generalized_Marker_Cortical_Hyperexcitability ~ AF_data$Weight_adjusted_total_AF_dose_g_per_kg, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### AF dose [100mg/kg]
log.out = glm(AF_data$Generalized_Marker_Cortical_Hyperexcitability ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### AF dose [10mg/kg]
log.out = glm(AF_data$Generalized_Marker_Cortical_Hyperexcitability ~ AF_data$Weight_adjusted_total_AF_dose_10mg_per_kg, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### CPB_yes_or_no
log.out = glm(AF_data$Generalized_Marker_Cortical_Hyperexcitability ~ AF_data$CPB_yes_or_no, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Total_CPB_Time_mins
log.out = glm(AF_data$Generalized_Marker_Cortical_Hyperexcitability ~ AF_data$Total_CPB_Time_mins, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Sex
log.out = glm(AF_data$Generalized_Marker_Cortical_Hyperexcitability ~ AF_data$Sex, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Age
log.out = glm(AF_data$Generalized_Marker_Cortical_Hyperexcitability ~ AF_data$Age, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Stroke_Prior_to_Admit
log.out = glm(AF_data$Generalized_Marker_Cortical_Hyperexcitability ~ AF_data$Stroke_Prior_to_Admit, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Acute_Hemorrhage_brain
log.out = glm(AF_data$Generalized_Marker_Cortical_Hyperexcitability ~ AF_data$Acute_Hemorrhage_brain, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Acute_Infarction_brain
log.out = glm(AF_data$Generalized_Marker_Cortical_Hyperexcitability ~ AF_data$Acute_Infarction_brain, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Seizure_or_Epilepsy_Hx
log.out = glm(AF_data$Generalized_Marker_Cortical_Hyperexcitability ~ AF_data$Seizure_or_Epilepsy_Hx, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Race
ch = AF_data$Generalized_Marker_Cortical_Hyperexcitability[-which(AF_data$Race == 'Other')]
rr = AF_data$Race[-which(AF_data$Race == 'Other')]
rr = factor(rr, levels=c('White', 'Black', 'Asian'))

log.out = glm(ch ~ rr, family = binomial (link=logit))

# get stats
exp(log.out$coeff)
summary(log.out)$coeff
exp(confint(log.out))

#######################################################################################################
### Univariate analyses (outcome: Lateralized_Marker_Cortical_Hyperexcitability)
#######################################################################################################

###### AF: yes or no
log.out = glm(AF_data$Lateralized_Marker_Cortical_Hyperexcitability ~ AF_data$Antifibrinolytic, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))


###### AF dose [g/kg]
log.out = glm(AF_data$Lateralized_Marker_Cortical_Hyperexcitability ~ AF_data$Weight_adjusted_total_AF_dose_g_per_kg, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### AF dose [100mg/kg]
log.out = glm(AF_data$Lateralized_Marker_Cortical_Hyperexcitability ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### AF dose [10mg/kg]
log.out = glm(AF_data$Lateralized_Marker_Cortical_Hyperexcitability ~ AF_data$Weight_adjusted_total_AF_dose_10mg_per_kg, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### CPB_yes_or_no
log.out = glm(AF_data$Lateralized_Marker_Cortical_Hyperexcitability ~ AF_data$CPB_yes_or_no, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Total_CPB_Time_mins
log.out = glm(AF_data$Lateralized_Marker_Cortical_Hyperexcitability ~ AF_data$Total_CPB_Time_mins, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Sex
log.out = glm(AF_data$Lateralized_Marker_Cortical_Hyperexcitability ~ AF_data$Sex, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Age
log.out = glm(AF_data$Lateralized_Marker_Cortical_Hyperexcitability ~ AF_data$Age, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Stroke_Prior_to_Admit
log.out = glm(AF_data$Lateralized_Marker_Cortical_Hyperexcitability ~ AF_data$Stroke_Prior_to_Admit, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Acute_Hemorrhage_brain
log.out = glm(AF_data$Lateralized_Marker_Cortical_Hyperexcitability ~ AF_data$Acute_Hemorrhage_brain, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Acute_Infarction_brain
log.out = glm(AF_data$Lateralized_Marker_Cortical_Hyperexcitability ~ AF_data$Acute_Infarction_brain, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Seizure_or_Epilepsy_Hx
log.out = glm(AF_data$Lateralized_Marker_Cortical_Hyperexcitability ~ AF_data$Seizure_or_Epilepsy_Hx, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Race
ch = AF_data$Lateralized_Marker_Cortical_Hyperexcitability[-which(AF_data$Race == 'Other')]
rr = AF_data$Race[-which(AF_data$Race == 'Other')]
rr = factor(rr, levels=c('White', 'Black', 'Asian'))

log.out = glm(ch ~ rr, family = binomial (link=logit))

# get stats
exp(log.out$coeff)
summary(log.out)$coeff
exp(confint(log.out))

#######################################################################################################
### Univariate analyses (outcome: Acute_Infarction_brain)
#######################################################################################################

###### AF: yes or no
log.out = glm(AF_data$Acute_Infarction_brain ~ AF_data$Antifibrinolytic, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))


###### AF dose [g/kg]
log.out = glm(AF_data$Acute_Infarction_brain ~ AF_data$Weight_adjusted_total_AF_dose_g_per_kg, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### AF dose [100mg/kg]
log.out = glm(AF_data$Acute_Infarction_brain ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### AF dose [10mg/kg]
log.out = glm(AF_data$Acute_Infarction_brain ~ AF_data$Weight_adjusted_total_AF_dose_10mg_per_kg, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### CPB_yes_or_no
log.out = glm(AF_data$Acute_Infarction_brain ~ AF_data$CPB_yes_or_no, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Total_CPB_Time_mins
log.out = glm(AF_data$Acute_Infarction_brain ~ AF_data$Total_CPB_Time_mins, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Sex
log.out = glm(AF_data$Acute_Infarction_brain ~ AF_data$Sex, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Age
log.out = glm(AF_data$Acute_Infarction_brain ~ AF_data$Age, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Stroke_Prior_to_Admit
log.out = glm(AF_data$Acute_Infarction_brain ~ AF_data$Stroke_Prior_to_Admit, family = binomial (link=logit))

# get stats
as.numeric(exp(log.out$coeff)[2])
summary(log.out)$coeff[2,4]
exp(confint(log.out))

###### Race
stroke = AF_data$Acute_Infarction_brain[-which(AF_data$Race == 'Other')]
rr = AF_data$Race[-which(AF_data$Race == 'Other')]
rr = factor(rr, levels=c('White', 'Black', 'Asian'))

log.out = glm(stroke ~ rr, family = binomial (link=logit))

# get stats
exp(log.out$coeff)
summary(log.out)$coeff
exp(confint(log.out))

#######################################################################################################
### Multivariate analyses (outcome: Electrographic_Seizure)
#######################################################################################################

log.out = glm(AF_data$Electrographic_Seizure ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg + AF_data$Total_CPB_Time_mins + AF_data$Acute_Infarction_brain, family = binomial (link=logit))

# get stats
exp(log.out$coeff)
summary(log.out)$coeff
exp(confint(log.out))

# > exp(log.out$coeff)
#                                        (Intercept)
#                                         0.04764651
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg
#                                         1.49635908
#                        AF_data$Total_CPB_Time_mins
#                                         1.00228622
#                     AF_data$Acute_Infarction_brain
#                                         2.59877857

# > summary(log.out)$coeff
#                                                        Estimate  Std. Error
# (Intercept)                                        -3.043945910 0.369698566
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg  0.403034879 0.184329682
# AF_data$Total_CPB_Time_mins                         0.002283614 0.003259187
# AF_data$Acute_Infarction_brain                      0.955041554 0.472747391
#                                                       z value     Pr(>|z|)
# (Intercept)                                        -8.2335886 1.816868e-16
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg  2.1864893 2.877983e-02
# AF_data$Total_CPB_Time_mins                         0.7006698 4.835091e-01
# AF_data$Acute_Infarction_brain                      2.0201942 4.336324e-02

# > exp(confint(log.out))
# Waiting for profiling to be done...
#                                                         2.5 %     97.5 %
# (Intercept)                                        0.02153716 0.09279836
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg 1.04416689 2.17358187
# AF_data$Total_CPB_Time_mins                        0.99550426 1.00849141
# AF_data$Acute_Infarction_brain                     1.00673164 6.52052138

# looks like the variable Total_CPB_Time_mins does not provide extra information that helps us predict the chance of seiuzre
# lets remove it and compare the new model to the old model (via the likelihod ratio test)

log.out_AF_CPB_infarct = glm(AF_data$Electrographic_Seizure ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg + AF_data$Total_CPB_Time_mins + AF_data$Acute_Infarction_brain, family = binomial (link=logit))
log.out_AF_infarct = glm(AF_data$Electrographic_Seizure ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg + AF_data$Acute_Infarction_brain, family = binomial (link=logit))

library(lmtest)
lrtest(log.out_AF_infarct, log.out_AF_CPB_infarct)
#   #Df  LogLik Df Chisq Pr(>Chisq)
# 1   3 -69.690
# 2   4 -69.454  1 0.473     0.4916

# so based on the likelihood ratio test, we can remove Total_CPB_Time_mins from the model

# let's try a different model, this time using a binary variable for CPB (yes/no)

log.out = glm(AF_data$Electrographic_Seizure ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg + AF_data$Total_CPB_Time_mins + AF_data$Acute_Infarction_brain, family = binomial (link=logit))

# get stats
exp(log.out$coeff)
summary(log.out)$coeff
exp(confint(log.out))

log.out = glm(AF_data$Electrographic_Seizure ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg + AF_data$CPB_yes_or_no + AF_data$Acute_Infarction_brain, family = binomial (link=logit))

# get stats
exp(log.out$coeff)
summary(log.out)$coeff
exp(confint(log.out))

# > exp(log.out$coeff)
#                                        (Intercept)
#                                         0.04414003
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg
#                                         1.43840878
#                              AF_data$CPB_yes_or_no
#                                         1.72889499
#                     AF_data$Acute_Infarction_brain
#                                         2.57002320
# > summary(log.out)$coeff
#                                                      Estimate Std. Error
# (Intercept)                                        -3.1203881  0.4100930
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg  0.3635375  0.2072151
# AF_data$CPB_yes_or_no                               0.5474825  0.6859386
# AF_data$Acute_Infarction_brain                      0.9439149  0.4716583
#                                                       z value     Pr(>|z|)
# (Intercept)                                        -7.6089766 2.762748e-14
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg  1.7543964 7.936260e-02
# AF_data$CPB_yes_or_no                               0.7981508 4.247830e-01
# AF_data$Acute_Infarction_brain                      2.0012687 4.536344e-02
# > exp(confint(log.out))
# Waiting for profiling to be done...
#                                                         2.5 %    97.5 %
# (Intercept)                                        0.01793319 0.0913696
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg 0.95886225 2.1841632
# AF_data$CPB_yes_or_no                              0.44411147 6.7554291
# AF_data$Acute_Infarction_brain                     0.99856809 6.4366953

# still doesn't add any new info to the model

log.out_AF_CPB_infarct = glm(AF_data$Electrographic_Seizure ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg + AF_data$CPB_yes_or_no + AF_data$Acute_Infarction_brain, family = binomial (link=logit))
log.out_AF_infarct = glm(AF_data$Electrographic_Seizure ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg + AF_data$Acute_Infarction_brain, family = binomial (link=logit))

lrtest(log.out_AF_infarct, log.out_AF_CPB_infarct)
#   #Df  LogLik Df Chisq Pr(>Chisq)
# 1   3 -69.690
# 2   4 -69.372  1 0.637     0.4248

# so our preliminarty main effects model is:
log.out_AF_infarct = glm(AF_data$Electrographic_Seizure ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg + AF_data$Acute_Infarction_brain, family = binomial (link=logit))

# check that Weight_adjusted_total_AF_dose_100mg_per_kg satisfies the assumption of linearity in the logit

q = quantile(AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg, probs=c(0,0.25,0.5,0.75,1))
#       0%      25%      50%      75%     100%
# 0.000000 0.000000 0.000000 2.033584 7.593052

zeroIndexes = which(AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg == 0)
#length(zeroIndexes) 138

p1 = mean(AF_data$Electrographic_Seizure[zeroIndexes])
p2 = mean( AF_data$Electrographic_Seizure[which(AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg > 0 & AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg <= q[4])] )
p3 = mean( AF_data$Electrographic_Seizure[which(AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg > q[4] & AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg <= q[5])] )

probs = c(p1,p2,p3)
logits = log(probs/(1-probs))

medians = c(0, 
	median(AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg[which(AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg > 0 & AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg <= q[4])]),
	median(AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg[which(AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg > q[4] & AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg <= q[5])])
)

pdf('/athena/masonlab/scratch/users/nai2008/Anesthesiology/SafavyniaLab/AF/AF_dose_logit_linearity_check.pdf')
plot(medians,logits,xlab='Medians of quartiles (weight adjusted total AF dose (100mg/kg)',ylab='Log-odds of outcome (seizure)')
dev.off()

# check for interaction
log.out_AF_infarct = glm(AF_data$Electrographic_Seizure ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg + AF_data$Acute_Infarction_brain + AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg*AF_data$Acute_Infarction_brain, family = binomial (link=logit))

# get stats
exp(log.out_AF_infarct$coeff)
summary(log.out_AF_infarct)$coeff
exp(confint(log.out_AF_infarct))

# > # get stats
# > exp(log.out_AF_infarct$coeff)
#                                                                       (Intercept)
#                                                                        0.05903991
#                                AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg
#                                                                        1.45421398
#                                                    AF_data$Acute_Infarction_brain
#                                                                        1.60185370
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg:AF_data$Acute_Infarction_brain
#                                                                        1.29184610
# > summary(log.out_AF_infarct)$coeff
#                                                                                     Estimate
# (Intercept)                                                                       -2.8295416
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg                                 0.3744655
# AF_data$Acute_Infarction_brain                                                     0.4711615
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg:AF_data$Acute_Infarction_brain  0.2560723
#                                                                                   Std. Error
# (Intercept)                                                                        0.3757306
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg                                 0.1963137
# AF_data$Acute_Infarction_brain                                                     0.7658137
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg:AF_data$Acute_Infarction_brain  0.3101855
#                                                                                      z value
# (Intercept)                                                                       -7.5307716
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg                                 1.9074856
# AF_data$Acute_Infarction_brain                                                     0.6152430
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg:AF_data$Acute_Infarction_brain  0.8255456
#                                                                                       Pr(>|z|)
# (Intercept)                                                                       5.044135e-14
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg                                5.645774e-02
# AF_data$Acute_Infarction_brain                                                    5.383943e-01
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg:AF_data$Acute_Infarction_brain 4.090619e-01
# > exp(confint(log.out_AF_infarct))
# Waiting for profiling to be done...
#                                                                                        2.5 %
# (Intercept)                                                                       0.02606516
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg                                0.97228761
# AF_data$Acute_Infarction_brain                                                    0.29856941
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg:AF_data$Acute_Infarction_brain 0.72532088
#                                                                                      97.5 %
# (Intercept)                                                                       0.1155116
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg                                2.1340569
# AF_data$Acute_Infarction_brain                                                    6.5425248
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg:AF_data$Acute_Infarction_brain 2.4893317

# no interaction b/w Amicar and infarct

# so our main main effects model is:
log.out_AF_infarct = glm(AF_data$Electrographic_Seizure ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg + AF_data$Acute_Infarction_brain, family = binomial (link=logit))

# get stats
exp(log.out_AF_infarct$coeff)
summary(log.out_AF_infarct)$coeff
exp(confint(log.out_AF_infarct))
# > exp(log.out_AF_infarct$coeff)
#                                        (Intercept)
#                                         0.05099238
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg
#                                         1.62169132
#                     AF_data$Acute_Infarction_brain
#                                         2.58237754
# > summary(log.out_AF_infarct)$coeff
#                                                      Estimate Std. Error
# (Intercept)                                        -2.9760791  0.3505037
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg  0.4834696  0.1439412
# AF_data$Acute_Infarction_brain                      0.9487105  0.4720719
#                                                      z value     Pr(>|z|)
# (Intercept)                                        -8.490863 2.051088e-17
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg  3.358800 7.828163e-04
# AF_data$Acute_Infarction_brain                      2.009674 4.446572e-02
# > exp(confint(log.out_AF_infarct))
# Waiting for profiling to be done...
#                                                         2.5 %    97.5 %
# (Intercept)                                        0.02408046 0.0961278
# AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg 1.23099037 2.1741302
# AF_data$Acute_Infarction_brain                     1.00067950 6.4647881

## mediation analysis
#         c
# ACA --------> Seizure

#       Stroke
#     ^       \
#   a/         \b
#   /      c'   \
# ACA ---------Seizure

install.packages("mediation")
library(mediation)

c_path = glm(Electrographic_Seizure ~ Antifibrinolytic, data = AF_data, family = binomial (link=logit))
exp(c_path$coeff)
summary(c_path)$coeff
exp(confint(c_path))

a_path = glm(Acute_Infarction_brain ~ Antifibrinolytic, data = AF_data, family = binomial (link=logit))
exp(a_path$coeff)
summary(a_path)$coeff
exp(confint(a_path))

b_path = glm(Electrographic_Seizure ~ Acute_Infarction_brain + Antifibrinolytic, data = AF_data, family = binomial (link=logit))
exp(b_path$coeff)
summary(b_path)$coeff
exp(confint(b_path))

m = mediate(a_path, b_path, sims = 5000, treat = 'Antifibrinolytic', mediator = 'Acute_Infarction_brain', boot = TRUE)
summary(m)
# Causal Mediation Analysis

# Nonparametric Bootstrap Confidence Intervals with the Percentile Method

#                          Estimate 95% CI Lower 95% CI Upper p-value
# ACME (control)            0.01955      0.00133         0.05   0.022 *
# ACME (treated)            0.03364      0.00273         0.07   0.022 *
# ADE (control)             0.06646     -0.00966         0.14   0.086 .
# ADE (treated)             0.08055     -0.01210         0.17   0.086 .
# Total Effect              0.10010      0.01228         0.18   0.022 *
# Prop. Mediated (control)  0.19534      0.00275         1.17   0.042 *
# Prop. Mediated (treated)  0.33606      0.01200         1.12   0.042 *
# ACME (average)            0.02660      0.00218         0.06   0.022 *
# ADE (average)             0.07351     -0.01072         0.15   0.086 .
# Prop. Mediated (average)  0.26570      0.00736         1.15   0.042 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Sample Size Used: 224

########### Run a model of seizure as a function of stroke, ACA dose, age, sex race
log.out_AF_infarct = glm(AF_data$Electrographic_Seizure ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg + AF_data$Acute_Infarction_brain + AF_data$Age + AF_data$Sex + AF_data$Race, family = binomial (link=logit))

# get stats
exp(log.out_AF_infarct$coeff)
summary(log.out_AF_infarct)$coeff
exp(confint(log.out_AF_infarct))


#######################################################################################################
### Multivariate analyses (outcome: any cortical hyperexcitability)
#######################################################################################################

#log.out = glm(AF_data$Markers_of_cortical_hyperexcitability ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg + AF_data$Acute_Infarction_brain + AF_data$Total_CPB_Time_mins + AF_data$Age + AF_data$Sex + AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg*AF_data$Acute_Infarction_brain, family = binomial (link=logit))

log.out = glm(AF_data$Markers_of_cortical_hyperexcitability ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg + AF_data$Acute_Infarction_brain + AF_data$Total_CPB_Time_mins + AF_data$Age + AF_data$Sex + AF_data$Race + AF_data$Acute_Hemorrhage_brain, family = binomial (link=logit))

# get stats
exp(log.out$coeff)
summary(log.out)$coeff
exp(confint(log.out))

#######################################################################################################
### Multivariate analyses (outcome: Acute_Infarction_brain)
#######################################################################################################

log.out = glm(AF_data$Acute_Infarction_brain ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg + AF_data$Total_CPB_Time_mins + AF_data$Age + AF_data$Stroke_Prior_to_Admit, family = binomial (link=logit))

# get stats
exp(log.out$coeff)
summary(log.out)$coeff
exp(confint(log.out))

#######################################################################################################
### Multivariate analyses (outcome: generalized markers of cortical hyperexcitability)
#######################################################################################################

log.out = glm(AF_data$Generalized_Marker_Cortical_Hyperexcitability ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg + AF_data$Total_CPB_Time_mins + AF_data$Sex + AF_data$Age, family = binomial (link=logit))

# get stats
exp(log.out$coeff)
summary(log.out)$coeff
exp(confint(log.out))

#######################################################################################################
### Multivariate analyses (outcome: lateralized markers of cortical hyperexcitability)
#######################################################################################################

log.out = glm(AF_data$Lateralized_Marker_Cortical_Hyperexcitability ~ AF_data$Weight_adjusted_total_AF_dose_100mg_per_kg + AF_data$Total_CPB_Time_mins + AF_data$Sex + AF_data$Age + AF_data$Acute_Infarction_brain, family = binomial (link=logit))

# get stats
exp(log.out$coeff)
summary(log.out)$coeff
exp(confint(log.out))

################ DRAFTS: 

log.out_AF_infarct = glm(AF_data$Electrographic_Seizure ~ AF_data$Antifibrinolytic + AF_data$Acute_Infarction_brain + AF_data$CPB_yes_or_no, family = binomial (link=logit))

# get stats
exp(log.out_AF_infarct$coeff)
summary(log.out_AF_infarct)$coeff
exp(confint(log.out_AF_infarct))
