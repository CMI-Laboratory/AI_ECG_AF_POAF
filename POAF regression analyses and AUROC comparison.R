setwd("G:/Workspace/Codes/Python/postopafib")

df<-read.csv('20250117_df_total5.csv')

df$prediction_score <- df$prediction_score * 10
df$agecat <- with(df, ifelse(age < 60, 1,
                             ifelse(age >= 60 & age <= 70, 2,
                                    ifelse(age > 70 & age <= 80, 3, 4))))
df$eGFR_below15 <- ifelse(df$eGFR < 15, 1, 0)
df$LVEF_below30 <- ifelse(df$EF < 30, 1, 0)
mva_full_glm  <- glm(new.AF  ~ as.factor(agecat) +as.factor(COPD) +as.factor(eGFR_below15) +as.factor(emergency) + as.factor(ECMOVADIABP_label) + 
                       as.factor(LVEF_below30) + as.factor(valve_surgery)
                     +prediction_score
                     , data=df,family="binomial",na.action = na.exclude)
summary(mva_full_glm)

conf_intervals <- confint(mva_full_glm, level = 0.95)

# Exponentiate coefficients and confidence intervals to get odds ratios
odds_ratios <- exp(coef(mva_full_glm))
odds_ratios_ci <- exp(conf_intervals)

# Combine results into a data frame
results <- data.frame(
  Variable = names(odds_ratios),
  Odds_Ratio = odds_ratios,
  CI_Lower = odds_ratios_ci[, 1],
  CI_Upper = odds_ratios_ci[, 2]
)

# Print the results
print(results)

df<-read.csv('20250117_df_total5.csv')

df$prediction_score <- df$prediction_score * 10
df$sex <- relevel(as.factor(df$sex), ref = "1")
df$WBC <- df$WBC * 0.01
df$Plt...1000. <- df$Plt...1000. * 0.01
head(df,1)
mva_full_glm  <- glm(new.AF  ~ 
                       age + height + weight + BSA + SBP + DBP + WBC + Hb + Plt...1000. +
                       Na + K + BUN + Cr + eGFR + OT + PT + Albumin + Glu + PT_sec +
                       PTT + EF + as.factor(valve_surgery) + 
                       as.factor(sex) + as.factor(emergency) + as.factor(HTN) + 
                       as.factor(CKD) + as.factor(DM) + as.factor(old.CVA) +
                       as.factor(liver.cirrhosis) + as.factor(CHF) +
                       as.factor(old.MI) + as.factor(COPD) +
                       as.factor(recentMI..3mon.) + as.factor(acuteMI..1WK.) + as.factor(BB) +
                       as.factor(CCB) + as.factor(ACEi) + as.factor(ARB) + as.factor(Statin) +
                       as.factor(PO.Nitrate) + as.factor(diuretics) + as.factor(Warfarin) +
                       as.factor(heparinization) + as.factor(IV.nitrate) + as.factor(NOAC) +
                       as.factor(ECMOVADIABP_label)   +prediction_score
                     , data=df,family="binomial")

mva_full_glm_step <- step(mva_full_glm,direction='backward')
summary(mva_full_glm_step)

conf_intervals <- confint(mva_full_glm_step, level = 0.95)

# Exponentiate coefficients and confidence intervals to get odds ratios
odds_ratios <- exp(coef(mva_full_glm_step))
odds_ratios_ci <- exp(conf_intervals)

# Combine results into a data frame
results <- data.frame(
  Variable = names(odds_ratios),
  Odds_Ratio = odds_ratios,
  CI_Lower = odds_ratios_ci[, 1],
  CI_Upper = odds_ratios_ci[, 2]
)

# Print the results
print(results)

df<-read.csv('20250119_df_total5.csv')
df$prediction_score <- df$prediction_score * 10
df$agecat <- with(df, ifelse(age < 60, 1,
                             ifelse(age >= 60 & age <= 70, 2,
                                    ifelse(age > 70 & age <= 80, 3, 4))))
df$eGFR_below15 <- ifelse(df$eGFR < 15, 1, 0)
df$LVEF_below30 <- ifelse(df$EF < 30, 1, 0)

df_M <- subset(df, sex == 0)
df_F <- subset(df, sex == 1)

mva_full_glm_M  <- glm(new.AF  ~ as.factor(agecat) +as.factor(COPD) +as.factor(eGFR_below15) +as.factor(emergency) + as.factor(ECMOVADIABP_label) + 
                       as.factor(LVEF_below30) + as.factor(valve_surgery)
                     +prediction_score
                     , data=df_M,family="binomial",na.action = na.exclude)

mva_full_glm_F  <- glm(new.AF  ~ as.factor(agecat) +as.factor(COPD) +as.factor(eGFR_below15) +as.factor(emergency) + as.factor(ECMOVADIABP_label) + 
                         as.factor(LVEF_below30) + as.factor(valve_surgery)
                       +prediction_score
                       , data=df_F,family="binomial",na.action = na.exclude)

summary(mva_full_glm_M)
summary(mva_full_glm_F)

# Extract the coefficients and standard errors for 'AI_prediction_average' from both models
coef_M <- summary(mva_full_glm_M)$coefficients["prediction_score", "Estimate"]
se_M <- summary(mva_full_glm_M)$coefficients["prediction_score", "Std. Error"]

coef_F <- summary(mva_full_glm_F)$coefficients["prediction_score", "Estimate"]
se_F <- summary(mva_full_glm_F)$coefficients["prediction_score", "Std. Error"]

exp(coef_M)
exp(coef_F)

confint(mva_full_glm_M, level=0.95)
confint(mva_full_glm_F, level=0.95)

# Apply Z-test to compare the log-hazard ratios
z_value_sex <- (coef_M - coef_F) / sqrt(se_M^2 + se_F^2)

# Calculate the p-value for the Z-test
p_value_sex <- 2 * pnorm(-abs(z_value_sex))

# Print Z-value and p-value for sex subgroup comparison
cat("Z-value for sex:", z_value_sex, "\n")
cat("P-value for sex:", p_value_sex, "\n")


# Create year_difference subgroup variable (0 for < 65, 1 for >= 65)
df$year_group <- ifelse(df$age < 65, 0, 1)

# Subset data by year group
df_year0 <- subset(df, year_group == 0)  # year_difference < 65
df_year1 <- subset(df, year_group == 1)  # year_difference >= 65

# Fit model for age < 65
mva_full_glm_year0  <- glm(new.AF  ~ as.factor(agecat) +as.factor(COPD) +as.factor(eGFR_below15) +as.factor(emergency) + as.factor(ECMOVADIABP_label) + 
                         as.factor(LVEF_below30) + as.factor(valve_surgery)
                       +prediction_score
                       , data=df_year0,family="binomial",na.action = na.exclude)

# Fit Cox model for age >=65
mva_full_glm_year1  <- glm(new.AF  ~ as.factor(agecat) +as.factor(COPD) +as.factor(eGFR_below15) +as.factor(emergency) + as.factor(ECMOVADIABP_label) + 
                         as.factor(LVEF_below30) + as.factor(valve_surgery)
                       +prediction_score
                       , data=df_year1,family="binomial",na.action = na.exclude)


summary(mva_full_glm_year0)
summary(mva_full_glm_year1)

# Extract the coefficients and standard errors for 'AI_prediction_average' from both models
coef_year0 <- summary(mva_full_glm_year0)$coefficients["prediction_score", "Estimate"]
se_year0 <- summary(mva_full_glm_year0)$coefficients["prediction_score", "Std. Error"]

coef_year1 <- summary(mva_full_glm_year1)$coefficients["prediction_score", "Estimate"]
se_year1 <- summary(mva_full_glm_year1)$coefficients["prediction_score", "Std. Error"]

exp(coef_year0)
exp(coef_year1)

confint(mva_full_glm_year0, level=0.95)
confint(mva_full_glm_year1, level=0.95)

# Apply Z-test to compare the log-hazard ratios
z_value_year <- (coef_year0 - coef_year1) / sqrt(se_year0^2 + se_year1^2)

# Calculate the p-value for the Z-test
p_value_year <- 2 * pnorm(-abs(z_value_year))

# Print Z-value and p-value for year_difference subgroup comparison
cat("Z-value for year_difference:", z_value_year, "\n")
cat("P-value for year_difference:", p_value_year, "\n")


### Subgroup analysis by surgery type (valve vs. non-valve)

df_valve0 <- subset(df, valve_surgery == 0) 
df_valve1 <- subset(df, valve_surgery == 1) 


mva_full_glm_valve0  <- glm(new.AF  ~ as.factor(agecat) +as.factor(COPD) +as.factor(eGFR_below15) +as.factor(emergency) + as.factor(ECMOVADIABP_label) + 
                             as.factor(LVEF_below30) 
                           +prediction_score
                           , data=df_valve0,family="binomial",na.action = na.exclude)

mva_full_glm_valve1  <- glm(new.AF  ~ as.factor(agecat) +as.factor(COPD) +as.factor(eGFR_below15) +as.factor(emergency) + as.factor(ECMOVADIABP_label) + 
                             as.factor(LVEF_below30) 
                           +prediction_score
                           , data=df_valve1,family="binomial",na.action = na.exclude)


summary(mva_full_glm_valve0)
summary(mva_full_glm_valve1)

# Extract the coefficients and standard errors for 'AI_prediction_average' from both models
coef_valve0 <- summary(mva_full_glm_valve0)$coefficients["prediction_score", "Estimate"]
se_valve0 <- summary(mva_full_glm_valve0)$coefficients["prediction_score", "Std. Error"]

coef_valve1 <- summary(mva_full_glm_valve1)$coefficients["prediction_score", "Estimate"]
se_valve1 <- summary(mva_full_glm_valve1)$coefficients["prediction_score", "Std. Error"]

exp(coef_valve0)
exp(coef_valve1)

confint(mva_full_glm_valve0, level=0.95)
confint(mva_full_glm_valve1, level=0.95)

# Apply Z-test to compare the log-hazard ratios
z_value_year <- (coef_year0 - coef_year1) / sqrt(se_year0^2 + se_year1^2)

# Calculate the p-value for the Z-test
p_value_year <- 2 * pnorm(-abs(z_value_year))

# Print Z-value and p-value for year_difference subgroup comparison
cat("Z-value for year_difference:", z_value_year, "\n")
cat("P-value for year_difference:", p_value_year, "\n")


df<-read.csv('20250117_df_total5.csv')
df$prediction_score <- df$prediction_score * 10
df$agecat <- with(df, ifelse(age < 60, 1,
                             ifelse(age >= 60 & age <= 70, 2,
                                    ifelse(age > 70 & age <= 80, 3, 4))))
df$eGFR_below15 <- ifelse(df$eGFR < 15, 1, 0)
df$LVEF_below30 <- ifelse(df$EF < 30, 1, 0)

df_M <- subset(df, sex == 0)
df_F <- subset(df, sex == 1)


mva_full_glm_M  <- glm(new.AF  ~ 
                       age + weight + DBP + WBC + Hb + Plt...1000. +
                       K + BUN + Cr + eGFR +PT + Albumin + PT_sec +
                       PTT + as.factor(valve_surgery) + 
                        as.factor(HTN) + 
                       as.factor(CKD) + as.factor(old.CVA) +
                       as.factor(CHF) +
                       as.factor(old.MI) + as.factor(COPD) +
                       as.factor(recentMI..3mon.) + as.factor(acuteMI..1WK.) + as.factor(BB) +
                       as.factor(CCB) + as.factor(ARB)  +
                       as.factor(Warfarin) +
                       as.factor(heparinization) + as.factor(IV.nitrate) + as.factor(NOAC) +
                       prediction_score
                     , data=df_M,family="binomial")

mva_full_glm_F  <- glm(new.AF  ~ 
                         age + weight + DBP + WBC + Hb + Plt...1000. +
                         K + BUN + Cr + eGFR +PT + Albumin + PT_sec +
                         PTT + as.factor(valve_surgery) + 
                          as.factor(HTN) + 
                         as.factor(CKD) + as.factor(old.CVA) +
                         as.factor(CHF) +
                         as.factor(old.MI) + as.factor(COPD) +
                         as.factor(recentMI..3mon.) + as.factor(acuteMI..1WK.) + as.factor(BB) +
                         as.factor(CCB) + as.factor(ARB)  +
                         as.factor(Warfarin) +
                         as.factor(heparinization) + as.factor(IV.nitrate) + as.factor(NOAC) +
                         prediction_score
                       , data=df_F,family="binomial")

summary(mva_full_glm_M)
summary(mva_full_glm_F)

# Extract the coefficients and standard errors for 'AI_prediction_average' from both models
coef_M <- summary(mva_full_glm_M)$coefficients["prediction_score", "Estimate"]
se_M <- summary(mva_full_glm_M)$coefficients["prediction_score", "Std. Error"]

coef_F <- summary(mva_full_glm_F)$coefficients["prediction_score", "Estimate"]
se_F <- summary(mva_full_glm_F)$coefficients["prediction_score", "Std. Error"]

exp(coef_M)
exp(coef_F)

confint(mva_full_glm_M, level=0.95)
confint(mva_full_glm_F, level=0.95)

# Apply Z-test to compare the log-hazard ratios
z_value_sex <- (coef_M - coef_F) / sqrt(se_M^2 + se_F^2)

# Calculate the p-value for the Z-test
p_value_sex <- 2 * pnorm(-abs(z_value_sex))

# Print Z-value and p-value for sex subgroup comparison
cat("Z-value for sex:", z_value_sex, "\n")
cat("P-value for sex:", p_value_sex, "\n")

### Subgroup analysis by year_difference (< 65 years vs. >= 65 years)

# Create year_difference subgroup variable (0 for < 65, 1 for >= 65)
df$year_group <- ifelse(df$age < 65, 0, 1)

# Subset data by year group
df_year0 <- subset(df, year_group == 0)  # year_difference < 65
df_year1 <- subset(df, year_group == 1)  # year_difference >= 65


mva_full_glm_year0  <- glm(new.AF  ~ 
                         age + weight + DBP + WBC + Hb + Plt...1000. +
                         K + BUN + Cr + eGFR +PT + Albumin + PT_sec +
                         PTT + as.factor(valve_surgery) + 
                        as.factor(sex) + as.factor(HTN) + 
                         as.factor(CKD) + as.factor(old.CVA) +
                         as.factor(CHF) +
                         as.factor(old.MI) + as.factor(COPD) +
                         as.factor(recentMI..3mon.) + as.factor(acuteMI..1WK.) + as.factor(BB) +
                         as.factor(CCB) + as.factor(ARB)  +
                         as.factor(Warfarin) +
                         as.factor(heparinization) + as.factor(IV.nitrate) + as.factor(NOAC) +
                         prediction_score
                       , data=df_year0,family="binomial")

mva_full_glm_year1  <- glm(new.AF  ~ 
                         age + weight + DBP + WBC + Hb + Plt...1000. +
                         K + BUN + Cr + eGFR +PT + Albumin + PT_sec +
                         PTT + as.factor(valve_surgery) + 
                         as.factor(sex) + as.factor(HTN) + 
                         as.factor(CKD) + as.factor(old.CVA) +
                         as.factor(CHF) +
                         as.factor(old.MI) + as.factor(COPD) +
                         as.factor(recentMI..3mon.) + as.factor(acuteMI..1WK.) + as.factor(BB) +
                         as.factor(CCB) + as.factor(ARB)  +
                         as.factor(Warfarin) +
                         as.factor(heparinization) + as.factor(IV.nitrate) + as.factor(NOAC) +
                         prediction_score
                       , data=df_year1,family="binomial")


summary(mva_full_glm_year0)
summary(mva_full_glm_year1)

# Extract the coefficients and standard errors for 'AI_prediction_average' from both models
coef_year0 <- summary(mva_full_glm_year0)$coefficients["prediction_score", "Estimate"]
se_year0 <- summary(mva_full_glm_year0)$coefficients["prediction_score", "Std. Error"]

coef_year1 <- summary(mva_full_glm_year1)$coefficients["prediction_score", "Estimate"]
se_year1 <- summary(mva_full_glm_year1)$coefficients["prediction_score", "Std. Error"]

exp(coef_year0)
exp(coef_year1)

confint(mva_full_glm_year0, level=0.95)
confint(mva_full_glm_year1, level=0.95)

# Apply Z-test to compare the log-hazard ratios
z_value_year <- (coef_year0 - coef_year1) / sqrt(se_year0^2 + se_year1^2)

# Calculate the p-value for the Z-test
p_value_year <- 2 * pnorm(-abs(z_value_year))

# Print Z-value and p-value for year_difference subgroup comparison
cat("Z-value for year_difference:", z_value_year, "\n")
cat("P-value for year_difference:", p_value_year, "\n")


### Subgroup analysis by surgery type (valve vs. non-valve)

df_valve0 <- subset(df, valve_surgery == 0) 
df_valve1 <- subset(df, valve_surgery == 1) 


mva_full_glm_valve0  <- glm(new.AF  ~ 
                             age + weight + DBP + WBC + Hb + Plt...1000. +
                             K + BUN + Cr + eGFR +PT + Albumin + PT_sec +
                             PTT + 
                             as.factor(sex) + as.factor(HTN) + 
                             as.factor(CKD) + as.factor(old.CVA) +
                             as.factor(CHF) +
                             as.factor(old.MI) + as.factor(COPD) +
                             as.factor(recentMI..3mon.) + as.factor(acuteMI..1WK.) + as.factor(BB) +
                             as.factor(CCB) + as.factor(ARB)  +
                             as.factor(Warfarin) +
                             as.factor(heparinization) + as.factor(IV.nitrate) + as.factor(NOAC) +
                             prediction_score
                           , data=df_valve0,family="binomial")

mva_full_glm_valve1  <- glm(new.AF  ~ 
                             age + weight + DBP + WBC + Hb + Plt...1000. +
                             K + BUN + Cr + eGFR +PT + Albumin + PT_sec +
                             PTT +  
                             as.factor(sex) + as.factor(HTN) + 
                             as.factor(CKD) + as.factor(old.CVA) +
                             as.factor(CHF) +
                             as.factor(old.MI) + as.factor(COPD) +
                             as.factor(recentMI..3mon.) + as.factor(acuteMI..1WK.) + as.factor(BB) +
                             as.factor(CCB) + as.factor(ARB)  +
                             as.factor(Warfarin) +
                             as.factor(heparinization) + as.factor(IV.nitrate) + as.factor(NOAC) +
                             prediction_score
                           , data=df_valve1,family="binomial")

summary(mva_full_glm_valve0)
summary(mva_full_glm_valve1)

# Extract the coefficients and standard errors for 'AI_prediction_average' from both models
coef_valve0 <- summary(mva_full_glm_valve0)$coefficients["prediction_score", "Estimate"]
se_valve0 <- summary(mva_full_glm_valve0)$coefficients["prediction_score", "Std. Error"]

coef_valve1 <- summary(mva_full_glm_valve1)$coefficients["prediction_score", "Estimate"]
se_valve1 <- summary(mva_full_glm_valve1)$coefficients["prediction_score", "Std. Error"]

exp(coef_valve0)
exp(coef_valve1)

confint(mva_full_glm_valve0, level=0.95)
confint(mva_full_glm_valve1, level=0.95)

# Apply Z-test to compare the log-hazard ratios
z_value_year <- (coef_year0 - coef_year1) / sqrt(se_year0^2 + se_year1^2)

# Calculate the p-value for the Z-test
p_value_year <- 2 * pnorm(-abs(z_value_year))

# Print Z-value and p-value for year_difference subgroup comparison
cat("Z-value for year_difference:", z_value_year, "\n")
cat("P-value for year_difference:", p_value_year, "\n")



setwd("G:/Workspace/Codes/Python/postopafib")

df<-read.csv('20250117_R_comparison.csv')

library(pROC)
roc0<-roc(df$new.AF,df$POAF.score)
roc1<-roc(df$new.AF,df$prediction_score)
roc2<-roc(df$new.AF,df$POAF_AI.score)
roc3<-roc(df$new.AF,df$NEWscore)

roc.test(roc0,roc1)
roc.test(roc1,roc2)
roc.test(roc2,roc3)
roc.test(roc1,roc3)

