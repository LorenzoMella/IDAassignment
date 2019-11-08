### IDA - Assignment
### Code for Question 2

## Preliminary configuration
library(tidyverse)

printf = function(frmt, ...) cat(sprintf(frmt, ...))

load("~/rprojects/ida/databp.Rdata")

# Hack to make a one-level tibble
bp_tb = tibble(logdose = databp$logdose,
               bloodp = databp$bloodp,
               recovtime = databp$recovtime,
               R = databp$R)

## (a) Complete Case Analysis

# Mean recovery time and std error of the mean
cca_stats = bp_tb %>% summarize(mean_rec_time = mean(recovtime, na.rm = TRUE),
                    mean_std_rec_time = sd(recovtime, na.rm = TRUE) / sqrt(sum(!is.na(recovtime))),
                    total_observed = sum(!is.na(recovtime)))

print(cca_stats)

# Pearson correlations
cca_cors = bp_tb %>%
  summarize(dose_rectime_cor = cor(logdose, recovtime, use = "complete.obs", method = "pearson"),
            bloodp_rectime_cor = cor(bloodp, recovtime, use = "complete.obs", method = "pearson"))

print(cca_cors)

# Elementary way
cor(bp_tb$logdose[bp_tb$R == 1], bp_tb$recovtime[bp_tb$R == 1], method = "pearson")
cor(bp_tb$bloodp[bp_tb$R == 1], bp_tb$recovtime[bp_tb$R == 1], method = "pearson")


## (b) Mean Imputation

# Create a new tibble with the imputed data
mi_bp_tb = bp_tb %>%
  # Replace the NA recovtime entries with their mean value (taken the CCA results in point (a))
  mutate_at(.vars = "recovtime", .funs = function(x) ifelse(is.na(x), cca_stats$mean_rec_time, x)) %>%
  # We don't need the indicator anymore
  select(-R) 

# Mean recovery time and std error of the mean
mi_stats = mi_bp_tb %>% summarize(mean_rec_time = mean(recovtime),
                                  mean_std_rec_time = sd(recovtime) / sqrt(sum(!is.na(recovtime))),
                                  total_observed = sum(!is.na(recovtime)))

print(mi_stats)

# Compute and print Pearson correlations
mi_cors = mi_bp_tb %>% summarize(dose_rectime_cor = cor(logdose, recovtime, method = "pearson"),
                                 bloodp_rectime_cor = cor(bloodp, recovtime, method = "pearson"))

print(mi_cors)

# Elementary way (to double-check)
cor(mi_bp_tb$logdose, mi_bp_tb$recovtime, method = "pearson")
cor(mi_bp_tb$bloodp, mi_bp_tb$recovtime, method = "pearson")

## (c) Mean Regression Imputation

# Perform linear regression
regression_data = lm(formula = recovtime ~ logdose + bloodp, data = bp_tb)
# Print formatted summary
summary(regression_data)

# Impute values as predictions using the linear model
# (all entries are predicted just for alignment purposes)
linreg_pred_recovtime = bp_tb %>%
  predict(regression_data, newdata = .)

linreg_bp_tb = bp_tb %>%
  # Add column (temporarily) with predicted data to the tibble
  add_column(pred_rec_time = linreg_pred_recovtime) %>%
  # Replace NA entries in recovtime with the predicted values
  mutate(recovtime = ifelse(R == 1, recovtime, pred_rec_time)) %>%
  # Remove the temporary column and the indicator
  select(-pred_rec_time, -R)

# Mean recovery time and std error of the mean
linreg_stats = linreg_bp_tb %>% summarize(mean_rec_time = mean(recovtime),
                                  mean_std_rec_time = sd(recovtime) / sqrt(sum(!is.na(recovtime))),
                                  total_observed = sum(!is.na(recovtime)))

print(linreg_stats)

# Pearson correlations
linreg_cors = linreg_bp_tb %>% summarize(dose_rectime_cor = cor(logdose, recovtime, method = "pearson"),
                                 bloodp_rectime_cor = cor(bloodp, recovtime, method = "pearson"))

print(linreg_cors)

# Elementary way
cor(linreg_bp_tb$logdose, linreg_bp_tb$recovtime, method = "pearson")
cor(linreg_bp_tb$bloodp, linreg_bp_tb$recovtime, method = "pearson")


## (d) Stochastic Regression Imputation

# Regression data inherited from point (c)
linreg_residual_sd = summary(regression_data)$sigma

# Impute values as in (c) but adding Gaussian noise
sri_pred_recovtime = linreg_pred_recovtime +
  rnorm(nrow(bp_tb), mean = 0, sd = linreg_residual_sd)

sri_bp_tb = bp_tb %>%
  add_column(pred_rec_time = sri_pred_recovtime) %>%
  mutate(recovtime = ifelse(R == 1, recovtime, pred_rec_time)) %>%
  select(-pred_rec_time, -R)

# Mean recovery time and std error of the mean
sri_stats = sri_bp_tb %>% summarize(mean_rec_time = mean(recovtime),
                                  mean_std_rec_time = sd(recovtime) / sqrt(sum(!is.na(recovtime))),
                                  total_observed = sum(!is.na(recovtime)))

print(sri_stats)

# Pearson correlations
sri_cors = sri_bp_tb %>% summarize(dose_rectime_cor = cor(logdose, recovtime, method = "pearson"),
                                 bloodp_rectime_cor = cor(bloodp, recovtime, method = "pearson"))

print(sri_cors)

# Elementary way
cor(sri_bp_tb$logdose, sri_bp_tb$recovtime, method = "pearson")
cor(sri_bp_tb$bloodp, sri_bp_tb$recovtime, method = "pearson")

# Diagnostics

# Linear plot of studentized residuals
sri_studentized_res = rstudent(regression_data)
plot(regression_data$fitted.values, sri_studentized_res,
     xlab = "Fitted values (LR)", ylab = "Studentised Residuals")

# QQ-Plot of studentized residuals (against the theoretical Student's t distribution)
sri_standardized_res = rstandard(regression_data)
sri_standardized_res_no_outliers = sri_standardized_res[sri_standardized_res < 2]
                                            
qqnorm(sri_standardized_res_no_outliers)
qqline(sri_standardized_res, col = 2)


## (e) Predictive Mean Matching
