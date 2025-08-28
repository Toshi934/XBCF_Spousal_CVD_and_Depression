####### Heterogeneous Treatment Effect Analysis Code #######
### Last Update: 2025/08/28 ###
rm(list = ls(all.names = TRUE))

#########################################################################
### Install packages
library(tidyverse)
library(gtsummary)
library(XBCF) # XBCF package can be installed from: https://github.com/socket778/XBCF
### Load data
df <- read.csv("spousal_CVD.csv")
#########################################################################

### Step1: data preparation ###
# List of variables
continuous <- c("age", "age.dep", "BMI", "total.cholesterol", "sysbp", "diabp", "HDL", "fbg", "eGFR")
categorical <- c("sex", "income", "hypertension", "CVD", "diabetes", "hypertension.s", "depression.s", "diabetes.s", "smoking", "drinking", "physical.activity", "antihypertensive")
exposure <- c("CVD.s")
outcome <- c("depression")

covariates <- c(continuous, categorical)

# Adjust the variable type
df[categorical] <- lapply(df[categorical], factor)

# Split data into 1) training sample and 2) test sample
set.seed(1)
df_train <- df[sample(nrow(df), 0.5 * nrow(df)),]
df_test <- df[!df$ID %in% df_train$ID,]

# Estimate propensity scores for training and test samples 
ps_model_train <- glm(exposure ~ age + age.dep + BMI + total.cholesterol + sysbp + diabp + HDL + fbg + eGFR + sex + income + hypertension + CVD + diabetes + hypertension.s + depression.s + diabetes.s + smoking + drinking + physical.activity + antihypertensive, 
                data=df_train, family=binomial(link = "logit"))
df_train$PS <- stats::predict(ps_model_train, type = "response")

ps_model_test <- glm(exposure ~ age + age.dep + BMI + total.cholesterol + sysbp + diabp + HDL + fbg + eGFR + sex + income + hypertension + CVD + diabetes + hypertension.s + depression.s + diabetes.s + smoking + drinking + physical.activity + antihypertensive, 
                      data=df_test, family=binomial(link = "logit"))
df_test$PS <- stats::predict(ps_model_test, type = "response")

# Adjust the data type into numeric 
df_train[categorical] <- lapply(df_train[categorical], numeric)
df_test[categorical] <- lapply(df_test[categorical], numeric)



### Step2: run accelerated Bayesian causal forest ###
# Run XBCF
xbcf <- XBCF(y = df_train[,outcome], # Outcome 
                z = df_train[,exposure], # Exposure 
                x_con = df_train[covariates], # Covariates 
                pihat = df_train$PS, # Propensity score 
                num_sweeps = 500, # Iterations for CATE estimation 
                burnin = 300, # Burn-in 
                pcat_con = length(categorical), # The number of categorical variables 
                n_trees_con = 300, # The number of trees for prognostic score
                n_trees_mod = 300, # The number of trees for treatment effect
                random_seed = 1) 

# Estimate CATE for training and test samples
df_train$tau <- predictTaus(xbcf, df_train[covariates])
df_test$tau <- predictTaus(xbcf, df_test[covariates])



### Step3: calibration of estimated CATE ###
# 1) Calibration plot
# Classify samples into CATE quartiles 
df_test <- df_test %>% mutate(ranking = ntile(tau, 4))
ranking <- df_test$ranking

# Average difference-in-means within each ranking
# Formula y ~ 0 + ranking + ranking:w 
fmla <- paste0(outcome, "~ 0 + ranking + ranking:", exposure)
ols.ate <- lm(fmla, data = transform(df_test, ranking=factor(ranking)))
ols.ate <- coeftest(ols.ate, vcov=vcovHC(ols.ate, type = "HC2"))
interact <- which(grepl(":", rownames(ols.ate)))
ols.ate <- data.frame("ols", paste0("Q", seq(4)), ols.ate[interact, 1:2])
rownames(ols.ate) <- NULL 
colnames(ols.ate) <- c("method", "ranking", "estimate", "std.err")

# Convert data type
df_test[categorical] <- lapply(df_test[categorical], factor)

# Compute Y hat
model_outcome <- glm(depression ~ age + age.dep + BMI + total.cholesterol + sysbp + diabp + HDL + fbg + eGFR + sex + income + hypertension + CVD + diabetes + hypertension.s + depression.s + diabetes.s + smoking + drinking + physical.activity + antihypertensive, 
                     data=df_test, family=binomial(link = "logit"))
df_test$Y_hat <- stats::predict(model_out, type = "response")

# Computing AIPW scores 
tau.hat <- df_test$tau
e.hat <- df_test$PS #P[W=1|X]
m.hat <- df_test$Y_hat #E[Y|X]

# Estimating mu.hat(X, 1) and mu.hat(X, 0) for observations
mu.hat.0 <- m.hat - e.hat * tau.hat
mu.hat.1 <- m.hat + (1 - e.hat) * tau.hat

# AIPW scores
aipw.scores <- tau.hat + df_test$CVD.s / e.hat * (df_test$depression - mu.hat.1) - (1-df_test$CVD.s) / (1-e.hat) * (df_test$depression - mu.hat.0)
ols <- lm(aipw.scores ~ 0 + factor(ranking))
forest.ate <- data.frame("aipw", paste0("Q", seq(4)), coeftest(ols, vcov = vcovHC(ols, "HC2"))[,1:2])
colnames(forest.ate) <- c("method", "ranking", "estimate", "std.err")
rownames(forest.ate) <- NULL

# Summarize results from OLS and AIPW 
res <- rbind(forest.ate, ols.ate)

# Plot the calibration results
ggplot(res) + 
  aes(x = ranking, y = estimate, group = method, color = method) + 
  geom_point(position = position_dodge(0.2)) + 
  geom_errorbar(aes(ymin = estimate -2 * std.err, ymax = estimate +2 * std.err), width = .2, position = position_dodge(0.2)) + 
  ylab("") + xlab("") + 
  ggtitle("Average CATE within each ranking (as defined by predicted CATE)") + 
  theme_minimal() + 
  theme(legend.position = "bottom", legend.title = element_blank()) + 
  ggsci::scale_color_jama()
gc()
gc()



# 2) Estimate risk difference and odds ratio of index individuals' depression by exposure and CATE ranking
# Linear regression for risk difference
model_rd <- lm(depression ~ CVD.s + CVD.s*ranking + age + age.dep + BMI + total.cholesterol + sysbp + diabp + HDL + fbg + eGFR + sex + income + hypertension + CVD + diabetes + hypertension.s + depression.s + diabetes.s + smoking + drinking + physical.activity + antihypertensive, 
                     data=df_test)

# Logistic regression for odds ratio
model_or <- lm(as.factor(depression) ~ CVD.s + CVD.s*ranking + age + age.dep + BMI + total.cholesterol + sysbp + diabp + HDL + fbg + eGFR + sex + income + hypertension + CVD + diabetes + hypertension.s + depression.s + diabetes.s + smoking + drinking + physical.activity + antihypertensive, 
               data=df_test, family=binomial(link = "logit"))

# Summarize results
summary(model_rd)
summary(model_or)



# 3) Comparisons of CATE subgroups
# Create a table to compare characteristics of each CATE subgroup 
df_test[c(categorical, "ranking")] %>% 
  tbl_summary(statistic = list(all_continuous() ~ "{mean} ({sd})"),
                          digits = all_continuous() ~ 2,
                          by = ranking)
  


