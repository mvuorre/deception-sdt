##############################################################################
# Title: Stop aggregating veracity judgements! The benefits and necessity of #
# proper data analyses; Decepticon II Conference                             #
# Author: Zloteanu, M.                                                       #
# Date: 18 Oct. 2022                                                         #
##############################################################################

# below is an example of how to analyse binary (lie/truth) data using a probit
# mixed effects model, which roughly translate to a Signal Detection Theory 
# (SDT) model without the need to re-code or compute new variables

#### setup ===================================================================
rm(list=ls(all=TRUE)) #clean workspace environment
options(mc.cores = parallel::detectCores()) #use all CPU cores
set.seed(123) #set seed for session
#set appropriate bayesian contrasts for factors with more than 2 levels.
options(contrasts = c('contr.equalprior', 'contr.poly')) 

#### Packages ================================================================
library(tidyverse) #suite of packages for handling data (contains ggplot2)
library(tidybayes) #suite of packages for brms objects
library(brms) #model building and analysis [main package]
library(emmeans) #for calculating post hoc comparisons
library(sjPlot) #useful plotting function
library(bayestestR) #clean and complete model output
library(performance) #built-in model fit tools
library(report) #automated results reporting

# note the below brms code assumes you installed cmdstanr to speed up the 
# compiling, see here https://mc-stan.org/cmdstanr/

#### Load Data ===============================================================
# run below code to import data file using easy-to-use pop-up window:
sdt.data <- read.csv(file.choose())

# if you don't have any data, you can use this package to test things out:
# library(sdtalt), and load data(confcontr)

# this script assumes that 'Participant' is your participant variable,
# 'Condition' is your between-subjects grouping variable,
# 'LieType' is your within-subjects variable
# 'Stimuli' is your trial variable
# 'Veracity' should reflect the video veracity; coded as truth=1 and lie=0
# 'Answer' should be your DV, as each answer participants give; truth=1, lie=0

# setup data: make sure model understands which are factors
# run each line of code based on what factors your dataframe has
sdt.data$Condition <-factor(sdt.data$Condition)
sdt.data$Condition <- relevel(sdt.data$Condition, "NT") #for 3+ levels, pick baseline
sdt.data$Stimuli <-factor(sdt.data$Stimuli)
sdt.data$Participant <-factor(sdt.data$Participant)
sdt.data$LieType <-factor(sdt.data$LieType)

# read data and verify correct
head(sdt.data)
summary(sdt.data)
glimpse(sdt.data)

################ SKIP SECTION; Just for illustration #########################
## Typical (Bad) Aggregate SDT analysis                                      #
# You can ignore this code and skip to the next section for the proper model #

#the below code will create a new variable "type" which classifies each Answer
#as either a hit, miss, correct rejection, or false alarm
sdt.agg <- sdt.data %>% 
  mutate(type = "hit",
         type = ifelse(Veracity==1 & Answer==0, "miss", type),
         type = ifelse(Veracity==0 & Answer==0, "cr", type),  
         type = ifelse(Veracity==0 & Answer==1, "fa", type))

#now we add up all corresponding answers to create one row per participant, as 
#in a typical ANOVA design analysis
sdt.agg <- sdt.agg %>% 
  group_by(Participant, type) %>% 
  summarise(count = n()) %>% 
  spread(type, count)
sdt.agg #this is a new dataframe for the aggregate data, you can delete after.

#now we can compute dprime, c, and bias (aka the conventional formula for c
#as found in Macmillan & Creelman, 2005)
sdt.agg <- sdt.agg %>% 
  mutate(hr = hit / (hit+miss),
         far = fa / (fa+cr),
         zhr = qnorm(hr),
         zfa = qnorm(far),
         dprime = zhr-zfa,
         c = -zfa,
         bias = -0.5*(zhr+zfa),
         model = "no pooling")

#check the new aggregate dataframe created
glimpse(sdt.agg)
summary(sdt.agg)

#run a model that estimates the average dprime for your sample
m.agg <- brm(dprime ~ 1,
             data = sdt.agg,
             family = gaussian(), #assumes a normal distribution for dprime
             backend = "cmdstanr",
             seed = 123)

#see dprime estimate
summary(m.agg)

#plot posterior of dprime;
m.agg %>% 
  gather_draws(b_Intercept) %>%
  ggplot(aes(y = .variable, x = .value)) +
  stat_halfeye()

#you will see the Mean value is (approx.) the same as the mixed effects models, 
#but the uncertainty is much narrower, as the aggregate model does not capture 
#(read: ignores) the uncertainty due to individual participants responding to
#each trial and the influence of each video (trial) on the overall result.

#same for bias (c) (plot not shown)
m.agg.bias <- brm(c ~ 1,
             data = sdt.agg,
             family = gaussian(),
             backend = "cmdstanr",
             seed = 123)

#summary of c aggregate analysis
summary(m.agg.bias)

#############################[END OF SKIP PART]###############################

#### Model Specification and Priors ==========================================
#now we begin to build the proper Probit Mixed Effects Model

#set some weakly regularizing priors
sdt_priors <- c(prior(normal(0,1), class = b),
                prior(student_t(3, 0, 1), class = sd),
                prior(lkj(1), class = cor))

#simple model with no factors other than video Veracity to classify accuracy; 
#the model estimates dprime and c for your sample while considering random 
#effects for Participants and Stimuli, and random slope for Veracity by PP
m0 <- brm(Answer ~ 1 + Veracity + (1 + Veracity | Participant) + (1 | Stimuli),  
          data = sdt.data, 
          family = bernoulli(link = probit), #change to "link = logit" for ORs
          warmup = 2000,
          iter = 4000,
          save_pars = save_pars(all = TRUE), 
          prior = sdt_priors,
          backend = "cmdstanr",
          init = 0,
          thin = 1,
          chains = 4,
          seed = 123,
          control = list(adapt_delta = 0.9))

# assess model convergence
plot(m0) #each parameter should show nice fuzzy caterpillars
mcmc_plot(m0, type = "acf_bar") #check autocorrelation, should go towards 0 fast
pp_check(m0) #how well does the model (y-rep) fit the data (y)

# see model coefficients
summary(m0)
#better table with additional metrics to report: pd, ROPE, etc.
#for details on what each means: https://easystats.github.io/bayestestR/
describe_posterior(m0)

#you can also use the automated reporting package "report" to give you some 
#guidance on how to interpret the parameters and model, but it is a bit verbose
report(m0) #takes some time to recompile

#make a ggplot (fancier) plot of the parameter posterior distributions
m0 %>%
  posterior_samples() %>%
  clean_names() %>%
  select(starts_with("b_")) %>%
  pivot_longer(cols = everything(), 
               names_to = "variable", 
               values_to = "value") %>%
  ggplot(data = .,
         mapping = aes(y = fct_rev(variable),
                       x = value)) +
  stat_halfeye(fill = "lightblue") + 
  theme(axis.title.y = element_blank()) + 
  theme_minimal()

# model fit
m0.R2 <- bayes_R2(m0) # Bayesian pseudo-R2 (Gelman et al., 2018)
#run the code at the bottom of this script to create the function bayes_R2_mz()
m0.R2_mz <- bayes_R2_mz(m0) #corrected pseudo-R2 for probit model (McKelvey & Zavoina, 1975)

#note, Bayesian R2 should not be interpreted as comparing variance explained
#between models, as would a frequentist R2 in a regression; they are useful only 
#for predicting how much future data would be explained by your current model 
#(i.e., a model with the same exact fixed factors as this one); (Gelman et al., 2018)

m0.loo <- loo(m0) #leave-one-out cross-validation; use to select best model
m0.waic <- waic(m0) #Widely Applicable Information Criterion (WAIC); loo is better
m0.kfold <- kfold(m0, K = 3) #if loo fails, perform exact K-fold cross-validation

#fixed effects for Condition (between subjects), Veracity, and their Interaction
#random effects for Participants and for Stimuli
#random slopes based on veracity (i.e., only within-subjects effects)
m1 <- brm(Answer ~ 1 + Condition*Veracity + (1 + Veracity | Participant) + (1 | Stimuli),  
          data = sdt.data, 
          family = bernoulli(link = probit),
          warmup = 2000,
          iter = 4000,
          save_pars = save_pars(all = TRUE), 
          prior = sdt_priors,
          backend = "cmdstanr",
          init = 0,
          thin = 1,
          chains = 4,
          seed = 123,
          control = list(adapt_delta = 0.9))

# model output
conditional_effects(m1)
describe_posterior(m1)

## pairwise comparisons [post hoc tests]
#differences in c (bias) between Condition levels
m1.pw.bias <- emmeans(m1, specs = revpairwise ~Condition)
describe_posterior(m1.pw.bias)

#difference in dprime (accuracy) between Condition levels
#aka a contrast of contrasts
m1.pw.accuracy <- emmeans (m1, ~ Veracity + Condition) #make the interaction

# dprime for each Condition level
m1.pw.each <- contrast (m1.pw.accuracy, "pairwise", by = "Condition")
describe_posterior(m1.pw.each)
plot(m1.pw.each)

# #fancier plot
# ame_cond <- m1 %>%
#   epred_draws(newdata = expand_grid(Condition = c("NT", "BT", "ERT"),
#                                     Veracity = c(0.5)), re_formula = NA)
# plotconditions <- ggplot(ame_cond,
#                          aes(x = .epred,
#                              y = "Mean criterion (c) value",
#                              fill = Condition, alpha=.5)) +
#                           stat_halfeye() + theme_minimal()
# plotconditions + xlab("c (bias)") #plot it
# 
# #while not very relevant, we can also plot the overall effect of Condition
# maineffectplot <- m1 %>%
#   emmeans(~ Condition,
#           at = list(Veracity = 0.5),
#           epred = FALSE, re_formula = NA) %>%
#   contrast(method = "pairwise") %>%
#   gather_emmeans_draws()
# 
# plotAME <- ggplot(plotidea,
#                   aes(x = .value, y = "Grand AME")) +
#   stat_halfeye() + theme_minimal()
# plotAME
# #Or both side by side
# # requires library(patchwork)
# (plotAME / plotconditions) +
#   plot_annotation(title = "Global grand mean") + xlab("c (bias)")

# difference in dprime between Conditions
m1.pw.diff <- contrast (m1.pw.accuracy, interaction = "pairwise") 
describe_posterior(m1.pw.diff) 
plot(m1.pw.diff)

# model fit
m1.R2_mz <- bayes_R2_mz(m1)
m1.loo <- loo(m1)

#fixed effects for LieType (repeated factor), Veracity, and their Interaction
#random effects for Participants and for Stimuli
#random slopes based on veracity and LieType
m2 <- brm(Answer ~ 1 + LieType*Veracity + (1 + LieType*Veracity | Participant) + (1 | Stimuli),  
          data = sdt.data, 
          family = bernoulli(link = probit),
          warmup = 2000,
          iter = 4000,
          save_pars = save_pars(all = TRUE), 
          prior = sdt_priors,
          backend = "cmdstanr",
          init = 0,
          thin = 1,
          chains = 4,
          seed = 123,
          control = list(adapt_delta = 0.9))

# model output
conditional_effects(m2)
describe_posterior(m2)

# pairwise comparisons
#bias of each LieType level and their difference
m2.pw.bias <- emmeans(m2, specs = revpairwise ~LieType)
describe_posterior(m2.pw.bias)

#accuracy of each LieType
m2.pw.accuracy <- emmeans (m2, ~ Veracity + LieType) #make the interaction

m2.pw.each <- contrast (m2.pw.accuracy, "pairwise", by = "LieType")
describe_posterior(m2.pw.each)# dprime for each LieType level

m2.pw.diff <- contrast (m2.pw.accuracy, interaction = "pairwise") 
describe_posterior(m2.pw.diff)# difference in dprime between LieType levels 

# model fit
m2.R2_mz <- bayes_R2_mz(m2)
m2.loo <- loo(m2)

# Full model:
#fixed effects for Condition (bw), LieType (ws), Veracity, and their Interaction
#random effects for Participants and for Stimuli
#random slopes based on veracity and LieType
m3 <- brm(Answer ~ 1 + LieType*Condition*Veracity + (1 + LieType*Veracity | Participant) + (1 | Stimuli),  
          data = sdt.data, 
          family = bernoulli(link = probit),
          warmup = 2000,
          iter = 4000,
          save_pars = save_pars(all = TRUE), 
          prior = sdt_priors,
          backend = "cmdstanr",
          init = 0,
          thin = 1,
          chains = 4,
          seed = 123,
          control = list(adapt_delta = 0.9))

# model output
conditional_effects(m3)
describe_posterior(m3)

# pairwise comparisons
#note, there are a lot of possible contrasts; run only hypothesised ones

##BIAS
#bias difference in within-subjects (ws) factor [akin to main effect]
m3.pw.bias.ws <- emmeans(m3, specs = revpairwise ~LieType)
describe_posterior(m3.pw.bias.ws)
#bias difference in between-subjects (bw) factor [akin to main effect]
m3.pw.bias.bw <- emmeans(m3, specs = revpairwise ~Condition)
describe_posterior(m3.pw.bias.bw)
#bias difference of ws factor at each bw level
m3.pw.bias.wsxbw <- emmeans(m3, specs = revpairwise ~LieType | Condition)
describe_posterior(m3.pw.bias.wsxbw)
#bias difference of bw factor at each ws level
m3.pw.bias.bwxws <- emmeans(m3, specs = revpairwise ~Condition | LieType)
describe_posterior(m3.pw.bias.bwxws)

##ACCURACY
#While I am sure there is a better way to do this, for clarity, I break down
#the interaction piece by piece

#create Veracity interaction (i.e., accuracy) for LieType, Condition, and their X
m3.pw.accuracy.2way.ws <- emmeans (m3, ~ Veracity + LieType)
m3.pw.accuracy.2way.bw <- emmeans (m3, ~ Veracity + Condition)
m3.pw.accuracy.3way <- emmeans (m3, ~ Veracity + LieType + Condition)

#accuracy of each LieType (ws) [ignoring interactions]
m3.pw.acc.ws <- contrast (m3.pw.accuracy.2way.ws, "pairwise", by = "LieType")
describe_posterior(m3.pw.acc.ws)
#accuracy difference of LieType levels [ignoring interactions]
m3.pw.acc.ws.diff <- contrast (m3.pw.accuracy.2way.ws, interaction = "pairwise") 
describe_posterior(m3.pw.acc.ws.diff) 

#accuracy of each Condition (bw) [ignoring interactions]
m3.pw.acc.bw <- contrast (m3.pw.accuracy.2way.bw, "pairwise", by = "Condition")
describe_posterior(m3.pw.acc.bw)
#accuracy difference between Condition groups [ignoring interactions]
m3.pw.acc.bw.diff <- contrast (m3.pw.accuracy.2way.bw, interaction = "pairwise")
describe_posterior(m3.pw.acc.bw.diff)

#accuracy for each LieType, at each Condition group
m3.pw.acc.ws.x <- contrast (m3.pw.accuracy.3way, "pairwise", by = c("LieType", "Condition"))
describe_posterior(m3.pw.acc.ws.x)
#accuracy difference of LieType levels, at each Condition group
m3.pw.acc.ws.x.diff <- contrast (m3.pw.accuracy.3way, interaction = "pairwise", by = "Condition")
describe_posterior(m3.pw.acc.ws.x.diff)

#accuracy for each Condition group, at each LieType level
m3.pw.acc.bw.x <- contrast (m3.pw.accuracy.3way, "pairwise", by = c("Condition", "LieType"))
describe_posterior(m3.pw.acc.bw.x) #note, this is just a reordered m3.pw.acc.ws.x
#accuracy difference of Condition groups, at each LieType level
m3.pw.acc.bw.x.diff <- contrast (m3.pw.accuracy.3way, interaction = "pairwise", by = "LieType")
describe_posterior(m3.pw.acc.bw.x.diff)

# The below comparison is nonsensical in my model, but may have some utility in
# other designs: it's a contrast of a contrast of a contrast; namely, it compares
# if the accuracy difference in classifying two levels of the WS factor is 
# different between two groups of the BW factor.
#difference in accuracy between LieType levels based on difference in Condition levels
m3.diff.of.diff.of.diff <- contrast (m3.pw.accuracy.3way, interaction = "pairwise")
describe_posterior(m3.diff.of.diff.of.diff)

# model fit
m3.R2_mz <- bayes_R2_mz(m3)
m3.loo <- loo(m3)

#### Compare models to find best one =========================================
#while "best" is relative (as theory-driven should take precedent), if you do 
#not know which factors to keep or discard, you can easily compare models using 
#'loo'; this will tell you which model best explains the data while considering 
#model complexity; thus it will identify the most parsimonious model. As a 
#rule-of-thumb a model is better only if the decrease in error is ~3x it's SE_diff

loo_models <- loo_compare(m0.loo, m1.loo, m2.loo, m3.loo)
loo_models #the output is sorted from "best" model to worst

perf_models <- compare_performance(m0, m1, m2, m3) #all-in-one fit checker (slow)
perf_models

comp_models <- compare_performance(m0, m1, m2, m3, rank = TRUE) #identifies best
comp_models

#this makes a radar plot of the models on all metrics; cool, but not very useful
plot(compare_performance(m0, m1, m2, m3, rank = TRUE))

#CONCLUSION ==================================================================
#there are many aspects of bayesian modelling that I do not cover here, such as:
#predicting new data or groups that do not appear in the current sample or 
#Bayesian Model Averaging, where we can combine the predictions of all models to
#estimate how much uncertainty there truly is if we didn't just pick "the best"
#model (see brms::posterior_average() and pp_average()), like what JASP does.

#-----If you found this script useful, please cite it from OSF---------------#

### References
# DeCarlo, L. T. (1998). Signal Detection Theory and Generalized Linear Models.
#   Psychological Methods 3 (2): 186–205.
# Macmillan, N. A., & Creelman C.D. (2005). Detection Theory: A User’s Guide. 
#   2nd ed. Mahwah, N.J: Lawrence Erlbaum Associates.
# Milne, A. J. (2020, June 19). Balance, evenness, and entropy in rhythm 
#   perception and aesthetic evaluation. https://doi.org/10.17605/OSF.IO/STNYG 
# Vuorre, M. (2017, October 9).Bayesian Estimation of Signal Detection Models. 
#   https://sometimesir.com/posts/2017-10-09-bayesian-estimation-of-signal-detection-theory-models 

### Acknowledgements
# Mattan S. Ben-Shachar - help with sorting out code and inspiration
# Andrew J. Milne - R2_mz custom function code creator
# Matti Vuorre - inspiration 

# Custom function ============================================================
# Bayesian version of McKelvey and Zavoina's pseudo-R2 for binary and ordinal
# models (McKelvey and Zavoina, 1975). This pseudo-R2 closely approximates the 
# R2 that would have been obtained if a linear model had have been run on 
# observations of the continuous latent variable underlying the discrete 
# responses.
# [RUN entire below code to create the bayes_R2_mz() function]
bayes_R2_mz <- function(fit, ...) {
  y_pred <- fitted(fit, scale = "linear", summary = FALSE, ...)
  var_fit <- apply(y_pred, 1, var)
  if (fit$formula$family$family == "cumulative" ||
      fit$formula$family$family == "bernoulli") {
    if (fit$formula$family$link == "probit" || 
        fit$formula$family$link == "probit_approx") {
      var_res <- 1
    }
    else if (fit$formula$family$link == "logit") {
      var_res <- pi^2 / 3 
    }
  } 
  else {
    sum_fit <- summary(fit)
    sig_res <- sum_fit$spec_pars["sigma","Estimate"]
    var_res <- sig_res^2
  } 
  R2_mz <- var_fit / (var_fit + var_res)
  return(
    data.frame(
      Estimate = mean(R2_mz), 
      Est.Error = sd(R2_mz), 
      Q2.5 = quantile(R2_mz, 0.025),
      Q97.5 = quantile(R2_mz, 0.975),
      row.names = "R2_mz"
    )
  )
  # return(R2_mz)
}
# ================================= END ======================================