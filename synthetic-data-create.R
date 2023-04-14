# Create synthetic dataset

# To run this script, ensure the raw data is in data-raw/SDT-Probit-JONB-v1.csv
# that's not checked to git because it can't be shared

# Libraries: synthpop for creating synthetic data
library(tidyverse)
library(brms)
library(synthpop)

# Load data
d <- read_csv("data-raw/SDT-Probit-JONB-v1.csv")

# Recode things
d <- d %>%
  mutate(
    Participant = fct_anon(factor(Participant)),
    Training = factor(
      Condition,
      levels = c("ERT", "BT", "NT"),
      labels = c("Emotion", "Bogus", "None")
    ),
    Stimulus = factor(str_split_i(Stimuli, pattern = "_", 1)),
    LieType = factor(
      LieType,
      levels = c("AFF", "EXP"),
      labels = c("Affective", "Experiential")
    ),
    Veracity = factor(Veracity, levels = 0:1, labels = c("False", "True")),
    Answer = Answer
  ) %>%
  select(Participant, Training, Stimulus, LieType, Veracity, Answer)

contrasts(d$Veracity) <- c(-0.5, 0.5)
d$Training <- fct_relevel(d$Training, "None")

fit <- brm(
  Answer ~ 1 + Veracity * Training * LieType + (1 + Veracity * LieType | Participant) + (1 | Stimulus),
  family = bernoulli(link = probit_approx),
  init = "0",
  data = d,
  backend = "cmdstanr",
  threads = 2,
  cores = 8
)
fit <- m2
summary(fit)
d$Answer <- posterior_predict(fit, ndraws = 1)[1,]
write_rds(d, "data/dataset-synthetic.rds", compress = "gz")
