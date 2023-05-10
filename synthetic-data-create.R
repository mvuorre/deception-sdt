# Create synthetic dataset

# To run this script, ensure the raw data is in data-raw/SDT-Probit-JONB-v1.csv
# that's not checked to git because it can't be shared

library(tidyverse)
library(brms)

# Load data
d <- read_csv("data-raw/SDT-Probit-JONB-v1.csv")

# Original data values
# - `Training` is a between-subjects variable indicating which training group the participant was in,
# - `LieType` is a within-subjects variable; what type of lie was presented on the trial
# - `Stimulus` indicates which stimulus was presented
# - `Veracity` indicates the stimulus veracity (True or False)
# - `Answer` is the DV, as each answer participants give; truth=1, lie=0
count(d, Veracity, Answer)

# Recode things
d <- d %>%
  mutate(
    Participant = fct_anon(factor(Participant), "p"),
    Training = factor(
      Condition,
      levels = c("NT", "BT", "ERT"),
      labels = c("None", "Bogus", "Emotion")
    ),
    LieType = factor(
      LieType,
      levels = c("AFF", "EXP"),
      labels = c("Affective", "Experiential")
    ),
    # Remove file extension from stimulus labels
    Stimulus = str_split_i(Stimuli, "\\.", 1),
    # Is stimulus a lie?
    isLie = factor(Veracity, levels = c(1, 0), labels = c("No", "Yes")),
    # Did participant respond with a "lie"?
    sayLie = if_else(Answer == 0, 1, 0)
  ) %>%
  select(Participant, Training, LieType, Stimulus, isLie, sayLie)
d
count(d, isLie, sayLie)
contrasts(d$isLie) <- c(-0.5, 0.5)

fit <- brm(
  sayLie ~ 1 + isLie * Training * LieType +
    (1 + isLie * LieType | Participant) +
    (1 | Stimulus),
  family = bernoulli(link = probit),
  init = "0",
  data = d,
  backend = "cmdstanr",
  threads = 2,
  cores = 8,
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.95)
)
summary(fit)

# Hacky way to simulate good dataset from the posterior
library(tidybayes)
library(posterior)
max_draw <- spread_draws(fit, lp__) %>%
  filter(lp__ == max(lp__)) %>%
  pull(.draw)
d$sayLie <- posterior_predict(fit, ndraws = 1, draw_ids = max_draw)[1,]
write_rds(d, "data/dataset-synthetic.rds", compress = "gz")
