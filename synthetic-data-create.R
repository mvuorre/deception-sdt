# Create synthetic dataset

# To run this script, ensure the raw data is in data-raw/SDT-Probit-JONB-v1.csv
# that's not checked to git because it can't be shared

# Libraries: synthpop for creating synthetic data
library(tidyverse)
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
    Answer = factor(Answer)
  ) %>%
  select(Participant, Training, Stimulus, LieType, Veracity, Answer)


# Create a synthetic dataset for each participant.
# Doing this by participant is much faster and avoids complexities.
# The only variable that needs to be synthesized is `Answer`,
# because all others are determined by experimenter.
d <- d %>%
  nest(.by = Participant) %>%
  mutate(
    data2 = map(
      data,
      ~syn(
        data = .x,
        method = c("constant", "", "", "", "cart")
      ) %>%
        .[["syn"]] %>%
        tibble()
    )
  )

# Write synthetic dataset to disk
dir.create("data", FALSE)
d %>%
  select(Participant, data2) %>%
  unnest(data2) %>%
  mutate(Answer = as.integer(Answer) - 1) %>%
  write_rds("data/dataset-synthetic.rds", compress = "gz")
