# Mission: unsupervised clustering analyses

# check and set working directory
getwd()

# call renv::restore() to install packages from lockfile
renv::restore()

# set up renv
renv::init()
# install core setup + utility packages
renv::install("tidyverse")
renv::install("here")
renv::install("usethis")
renv::install("pak")
renv::install("devtools")
renv::install("ddauber/r4np")
renv::install("patchwork")
renv::install("tictoc")
renv::install("brooke-watson/BRRR")
renv::install("tsibble")

# machine learning
renv::install("YuHuiDeakin/rabc")
renv::install("theft")
renv::install("theftdlc")

# load packages
library(tidyverse)
library(here)
library(usethis)
library(pak)
library(devtools)
library(r4np)
library(patchwork)
library(tictoc)
library(BRRR)
library(tsibble)

# machine learning
library(rabc)
library(theft)
library(theftdlc)


# save packages to lockfile
renv::snapshot()

# ensure variables don't save between sessions
usethis::use_blank_slate()

# set up folder system using r4np package - if not already
# ensure r4np is loaded before running the below code (hashed out)
# r4np::create_project_folder()

# set up git
# usethis::use_git()
# usethis::use_github() # use this code to create GitHub repo, but I already did this before manually

# load helper functions
source("02_r_scripts/utils.R")

# ensure variables don't save between sessions
usethis::use_blank_slate()

# section 1: load raw ACC data and format for analysis -----
## vulture example - Vaadia et al 2025
# this data is already formatted into segments.
# we want to reformat it into a tsibble,
# common format so we can easily do all subsequent analyses

# load data
## Vultures training dataset Vaadia et al 2025 ------------
vultures <- vroom::vroom(here("00_raw_data/vulture_training_dataset.csv"))

# explore data
glimpse(vultures)
head(vultures)
n_distinct(vultures$device_id)

# data structure:
# bout id - each unique segment
# device id - individual id
# harness - device attachment style
# observed_beh - labelled behaviour
# acc_x - rows 1:100, for a given segment. Sampling at 20Hz, segment 5s long
# acc_y - same as above
# acc_z - same as above
# then made their own features based on these acc measurements -
# can be another feature set for us

# time is not explicit in this dataset, is assumed.
# We either need a timestamp or an index to assign as our time index for the tsibble

fs <- 20 # Hz; 100 samples per bout = 5 s
vultures_long <- vultures |>
  select(
    id = device_id,
    bout = bout_id,
    attachment = harness,
    label = observed_beh,
    matches("^acc_[xyz]_\\d+$")
  ) |>
  pivot_longer(
    cols = matches("^acc_[xyz]_\\d+$"),
    names_to = c(".value", "sample_in_bout"),
    names_pattern = "(acc_[xyz])_(\\d+)"
  ) |>
  mutate(
    id = as.factor(id),
    sample_in_bout = as.integer(sample_in_bout),
    t_in_bout = (sample_in_bout - 1) / fs
  ) |>
  relocate(
    id,
    bout,
    sample_in_bout,
    t_in_bout,
    acc_x,
    acc_y,
    acc_z,
    attachment,
    label
  )

vultures_tsbl <- as_tsibble(
  vultures_long,
  key = c(id, bout),
  index = sample_in_bout
)

# print our tsibble
vultures_tsbl

# check the data is ok
nrow(vultures_tsbl) == nrow(vultures) * 100 # 578,300 expected
any(has_gaps(vultures_tsbl)$.gaps) # should be FALSE
count(vultures, device_id, bout_id) |> filter(n > 1) # should return 0 rows
# if this final check returns rows, bout_id is not unique within a device
# as_tsibble() will then erorr on duplicate key_index pairs, and bout will need to be made unique first
# e.g. with consecutive_id() per device

# if we want seconds within a bout for plotting or spectral work,
# add this in after building the tsibble.
# keep sample_in_bout as the index
# vultures_tsbl <- mutate(vultures_tsbl, t_in_bout = (sample_in_bout - 1) / fs)

# data exploration
n_distinct(vultures_tsbl$id)
# how many bouts per ind?
vultures_tsbl %>%
  as_tibble() %>%
  group_by(id) %>%
  summarise(n_bouts = n_distinct(bout)) %>%
  arrange(desc(n_bouts))

# how many behaviours & bouts per ind?
vultures_tsbl %>%
  as_tibble() %>%
  group_by(id, label) %>%
  summarise(n_bouts = n_distinct(bout)) %>%
  arrange(desc(n_bouts))
# View()

# how many behaviour bouts per ind - showing all behaviours per ind?
vultures_tsbl %>%
  as_tibble() %>%
  group_by(id, label) %>%
  summarise(n_bouts = n_distinct(bout)) %>%
  arrange(desc(n_bouts)) %>%
  pivot_wider(names_from = label, values_from = n_bouts)
# %>%
# View()

# stacked bar chart + total bouts per ind underneath
vulture_bouts <- vultures_tsbl |>
  as_tibble() |>
  distinct(id, bout, label) # one row per bout

id_order_vultures <- vulture_bouts |> count(id, sort = TRUE) |> pull(id) # most bouts first

p_n_vult <- vulture_bouts |>
  count(id, name = "n_bouts") |>
  ggplot(aes(factor(id, id_order_vultures), n_bouts)) +
  geom_col(fill = "grey40") +
  labs(x = NULL, y = "Bouts") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

p_prop_vult <- vulture_bouts |>
  count(id, label) |>
  ggplot(aes(factor(id, id_order_vultures), n, fill = label)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "ID", y = "Share of bouts", fill = "Behaviour") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

p_n_vult / p_prop_vult + plot_layout(heights = c(1, 2), guides = "collect")


# vizz data with plot_signal
# 3 functions
# 1) time mode: simple data streams over time, stitched bouts together if so
# 2) behaviour mode: plot all behaviours for a given id
# 3) individuals mode: plot all individuals for a given behaviour
# can specify however many data streams you want!
# Really easy

plot_signal(vultures_tsbl, acc_x, mode = "time", id = "202367", fs = 20)
plot_signal(vultures_tsbl, c(acc_x, acc_y, acc_z), id = "202367", fs = 20)
plot_signal(
  vultures_tsbl,
  c(acc_x, acc_y, acc_z),
  mode = "behaviours",
  id = "202367",
  fs = 20
)
plot_signal(
  vultures_tsbl,
  c(acc_x, acc_y, acc_z),
  mode = "behaviours",
  id = "202378",
  fs = 20
)

plot_signal(
  vultures_tsbl,
  acc_x:acc_z,
  mode = "individuals",
  behaviour = "Standing",
  fs = 20
)
plot_signal(
  vultures_tsbl,
  acc_x:acc_z,
  mode = "individuals",
  behaviour = "Flapping",
  fs = 20
)
plot_signal(
  vultures_tsbl,
  acc_x:acc_z,
  mode = "individuals",
  behaviour = "Soaring",
  fs = 20
)

# feature calculation -----
tic()
vulture_features <- calc_features(
  vultures_tsbl,
  c(acc_x, acc_y, acc_z),
  feature_set = c("catch22", "rabc_time2", "rabc_freq"),
  winlen_dba_s = 1, # same as original winlen_dba = 21
  fs = 20,
  catch24 = TRUE
)
toc()
BRRR::skrrrahh("soulja")


calc_features(
  vultures_tsbl,
  c(acc_x, acc_y, acc_z),
  feature_set = c("catch22", "rabc_time2", "rabc_freq"),
  winlen_dba_s = 1,
  fs = 20
)
calc_features(
  vultures_tsbl,
  acc_x,
  feature_set = "catch22",
  features = list(mean = mean, sd = sd)
)
calc_features(
  vultures_tsbl,
  c(acc_x, acc_y, acc_z),
  feature_set = NULL,
  features = list(mean = mean, sd = sd)
)
calc_features(
  vultures_tsbl,
  c(acc_x, acc_y, acc_z),
  feature_set = NULL,
  features = list(mean = mean, range = \(v) max(v) - min(v))
)

calc_features(
  vultures_tsbl,
  c(acc_x, acc_y, acc_z),
  feature_set = NULL,
  winlen_dba_s = 1,
  fs = 20,
  features = list(
    mean = mean,
    odba = cross_feature(
      c("acc_x", "acc_y", "acc_z"),
      function(x, y, z, winlen) {
        sba <- function(v) roll_static(v, winlen)
        mean(
          abs(x - sba(x)) + abs(y - sba(y)) + abs(z - sba(z)),
          na.rm = TRUE
        )
      },
      window = TRUE
    )
  )
)

v2 <- calc_features(
  vultures_tsbl,
  c(acc_x, acc_y, acc_z),
  feature_set = c("catch22", "rabc_time2", "rabc_freq"),
  features = list(
    mean = mean,
    odba = cross_feature(
      c("acc_x", "acc_y", "acc_z"),
      function(x, y, z, winlen) {
        sba <- function(v) roll_static(v, winlen)
        mean(
          abs(x - sba(x)) + abs(y - sba(y)) + abs(z - sba(z)),
          na.rm = TRUE
        )
      },
      window = TRUE
    )
  ),
  winlen_dba_s = 1, # same as original winlen_dba = 21
  fs = 20,
  catch24 = TRUE
)

glimpse(vulture_features)
head(vulture_features)
# decide whether to set winlen_dba as based on segment sample length (0.05s for vultures)
# or based on seconds.
# for now have gone with seconds, but can change back if needed

# winlen sensitvity - figure out wtf this is doing
# res <- winlen_sensitivity(vultures_tsbl, fs = 20)
res <- winlen_sensitivity(
  vultures_tsbl,
  fs = 20,
  winlen_s = c(0.25, 0.5, 0.75, 1, 1.5, 2, 3),
  ref_s = 1
)
res$overview # one row per window
res$plots$separation # also $odba_by_label, $stability_cost

# plot dimension_reduction ------
#' calc_features(vultures_tsbl, c(acc_x, acc_y, acc_z), feature_set = "catch22") |>
#'   project_features(norm_method = "RobustSigmoid", unit_int = TRUE,
#'                    low_dim_method = "PCA") |>
#'   plot()

v2 %>%
  project_features(
    norm_method = "zScore",
    unit_int = TRUE,
    low_dim_method = "UMAP",
    seed = 12
  ) %>%
  plot()

v2 %>%
  project_features(
    norm_method = "zScore",
    feature_set = "rabc_time2",
    unit_int = TRUE,
    low_dim_method = "UMAP",
    seed = 12
  ) %>%
  plot()

v2 %>%
  project_features(
    norm_method = "zScore",
    feature_set = "rabc_time2",
    unit_int = TRUE,
    low_dim_method = "PCA",
    seed = 12
  ) %>%
  plot()

v2 %>%
  project_features(
    norm_method = "zScore",
    feature_set = "rabc_time2",
    unit_int = TRUE,
    low_dim_method = "tSNE",
    seed = 12
  ) %>%
  plot()

# make object to inspect
v2_proj <- v2 %>%
  project_features(
    norm_method = "zScore",
    feature_set = "rabc_time2",
    unit_int = TRUE,
    low_dim_method = "PCA",
    seed = 12
  )

# plot dimension reduction for a specific behaviour only - e.g. flapping
v2 %>%
  #   mutate(label = if_else(label == "Flapping", "Flapping", "Other")) %>%
  project_features(
    norm_method = "zScore",
    feature_set = "rabc_time2",
    unit_int = TRUE,
    low_dim_method = "UMAP",
    seed = 12
  ) %>%
  plot()

# to do:
# create function that directly calculates outputs of rabc df_time and df_freq
# so I can plot a UMAP in the shiny app
# and then eventually build my own replacement function
# / the project function from theft, but adapted to my data/specified feature sets!
# main mission now is to vizz how diff features affect UMAP output,
# and start unsupervised feature selection!

# import another dataset - cats ----
cats_raw <- vroom::vroom(here(
  "00_raw_data/Dunford_et_al._Cats_calibrated_data.csv"
))
glimpse(cats_raw)

cats <- cats_raw %>%
  select(
    id = ID,
    time = Time,
    acc_x = AccX,
    acc_y = AccY,
    acc_z = AccZ,
    label = Behaviour
  ) %>%
  mutate(id = as.factor(id), attachment = "collar") %>%
  relocate(attachment, .before = label) %>%
  arrange(id, time)

# now need to separate into bouts
# cats |> add_segments(min_gap_s = 5, label_col = "label")
# cats |> make_bouts(min_gap_s = 5, nrow = 100)   # label ignored entirely

cats_tsbl <- cats |>
  add_segments(min_gap_s = 5, label_col = "label") |>
  make_bouts(min_gap_s = 5, duration_s = 2, fs = 40)

# run workflow on these to check feature calculation and UMAP works

# then can go onto clustering or some other thing to test & select features
