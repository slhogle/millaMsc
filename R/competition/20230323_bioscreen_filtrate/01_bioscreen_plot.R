library(here)
library(tidyverse)
library(stringr)
library(stringi)
library(lubridate)
library(fs)
library(fishualize)

# Read data ---------------------------------------------------------------

# OD data
od <-
  read_table(
    here::here(
      "data_raw",
      "competition",
      "20230323_bioscreen_filtrate",
      "Iina_EXP_Nutrient-spike_recoded.txt"
    ),
    skip = 2,
    locale = locale(
      encoding = "UTF-8",
      decimal_mark = ",",
      grouping_mark = "."
    )
  )
  

# growth summary data from amiga
mapdir <- here::here("data", "competition", "20230217_bioscreen_prelim_run")
files <- list.files(mapdir, full.names = TRUE, pattern="_summary.txt")
files <- set_names(files,
                   str_extract(files,
                               regex("(?<=[/])([^/]+)(?=\\.[^.]+)")))

am <- map_df(
  files,
  read_tsv,
  comment = "#")


# metadata
mdf <- read_tsv(here::here(
  "data_raw",
  "competition",
  "20230323_bioscreen_filtrate",
  "bioscreen_metadata.tsv"
))


# Format data -------------------------------------------------------------

od_fmt <- od |>
  mutate(Time = lubridate::period_to_seconds(lubridate::hms(Time))) |>
  dplyr::select(-Blank) |>
  pivot_longer(-Time, values_to = "OD600", names_to = "bioscreen_well") |>
  drop_na() |>
  mutate(bioscreen_well = as.numeric(bioscreen_well)) |>
  left_join(mdf) |>
  mutate(KB_spike_conc_perc_lab = paste0(KB_spike_conc_perc, "% KB spike")) |>
  mutate(KB_spike_conc_perc_lab = factor(
    KB_spike_conc_perc_lab,
    levels = c(
      "0% KB spike",
      "2.5% KB spike",
      "5% KB spike",
      "7.5% KB spike",
      "10% KB spike"
    )
  ))

amdf <- mdf |>
  mutate(Plate_ID = paste0("plate_", plate)) |>
  separate(well, into = c("a", "b"), sep = "[A-Z]", remove=FALSE) |>
  mutate(a = str_sub(well, start=1, end=1),
         b = as.numeric(b)) |>
  mutate(Well = paste0(a, b)) |>
  dplyr::select(Plate_ID, Well, strep_conc:last_col())
  
am_fmt <- left_join(am, amdf, by = join_by(Well, Plate_ID)) |>
  dplyr::select(Plate_ID, Well, strep_conc, strainID, evolution, k_lin, gr, td, lagC, t_gr, t_k) |>
  mutate(lagC = ifelse(gr < 0.05, NA_real_, lagC),
         t_gr = ifelse(gr < 0.05, NA_real_, t_gr),
         t_k = ifelse(gr < 0.05, NA_real_, t_k),
         td = ifelse(gr < 0.05, NA_real_, td)) |>
  pivot_longer(k_lin:last_col(), names_to = "parameter") |>
  mutate(parameter = factor(parameter, levels = c("gr", "td", "k_lin",
                                                  "lagC", "t_gr", "t_k"),
                            labels = c("Growth rate (1/hr)",
                                       "Doubling time (hr)",
                                       "Carrying capactiy (OD units)",
                                       "Lag time (hr)",
                                       "Time at max growth rate (hr)",
                                       "Time at carrying capacity (hr)")))

# Plot data ---------------------------------------------------------------

pgc <- ggplot(od_fmt) +
  geom_vline(xintercept = 24*60*60, linetype = "dashed", color = "grey80") +
  geom_vline(xintercept = 30*60*60, linetype = "dashed", color = "grey80") +
  geom_line(aes(
    x = Time,
    y = OD600,
    color = evolution,
    group = bioscreen_well
  )) +
  scale_color_fish(option = "Synchiropus_splendidus", 
                   discrete = TRUE, direction = -1) +
  scale_x_continuous(breaks = seq(0, 48*60*60, 6*60*60), 
                     labels = seq(0, 48, 6)) +
  facet_grid(KB_spike_conc_perc_lab ~ strainID) +
  labs(x="Elapsed time (hours)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.placement = 'outside',
    strip.background = element_blank()
  )

ggsave(here::here("figs", "competition", "20230323_bioscreen_filtrate", "01_growthcurve_filtrate.png"), 
       pgc, 
       width=6, 
       height=8, 
       units="in", 
       device="png", 
       #dpi=320
       )

pgc2 <- ggplot(od_fmt) +
  geom_vline(xintercept = 24*60*60, linetype = "dashed", color = "grey80") +
  geom_vline(xintercept = 30*60*60, linetype = "dashed", color = "grey80") +
  geom_line(aes(
    x = Time,
    y = OD600,
    color = KB_spike_conc_perc_lab,
    group = bioscreen_well
  )) +
  scale_color_viridis_d() + 
  scale_x_continuous(breaks = seq(0, 48*60*60, 6*60*60), 
                     labels = seq(0, 48, 6)) +
  facet_grid(evolution ~ strainID,
             labeller = label_parsed) +
  labs(x="Elapsed time (hours)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.placement = 'outside',
    strip.background = element_blank()
  )

ggsave(here::here("figs", "competition", "20230323_bioscreen_filtrate", "02_growthcurve_filtrate.png"), 
       pgc2, 
       width=7, 
       height=5, 
       units="in", 
       device="png", 
       #dpi=320
)
