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
      "20230217_bioscreen_prelim_run",
      "bioscreen_out_recoded.txt"
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
  "20230217_bioscreen_prelim_run",
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
  mutate(strep_conc_lab = paste0(strep_conc, "~mu*g/ml")) |>
  mutate(strep_conc_lab = factor(
    strep_conc_lab,
    levels = c(
      "0~mu*g/ml",
      "1~mu*g/ml",
      "2.5~mu*g/ml",
      "5~mu*g/ml",
      "12.5~mu*g/ml",
      "25~mu*g/ml",
      "50~mu*g/ml",
      "100~mu*g/ml",
      "500~mu*g/ml",
      "5000~mu*g/ml"
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
  facet_grid(strep_conc_lab ~ strainID,
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

ggsave(here::here("figs", "competition", "20230217_bioscreen_prelim_run", "01_growthcurve_strep.png"), 
       pgc, 
       width=6, 
       height=12, 
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
    color = strep_conc_lab,
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

ggsave(here::here("figs", "competition", "20230217_bioscreen_prelim_run", "02_growthcurve_strep.png"), 
       pgc2, 
       width=7, 
       height=5, 
       units="in", 
       device="png", 
       #dpi=320
)
# Growth rate plot --------------------------------------------------------


pgr <- am_fmt |>
  group_by(strainID, evolution, strep_conc, parameter) |>
  summarize(mn = mean(value), sdv = sd(value)) |>
  ggplot() +
  geom_pointrange(aes(x=strep_conc, y=mn, ymin = mn-sdv, ymax = mn+sdv, color = evolution)) + 
  geom_line(aes(x=strep_conc, y=mn, group= evolution, color=evolution)) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 5, 50, 500, 5000)) + 
  scale_color_fish(option = "Synchiropus_splendidus", 
                   discrete = TRUE, direction = -1) +
  facet_grid(parameter ~ strainID, scales="free_y",
             labeller = labeller(parameter = label_wrap_gen(width = 20))) +
  labs(x=expression("Streptomycin concentration"~mu*"g/ml")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.placement = 'outside',
    strip.background = element_blank()
  )

ggsave(here::here("figs", "competition", "20230217_bioscreen_prelim_run", "03_growthrate_strep.png"), 
       pgr, 
       width=6, 
       height=10, 
       units="in", 
       device="png", 
       #dpi=320
)

