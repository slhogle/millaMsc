library(here)
library(tidyverse)
library(Polychrome)
library(withr)

counts <- read_tsv(here::here("data", "amplicon", "corrected_species_counts.tsv")) |>
  # replace dash in strainID w/ underscore because it plays nicer w/ R
  mutate(strainID = str_replace(strainID, "-", "_"))
metadf <- read_tsv(here::here("data", "amplicon", "amplicon_metadata.tsv")) 

counts_f <- left_join(metadf, counts) |>
  filter(!is.na(strainID)) |>
  mutate(treatment = factor(
    treatment,
    levels = c("bact", "bact_strep"),
    labels = c(
      "Bacteria alone",
      "Bacteria + streptomycin")
  ))

strain.order <- counts_f %>%
  group_by(strainID) %>%
  summarize(tot=mean(f_correct)) %>%
  arrange(tot) %>%
  # pull converts a tibble dataframe to a vector
  pull(strainID)

with_seed(234523,
          firstpal <- unname(createPalette(23, c("#F3874AFF", "#FCD125FF"), M=5000)))

names(firstpal) <- strain.order

swatch(firstpal)

myartheme <- function(...){
  theme_bw() + theme(
    panel.spacing.x = unit(0.25,"line"),
    strip.placement = 'outside',
    strip.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    #axis.text.x = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    legend.title = element_blank(),
    legend.position="bottom",
    legend.background = element_blank(),
    legend.key = element_blank(),
    ...)
}

area_plot <- counts_f |>
  mutate(strainID=factor(strainID, levels=strain.order)) |>
  ggplot() +
  geom_area(aes(x=week, y=f_correct, fill=strainID),
            color="black", size=0.1) +
  facet_grid(treatment ~ replicate) +
  scale_fill_manual(values = firstpal, ) + 
  scale_y_continuous(limits = c(0,1), expand = c(0, 0), labels = scales::percent) +
  scale_x_continuous(limits = c(0,80), breaks = c(0, 20, 40, 60, 80)) +
  labs(x="Sampling time (weeks)", y="Species abundance", fill="", title="") + 
  myartheme()

ggsave(
  here::here("figs", "amplicon", "species_comp_area.svg"),
  area_plot,
  width = 13,
  height = 7,
  units = "in",
  device = "svg"
)

richness <- function(count_vec){
  sum(count_vec > 0)
}

shannon <- function(f){
  # log of zero is undefined
  -sum(f[f>0] * log(f[f>0]))
}

diversity <- counts_f |>
  group_by(sample) |>
  summarize(richness = richness(count_correct),
            shannon = shannon(f_correct)) |>
  ungroup() |>
  left_join(metadf)

diversity |>
  # to get alpha and shannon into 'long' format
  pivot_longer(cols = c("richness", "shannon"), names_to="index_name", values_to="index_val") |>
  mutate(treatment = factor(
    treatment,
    levels = c("bact", "bact_strep"),
    labels = c(
      "Bacteria alone",
      "Bacteria + streptomycin")
  )) |>
  ggplot() +
  geom_line(aes(x=week, y=index_val, group=interaction(index_name, treatment, replicate))) +
  geom_point(aes(x=week, y=index_val)) +
  facet_grid(index_name ~ treatment, scales="free") +
  labs(x="Sampling time (weeks)", y="Diversity", fill="", title="") + 
  myartheme()