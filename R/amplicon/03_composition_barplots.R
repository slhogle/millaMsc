library(here)
library(tidyverse)
library(Polychrome)
library(withr)
library(patchwork)

source(here::here("R", "utils_generic.R"))

# Read data ---------------------------------------------------------------
counts <- read_rds(here::here("data", "amplicon", "20230308_YSK_btk", "corrected_species_counts.rds")) |>
  mutate(strainID = str_replace(strainID, "-", "_"))
metadf <- read_tsv(here::here("data", "amplicon", "20230308_YSK_btk", "amplicon_metadata.tsv")) 

counts_f <- left_join(metadf, counts) |>
  group_by(sample) |>
  mutate(f = count_correct / sum(count_correct)) |>
  ungroup() |>
  filter(!is.na(strainID)) |>
  mutate(treatment = factor(
    treatment,
    levels = c("bact", "bact_strep",
               "bact_pred", "bact_pred_strep"),
    labels = c(
      "Bacteria alone",
      "Bacteria + streptomycin",
      "Bacteria + ciliate",
      "Bacteria + streptomycin + ciliate"
    )
  ))
strain.order <- counts_f %>%
  group_by(strainID) %>%
  summarize(tot=mean(f)) %>%
  arrange(tot) %>%
  pull(strainID)

# colors
with_seed(234523,
          firstpal <- unname(createPalette(23, c("#F3874AFF", "#FCD125FF"), M=5000)))

names(firstpal) <- strain.order

swatch(firstpal)

# Plot --------------------------------------------------------------------

myartheme <- function(...){
  theme(
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
    legend.background = element_blank(),
    legend.key = element_blank(),
    ...)
}

area_plot <- counts_f |>
  mutate(strainID=factor(strainID, levels=strain.order)) |>
  ggplot() +
  geom_area(aes(x=week, y=f, fill=strainID),
            color="black", size=0.1) +
  facet_grid(replicate ~ treatment) +
  scale_fill_manual(values = firstpal, ) + 
  scale_y_continuous(limits = c(0,1), expand = c(0, 0), labels = scales::percent) +
  scale_x_continuous(limits = c(0,80), breaks = c(0, 20, 40, 60, 80)) +
  labs(x="Sampling time (weeks)", y="Sp. abundance", fill="", title="") + 
  theme_bw() + 
  myartheme()

ggsave(here::here("figs", "amplicon", "20230308_YSK_btk", "species_comp_area.png"), area_plot, width=13, height=7, units="in",
       device="png", dpi=320)

ggsave(here::here("figs", "amplicon", "20230308_YSK_btk", "species_comp_area.svg"), area_plot, width=11, height=7, units="in",
       device="svg")

