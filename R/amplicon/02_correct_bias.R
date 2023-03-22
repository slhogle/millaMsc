library(here)
library(tidyverse)
#pak::pkg_install("mikemc/metacal")
library(metacal)
library(Polychrome)
library(withr)

# Read data ---------------------------------------------------------------
counts    <- readr::read_tsv(here::here("data", "amplicon", "species_counts.tsv"))
expdesign <- readr::read_tsv(here::here("data", "amplicon", "amplicon_metadata.tsv"))
pos_ctrl  <- readr::read_tsv(here::here("data", "amplicon", "pos_control.tsv")) |>
  mutate(ra=pos_control_cells_ml/sum(pos_control_cells_ml))

# Prepare -----------------------------------------------------------------
# We will use the metacal package for estimating bias and performing 
# calibration in the special case where the bias of all the taxa of interest 
# can be directly measured from the control sample. Since we know the exact 
# concentration of cells added in the beginning of the experiment we can 
# estimate a corrected/unbiased relative abundance.

# Make a matrix of observed counts
observed_mat <- left_join(counts, expdesign) |> 
  filter(str_detect(sample, "Community_control")) |>
  left_join(pos_ctrl) |>
  select(sample, strainID, count) |>
  pivot_wider(names_from="strainID", values_from="count") |>
  column_to_rownames(var="sample") |>
  as.matrix()

# Make a matrix of true proportions
actual_mat <- left_join(counts, expdesign) |> 
  filter(str_detect(sample, "Community_control")) |>
  left_join(pos_ctrl) |>
  select(sample, strainID, ra) |>
  pivot_wider(names_from="strainID", values_from="ra") |>
  column_to_rownames(var="sample") |>
  as.matrix()

# Bias estimation ---------------------------------------------------------
with_seed(12378, mc_fit <- estimate_bias(observed_mat, actual_mat, 1, boot=TRUE) |> print())
            
bias <- coef(mc_fit)
print(bias)
mc_fit.summary <- summary(mc_fit)
coef_tb <- mc_fit.summary[['coefficients']]

# Plot bias estimation
bias_est <- coef_tb |>
  mutate(taxon = fct_reorder(taxon, estimate)) |>
  ggplot(aes(taxon, estimate, 
             ymin = estimate / gm_se^2, ymax = estimate * gm_se^2)) +
  geom_hline(yintercept = 1, color = "grey") +
  geom_pointrange() +
  scale_y_log10() +
  coord_flip() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(
  here::here("figs", "amplicon", "metalcal_bias_est.svg"),
  plot = bias_est,
  device = "svg",
  units = "cm",
  width = 17.8,
  height = 12.8
)

# Plot model fit
mypal <- with_seed(14123, unname(createPalette(23, c("#F3874AFF", "#FCD125FF"), M=5000)))

a <- as_tibble(fitted(mc_fit)) |>
  pivot_longer(everything())|> 
  mutate(type="Fitted")

b <- as_tibble(actual_mat) |>
  pivot_longer(everything())|> 
  mutate(type="Actual")

c <- as_tibble(observed_mat/rowSums(observed_mat)) |>
  pivot_longer(everything(),  values_to="observed")

metacal_fit_plot <- bind_rows(a,b) |>
  left_join(c) |>
  ggplot(aes(x=observed, y=value, color = name)) +
  geom_abline(color = "darkgrey") +
  geom_point(show.legend = T) +
  scale_color_manual(values=mypal) +
  facet_grid(~type) + 
  scale_x_sqrt() + 
  scale_y_sqrt() + 
  coord_fixed() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(
  here::here("figs", "amplicon", "metalcal_fit.svg"),
  plot = metacal_fit_plot,
  device = "svg",
  units = "cm",
  width = 17.8,
  height = 17.8
)

# Calibrate ---------------------------------------------------------------

# Make a matrix of observed counts
ps_mat <- left_join(counts, expdesign) |> 
  left_join(pos_ctrl) |>
  select(sample, strainID, count) |> 
  pivot_wider(names_from="strainID", values_from="count") |>
  column_to_rownames(var="sample") |>
  as.matrix()

# run the calibrate function
with_seed(435761, ps_cal <- calibrate(ps_mat, bias, margin=1))

counts_final <- rownames_to_column(as.data.frame(ps_cal)) |>
  rename(sample=rowname) |>
  pivot_longer(`HAMBI-0006`:last_col(), names_to="strainID", values_to="f_correct") |>
  left_join(counts) |>
  group_by(sample) |>
  mutate(total=sum(count)) |>
  ungroup() |>
  mutate(count_correct=total*f_correct,
         f_obs=count/total)

with_seed(435761,
          samples <- counts |>
            select(sample) |>
            distinct() |>
            sample_n(10) |>
            pull())

metalcal_abund_correct <- counts_final |>
  filter(sample %in% samples) |>
  dplyr::select(sample, strainID, f_correct, f_obs) |>
  pivot_longer(cols=c(f_correct, f_obs), names_to="type", values_to="prop") |>
  ggplot(aes(x=type, prop, fill = strainID)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values=mypal) +
  facet_wrap(~sample, ncol = 5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(
  here::here("figs", "amplicon", "metalcal_abund_correct.svg"),
  plot = metalcal_abund_correct,
  device = "svg",
  units = "cm",
  width = 20.8,
  height = 10.8
)

# Output ------------------------------------------------------------------
# What I am doing here is multiplying each species relative abundance by the 
# total number of reads per sample so you get corrected counts. Counts are 
# necessary for DEseq normalization etc...

# Also make sure to remove the negative contols
counts_final_fmt <- counts_final |>
  mutate(count_correct = round(count_correct)) |>
  filter(!str_detect(sample, "Community_control")) |>
  dplyr::select(sample, strainID, genus, species, count, count_correct, f_correct)
  
write_tsv(counts_final_fmt, here::here("data", "amplicon", 
                                       "corrected_species_counts.tsv"))
