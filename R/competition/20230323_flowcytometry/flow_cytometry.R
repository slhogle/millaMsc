library(here)
library(tidyverse)


# Read conc data from guava -----------------------------------------------

# these data were gated using the GuavaSoft software so we can plot concentrations
# without further preprocessing

cyt <-
  read_tsv(
    here::here(
      "data_raw",
      "competition",
      "20230323_flowcytometry",
      "2023-03-22-SLH_processed.tsv"
    ))

cp <- cyt |>
  pivot_longer(cols = -(Well:dilution), names_to = "gate", values_to = "concentration") |>
  mutate(gate_short = str_extract(gate, "([^:alpha:]+)(?=\\.[^.]+)"),
         conc_corr = concentration*dilution) |>
  separate(`Sample ID`, into = c("sp", "b")) |>
  mutate(lab = paste(sp, dilution, syber_stained)) |>
  #mutate(conc_corr1 = ifelse(Well %in% c("A10", "A11", "A12"), concentration*5000, conc_corr)) |>
  ggplot() +
  geom_point(aes(x=Well, y=conc_corr, color=gate_short)) +
  annotate("text", x = 2, y = 1e06, label = "100 dilution") +
  annotate("text", x = 5, y = 1e07, label = "500 dilution") +
  annotate("text", x = 8, y = 1e08, label = "1000 dilution") +
  annotate("text", x = 11, y = 1e05, label = "Unstained") +
  annotate("text", x = 13.5, y = 1e06, label = "Filtrate") +
  scale_y_log10() +
  theme_bw()


ggsave(here::here("figs", "competition", "20230323_flowcytometry", "01_flow_cyt_test.png"), 
       cp, 
       width=7, 
       height=5, 
       units="in", 
       device="png", 
       #dpi=320
)

# https://med.virginia.edu/flow-cytometry-facility/resources/r-script/

# https://rstudio-pubs-static.s3.amazonaws.com/817768_fec70ef38efb4da395c8b9b9547e345d.html

# https://jchellmuth.com/posts/FACS-with-R/