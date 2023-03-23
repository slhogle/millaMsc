library(here)
library(tidyverse)
source(here::here("R", "utils_generic.R"))

untar(tarfile=here::here("data_raw", "amplicon", "bbmap_rpkm.tar.gz"), 
      exdir=here::here("data_raw", "amplicon"))

mapdir <- here::here("data_raw", "amplicon", "bbmap_rpkm")
filenames <- list.files(mapdir, full.names = TRUE, pattern = "*.rpkm")

files <- set_names(filenames, 
                   str_extract(filenames,
                     regex("(?<=[/])([^/]+)(?=\\.[^.]+)")
                   ))

counts <- map_df(
  files,
  read_tsv,
  comment = "#",
  col_names = c(
    "strainID",
    "Length",
    "Bases",
    "Coverage",
    "count",
    "RPKM",
    "Frags",
    "FPKM"
  ),
  .id = "sample"
) |>
  left_join(tax, by = join_by(strainID)) |>
  select(sample, strainID, genus, species, count)

counts_wide <- counts |>
  # build group index
  group_by_at(vars(-count)) |>
  mutate(row_id = 1:n()) |> 
  ungroup() |>
  pivot_wider(names_from = sample, values_from = count) |> 
  # drop the index
  select(-row_id) |> 
  drop_na()

write_tsv(counts, here("data", "amplicon", "species_counts.tsv"))
write_tsv(counts_wide, here("data", "amplicon", "species_counts_wide.tsv"))

# Clean up ----------------------------------------------------------------
unlink("data_raw/amplicon/bbmap_rpkm", recursive = T)
