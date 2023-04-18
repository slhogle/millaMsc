library(here)
library(tidyverse)
library(stringr)
library(stringi)
library(lubridate)
library(fs)

# Deal with encoding issues -----------------------------------------------

# Note: in this case I needed to manually convert the file from UTF-16LE to UTF-8
# using iconv command line tool
# iconv -f utf-16le -t utf-8 bioscreen_out.txt > bioscreen_out_recoded.txt

guess_encoding(
  here::here(
    "data_raw",
    "competition",
    "20230217_bioscreen_prelim_run",
    "bioscreen_out.txt"
  )
)
guess_encoding(
  here::here(
    "data_raw",
    "competition",
    "20230217_bioscreen_prelim_run",
    "bioscreen_out_recoded.txt"
  )
)

# Read data ---------------------------------------------------------------

df <-
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
  ) |>
  mutate(Time = lubridate::period_to_seconds(lubridate::hms(Time)))

mdf <- read_tsv(here::here(
  "data_raw",
  "competition",
  "20230217_bioscreen_prelim_run",
  "bioscreen_metadata.tsv"
))

# Format data -------------------------------------------------------------

# Setup well list for 96 well plate. Because these are bioscreen runs it is a 
# bit of a pain because bioscreen uses weird 100 well format plates.

nsamps <- 96
combinations <-  expand.grid(stri_pad_left(1:12, pad="0", width=2), LETTERS[1:8])
wells <- paste0(combinations[1:nsamps, 2], combinations[1:nsamps, 1])

plate_no <- c(paste0("plate_", rep(1, 96)),
              paste0("plate_", rep(2, 96)),
              paste0("plate_", rep(3, 8))
)

well_id <- c(wells, wells, wells[1:8])

df_fmt <- df |>
  dplyr::select(-Blank) |>
  pivot_longer(-Time, values_to = "OD600") |>
  pivot_wider(names_from = Time, values_from = OD600) |>
  drop_na() |>
  mutate(plate_no = plate_no,
         well = well_id) |>
  relocate(well, plate_no) |>
  dplyr::select(-name) |>
  as_tibble()

# Write OD data -----------------------------------------------------------

# make a directory to write to
data_write_dir <- here::here(
  "data_raw",
  "competition",
  "20230217_bioscreen_prelim_run",
  "amiga",
  "data"
)

dir_create(
  data_write_dir,
  mode = "u=rwx,go=rx",
  recurse = TRUE
)

# write the OD data to amiga format
df_fmt |>
  group_by(plate_no) |>
  nest() |>
  mutate(file_out = paste0(plate_no, ".tsv"),
         file_out_path = here::here(data_write_dir, file_out)) |>
  transpose() |>
  walk(\(l) write_tsv(l$data, l$file_out_path))


# Write metadata ----------------------------------------------------------

# make a directory to write to
map_write_dir <- here::here(
  "data_raw",
  "competition",
  "20230217_bioscreen_prelim_run",
  "amiga",
  "mapping"
)

dir_create(
  map_write_dir,
  mode = "u=rwx,go=rx",
  recurse = TRUE
)

# write the metadata to amiga format
mdf |>
  dplyr::select(plate:last_col()) |>
  group_by(plate) |>
  nest() |>
  mutate(file_out = paste0("plate_", plate, ".map"),
         file_out_path = here::here(map_write_dir, file_out)) |>
  transpose() |>
  walk(\(l) write_tsv(l$data, l$file_out_path))
