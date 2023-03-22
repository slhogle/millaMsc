tax <- tibble::tribble(
  ~strainID,             ~genus,        ~species,
  "HAMBI-0006",      "Pseudomonas",        "putida",
  "HAMBI-0097",    "Acinetobacter",       "lwoffii",
  "HAMBI-0105",    "Agrobacterium",   "tumefaciens",
  "HAMBI-0262",    "Brevundimonas",       "bullata",
  "HAMBI-0403",        "Comamonas",  "testosteroni",
  "HAMBI-1279",           "Hafnia",         "alvei",
  "HAMBI-1287",      "Citrobacter",        "koseri",
  "HAMBI-1292",       "Morganella",      "morganii",
  "HAMBI-1299",         "Kluyvera",    "intermedia",
  "HAMBI-1842",      "Sphingobium",    "yanoikuyae",
  "HAMBI-1896", "Sphingobacterium",  "spiritivorum",
  "HAMBI-1972",        "Aeromonas",        "caviae",
  "HAMBI-1977",      "Pseudomonas",  "chlororaphis",
  "HAMBI-1988",     "Chitinophaga",        "sancti",
  "HAMBI-2159", "Paraburkholderia",   "caryophylli",
  "HAMBI-2160",       "Bordetella",         "avium",
  "HAMBI-2164",      "Cupriavidus",       "necator",
  "HAMBI-2443",       "Paracoccus", "denitrificans",
  "HAMBI-2494", "Paraburkholderia",   "kururiensis",
  "HAMBI-2659", "Stenotrophomonas",   "maltophilia",
  "HAMBI-2792",        "Moraxella",         "canis",
  "HAMBI-3031",         "Niabella",  "yanshanensis",
  "HAMBI-3237",       "Microvirga",   "lotononidis"
)

# ggplot theme
mytheme = function() {
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank()
  )
}

# opposite of %in% fuction
`%nin%` = Negate(`%in%`)

# logit transform
logit = function(x){
  log(x/(1-x))
}

minnz = function(V) {
  # Calculates the smallest value of the vector except for 0 (non-zero minumum)
  # Argument: vector
  C <- NULL        # prepare
  k <- length(V)   # count to
  for (i in 1:k) { # check all
    if ((V[i] == 0) == FALSE) (C[i] <- V[i]) else (C[i] <- 9999919) # if V[i] is not 0, add it to C
  }
  m <- min(C)               # minimum of V, not counting 0
  if (max(V) == 1) (m <- 1) # fix for binary vectors (0,1)
  if (m == 9999919) (warning("Error: Minimum calculation failed."))  # warning because of hard-coded replacement
  return(m)
}

quibble95 = function(x, q = c(0.025, 0.5, 0.975)) {
  tibble(x = quantile(x, q), quantile = c("q2.5", "q50", "q97.5"))
}


