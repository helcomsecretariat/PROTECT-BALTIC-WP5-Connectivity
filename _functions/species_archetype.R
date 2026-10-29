# Species -> archetype matrix and PLD (NA if not passive)
species_archetype <- function() {
  sp <- read.csv("WP5_species_list.csv")
  aa <- read.csv("WP5_archetype_assignments.csv")
  aa <- aa[!duplicated(aa$archetype),
           c("archetype", "archetype_matrix", "drift_duration_average", "season")]
  m <- merge(sp[, c("scientific_name", "archetype")], aa, by = "archetype", all.x = TRUE)
  names(m)[names(m) == "drift_duration_average"] <- "pld"
  m
}
