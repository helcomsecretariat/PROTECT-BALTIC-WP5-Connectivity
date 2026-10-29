library(dplyr)
library(tidyr)

spec_df <- read.csv("WP5_species_list.csv")
spec_df <- select(spec_df, scientific_name)
spec_df$do_run <- 1
spec_df$run_started <- 0

# Run species with overall_min_score >= 2
df <- read.csv("inputs/WP3_species_review.csv")
df <- select(df, scientific_name, overall_min_score)
spec_df <- left_join(spec_df, df, by = "scientific_name")

spec_df$do_run <- 0
spec_df$do_run[spec_df$overall_min_score >= 2] <- 1

# subset
#s <- seq(1, 630, by = 10)
#s <- c(30, 32, 67, 463, 500)
#spec_df$do_run <- 0
#spec_df$do_run[s] <- 1

write.csv(spec_df, "run_log.csv", row.names = FALSE)

# Create and clear the multi-session claim folders
dir.create("_run_log")

dir.create("_run_log/matrices")
dir.create("_run_log/matrices/started")
dir.create("_run_log/matrices/finished")

dir.create("_run_log/archetypes")
dir.create("_run_log/archetypes/started")
dir.create("_run_log/archetypes/finished")

dir.create("_run_log/strength")
dir.create("_run_log/strength/started")
dir.create("_run_log/strength/finished")

unlink("_run_log/matrices/started/*")
unlink("_run_log/matrices/finished/*")

unlink("_run_log/archetypes/started/*", recursive = TRUE)
unlink("_run_log/archetypes/finished/*", recursive = TRUE)

unlink("_run_log/strength/started/*", recursive = TRUE)
unlink("_run_log/strength/finished/*", recursive = TRUE)
