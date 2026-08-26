library(here)
library(tidyverse)
library(multilevel)

# load data 
## to request the data, please contact the first author 
esm1 <- readRDS(here("data", "pace_esm_data_s1.rds"))
esm2 <- readRDS(here("data", "pace_esm_data_s2.rds"))
dat <- bind_rows(esm1, esm2)


baseline1 <- readRDS(here("data", "pace_baseline_data_s1.rds")) 
baseline2 <- readRDS(here("data", "pace_baseline_data_s2.rds"))
baseline <- bind_rows(baseline1, baseline2)


followup1 <- readRDS(here("data", "pace_followup_data_s1.rds")) 
followup2 <- readRDS(here("data", "pace_followup_data_s2.rds"))
followup <- bind_rows(followup1, followup2)


# Broad overview ----------------------------------------------------------
dat_c <- dat %>% # only complete ESM surveys
  filter(!is.na(irri))

ids <- distinct(dat_c, session) %>% 
  pull(session)

length(ids)

# to get the number of incomplete surveys, only keep participants with at least one complete survey
dat <- dat %>% 
  filter(session %in% ids)

nrow(dat) # all ESM observations including incomplete responses


# Table with completion rates ---------------------------------------------


# Count ESM observations per participant
esm_counts <- dat_c %>%
  group_by(session, sample, condition) %>%
  summarise(n_obs = n(), .groups = "drop")

# Create follow-up completion indicators
fudavis_sessions <- followup %>% filter(sample == "davis") %>% pull(session) 
fuubc_sessions <- followup %>% filter(sample == "ubc") %>% pull(session)
fude_sessions <- followup %>% filter(sample == "de") %>% pull(session)

# Add follow-up completion
esm_counts <- esm_counts %>%
  mutate(
    completed_fu = case_when(
      sample == "davis" ~ session %in% fudavis_sessions,
      sample == "ubc" ~ session %in% fuubc_sessions,
      sample == "de" ~ session %in% fude_sessions
    )
  )

# Function to compute completion metrics per site and condition
get_metrics <- function(data, site_label, sample_filter) {
  d <- data %>% filter(sample %in% sample_filter)
  
  tibble(
    Site = site_label,
    Metric = c("Participants started ESM",
               "Completed ≥25 ESM surveys",
               "Completed ≥60 ESM surveys",
               "Completed entire study"),
    A = c(
      sum(d$condition == "A"),
      sum(d$condition == "A" & d$n_obs >= 25),
      sum(d$condition == "A" & d$n_obs >= 60),
      sum(d$condition == "A" & d$completed_fu)
    ),
    B = c(
      sum(d$condition == "B"),
      sum(d$condition == "B" & d$n_obs >= 25),
      sum(d$condition == "B" & d$n_obs >= 60),
      sum(d$condition == "B" & d$completed_fu)
    )
  ) %>%
    mutate(
      Total = A + B,
      pct_starters = Total / Total[1] * 100,
      pct_starters = ifelse(Metric == "Participants started ESM", NA, round(pct_starters, 2))
    )
}

# Compute for each site
table_davis <- get_metrics(esm_counts, "UC Davis", "davis")
table_ubc <- get_metrics(esm_counts, "UBC", "ubc")
table_de <- get_metrics(esm_counts, "HU Berlin & U Siegen", "de")

# Combine
table_final <- bind_rows(table_davis, table_ubc, table_de)
writexl::write_xlsx(table_final, "completion_rates_table.xlsx")

# Figure completion rates -------------------------------------------------


tmp <- dat_c %>%
  group_by(session, condition) %>%
  summarise(completed_surveys = n(), .groups = "drop") %>%
  mutate(
    condition = factor(condition),
    completed_surveys = as.integer(completed_surveys)
  ) %>%
  count(condition, completed_surveys, name = "n_participants")


figure_3 <- ggplot(
  tmp,
  aes(x = factor(completed_surveys), y = n_participants)
) +
  geom_col(fill = "grey70", color = "grey70") +
  facet_wrap(
    ~ condition,
    nrow = 1,
    labeller = as_labeller(c(
      A = "Condition A – \nHigh Frequency",
      B = "Condition B – \nHigh Duration"
    ))
  ) +
  labs(
    x = "Number of completed surveys per participant",
    y = "Number of participants"
  ) +
  scale_x_discrete(breaks = c("25", "50", "75", "100"))+
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80))+
  theme_classic(base_size = 18) +
  theme(
    panel.spacing.x = unit(3, "lines"),
    axis.title = element_text(face = "bold", size = 14),
    axis.text  = element_text(size = 14),
    # strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(
      face = "bold",
      size = 20
    ))


ggsave("Figure_3.png", figure_3,
       width = 8, height = 4,   # larger = better quality
       dpi = 300)                # 300+ for publication
 
figure_3



# Sample description ------------------------------------------------------

# only keep final ids
ids_esm <- distinct(dat_c, session) %>% pull(session)
baseline <- baseline %>% 
  filter(session %in% ids_esm)

# First look at unique major values as specified in an open text-format
baseline %>%
  filter(sample %in% c("davis", "ubc")) %>%
  count(major) %>%
  arrange(desc(n)) %>%
  print(n = 50)

# Check remaining rows for any psychology variants we might have missed
baseline %>%
  filter(sample %in% c("davis", "ubc")) %>%
  filter(str_detect(major, regex("psy", ignore_case = TRUE))) %>%
  count(major) %>%
  arrange(desc(n)) %>%
  print(n = Inf)


#Count all entries as psychology major that have "psy" in it in North American universities
#in Germany, participants were simply asked whether they study psychology (1) or not (2)
baseline <- baseline %>%
  mutate(is_psychology = case_when(
    sample %in% c("davis", "ubc") ~ str_detect(major, regex("psy", ignore_case = TRUE)),
    sample == "de" ~ psychology == 1, 
    TRUE ~ NA
  ))


# Next see all unique values that were self-specified for ethnicity
baseline %>%
  filter(sample %in% c("davis", "ubc")) %>%
  count(ethnicity) %>%
  arrange(desc(n)) %>%
  print(n = Inf)


# Inspect the frequencies of each category
## Based on that overview, it seems useful to create new categories for (a) Hispanic/Latino/a, (b) Middle Eastern/Arab, and (c) Multiracial 
## Categories 3 (Native American), 4 (Alaska Native), and 6 (Native Hawaiian or Other Pacific Islander) were so rare that they were merged into "Another Ethnic Identity" Category 

baseline %>%
  filter(sample %in% c("davis", "ubc")) %>%
  mutate(ethnicity_recode = case_when(
    # Asian
    ethnicity == "5" ~ "Asian",
    ethnicity %in% c("South Asian", "south asian", "Chinese", "chinese and white",
                     "Korean", "korean", "Filipino-American", "Laos", 
                     "Southeast Asian", "southeast asian", "east asian",
                     "East Asian") ~ "Asian",
    
    # White/Caucasian
    ethnicity == "1" ~ "White/Caucasian",
    ethnicity %in% c("Brazilian; white") ~ "White/Caucasian",
    
    # Black/African American
    ethnicity == "2" ~ "Black/African American",
    ethnicity %in% c("East African") ~ "Black/African American",
    
    # Hispanic/Latino/a
    ethnicity %in% c("Hispanic", "hispanic", "Latina", "latina", "Latino", 
                     "latino", "Mexican", "mexican", "Chicana", 
                     "Hispanic or Latinx", "Hispanic/Latino", "Latino/Hispanic",
                     "Hispanic/Latina", "Hispanic/latina", "hispanic/latina",
                     "Latine", "Latino/a", "Mestizo", "Mexican American",
                     "Latin American", "Latin American (Mexican)", 
                     "European Latin American", "hispanic mexican",
                     "mixed race/hispanic", "Hispanic/Latinx") ~ "Hispanic/Latino/a",
    
    # Middle Eastern/Arab
    ethnicity %in% c("Middle Eastern", "middle eastern", "Persian", "Afghan",
                     "Egyptian", "arab", "arabic", "Arab", "west asian",
                     "middle eastern/arab", "Middle eastern/arab") ~ "Middle Eastern/Arab",
    
    # Multiracial
    ethnicity %in% c("Mixed", "mixed", "Biracial (Mexican + Hmong)", 
                     "Biracial: Asian and White", "Caucasian and Asian",
                     "Caucasian and Latina", "Mixed black and white",
                     "Mixed; White and African American", "Multiracial",
                     "multiracial", "Southeast Asian and Caucasian",
                     "White and East Asian", "White/East Asian", "White/Korean",
                     "asian and white", "mixed white and chinese",
                     "mixed: white; east asian", "white and asian",
                     "Mexican - White/Native American", "arab and white",
                     "chinese and white") ~ "Multiracial",
    
    # Prefer not to say
    ethnicity == "7" ~ "Prefer not to say",
    
    # Another Ethnic Identity
    ethnicity %in% c("3", "4", "6", "indgenous", "Indian") ~ "Another Ethnic Identity",
    
    TRUE ~ "CHECK - not coded"
  ))  %>% 
  count(ethnicity_recode) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  arrange(desc(n))

#Middle Eastern/Arab is also pretty rare, so merge that with the "Other" category as well 

# save the new categories in the dataset 
baseline <- baseline %>%
  mutate(ethnicity_recode = case_when(
    sample == "de" ~ NA_character_, # ethnicity was not assessed at the German universities
    ethnicity == "5" ~ "Asian",
    ethnicity %in% c("South Asian", "south asian", "Chinese", "chinese and white",
                     "Korean", "korean", "Filipino-American", "Laos", 
                     "Southeast Asian", "southeast asian", "east asian",
                     "East Asian") ~ "Asian",
    ethnicity == "1" ~ "White/Caucasian",
    ethnicity %in% c("Brazilian; white") ~ "White/Caucasian",
    ethnicity == "2" ~ "Black/African American",
    ethnicity %in% c("East African") ~ "Black/African American",
    ethnicity %in% c("Hispanic", "hispanic", "Latina", "latina", "Latino", 
                     "latino", "Mexican", "mexican", "Chicana", 
                     "Hispanic or Latinx", "Hispanic/Latino", "Latino/Hispanic",
                     "Hispanic/Latina", "Hispanic/latina", "hispanic/latina",
                     "Latine", "Latino/a", "Mestizo", "Mexican American",
                     "Latin American", "Latin American (Mexican)", 
                     "European Latin American", "hispanic mexican",
                     "mixed race/hispanic", "Hispanic/Latinx") ~ "Hispanic/Latino/a",
    ethnicity %in% c("Mixed", "mixed", "Biracial (Mexican + Hmong)", 
                     "Biracial: Asian and White", "Caucasian and Asian",
                     "Caucasian and Latina", "Mixed black and white",
                     "Mixed; White and African American", "Multiracial",
                     "multiracial", "Southeast Asian and Caucasian",
                     "White and East Asian", "White/East Asian", "White/Korean",
                     "asian and white", "mixed white and chinese",
                     "mixed: white; east asian", "white and asian",
                     "Mexican - White/Native American", "arab and white",
                     "chinese and white") ~ "Multiracial",
    ethnicity %in% c("Middle Eastern", "middle eastern", "Persian", "Afghan",
                     "Egyptian", "arab", "arabic", "Arab", "west asian",
                     "middle eastern/arab", "Middle eastern/arab") ~ "Another Ethnic Identity",
    ethnicity == "7" ~ "Prefer not to say",
    ethnicity %in% c("3", "4", "6", "indgenous", "Indian") ~ "Another Ethnic Identity",
    TRUE ~ "CHECK - not coded"
  )) %>%
  mutate(ethnicity_recode = factor(ethnicity_recode,
                                   levels = c("Asian", "White/Caucasian", "Hispanic/Latino/a",
                                              "Black/African American", "Multiracial",
                                              "Another Ethnic Identity", "Prefer not to say")))

# Verify
baseline %>%
  count(sample, ethnicity_recode) %>%
  print(n = Inf)


### Create the final table

get_site_stats <- function(data, site) {
  d <- data %>% filter(sample == site)
  n <- nrow(d)
  
  # Sample composition (DE only)
  if (site == "de") {
    hu <- round(mean(d$uni == 1, na.rm = TRUE) * 100, 1)
    siegen <- round(mean(d$uni == 2, na.rm = TRUE) * 100, 1)
    other <- round(mean(d$uni == 3, na.rm = TRUE) * 100, 1)
    comp <- paste0("HU Berlin: ", hu, "%, U Siegen: ", siegen, "%, Other: ", other, "%")
  } else {
    comp <- "-"
  }
  
  # Age
  age_mean_sd <- paste0(round(mean(d$age, na.rm = TRUE), 1), 
                        " (", round(sd(d$age, na.rm = TRUE), 1), ")")
  
  # Gender helper
  gender_fn <- function(val) {
    count <- sum(d$gender == val, na.rm = TRUE)
    pct <- round(count / n * 100, 1)
    paste0(count, " (", pct, "%)")
  }
  
  # Psychology
  psych_count <- sum(d$is_psychology, na.rm = TRUE)
  psych_pct <- round(psych_count / n * 100, 1)
  psych <- paste0(psych_count, " (", psych_pct, "%)")
  
  # Ethnicity
  if (site == "de") {
    eth <- rep("-", 7)
  } else {
    eth_fn <- function(cat) {
      count <- sum(d$ethnicity_recode == cat, na.rm = TRUE)
      pct <- round(count / n * 100, 1)
      paste0(count, " (", pct, "%)")
    }
    eth <- unname(sapply(c("Asian", "White/Caucasian", "Hispanic/Latino/a",
                           "Black/African American", "Multiracial",
                           "Another Ethnic Identity", "Prefer not to say"), eth_fn))
  }
  
  # Return as named vector
  c(
    "Sample composition" = comp,
    "N" = as.character(n),
    "Mean age (SD)" = age_mean_sd,
    "Gender" = "",
    "Women" = gender_fn(1),
    "Men" = gender_fn(2),
    "Non-binary" = gender_fn(3),
    "Gender: Prefer not to say" = gender_fn(4),
    "Self-specified" = gender_fn(5),
    "Psychology major" = psych,
    "Ethnicity" = "",
    "Asian" = eth[1],
    "White/Caucasian" = eth[2],
    "Hispanic/Latino/a" = eth[3],
    "Black/African American" = eth[4],
    "Multiracial" = eth[5],
    "Another Ethnic Identity" = eth[6],
    "Ethnicity: Prefer not to say" = eth[7]
  )
}

# Compute for each site
davis_stats <- get_site_stats(baseline, "davis")
ubc_stats <- get_site_stats(baseline, "ubc")
de_stats <- get_site_stats(baseline, "de")

# Combine into table
table_final <- tibble(
  Variable = names(davis_stats),
  `UC Davis` = unname(davis_stats),
  `UBC` = unname(ubc_stats),
  `HU Berlin & U Siegen` = unname(de_stats)
)

print(table_final)

writexl::write_xlsx(table_final, "sample_composition.xlsx")



# Descriptive statistics of key ESM variables -----------------------------

# Define variables and their labels
var_info <- tibble(
  construct = c(
    rep("Situation Characteristics", 4),
    rep("Personality States", 15),
    rep("Emotion States", 2),
    rep("Biological States", 2)
  ),
  var = c(
    "duty", "sociality", "negativity", "typicality",
    "soci", "asse", "acti", "trus", "kind", "poli", "crea", "curi", "arti", "lazy", "orga", "reli", "irri", "rela", "care",
    "happ", "stres",
    "tire", "hung"
  ),
  label = c(
    "duty", "sociality", "negativity", "typicality",
    "sociable", "assertive", "active", "trustful", "kind", "polite", "creative", "curious", "artistic", "lazy", "organized", "reliable", "irritable", "relaxed", "carefree",
    "happy", "stressed",
    "tired", "hungry"
  )
)

compute_stats <- function(var_name, data) {
  d <- data %>%
    dplyr::select(session, value = dplyr::all_of(var_name)) %>%
    dplyr::filter(!is.na(value))
  
  # Overall mean and range
  overall_mean <- round(mean(d$value), 2)
  overall_min <- round(min(d$value), 2)
  overall_max <- round(max(d$value), 2)
  overall_range <- paste0(overall_min, "-", overall_max)
  
  # Between-person SD (SD of person means)
  person_means <- d %>%
    dplyr::group_by(session) %>%
    dplyr::summarise(person_mean = mean(value), .groups = "drop")
  between_sd <- round(sd(person_means$person_mean), 2)
  
  # Within-person SD
  person_within_sd <- d %>%
    dplyr::group_by(session) %>%
    dplyr::summarise(within_sd = sd(value), .groups = "drop") %>%
    dplyr::filter(!is.na(within_sd))
  
  within_sd_mean <- round(mean(person_within_sd$within_sd), 2)
  within_sd_min <- round(min(person_within_sd$within_sd), 2)
  within_sd_max <- round(max(person_within_sd$within_sd), 2)
  within_sd_range <- paste0(within_sd_min, "-", within_sd_max)
  
  tibble(
    mean = overall_mean,
    range = overall_range,
    between_sd = between_sd,
    within_sd_mean = within_sd_mean,
    within_sd_range = within_sd_range
  )
}

# Apply to all variables
desc_table <- var_info %>%
  dplyr::mutate(stats = purrr::map(var, ~compute_stats(.x, dat_c))) %>%
  tidyr::unnest(stats)



writexl::write_xlsx(desc_table, "esm_descriptives.xlsx")
