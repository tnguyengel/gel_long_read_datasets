# Retrieve all unsolved Rare Disease participants from the 100,000 Genomes Project.
# Exclude any withdrawn participants.
# Do not filter based on whether participant is in a trio or not.


library(tidyverse)
library(Rlabkey)



# Max rows to retrieve in a single query to labkey. Labkey can only handle 1M rows.
# If you need more than that, then you need to page the results.
MAX_ROWS <- 1e6  # 1000000000
LABKEY_URL <- "https://labkey-embassy.gel.zone/labkey/"
# We are using v17 data release of the NGRL, which includes 100KG Project data as well as some GMS data
LABKEY_FOLDER <- "/main_programme/main-programme_v17_2023-03-30"

# Max approximate FST distance from individual to UK Biobank reference population for ancestry assignment
FST_THRESH <- 0.002

# Publicly available ancestry inference for 100KG data release v17 using UK Biobank as reference populations
ANCESTRY_CSV <- "/gel_data_resources/gel_diverse_data/100k/310523/100k-main_programme_v17_ukbb_worldwide_populations_310523.csv.gz"


###########  Query Labkey

labkey.setDefaults(baseUrl = LABKEY_URL)


participant_sql <- "SELECT 
  p.participant_id, 
  p.participant_type,
  p.year_of_birth,
  p.programme,
  p.handling_gmc_trust,
  p.registration_date,
  p.rare_diseases_family_id,
  p.participant_ethnic_category,
  p.participant_phenotypic_sex,
  p.participant_karyotypic_sex,
  p.participant_pipeline_status 
  FROM participant as p 
  WHERE p.participant_pipeline_status != 'Withdrawn'"

participant_df <- labkey.executeSql(
    maxRows=MAX_ROWS,
    folderPath=LABKEY_FOLDER,
    schemaName="lists",
    sql=participant_sql,
    colNameOpt="rname") %>% 
    as_tibble() %>% 
    mutate(participant_id = as.character(participant_id),
           rare_diseases_family_id = as.character(rare_diseases_family_id))


rd_analysis_sql <- "SELECT 
  r.participant_id,
  r.plate_key,
  r.rare_diseases_family_id,
  r.reported_karyotypic_sex,
  r.inferred_sex_karyotype,
  r.biological_relationship_to_proband 
  from rare_disease_analysis r"

rd_analysis_df <- labkey.executeSql(
    maxRows=MAX_ROWS,
    folderPath=LABKEY_FOLDER,
    schemaName="lists",
    sql=rd_analysis_sql,
    colNameOpt="rname") %>% 
    as_tibble() %>% 
    mutate(participant_id = as.character(participant_id),
           rare_diseases_family_id = as.character(rare_diseases_family_id))

rd_interpreted_sql <- "SELECT 
  r.participant_id,
  r.normalised_specific_disease_proband, 
  r.affection_status,
  r.case_solved_family 
  from rare_disease_interpreted r"

rd_interpreted_df <- labkey.executeSql(
    maxRows=MAX_ROWS,
    folderPath=LABKEY_FOLDER,
    schemaName="lists",
    sql=rd_interpreted_sql,
    colNameOpt="rname") %>% 
    as_tibble() %>% 
    mutate(participant_id = as.character(participant_id))

rd_disease_sql <- "SELECT 
      d.participant_id, 
      d.normalised_specific_disease,
      d.normalised_disease_sub_group,
      d.normalised_disease_group
      FROM rare_diseases_participant_disease as d"

rd_disease_df <- labkey.executeSql(
    maxRows=MAX_ROWS,
    folderPath=LABKEY_FOLDER,
    schemaName="lists",
    sql=rd_disease_sql,
    colNameOpt="rname") %>% 
    as_tibble() %>% 
    mutate(participant_id = as.character(participant_id))


rd_family_sql <- "SELECT 
  r.rare_diseases_family_id,
  r.family_group_type,
  r.family_medical_review_qc_state_code 
  FROM rare_diseases_family r"


rd_family_df <- labkey.executeSql(
    maxRows=MAX_ROWS,
    folderPath=LABKEY_FOLDER,
    schemaName="lists",
    sql=rd_family_sql,
    colNameOpt="rname") %>% 
    as_tibble() %>% 
    mutate(rare_diseases_family_id = as.character(rare_diseases_family_id))


###########  Assign ancestries to the 100KG participants based on genetic similarity to UK Biobank reference populations

# See https://re-docs.genomicsengland.co.uk/gen_sim/
#
# We only assign an ancestry if the participant is within the approximate FST threshold to the reference population.
# FST is approximated here using normalized Euclidean distance to the UK Biobank reference populations, as per Prive et al. 2022.
# A participant can only be assigned one ancestry.  This means that participants with mixed ancestry 
# will be Unassigned.
# We assign the UK reference population that is closest in approximate FST to the participant.

ancestry_df <- read_csv(ANCESTRY_CSV) %>%
  mutate(`Participant ID` = as.character(`Participant ID`))

ancestry_col_names <- grep('Approximate FST to', names(ancestry_df), value = TRUE)

ancestry_df <- ancestry_df %>%
  rowwise %>%
  mutate(min_fst = min(c_across(ancestry_col_names))) %>%
  # toString() will add a comma when there are multiple entries
  mutate(min_ancestry = toString(ancestry_col_names[which(c_across(ancestry_col_names) == min_fst)])) %>%
  ungroup() %>%
  mutate(participant_genetic_category = case_when(
     min_fst < FST_THRESH & str_detect(min_ancestry, ",") == F ~ str_replace(min_ancestry, "Approximate FST to ", ""),
     T                                                          ~ "Unassigned"
  ))


curr_year <- as.integer(format(Sys.Date(), "%Y"))

unsolved_rd_df <- participant_df %>%
  inner_join(rd_analysis_df, by = c("participant_id", "rare_diseases_family_id")) %>%
  left_join(rd_interpreted_df, by = c("participant_id")) %>%
  left_join(rd_family_df, by = c("rare_diseases_family_id")) %>%
  left_join(ancestry_df, by=c("participant_id" = "Participant ID", "plate_key" = "Platekey")) %>%
  filter(case_solved_family == "no") %>%
  mutate(age  = curr_year - year_of_birth) 


write.table(unsolved_rd_df %>% select(
              participant_id,
              plate_key,
              handling_gmc_trust,
              rare_diseases_family_id,
              biological_relationship_to_proband, 
              age,
              inferred_sex_karyotype,
              normalised_specific_disease_proband,
              affection_status,
              case_solved_family,
              participant_ethnic_category,
              participant_genetic_category),
  file = "./unsolved_rd_cases_meta_v17.tsv", row.names=FALSE, sep="\t", quote=F)
