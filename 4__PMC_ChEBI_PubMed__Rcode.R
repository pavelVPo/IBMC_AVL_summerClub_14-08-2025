library(tidyverse)

# Prepare
in_path <- ".../summer_club/types_of_studies_14-08_2025/data/"
out_path <- in_path
chebi_pmid <- read_tsv(str_glue("{in_path}chebi_pmids__11-08-25.tsv")) |> select(pmid) |> distinct()
obs_pmid <- read_tsv(str_glue("{in_path}obs_ids__12-08-2025.tsv"))
meta_pmid <- read_tsv(str_glue("{in_path}meta_ids__12-08-2025.tsv"))

# Intersection
chebi_obs_pmid <- chebi_pmid |> inner_join(obs_pmid)
chebi_meta_pmid <- chebi_pmid |> inner_join(meta_pmid)

# Save the results
write_tsv(chebi_obs_pmid, str_glue("{out_path}obs_chebi_pmids__12-08-25.tsv"))
write_tsv(chebi_meta_pmid, str_glue("{out_path}meta_chebi_pmids__12-08-25.tsv"))