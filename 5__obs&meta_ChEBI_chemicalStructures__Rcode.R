library(tidyverse)
library(europepmc)
library(xml2)

## Prepare
# Download the ChEBI data: structures.tsv.gz FROM https://ftp.ebi.ac.uk/pub/databases/chebi-2/flat_files/
# Paths
in_path <- ".../summer_club/types_of_studies_14-08_2025/data/"
out_path <- ".../summer_club/types_of_studies_14-08_2025/output/"
# Read
obs_pmid <- read_tsv(str_glue("{in_path}obs_chebi_pmids__12-08-25.tsv")) |> mutate(category_ = 'observational study')
meta_pmid <- read_tsv(str_glue("{in_path}meta_chebi_pmids__12-08-25.tsv")) |> mutate(category_ = 'meta-analysis')
structures_all <- read_tsv(str_glue("{in_path}chebi_structures.tsv")) |> select(compound_id, smiles) |>
					mutate(compound_id = as.character(compound_id))
# Wide form
pmids <- bind_rows(obs_pmid, meta_pmid) |> group_by(pmid) |>
						mutate(category = str_c(category_, sep = " & ", collapse = " & ")) |>
						ungroup() |>
						select(pmid, category) |>
						distinct() |>
						mutate(chebid = NA)


## Get the ChEBI IDs

# Using europepmc for R, 
# epmc_db_count(ext_id = '11958740', data_src = "med")
# Ошибка в if (doc$hitCount == 0) { :аргумент нулевой длины

# Using Europe PMC annotations API, SEE: https://europepmc.org/AnnotationsApi
for (i in seq(1:nrow(pmids))) {
	chebi_vec <- NA
	# Get annotations
	annotations <- read_xml(str_glue("https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByArticleIds?articleIds=MED%3A{as.character(pmids[i,1])}&format=XML")) |>
						as_list()
	# Create list to store ChEBI IDs
	chebi_vec <- rep(NA, annotations[["List"]][["item"]][["annotations"]] |> length())
	for (k in seq(1 : length(chebi_vec))) {
		chebi_id <- str_match(annotations[["List"]][["item"]][["annotations"]][[k]][["tags"]][[1]][["uri"]], "CHEBI_.*")
		if (!is.na(chebi_id)) {
			chebi_id <- str_replace(chebi_id, "CHEBI_", "") |> str_trim()
			chebi_vec[[k]] <- chebi_id
		} else {
			next
		}
	}
	chebi_vec <- chebi_vec[!is.na(chebi_vec)] |> unique()
	pmids[i,3] <- str_c(chebi_vec, collapse = "_")
}

## Transform the data
pmids_chebids <- pmids |> separate_longer_delim(chebid, delim="_") |> filter(chebid != "")
write_tsv(pmids_chebids, str_glue("{out_path}pmIDs_chebIDs__12-08-25.tsv"))
ids <- pmids_chebids |> select(chebid, category) |>
								  distinct() |>
								  mutate(category = if_else(category == "observational study", "obs", "meta")) |>
								  group_by(chebid) |>
								  summarise(category = c(category) |> unique() |> str_c(collapse = " & "))

## Add the chemical structures
structures_studies <- ids |> inner_join(structures_all, by = c("chebid" = "compound_id")) |> filter(!is.na(smiles)) |>
							select(chebid, smiles, category)
structures_fromStudies <- structures_studies |> select(chebid, smiles)

## Save the results
write_tsv(structures_studies, str_glue("{out_path}structures-&-studies__12-08-25.tsv"))
write_tsv(structures_fromStudies, str_glue("{out_path}structures_fromStudies__12-08-25.tsv"))