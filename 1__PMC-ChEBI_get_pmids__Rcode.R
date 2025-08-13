library(tidyverse)
library(europmc)
library(jsonlite)


# europepmc facilitates access to the Europe PMC RESTful Web Service. The client furthermore supports the Europe PMC Annotations API to retrieve text-mined concepts and terms per article.
# Europe PMC covers life science literature and gives access to open access full texts. Europe PMC ingests all PubMed content and extends its index with other literature and patent sources.
# https://docs.ropensci.org/europepmc/

## Get the PubMed IDs of articles having associated ChEBI IDs
number_of_articles <- europepmc::epmc_hits("((ACCESSION_TYPE:'chebi') OR (HAS_CHEBI:y))") 
pmc_med_ids <- europepmc::epmc_search(query = "((ACCESSION_TYPE:'chebi') OR (HAS_CHEBI:y))", limit = number_of_articles)

# Filter the results
pmc_med_ids_processed <- pmc_med_ids |> filter(!is.na(pmid))

# Export the results
write_tsv(pmc_med_ids_processed, ".../summer_club/data/chebi_pmids__11-08-25.tsv")

## Get the ChEBI IDs of the interesting papers
# epmc_db(ext_id="40357648", data_src="med", db="CHEBI")  ---> NOT WORKING

# Use the other, fully official way: https://europepmc.org/annotationsapi#!/annotations45api45controller/getAnnotationsArticlesByProviderUsingGET