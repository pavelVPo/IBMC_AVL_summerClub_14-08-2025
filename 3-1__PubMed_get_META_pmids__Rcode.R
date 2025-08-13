library(tidyverse)
library(rentrez)

# Key
api_key = "af82cca6d1bd01f71dc4319bf810722a2708"
#C72C41 In order not to overload the E-utility servers, NCBI recommends that users post
#C72C41 no more than three URL requests per second and
#C72C41 limit large jobs to either weekends or between 9:00 PM and 5:00 AM Eastern time during weekdays.
#C72C41 Failure to comply with this policy may result in an IP address being blocked from accessing NCBI.
# Output path
out_path <- ".../summer_club/types_of_studies_14-08_2025/data/"

## Meta-Analyses
# Initial search, SEE: https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html#searching-databases-entrez_search
meta_r_search <- entrez_search(db="pubmed", term="(study) AND (meta-analysis[Publication Type])", retmax = 1,
								api_key = "af82cca6d1bd01f71dc4319bf810722a2708", use_history = FALSE)
# URL with the data (PMIDs)
# (study) AND (metaervational Study[Publication Type]) -> term=%28study%29+AND+%28metaervational+Study%5BPublication+Type%5D%29
# (study) AND (Meta-Analysis[Publication Type])		  -> term=%28study%29+AND+%28Meta-Analysis%5BPublication+Type%5D%29
meta_base_url <- "https://pubmed.ncbi.nlm.nih.gov/?term=%28study%29+AND+%28meta-analysis%5BPublication+Type%5D%29&api_key=af82cca6d1bd01f71dc4319bf810722a2708"
meta_first_url <- str_glue("{meta_base_url}&sort=pubdate&sort_order=asc&size=10")
meta_last_url <- str_glue("{meta_base_url}&sort=pubdate&sort_order=desc&size=10")
# Number of studies
total_count <- meta_r_search$count |> as.integer()
meta_ids <- rep(NA, total_count)
# Get the year of the first publication
search_rslt <- read_html(meta_first_url)
first_rslt <- xml_find_first(search_rslt, ".//div[contains(@class, 'search-results-chunk results-chunk')]")
first_id <- xml_attr(first_rslt, "data-chunk-ids") |> str_replace(",.*", "") |> as.integer()
first_year <- entrez_summary(db = "pubmed", id = first_id, retmax = 1, api_key = "af82cca6d1bd01f71dc4319bf810722a2708")$pubdate |> str_replace(" .*", '')
# Use date intervals to search further retriving the data for each month
year = as.integer(first_year)
month = 0
current_position = 0
repeat {
	# Prepare
	monthly_ids <- "empty"
	# Last month search
	if (month == 12) {
		monthly_meta_r_search <- tryCatch({ entrez_search(db="pubmed", term=str_glue("(study) AND (('{year}/{month}/1'[Date - Publication] : '{year}/{month}/31'[Date - Publication])) AND (meta-analysis[Publication Type])"),
														retmax = 9999, api_key = "af82cca6d1bd01f71dc4319bf810722a2708", use_history = FALSE) },
										warning = function(w) { "empty" },
										error = function(e) { "empty" })
	}
	# General search
	if (month != 12) {
		monthly_meta_r_search <- tryCatch({ entrez_search(db="pubmed", term=str_glue("(study) AND (('{year}/{month}/1'[Date - Publication] : '{year}/{month+1}/1'[Date - Publication])) AND (meta-analysis[Publication Type])"),
														retmax = 9999, api_key = "af82cca6d1bd01f71dc4319bf810722a2708", use_history = FALSE) },
										warning = function(w) { "empty" },
										error = function(e) { "empty" })
	}
	# Collect the IDs
	if (is.list(monthly_meta_r_search)) {
		monthly_ids <- monthly_meta_r_search$ids
		if (length(monthly_ids) > 0) {
			for (i in seq(1:length(monthly_ids))) {
				meta_ids[current_position + 1] <- monthly_ids[i]
				current_position <- current_position + 1
			}
		}
	}
	# Exit
	if (year > 2025) {break}
	# Change month
	month = month + 1
	# Change year
	if (month == 13) {
		year <- year + 1
		month <- 1}
}
# Process and save the results
meta_ids_t <- tibble(pmid = meta_ids) |> filter(!is.na(meta_ids)) |> distinct()
write_tsv(meta_ids_t, str_glue("{out_path}meta_ids__12-08-2025.tsv"))