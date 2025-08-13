library(tidyverse)
library(umap)
library(rJava)
library(rcdklibs)
library(rcdk)
library(fingerprint)
library(praznik)
library(mirt)

# Input
path <- ".../summer_club/types_of_studies_14-08_2025/output/"
data <- read_tsv(str_glue("{path}structures-&-studies__12-08-25.tsv")) |> mutate(parseable = "no") |> rename(type = category)

# Check SMILES for parseability
maccs_fps <- rep(NA, nrow(data)) |> as.list()
for (i in seq(1:nrow(data))) {
	mol <- parse.smiles(data[i,2] |> pull())
	if(!is.null(mol[[1]])) {
		data[i,4] <- "yes"
		maccs_fp <- get.fingerprint(mol[[1]], type = 'maccs')
		maccs_fps[[i]] <- maccs_fp
	}
}

# Filter to avoid the non-parseable structures
data <- data |> filter(parseable != "no")
maccs_fps <- maccs_fps[!is.na(maccs_fps)]

# FPs to matrix
maccs_matrix <- fp.to.matrix(maccs_fps)
maccs_int <- as.integer(maccs_matrix)
maccs_matrix_int <- matrix(maccs_int, nrow = nrow(maccs_matrix), ncol = ncol(maccs_matrix))

# ADD FPs to the data
data_fps <- bind_cols(data, maccs_matrix_int) |>
					mutate(y = row_number()) |>
					rename_with(~ paste0('maccs_', seq_along(.x)), starts_with('...')) |>
					select(-parseable)

# Select 11 features
bits_selection <- MRMR( data_fps |> select(-chebid, -smiles, -type, -y), data_fps |> pull(y), 11 )
bits_selected  <- names(bits_selection$selection)
data_fps_selected <- data_fps |> select(chebid, smiles, type, any_of(bits_selected)) |>
						group_by(maccs_139, maccs_145, maccs_151, maccs_147, maccs_148, maccs_150, maccs_132, maccs_112, maccs_146, maccs_153, maccs_137) |>
						summarise(type = c(type) |> unique() |> str_c(collapse = " & "),
									chebid = c(chebid) |> unique() |> str_c(collapse = " & ")) |>
						rowwise() |>
						mutate(type = type |> str_split_1(" & ") |> str_sort() |> unique() |> str_c(collapse = " & ")) |>
						ungroup()

# UMAP the selected bits and plot
set.seed(23)
data_bits <- data_fps_selected |> select( any_of(bits_selected) )
data_labels   <- data_fps_selected |> select(chebid, type)
umap <- umap(data_bits, n_neighbors = 11, metric = "manhattan", preserve.seed = TRUE)
coordinates	<- umap$layout |> bind_cols(data_labels) |>
							  rename(coord_one = `...1`, coord_two = `...2`) |>
							  mutate(order = NA) |>
							  mutate(order = if_else(type == "unknown", 1, order)) |>
							  mutate(order = if_else(type == "meta & obs", 2, order)) |>
							  mutate(order = if_else(type == "obs", 3, order)) |>
							  mutate(order = if_else(type == "meta", 4, order)) |>
							  mutate(type = fct_reorder(type, order))

umap_plot  <- ggplot(coordinates |> arrange(order), aes(x = coord_one, y = coord_two, color = type)) +
															 geom_point() +
															 scale_color_manual(values = c('#ffb588', '#AB6C82', '#475C7A')) +
															 coord_fixed(ratio = 1, xlim = c(-4, 4), ylim = c(-4, 4)) + 
															 theme_minimal() +
															 theme(plot.background = element_rect(fill = "white", color = "white"))
umap_plot
ggsave(str_glue("{path}basicUMAP.png"), width = 6, units = "in", dpi = 300)

# UMAP the whole space of selected bits and plot
set.seed(43)
#Enumerate feature vectors, target dimensionality is 11
featureSpace_11_raw <- thetaComb(theta = c(0L,1L), nfact = 11, intercept = FALSE)
colnames(featureSpace_11_raw) <- names(bits_selection$selection)
featureSpace_11 <- as_tibble(featureSpace_11_raw) |>
						mutate(chebid = str_glue("no_{row_number()}"),
								type = "unknown",
						full_description = str_c(maccs_139, maccs_145, maccs_151, maccs_147, maccs_148, maccs_150, maccs_132, maccs_112, maccs_146, maccs_153, maccs_137, sep = "") )
data_fps_selected <- data_fps_selected |> mutate(full_description = str_c(maccs_139, maccs_145, maccs_151, maccs_147, maccs_148, maccs_150, maccs_132, maccs_112, maccs_146, maccs_153, maccs_137, sep = ""))
space_points <- featureSpace_11 |> anti_join(data_fps_selected, by = "full_description")
# |> sample_frac(size = .3)
all_points <- bind_rows(data_fps_selected, space_points)
#UMAP
extended_data_bits <- all_points |> select( any_of(bits_selected) )
extended_data_labels   <- all_points |> select(chebid, type)
extended_umap <- umap(extended_data_bits, n_neighbors = 11, metric = "manhattan", preserve.seed = TRUE)
extended_coordinates	<- extended_umap$layout |> bind_cols(extended_data_labels) |> rename(coord_one = `...1`, coord_two = `...2`) |>
							  mutate(order = NA) |>
							  mutate(order = if_else(type == "unknown", 1, order)) |>
							  mutate(order = if_else(type == "meta & obs", 2, order)) |>
							  mutate(order = if_else(type == "obs", 3, order)) |>
							  mutate(order = if_else(type == "meta", 4, order)) |>
							  mutate(type = fct_reorder(type, order))
extended_umap_plot  <- ggplot(extended_coordinates |> arrange(order), aes(x = coord_one, y = coord_two, color = type)) +
															 geom_point() +
															 scale_color_manual(values = c('#CBCBD6FF', '#ffb588', '#AB6C82', '#475C7A')) +
															 coord_fixed(ratio = 1, xlim = c(-4, 4), ylim = c(-4, 4)) + 
															 theme_minimal() +
															 theme(plot.background = element_rect(fill = "white", color = "white"))
extended_umap_plot
ggsave(str_glue("{path}extendedUMAP.png"), width = 6, units = "in", dpi = 300)

# UMAP the actual points along with the random fraction of space points and plot
set.seed(63)
space_points_frac <- featureSpace_11 |> anti_join(data_fps_selected, by = "full_description") |>
									sample_frac(size = .3)
# |> sample_frac(size = .3)
many_points <- bind_rows(data_fps_selected, space_points_frac)
#UMAP
frac_data_bits <- many_points |> select( any_of(bits_selected) )
frac_data_labels   <- many_points |> select(chebid, type)
frac_umap <- umap(frac_data_bits, n_neighbors = 11, metric = "manhattan", preserve.seed = TRUE)
frac_coordinates	<- frac_umap$layout |> bind_cols(frac_data_labels) |> rename(coord_one = `...1`, coord_two = `...2`) |>
							  mutate(order = NA) |>
							  mutate(order = if_else(type == "unknown", 1, order)) |>
							  mutate(order = if_else(type == "meta & obs", 2, order)) |>
							  mutate(order = if_else(type == "obs", 3, order)) |>
							  mutate(order = if_else(type == "meta", 4, order)) |>
							  mutate(type = fct_reorder(type, order))
frac_umap_plot  <- ggplot(frac_coordinates |> arrange(order), aes(x = coord_one, y = coord_two, color = type)) +
															 geom_point() +
															 scale_color_manual(values = c('#CBCBD6FF', '#ffb588', '#AB6C82', '#475C7A')) +
															 coord_fixed(ratio = 1, xlim = c(-4, 4), ylim = c(-4, 4)) + 
															 theme_minimal() +
															 theme(plot.background = element_rect(fill = "white", color = "white"))
frac_umap_plot
ggsave(str_glue("{path}fracUMAP.png"), width = 6, units = "in", dpi = 300)