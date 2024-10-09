library('tidyverse')
library('RColorBrewer')

#' Read the expression data "csv" file as a dataframe, not tibble
#'
#' @param filename (str): the path of the file to read
#' @param delimiter (str): generalize the function so it can read in data with
#'   your choice of delimiter
#'
#' @return A dataframe containing the example intensity data with rows as probes
#'   and columns as samples
#' @export
#'
#' @examples
read_data <- function(intensity_data, delimiter) {
  data <- read.table(file = intensity_data, sep = delimiter, header = TRUE)  
  
  return(data)
}

#' Define a function to calculate the proportion of variance explained by each PC
#'
#' @param pca_results (obj): the results returned by `prcomp()`
#'
#' @return A vector containing the values of the variance explained by each PC
#' @export
#'
#' @examples
calculate_variance_explained <- function(pca_results) {
  #extract st deviations and square to get variances
  variances <- pca_results$sdev^2
  
  #calculate the proportion of variance explained by each PC
  variance_proportion <- variances / sum(variances)
  
  return(variance_proportion)
}

#' Define a function that takes in the variance values and the PCA results to
#' make a tibble with PC names, variance explained by each PC, and the
#' cumulative sum of variance explained. These columns should be named 
#' "principal_components", "variance_explained", and "cumulative", respectively.
#' 
#'
#'
#' @param pca_ve (vector): the vector generated in the previous function with
#'   the variance explained values
#' @param pca_results (object): the results returned by `prcomp()`
#' @return A tibble that contains the names of the PCs, the individual variance
#'   explained, and the cumulative variance explained with names described above
#' @export
#' @examples 
make_variance_tibble <- function(pca_ve, pca_results) {
  
  #create PCs to use as factors for tibble
  len_components <- length(pca_ve)
  p_components <- paste0("PC", 1:len_components)
  
  #cumulative sum of variance explained
  cumulative_variance <- cumsum(pca_ve)
  
  #create tibble, use PCs as factors
  variance_tibble <- tibble(
    variance_explained = pca_ve,
    principal_components = factor(p_components),
    cumulative = cumulative_variance
  )
  
  return(variance_tibble)
}



#' Define a function to create a biplot of PC1 vs. PC2 labeled by
#' SixSubTypesClassification
#'
#' @param metadata (str): The path to the proj_metadata.csv file
#' @param pca_results (obj): The results returned by `prcomp()`
#'
#' @return A ggplot consisting of a scatter plot of PC1 vs PC2 labeled by
#'   SixSubTypesClassification found in the metadata
#' @export
#'
#' @examples
make_biplot <- function(metadata_path, pca_results) {
    
  metadata <- read_csv(metadata_path)
  
  pca_scores <- as_tibble(pca_results$x[, 1:2])  # PC1 and PC2
  pca_scores <- pca_scores %>%
    mutate(sample = rownames(pca_results$x))

  names(pca_scores) [names(pca_scores) == 'sample'] <- 'geo_accession'
  
  merged_data <- inner_join(pca_scores, metadata, by = 'geo_accession')

  biplot <- merged_data %>% ggplot() +
    geom_point(aes(x = PC1, y = PC2, color = SixSubtypesClassification)) +
    theme_minimal()

  return(biplot)
}

  

#' Define a function to return a list of probeids filtered by signifiance
#'
#' @param diff_exp_tibble (tibble): A tibble containing the differential expression results
#' @param fdr_threshold (float): an appropriate FDR threshold, we will use a
#'   value of .01. This is the column "padj" in the tibble.
#'
#' @return A list with the names of the probeids passing the fdr_threshold
#' @export
#'
#' @examples
list_significant_probes <- function(diff_exp_tibble, fdr_threshold) {
  # Filter the tibble for rows where padj is less than the FDR threshold
  significant_probes <- diff_exp_tibble %>%
    filter(padj < fdr_threshold) %>%
    select(probeid) # Select only the probeid column
  
  # Return the list of significant probe IDs
  return(significant_probes$probeid)
}

#' Define a function that uses the list of significant probeids to return a
#' matrix with the intensity values for only those probeids.
#' @param intensity (dataframe): The dataframe of intensity data generated in
#'   part 1
#' @param sig_ids_list (list/vector): The list of differentially expressed
#'   probes generated in part 6
#'
#' @return A `matrix()` of the probe intensities for probes in the list of
#'   significant probes by FDR determined in the previous function.
#'
#' @export
#'
#' @examples
return_de_intensity <- function(intensity, sig_ids_list) {
  # Subset the intensity matrix/data frame by matching the row names (probe IDs) with sig_ids_list
  de_intensity <- intensity[rownames(intensity) %in% sig_ids_list, ]
  
  # Convert the result to a matrix
  de_intensity_matrix <- as.matrix(de_intensity)
  
  # Return the subsetted matrix
  return(de_intensity_matrix)
}

#' Define a function that takes the intensity values for significant probes and
#' creates a color-blind friendly heatmap
#'
#' @param de_intensity (matrix): The matrix of intensity values for significant
#'   differentially expressed probes returned in part 7
#' @param num_colors (int): The number of colors in a specificed RColorBrewer
#'   palette
#' @param palette (str): The name of the chosen RColorBrewer palette
#'
#' @return A heatmap displaying the intensity values for the differentially
#'   expressed probes
#' @export
#'
#' @examples
plot_heatmap <- function(de_intensity, num_colors, palette) {
  # Generate the color palette using RColorBrewer
  color_palette <- brewer.pal(num_colors, palette)
  
  # Plot the heatmap using the intensity matrix and the specified color palette
  heatmap(de_intensity, col = color_palette, scale = "row")
  
  mtext(caption_8, side = 1, line = 5, outer = FALSE, cex = 1.2)
}

