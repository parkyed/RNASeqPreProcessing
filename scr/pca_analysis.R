# script to implement pca analysis and plots based on annotation, gene signatures, and included samples
# flexibility to colour plots based on any phenotype variable, and to reduce the set of genes and samples used
# enables sex check, and classification check with alternative gene signatures


# function to filter the dataset based on a specified gene signature, set of examples and perform the PCA analysis

calc_pca_results_gene_sig <-  function(log2_counts, targets, gene_superset, g_signature, c_inc, n_pc){

  #' Calculate principal component eigenvalues
  #' 
  #' This function wraps the prcomp function giving flexibility to filter the input data set by a specific gene signature, and set of conditions to include
  #'
  #' @param log2_counts the matrix of counts, either vst or rlog transformed
  #' @param targets the targets matrix
  #' @param gene_superset the superset of genes to filter from, e.g. all genes, or just Y chromosome genes
  #' @param g_signature the gene names of the genes we want to filter on
  #' @param c_inc the conditions to include, e.g. N-LOS, P-LOS etc.
  #' @param n_pc the number of principal components to calculate

  # selected samples to compare based on the included conditions
  pca_selected_samples <- dplyr::filter(targets, condition %in% c_inc)$analysisID
  
  # filter genes list by selected gene names to give subset of ensembl IDs
  pca_gene_signature <- dplyr::filter(gene_superset, Gene.name %in% g_signature)$Gene.stable.ID
  
  ## filter the log2 counts dataset by selected genes and samples
  pca_gs_counts <- as.data.frame(log2_counts)
  pca_gs_counts <- pca_gs_counts %>% filter(row.names(pca_gs_counts) %in% pca_gene_signature) %>% dplyr::select(all_of(pca_selected_samples))
  
  ## calculate principle components, using the transpose of the counts matrix.
  gs_pca_results <- prcomp(t(pca_gs_counts), center = TRUE, scale = FALSE, rank = n_pc)
  
  ## check PCA output matrix has samples in the same order as the original targets file
  if (!identical(rownames(gs_pca_results$x), pca_selected_samples)) { stop() }
  
  return(gs_pca_results)
}


# function to create scree plot

pca_scree_plot <- function(pca_results){
  
  #' Generate a scree plot showing the % of variation in dataset explained by each principal component
  #'
  #' @param pca_results the matrix of eigen values output from the PCA cacluation
  
  ## compute the variance explained by each PC for the n PCs calculated
  pca_var_explained <- (pca_results$sdev^2/sum(pca_results$sdev^2))[1:ncol(pca_results$x)]
  
  ## Scree plot of the proportion of variance explained by each PC
  scree_plot <- data.frame('principal_component'=colnames(pca_results$x), 'variance_explained'= pca_var_explained) %>%
    ggplot(aes(x=principal_component, y=pca_var_explained)) +
    geom_bar(stat='identity') +
    ggtitle("Scree Plot - % Variance Explained")
  
  return(scree_plot)
  
}

## annotated PCA plot function
pca_scatter_annotated <- function(pca_results, targets, annotation) {
  
  #' Generate a 2D scatter plot of the first two PCs, annotated with a given column in the phenotype data
  #'
  #' @param pca_results the matrix of eigen values output from the PCA cacluation
  #' @param targets the targets matrix
  #' @param annotation The the name of the column in the targets dataframe used to annotate plot
  
  annotation <- enquo(annotation)  
  
  ## compute the variance explained by each PC for the n PCs calculated
  pca_var_explained <- (pca_results$sdev^2/sum(pca_results$sdev^2))
  
  ## Scatter plots of pairs of principle components
  pca_results_merged <- pca_results$x %>%
    as_tibble(rownames = "analysisID") %>%
    merge(targets, by="analysisID",  sort=FALSE)    # may need to left join this...
    
  pca_plot <- ggplot(pca_results_merged, aes(x=PC1,y=PC2)) +
    geom_point(aes(color=(!!annotation)),size=4) +
    theme_bw(base_size=32) +
    labs(x=paste0("PC1: ",round(pca_var_explained[1]*100,1),"%"),
         y=paste0("PC2: ",round(pca_var_explained[2]*100,1),"%")) +
    theme(legend.position="bottom",
          legend.margin=ggplot2::margin(-5,-5,-5,-5),
          legend.box.margin=ggplot2::margin(-5,-5,-5,-5),
          axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          plot.title=element_text(size=14),
          legend.text=element_text(size=10),
          legend.title=element_text(size=10)) +
    # scale_colour_discrete(drop=TRUE, limits = levels(pca_results_merged$annotation)) + # need to figure out how to ensure label colours are consistent
    # the colours need to be defined by the factors of the annotation column, rather than explicitly specified
    # maybe use limits=levels(get(names(xxx))
    # scale_colour_manual(values = plot_colours, name = "Patient Condition:") +
    geom_text_repel(aes(label = analysisID),
                    size = 3,
                    box.padding   = 0.15,
                    point.padding = 0.25,
                    segment.color = 'grey50')
  
  # excluded from function
  # stat_ellipse(level = 0.85) +
  
  return(pca_plot)
}

## annotated PCA plot function, with flexible labelling of the points

pcaScatterAnnoLabel <- function(pca_results, targets, annotation, label) {
  
  #' Generate a 2D scatter plot of the first two PCs, annotated with a given column in the phenotype data
  #'
  #' @param pca_results the matrix of eigen values output from the PCA cacluation
  #' @param targets the targets matrix
  #' @param annotation The the name of the column in the targets dataframe used to annotate plot
  
  annotation <- enquo(annotation)  
  label <- enquo(label)  
  
  ## compute the variance explained by each PC for the n PCs calculated
  pca_var_explained <- (pca_results$sdev^2/sum(pca_results$sdev^2))
  
  ## Scatter plots of pairs of principle components
  pca_results_merged <- pca_results$x %>%
    as_tibble(rownames = "analysisID") %>%
    merge(targets, by="analysisID",  sort=FALSE)    # may need to left join this...
  
  pca_plot <- ggplot(pca_results_merged, aes(x=PC1,y=PC2)) +
    geom_point(aes(color=(!!annotation)),size=4) +
    theme_bw(base_size=32) +
    labs(x=paste0("PC1: ",round(pca_var_explained[1]*100,1),"%"),
         y=paste0("PC2: ",round(pca_var_explained[2]*100,1),"%")) +
    theme(legend.position="bottom",
          legend.margin=margin(-5,-5,-5,-5),
          legend.box.margin=margin(-5,-5,-5,-5),
          axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          plot.title=element_text(size=14),
          legend.text=element_text(size=10),
          legend.title=element_text(size=10)) +
    # scale_colour_discrete(drop=TRUE, limits = levels(pca_results_merged$annotation)) + # need to figure out how to ensure label colours are consistent
    # the colours need to be defined by the factors of the annotation column, rather than explicitly specified
    # maybe use limits=levels(get(names(xxx))
    # scale_colour_manual(values = plot_colours, name = "Patient Condition:") +
    geom_text_repel(aes(label = (!!label)),
                    size = 3,
                    box.padding   = 0.15,
                    point.padding = 0.25,
                    segment.color = 'grey50')
  
  # excluded from function
  # stat_ellipse(level = 0.85) +
  
  return(pca_plot)
}


#plot_colours <- c("C-LOS" = "red", "P-LOS" = "blue", "N-LOS" = "darkgreen")

