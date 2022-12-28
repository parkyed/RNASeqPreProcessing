
# Create a box plot over all samples, density plots for each class, calculate the K-S statistic for each samples and a horizontal bar
# chart plotting k-s statistic by sample. Identify outliers based on k-s statistic using Tukey's method.


box_density_k_s <- function(logNormCount, targets, annotation){
  
  #' Generate a box plot, density plots, calculate K-S statistic for each sample and show in a horizontal bar. Identify outliers based on k-s.
  #'
  #' @param logNormCount The matrix of log2 transformed, DESeq2 median-ratio normalised counts
  #' @param targets The targets dataframe
  #' @param annotation The factor in the targets dataframe used to annotate the chart
  
  # filter targets file based on the column name in annotation argument
  annotation <- enquo(annotation)
  targets.filtered <- dplyr::select(targets, analysisID, (!!annotation))
  if (!identical(colnames(logNormCount), targets.filtered$analysisID)) { stop() }
  
  # melt dataframe or normalised counts, add the annotation column
  logNormCount.melt <- as.data.frame(logNormCount) %>%               
    tibble::rownames_to_column("ensemblID")  %>%      
    melt(id.vars = 'ensemblID', value.name='count', variable.name=c('analysisID')) %>%  
    left_join(targets.filtered, by = 'analysisID')
  
  # boxplot of the distribution of each sample
  dist_box_plot <- ggplot(logNormCount.melt, aes(x=reorder(analysisID, count, FUN = median), y = count, fill=get(names(logNormCount.melt)[4]))) + 
    geom_boxplot() + 
    ggtitle("Box Plot of Distribution of log2 normalised counts")+ 
    xlab("sampleID") + 
    coord_flip() +
    theme(legend.title = element_blank()) +
    ylab("log[2](count + 1)")
  
  # density plot of the distribution of each sample, split out by annotation
  dist_density_plot <- ggplot(logNormCount.melt, aes(x = count, fill=analysisID)) + 
    geom_density(alpha = 0.2, size = 1.25) + 
    ggtitle("Density Plot of log2 normalised counts") + 
    xlab("counts") + 
    ylab("density") + 
    facet_wrap(~ get(names(logNormCount.melt)[4])) + 
    theme(legend.position="none")
  
  # pool all the log counts
  pooled_logNormCount = c(logNormCount)                                             
  
  # loop over every column in the log2Normalised counts and perform k-s test
  # matricise this using purrr - get rid of the for loop!!
  k_stats <- character()
  for (i in 1:nrow(targets.filtered)){
    k_stat <- ks.test(logNormCount[,i], pooled_logNormCount, exact = FALSE)
    k_stats[i] <- k_stat$statistic
  }
  targets.filtered$k_stat <- as.numeric(k_stats)
  
  # identify outlier values that lie more than 1.5IQR beyond 25th or 75th percentiles
  ks_outlier_threshold <- boxplot.stats(targets.filtered$k_stat)$stats[5]
  dist_outliers <- filter(targets.filtered, k_stat > ks_outlier_threshold)
  
  # bar chart of D statistic for each sample
  dist_outlier_plot <- ggplot(targets.filtered, aes(x=reorder(analysisID, k_stat), y = k_stat, fill=get(names(logNormCount.melt)[4]))) +
    geom_bar(stat='identity') +
    ggtitle("Bar plot of k-s statistic by sample") + 
    xlab("sampleID") + 
    ylab("k-s statistic") + 
    coord_flip() +
    theme(legend.title = element_blank()) +
    geom_hline(yintercept = ks_outlier_threshold)
  
  output <- list('box_plot'= dist_box_plot, 'density_plot' = dist_density_plot, 'outliers' = dist_outliers, 'k_stat_bar' = dist_outlier_plot)
  
  return(output)
}