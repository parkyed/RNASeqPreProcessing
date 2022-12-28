# MA plots and Hoeffding's D-Statistic


# The MAcal_func function computes M and A matrices, while use the intensity of 20000 randomly selected probes
# Input is matrix of log2 normalised counts

MAcal_func <- function(x) { # matrix (row: ensembl gene IDs, col: array (samples)
  medArray = rowMedians(x, na.rm=TRUE)
  M =  x - medArray
  A = (x + medArray)/2
  subsample=20000
  if(nrow(M)>subsample) {
    set.seed(12345)
    sel = sample(nrow(M), subsample)
    sM = M[sel, ]
    sA = A[sel, ]
  } else {
    sM = M
    sA = A
  }
  output <- list(M=sM,A=sA) # return a list with M and A data matrices
  
  return(output)
}



# function to calculate D- statistic
d_stat_function <- function(logNormCount, targets, annotation, threshold) {
  
  # Compute the M and A matrices and assign to variables
  ma_list <- MAcal_func(logNormCount)
  M=ma_list$M
  A=ma_list$A
  
  # filter targets file based on the column name in annotation argument
  annotation <- enquo(annotation)
  ma_dstats_df <- dplyr::select(targets, analysisID, (!!annotation))
  if (!identical(colnames(M), ma_dstats_df$analysisID)) { stop() }
  
  # Compute the Hoeffding's statistic (Da) statistics for outlier detection
  Dstats =  sapply(1:ncol(M), function(x){hoeffd(A[,x], M[,x])$D[1,2]})
  ma_dstats_df$d_stat <- Dstats
  
  # identify outlier values that lie more than specified outlier threshold
  dstat_outlier_threshold <-  threshold
  ma_dstat_outliers <- filter(ma_dstats_df, d_stat > dstat_outlier_threshold)
  
  # bar chart of D statistic for each sample
  d_stat_barplot <- ggplot(ma_dstats_df, aes(x=reorder(analysisID, d_stat), y = d_stat, fill=get(names(ma_dstats_df)[2]))) +
    geom_bar(stat='identity') + coord_flip() +
    ggtitle("Hoeffding's d-statistic on joint dist. of M and A values") + 
    xlab("sampleID") + 
    ylab("Hoeffding's d-statistic") + 
    theme(legend.title = element_blank()) +
    geom_hline(yintercept = dstat_outlier_threshold)
  
  output <- list('d_stats_df' = ma_dstats_df, 'outliers' = ma_dstat_outliers, 'd_stat_bar' = d_stat_barplot)
  
  return(output)
  
}

# plot the MA plot for selected samples based on their d-statistic (top x, bottom x)
# select arrays with top x highest and bottom x lowest d-stat, and filter M and A matrices

ma_plots_function <- function(logNormCount, ma_dstats_dataframe, top_x, bottom_x) {
  
  # Compute the M and A matrices and assign to variables
  ma_list <- MAcal_func(logNormCount)
  M=ma_list$M
  A=ma_list$A
  
  # select the top x and bottom x rows of the M and A
  colnames(ma_dstats_dataframe)[2] <- "annotation"
  ma_stats_ordered <- ma_dstats_dataframe[order(-ma_dstats_dataframe$d_stat),]
  ma_stats_id_sel <- ma_stats_ordered[c(1:top_x,(ncol(M)-(bottom_x-1)):ncol(M)),]
  ma_stats_id_sel <- mutate(ma_stats_id_sel, group = c(rep('top', times=top_x), rep('bottom', times=bottom_x)))
  M_sel <- M[,ma_stats_id_sel$analysisID]
  A_sel <- A[,ma_stats_id_sel$analysisID]
  
  # create data frame for the plots
  ma_plot_df <- data.frame(
    sample_id=rep(ma_stats_id_sel$analysisID,each=nrow(M_sel)),
    annotation = rep(ma_stats_id_sel$annotation,each=nrow(M_sel)),
    group = factor(rep(ma_stats_id_sel$group,each=nrow(M_sel)), levels = c('top', 'bottom')),
    M=as.numeric(M_sel),
    A=as.numeric(A_sel)
  )
  
  # plots
  ma_plots <- ggplot(ma_plot_df,aes(x=A, y=M, color=annotation)) +
    geom_point(alpha=0.1) +
    theme_bw() +
    facet_wrap(~group + sample_id, ncol=5)
  
  return(ma_plots)
}



# gives MA plots of a subset of samples in a dataset (subset could equal entire set)
MAPlots <- function(exp.data, sample.set, n.col){
  
  # calculate M and A values
  m.a.values <- MAcal_func(exp.data)
  
  M.sel <- m.a.values$M[,sample.set]
  A.sel <- m.a.values$A[,sample.set]
  
  # create data frame for the plots
  ma.df <- data.frame(
    sample_id=rep(sample.set,each=nrow(M.sel)),
    M=as.numeric(M.sel),
    A=as.numeric(A.sel)
  )
  # MA plots
  ma.plots <- ggplot(ma.df,aes(x=A, y=M)) +
    geom_point(alpha=0.05) +
    theme_bw() +
    facet_wrap(~sample_id, ncol=n.col)
  
  return(ma.plots)
}
