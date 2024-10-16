# Function that computes geometric mean
GmMean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
}

# Clr transformation
ClrFunction = function(x, base=2){
  x <- log((x / GmMean(x)), base)
  fin <- x[is.finite(x)]
  fin <- fin[!is.na(fin)]
  x[!is.finite(x) | is.na(x)] <- min(fin)
  return(x)
}

# Function to plot a PCOA for a given set of axes and save in a given save directory. 
# "Type" indicates the type of pcoa (ex "Euclidean", "Unifrac")
plot_pcoa <- function(po_in, pcoa_in, axes = 1:2, fill_var = "Treatment_Cycle", colorscale = NULL, 
                      type = "PCoA", savedir, suffix = "", connect = FALSE, ptcolor = NULL, rmlegend = FALSE, 
                      h = 5, w = 5.5, title = TRUE, connectmouse = FALSE, diffshape = FALSE, returnplot = FALSE,
                      facetDayPt = FALSE, facetDay = FALSE, ncolDay = 6, facetPt = FALSE, facetTxCat = FALSE, mousetxshape = FALSE,
                      sictxshape = FALSE){
  
  if (!diffshape){
  p <- plot_ordination(po_in, pcoa_in, axes = axes) +
    geom_point(aes(fill = !!sym(fill_var)), size = 3, stroke=0, pch = 21) + # size = 'Day' if interested
    scale_fill_manual(values = colorscale, name = fill_var) +
    theme_pubr() 
  } else if (mousetxshape) {
    p <- plot_ordination(po_in, pcoa_in, axes = axes, shape = "Donor") + 
      geom_point(size = 2, aes(shape = Donor, color = !!sym(fill_var))) + # size = 'Day' if interested
      scale_color_manual(values = colorscale, name = "Treatment") +
      theme_pubr() + 
      theme(legend.position = "right")
  } else if (sictxshape) {
    p <- plot_ordination(po_in, pcoa_in, axes = axes, shape = "Patient_ID") + 
      geom_point(size = 2, aes(shape = Patient_ID, color = !!sym(fill_var))) + # size = 'Day' if interested
      scale_color_manual(values = colorscale, name = "Treatment") +
      theme_pubr() + 
      theme(legend.position = "right")
  } else {
    p <- plot_ordination(po_in, pcoa_in, axes = axes) +
      geom_point(aes(shape = "Patient_ID", color = !!sym(fill_var))) + # size = 'Day' if interested
      scale_color_manual(values = colorscale, name = fill_var) +
      theme_pubr() 
  }
  
  axes_str = paste0(axes, collapse = "")
  fn <- paste0(savedir, "Beta_", type, "_PCA_axes", axes_str, suffix, ".pdf")
  if (connect){
    p <- p + geom_path(aes(group=Patient_ID, colour=Patient_ID), alpha=0.4, arrow = arrow(length = unit(0.3, "cm"))) + 
      scale_colour_manual(values = ptcolor) + 
      guides(colour = "none")
    fn <- paste0(savedir, "Beta_", type, "_PCA_axes", axes_str, suffix, "_ConnectedPt.pdf")
  }
  if (connectmouse){
    newtab = data.table(p$data)
    setorder(newtab, Day)
    p$data <- newtab
    
    p <- p + geom_path(aes(group=MouseID, colour=MouseID), alpha=0.4, arrow = arrow(length = unit(0.3, "cm"))) + 
      scale_colour_manual(values = ptcolor) + 
      guides(colour = "none")
    fn <- paste0(savedir, "Beta_", type, "_PCA_axes", axes_str, suffix, "_ConnectedPt.pdf")
  }
  if (title){
    p <- p + ggtitle(paste0(type, ' PCoA')) 
  }
  if (rmlegend){
    p <- p + theme(legend.position = 'none')
  }
  if (facetDayPt){
    p <- p + facet_grid(Donor ~ Day, scales = 'free')
  }
  if (facetDay){
    p <- p + facet_wrap( ~ Day, ncol = ncolDay)
  }
  if (facetPt){
    p <- p + facet_wrap( ~ Patient_ID, scales = "free")
  }
  if (facetTxCat){
      p <- p + facet_wrap( ~ TreatmentCyclev1)
  }
  if (returnplot){
    return(p)
  } else{
    fn <- gsub("-", "", fn)
    ggsave(filename = fn, plot = p, useDingbats=F, height=h, width=w)
  }
}

## Same as above, but with baseline normalization
# Function to plot a PCOA for a given set of axes and save in a given save directory. 
# "Type" indicates the type of pcoa (ex "Euclidean", "Unifrac")
plot_pcoa_v2 <- function(po_in, pcoa_in, axes = 1:2, fill_var = "Treatment_Cycle", colorscale = NULL, 
                         type = "PCoA", savedir, suffix = "", connect = FALSE, ptcolor = NULL, rmlegend = FALSE, 
                         returnplot = FALSE){
  dtest <- plot_ordination(po_in, pcoa_in, axes = axes, justDF = TRUE)

  ax1name <- paste0("PC", axes[1])
  ax2name <- paste0("PC", axes[2])
  if (type %in% c("Bray", "Bray-Curtis")){
    ax1name <- paste0("Axis.", axes[1])
    ax2name <- paste0("Axis.", axes[2])
  }
  dtest$PC1 <- dtest[[ax1name]]
    dtest$PC2 <- dtest[[ax2name]]
    
    rownams <- rownames(dtest)
    nams <- gsub(".*GO-", "", rownams)
    nams <- gsub("-.*", "", nams)
    nams <- gsub(".*Pt", "", nams)
    dtest$Patient_ID <- as.factor(nams)
    
    dtest$Treatment_Cycle <- gsub(".*-", "", rownams)
    dtest$Treatment_Cycle[dtest$Treatment_Cycle == "C1D7"] <- "C1Mid"
    dtest$Treatment_Cycle <- as.factor(dtest$Treatment_Cycle)
  
  dtest2 <- dtest %>% group_by(Patient_ID) %>%
    mutate(PC1_base = PC1[which.min(Treatment_Cycle)]) %>%
    group_by(Patient_ID) %>% 
    mutate(PC2_base = PC2[which.min(Treatment_Cycle)]) 
  dtest2 <- dtest2 %>% mutate(PC1_norm = PC1 - PC1_base) %>% 
    mutate(PC2_norm = PC2 - PC2_base)
  
  p <- ggplot(dtest2, aes(x = PC1_norm, y = PC2_norm)) + 
    geom_point(aes(fill = !!sym(fill_var)), size=3, stroke=0, pch = 21) +
    scale_fill_manual(values = colorscale, name = fill_var) +
    xlab(paste0("Normalized PC", axes[1])) + 
    ylab(paste0("Normalized PC", axes[2])) + 
    theme_pubr() + 
    ggtitle(paste0(type, ' PCoA')) 
  
  axes_str = paste0(axes, collapse = "")
  fn <- paste0(savedir, type, "_PCA_axes", axes_str, suffix, "_BaselineNrml.pdf")
  if (connect){
    p <- p + geom_path(aes(group=Patient_ID, colour=Patient_ID), alpha=0.4, arrow = arrow(length = unit(0.3, "cm"))) + 
      scale_colour_manual(values = ptcolor) + 
      guides(colour = "none")
    fn <- paste0(savedir, type, "_PCA_axes", axes_str, suffix, "_ConnectedPt_BaselineNrml.pdf")
  }
  if (rmlegend){
    p <- p + theme(legend.position = 'none')
  }
  
  if (returnplot){
    return(p)
  } else{
    fn <- gsub("-", "", fn)
    ggsave(filename = fn, plot = p, useDingbats=F, height=5, width=5.5)
  }
}


# Function to plot taxa. w = width, h = height
# Can set "mode" to either "normal" (standard axes) or "logplot" (scale_y_log10())
# Plots only taxa that are nonzero in >1/2 of samples unless a specific taxa is selected to plot with toi
plot_taxa <- function(TaxaSummary, taxa_level = "Phylum", taxa_letter = "p", metadata_in, savedir, w = 6, h = 5, 
                      toi = NULL, suffix = "", mode = "normal", clr = FALSE, xvar = "Treatment_Cycle", 
                      xlabstr = 'Treatment cycle', facettx = FALSE){
  dat <- TaxaSummary[[taxa_level]] # Extract taxa level data
  # Fix row names
  pattern = paste0(".*", taxa_letter, "__")
  dat <- dat %>% as.data.frame()
  dat$NewNames <- gsub(pattern, "", rownames(dat))
  dat$NewNames <- case_when(!duplicated(dat$NewNames) ~ dat$NewNames,
                            TRUE ~ rownames(dat))
  rownames(dat) <- dat$NewNames 
  dat <- dat %>% select(-c('NewNames')) 
  if (clr){ # clr normalize (optional)
    newmin <- 0.65*min(dat[dat > 0])
    dat[dat == 0] <- newmin; dat[is.na(dat)] <- newmin
    dat <- clr(t(dat))
  }
  else{
    dat <- t(dat)/rowSums(t(dat)) # Normalize by abundance
  }
  
  # Filter for taxa of interest; if not filtering, remove all taxa that are zero in > half of samples
  if (length(toi) == 0){
    dat <- dat[, names(colSums(dat > 0)[colSums(dat > 0) > length(rownames(dat))/2])] # Filter out taxa with low reads (nonzero in at least one half of all samples)
  }
  dat <- melt(dat) %>% as.data.frame() # Melt
  
  # Fix names
  colnames(dat) <- c('Sample_ID', 'Taxa', 'Value')
  df_taxa <- left_join(dat, metadata_in)
  
  # Filter
  if (length(toi) > 0){
    df_taxa <- df_taxa %>% filter(Taxa %in% toi)
  }
  
  minval2 = 0
  if (mode == "logplot"){
    minval2 = min(df_taxa$Value[df_taxa$Value > 0]) # get second smallest value in dataframe for log plot
  }
  p <- ggplot(df_taxa, aes(x = .data[[xvar]], y = Value + minval2, fill = Taxa)) + 
    geom_jitter(alpha = 0.2, width = 0.2) + 
    stat_summary(geom='ribbon', 
                 fun.data = mean_se, 
                 #fun.args=list(conf.int=0.95),
                 alpha = 0.5,
                 aes(fill = Taxa, group = Taxa)) + 
    theme_pubr() + 
    ylab('Relative abundance (proportion)') + 
    xlab(xlabstr) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(legend.position = 'none')
  if (mode == "logplot"){ 
    p <- p + scale_y_log10() # log y axis
  }
  if (clr){
    p <- p + ylab('Abundance (CLR)')
  }
  if (facettx){
    p <- p + facet_grid(Treatment ~ Taxa, scales ="free")
  } else {
    p <- p + facet_wrap(~ Taxa, scales ="free")
  }
  fn <- paste0(savedir, taxa_level, "VsTime", suffix, ".pdf")
  ggsave(plot = p, filename = fn, width = w, height = h, units = "in")
}

# Similar to above, but for coloring based on Tx
plot_taxa_tx <- function(TaxaSummary, taxa_level = "Phylum", taxa_letter = "p", metadata_in, savedir, w = 6, h = 5, 
                      toi = NULL, suffix = "", mode = "normal", clr = FALSE){
  dat <- TaxaSummary[[taxa_level]] # Extract taxa level data
  # Fix row names
  pattern = paste0(".*", taxa_letter, "__")
  dat <- dat %>% as.data.frame()
  dat$NewNames <- gsub(pattern, "", rownames(dat))
  dat$NewNames <- case_when(!duplicated(dat$NewNames) ~ dat$NewNames,
                            TRUE ~ rownames(dat))
  rownames(dat) <- dat$NewNames 
  dat <- dat %>% select(-c('NewNames')) 
  if (clr){ # clr normalize (optional)
    newmin <- 0.65*min(dat[dat > 0])
    dat[dat == 0] <- newmin; dat[is.na(dat)] <- newmin
    dat <- clr(t(dat))
  }
  else{
    dat <- t(dat)/rowSums(t(dat)) # Normalize by abundance
  }
  
  # Filter for taxa of interest; if not filtering, remove all taxa that are zero in > half of samples
  if (length(toi) == 0){
    dat <- dat[, names(colSums(dat > 0)[colSums(dat > 0) > length(rownames(dat))/2])] # Filter out taxa with low reads (nonzero in at least one half of all samples)
  }
  dat <- melt(dat) %>% as.data.frame() # Melt
  
  # Fix names
  colnames(dat) <- c('Sample_ID', 'Taxa', 'Value')
  df_taxa <- left_join(dat, metadata_in)
  
  # Filter
  if (length(toi) > 0){
    df_taxa <- df_taxa %>% filter(Taxa %in% toi)
  }
  
  minval2 = 0
  if (mode == "logplot"){
    minval2 = min(df_taxa$Value[df_taxa$Value > 0]) # get second smallest value in dataframe for log plot
  }
  p <- ggplot(df_taxa, aes(x = Treatment_Cycle, y = Value + minval2, fill = Treatment)) + 
    geom_jitter(aes(color = Treatment), alpha = 0.2, width = 0.2) + 
    #geom_line(aes(group=Patient_ID, color = Treatment), alpha = 0.2) +
    stat_summary(geom='ribbon', 
                 fun.data = mean_se, 
                 #fun.args=list(conf.int=0.95),
                 alpha = 0.5,
                 aes(fill = Treatment, group = Treatment)) + 
    facet_grid(Treatment ~ Taxa, scales ="free") + 
    theme_pubr() + 
    ylab('Relative abundance (proportion)') + 
    xlab('Treatment cycle') + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values = treatment.colours) + 
    scale_color_manual(values = treatment.colours) + 
    theme(legend.position = 'none')
  if (mode == "logplot"){ 
    p <- p + scale_y_log10() # log y axis
  }
  if (clr){
    p <- p + ylab('Abundance (CLR)')
  }
  fn <- paste0(savedir, taxa_level, "VsTime", suffix, "_Treatment.pdf")
  ggsave(plot = p, filename = fn, width = w, height = h, units = "in")
}

# Version 2
# Similar to above, but for coloring based on Tx
plot_taxa_tx_v2 <- function(TaxaSummary, taxa_level = "Phylum", taxa_letter = "p", metadata_in, savedir, w = 6, h = 5, 
                         toi = NULL, suffix = "", mode = "normal", clr = FALSE, pt.colors = NULL){
  dat <- TaxaSummary[[taxa_level]] # Extract taxa level data
  # Fix row names
  pattern = paste0(".*", taxa_letter, "__")
  dat <- dat %>% as.data.frame()
  dat$NewNames <- gsub(pattern, "", rownames(dat))
  dat$NewNames <- case_when(!duplicated(dat$NewNames) ~ dat$NewNames,
                            TRUE ~ rownames(dat))
  rownames(dat) <- dat$NewNames 
  dat <- dat %>% select(-c('NewNames')) 
  if (clr){ # clr normalize (optional)
    newmin <- 0.65*min(dat[dat > 0])
    dat[dat == 0] <- newmin; dat[is.na(dat)] <- newmin
    dat <- clr(t(dat))
  }
  else{
    dat <- t(dat)/rowSums(t(dat)) # Normalize by abundance
  }
  
  # Filter for taxa of interest; if not filtering, remove all taxa that are zero in > half of samples
  if (length(toi) == 0){
    dat <- dat[, names(colSums(dat > 0)[colSums(dat > 0) > length(rownames(dat))/2])] # Filter out taxa with low reads (nonzero in at least one half of all samples)
  }
  dat <- melt(dat) %>% as.data.frame() # Melt
  
  # Fix names
  colnames(dat) <- c('Sample_ID', 'Taxa', 'Value')
  df_taxa <- left_join(dat, metadata_in)
  
  # Filter
  if (length(toi) > 0){
    df_taxa <- df_taxa %>% filter(Taxa %in% toi)
  }
  
  minval2 = 0
  multiplier = 1
  if (mode == "logplot"){
    minval2 = min(df_taxa$Value[df_taxa$Value > 0]) # get second smallest value in dataframe for log plot
    multiplier = 100
  }
  p <- ggplot(df_taxa, aes(x = Treatment_Cycle, y = multiplier*(Value + minval2), fill = Treatment)) + 
    #geom_jitter(aes(color = Treatment), alpha = 0.2, width = 0.2) + 
    geom_line(aes(group=Patient_ID, color = Patient_ID), alpha = 0.2) +
    stat_summary(geom='ribbon', 
                 fun.data = mean_se, 
                 #fun.args=list(conf.int=0.95),
                 alpha = 0.5,
                 aes(fill = Treatment, group = Treatment)) + 
    facet_grid(Treatment ~ Taxa, scales ="free") + 
    theme_pubr() + 
    ylab('Relative abundance (percent)') + 
    xlab('Treatment cycle') + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values = treatment.colours) + 
    scale_color_manual(values = pt.colors) + 
    theme(legend.position = 'none')
  if (mode == "logplot"){ 
    p <- p + scale_y_log10() # log y axis
  }
  if (clr){
    p <- p + ylab('Abundance (CLR)')
  }
  fn <- paste0(savedir, taxa_level, "VsTime", suffix, "_Treatment.pdf")
  ggsave(plot = p, filename = fn, width = w, height = h, units = "in")
}

# Function to fit a linear model to taxa-level data
# Inputs: TaxaSummary dataframe, taxa_level (ex Phylum), corresponding first letter in lowercase (ex "p"),
# and metadata 
# savedir - directory to save resulting dataframes
# fitbatch - BOOLEAN whether to include the factor Batch_16S in the model
fit_lmer_model <- function(TaxaSummary, taxa_level = "Phylum", taxa_letter = "p", metadata_in, savedir, 
                           suffix = "All", Treatments = c("CAP", "TAS102", "CAP+IO"), fit_batch = TRUE){
  dat <- TaxaSummary[[taxa_level]]
  # Fix row names, including removing duplicate names
  pattern = paste0(".*", taxa_letter, "__")
  dat <- dat %>% as.data.frame
  dat$NewNames <- gsub(pattern, "", rownames(dat))
  dat$NewNames <- case_when(!duplicated(dat$NewNames) ~ dat$NewNames,
                            TRUE ~ rownames(dat))
  rownames(dat) <- dat$NewNames 
  dat <- dat %>% select(-c('NewNames')) 
  dat <- t(dat) # Transpose
  lowtax <- colnames(dat[, names(colSums(dat > 0)[colSums(dat > 0) < length(rownames(dat))/2])]) # save taxa that are present in <1/2 samples
  vlowtax <- colnames(dat[, names(colSums(dat > 0)[colSums(dat > 0) < length(rownames(dat))/10])]) # save taxa that are present in <1/10 samples
  # replace zeros with 0.65 the minimum value (i.e. 0.65 the observed detection limit), then transform
  newmin <- 0.65*min(dat[dat > 0 & is.finite(dat)])
  dat[dat == 0] <- newmin; dat[is.na(dat)] <- newmin
  dat <- clr(dat)
  dat <- melt(dat) %>% as.data.frame() # Melt
  colnames(dat) <- c('Sample_ID', 'Taxa', 'Abundance_CLR')
  df_taxa <- left_join(dat, metadata_in)
  df_taxa <- df_taxa %>% filter(Treatment %in% c(Treatments))
  # df_taxa <- df_taxa %>% filter(Taxa != vlowtax) # exclude extremely rare taxa
  
  # Create a dataframes for model results
  df_cat <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(df_cat) <- c('Taxa', 'Estimate', 'StdError', 'tValue')
  
  df_prepost <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(df_prepost) <- c('Taxa', 'Estimate', 'StdError', 'tValue')
  
  # For each taxa member, add it to the dataframe and fit a model
  all_taxa <- unique(df_taxa$Taxa)
  ntax <- length(all_taxa)
  nsamp <- length(df_taxa$Sample_ID)
  for (i in 1:ntax){
    tax <- all_taxa[i]
    tax_df <- df_taxa %>% filter(Taxa == tax)
    tax_df$Abundance_CLR[tax_df$Abundance_CLR == -Inf] <- min(tax_df$Abundance_CLR[tax_df$Abundance_CLR != -Inf])
    
    # Fit treatment_cycle_cat
    if (fit_batch){
      mod <- lmer(Abundance_CLR ~ Treatment_Cycle_Cat + (1 | Batch_16S) + (1|Patient_ID), data = tax_df)
      
    } else {
      mod <- lmer(Abundance_CLR ~ Treatment_Cycle_Cat + (1|Patient_ID), data = tax_df)
    }
    coef <- summary(mod)$coefficients %>% as.data.frame() 
    coef$Treatment_Cycle_Cat <- rownames(coef)
    rownames(coef) <- NULL
    coef$Taxa <- rep(tax, length(coef$Estimate))
    df_cat <- rbind(df_cat, coef)
    
    # fit treatment_cycle_prepost
    if (fit_batch){
      mod <- lmer(Abundance_CLR ~ Treatment_Cycle_PrePost + (1 | Batch_16S) + (1|Patient_ID), data = tax_df)
    } else {
      mod <- lmer(Abundance_CLR ~ Treatment_Cycle_PrePost + (1|Patient_ID), data = tax_df)
    }
    coef <- summary(mod)$coefficients %>% as.data.frame() 
    coef$Treatment_Cycle_PrePost <- rownames(coef)
    rownames(coef) <- NULL
    coef$Taxa <- rep(tax, length(coef$Estimate))
    df_prepost <- rbind(df_prepost, coef)
  }
  
  # Get p values and adjust for multiple hypothesis testing
  df_cat_final <- df_cat %>% filter(Treatment_Cycle_Cat != "(Intercept)")
  df_prepost_final <- df_prepost %>% filter(Treatment_Cycle_PrePost != "(Intercept)")
  
  df_cat_final$pValue <- 2*pt(q = abs(df_cat_final$`t value`), df = nsamp - 1, lower=FALSE)
  df_prepost_final$pValue <- 2*pt(q = abs(df_prepost_final$`t value`), df = nsamp - 1, lower=FALSE)
  
  df_cat_final$FDR <- p.adjust(df_cat_final$pValue, method = "BH")
  df_prepost_final$FDR <- p.adjust(df_prepost_final$pValue, method = "BH")
  
  # Sort and save dataframes
  df_cat_final <- arrange(df_cat_final, pValue)
  df_prepost_final <- arrange(df_prepost_final, pValue)
  
  fn <- paste0(savedir, taxa_level, "_PreVsEarlyVsLate_", suffix, ".csv")
  write.csv(df_cat_final, file = fn)
  
  fn <- paste0(savedir, taxa_level, "_PreVsPost_", suffix, ".csv")
  write.csv(df_prepost_final, file = fn)
}

## Function to fit all taxa-level models, for various taxa levels
fit_all_models <- function(TaxaSummary, metadata_in, savedir, suffix, Treatments = c("CAP", "CAP+IO", "TAS102"), fit_batch = TRUE){
  fit_lmer_model(TaxaSummary = TaxaSummary, taxa_level = "Phylum", taxa_letter = "p", metadata_in = metadata_in, savedir = savedir, suffix = suffix, Treatments = Treatments, fit_batch = fit_batch)
  fit_lmer_model(TaxaSummary = TaxaSummary, taxa_level = "Class", taxa_letter = "c", metadata_in = metadata_in, savedir = savedir, suffix = suffix, Treatments = Treatments, fit_batch = fit_batch)
  fit_lmer_model(TaxaSummary = TaxaSummary, taxa_level = "Order", taxa_letter = "o", metadata_in = metadata_in, savedir = savedir, suffix = suffix, Treatments = Treatments, fit_batch = fit_batch)
  fit_lmer_model(TaxaSummary = TaxaSummary, taxa_level = "Family", taxa_letter = "f", metadata_in = metadata_in, savedir = savedir, suffix = suffix, Treatments = Treatments, fit_batch = fit_batch)
  fit_lmer_model(TaxaSummary = TaxaSummary, taxa_level = "Genus", taxa_letter = "g", metadata_in = metadata_in, savedir = savedir, suffix = suffix, Treatments = Treatments, fit_batch = fit_batch)
  fit_lmer_model(TaxaSummary = TaxaSummary, taxa_level = "Species", taxa_letter = "s", metadata_in = metadata_in, savedir = savedir, suffix = suffix, Treatments = Treatments, fit_batch = fit_batch)
}

## Function to plot volcanoplot. Inputs: 
# tabledir - directory table is in
# level_in - either "Phylum", "Class", "Order", "Family", "Genus", "Species". Default Genus
# analysis_mode - either "_PreVsPost_" or "_PreVsEarlyVsLate_". Default "_PreVsPost_"
# subset_in - either "All", "CAP", "CAP+IO", "TAS102". Default "CAP"
# savedir - directory to save resulting plot in 
# fdrcut - FDR cutoff. Default 0.25
# estcut - estimate cutoff. Default 0.25. 
# hsize = height in inches. Default 6
# wsize = width in inches. Default 6
plot_volcano <- function(tabledir, savedir, level_in = 'Genus', analysis_mode = "_PreVsPost_", subset_in = "CAP",
                         fdrcut = 0.2, estcut = 0.3, hsize = 6, wsize = 6, label_plot = TRUE, suffix = ""){
  dat <- read.csv(paste0(tabledir, level_in, analysis_mode, subset_in, ".csv")) %>% as.data.frame() # dataframe with at least the following columns: FDR, Estimate, Taxa
  
  volc <- dat %>% 
    filter(!is.na(FDR)) %>% 
    mutate(color = case_when(
      Estimate > estcut & FDR < fdrcut ~ "Increased", 
      Estimate < -estcut & FDR < fdrcut ~ "Decreased", 
      FDR > fdrcut ~ "Nonsignificant"
    ))
  
  p <- ggplot(volc, aes(x = Estimate, y = -log10(FDR))) + 
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = -log10(fdrcut), color="grey", linetype=2) +
    geom_vline(xintercept = -estcut, color="grey", linetype=2) +
    geom_vline(xintercept = estcut, color="grey", linetype=2) +
    xlab("Estimate") +
    ylab(expression(-log[10]("FDR"))) +
    #gghighlight(color!="nonsignificant", label_key=Genus)
    geom_point(data=subset(volc, color=="Increased"), color="#669DB3FF") +
    geom_point(data=subset(volc, color=="Decreased"), color="#FF4F58FF") + 
    theme_pubr()
  
  if (label_plot == TRUE){
    p <- p + geom_text_repel(data=subset(volc, color %in% c("Decreased", "Increased")), aes(label=Taxa))
  }
  
  fn <- paste0(savedir, "Volcano", level_in, analysis_mode, subset_in, suffix, ".pdf")
  ggsave(plot = p, filename = fn, useDingbats=F, height=hsize, width=wsize)
}


## For ggplot color palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



# Function to normalize a dataframe by its columns (i.e. all columns sum to 1)
colNorm <- function(m){
  sweep(m,2,colSums(m),`/`)
}

# Filter taxa summary to only include taxa present in at least "samplecut" samples, with at least "countcut" counts/sample
filterTaxaSummary <- function(TaxaSummary, countcut = 10, samplecut = 3){
  TaxNew <- TaxaSummary
  sublevels = c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  for (sublevel in sublevels){
    taxsub <- TaxNew[[sublevel]]
    taxsub$Sum <- rowSums(taxsub >= countcut)
    taxsub <- taxsub %>% filter(Sum >= samplecut) %>%
      select(-Sum)
    TaxNew[[sublevel]] <- taxsub
  }
  return(TaxNew)
}

# Filter ASVs
filterSVTableData <- function(SVTableData, countcut = 10, samplecut = 3){
  taxsub <- SVTableData %>% as.data.frame()
  taxsub$Sum <- rowSums(taxsub >= countcut)
  taxsub <- taxsub %>% filter(Sum >= samplecut) %>%
    select(-Sum)
  return(taxsub)
}



## HUMAN
plotTaxaTimecourse_Human <- function(txSubset = c("CAP"), taxLevel = "Genus", tableDir = paste0(tabledir, "TaxaTables/"), 
                                     taxName = "d__Bacteria;  p__Firmicutes;  c__Clostridia;  o__Oscillospirales;  f__Ruminococcaceae;  g__Subdoligranulum", 
                                     separateTreatment = TRUE, rmZeros = TRUE, fillColor = 'purple',
                                     treatment.colours = NULL, patient.colours = NULL, spaghetti = TRUE, ncols = 5,
                                     TreatmentWrap = FALSE, returnDF = FALSE){
  if (length(taxName) == 0){
    return(NULL) # break out of loop if no taxa name is available
  }
  taxNameNew <- gsub(";  ", "...", taxName)
  taxNameNew <- gsub("\\|", "...", taxNameNew)
  taxNameNew <- gsub("; ", "..", taxNameNew)
  taxNameNew <- gsub("-", ".", taxNameNew)
  taxNameNew <- gsub("\\[", ".", taxNameNew)
  taxNameNew <- gsub("\\]", ".", taxNameNew)
  
  suffix = paste0(txSubset, collapse = "")
  fn <- paste0(tableDir, "TaxaSummary", taxLevel, "CLR_", suffix, ".csv")
  df <- read.csv(fn) # read in
  if (returnDF){
    return(df)
  }
  
  df <- df %>% filter(Treatment %in% txSubset) %>% select(all_of(taxNameNew), Treatment_Cycle, Treatment, Patient_ID)
  df <- melt(df, id.vars=c('Treatment', 'Treatment_Cycle', 'Patient_ID'), var='Taxa')
  colnames(df)[colnames(df) == 'value'] <- 'Relative_Abundance' 
  if (rmZeros){
    df_new = data.frame(matrix(ncol = length(colnames(df)), nrow = 0))
    colnames(df_new) <- colnames(df)
    for (name in taxNameNew){
      df_sub <- df %>% filter(Taxa == name)
      minval = min(as.numeric(df_sub$Taxa), na.rm = TRUE); epsilon = 0.01
      pids <- df_sub %>% group_by(Patient_ID) %>% dplyr::summarise(M = mean(Relative_Abundance, na.rm = TRUE)) %>% as.data.frame() %>% filter(M > minval + epsilon) %>% select(Patient_ID) %>% pull()
      df_sub <- df_sub %>% filter(Patient_ID %in% pids) 
      df_new = rbind(df_new, df_sub)
    }
    df = df_new
  }
  
  df$Patient_ID <- factor(df$Patient_ID) # set Patient_ID to factor for coloring
  df$Taxa <- gsub(".*__", "", df$Taxa) # rename taxa for simplicity in labeling and viewing
  p <- ggplot(df, aes(x = Treatment_Cycle, y = Relative_Abundance)) + 
    theme_pubr() +
    ylab('Relative abundance (CLR)') + 
    xlab('Treatment cycle') + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  if (spaghetti){
    p <- p + geom_line(aes(group=Patient_ID, color=Patient_ID), alpha = 0.2)  
    # + scale_color_manual(values = patient.colours)
  } else if (separateTreatment){
    p <- p + stat_summary(geom='line', fun = mean, alpha = 0.5, aes(color = Treatment, group = Treatment)) 
  } else {
    p <- p + stat_summary(geom='line', fun = mean, alpha = 0.5, aes(group = 'all'), color = fillColor)  
  }
  
  if (separateTreatment){
    p <- p + stat_summary(geom='ribbon', fun.data = mean_se, alpha = 0.5, aes(fill = Treatment, group = Treatment)) +
      scale_fill_manual(values = treatment.colours) 
  } else {
    p <- p + stat_summary(geom='ribbon', fun.data = mean_se, alpha = 0.5, aes(group = 'all'), fill = fillColor)  
  }
  
  if (!TreatmentWrap){
    p <- p + facet_wrap(~Taxa, scales = 'free_y', ncol = ncols)
  } else {
    p <- p + facet_grid(Treatment ~ Taxa, scales = 'free')
  }
  
  p <- p + theme(legend.position = 'none')
  return(p)
}

## MOUSE
plotTaxaTimecourse_Mice <- function(txSubset = c("Cap"), taxLevel = "Genus", tableDir = paste0(tabledir, "TaxaTables/"), 
                                    taxName = "d__Bacteria;  p__Firmicutes;  c__Clostridia;  o__Oscillospirales;  f__Ruminococcaceae;  g__Subdoligranulum", 
                                    rmZeros = TRUE, fillColor = 'purple',
                                    treatment.colours = NULL, patient.colours = NULL, spaghetti = TRUE, ncols = 5,
                                    TreatmentWrap = FALSE){
  if (length(taxName) == 0){
    return(NULL) # break out of loop if no taxa name is available
  }
  taxNameNew <- gsub(";  ", "...", taxName)
  taxNameNew <- gsub("; ", "..", taxNameNew)
  taxNameNew <- gsub("-", ".", taxNameNew)
  taxNameNew <- gsub("\\[", ".", taxNameNew)
  taxNameNew <- gsub("\\]", ".", taxNameNew)
  taxNameNew <- gsub(".*\\__", "",taxNameNew)
  
  fn <- paste0(tableDir, "TaxaSummary", taxLevel, "CLR.csv")
  df <- read.csv(fn) # read in
  df <- df %>% filter(Group %in% txSubset) %>% 
    filter(TissueType2 == "Stool") %>%
    filter(Taxa %in% taxNameNew) %>% 
    select(Group, Taxa, TreatmentCycle, Day, MouseID, Abundance_CLR)
  colnames(df)[colnames(df) == 'Abundance_CLR'] <- 'Relative_Abundance' 
  if (rmZeros){
    df_new = data.frame(matrix(ncol = length(colnames(df)), nrow = 0))
    colnames(df_new) <- colnames(df)
    for (name in taxNameNew){
      df_sub <- df %>% filter(Taxa == name)
      minval = min(as.numeric(df_sub$Relative_Abundance), na.rm = TRUE); epsilon = 0.01
      pids <- df_sub %>% group_by(MouseID) %>% dplyr::summarise(M = mean(Relative_Abundance, na.rm = TRUE)) %>% as.data.frame() %>% filter(M > minval + epsilon) %>% select(MouseID) %>% pull()
      df_sub <- df_sub %>% filter(MouseID %in% pids) 
      df_new = rbind(df_new, df_sub)
    }
    df = df_new
  }

  
  df$MouseID <- factor(df$MouseID) # set Patient_ID to factor for coloring
  p <- ggplot(df, aes(x = Day, y = Relative_Abundance)) + 
    theme_pubr() +
    ylab('Relative abundance (CLR)') + 
    xlab('Treatment duration (Days)')
  
  if (spaghetti){
    p <- p + geom_line(aes(group=MouseID, color=Group), alpha = 0.2)
  } else {
    p <- p + stat_summary(geom='line', fun = mean, alpha = 0.5, aes(color = Group))
    #p <- p + stat_summary(geom='line', fun = mean, alpha = 0.5, aes(group = 'all'), color = fillColor)  
  }
  
  p <- p + stat_summary(geom='ribbon', fun.data = mean_se, alpha = 0.5, aes(fill = Group))
  #p <- p + stat_summary(geom='ribbon', fun.data = mean_se, alpha = 0.5, aes(group = 'all'), fill = fillColor)  
  
  if (length(treatment.colours) > 0){
    p <- p + scale_color_manual(values = treatment.colours) + scale_fill_manual(values = treatment.colours)
  }
  
  if (length(unique(df$Taxa)) > 1){
    p <- p + facet_wrap(~Taxa, scales = 'free_y', ncol = ncols)
  }
  
  #p <- p + theme(legend.position = 'none')
  return(p)
}

## SIC
plotTaxaTimecourse_SIC <- function(txSubset = c("CAP"), taxLevel = "Genus", tableDir = paste0(tabledir, "TaxaTables/"), 
                                   taxName = "d__Bacteria;  p__Firmicutes;  c__Clostridia;  o__Oscillospirales;  f__Ruminococcaceae;  g__Subdoligranulum", 
                                   rmZeros = TRUE, fillColor = 'purple',
                                   treatment.colours = NULL, patient.colours = NULL, spaghetti = TRUE, ncols = 5,
                                   TreatmentWrap = FALSE){
  if (length(taxName) == 0){
    return(NULL) # break out of loop if no taxa name is available
  }
  taxNameNew <- gsub(";  ", "...", taxName)
  taxNameNew <- gsub("; ", "..", taxNameNew)
  taxNameNew <- gsub("-", ".", taxNameNew)
  taxNameNew <- gsub("\\[", ".", taxNameNew)
  taxNameNew <- gsub("\\]", ".", taxNameNew)
  taxNameNew <- gsub(".*\\__", "",taxNameNew)
  
  
  fn <- paste0(tableDir, "TaxaSummary", taxLevel, "CLR.csv")
  df <- read.csv(fn) # read in
  df <- df %>% filter(Treatment %in% txSubset) %>% 
    filter(Taxa %in% taxNameNew) %>% 
    select(Treatment, Taxa, Passage, Patient_ID, Abundance_CLR)
  colnames(df)[colnames(df) == 'Abundance_CLR'] <- 'Relative_Abundance' 
  if (rmZeros){
    df_new = data.frame(matrix(ncol = length(colnames(df)), nrow = 0))
    colnames(df_new) <- colnames(df)
    for (name in taxNameNew){
      df_sub <- df %>% filter(Taxa == name)
      minval = min(as.numeric(df_sub$Relative_Abundance), na.rm = TRUE); epsilon = 0.01
      pids <- df_sub %>% group_by(Patient_ID) %>% dplyr::summarise(M = mean(Relative_Abundance, na.rm = TRUE)) %>% as.data.frame() %>% filter(M > minval + epsilon) %>% select(Patient_ID) %>% pull()
      df_sub <- df_sub %>% filter(Patient_ID %in% pids) 
      df_new = rbind(df_new, df_sub)
    }
    df = df_new
  }
  
  df$Patient_ID <- factor(df$Patient_ID) # set Patient_ID to factor for coloring
  p <- ggplot(df, aes(x = Passage, y = Relative_Abundance)) + 
    theme_pubr() +
    ylab('Relative abundance (CLR)') + 
    xlab('Passage (#)')
  
  if (spaghetti){
    p <- p + geom_line(aes(group=Patient_ID, color=Patient_ID), alpha = 0.2)
  } else {
    p <- p + stat_summary(geom='line', fun = mean, alpha = 0.5, aes(color = Treatment))
  }
  
  p <- p + stat_summary(geom='ribbon', fun.data = mean_se, alpha = 0.5, aes(fill = Treatment))
  p <- p + facet_wrap(~Taxa, scales = 'free_y', ncol = ncols)
  colors <- c("grey", "orange", "red")
  p <- p + scale_colour_manual(values = colors)
  p <- p + scale_fill_manual(values = colors)
  #p <- p + theme(legend.position = 'none')
  p <- p + xlim(c(1,3))
  
  # CREATE NEW PLOT
  df$Passage <- factor(df$Passage)
  df <- df %>% filter(Passage %in% c("1", "3", "4", "6"))
  df$Treatment <- factor(df$Treatment, levels = c("Veh", "5FU", "CAP"), ordered = TRUE)
  p <- ggplot(df, aes(x = Passage, y = Relative_Abundance, color = Treatment)) + 
    theme_pubr() +
    geom_boxplot(outlier.shape = NA) + 
    ylab('Relative abundance (CLR)') + 
    xlab('Passage (#)')
  p <- p + scale_colour_manual(values = colors)
  p <- p + stat_compare_means(method = "kruskal.test")
  # p <- p + scale_fill_manual(values = colors)
  
  
  return(p)
}
