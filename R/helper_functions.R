#Helper functions
options(stringsAsFactors = FALSE)

##################################
### QC and normalisation steps ###
##################################

QC_filter <- function(
  scData.pass,
  exprs_data = "counts", 
  detect_limit = 4, 
  min_cells = 2, 
  min_batches = 1,
  plot_title = "Remove low-quality cells and features"
){
  
  # Calculate which features are expressed using the supplied thresholds (returns a logical matrix)
  is_exprs(scData.pass) <- calcIsExprs(scData.pass, lowerDetectionLimit = detect_limit, exprs_data = exprs_data)
  
  # Identify features that are not expressed in the required number or batches or cells
  batch <- pData(scData.pass)$batchvar
  qc.feat0 <- sapply(levels(batch), function(b) rowSums(is_exprs(scData.pass)[, batch==b]))
  qc.feat <- data.frame(
    fail_n_batches = apply(qc.feat0, 1, function(x) sum(x > 0) < min_batches),
    fail_n_cells = rowSums(qc.feat0) < min_cells
  )
  
  # Identify features that pass both cell and batch number tests
  qc.feat$pass <- !(qc.feat$fail_n_batches | qc.feat$fail_n_cells)
  
  # Plot the results
  vennDiagram(qc.feat, names=c(paste("<", min_batches,"batches"), paste("<", min_cells,"cells"), "pass"), 
              mar=c(6,6,4,4), cex = c(1,1,0.7), circle.col = c("firebrick","firebrick","forestgreen"))
  title(main = plot_title)
  legend("bottomleft", paste0("features [", length(qc.feat$pass), " to ", sum(qc.feat$pass), "]"), bty = "n", xjust = 0.5)
  
  # Subset and return the SCESet
  scData.pass <- scData.pass[qc.feat$pass,]
  scData.pass
}

calc_cell_RLE <- function(expr_mat, spikes = NULL) {
  # Local function to find per-gene RLE values. Return NA for genes with zero median expression.
  RLE_gene <- function(x) {
    if (median(unlist(x)) > 0) {
      log( (x + 1) / ( median(unlist(x)) + 1 ) ) / log(2)
    } else {
      rep(NA, times = length(x))
    }
  }
  
  # Calculate per-gene RLE values
  # Optionally remove spike-ins
  if(!is.null(spikes)) {
    RLE_matrix <- t(apply(expr_mat[-spikes, ], 1, RLE_gene))
  } else {
    RLE_matrix <- t(apply(expr_mat, 1, RLE_gene))
  }
  
  # Find median gene RLE value for each cell and return
  cell_RLE <- apply(RLE_matrix, 2, median, na.rm = T)
  cell_RLE
}

norm_plotRLE <- function(
  scData, 
  assays = c("exprs", "norm_exprs", "ruv_exprs"), 
  by_batch = TRUE,
  use_spikes = FALSE
) {
  
  # Option to use spike ins or endogenous genes
  f_set <- isSpike(scData)
  if(!use_spikes) f_set = !f_set
  
  # Option to facet by batch
  if(by_batch) {
    batch <- pData(scData)$batchvar
  } else {
    batch <- factor(rep("all cells", ncol(scData)))
  }
  
  rle <- list()
  
  for(assay in assays) {
    
    # Convert exprs values back to integer counts
    counts_mat <- round((2^assayData(scData)[[assay]])-1)
    
    # Calculate cell RLE values
    rle[[assay]] <- data.frame(
      exprs_type = rep(assay, ncol(counts_mat)),
      batch = batch,
      cell_rle = calc_cell_RLE(counts_mat[f_set,])
    )
  }
  
  # Reformat data for ggplot
  df_rle <- rbind.fill(rle)
  df_rle$exprs_type <- factor(df_rle$exprs_type, levels = assays)
  
  # Return plot
  ggplot(data = df_rle, 
         aes(x = exprs_type, y = cell_rle)) +
    geom_boxplot(aes(fill = exprs_type), alpha = 0.5) +
    facet_wrap( ~ batch) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_discrete() +
    xlab("normalisation method") +
    ylab("median cell RLE value") +
    theme_bw() +
    theme(legend.position = "none")
}

norm_var_features <- function(
  scData,
  assayName,
  span_endog = 0.1,
  span_spike = 0.3,
  family = "symmetric",
  use_endog = TRUE,
  plot_title = NULL
){
  
  # Scran method - total variance
  fit_endog <- trendVar(scData, assay = assayName, use.spikes = FALSE, span = span_endog, family = family)
  decomp_endog <- decomposeVar(scData, fit_endog, assay = assayName, get.spikes = FALSE)
  
  # Define biologically variable genes
  bv_genes <- (decomp_endog$bio > 0)
  
  # Draw plot
  plot(decomp_endog$mean[bv_genes], decomp_endog$total[bv_genes], xlab="Mean log-expression", ylab="Variance", pch = ".", col = "forestgreen")
  points(decomp_endog$mean[!bv_genes], decomp_endog$total[!bv_genes], pch = ".", col = "firebrick")
  o <- order(decomp_endog$mean)
  lines(decomp_endog$mean[o], decomp_endog$tech[o], col="blue", lwd=2)
  legend("topright", c("bvg","loess fit","non-bvg"), lty=c(NA,1,NA), lwd=2.5, pch=c(16,NA,16), col=c("forestgreen","blue","firebrick"), bty = "n", inset = 0.05)
  title(plot_title)
  
  # Exclude spike ins
  bv_genes[isSpike(scData)] <- FALSE
  
  # Return logical vector
  bv_genes
}

norm_RUVg <- function(scData,
                      ctrl_f = fData(scData)[["is_feature_control_ERCC"]],
                      assay_in = "norm_exprs",
                      assay_out = "ruv_exprs",
                      k = 1,
                      isLog = TRUE
) {
  
  #Batch effect adjustment supplied control feature set
  ruv_exprs <- RUVSeq::RUVg(
    assayData(scData)[[assay_in]],
    ctrl_f,
    k = k, isLog = isLog
  )$normalizedCounts
  
  #Correct negative expression values back to zero
  ruv_exprs[ruv_exprs < 0] <- 0
  
  #Add to SCESet and return it
  assayDataElement(scData, assay_out) <- ruv_exprs
  
  scData
}

#################################
### Analysis helper functions ###
#################################

.datTraits <- function() {
  traitData <- pData(scData)[incSamples,names(traits)]
  
  # convert data types, split factors
  allTraits <- sapply(1:length(traits), function(i) {
    if(is.numeric(traitData[,i]) | is.integer(traitData[,i])) {
      val_out <- traitData[,i]
    } else if(is.logical(traitData[,i]) | (is.factor(traitData[,i]) & !traits[[i]])) {
      val_out <- as.integer(traitData[,i])
    } else {
      traitData[,i] <- factor(traitData[,i], levels = unique(traitData[,i]))
      val_out <- sapply(levels(traitData[,i]), function(lvl) {
        as.integer(traitData[,i] == lvl)
      })
    }
    val_out
  })
  
  names(allTraits) <- names(traits)
  traitnames <- unlist(sapply(1:length(allTraits), function(x) {
    paste0(names(allTraits[x]),
           if(is.matrix(allTraits[[x]])) { paste0("_",colnames(allTraits[[x]])) })
  }))
  allTraits <- do.call(cbind, allTraits)
  
  # Form a data frame analogous to expression data that will hold the clinical traits.
  row.names(allTraits) <- sampleNames(scData)[incSamples]
  colnames(allTraits) <- traitnames
  data.frame(allTraits)  
}

.orderCells <- function() {
  if(all(unique(dynamicCells) %in% clust_order)) {
    ord_cells <- dynamicCells
    for(i in 1:length(clust_order)) {
      ord_cells[dynamicCells == clust_order[i]] <- i
    }
  } else {
    mod_col_ord <- rev(substring(cl_hr$labels[cl_hr$order],3))
    ord_cells <- cell_cols[cl_hc$order]
    for(i in 1:ncol(iMEs)) {
      ord_cells[ord_cells == mod_col_ord[i]] <- i
    }
    ord_cells[ord_cells == "grey" | is.na(ord_cells)] <- 0
  }
  
  ord_cells <- ord_cells[order(ord_cells, decreasing = FALSE)]
  hcl_rotate <- rotate(cl_hc, names(ord_cells))
  hcl_rotate
}

switchIDs <- function(ids, id_df = fData(scData), col_in = "feature_symbol", col_out = "feature_id") {
  out_df <- id_df[which(id_df[[col_in]] %in% ids ),c(col_in,col_out)]
  row.names(out_df) <- out_df[[col_in]]
  out_df <- out_df[ids,]
  return(out_df[[col_out]])
}

RswitchIDs <- function(ids, id_df = fData(scData), col_in = "feature_id", col_out = "feature_symbol") {
  return(switchIDs(ids,id_df,col_in,col_out))
}

.assign_celltype <- function() {
  if(use_ct_names) {
    mark_mm <- sapply(names(ct_markers), function(x) {
      mm <- geneModuleMembership[switchIDs(x),]
      best_mm <- max(mm)
      names(best_mm) <- substring(colnames(geneModuleMembership)[which(mm %in% best_mm)],3)
      best_mm
    })
    ct_df <- data.frame(ct_name = unlist(ct_markers), best_mm = mark_mm, mm_col = gsub(".*[.]","",names(mark_mm)))
    ct_df <- ct_df[order(ct_df$best_mm, decreasing = TRUE),]
    ct_df <- ct_df[!duplicated(ct_df$mm_col),]
    cell_type <- c(ct_df$ct_name,"unclassified")
    names(cell_type) <- c(ct_df$mm_col, "grey")
  } else {
    cell_type <- unique(moduleColors)
    names(cell_type) <- unique(moduleColors)
  }
  
  cell_type
}

.getAUC <- function(gene, labels) {
  ranked <- rank(gene);
  ms <- aggregate(ranked~unlist(labels),FUN=mean); #Get average score for each cluster
  posgroup <- as.character(unlist(ms[which(ms[,2]==max(ms[,2])),1])); #Get cluster with highest average score
  if (length(posgroup) > 1) {return (c(-1,-1,-1))} # Return negatives if there is a tie for cluster with highest average score (by definition this is not cluster specific)
  
  # Create 1/0 vector of truths for predictions, cluster with highest average score vs everything else
  truth <- labels == posgroup
  
  #Make predictions & get auc using RCOR package.
  pred <- ROCR::prediction(ranked,as.numeric(truth))
  val <- unlist(ROCR::performance(pred,"auc")@y.values)
  pval <- wilcox.test(gene[truth],gene[!truth])$p.value
  if (!exists("pval")) {pval <- 1}
  
  return(c(val,posgroup,pval))
}

.getMarkers <- function(expr_mat, labels) {
  if (length(labels) != length(expr_mat[1,])) {
    stop("Length of labels does not match number of cells.")
  }
  aucs <- apply(expr_mat,1,.getAUC,labels=labels)
  auc_df <- data.frame(matrix(unlist(aucs), ncol=3, byrow=TRUE))
  rownames(auc_df) <- rownames(expr_mat)
  colnames(auc_df) <- c("AUC","Group", "pval")
  auc_df$Group <- as.character(auc_df$Group)
  if(sum(auc_df$Group == "-1") > 0) {
    auc_df$Group[auc_df$Group == "-1"] <- "Ambiguous";
  }
  auc_df[,1] <- as.numeric(as.character(auc_df[,1]))
  auc_df[,3] <- as.numeric(as.character(auc_df[,3]))
  auc_df <- auc_df[auc_df[,1] > 0,]
  auc_df <- auc_df[order(-auc_df$AUC),]
  return(auc_df);
}

not_grey <- function(x) {
  x[x!="grey"]
}

###################################
### Analysis plotting functions ###
###################################

plotMEmerge <- function() {
  # Calculate eigengenes
  MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
  MEs <- MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss <- 1-cor(MEs);
  # Cluster module eigengenes
  METree <- hclust(as.dist(MEDiss), method = "average");
  # Plot the result
  par(mar=c(8,4,4,2))
  plot(as.dendrogram(METree), main = "Clustering of module eigengenes",
       xlab = "", sub = "", cex = 0.7)
  
  # Plot the cut line into the dendrogram
  abline(h=merge_height, col = "red")
}

plotCellDendro <- function() {
  d_plot <- rcl_hc
  d_plot$labels <- rep(".", length(d_plot$labels))
  d_plot$labels[!duplicated(dynamicCells)] <- dynamicCells[!duplicated(dynamicCells)]
  d_plot_col <- d_plot %>% as.dendrogram
  
  d_plot_col %>%  
    set("branches_k_color", k = length(unique(dynamicCells))) %>%
    set("branches_lwd", 2) %>% plot
  
  title("Cell clustering dendrogram")
  colored_bars(cell_cols, dend = d_plot)
}

plotCellClusters <- function() {
  rcl_hc_col <- rcl_hc %>% 
    as.dendrogram %>%
    set("branches_lwd", 1)
  
  heatmap.2(1-as.matrix(cl_dist), Rowv=rev(rcl_hc_col), Colv=rcl_hc_col, labCol = NA, labRow = NA, 
            col = viridis(75,option = "plasma"), scale="none", density.info="none", 
            ColSideColors = cell_cols,
            RowSideColors = cell_cols,
            trace = "none", margins = c(5, 5), cexRow = 0.5, cexCol = 0.5,
            breaks = seq(0,1,length.out = 76),
            keysize = 1.2, key.xlab = NA, key.ylab = NA, 
            key.xtickfun = function() {return(list(at=c(0,0.5,1),labels=c(0,0.5,1)))},
            xlab = paste(nrow(datExpr),"single cells"),
            ylab = paste(nrow(datExpr),"single cells")
  )
}

plotMEbars <- function() {
  df_brk <- rep(0,nrow(iMEs))
  df_brk[row.names(iMEs) %in% rcl_hc$labels[rcl_hc$order][match(unique(dynamicCells), dynamicCells[rcl_hc$order])]] <- 1
  df_MEs <- iMEs %>% cbind(cell = row.names(iMEs)) %>% cbind(brk = df_brk) %>% melt(id.vars = c("brk","cell"))
  df_MEs$variable <- df_MEs$variable %>% substring(3)
  df_MEs$cell <- factor(df_MEs$cell, levels = rcl_hc$labels[rcl_hc$order])
  ggplot(data = df_MEs, aes(x = cell, y = value, fill = variable, facet_x = factor(1))) + 
    geom_bar(colour = "black",stat = "identity", position = "identity", size = 0.3) + 
    scale_fill_identity() + 
    facet_grid(variable ~ .) +
    geom_vline(data=df_MEs[df_MEs$brk == 1,],aes(xintercept = as.numeric(cell) - 0.5), linetype = "dashed") +
    ggtitle("Clustered module eigengene expression") +
    theme_bw() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
}

plotPCAclusters <- function(ndims = 2, ... ) {
  scPCA <- scData[incFeatures,incSamples]
  pData(scPCA)$cluster <- cell_cols
  df_out <- pcs$x
  rownames(df_out) <- sampleNames(scPCA)
  reducedDimension(scPCA) <- df_out
  
  plotReducedDim(scPCA, ncomponents = ndims, colour_by = "cluster", ... ) + scale_fill_identity()
}

plotTSNE <- function() {
  col_scale <- names(subset_names)
  names(col_scale) <- subset_names
  
  ggplot(data = data.frame(cell_tSNE$Y, col = subset_names[cell_cols], detected_genes = pData(scData)$detected_genes[incSamples]),
         aes(x = X1, y = X2, fill = col, size = detected_genes)) + 
    geom_point(colour = "black", pch = 21, alpha = 0.7) + 
    xlab("dimension 1") + ylab("dimension 2") +
    scale_fill_manual(name = "subset", values = col_scale, guide = guide_legend(override.aes = list(size = 6))) + 
    scale_size_continuous(name = "detected \ngenes", limits=c(0,NA)) +
    theme_classic() + theme(axis.ticks = element_blank(), axis.text = element_blank())
}

plotGeneViolins <- function(
  geneNames,
  cex_col = 8,
  cex_row = 10,
  cex_y = 12,
  cex_ylab = 10,
  angle_col = 45,
  angle_row = 45,
  inc_ylab = TRUE,
  inc_legend = TRUE,
  leg_below = TRUE
) {
  # Get data
  geneIDs <- switchIDs(geneNames)
  
  # Remove NAs
  geneNames <- geneNames[!is.na(geneIDs)]
  geneIDs <- geneIDs[!is.na(geneIDs)]
  
  plot_df <- data.frame(
    cell_type = factor(subset_names[cell_cols], levels = c(unlist(ct_markers),"unclassified")), 
    group_fill = cell_cols,
    exprs = t(assayData(scData)[[ROC_assay]])[incSamples,geneIDs], 
    stringsAsFactors = FALSE)
  colnames(plot_df)[3:ncol(plot_df)] <- geneNames
  
  # Fill colours by mean expression (scaled)
  fill_stat <- sapply(levels(as.factor(cell_cols)), function(x) colMeans(as.matrix(plot_df[cell_cols == x,3:ncol(plot_df)])))
  if(!is.matrix(fill_stat)) fill_stat <- t(as.matrix(fill_stat))
  fill_max <- max(t(plot_df[,3:ncol(plot_df)]))
  fill_stat <- (fill_stat/rowMaxs(fill_stat)) * fill_max
  row.names(fill_stat) <- geneNames
  
  plot_df <- reshape2::melt(plot_df, id.vars = c("cell_type","group_fill"))
  plot_df$fill <- sapply(1:nrow(plot_df), function(x) fill_stat[plot_df$variable[x],plot_df$group_fill[x]])
  
  vio_plot <- ggplot(data = plot_df, aes(x=factor(1),y=value)) + 
    geom_violin(aes(fill=fill), alpha = 0.7, show.legend = inc_legend) +
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                 geom = "crossbar", width = 0.5) +
    ylim(c(0,NA)) +
    ylab("log2 ( normalised cpm + 1 )") +
    facet_grid(variable ~ cell_type, scales = "free_y")
  
  if(inc_ylab) {
    vio_plot <- vio_plot + theme_bw() + 
      theme(axis.text.y=element_text(size = cex_ylab))
  } else {
    vio_plot <- vio_plot + theme_classic() + 
      theme(axis.text.y=element_blank(), axis.ticks.y = element_blank())
  }
  
  vio_plot <- vio_plot +
    theme(
      strip.text.x = element_text(size = cex_col, angle = angle_col),
      strip.text.y = element_text(size = cex_row, angle = angle_row),
      strip.background = element_blank(),
      axis.ticks.x=element_blank(), 
      axis.line.x=element_blank(), 
      axis.line.y=element_blank(),
      axis.title.x=element_blank(), 
      axis.title.y=element_text(size = cex_y, margin = margin(0,30,0,0)),
      axis.text.x=element_blank(),
      legend.position = ifelse(leg_below, "bottom", "right")
    )
  
  vio_plot + 
    scale_fill_viridis(
      "", begin = 0.2, end = 0.8, option = "inferno", 
      limits = c(0,fill_max), breaks = c(0,fill_max),
      labels = c(0,"Maximum")
    )
}

tidylabel <- function(lab) {
  sapply(lab, function(x) {
    x <- gsub('(.{1,35})(\\s|$)', '\\1\n', x)
    x <- substring(x,1,nchar(x) - 1)
    x
  })
}

###############################
### Export results to excel ###
###############################

saveGeneLists <- function(clusts = not_grey(modNames)) {
  wb<-createWorkbook(type="xlsx")
  TITLE_STYLE <- CellStyle(wb)+ Font(wb, isBold=TRUE)
  
  for(i in 1:length(clusts)) {
    df_out <- ROCmark[ROCmark$Group == clusts[i],]
    df_out <- df_out[order(df_out$AUC, decreasing = TRUE),]
    df_out$AUC <- signif(df_out$AUC, digits = 2)
    df_out$pval <- signif(df_out$pval, digits = 3)
    sheet <- createSheet(wb, sheetName = gsub("/","_",subset_names[clusts[i]]))
    addDataFrame(subset(df_out, select = -Group), sheet, colnamesStyle = TITLE_STYLE)
    autoSizeColumn(sheet, 1:(ncol(df_out)-1))
  }
  
  saveWorkbook(wb, paste0("output/",params$project,"_genelist.xlsx"))
}

saveEnrichLists <- function(en_obj, suffix) {
  wb<-createWorkbook(type="xlsx")
  TITLE_STYLE <- CellStyle(wb)+ Font(wb, isBold=TRUE)
  
  enRes <- en_obj@compareClusterResult
  
  for(clust in unique(enRes$Cluster)) {
    df_out <- enRes[enRes$Cluster == clust,]
    df_out$pvalue <- signif(df_out$pvalue, digits = 3)
    df_out$p.adjust <- signif(df_out$p.adjust, digits = 3)
    df_out$qvalue <- signif(df_out$qvalue, digits = 3)
    sheet <- createSheet(wb, sheetName = gsub("/","_",clust))
    addDataFrame(df_out[,c(1:8,10,9)], sheet, colnamesStyle = TITLE_STYLE)
    autoSizeColumn(sheet, 1:9)
  }
  
  saveWorkbook(wb, paste0("output/",params$project, "_", suffix, "_", en_obj@fun, ".xlsx"))
}