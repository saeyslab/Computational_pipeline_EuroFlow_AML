correct_time <- function(ff){
  startTime <- ff@description$`$BTIM`
  endTime <- ff@description$`$ETIM`
  
  startTime <- as.POSIXlt(startTime, tryFormats = "%H:%M:%OS")
  endTime <- as.POSIXlt(endTime, tryFormats = "%H:%M:%OS")
  
  diffTime <- endTime - startTime
  
  return(diffTime)
}

AddClusterlabels <- function(p, 
                             clustering = NULL, 
                             whichGeom = "shadow", 
                             color = "cluster", ...){
  
  df <- ggplot2::ggplot_build(p)$plot$data
  if (is.null(clustering)) clustering <- df[, 3]
  if (!is.factor(clustering)) clustering <- factor(clustering)
  clusters <- levels(clustering)
  clusterCenters <- lapply(clusters, function(cl){
    clusterData <- df[which(clustering == cl), 1:2, drop = FALSE]
    if (nrow(clusterData) == 1) return(clusterData[1, ])
    medians <- apply(clusterData, 2, mean)
    return(medians)
  })

  clusterCenters <- Reduce(rbind, clusterCenters)
  clusterCenters <- data.frame(x = clusterCenters[, 1],
                               y = clusterCenters[, 2],
                               labels = clusters)
  
  
  if (whichGeom == "shadow"){
    if (color == "cluster"){
      p <- p + shadowtext::geom_shadowtext(data = clusterCenters,
                                           aes(x=x, y=y, label = labels, color = clusters),
                                           ...)
    } else {
      p <- p + shadowtext::geom_shadowtext(data = clusterCenters,
                                           aes(x=x, y=y, label = labels),
                                           col = color, ...)
    }

  } else if (whichGeom == "label"){
    if (color == "cluster"){
      p <- p + ggplot2::geom_label(data = clusterCenters,
                                   aes(x=x, y=y, label = labels, color = clusters),
                                   ...)
    } else {
      p <- p + ggplot2::geom_label(data = clusterCenters,
                                   aes(x=x, y=y, label = labels),
                                   col = color,...)
    } 
  } else if (whichGeom == "text"){
    if (color == "cluster"){
      p <- p + ggplot2::geom_text(data = clusterCenters,
                                  aes(x=x, y=y, label = labels, color = clusters),
                                  ...)
    } else {
      p <- p + ggplot2::geom_text(data = clusterCenters,
                                  aes(x=x, y=y, label = labels),
                                  col = color, ...)
    }
  } else {
    stop("whichGeom should be \"shadow\", \"label\" or \"text\"")
  }
    

  return(p)
}

minMaxScaling <- function(x, new_min, new_max){
  x <- new_min + ((x - min(x)) * (new_max - new_min)) / (max(x) - min(x))
  return(x)
}
transformCoordinates <- function(coordinates, transform, angle = 90){
  
  if (transform == "rotate"){
    angle <- angle * (pi / 180)
    M <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), nrow = 2)
  } else if (transform == "mirror h"){
    M <- matrix(c(1,0, 0, -1), nrow = 2)
  } else if (transform == "mirror v"){
    M <- matrix(c(-1,0, 0, 1), nrow = 2)
  } else {
   stop("transform should be \"rotate\", \"mirror h\" or \"mirror v\"") 
  }
  
  transformed_point <- function(coords){
    return(as.vector(M %*% coords))
  }
  
  coordinates <- t(apply(coordinates, 1, transformed_point))
  
  return(coordinates)
}

plot_relative_expression <- function(trajectory, 
                                     expression_source,
                                     scale = dynutils::scale_quantile,
                                     markers_oi = NULL, 
                                     stacked = F,
                                     file_dir = NULL,
                                     file_name= "Relative_expression.png"){
  
  linearised <- dynplot:::linearise_cells(trajectory = trajectory, 
                                          equal_cell_width = TRUE, 
                                          margin = 0.02)
  milestones_network <- linearised$milestone_network 
  expression <- get_expression(trajectory, expression_source)
  
  if (dynutils::is_sparse(expression)) {
    expression <- as.matrix(expression)
  }
  if (is.function(scale)) {
    expression <- scale(expression)
  } else if (is.logical(scale) && scale) {
    expression <- dynutils::scale_quantile(expression)
  }
  #smooth per part
  
  expression <- as.data.frame(expression)
  expression <- tibble::rownames_to_column(expression, "cell_id")
  df <- left_join(expression, linearised$progressions[, c("cell_id", 
                                                          "cumpercentage")],
                  by = "cell_id")
  df <- df[order(df$cumpercentage), ]
  
  
  
  
  milestones_network_order <- as.vector(t(milestones_network[, c("cumstart", "cumend")]))
  milestones <- milestones_network_order[ave(milestones_network_order, 
                                             milestones_network_order, 
                                             FUN = length) == 1]
  milestones_sewn <- as.data.frame(matrix(milestones, ncol = 2, byrow = T))
  colnames(milestones_sewn) <- c("cumstart", "cumend")
  
  
  smoothed_expression <- data.frame()
  for (i in 1:nrow(milestones_sewn)){
    begin <- milestones_sewn$cumstart[i]
    end <- milestones_sewn$cumend[i]
    chunk_df <- df %>% filter(cumpercentage >= begin, cumpercentage <= end)
    
    smoothed_expression_chunk <- apply(chunk_df[, -which(names(chunk_df) %in% 
                                                           c("cell_id", 
                                                             "cumpercentage"))], 
                                       2, 
                                       smooth.spline)
    smoothed_expression_chunk <- as.data.frame(sapply(smoothed_expression_chunk, 
                                                      "[", 
                                                      2))
    smoothed_expression_chunk <- cbind(smoothed_expression_chunk,
                                       cumpercentage = chunk_df$cumpercentage)
    smoothed_expression <- rbind(smoothed_expression, smoothed_expression_chunk)
    
    if (i > 1){
      na_begin <- milestones_sewn$cumend[i - 1] + 1
      na_end <- milestones_sewn$cumstart[i] - 1
      missing_data <- data.frame(cumpercentage = na_begin:na_end)
      smoothed_expression <- full_join(smoothed_expression, 
                                       missing_data,
                                       by = "cumpercentage")
    }
  }
  
  
  colnames(smoothed_expression) <- gsub("(.*)\\.y", "\\1", 
                                        colnames(smoothed_expression))
  
  molten <- reshape2::melt(smoothed_expression, 
                           id.vars = "cumpercentage",
                           variable.name = "Markers",
                           value.name = "Relative expression")
  
  
  if (!is.null(markers_oi)){
    molten <- molten %>% filter(Markers %in% markers_oi)
  }
  
  
  if (stacked == F){
    p <- ggplot(molten) + 
      geom_line(aes(x = cumpercentage, 
                    y = `Relative expression`, 
                    col = Markers), size = 1.2) + 
      xlab("Pseudotime") +
      geom_segment(data = milestones_network, aes(x = cumstart, 
                                                  xend = cumend, 
                                                  y = 0, 
                                                  yend = 0)) + 
      geom_text(data = milestones_network, aes(x = 0, 
                                               y = -0.05, 
                                               label = from[1])) + 
      geom_text(data = milestones_network, aes(x = cumend, y = -0.05, label = to)) + 
      theme_classic() +
      theme(axis.line = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank())
    
  } else {
    p <- ggplot(molten) + 
      geom_line(aes(x = cumpercentage, y = `Relative expression`, col = Markers), 
                size = 1.2) + 
      xlab("Pseudotime") + 
      scale_y_continuous(limits = c(-0.5, 1)) +
      geom_segment(data = milestones_network, aes(x = cumstart, 
                                                  xend = cumend, 
                                                  y = 0, 
                                                  yend = 0)) + 
      geom_text(data = milestones_network, aes(x = 0, 
                                               y = -0.05, 
                                               label = milestones_network$from[1])) + 
      geom_text(data = milestones_network, aes(x = cumend, y = -0.05, label = to)) + 
      theme_classic() + 
      theme(axis.line = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(), 
            legend.position = "none",
            strip.background = element_blank()) + 
      facet_wrap(~Markers, ncol = 1, scales = "free", strip.position = "left") 
    
  }
  
  
  if (is.null(file_dir)){
    print(p)
  } else {
    if (stacked == F){
      ggsave(filename = file_name, 
             plot = p, 
             device = "png", 
             path = file_dir, 
             width = 30, 
             height =  20, 
             units = "cm")
    } else {
      ggsave(filename = file_name,
             plot = p,
             device = "png",
             path = file_dir,
             width = 15,
             height = 5 * length(unique(molten$Markers)),
             units = "cm",
             limitsize = F)
    }
  }
}



spadeDownSampling <- function (ff, cols = NULL, kernel_mult = 5, apprx_mult = 1.5, 
                               med_samples = 2000, comp = TRUE, exclude_pctile = 0.01, 
                               target_pctile = NULL, target_number = round(nrow(ff) * 0.1), 
                               seed = 1) 
{
  #----Calculate Density----
  if (nrow(ff) <= target_number){
    message("Not downsampled")
    return(list(ffSubs = ff))
  }
  in_fcs <- ff
  in_data <- exprs(in_fcs)
  params <- parameters(in_fcs)
  pd <- pData(params)
  
  if (is.null(cols)) {
    cols <- as.vector(pd$name)
  }
  idxs <- match(cols, pd$name)
  if (any(is.na(idxs))) {
    stop("Invalid column specifier")
  }
  message("Calculating density")
  density <- SPADE.density(in_data[, idxs], kernel_mult = kernel_mult, 
                           apprx_mult = apprx_mult, med_samples = med_samples)
  
  channel_number <- ncol(in_fcs) + 1
  channel_id <- paste("$P", channel_number, sep = "")
  channel_name <- "density"
  channel_range <- max(density) + 1
  plist <- matrix(c(channel_name, channel_name, channel_range, 
                    0, channel_range - 1))
  rownames(plist) <- c("name", "desc", "range", "minRange", 
                       "maxRange")
  colnames(plist) <- c(channel_id)
  pd <- rbind(pd, t(plist))
  pData(params) <- pd
  
  out_data <- cbind(in_data, density = density)
  out_frame <- flowFrame(out_data, params, description = description(in_fcs))
  keyval <- list()
  keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
  keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
  keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
  keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
  keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
  keyword(out_frame) <- keyval
  
  in_fcs <- out_frame
  
  #----DownSample----
  message("Downsampling")
  in_data <- exprs(in_fcs)
  cell_ids <- 1:nrow(in_fcs)
  params <- parameters(in_fcs)
  pd <- pData(params)
  d_idx <- match("density", pd$name)
  if (is.na(d_idx)) {
    stop("No density parameter in FCS file")
  }
  boundary <- quantile(in_data[, d_idx], c(exclude_pctile, 
                                           target_pctile), names = FALSE)
  out_data <- subset(in_data, in_data[, d_idx] > boundary[1])
  cellsKept <- cell_ids[which(in_data[, d_idx] > boundary[1])]
  
  cat("  Targeting ", target_number, " events.\n")
  density <- out_data[, d_idx]
  if (target_number < nrow(out_data)) {
    density_s <- sort(density)
    cdf <- rev(cumsum(1/rev(density_s)))
    boundary <- target_number/cdf[1]
    if (boundary > density_s[1]) {
      targets <- (target_number - 1:length(density_s))/cdf
      boundary <- targets[which.min(targets - density_s > 
                                      0)]
    }
    set.seed(seed)
    random_uni_distr <- runif(length(density))
    out_data <- subset(out_data, boundary/density > random_uni_distr)
    cellsKept <- cellsKept[which(boundary/density > random_uni_distr)]
  } else if (target_number > nrow(out_data)) {
    stop("More events requested than present in file")
  }
  cellsRemoved <- setdiff(1:nrow(in_data), cellsKept)
  out_frame <- flowFrame(out_data, params, description = description(in_fcs))
  return(list(ffSubs = out_frame, cellsRemoved = cellsRemoved))
}

check_deGate <- function(ff, maxX = NULL, maxY = NULL, markers = colnames(ff), 
                         custom = NULL, removeZeros = FALSE, sortMarkers = FALSE,
                         plotFile = "./density_deGate.png", ...){

  markers <- GetMarkers(ff, markers)
  channels <- GetChannels(ff, markers)
  if (sortMarkers) markers <- sortMarkers(markers)
  plotList <- list()
  resultList <- list()
  for(marker in markers){
    ch <- GetChannels(ff, marker)
    if (removeZeros){
      ff_tmp <- ff[-which(exprs(ff)[, ch] == 0), ]
    } else {
      ff_tmp <- ff
    }
    p <- ggplot(data.frame(counts = exprs(ff_tmp)[, ch])) +
      geom_density(aes(x = counts))
    if (!is.null(maxX)) p <- p + xlim(0, maxX)
    if (!is.null(maxY)) p <- p + ylim(0, maxY)
    if (marker %in% names(custom)){
      cut <- custom[marker]
      p <- p +
        geom_vline(xintercept = cut, color = "deepskyblue3")
      resultList[[ch]] <- cut
    } else {
      tryCatch({
        cut <- deGate(ff_tmp, ch, ...)
        resultList[[marker]] <- cut
        p <- p +
          geom_vline(xintercept = cut, color = "red")
      }, error = function(e){
        message(paste0(marker,": no deGate"))
        resultList[[marker]] <- e
      }, warning = function(w){
        message(paste0(marker,": no deGate"))
        resultList[[marker]] <- w})
    }    
    p <- p + 
      theme_minimal() +
      ggtitle(marker)
    plotList[[marker]] <- p
  }
  
  if (!is.null(plotFile)){
    p <- ggpubr::annotate_figure(ggpubr::ggarrange(plotlist = plotList), top = ggpubr::text_grob(ff@description$FIL))
    png(plotFile, width = round(sqrt(length(markers))) * 300, height = round(sqrt(length(markers))) * 300)
    print(p)
    dev.off()
  } 
  resultList[["plot"]] <- plotList
  return(resultList)
}

FlowSOMCDD <- function(ff, clustering = NULL, markers = NULL, xdim = NULL, ydim = NULL, nClus = NULL, nCells, seed = 1){
  cat("FlowSOM\n")
  if (nrow(ff) <= nCells){
    message("Not downsampled")
    return(list(ffSubSampled = ff))
  }
  if (is.null(clustering)){
    fsom <- suppressWarnings(FlowSOM(ff,
                                     scale = FALSE,
                                     colsToUse = markers,
                                     xdim = xdim, ydim = ydim,
                                     nClus = nClus,
                                     seed = seed))
    metaclusters <- GetMetaclusters(fsom)
    pre_nCells <- nrow(fsom$FlowSOM$data)
  } else {
    pre_nCells <- length(clustering)
    metaclusters <- clustering
    fsom <- NULL
  }
  cat("Determine cells to remove\n")
  nCellsToRemove <- pre_nCells - nCells
  countMetaclusters <- data.frame(mcl = metaclusters) %>% 
    count(mcl) %>% 
    arrange(desc(n)) %>% 
    as.data.frame()
  differences <- countMetaclusters[1:(nrow(countMetaclusters) - 1), 2] - 
    countMetaclusters[2:nrow(countMetaclusters), 2]
  if (differences[1] < nCellsToRemove){
    iClusterStop <- which(cumsum(differences * 1:(length(differences))) <= nCellsToRemove)
    iClusterStop <- iClusterStop[length(iClusterStop)] 
    nCellsPerCluster <- c()
    for (i in seq_len(iClusterStop)){
      nCellsPerCluster <- c(nCellsPerCluster, sum(differences[i:iClusterStop]))
    }
    indexToBeRemoved <- c()
    for (i in seq_len(length(nCellsPerCluster))){
      cl <- countMetaclusters[i, 1]
      mclData <- which(metaclusters == cl)
      set.seed(seed)
      indices <- sample(x = mclData, size = nCellsPerCluster[i])
      indexToBeRemoved <- c(indexToBeRemoved, indices)
    }
    nStillToRemove <- nCellsToRemove - sum(nCellsPerCluster)
    leftoverRows <- seq.int(1,nrow(ff))[-indexToBeRemoved]
    stillToRemove <- sample(leftoverRows, nStillToRemove)
    indexToBeRemoved <- c(indexToBeRemoved, stillToRemove)
  } else {
    set.seed(seed)
    indexToBeRemoved <- sample(which(metaclusters == countMetaclusters[1, 1]), nCellsToRemove)
  }
  cat("Downsampling\n")
  ffSub <- ff[-indexToBeRemoved]
  res <- list(ffSubSampled = ffSub, cellsRemoved = sort(indexToBeRemoved), fsom = fsom)
  return(res)
} 


splitPlotStars <- function(fsom, 
                           markersPerTree, 
                           colorPalette = FlowSOM_colors,
                           file = NULL,
                           ...){
  metaclustering <- fsom$metaclustering
  fsom <- UpdateFlowSOM(fsom)
  markers <- unlist(markersPerTree)
  l <- length(markers)
  colors <- colorPalette(l)
  plotTrees <- list()
  plotStarLegends <- list()
  for (subTree in seq_along(markersPerTree)){
    markersSubtree <-  markersPerTree[[subTree]]
    colorsSub <- colors[1:length(markersSubtree)]
    p <- PlotStars(fsom, 
                   markers = markersSubtree,
                   colorPalette = colorsSub,
                   list_insteadof_ggarrange = TRUE,
                   ...)
    plotTrees[[subTree]] <- p$tree
    plotStarLegends[[subTree]] <- p$starLegend
    colors <- colors[-c(1:length(markersSubtree))]
  }
  legends <- ggarrange(plotlist = plotStarLegends, nrow = 1)
  trees <- ggarrange(plotlist = plotTrees, nrow = 1)
  legendsAndTrees <- ggarrange(legends, trees, nrow = 2, heights = c(1, 5))
  res <- ggarrange(legendsAndTrees, p[[3]], widths = c(15, 1), nrow = 1)
  if (!is.null(file)){
    pdf(file, width = length(markersPerTree) * 10, height = 15)
    print(res)
    dev.off()
  }
  return(res)
}


aggregatemaker <- function (fileNames, cTotal, writeOutput = FALSE, outputFile = "aggregate.fcs", 
                            writeMeta = FALSE, keepOrder = FALSE, verbose = FALSE) 
{
  nFiles <- length(fileNames)
  cFile <- ceiling(cTotal/nFiles)
  flowFrame <- NULL
  for (i in seq_len(nFiles)) {
    if (verbose) {
      message("Reading ", fileNames[i])
    }
    f <- flowCore::read.FCS(fileNames[i])
    res <- FlowSOMCDD(f, colnames(ff@exprs), xdim = 10, ydim = 10, nClus = 20, nCells = cFile)
    c <- seq(from = 1, to = nrow(f))[-res$cellsRemoved]
    if (keepOrder) 
      c <- sort(c)
    if (writeMeta) {
      utils::write.table(c, paste(gsub("[^/]*$", "", outputFile), 
                                  gsub("\\.[^.]*$", "", gsub(".*/", "", fileNames[i])), 
                                  "_selected_", gsub("\\.[^.]*$", "", gsub(".*/", 
                                                                           "", outputFile)), ".txt", sep = ""))
    }
    m <- matrix(rep(i, min(nrow(f), cFile)))
    m2 <- m + stats::rnorm(length(m), 0, 0.1)
    m <- cbind(m, m2)
    colnames(m) <- c("File", "File_scattered")
    prev_agg <- length(grep("File[0-9]*$", colnames(f)))
    if (prev_agg > 0) {
      colnames(m) <- paste0(colnames(m), prev_agg + 1)
    }
    f <- flowCore::cbind2(f[c, ], m)
    if (is.null(flowFrame)) {
      flowFrame <- f
      flowFrame@description$`$FIL` <- gsub(".*/", "", 
                                           outputFile)
      flowFrame@description$FILENAME <- gsub(".*/", "", 
                                             outputFile)
    }
    else {
      flowCore::exprs(flowFrame) <- rbind(flowCore::exprs(flowFrame), 
                                          flowCore::exprs(f)[, flowCore::colnames(flowCore::exprs(flowFrame))])
    }
  }
  flowFrame@description[[paste("flowCore_$P", ncol(flowFrame) - 
                                 1, "Rmin", sep = "")]] <- 0
  flowFrame@description[[paste("flowCore_$P", ncol(flowFrame) - 
                                 1, "Rmax", sep = "")]] <- nFiles + 1
  flowFrame@description[[paste("flowCore_$P", ncol(flowFrame), 
                               "Rmin", sep = "")]] <- 0
  flowFrame@description[[paste("flowCore_$P", ncol(flowFrame), 
                               "Rmax", sep = "")]] <- nFiles + 1
  flowFrame@description[[paste("$P", ncol(flowFrame) - 1, 
                               "B", sep = "")]] <- 32
  flowFrame@description[[paste("$P", ncol(flowFrame), "B", 
                               sep = "")]] <- 32
  flowFrame@description$FIL <- gsub(".*/", "", outputFile)
  if (writeOutput) {
    flowCore::write.FCS(flowFrame, filename = outputFile)
  }
  flowFrame
}

parse_flowjo <- function(files,
                         wsp_file,
                         group = "All Samples",
                         plot = FALSE) {
  wsp <- flowWorkspace::openWorkspace(wsp_file)
  o <- capture.output(
    gates <- suppressMessages(flowWorkspace::parseWorkspace(wsp, group))
  )
  files_in_wsp <- gates@data@origSampleVector
  counts <- as.numeric(gsub(".*_([0-9]*)$", "\\1", files_in_wsp))
  result <- list()
  for(file in files){
    print(paste0("Processing ", file))
    file_id <- grep(gsub(".*/", "", file), files_in_wsp)
    if(length(file_id) == 0) {stop("File not found. Files available: ",
                                   gsub("_[0-9]*$", "\n", files_in_wsp))}
    gate_names <- flowWorkspace::getNodes(gates[[file_id]], path = "auto")
    gatingMatrix <- matrix(FALSE,
                           nrow = counts[file_id],
                           ncol = length(gate_names),
                           dimnames = list(NULL, gate_names))
    for (gate in gate_names) {
      gatingMatrix[, gate] <- flowWorkspace::getIndiceMat(gates[[file_id]],
                                                          gate)
    }
    ff <- flowWorkspace::getData(gates[[file_id]], "root")
    ff@exprs[, "Time"] <- ff@exprs[, "Time"] * 100
    result[[file]] <- list("flowFrame" = ff,
                           "gates" = gatingMatrix)
    
    if (plot) {
      flowWorkspace::plot(gates[[file_id]])
    }
  }
  if (length(files) == 1){
    result <- result[[1]]
  } else {
    result <- list(flowSet = flowCore::flowSet(lapply(result, function(x) x$flowFrame)),
                   gates = lapply(result, function(x) x$gates))
  }
  return(result)
}

manual_vector <- function(manual_matrix, cell_types){
  
  if(is.list(manual_matrix)){ manual_matrix <- do.call(rbind, manual_matrix) }
  
  manual <- rep("Unknown",nrow(manual_matrix))
  for(cellType in cell_types){
    manual[manual_matrix[,cellType]] <- cellType
  }
  manual <- factor(manual, levels=c("Unknown",cell_types))
  return(manual)
}

label_metaclusters <- function(fsom, manual_labels){
  counts <- as.matrix(table(FlowSOM::GetMetaclusters(fsom),
                            manual_labels))
  metacluster_names <- apply(counts,
                             1,
                             function(x) colnames(counts)[which.max(x)])
  metacluster_names <- number_duplicates(metacluster_names)
  return(metacluster_names)
}

label_clusters <- function(fsom, manual_labels){
  counts <- as.matrix(table(FlowSOM::GetClusters(fsom),
                            manual_labels))
  cluster_names <- apply(counts,
                         1,
                         function(x) colnames(counts)[which.max(x)])
  return(cluster_names)
}

number_duplicates <- function(x){
  counts <- table(x)
  for (value in names(counts)) {
    if (counts[value] > 1) {
      x[which(x == value)] <- paste0(value, "_", seq_len(counts[value]))
    }
  }
  return(x)
}


preprocess <- function(ff, 
                       transformList = NULL, 
                       QC_algorithm = "PeacoQC", 
                       QC_file = FALSE,
                       linear_rescale = FALSE,
                       plot_densities = FALSE,
                       ...)
{
  
  cat(paste("Processing file:", ff@description$FILENAME, "\n"))
  cat("-----------------------------------------\n")
  #----Remove Margins----
  indices_margins <- PeacoQC::RemoveMargins(ff, colnames(ff), output = "full")$indices_margins
  selection_inrange <- !seq_len(nrow(ff)) %in% indices_margins 
  selection <- selection_inrange
  perc_remove_margins <- (sum(!selection_inrange) / length(selection_inrange)) * 100
  cat(paste0("Remove margins removed ", round(perc_remove_margins, 2), "% of cells.\n\n"))
  ff_tmp1 <- ff[selection_inrange, ]
  #----Compensation----
  ff <- compensate(ff, ff@description$SPILL)
  
  #----Transformation----
  if (is.null(transformList)){
    transformList <- estimateLogicle(ff, channels = colnames(ff@description$SPILL))
  }
  ff <- transform(ff, transformList)

  #----FlowAI----
  if (QC_algorithm == "FlowAI"){
    ff_qc <- flow_auto_qc(FlowAI::ff,
                          output = 2,
                          html_report = FALSE,
                          mini_report = FALSE,
                          fcs_QC = FALSE,
                          folder_results = FALSE)
    selection_quality <- ff_qc@exprs[, "QCvector"] < 10000
    perc_selection_q <- sum(!selection_quality) / length(selection_quality) * 100
  }
  
  #----PeacoQC----
  if (QC_algorithm == "PeacoQC"){
    ff_qc <- suppressWarnings(PeacoQC::PeacoQC(ff = ff[selection_inrange, ],
                                               channels = colnames(ff)[-which(colnames(ff) == "Time")],
                                               plot = T,
                                               save_fcs = F,
                                               report = T,
                                               output_directory = ".",
                                               ...))
    selection_quality <- ff_qc$GoodCells
    perc_selection_q <- ff_qc$PercentageRemoved
  }
  cat(paste0("\n", QC_algorithm, " removed ", round(perc_selection_q, 2), "% of cells.\n"))
  selection[selection] <- selection_quality
  ff_tmp2 <- ff_tmp1[selection_quality, ]
  #----Debris----
  FSC_cuts <- deGate(ff[selection, ], "FSC-A", all.cuts = TRUE,
                     upper = FALSE)
  FSC_cut <- FSC_cuts[which.min(abs(FSC_cuts - 35000))]
  if(FSC_cut > 65000) FSC_cut <- 35000
  selection_debris <- no_debris(ff[selection, ], FSC_cut)
  perc_selection_debris <- (sum(!selection_debris) / length(selection)) * 100
  cat(paste0("File has ", round(perc_selection_debris, 2), "% of debris.\n"))
  selection[selection] <- selection_debris
  ff_tmp3 <- ff_tmp2[selection_debris, ]
  #----Doublets----
  selection_single <- is_single(ff[selection, ], nMAD = 5)
  perc_selection_doublets <- (sum(!selection_single) / length(selection)) * 100
  cat(paste0("File has ", round(perc_selection_doublets, 2), "% of doublets.\n"))
  tmp_selection <- selection
  selection[selection] <- selection_single
  #----Plot densities----
  if (plot_densities){
    cat("Plotting densities\n")
    png(paste0("./Densities/", ff@description$`$FIL`, ".png"), width = 1200, height = 300)
    par(mfrow = c(1, 4))
    plotDens(ff, c("FSC-A", "SSC-A"), main = "Margin events")
    points(ff@exprs[, c("FSC-A","SSC-A")][!selection_inrange,], col = "red")
    plotDens(ff_tmp1, c("FSC-A", "SSC-A"), main = "PeacoQC")
    points(ff_tmp1@exprs[, c("FSC-A","SSC-A")][!selection_quality,], col = "red")
    plotDens(ff_tmp2, c("FSC-A", "SSC-A"), main = "Debris")
    points(ff_tmp2@exprs[, c("FSC-A","SSC-A")][!selection_debris,], col = "red")
    plotDens(ff_tmp3, c("FSC-A", "FSC-H"), main = "Doublets")
    points(ff_tmp3@exprs[, c("FSC-A","FSC-H")][!selection_single,], col = "red")
    dev.off()
  }
  
  ff <- ff[selection, ]
  perc_total <- (sum(!selection) / length(selection)) * 100

  cat(paste0("\nIn total, ", round(perc_total, 2), "% was removed.\n"))
  cat("-----------------------------------------\n\n")
  
  #----Linear rescaling----
  if (linear_rescale == T){
    rescale_linear <- function(x){(x - min(x)) / (max(x) - min(x)) * 5} 
    #maybe here it's better to look for quantiles instead of min-max
    for (ch in c("FSC-A", "FSC-H", "SSC-A")) {
      ff@exprs[, ch] <- rescale_linear(ff@exprs[, ch])
    }
  }
  
  #----Print QC_file----
  if (QC_file == T){
    if (!file.exists("Quality_control.txt")){
      write("Quality control\n------------------------\n", "Quality_control.txt", 
            append = T)
      write(paste0("Filename\tMargins\t", QC_algorithm, "\tDebris\tDoublet\tTotal"), 
            file = "Quality_control.txt",
            append = T)
    }
    write(paste0(ff@description$FILENAME, "\t", 
                 round(perc_remove_margins, 2), "\t", 
                 round(perc_selection_q, 2), "\t",
                 round(perc_selection_debris, 2), "\t",
                 round(perc_selection_doublets, 2), "\t",
                 round(perc_total, 2)),
          file = "Quality_control.txt",
          append = T)
  }
  
  #----Return----
  return(list(FinalFF = ff, Selection = selection))
}

mean_jitter <- function(x){
  missing <- is.na(x)
  x[missing] <- mean(x, na.rm = T) + 
    rnorm(length(x[missing]), mean = 0, sd = 0.01)
  return(x)
}

multiUnivariateCox <- function(data, #all the data including covariates if not NULL
                               covariates = NULL, #colname of covariates
                               time, #Survival time
                               event){ # Survival status
  if (!time %in% colnames(data)) stop(paste0(time, " is not a column in data.")) 
  if (!event %in% colnames(data)) stop(paste0(event, " is not a column in data.")) 
  for (covariate in covariates){
    if (!covariate %in% colnames(data)) stop(paste0(covariate, 
                                                    " is not a column in data."))
  }
  data[, event] <- as.numeric(as.character(data[, event]))
  
  features <- colnames(data)[which(!colnames(data) %in% c(time, event, covariates))]
  
  features_form <- sapply(features, function(x){
    if (!is.null(covariates)){
      x <- paste0(x, " + ", paste0(covariates, collapse = " + "))
    }
    as.formula(paste("Surv(", time, ",", event, ") ~", x))
  })
  resList <- list()
  cox_models <- lapply(features_form, function(x) {survival::coxph(x, data)})
  resList[["individual_models"]] <- cox_models
  cox_results <- lapply(cox_models,
                        function(x){ 
                          x <- summary(x)
                          p_value <- signif(x$coefficients[1, 5], digits = 2)
                          wald_test <- signif(x$wald["test"], digits = 2)
                          beta <- signif(x$coef[1], digits = 2)
                          res <- c(beta,
                                   wald_test, 
                                   p_value)
                          names(res) <- c("Beta", 
                                          "Wald_test", 
                                          "p_value")
                          return(res)
                        })
  cox_results <- as.data.frame(do.call(rbind, cox_results))
  cox_results$p_BH_adj <- p.adjust(cox_results$p_value, method = "BH") 
  cox_results <- cox_results %>% dplyr::arrange(p_BH_adj)
  resList[["summary"]] <- cox_results
  return(resList)
}


rotate_fsom <- function(fsom, angle){
  rotation_matrix <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), nrow = 2)
  fsom$MST$l <- t(apply(fsom$MST$l, 1, function(x) rotation_matrix %*% x))  
  return(fsom)
}

plot_survivalplots_median <- function(data, 
                                      features_to_plot,
                                      ncols,
                                      plot_dir = NULL)
{
  medians <- apply(data[, features_to_plot], 2, function(x){
    median <- median(x, na.rm = T)
    x <- x > median
    x <- gsub(TRUE, "Above Median", x)
    x <- gsub(FALSE, "Under Median", x)
  })
  
  data[, features_to_plot] <- medians
  plot_list <- list()
  for (feature in features_to_plot){
    f <- as.formula(paste("Surv(Survival_time, Death_or_alive) ~", feature))
    fit <- surv_fit(f, data = data)
    
    p <- ggsurvplot(fit,
                    data = data, 
                    pval = T,
                    title = feature,
                    subtitle = " ")
    
    plot_list[[feature]] <- p
  }
  
  nrows_plot <- ceiling(length(features_to_plot ) / ncols)
  
  if (is.null(plot_dir)){
    arrange_ggsurvplots(plot_list, print = T, ncol = ncols, nrow = nrows_plot)
  } else {
    png(plot_dir, width = ncols * 400, height = nrows_plot * 350)
    arrange_ggsurvplots(plot_list, print = T, ncol = ncols, nrow = nrows_plot)
    dev.off()
  }
}


ff_density <- function(ff, maxX = NULL, markers = NULL, custom = NULL){
  if (!is.null(markers)){
    markers <- get_channels(ff, markers)
  } else {
    markers <- colnames(ff)
  }
  plotList <- list()
  markers <- sortMarkers(markers)
  i <- 1
  max <- length(markers)
  for(marker in markers){
    cat(i, "of ", max, "\r")
    i <- i + 1
    ch <- get_markers(ff, marker)
    p <- ggplot(data.frame(counts = ff@exprs[,marker])) +
      geom_density(aes(x = counts)) +
      geom_rug(aes(x = counts), col = "red", alpha = 0.3) 
    if (!is.null(maxX)){
      p <- p + xlim(0, maxX)
    }
    p <- p + 
      theme_minimal() +
      theme(legend.position = "none") +
      ggtitle(ch)
    plotList[[ch]] <- p
  }
  return(plotList)
}

sortMarkers <- function(markers){
  markers <- sort(markers)
  whichCD <- grep("^CD", markers)
  CD_markers <- markers[whichCD]
  CD_markers <- as.numeric(gsub("CD([0-9]*\\.?[0-9]*).*", "\\1", CD_markers))
  order_CD_markers <- order(CD_markers)
  markers[grep("^CD", markers)] <- markers[whichCD][order_CD_markers]
  return(markers)
}

jaccard.table <- function(clustering_1, 
                          clustering_2){
  
  if(length(clustering_1) != length(clustering_2)){
    stop("The two cluster results should have the same length")
  }
  if(!is.factor(clustering_1)) clustering_1 <- factor(clustering_1)
  if(!is.factor(clustering_2)) clustering_2 <- factor(clustering_2)
  
  jaccard <- matrix(NA,
                    nrow = length(levels(clustering_1)),
                    ncol = length(levels(clustering_2)),
                    dimnames = list(levels(clustering_1),
                                    levels(clustering_2)))
  for (cluster1 in levels(clustering_1)){
    for (cluster2 in levels(clustering_2)){
      intersection <- sum(clustering_1 == cluster1 & 
                            clustering_2 == cluster2)
      union <- sum(clustering_1 == cluster1 | 
                     clustering_2 == cluster2)
      jaccard[cluster1, cluster2] <-  intersection / union
    }
  }
  return(jaccard)
}  

expression_matrix_degate <- function(ff, cl, st){
  degates <- check_deGate(ff, removeZeros = T)
  expression <- lapply(names(degates), function(degate){
    if (is.infinite(degates[[degate]])) degates[[degate]] <- max(ff[, degate]@exprs)
    exprMarker <- ff[, degate]@exprs > degates[[degate]]
    return(exprMarker)
  })
  expression <- as.data.frame(Reduce(cbind, expression))
  resultMatrix <- matrix(0, ncol = ncol(expression), nrow = ncol(expression),
                         dimnames = list(colnames(expression), 
                                         colnames(expression)))
  for (i in seq_len(nrow(expression))){
    row <- expression[i, ]
    cellMarker <- colnames(expression)[which(row == TRUE)]
    combinations <- expand.grid(cellMarker, cellMarker)
    for (j in seq_len(nrow(combinations))){
      x <- as.character(combinations[j, 1])
      y <- as.character(combinations[j, 2])
      resultMatrix[x, y] <- resultMatrix[x, y]  + 1
    }
  }
  return(resultMatrix)
}

co_expression_matrix <- function(expression_matrix){
  nCols <- ncol(expression_matrix)
  resultMatrix <- matrix(0, ncol = nCols, nrow = nCols,
                         dimnames = list(colnames(expression_matrix),
                                         colnames(expression_matrix)))
  for (row in seq_len(nrow(expression_matrix))){
    co_expressed <- which(expression_matrix[row, ] == TRUE)
    resultMatrix[co_expressed, co_expressed] <- resultMatrix[co_expressed, co_expressed] + 1
  }
  return(resultMatrix)
}


merge.png.pdf <- function(plotList, pdfName, pdfHeight = 8, pdfWidth = 8, 
                          pngHeight = 600, pngWidth = 600, deletePngFiles = FALSE) {
  tmp_dir_name <- paste0("./tmp_", gsub(" |:|-","", Sys.time()))
  dir.create(tmp_dir_name)
  for (pic_i in seq_along(plotList)){
    png(paste0(tmp_dir_name, "/", pic_i, ".png"), width = pngWidth, height = pngHeight)
    print(plotList[[pic_i]])
    dev.off()
  }
  pngFiles <- paste0(tmp_dir_name, "/", list.files(path = tmp_dir_name, pattern = ".png$"))
  pdf(pdfName, height = pdfHeight, width = pdfWidth)
  n <- length(pngFiles)
  for( i in 1:n) {
    pngFile <- pngFiles[i]
    pngRaster <- png::readPNG(pngFile)
    grid::grid.raster(pngRaster, width=grid::unit(0.8, "npc"), height= grid::unit(0.8, "npc"))
    if (i < n) plot.new()
  }
  dev.off()
  if (deletePngFiles) {
    unlink(pngFiles)
    unlink(tmp_dir_name, recursive = T)
  }
}


hierarchicalFeatures <- function(fsom, percentages, thresholds = NULL,
                                 hierarchical_clustering = NULL, ...){
  #' @param fsom        A FlowSOM object
  #' @param percentages A dataframe with the cluster (columns) percentages per 
  #'                    sample (row)
  #' @param thresholds  If a percentage is lower than this threshold, the ratio
  #'                    will not be calculated            
  #' @param clustering  An optional hclust result                       
 
  if (is.null(hierarchical_clustering)){
    hierarchical_clustering <- hclust(dist(fsom$map$codes), ...)
    hierarchical_clustering$labels <- paste0("Cl", seq_len(NClusters(fsom)))
  }
  merge_groups <- list()
  
  # Combine groups hierarchically
  for (merge_id in 1:nrow(hierarchical_clustering$merge)){
    merge_groups[[merge_id]] <- list()
    groups <- hierarchical_clustering$merge[merge_id,]
    for(side in 1:2){
      if(groups[side] < 0){
        merge_groups[[merge_id]][[side]] <- -groups[side]
      } else {
        merge_groups[[merge_id]][[side]] <- c(merge_groups[[groups[side]]][[1]], 
                                              merge_groups[[groups[side]]][[2]])
      }
    }
  }
  
  # Prepare empty matrices
  ratios <- matrix(NA,
                   nrow = nrow(percentages),
                   ncol = length(merge_groups),
                   dimnames = list(rownames(percentages), 
                                   seq_len(length(merge_groups))))
  parent_clusters <- matrix(NA,
                            nrow = nrow(percentages),
                            ncol = length(merge_groups),
                            dimnames = list(rownames(percentages), 
                                            seq_len(length(merge_groups))))
  
  perc_of_parents <- matrix(NA,
                            nrow = nrow(percentages),
                            ncol = length(merge_groups),
                            dimnames = list(rownames(percentages), 
                                            seq_len(length(merge_groups))))
  
  for (merge_id in seq_along(merge_groups)){
    namesSide1 <- hierarchical_clustering$labels[merge_groups[[merge_id]][[1]]]
    namesSide2 <- hierarchical_clustering$labels[merge_groups[[merge_id]][[2]]]
    colnames(ratios)[merge_id] <- 
      paste("(", paste(namesSide1, collapse = " + "), ") / (",
            paste(namesSide2, collapse = " + "), ")")
    colnames(parent_clusters)[merge_id] <- 
      paste(c(namesSide1, namesSide2), collapse = " + ")
    colnames(perc_of_parents)[merge_id] <- 
      paste("(", paste(namesSide1, collapse = " + "), ") / (",
            paste(c(namesSide1, namesSide2), collapse = " + "), ")")
    
    # Calculate sums of the sides
    side1 <- rowSums(percentages[, merge_groups[[merge_id]][[1]], drop = FALSE]) + 1e-100
    side2 <- rowSums(percentages[, merge_groups[[merge_id]][[2]], drop = FALSE]) + 1e-100
  
    # Calculate sums
    parent_clusters[, merge_id] <- side1 + side2
    
    # Calculate percentage of parents
    perc_of_parents_tmp <- side1 / (side1 + side2)
    if (!is.null(thresholds)){
      belowThreshold <- (side2 + side1) < thresholds
      perc_of_parents_tmp[belowThreshold] <- NA
    }
    perc_of_parents[, merge_id] <- perc_of_parents_tmp
    
    # Calculate ratio's
    ratio_tmp <- side1 / side2
    if (!is.null(thresholds)){
      belowThreshold <- side1 < thresholds | side2 < thresholds
      ratio_tmp[belowThreshold] <- NA
    }
    ratios[, merge_id] <- ratio_tmp
  }
  return(list("ratios" = ratios, "parent_clusters" = parent_clusters,
              "perc_of_parents" = perc_of_parents, 
              "hierarchical_clustering" = hierarchical_clustering))
}

PlotClusterDendrogram <- function(hierarchical_clustering, #hclust
                                  values, #-log10p values or other values for the nodes
                                  values_title, #title of legend
                                  leaf_col, #colors for the leaves
                                  leaf_labels = NULL, #labels for leaves
                                  main){
  
  values_col <- rev(RColorBrewer::brewer.pal(8, "RdYlBu"))[
    cut(values, 
        breaks = 7,
        labels = FALSE)]
  
  
  get_h_properties <- function(x){
    if(!is.null(attr(x, "leaf"))) {
      return( list(height = attr(x, "height")))
    } else {
      return( list(height = attr(x, "height"),
                   child1 = get_h_properties(x[[1]]),
                   child2 = get_h_properties(x[[2]])))
    }
  }
  dend <- as.dendrogram(hierarchical_clustering)
  
  heights <- get_h_properties(dend)
  heights <- unlist(heights)
  heights_o <- order(heights)
  nodes_col <- rep("black", length(values)*2 + 1)
  nodes_col[heights_o[(length(values) + 2):(length(values) * 2 + 1)]] <- values_col[seq_len(length(values))]
  
  layout(matrix(c(1,rep(2, 3),
                  3, rep(2, 3)), nrow = 2, byrow = TRUE))
  plot(values,
       col = values_col,
       pch = 19,
       ylab = values_title,
       xlab = "",
       xaxt = "n")
  
  plot(dend %>% 
         set("nodes_pch", 19) %>%   
         set("nodes_col", nodes_col) %>% 
         set("labels_cex", 0.5) %>%
         set("leaves_col", leaf_col[hierarchical_clustering$order]), 
       horiz = TRUE,
       main = main)
  
}

PrepareDendrogram <- function(features, #results from ParseSumsRatios
                              groupVector) #a vector containing which row fits in which group
{
  if (!is.factor(groupVector)) groupVector <- factor(groupVector)
  if (nrow(features) != length(groupVector)) {
    stop("Warning groupVector should be as long as the number of rows in features")
  }
  groups <- levels(groupVector)
  
  p_values <- apply(features, 2, function(feature){
    values_group1 <- feature[groupVector == groups[1]]
    values_group2 <- feature[groupVector == groups[2]]
    if (! (all(is.na(values_group1)) | all(is.na(values_group2)))) {
      test_res <- wilcox.test(values_group1, values_group2, exact = FALSE)
      test_res$p.value
    } else {
      1
    }
  })
  
  fold_changes <- apply(features, 2, function(feature){
    medians <- c(median(feature[groupVector == groups[1]], na.rm = TRUE) + 1e-100,
                 median(feature[groupVector == groups[2]], na.rm = TRUE) + 1e-100)
    fold_change <- max(medians) / min(medians) * (-1) ^ which.max(medians)
    if (length(fold_change) == 0) fold_change <- NA
    return(fold_change)
  })
  
  group_name <- paste0(groups[1], " vs ", groups[2])
  to_plot <- data.frame("feature" = colnames(features),
                        "pvalue" = p_values,
                        "-log10p" = -log10(p_values),
                        "foldchange" = fold_changes,
                        "log10foldchange" = sign(fold_changes)*log10(abs(fold_changes)), #Double check
                        "group" = group_name,
                        "subset" = rep("feature", ncol(features)),
                        check.names = FALSE)
  
  dendrogram_values <- list()
  dendrogram_values[[paste0(group_name, " feature", " -log10(p)")]] <- -log10(p_values[to_plot$subset == "feature"])
  dendrogram_values[[paste0(group_name, " feature", " log10foldchange")]] <- sign(fold_changes[to_plot$subset == "feature"]) *
    log10(abs(fold_changes[to_plot$subset == "feature"]))
  dendro_outlier <- dendrogram_values[[paste0(group_name, " feature", " log10foldchange")]] > 3
  dendrogram_values[[paste0(group_name, " feature", " log10foldchange")]][dendro_outlier] <- 3
  
  dendro_outlier <- dendrogram_values[[paste0(group_name, " feature", " log10foldchange")]] < -3
  dendrogram_values[[paste0(group_name, " feature", " log10foldchange")]][dendro_outlier] <- -3
  return(dendrogram_values)
}

plotColor <- function(colorVector){
  plotdf <- data.frame(x = colorVector, y = 1)
  plotdf$x <- factor(plotdf$x, levels = plotdf$x)
  ggplot2::ggplot(plotdf) +
    geom_bar(aes(x = x, y = y, fill = x), color = "black", stat = "identity") +
    scale_fill_manual(values = colorVector) +
    theme_minimal() + 
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank()) +
  coord_fixed() +
    ggtitle("Colors")
}

ggplotAsList <- function(p){
  tree <- p + theme(legend.position = "none")
  l <- get_legend(p)
  l1 <- as_ggplot(l$grobs[[1]])
  l2 <- as_ggplot(l$grobs[[2]])
  return(list("tree" = tree, "l1" = l1, "l2" = l2))
}

calculatePercentiles <- function(FCSfiles, percentiles, channels = NULL){
  if (is.null(channels)) channels <- colnames(exprs(read.FCS(FCSfiles[1])))
  ff <- read.FCS(FCSfiles[1])
  channels <- FlowSOM::GetChannels(ff, channels)
  names(channels) <- channels
  allPercentiles <- lapply(FCSfiles, function(fcsFile){
    ff_tmp <- read.FCS(fcsFile)
    exprs_mat <- exprs(ff_tmp)[, channels, drop = FALSE]
    percentilesPerFile <- lapply(channels, function(channel){
      return(quantile(exprs_mat[, channel], percentiles))
    })
  })
  results <- lapply(channels, function(channel){
    pcts_tmp <- do.call(rbind, lapply(allPercentiles, function(percentilePerFile){
      percentilePerFile[[channel]]
    }))
    return(apply(pcts_tmp, 2, median))
  })
  return(results)
}


doubletFinder2 <- function(ff, extraSpace = 0.02, channels = c("FSC-A", "FSC-H"), recursive = TRUE){
  pc <- prcomp(scale(exprs(ff)[, channels]))
  minY <- quantile(pc$x[, 2], 0.01)
  maxY <- quantile(pc$x[, 2], 0.99)
  if (abs(minY) > abs(maxY)) {
    cutoff <- -maxY
    isSinglets <- pc$x[, 2] > cutoff  - extraSpace
  } else {
    cutoff <- -minY
    isSinglets <- pc$x[, 2] < cutoff + extraSpace
  }
  if (recursive == TRUE){
    isSinglets2 <- doubletFinder(ff[isSinglets, ], recursive = FALSE)
    isSinglets[isSinglets] <- isSinglets2
  } 
  # pc_semiclean <- prcomp(scale(exprs(ff)[isSinglets, c("FSC-A", "FSC-H")]))
  # minY <- quantile(pc_semiclean$x[, 2], 0.01)
  # maxY <- quantile(pc_semiclean$x[, 2], 0.99)
  # if (abs(minY) > abs(maxY)) {
  #   cutoff <- -maxY
  #   isSinglets2 <- pc_semiclean$x[, 2] > cutoff - extraSpace
  # } else {
  #   cutoff <- -minY
  #   isSinglets2 <- pc_semiclean$x[, 2] < cutoff + extraSpace
  # }
  
  return(isSinglets) 
}



