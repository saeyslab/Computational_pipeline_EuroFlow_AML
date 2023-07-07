############################################################################## |
# LIBRARIES
############################################################################### |
library(flowCore)
library(FlowSOM)
library(Rtsne)
library(dplyr)
library(flowDensity)
library(ggplot2)
library(PeacoQC)
library(caret)
library(tidyr)
library(openxlsx)
library(randomForest)
library(caret)
library(stats)
library(parallel)
source("./3_helperFunctions.R")


############################################################################### /
# SCRIPT
############################################################################### /
#----------------------------------Define variables-----------------------------
seed <- 2022
set.seed(seed)

input_metadata <- "./full_database_02082022.xlsx"

output_norm_FCS_files <- "../../../group/irc/personal/artuurc/blasts/Normalized_fcs/"
output_RDS <- "../../../group/irc/personal/artuurc/blasts/Results/"
mc.cores <- 24

#--------------------------------------Functions--------------------------------
GetFeaturesPar <- function(fsom,
                           files,
                           level = c("clusters", "metaclusters"),
                           type = "counts",
                           MFI = NULL,
                           positive_cutoffs = NULL,
                           filenames = NULL,
                           silent = FALSE,
                           mc.cores = parallel::detectCores()) {
  nclus <- NClusters(fsom)
  nfiles <- length(files)
  i <- 0

  #----Warnings----
  if (!is.null(filenames) & length(filenames) != nfiles) {
    stop("Filenames vector should have same length as files vector.")
  }

  if (sum(level %in% c("clusters", "metaclusters")) != length(level)) {
    stop("level should be \"clusters\" and/or \"metaclusters\".")
  }

  if (sum(type %in% c("counts", "percentages", "MFIs", "percentages_positive")) != length(type)) {
    stop("level should be \"counts\", \"percentages\", \"MFIs\" and/or \"percentages_positive\".")
  }

  if ("MFIs" %in% type & is.null(MFI)) {
    stop("Please provide channel names for MFI calculation")
  }

  if ("percentages_positive" %in% type & is.null(positive_cutoffs)) {
    stop("Please provide positive cutoffs for percentage-positive calculation")
  }

  if (!is.null(positive_cutoffs) & length(names(positive_cutoffs)) == 0) {
    stop("positive_cutoffs must be a named vector")
  }

  matrices <- list()

  #----Prepare variables----
  if (is.null(filenames)) {
    if (is.character(files)) {
      filenames <- files
    } else {
      filenames <- as.character(seq_len(length(files)))
    }
  }

  if ("MFIs" %in% type) {
    nmetaclus <- NMetaclusters(fsom)
    MFI <- GetChannels(fsom, MFI)
    nmarker <- length(MFI)
  }

  if ("percentages_positive" %in% type) {
    nmetaclus <- NMetaclusters(fsom)
    positive_cutoffs <- unlist(positive_cutoffs)
    perc_pos <- GetChannels(fsom, names(positive_cutoffs))
    nmarker <- length(perc_pos)
  }
  #----Loop over files----
  featuresPerPatient <- parallel::mclapply(seq_along(files), function(fileI) {
    file <- files[[fileI]]
    fileName <- filenames[fileI]
    i <- i + 1
    if (isFALSE(silent)) {
      message(paste0("Mapping file ", i, " of ", nfiles, "."))
    }
    resList <- list()
    fsom_tmp <- suppressWarnings(FlowSOM::NewData(
      fsom = fsom,
      input = file,
      silent = silent
    ))
    C_counts <- C_outliers <- matrix(0,
      nrow = 1, ncol = FlowSOM::NClusters(fsom_tmp) + 1,
      dimnames = list(NULL, c("fileName", paste0("C", seq_len(FlowSOM::NClusters(fsom_tmp)))))
    )
    C_counts[1, 1] <- C_outliers[1, 1] <- fileName
    counts_t <- table(FlowSOM::GetClusters(fsom_tmp))
    C_counts[1, paste0("C", names(counts_t))] <- counts_t
    outliers_t <- fsom_tmp$outliers[, "Number_of_outliers", drop = FALSE]
    resList[["C_counts"]] <- C_counts
    if (nrow(outliers_t) != 0) {
      C_outliers[, paste0("C", rownames(outliers_t))] <-
        outliers_t$Number_of_outliers
      resList[["C_outliers"]] <- C_outliers
    }

    if ("MFIs" %in% type) {
      if ("clusters" %in% level) {
        C_MFIs <- as.vector(t(FlowSOM::GetClusterMFIs(fsom_tmp)[, MFI]))
        names(C_MFIs) <- paste0(
          rep(paste0("C", seq_len(nclus)),
            each = nmarker
          ),
          " ", fsom$prettyColnames[MFI]
        )
        resList[["C_MFIs"]] <- c(fileName, C_MFIs)
      }
      if ("metaclusters" %in% level) {
        MC_MFIs <- as.vector(t(FlowSOM::GetMetaclusterMFIs(fsom_tmp)[, MFI]))
        names(MC_MFIs) <- paste0(
          rep(
            paste0(
              "MC",
              seq_len(nmetaclus)
            ),
            each = nmarker
          ),
          " ", fsom$prettyColnames[MFI]
        )
        resList[["MC_MFIs"]] <- c(fileName, MC_MFIs)
      }
    }

    if ("percentages_positive" %in% type) {
      if ("clusters" %in% level) {
        C_perc_pos <- as.vector(t(FlowSOM::GetClusterPercentagesPositive(fsom_tmp, cutoffs = positive_cutoffs)))
        names(C_perc_pos) <- paste0(
          rep(paste0("C", seq_len(nclus)),
            each = nmarker
          ),
          " ", fsom$prettyColnames[perc_pos]
        )

        resList[["C_perc_pos"]] <- c(fileName, C_perc_pos)
      }
      if ("metaclusters" %in% level) {
        MC_perc_pos <- as.vector(t(FlowSOM::GetMetaclusterPercentagesPositive(fsom_tmp, cutoffs = positive_cutoffs)))
        names(MC_perc_pos) <- paste0(
          rep(
            paste0(
              "MC",
              seq_len(nmetaclus)
            ),
            each = nmarker
          ),
          " ", fsom$prettyColnames[perc_pos]
        )

        resList[["MC_perc_pos"]] <- c(fileName, MC_perc_pos)
      }
    }
    return(resList)
  }, mc.cores = mc.cores)

  results <- lapply(names(featuresPerPatient[[1]]), function(m) {
    df <- lapply(featuresPerPatient, function(p) {
      p[[m]]
    })
    df <- do.call(rbind, df) %>%
      as.data.frame()

    rownames_df <- df[, 1]
    df[, 1] <- NULL
    df <- apply(df, 2, as.numeric)
    rownames(df) <- rownames_df
    return(df)
  })
  names(results) <- names(featuresPerPatient[[1]])

  C_counts <- results$C_counts
  C_outliers <- results$C_outliers

  if ("MFIs" %in% type) {
    C_MFIs <- results$C_MFIs
    MC_MFIs <- results$MC_MFIs
  }

  if ("percentages_positive" %in% type) {
    C_perc_pos <- results$C_perc_pos
    MC_perc_pos <- results$MC_perc_pos
  }

  #----Add matrices to list----
  if ("clusters" %in% level) {
    if ("counts" %in% type) {
      C_counts_tmp <- C_counts
      attr(C_counts_tmp, "outliers") <- C_outliers
      matrices[["cluster_counts"]] <- C_counts_tmp
    }
    if ("percentages" %in% type) {
      C_pctgs <- prop.table(C_counts, margin = 1)
      colnames(C_pctgs) <- paste0("%", colnames(C_pctgs))
      attr(C_pctgs, "outliers") <- C_outliers
      matrices[["cluster_percentages"]] <- C_pctgs
    }
    if ("MFIs" %in% type) {
      matrices[["cluster_MFIs"]] <- C_MFIs
    }
    if ("percentages_positive" %in% type) {
      matrices[["cluster_percentages_positive"]] <- C_perc_pos
    }
  }

  if ("metaclusters" %in% level) {
    MC_counts <- t(apply(
      C_counts,
      1,
      function(x) {
        tapply(x, fsom$metaclustering, sum)
      }
    ))
    MC_counts[is.na(MC_counts)] <- 0
    colnames(MC_counts) <- paste0("MC", colnames(MC_counts))

    if ("counts" %in% type) {
      matrices[["metacluster_counts"]] <- MC_counts
    }
    if ("percentages" %in% type) {
      MC_pctgs <- prop.table(MC_counts, margin = 1)
      colnames(MC_pctgs) <- paste0("%", colnames(MC_pctgs))
      matrices[["metacluster_percentages"]] <- MC_pctgs
    }
    if ("MFIs" %in% type) {
      matrices[["metacluster_MFIs"]] <- MC_MFIs
    }
    if ("percentages_positive" %in% type) {
      matrices[["metacluster_percentages_positive"]] <- MC_perc_pos
    }
  }

  return(matrices)
}

getPositiveMFIs <- function(MFIs, metaclusterLabels, mc.cores = parallel::detectCores()) {
  positiveMFIs <- parallel::mclapply(metaclusterLabels, function(label) {
    MCL <- paste0("MC", stringr::str_extract(label, "[0-9]*"))
    label_split <- unlist(stringr::str_split(label, ", "))
    non_scatters_labels <- label_split[!grepl("FSC|SSC", label_split)]
    pos_labels <- non_scatters_labels[grepl("\\+", non_scatters_labels)]
    if (length(pos_labels) != 0) {
      pos_labels <- substr(pos_labels, 0, nchar(pos_labels) - 1)
      pos_labels <- paste0(MCL, " ", pos_labels)
      return(pos_labels)
    } else {
      return(NULL)
    }
  }, mc.cores = mc.cores)
  positiveMFIs <- do.call(c, positiveMFIs)
  cleanedColnames <- gsub(" <.*>", "", colnames(MFIs))
  whichCol <- c(
    match(positiveMFIs, cleanedColnames),
    grep("FSC|SSC", cleanedColnames)
  )
  resultMFIs <- MFIs[, whichCol]
  return(resultMFIs)
}

#--------------------------------Read metadata----------------------------------
metadata <- openxlsx::read.xlsx(input_metadata,
  check.names = FALSE,
  sep.names = "_"
)

for (x in colnames(metadata)[grepl("date|Date", colnames(metadata))]) {
  y <- convDate(metadata[, x])
  metadata[, x] <- y
}

non_intensive_chemotherapy_patients <- metadata %>%
  filter(treatment_modality != "intensive chemotherapy") %>%
  pull(PatientID)
noAML_patients <- c(
  "1861", "1870", "1873", "1901", "1903", "1914", "1919",
  "1922", "1923", "1934", "1939", "1943", "1947", "0011"
)
relapse_patients <- c("1925")
noInfo_patients <- c(
  "1908", "1910", "1911", "1912", "1928", "1942", "1946",
  "1734"
)
tooYoung_patients <- c("1866", "1868", "1869", "1605", "1906")
secAML_or_otherDiagnosis_patients <- c(
  "1616", "1723", "1812", "1904", "1740",
  "1711", "1879"
)
strangeFlowData_patient <- c("1951")
patientsToRemove <- c(
  noAML_patients, noInfo_patients, tooYoung_patients,
  secAML_or_otherDiagnosis_patients, strangeFlowData_patient,
  relapse_patients, non_intensive_chemotherapy_patients
)
metadata <- metadata %>% filter(!PatientID %in% patientsToRemove)
#------------------------Parse normalized files---------------------------------
fcs_files_norm <- paste0(
  output_norm_FCS_files,
  list.files(output_norm_FCS_files, pattern = ".fcs$")
)

samples_norm <- data.frame(
  File = fcs_files_norm,
  Panel = stringr::str_match(fcs_files_norm, "ALOT|AML")[, 1],
  PatientID = stringr::str_match(fcs_files_norm, "([0-9]+-?[0-9])")[, 2],
  Tube = as.numeric(stringr::str_match(fcs_files_norm, "AML([1-7]*)")[, 2]),
  stringsAsFactors = FALSE
) %>%
  filter(!PatientID %in% patientsToRemove)

samples_norm$Tube[is.na(samples_norm$Tube)] <- 0
tmp <- factor(substr(samples_norm$PatientID, 1, 2))
levels(tmp) <- c("2015", "2016", "2017", "2018", "2019")
samples_norm$Year <- tmp
samples_norm <- samples_norm %>% # filter out because double file
  dplyr::filter(!grepl("1747_ALOT_PBO_norm.fcs", File)) %>%
  inner_join(metadata[, c("PatientID", "date_flow_acquisition_ALOT")], by = "PatientID") %>%
  arrange(date_flow_acquisition_ALOT)

#-------FlowSOM and feature extraction in five-fold-cross validation---------
gate_values_all <- readRDS("./gate_values_all_updated.rds")
set.seed(seed)
all_features <- list()
for (tube in seq(0, 6)) {
  print(tube)
  samples_tube <- samples_norm %>%
    dplyr::filter(Tube == tube) %>%
    pull(File)
  saveRDS(samples_tube, paste0(output_RDS, "samples_order_tube_", tube, ".RDS"))
  ff_agg <- AggregateFlowFrames(samples_tube,
    cTotal = 3000000,
    keepOrder = TRUE, writeOutput = TRUE,
    silent = TRUE,
    outputFile = paste0(
      output_RDS, "agg_fcs_all_tube_",
      tube, ".fcs"
    )
  )
  fsom <- FlowSOM(ff_agg,
    colsToUse = colnames(ff_agg)[
      which(!colnames(ff_agg) %in% c(
        "Time",
        "File",
        "File_scattered",
        "Original_ID",
        "FSC-H"
      ))
    ],
    scale = FALSE,
    xdim = 20,
    ydim = 20,
    nClus = 30,
    seed = seed
  )
  # Split up metaclusters based on marker cutoffs
  C_centers <- GetClusterMFIs(fsom = fsom, prettyColnames = TRUE)

  gating_matrix <- matrix(
    data = NA,
    nrow = nrow(C_centers), ncol = length(fsom$map$colsUsed),
    dimnames = list(
      1:NClusters(fsom),
      colnames(C_centers)[match(
        fsom$map$colsUsed,
        colnames(fsom$data)
      )]
    )
  )
  gate_values <- gate_values_all[[as.character(tube)]]

  for (row in rownames(gating_matrix)) {
    for (col in colnames(gating_matrix)) {
      gating_matrix[row, col] <- C_centers[row, col] > gate_values[[col]]
    }
  }

  MC_labels <- c()
  for (row in rownames(gating_matrix)) {
    label <- ""
    for (col in colnames(gating_matrix)) {
      marker <- sub("(.*) <.*", "\\1", col)
      label <- paste0(label, ifelse(test = gating_matrix[row, col],
        yes = paste0(marker, "+, "),
        no = paste0(marker, "-, ")
      ))
    }
    MC_labels <- c(MC_labels, sub("(.*), $", "\\1", label))
  }

  MC_labels <- paste0(match(MC_labels, sort(unique(MC_labels))), ": ", MC_labels)

  sorted_MClabels <- gtools::mixedsort(unique(MC_labels))
  fsom_fullMC <- UpdateMetaclusters(
    fsom = fsom,
    clusterAssignment = MC_labels,
    levelOrder = sorted_MClabels
  )
  MC_df <- as.data.frame(sorted_MClabels)
  openxlsx::write.xlsx(MC_df, file = paste0(output_RDS, "MC_fullnames_tube_", tube, ".xlsx"))

  MC_labels_short <- sub("([0-9]*).*", "\\1", MC_labels)
  fsom_short <- UpdateMetaclusters(
    fsom = fsom,
    clusterAssignment = MC_labels_short,
    levelOrder = gtools::mixedsort(unique(MC_labels_short))
  )
  saveRDS(fsom_fullMC, paste0(
    output_RDS, "fsom_agg_toMapFiles_all_fullNames_tube_",
    tube, ".RDS"
  ))

  saveRDS(fsom_short, paste0(
    output_RDS, "fsom_agg_toMapFiles_all_shortNames_tube_",
    tube, ".RDS"
  ))
  FlowSOMmary(fsom_short, plotFile = paste0(
    output_RDS, "FlowSOMmary_all_tube_",
    tube, ".pdf"
  ))
  # Feature extraction
  fsom_agg <- fsom_short
  features <- GetFeaturesPar(fsom_agg,
    files = samples_tube,
    level = c("metaclusters"),
    type = c("percentages"), silent = TRUE,
    mc.cores = mc.cores
  )

  saveRDS(features, paste0(
    output_RDS, "unparsed_tube_features_all_tube_",
    tube, ".RDS"
  ))

  ## Parse features tube
  features_tube <- as.data.frame(features$metacluster_percentages)
  tmp <- gsub(" \\+ ", "p", colnames(features_tube))
  tmp <- gsub("\\Q ) / ( \\E", "d", tmp)
  tmp <- gsub("%", "X", tmp)
  tmp <- gsub("Cl", "", tmp)
  tmp <- gsub("\\( (.*) \\)", "\\1", tmp)
  tmp <- gsub(" <.*>", "", tmp)
  tmp <- gsub(" |-", "_", tmp)


  colnames(features_tube) <- paste0("Tube_", tube, "_", tmp)
  features_tube$PatientID <- stringr::str_extract(
    rownames(features_tube), "[0-9]{4}-?[0-9]?"
  )
  saveRDS(features_tube, paste0(output_RDS, "parsed_tube_features_all_tube_", tube, ".RDS"))
  all_features[[as.character(tube)]] <- features_tube
}

# Parse all features
all_features <- all_features %>%
  Reduce(function(dtf1, dtf2) dplyr::full_join(dtf1, dtf2, by = "PatientID"), .)
saveRDS(all_features, paste0(output_RDS, "all_features_all.RDS"))

print("Done")