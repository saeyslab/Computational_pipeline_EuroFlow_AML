############################################################################### /
# LIBRARIES
############################################################################### /
library(flowCore)
library(FlowSOM)
library(dplyr)
library(ggplot2)
library(treeshap)
library(shapviz)
library(caret)
library(randomForest)
library(xgboost)
library(ggpubr)
source("./3_helperFunctions.R")

############################################################################### /
# SCRIPT
############################################################################### /
#----------------------------------Define variables-----------------------------
seed <- 2022
set.seed(seed)

input_metadata <- "./Excel_metadata/"
input_RDS <- "./blasts/XGBoost/RDS/"


output_norm_FCS_files <- "./blasts/Normalized_fcs/"
output_downstream_analysis <- "./blasts/SHAP1/"

mc.cores <- 1

#--------------------------------Read metadata----------------------------------
metadata <- openxlsx::read.xlsx(paste0(input_metadata, "full_database_02082022.xlsx"),
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
  list.files(output_norm_FCS_files, pattern = ".*fcs$")
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
  dplyr::filter(File != "D:/PhD/EuroFlow/Normalized_fcs/1747_ALOT_PBO_norm.fcs") %>%
  inner_join(metadata[, c("PatientID", "date_flow_acquisition_ALOT")], by = "PatientID") %>%
  arrange(date_flow_acquisition_ALOT)

#------------------Split patientIDs into outcome of interest--------------------
groupsList <- list(
  "A2_vs_D2" = list( # alive 2y after diagnosis vs not alive 2y after diagnosis
    "group1" = metadata %>% dplyr::filter(Date_of_last_contact - date_of_diagnosis >= 365 * 2 &
      (survival_status_at_last_contact == "alive" |
        (survival_status_at_last_contact == "died" &
          date_of_death - date_of_diagnosis >= 365 * 2))) %>%
      pull(PatientID),
    "group2" = metadata %>% dplyr::filter(Date_of_last_contact - date_of_diagnosis < 365 * 2 &
      cause_of_death %in% c(
        "death from AML progressive disease",
        "death in aplasie",
        "death with (R/R) AML but from another cause"
      )) %>%
      pull(PatientID)
  ),
  "favorable_vs_poor" = list( # alive 2y after diagnosis vs not alive 2y after diagnosis
    "group1" = metadata %>% dplyr::filter(ELN_risk_classification_at_diagnosis == "favorable") %>%
      pull(PatientID),
    "group2" = metadata %>% dplyr::filter(ELN_risk_classification_at_diagnosis == "poor") %>%
      pull(PatientID)
  ),
  "NPM1_vs_noNPM1" = list(
    # Relapse versus no relapse
    "group1" = metadata %>% dplyr::filter(NPM1 == "present") %>%
      pull(PatientID),
    "group2" = metadata %>% dplyr::filter(NPM1 == "absent") %>%
      pull(PatientID)
  ),
  "inv16_vs_noInv16" = list(
    # alive 2y after diagnosis vs not alive 2y after diagnosis
    "group1" = metadata %>% dplyr::filter(inv16_at_diagnosis == "100") %>%
      pull(PatientID),
    "group2" = metadata %>% dplyr::filter(inv16_at_diagnosis == "absent") %>%
      pull(PatientID)
  )
)
#----------------------Train best models and SHAP analysis----------------------
if (!dir.exists(output_downstream_analysis)) dir.create(output_downstream_analysis)
# Parse metadata
metadata_lim_tmp <- metadata[, c(
  "PatientID", "gender",
  "ELN_risk_classification_at_diagnosis",
  "age_at_diagnosis", "WBC_count_at_diagnosis",
  "peripheral_blood_blast_percentage",
  "bone_marrow_blasts_percentage"
)]

metadata_lim_tmp$age2 <- metadata_lim_tmp$age^2
metadata_lim_tmp[, "gender"] <- as.numeric(as.character(
  factor(metadata_lim_tmp[, "gender"],
    levels = c("female", "male"),
    labels = c("0", "1")
  )
))

metadata_lim_tmp[, "ELN_risk_classification_at_diagnosis"] <- as.numeric(as.character(
  factor(metadata_lim_tmp[, "ELN_risk_classification_at_diagnosis"],
    levels = c("poor", "intermediate", "favorable"),
    labels = c("0", "1", "2")
  )
))

reset <- TRUE
set.seed(seed)
plotlist <- list()
for (compareGroups in names(groupsList)) {
  title <- if (compareGroups == "A2_vs_D2") {
    "survival status after 2y"
  } else if (compareGroups == "favorable_vs_poor") {
    "ELN"
  } else if (compareGroups == "NPM1_vs_noNPM1") {
    "NPM1"
  } else if (compareGroups == "inv16_vs_noInv16") {
    "inv(16)"
  }
  if (compareGroups == "favorable_vs_poor") {
    metadata_lim <- metadata_lim_tmp %>%
      select(-ELN_risk_classification_at_diagnosis)
  } else {
    metadata_lim <- metadata_lim_tmp
  }
  message(compareGroups)
  if (!dir.exists(paste0(output_downstream_analysis, compareGroups))) {
    dir.create(paste0(output_downstream_analysis, compareGroups))
  }
  # Initialize
  namesOfGroups <- unlist(strsplit(compareGroups, "_"))
  if (length(namesOfGroups) == 4) {
    nameOfGroup1 <- paste0(namesOfGroups[1], "_", namesOfGroups[4])
    nameOfGroup2 <- paste0(namesOfGroups[3], "_", namesOfGroups[4])
  } else {
    nameOfGroup1 <- namesOfGroups[1]
    nameOfGroup2 <- namesOfGroups[3]
  }
  wronglyPredicted <- readRDS(paste0(
    input_RDS, "wronglyPredicted_",
    compareGroups, ".RDS"
  ))
  # Determine best model and best average models
  SHAPModels <- readRDS(paste0(
    input_RDS, "QualityMetrics_",
    compareGroups, ".RDS"
  )) %>%
    filter(name == "balanced_accuracies") %>%
    arrange(desc(value)) %>%
    dplyr::slice(1) %>%
    select(featureSet, value, groups)

  # Loop over best (average) models

  featureSet <- featureSet_tmp <- SHAPModels$featureSet
  model <- "GBM"
  wronglyPredicted_OI <- do.call(c, lapply(wronglyPredicted, function(x) {
    names(x[[featureSet]])
  }))
  dir_tmp <- paste0(
    output_downstream_analysis, compareGroups, "/", featureSet,
    "_", model, "/"
  )
  print(basename(dir_tmp))
  if (!dir.exists(dir_tmp)) {
    dir.create(dir_tmp)
  }

  # Parse all features
  all_features <- all_features_bu <- readRDS("./blasts/Heatmaps_all_patients/all_features_all.RDS")
  patientids <- rownames(all_features) <- all_features$PatientID
  all_features$PatientID <- NULL

  tube <- stringr::str_extract(featureSet, "tube_[0-9]")
  if (!is.na(tube)) {
    all_features <- all_features[, grepl(tube, colnames(all_features), ignore.case = TRUE)]
    featureSet <- substr(featureSet, 8, nchar(as.character(featureSet)))
  }

  if (!file.exists(paste0(dir_tmp, "all_features.csv")) | reset) {
    # Select best features
    if (featureSet == "only_metadata") {
      all_features <- metadata_lim
      patientids <- metadata_lim$PatientID
    } else {
      if (featureSet %in% c("significant_3foldchange", "uncorrelated")) {
        tmp <- getSign3FoldchangeFeatures(all_features, groupsList[[compareGroups]])
        if (featureSet == "uncorrelated") {
          all_features <- getFloremiFeatures(all_features,
            order = rownames(tmp$pValues),
            maxCor = 0.2
          )
        } else {
          all_features <- tmp$features
        }
      } else if (featureSet == "MCL_pctgs") {
        all_features <- all_features[
          , grepl("XMC", colnames(all_features))
        ]
      }
      all_features$PatientID <- patientids
      all_features <- full_join(all_features, metadata_lim, by = "PatientID")
    }

    # Prepare features for model training
    groupMetadata <- rbind(
      data.frame(PatientID = groupsList[[compareGroups]]$group1, groups = "group1"),
      data.frame(PatientID = groupsList[[compareGroups]]$group2, groups = "group2")
    )
    all_features <- inner_join(groupMetadata, all_features, by = "PatientID")
    rownames(all_features) <- all_features$PatientID
    all_features$PatientID <- NULL
    all_features$groups <- factor(all_features$groups)
    write.csv(all_features, paste0(dir_tmp, "all_features.csv"))
  } else {
    all_features <- read.csv(paste0(dir_tmp, "all_features.csv"),
      stringsAsFactors = FALSE,
      colClasses = c(
        "groups" = "character",
        "X" = "character"
      ),
      row.names = 1
    )
  }

  colnames(all_features) <- gsub("\\+", "_", colnames(all_features))

  # Train model
  set.seed(seed)
  groups_GBM <- all_features$groups
  all_features_GBM <- all_features
  groups_GBM <- as.numeric(groups_GBM == "group1")
  all_features_GBM$groups <- NULL

  traindata <- xgboost::xgb.DMatrix(data = as.matrix(all_features_GBM), label = groups_GBM)
  resultModel <- suppressWarnings(xgboost::xgb.train(
    data = traindata,
    params = list(
      subsample = 0.8,
      eta = 0.05
    ),
    nrounds = 50,
    verbose = FALSE,
    objective = "binary:logistic",
    eval_metric = "logloss",
    nthread = 1
  ))
  saveRDS(resultModel, paste0(dir_tmp, featureSet_tmp, "_", model, ".RDS"))

  # SHAP analysis
  toShow <- 5
  groups <- all_features$groups
  all_features$groups <- NULL
  unified <- xgboost.unify(resultModel, all_features)
  base_value <- mean(predict(resultModel, data.matrix(all_features)))
  treeshap_values <- treeshap(unified, as.data.frame(all_features), verbose = 0)
  order <- treeshap_values$shaps %>%
    mutate(across(everything(), abs)) %>%
    summarize(across(everything(), mean)) %>%
    as.character() %>%
    order(decreasing = TRUE)
  breaks <- colnames(all_features[, order])[1:toShow]
  labels <- gsub("_", " ", breaks)
  labels <- gsub(" at diagnosis|Tube |X", "", labels)
  shp <- shapviz::shapviz(treeshap_values)
  p <- sv_importance(shp, kind = "beeswarm", max_display = toShow, size = 0.5) +
    theme_classic() +
    xlim(-2.5, 2.5) +
    scale_y_discrete(breaks = breaks, label = labels) +
    ggtitle(title) +
    theme(
      legend.position = "None",
      axis.text.x = element_text(size = 5),
      axis.text.y = element_text(size = 5),
      axis.line.x = element_line(linewidth = 0.3),
      axis.line.y = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      axis.title = element_text(size = 5),
      title = element_text(size = 7)
    )
  plotlist[[compareGroups]] <- p
  write.csv(treeshap_values$shaps, paste0(dir_tmp, "shap_values_group1.csv"), row.names = FALSE)
  write.csv(rownames(all_features), paste0(dir_tmp, "patientIDs_group1.csv"))
  saveRDS(treeshap_values, paste0(dir_tmp, "shap_values_group1.RDS"))
  tmp_metadata <- metadata_lim
}
pdf("./blasts/SHAP/SHAP.pdf", width = 6.88, height = 4)
ggarrange(plotlist = plotlist, common.legend = TRUE, legend = "right", align = "hv")
dev.off()

