############################################################################### |
# LIBRARIES
############################################################################### |
library(flowCore)
library(FlowSOM)
library(xgboost)
library(dplyr)
library(flowDensity)
library(ggplot2)
library(ggpubr)
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

input_norm_FCS_files <- "../../../group/irc/personal/artuurc/blasts/Normalized_fcs/"

output_cross_validation <- "../../../group/irc/personal/artuurc/blasts/Results/CV/"
output_RDS <- "../../../group/irc/personal/artuurc/blasts/RDS/"

mc.cores <- 24

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
  input_norm_FCS_files,
  list.files(input_norm_FCS_files, pattern = ".*fcs$")
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

#---------FlowSOM and feature extraction in five-fold-cross validation----------
patientIDs <- unique(samples_norm$PatientID)

gate_values_all <- readRDS("C:/Users/artuurc/Desktop/gate_values_all_updated.rds")

set.seed(seed)

handle <- qsub_lapply(
  X = patientIDs,
  qsub_environment = c("samples_norm", "getPositiveMFIs", "gate_values_all", "output_RDS", "seed"),
  qsub_config = qsub_config,
  FUN = function(part) {
    # Initialization
    library(FlowSOM)
    library(flowCore)
    library(dplyr)
    library(parallel)
    set.seed(seed)
    # Parse training and testing dataset
    df_train <- samples_norm %>%
      dplyr::filter(!PatientID %in% part)
    df_test <- samples_norm %>%
      dplyr::filter(PatientID %in% part)

    # Model building
    all_features <- list()
    for (tube in seq(0, 6)) {
      samples_train_tube <- df_train %>%
        dplyr::filter(Tube == tube) %>%
        pull(File)
      samples_test_tube <- df_test %>%
        dplyr::filter(Tube == tube) %>%
        pull(File)

      ff_agg <- AggregateFlowFrames(samples_train_tube,
        cTotal = 3000000,
        keepOrder = TRUE, writeOutput = FALSE,
        silent = TRUE,
        outputFile = paste0(
          output_RDS, "agg_fcs_cv_part_",
          part, "_tube_", tube, ".fcs"
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

      openxlsx::write.xlsx(MC_df, file = paste0(
        output_RDS, "MC_fullnames_tube_",
        tube, "part_", part, ".xlsx"
      ))
      saveRDS(fsom_fullMC, paste0(
        output_RDS, "fsom_agg_toMapFiles_cv_",
        part, "_fullNames_tube_", tube, ".RDS"
      ))

      MC_labels_short <- sub("([0-9]*).*", "\\1", MC_labels)
      fsom_short <- UpdateMetaclusters(
        fsom = fsom,
        clusterAssignment = MC_labels_short,
        levelOrder = gtools::mixedsort(unique(MC_labels_short))
      )
      # Feature extraction
      fsom_agg <- fsom_short
      features <- GetFeatures(fsom_agg,
        files = c(samples_train_tube, samples_test_tube),
        level = c("metaclusters"),
        type = c("percentages"), silent = TRUE
      )
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

      all_features[[as.character(tube)]] <- features_tube
    }
    print(paste0(output_RDS, "all_features_", part, ".RDS"))
    # Parse all features
    all_features <- all_features %>%
      Reduce(function(dtf1, dtf2) dplyr::full_join(dtf1, dtf2, by = "PatientID"), .)
    saveRDS(all_features, paste0(output_RDS, "all_features_", part, ".RDS"))
  }
)

#------------------Split patientIDs into outcome of interest--------------------
groupsList <- list(
  "relapse_vs_noRelapse_noTransplant" = list( # Relapse versus no relapse
    "group1" = metadata %>% dplyr::filter(grepl("relapse", relapse_.morphological.) &
      (transplant == "no" | date_of_relapse < date_of_transplantation)) %>%
      pull(PatientID),
    "group2" = metadata %>% dplyr::filter(relapse_.morphological. == "no" & transplant == "no") %>%
      pull(PatientID)
  ),
  "relapse_vs_noRelapse_transplant" = list( # alive 2y after diagnosis vs not alive 2y after diagnosis
    "group1" = metadata %>% dplyr::filter(grepl("relapse", relapse_.morphological.) &
      (transplant == "yes" & date_of_relapse > date_of_transplantation)) %>%
      pull(PatientID),
    "group2" = metadata %>% dplyr::filter(transplant == "yes" &
      relapse_.morphological. == "no") %>%
      pull(PatientID)
  ),
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
    "group1" = metadata %>% dplyr::filter(
      ELN_risk_classification_at_diagnosis == "favorable" &
        NPM1 == "present"
    ) %>%
      pull(PatientID),
    "group2" = metadata %>% dplyr::filter(ELN_risk_classification_at_diagnosis == "favorable" &
      NPM1 == "absent") %>%
      pull(PatientID)
  ),
  "inv16_vs_noInv16" = list(
    # alive 2y after diagnosis vs not alive 2y after diagnosis
    "group1" = metadata %>% dplyr::filter(
      ELN_risk_classification_at_diagnosis == "favorable" &
        inv16_at_diagnosis == "100"
    ) %>%
      pull(PatientID),
    "group2" = metadata %>% dplyr::filter(
      ELN_risk_classification_at_diagnosis == "favorable" &
        inv16_at_diagnosis == "absent"
    ) %>%
      pull(PatientID)
  )
)

#-----------Train and test models in five fold cross validation-----------------
perTube_and_forAll <- c(
  "all", "significant_3foldchange", "uncorrelated",
  "MCL_pctgs", "MCL_MFIs"
)
order_feature_selection <- c(
  "only_metadata",
  perTube_and_forAll,
  paste0(
    rep(paste0("tube_", seq(0, 6)), each = 5), "_",
    perTube_and_forAll
  )
)

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



# Initialize cross validation
plotList <- list()
NA_features_list <- list()
set.seed(seed)
for (compareGroups in names(groupsList)) {
  if (compareGroups == "favorable_vs_poor") {
    metadata_lim <- metadata_lim_tmp %>%
      select(-ELN_risk_classification_at_diagnosis)
  } else {
    metadata_lim <- metadata_lim_tmp
  }
  message(compareGroups)
  groupMetadata <- data.frame(
    PatientID = c(groupsList[[compareGroups]]$group1, groupsList[[compareGroups]]$group2),
    groups = factor(
      c(
        rep("group1", length(groupsList[[compareGroups]]$group1)),
        rep("group2", length(groupsList[[compareGroups]]$group2))
      ),
      levels = c("group1", "group2")
    )
  ) %>% filter(!is.na(groups))
  message("--Train and validate models")
  samples_norm_tmp <- samples_norm %>%
    dplyr::filter(PatientID %in% groupMetadata$PatientID)
  resultList <- mclapply(unique(samples_norm_tmp$PatientID), function(part) {
    print(part)
    all_features <- readRDS(paste0(output_RDS, "all_features_", part, ".RDS"))

    colnames(all_features)[grepl("\\+", colnames(all_features))] <- gsub(
      "\\+", "_", colnames(all_features)[grepl("\\+", colnames(all_features))]
    )

    df_train <- samples_norm_tmp %>%
      filter(!PatientID %in% part)
    df_test <- samples_norm_tmp %>%
      filter(PatientID %in% part)

    # Split up features in training and testing data set
    all_features <- all_features %>%
      filter(PatientID %in% groupMetadata$PatientID)

    all_features_test <- all_features %>%
      filter(PatientID %in% df_test$PatientID)

    all_features_train <- all_features %>%
      filter(PatientID %in% df_train$PatientID)

    # Feature Selection
    patientID_train <- all_features_train$PatientID
    patientID_test <- all_features_test$PatientID
    all_features_train$PatientID <- all_features_test$PatientID <- NULL
    selected_features_train <- list()
    selected_features_test <- list()

    ## All
    ### Train
    selected_features_train[["all"]] <- all_features_train
    ### Test
    selected_features_test[["all"]] <- all_features_test[, colnames(all_features_train)]

    ## Calculate statistically relevant features with 3x fold change
    ### Train
    rownames(all_features_train) <- patientID_train
    tmp <- getSign3FoldchangeFeatures(all_features_train, groupsList[[compareGroups]])
    selected_features_train[["significant_3foldchange"]] <- tmp$features

    ### Test
    selected_features_test[["significant_3foldchange"]] <- all_features_test[, colnames(tmp$features)]

    ## Select uncorrelated features
    ### Train
    selected_features_train[["uncorrelated"]] <- getFloremiFeatures(all_features_train,
      order = rownames(tmp$pValues),
      maxCor = 0.2
    )
    all_features_train <- all_features_train[, rownames(tmp$pValues)]

    ### Test
    selected_features_test[["uncorrelated"]] <- all_features_test[
      , colnames(selected_features_train[["uncorrelated"]])
    ]

    ## Select only MCL percentages
    ### Train
    selected_features_train[["MCL_pctgs"]] <- all_features_train[
      , grepl("XMC", colnames(all_features_train))
    ]

    ### Test
    selected_features_test[["MCL_pctgs"]] <- all_features_test[
      , grepl("XMC", colnames(all_features_test))
    ]

    ## Select only metacluster MFIs
    ### Train
    selected_features_train[["MCL_MFIs"]] <- all_features_train[
      grepl("_MC[0-9]", colnames(all_features_train))
    ]

    ### Test
    selected_features_test[["MCL_MFIs"]] <- all_features_test[
      grepl("_MC[0-9]", colnames(all_features_test))
    ]

    ## Select features per tube
    for (tube in seq(0, 6)) {
      print(tube)
      ### Select all features per tube
      #### Train
      all_features_train_tube <- all_features_train[
        ,
        grepl(paste0("Tube_", tube), colnames(all_features_train))
      ]
      selected_features_train[[paste0("tube_", tube, "_all")]] <- all_features_train_tube

      #### Test
      selected_features_test[[paste0("tube_", tube, "_all")]] <- all_features_test[
        , colnames(all_features_train_tube)
      ]

      ### Calculate statistically relevant features with 3x fold change per tube
      #### Train
      tmp <- getSign3FoldchangeFeatures(all_features_train_tube, groupsList[[compareGroups]])
      selected_features_train[[
        paste0("tube_", tube, "_significant_3foldchange")
      ]] <- tmp$features

      #### Test
      selected_features_test[[
        paste0("tube_", tube, "_significant_3foldchange")
      ]] <- all_features_test[
        , colnames(selected_features_train[[paste0("tube_", tube, "_significant_3foldchange")]])
      ]

      ### Select uncorrelated features per tube
      #### Train
      selected_features_train[[
        paste0("tube_", tube, "_uncorrelated")
      ]] <- getFloremiFeatures(
        all_features_train_tube,
        order = rownames(tmp$pValues),
        maxCor = 0.2
      )

      #### Test
      selected_features_test[[paste0("tube_", tube, "_uncorrelated")]] <- all_features_test[
        , colnames(selected_features_train[[paste0("tube_", tube, "_uncorrelated")]])
      ]

      ### Select only MCL percentages per tube
      #### Train
      selected_features_train[[paste0("tube_", tube, "_MCL_pctgs")]] <- all_features_train_tube[
        , grepl("XMC", colnames(all_features_train_tube))
      ]

      #### Test
      selected_features_test[[paste0("tube_", tube, "_MCL_pctgs")]] <- all_features_test[
        , colnames(selected_features_train[[paste0("tube_", tube, "_MCL_pctgs")]])
      ]

      ### Select only metacluster MFIs
      #### Train
      selected_features_train[[paste0("tube_", tube, "_MCL_MFIs")]] <- all_features_train_tube[
        , grepl("_MC[0-9]", colnames(all_features_train_tube))
      ]

      #### Test
      selected_features_test[[paste0("tube_", tube, "_MCL_MFIs")]] <- all_features_test[
        , colnames(selected_features_train[[paste0("tube_", tube, "_MCL_MFIs")]])
      ]
    }
    saveRDS(selected_features_train, paste0(
      output_RDS, "selected_features_train_",
      compareGroups, "_part_", part, ".RDS"
    ))
    saveRDS(selected_features_test, paste0(
      output_RDS, "selected_features_test_",
      compareGroups, "_part_", part, ".RDS"
    ))

    # Machine learning models
    ## Merge data with metadata
    full_data_train <- lapply(selected_features_train, function(x) {
      x$PatientID <- patientID_train
      x <- inner_join(metadata_lim, x, by = "PatientID")
      x <- inner_join(x, groupMetadata, by = "PatientID")
      return(x)
    })
    full_data_train$only_metadata <- inner_join(metadata_lim[metadata_lim$PatientID %in% patientID_train, ],
      groupMetadata,
      by = "PatientID"
    )

    full_data_test <- lapply(selected_features_test, function(x) {
      x$PatientID <- patientID_test
      x <- inner_join(metadata_lim, x, by = "PatientID")
      x <- inner_join(x, groupMetadata, by = "PatientID")
      return(x)
    })
    full_data_test$only_metadata <- inner_join(metadata_lim[metadata_lim$PatientID %in% patientID_test, ],
      groupMetadata,
      by = "PatientID"
    )
    ## Intialize lists
    featureSet_models <- wronglyPredicted_FS <- resList <- resList_featureSet <- list()
    for (featureSet in names(full_data_train)) {
      # Train ML models
      ## Intitialize lists
      selectedFeatures <- full_data_train[[featureSet]]
      patientID_train <- selectedFeatures$PatientID
      selectedFeatures$PatientID <- NULL
      if (ncol(selectedFeatures) <= 2) next

      ## Gradient Boosting Machines
      groups_GBM <- selectedFeatures$groups
      selectedFeatures_GBM <- selectedFeatures
      groups_GBM <- as.numeric(groups_GBM == "group1")
      selectedFeatures_GBM$groups <- NULL

      traindata <- xgboost::xgb.DMatrix(data = as.matrix(selectedFeatures_GBM), label = groups_GBM)
      gbm <- suppressWarnings(xgboost::xgb.train(
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

      featureSet_models[[featureSet]] <- gbm

      # Validate models on test set
      selectedFeatures_test <- full_data_test[[featureSet]]
      groups <- selectedFeatures_test$groups
      ref <- factor(groups, levels = c("group1", "group2"))
      names(ref) <- rownames(selectedFeatures_test) <- selectedFeatures_test$PatientID
      selectedFeatures_test$groups <- selectedFeatures_test$PatientID <- NULL

      ## Predict GBM
      testdata <- xgboost::xgb.DMatrix(data = as.matrix(selectedFeatures_test[, gbm$feature_names]))
      predictGBM <- predict(
        object = gbm,
        newdata = testdata
      )
      predictGBM <- ifelse(predictGBM > 0.5, "group1", "group2")
      names(predictGBM) <- rownames(selectedFeatures_test)

      ## Parse Results
      resList[[featureSet]] <- predictGBM
      wronglyPredicted_FS[[featureSet]] <- predictGBM[predictGBM != ref]
    }
    return(list(
      patientID = part,
      wronglyPredicted = wronglyPredicted_FS,
      predictions = resList,
      models = featureSet_models
    ))
  }, mc.cores = mc.cores)

  # Downstream analysis
  message("--Downstream analysis")

  patientIDs_resultList <- lapply(resultList, function(x) x[["patientID"]])
  listOfWronglyPredicted <- lapply(resultList, function(x) x[["wronglyPredicted"]])
  listOfModels <- lapply(resultList, function(x) x[["models"]])
  predictions_tmp <- lapply(resultList, function(x) x[["predictions"]])
  names(listOfWronglyPredicted) <- names(listOfModels) <- names(predictions_tmp) <- do.call(c, patientIDs_resultList)

  resList <- lapply(names(predictions_tmp[[1]]), function(x) {
    prediction <- lapply(names(predictions_tmp), function(y) {
      return(predictions_tmp[[y]][[x]])
    })
    prediction <- do.call(c, prediction)
    return(prediction)
  })
  names(resList) <- names(predictions_tmp[[1]])
  saveRDS(listOfWronglyPredicted, paste0(
    output_RDS, "wronglyPredicted_",
    compareGroups, ".RDS"
  ))
  saveRDS(listOfModels, paste0(
    output_RDS, "listOfModels_",
    compareGroups, ".RDS"
  ))
  saveRDS(resList, paste0(
    output_RDS, "listOfResults_",
    compareGroups, ".RDS"
  ))

  evaluationMetrics <- mclapply(names(resList), function(featureSet) {
    ## BA on all years together
    all_predictions <- resList[[featureSet]]
    df <- data.frame(
      PatientID = names(all_predictions),
      prediction = all_predictions
    )
    joined_dataset <- inner_join(groupMetadata, df, by = "PatientID")
    cm <- caret::confusionMatrix(
      reference = factor(joined_dataset$groups,
        levels = c("group1", "group2")
      ),
      data = factor(joined_dataset$prediction,
        levels = c("group1", "group2")
      )
    )

    return(data.frame(
      featureSet = featureSet,
      balanced_accuracies = cm[["byClass"]][["Balanced Accuracy"]],
      accuracies = cm[["overall"]][["Accuracy"]],
      precisions = cm[["byClass"]][["Pos Pred Value"]],
      recalls = cm[["byClass"]][["Sensitivity"]],
      F1 = cm[["byClass"]][["F1"]]
    ))
  }, mc.cores = mc.cores)

  ## Parse scores into ggplots
  scores <- do.call(rbind, evaluationMetrics) %>%
    as.data.frame() %>%
    tidyr::pivot_longer(-featureSet) %>%
    mutate(groups = ifelse(grepl("tube", featureSet),
      substr(featureSet, 1, 6),
      "X"
    )) %>%
    mutate(featureSet = factor(featureSet, levels = order_feature_selection))
  bestScores <- scores %>%
    filter(name == "balanced_accuracies") %>%
    arrange(desc(value)) %>%
    slice(1) %>%
    pull(value)
  p <- ggplot(scores) +
    geom_point(aes(
      x = featureSet,
      y = value,
      col = groups
    ), size = 2) +
    geom_hline(aes(yintercept = 0.5), linetype = "dashed") +
    xlab("") +
    geom_hline(yintercept = bestScores, color = "red", linetype = "dotted") +
    ggtitle(compareGroups) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    facet_wrap(~name, ncol = 1)
  saveRDS(scores, paste0(output_RDS, "/", "QualityMetrics_", compareGroups, ".RDS"))

  pdf(paste0(output_cross_validation, compareGroups, "_scores_per_model.pdf"),
    width = 10, height = 20
  )
  print(p)
  dev.off()
}
