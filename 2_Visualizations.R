############################################################################### /
# LIBRARIES
############################################################################### /
library(flowCore)
library(FlowSOM)
library(dplyr)
library(ggplot2)
library(fastshap)
library(caret)
library(randomForest)
library(xgboost)
library(ComplexHeatmap)
library(tidycmprsk)
library(ggpubr)
library(ggsurvfit)
source("./3_helperFunctions.R")

############################################################################### /
# SCRIPT
############################################################################### /
#----------------------------------Define variables-----------------------------
seed <- 2022
set.seed(seed)

input_metadata <- "./Excel_metadata/"
input_norm_FCS_files <- "./blasts/Normalized_fcs/"
input_RDS <- "./blasts/Heatmaps_all_patients/"
output_heatmaps <- "./blasts/Heatmaps_all_patients/"

#--------------------------------Read metadata----------------------------------
metadata <-
  openxlsx::read.xlsx(
    paste0(input_metadata, "full_database_02082022.xlsx"),
    check.names = FALSE,
    sep.names = "_"
  )

for (x in colnames(metadata)[grepl("date|Date", colnames(metadata))]) {
  y <- convDate(metadata[, x])
  metadata[, x] <- y
}

non_intensive_chemotherapy_patients <- metadata %>%
  dplyr::filter(treatment_modality != "intensive chemotherapy") %>%
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
plots <- list(
  "0" = list(
    c("CD45", "SSC-A"),
    c("CD3", "cyCD3"), c("CD19", "cyCD79a"),
    c("CD34", "CD7"), c("FSC-A", "cyMPO")
  ),
  "1" = list(
    c("CD45", "SSC-A"),
    c("CD34", "CD117"), c("HLA-DR", "CD11b"),
    c("CD16", "CD13"), c("CD11b", "CD13"),
    c("CD16", "CD10")
  ),
  "2" = list(
    c("CD45", "SSC-A"),
    c("CD34", "CD117"), c("CD14", "CD64"),
    c("HLA-DR", "CD64"), c("CD14", "CD35"),
    c("CD300e", "CD35")
  ),
  "3" = list(
    c("CD45", "SSC-A"),
    c("CD34", "CD117"), c("HLA-DR", "CD33"),
    c("CD36", "CD71"), c("CD36", "CD105"),
    c("CD71", "CD105")
  ),
  "4" = list(
    c("CD45", "SSC-A"),
    c("CD34", "CD117"), c("CD34", "nuTdT"),
    c("CD7", "CD56"), c("CD7", "CD19"),
    c("CD117", "HLA-DR")
  ),
  "5" = list(
    c("CD45", "SSC-A"),
    c("CD34", "CD117"), c("CD34", "CD38"),
    c("CD22", "CD38"), c("CD15", "NG2"),
    c("CD117", "HLA-DR")
  ),
  "6" = list(
    c("CD45", "SSC-A"),
    c("CD34", "CD117"), c("CD123", "HLA-DR"),
    c("CD42a+CD61", "CD203c"), c("CD4", "CD123")
  )
)

subgroups_AML <- lapply(unique(metadata$AML_with_recurrent_genetic_abnormalities), function(x) {
  res <- na.omit(metadata$PatientID[metadata$AML_with_recurrent_genetic_abnormalities == x])
  if (is.na(x)) {
    res <- metadata$PatientID[is.na(metadata$AML_with_recurrent_genetic_abnormalities)]
  }
  return(res)
})
names(subgroups_AML) <- paste(unique(metadata$AML_with_recurrent_genetic_abnormalities))
subgroups_AML$`mutated NPM1` <- subgroups_AML$`NA` <- NULL
subgroups_AML <- c(subgroups_AML, list(
  "NPM1_A" = c( # CD117+HLADR+CD34+
    "1927", "1926", "19"
  ),
  "NPM1_B" = c( # CD117-HLADR-CD34- and CD117+HLADR-CD34-
    "1885", "0032", "1863", "1930"
  ),
  "NPM1_C" = c( # CD117-HLADR+CD34-
    "1741", "1949", "0026", "0008", "1730", "1708", "1710", "1802", "1924"
  ),
  "NPM1_D" = c( # CD117+HLADR-CD34-
    "0022", "1627", "0017", "0020", "1712", "1742", "1805", "1808", "1862", "0039",
    "1882", "1612"
  ),
  "NPM1_E" = c( # CD117+HLADR+CD34-
    "0044", "1623", "1940", "1929", "1945", "1722", "1721", "1617", "0023", "0041"
  ),
  "other_A" = c( # CD117+HLADR+CD34+
    "1875", "0028", "1815", "1602", "1909", "1871", "0021", "1878", "1804", "1731",
    "1749", "1626", "1735", "1816", "1714", "1738", "1952", "1743", "0035", "0031",
    "0025", "1603", "1750", "0029", "1705", "1739"
  ),
  "other_B" = c( # CD117-HLADR+CD34-
    "1716", "1877", "1936", "1813", "0040", "0027", "0037"
  ),
  "other_C" = c( # CD117+HLADR-CD34-
    "1752", "0030", "1864"
  ),
  "other_D" = c( # CD117+HLADR+CD34-
    "0033", "1707", "1745", "1953", "1608", "1881", "1609", "1916"
  )
))

subgroups_AML <- do.call(rbind, lapply(names(subgroups_AML), function(x) {
  data.frame("group" = x, "PatientID" = unlist(subgroups_AML[x]))
})) %>% as.data.frame()

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

#------------------------------Heatmap all features-----------------------------
all_features <- readRDS(paste0(input_RDS, "all_features_all.RDS"))
rownames(all_features) <- all_features$PatientID
all_features_pctgs <- all_features[, grep("XMC", colnames(all_features))]

hm_data <- t(as.matrix(all_features_pctgs))
rownames(hm_data) <- gsub("Tube_|X", "", rownames(hm_data))
emptyCols <- grep("1622|1801|1872|1905|0019|1917|1905|1754|1604", colnames(hm_data))
hm_data_tocluster <- hm_data[, -emptyCols]
clustering_cols <- hclust(dist(t(hm_data_tocluster)))$order
clustering_rows <- rownames(hm_data_tocluster)[hclust(dist(hm_data_tocluster))$order]
hm_data <- cbind(hm_data_tocluster[clustering_rows, clustering_cols], hm_data[clustering_rows, emptyCols])
col_fun <- circlize::colorRamp2(c(0, 0.25, 0.5, 1), c("white", "#79C0CE", "#054D48", "firebrick"))

metadata_to_show <- c(
  "PatientID", "NPM1", "FLT3ITD", "FLT3TKD", "CEBPalfa", "ASXL1", "RUNX1", "SF3B1", "SRSF2", "U2AF1",
  "DNMT3A", "IDH1", "IDH2", "TP53", "KIT", "TET2", "NRAS", "MLL.PTD", "BCR.ABL_at_diagnosis",
  "MLL.AF9_at_diagnosis", "pml_rara_at_diagnosis", "inv16_at_diagnosis", "AML1.ETO_at_diagnosis",
  "alive_2y_after_diagnosis", "transplant", "relapse_HSCT", "ELN_risk_classification_at_diagnosis",
  "cytogenetics", "WHO_classification", "AML_with_recurrent_genetic_abnormalities"
)
relapse_HSCT <- list(
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
  )
)

relapse_df <- rbind(
  data.frame(PatientID = relapse_HSCT$relapse_vs_noRelapse_noTransplant$group1, relapse_HSCT = "relapse after chemotherapy"),
  data.frame(PatientID = c(relapse_HSCT$relapse_vs_noRelapse_noTransplant$group2, relapse_HSCT$relapse_vs_noRelapse_transplant$group2), relapse_HSCT = "no relapse"),
  data.frame(PatientID = relapse_HSCT$relapse_vs_noRelapse_transplant$group1, relapse_HSCT = "relapse after transplant"),
  data.frame(PatientID = metadata %>% filter(relapse_.morphological. == "primary refractory") %>% pull(PatientID), relapse_HSCT = "primary refractory")
)


df_ha <- metadata %>%
  left_join(relapse_df, by = "PatientID") %>%
  select(all_of(metadata_to_show)) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("PatientID") %>%
  mutate(across(everything(), function(x) ifelse(x == "not available", NA, x))) %>%
  mutate(pml_rara_at_diagnosis = ifelse(pml_rara_at_diagnosis %in% c("100", "9.4E-2"), "+", pml_rara_at_diagnosis)) %>%
  mutate(inv16_at_diagnosis = ifelse(inv16_at_diagnosis == "100", "+", inv16_at_diagnosis)) %>%
  mutate(FLT3ITD = ifelse(grepl("present", FLT3ITD), "+", FLT3ITD)) %>%
  mutate(FLT3ITD = ifelse(FLT3ITD == "other", NA, FLT3ITD)) %>%
  mutate(FLT3TKD = ifelse(grepl("present", FLT3TKD), "+", FLT3TKD)) %>%
  mutate(CEBPalfa = ifelse(grepl("present", CEBPalfa), "+", CEBPalfa)) %>%
  mutate(BCR.ABL_at_diagnosis = ifelse(grepl("100", BCR.ABL_at_diagnosis), "+", BCR.ABL_at_diagnosis)) %>%
  mutate(MLL.AF9_at_diagnosis = ifelse(MLL.AF9_at_diagnosis == 100, "+", MLL.AF9_at_diagnosis)) %>%
  mutate(across(everything(), function(x) ifelse(x %in% c("present", "yes"), "+", x))) %>%
  mutate(across(everything(), function(x) ifelse(x %in% c("absent", "no"), "-", x)))

colors <- c(
  "+" = "#FFF599", "-" = "white",
  "inv(16), t(16;16)" = "#0074D9", "normal karyotype" = "#F012BE",
  "other" = "#B10DC9", "PML-RARA = t(15;17)" = "#85144b",
  "t(8;21) = RUNX1-RUNXT (AML1-ETO)" = "#3D9970", "inv(16)" = "#FFC300",
  "mutated NPM1" = "#00C2CB", "mutated RUNX1" = "#4B0082",
  "t(15;17)" = "#FF6E7F", "t(8;21)" = "#8B4513", "BCR-ABL t(9;22)" = "#FF4136",
  "RPN1-EVI1 = inv(3), t(3;3)" = "#39CCCC", "del7, del7q" = "#2ECC40",
  "complex karyotype" = "#FF851B", "monosomal karyotype" = "#7FDBFF",
  "t(v;11q23) MLL rearrang, excl t(9;11)" = "#AAAAAA", "del5, del5q" = "#FFDC00",
  "MLL = t(9;11) (p21-22;q23)" = "#001f3f", "BCR-ABL1" = "#FF5733",
  "inv(3)" = "#DAF7A6", "favorable" = "#7CEA9C", "intermediate" = "#F1BF98",
  "poor" = "#E84855", "with recurrent genetic abnormalities" = "#623B5A", "therapy-related" = "#BA9593",
  "with MDS-related changes" = "#89608E", "ambiguous lineage" = "#C8FFBE", "NOS" = "#EDFFAB"
)

cols <- lapply(metadata_to_show, function(x) {
  if (!x %in% c("relapse_HSCT", "alive_2y_after_diagnosis", "transplant")) {
    return(colors)
  } else {
    return(c(
      "+" = "#CDFCC5", "-" = "white",
      "relapse after chemotherapy" = "#1C3144", "primary refractory" = "#D00000", "relapse after transplant" = "#FFBA08", "no relapse" = "#FFFFFF"
    ))
  }
})
names(cols) <- metadata_to_show
size <- 0.15
annotation_labels <- gsub("_", " ", colnames(df_ha))
annotation_labels <- gsub("at diagnosis", " ", annotation_labels)
annotation_labels <- gsub("alive 2y after diagnosis", "survival status 2y", annotation_labels)
annotation_labels <- gsub("AML1.ETO  ", "AML1-ETO", annotation_labels)
annotation_labels <- gsub("inv16  ", "inv(16)", annotation_labels)
annotation_labels <- gsub("pml rara  ", "PML-RARA", annotation_labels)
annotation_labels <- gsub("MLL.AF9  ", "MLL-AF9", annotation_labels)
annotation_labels <- gsub("BCR.ABL  ", "BCR-ABL1", annotation_labels)
annotation_labels <- gsub("MLL.PTD", "MLL-PTD", annotation_labels)
annotation_labels <- gsub("CEBPalfa", "CEBPA", annotation_labels)
annotation_labels <- gsub("FLT3TKD", "FLT3-TKD", annotation_labels)
annotation_labels <- gsub("FLT3ITD", "FLT3-ITD", annotation_labels)
annotation_labels <- gsub("AML with recurrent genetic abnormalities", "genetic abnormalities", annotation_labels)
ha <- columnAnnotation(
  df = df_ha[colnames(hm_data), ],
  col = cols,
  gp = gpar(col = "grey80", lwd = 0.5),
  annotation_name_gp = gpar(fontsize = 5),
  annotation_label = annotation_labels,
  annotation_legend_param = list(
    "NPM1" = list(title = "Genetic abnormalities"),
    "transplant" = list(title = "Outcome")
  ),
  show_legend = c(TRUE, rep(FALSE, 22), rep(TRUE, 5)),
  simple_anno_size = unit(size, "cm"),
  na_col = "grey90"
)


column_split <- metadata[, c("PatientID", "AML_with_recurrent_genetic_abnormalities")]
column_split <- tibble::deframe(column_split)[colnames(hm_data)]
column_split[is.na(column_split)] <- "other"


rn_hm <- gsub("_", " ", rownames(hm_data))
MCs <- stringr::str_extract(rn_hm, "(?<=MC)[0-9]*")
rn_hm[nchar(MCs) == 1] <- gsub("MC([0-9])", "MC0\\1", rn_hm[nchar(MCs) == 1])
hm1 <- Heatmap(hm_data,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_fun,
  bottom_annotation = ha,
  column_split = column_split,
  row_labels = rn_hm,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  column_title_gp = gpar(fontsize = 6),
  na_col = "grey90",
  heatmap_legend_param = list(title = "Metacluster percentage"),
  rect_gp = gpar(col = "grey80", lwd = 0.5),
  height = unit(nrow(hm_data) * size, "cm"),
  width = unit(ncol(hm_data) * size, "cm")
)



expression_matrices_tmp <- lapply(seq(0, 6), function(tube) {
  features_tube <- colnames(all_features_pctgs)[grepl(paste0("Tube_", tube), colnames(all_features_pctgs))]
  fsom_full_names <- readRDS(paste0(input_RDS, "fsom_agg_toMapFiles_all_fullNames_tube_", tube, ".RDS"))
  all_markers <- GetMarkers(fsom_full_names, fsom_full_names$map$colsUsed)
  expression_matrix_tube <- matrix(do.call(rbind, lapply(features_tube, function(f) {
    MC <- as.numeric(gsub("Tube_[0-9]_XMC", "", f))
    title <- levels(fsom_full_names$metaclustering)[MC]
    marker_profile <- unlist(stringr::str_split(gsub("[0-9]*: ", "", title), ", "))
    return(substr(marker_profile, nchar(marker_profile), nchar(marker_profile)))
  })), nrow = length(features_tube), ncol = length(all_markers), dimnames = list(features_tube, all_markers))
  return(as.data.frame(expression_matrix_tube))
})

stats <- lapply(names(groupsList), function(group) {
  if (group == "relapse_vs_noRelapse_noTransplant") {
    nameOfGroup1 <- "relapse_after_chemotherapy"
    nameOfGroup2 <- "no_relapse_after_chemotherapy"
  } else if (group == "relapse_vs_noRelapse_transplant") {
    nameOfGroup1 <- "relapse_after_transplant"
    nameOfGroup2 <- "no_relapse_after_transplant"
  } else {
    namesOfGroups <- unlist(strsplit(group, "_"))
    nameOfGroup1 <- namesOfGroups[1]
    nameOfGroup2 <- namesOfGroups[3]
  }
  stats <- GroupStats(t(hm_data), groups = groupsList[[group]]) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature") %>%
    select(feature, `p values`) %>%
    mutate(significance = case_when(
      `p values` < 0.001 ~ "***",
      `p values` < 0.01 ~ "**",
      `p values` < 0.05 ~ "*",
      TRUE ~ "ns"
    )) %>%
    select(feature, significance)
  colnames(stats) <- c("feature", nameOfGroup1)
  return(stats)
})
stats <- stats %>%
  Reduce(function(dtf1, dtf2) full_join(dtf1, dtf2, by = "feature"), .)
rownames(stats) <- stats$feature
stats$feature <- NULL

col_labels <- c(
  "relapse after chemotherapy", "relapse after transplant", "survival status 2y", "ELN",
  "NPM1", "inv(16)"
)
rn_hm_sign <- gsub("_", " ", rownames(hm_data))
MCs <- stringr::str_extract(rn_hm_sign, "(?<=MC)[0-9]*")
rn_hm_sign[nchar(MCs) == 1] <- gsub("MC([0-9])", "MC0\\1", rn_hm_sign[nchar(MCs) == 1])
hm_sign <- Heatmap(as.matrix(stats),
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  col = c("*" = "#FFE5F5", "**" = "#FF99D8", "***" = "#FF0AA1", "ns" = "white"),
  na_col = "grey90",
  row_labels = rn_hm_sign,
  column_labels = col_labels,
  heatmap_legend_param = list(
    title = "Significance",
    at = c("*", "**", "***", "ns"),
    labels = c("p < 0.05", "p < 0.01", "p < 0.001", "not significant")
  ),
  rect_gp = gpar(col = "grey80", lwd = 0.5),
  row_names_side = "left",
  height = unit(nrow(stats) * size, "cm"),
  width = unit(ncol(stats) * size, "cm")
)


colors <- c("-" = "white", "+" = "linen")
expression_matrices <- bind_rows(expression_matrices_tmp)
rownames(expression_matrices) <- gsub("Tube_|X", "", rownames(expression_matrices))
expression_matrices <- expression_matrices[, -grep("FSC|SSC|CD45", colnames(expression_matrices))]
hm2_data <- as.matrix(expression_matrices[clustering_rows, c(
  "CD117", "HLA-DR", "CD34", "cyMPO", "cyCD79a",
  "CD3", "cyCD3", "CD19", "CD7", "CD16", "CD13",
  "CD11b", "CD10", "CD35", "CD64", "CD300e", "CD14",
  "CD36", "CD105", "CD33", "CD71", "nuTdT", "CD56",
  "CD15", "NG2", "CD22", "CD38", "CD42a+CD61",
  "CD203c", "CD123", "CD4"
)])

hm2 <- Heatmap(hm2_data,
  col = colors,
  na_col = "white",
  rect_gp = gpar(col = "grey80", lwd = 0.5),
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  heatmap_legend_param = list(title = "Expression"),
  row_names_side = "left",
  height = unit(nrow(hm2_data) * size, "cm"),
  width = unit(ncol(hm2_data) * size, "cm"),
  layer_fun = function(j, i, x, y, width, height, fill) {
    v <- pindex(hm2_data, i, j)
    l <- !is.na(v)
    grid.text(v[l], x[l], y[l], gp = gpar(fontsize = 4))
  }
)


pdf(paste0(output_heatmaps, "all_features_all_patients.pdf"), width = 15, height = 23)
print(hm_sign + hm2 + hm1)
dev.off()

for (group in names(groupsList)) {
  namesOfGroups <- unlist(strsplit(group, "_"))
  nameOfGroup1 <- namesOfGroups[1]
  nameOfGroup2 <- namesOfGroups[3]
  sign_features <- rownames(stats[grepl("\\*", stats[, nameOfGroup1]), ])
  sign_features_data <- hm_data[sign_features, ]
  pl <- lapply(rownames(sign_features_data), function(x) {
    plotdf <- rbind(
      data.frame(MCL_pctg = sign_features_data[x, groupsList$NPM1_vs_noNPM1$group1], group = nameOfGroup1),
      data.frame(MCL_pctg = sign_features_data[x, groupsList$NPM1_vs_noNPM1$group2], group = nameOfGroup2)
    )
    p <- ggplot(data = plotdf, aes(x = group, y = MCL_pctg)) +
      geom_boxplot(outlier.shape = NULL) +
      geom_point(position = ggbeeswarm::position_quasirandom()) +
      theme_minimal() +
      ggtitle(x)
    return(p)
  })
  pdf(paste0(output_heatmaps, "boxplots_", group, ".pdf"), width = 5 * ceiling(sqrt(length(pl))), height = 5 * ceiling(sqrt(length(pl))))
  print(ggarrange(plotlist = pl))
  dev.off()
}

#---------------------------Heatmap all comparisons-----------------------------
# Heatmap
size <- 0.15
# --Filter rows
emptyCols <- grep("1622|1801|1872|1905|0019|1917|1905|1754|1604", colnames(hm_data))
hm_data <- filter_heatmap(hm_data, cutoff = 0.5, nPatient = 3)

# --Cluster column
hm_data_tmp <- hm_data[, -emptyCols]
hm_data_tmp <- hm_data_tmp[, hclust(dist(t(hm_data_tmp)))$order]
hm_data <- cbind(hm_data_tmp, hm_data[, emptyCols])

# --column annotation
annotation_labels <- gsub("_", " ", colnames(df_ha[colnames(hm_data), ]))
annotation_labels <- gsub("at diagnosis", " ", annotation_labels)
annotation_labels <- gsub("alive 2y after diagnosis", "survival status 2y", annotation_labels)
annotation_labels <- gsub("AML1.ETO  ", "AML1-ETO", annotation_labels)
annotation_labels <- gsub("inv16  ", "inv(16)", annotation_labels)
annotation_labels <- gsub("pml rara  ", "PML-RARA", annotation_labels)
annotation_labels <- gsub("MLL.AF9  ", "MLL-AF9", annotation_labels)
annotation_labels <- gsub("BCR.ABL  ", "BCR-ABL1", annotation_labels)
annotation_labels <- gsub("MLL.PTD", "MLL-PTD", annotation_labels)
annotation_labels <- gsub("CEBPalfa", "CEBPA", annotation_labels)
annotation_labels <- gsub("FLT3TKD", "FLT3-TKD", annotation_labels)
annotation_labels <- gsub("FLT3ITD", "FLT3-ITD", annotation_labels)
annotation_labels <- gsub("AML with recurrent genetic abnormalities", "genetic abnormalities", annotation_labels)
ha <- columnAnnotation(
  df = df_ha[colnames(hm_data), ],
  col = cols,
  gp = gpar(col = "grey80", lwd = 0.5),
  annotation_name_gp = gpar(fontsize = 5),
  annotation_legend_param = list(
    "NPM1" = list(title = "Genetic abnormalities"),
    "transplant" = list(title = "Outcome")
  ),
  annotation_label = annotation_labels,
  show_legend = c(TRUE, rep(FALSE, 22), rep(TRUE, 5)),
  simple_anno_size = unit(size, "cm"),
  na_col = "grey90"
)
# --Feature heatmap
column_split <- metadata[, c("PatientID", "AML_with_recurrent_genetic_abnormalities")]
column_split <- tibble::deframe(column_split)[colnames(hm_data)]
column_split[is.na(column_split)] <- "other"

rn_hm <- gsub("_", " ", rownames(hm_data))
MCs <- stringr::str_extract(rn_hm, "(?<=MC)[0-9]*")
rn_hm[nchar(MCs) == 1] <- gsub("MC([0-9])", "MC0\\1", rn_hm[nchar(MCs) == 1])
hm1_subset <- Heatmap(hm_data,
  cluster_rows = TRUE,
  bottom_annotation = ha,
  cluster_columns = FALSE,
  column_split = column_split,
  row_labels = rn_hm,
  rect_gp = gpar(col = "grey80", lwd = 0.5),
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  column_title_gp = gpar(fontsize = 6),
  heatmap_legend_param = list("title" = "Metacluster percentage"),
  na_col = "grey90",
  col = col_fun,
  height = unit(nrow(hm_data) * size, "cm"),
  width = unit(ncol(hm_data) * size, "cm")
)

# -- Significance heatmap
col_labels <- c(
  "relapse after chemotherapy", "relapse after transplant", "survival status 2y", "ELN",
  "NPM1", "inv(16)"
)
rn_hm_sign <- gsub("_", " ", rownames(stats[rownames(hm_data), ]))
MCs <- stringr::str_extract(rn_hm_sign, "(?<=MC)[0-9]*")
rn_hm_sign[nchar(MCs) == 1] <- gsub("MC([0-9])", "MC0\\1", rn_hm_sign[nchar(MCs) == 1])
hm_sign_subset <- Heatmap(as.matrix(stats[rownames(hm_data), ]),
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  col = c("*" = "#FFE5F5", "**" = "#FF99D8", "***" = "#FF0AA1", "ns" = "white"),
  na_col = "grey90",
  heatmap_legend_param = list(
    title = "Significance",
    at = c("*", "**", "***", "ns"),
    labels = c("p < 0.05", "p < 0.01", "p < 0.001", "not significant")
  ),
  rect_gp = gpar(col = "grey80", lwd = 0.5),
  row_names_side = "left",
  row_labels = rn_hm_sign,
  column_labels = col_labels,
  height = unit(nrow(stats[rownames(hm_data), ]) * size, "cm"),
  width = unit(ncol(stats[rownames(hm_data), ]) * size, "cm")
)

# --Expression heatmap
hm2_data <- hm2_data[rownames(hm_data), ]
colors <- c("-" = "white", "+" = "linen")
hm2_subset <- Heatmap(hm2_data,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  col = colors,
  na_col = "white",
  heatmap_legend_param = list("title" = "Expression"),
  rect_gp = gpar(col = "grey80", lwd = 0.5),
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  row_names_side = "left",
  height = unit(nrow(hm2_data) * size, "cm"),
  width = unit(ncol(hm2_data) * size, "cm"),
  layer_fun = function(j, i, x, y, width, height, fill) {
    v <- pindex(hm2_data, i, j)
    l <- !is.na(v)
    grid.text(v[l], x[l], y[l], gp = gpar(fontsize = 4))
  }
)


pdf(
  paste0(output_heatmaps, "all_comparisons_all_patients_filtered.pdf"),
  width = 20,
  height = 6
)
print(hm_sign_subset + hm2_subset + hm1_subset)
dev.off()

#-------------------------Volcanoplot + CM--------------------------------------
all_features <- readRDS(paste0(input_RDS, "all_features_all.RDS"))
rownames(all_features) <- all_features$PatientID
all_features_pctgs <- all_features[, grep("XMC", colnames(all_features))]
all_features_pctgs_filtered <- filter_heatmap(t(all_features_pctgs), cutoff = 0.5, nPatient = 3) %>% rownames()

groupsList <- groupsList[c(2, 1, 3, 4, 6, 5)]

plots <- lapply(names(groupsList), function(group) {
  if (group == "relapse_vs_noRelapse_noTransplant") {
    title <- "relapse after chemotherapy"
  } else if (group == "relapse_vs_noRelapse_transplant") {
    title <- "relapse after transplant"
  } else if (group == "A2_vs_D2") {
    title <- "survival status 2y"
  } else if (group == "favorable_vs_poor") {
    title <- "ELN"
  } else if (group == "NPM1_vs_noNPM1") {
    title <- "NPM1"
  } else if (group == "inv16_vs_noInv16") {
    title <- "inv(16)"
  }
  namesOfGroups <- unlist(strsplit(group, "_"))
  nameOfGroup1 <- namesOfGroups[1]
  nameOfGroup2 <- namesOfGroups[3]
  # Volcano plot
  stats <- GroupStats(all_features_pctgs, groups = groupsList[[group]]) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    mutate(sign = ifelse(`-log10 p values` > -log10(0.05), "p_sign", "not_sign")) %>%
    mutate(p_adj_sign = ifelse(`adjusted p values` < 0.05, "p_adj_sign", sign)) %>%
    mutate(diff_exp = ifelse(`medians group1` > `medians group2`, paste0(p_adj_sign, "_1>2"), paste0(p_adj_sign, "_1<2")))
  medianWithoutNA <- function(x) {
    median(x[which(!is.na(x))])
  }

  fix_values <- function(FCs, max) {
    FCs[is.nan(FCs)] <- 0
    FCs[is.infinite(FCs)] <- max * sign(FCs[is.infinite(FCs)])
    return(FCs)
  }
  medians_group1 <- apply(all_features_pctgs[rownames(all_features_pctgs) %in% groupsList[[group]]$group1, ], 2, medianWithoutNA)
  medians_group2 <- apply(all_features_pctgs[rownames(all_features_pctgs) %in% groupsList[[group]]$group2, ], 2, medianWithoutNA)
  FCs <- pmax(medians_group2, medians_group1) / pmin(medians_group2, medians_group1) * (-1)^sapply(1:ncol(all_features_pctgs), function(i) which.max(c(medians_group1[i], medians_group2[i])))
  logFCs <- log(abs(FCs)) * sign(FCs)
  max <- max(abs(FCs[is.finite(FCs)]), na.rm = FALSE)
  logmax <- max(abs(logFCs[is.finite(logFCs)]), na.rm = FALSE)
  FCs <- fix_values(FCs, max)
  logFCs <- fix_values(logFCs, logmax)

  stats$`fold changes` <- FCs
  stats$`log10 fold changes` <- logFCs
  stats$short_rowname <- gsub("Tube_|X", "", stats$rowname)
  stats_text <- stats %>%
    filter(rowname %in% all_features_pctgs_filtered &
      `p values` <= 0.05)
  return(stats)
  stats$in_heatmap <- stats$rowname %in% stats_text$rowname
  vp <- ggplot(stats, aes(x = -`log10 fold changes`, y = `-log10 p values`, col = diff_exp)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.3) +
    ggrepel::geom_text_repel(data = stats_text, aes(label = short_rowname), size = 2, max.overlaps = 40, segment.size = 0.1) +
    geom_point(aes(shape = in_heatmap), size = 0.5) +
    scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1)) +
    scale_color_manual(values = c(
      "p_sign_1<2" = "#FF9B71", "p_sign_1>2" = "#61E294",
      "p_adj_sign_1<2" = "firebrick", "p_adj_sign_1>2" = "#054D48",
      "not_sign_1>2" = "black", "not_sign_1<2" = "black"
    )) +
    ylim(c(0, 6)) +
    theme_classic() +
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
  return(vp)
})

pdf("./blasts/Volcano_plots/volcano_plot.pdf", width = 6.88, height = 4)
ggarrange(plotlist = plots)
dev.off()
#---------------------------Relapse survival analysis---------------------------
all_features <- readRDS(paste0(input_RDS, "all_features_all.RDS"))
all_features_MCpctgs <- all_features[, grepl("XMC", colnames(all_features))]
all_features_MCpctgs_median <- as.data.frame(apply(all_features_MCpctgs, 2, function(x) {
  med_x <- median(x, na.rm = T)
  return(ifelse(x > med_x, "above_equal_median", "below_median"))
}))
rownames(all_features_MCpctgs_median) <- all_features_MCpctgs_median$PatientID <- rownames(all_features_MCpctgs) <- all_features$PatientID

data <- metadata %>%
  select(
    PatientID, date_of_diagnosis, cause_of_death,
    date_of_relapse, relapse_.morphological.,
    age_at_diagnosis, transplant, Date_of_last_contact,
    date_of_transplantation
  ) %>%
  filter(!is.na(relapse_.morphological.)) %>%
  mutate(transplant_before_relapse = ifelse(date_of_transplantation < date_of_relapse,
    "yes", "no"
  )) %>%
  mutate(transplant_before_relapse = ifelse(relapse_.morphological. == "no" & transplant == "yes",
    "yes",
    transplant_before_relapse
  )) %>%
  mutate(transplant_before_relapse = ifelse(relapse_.morphological. == "primary refractory" & transplant == "yes",
    "no",
    transplant_before_relapse
  )) %>%
  mutate(transplant_before_relapse = ifelse(is.na(transplant_before_relapse), "no", transplant_before_relapse)) %>%
  mutate(time_to_relapse = as.numeric(date_of_relapse - date_of_diagnosis)) %>%
  mutate(relapse_status = case_when(
    relapse_.morphological. == "no" ~ 0,
    grepl("relapse|primary", relapse_.morphological.) ~ 1
  )) %>%
  mutate(time_to_relapse = ifelse(is.na(time_to_relapse) & relapse_status == 0,
    Date_of_last_contact - date_of_diagnosis,
    time_to_relapse
  )) %>%
  mutate(time_to_relapse = ifelse(relapse_.morphological. == "primary refractory", 0, time_to_relapse)) %>%
  select(PatientID, relapse_status, time_to_relapse, age_at_diagnosis, transplant) %>%
  full_join(all_features_MCpctgs_median, by = "PatientID") %>%
  tibble::column_to_rownames("PatientID") %>%
  select(-age_at_diagnosis)


results <- multiUnivariateLogRank(data,
  time = "time_to_relapse",
  event = "relapse_status",
  covariates = "transplant"
)

sign_features <- results$summary %>%
  filter(p_value < 0.05) %>%
  rownames()

pl <- list()
for (t in seq(0, 6)) {
  tube_sign_features <- sign_features[grepl(paste0("Tube_", t), sign_features)]
  if (length(tube_sign_features) != 0) {
    fsom <- readRDS(paste0(input_RDS, "fsom_agg_toMapFiles_all_fullNames_tube_", t, ".RDS"))
    channelpairs <- plots[[as.character(t)]]
    for (feature in tube_sign_features) {
      tube <- stringr::str_extract(feature, "Tube_[0-9]")
      MC <- as.numeric(stringr::str_extract(feature, "(?<=MC)[0-9]*"))
      phenotype <- gsub(".*: ", "", levels(fsom$metaclustering)[MC])
      phenotype <- gsub("FSC-A\\+, |SSC-A\\+, |, CD45\\+", "", phenotype)
      f <- as.formula(paste("survival::Surv(survivalTime, survivalStatus) ~ transplant + ", feature))
      df <- data.frame(
        feature = data[, feature], survivalTime = data[, "time_to_relapse"],
        survivalStatus = data[, "relapse_status"], transplant = data[, "transplant"]
      )
      df <- na.omit(df)
      colnames(df)[1] <- feature

      # Survival plot
      fit <- survminer::surv_fit(f, data = df)
      p <- survminer::ggsurvplot(fit,
        data = df,
        palette = c("dodgerblue4", "dodgerblue", "orange", "gold"),
        pval = T,
        conf.int = F,
        fun = "event",
        break.time.by = 365,
        risk.table = T,
        title = paste0(feature, ": ", fsom$metaclustering[MC]),
        legend = "right"
      )
      order_at_2y <- p[["data.survplot"]] %>%
        mutate(time_close_to_2y = abs(time - 365 * 2)) %>%
        group_by(strata) %>%
        arrange(time_close_to_2y) %>%
        dplyr::slice(1) %>%
        pull(surv) %>%
        order()
      pl[[feature]] <- list(plot = p, order = paste0(order_at_2y, collapse = ""))
      p_legend <- as_ggplot(ggpubr::get_legend(p$plot))
      p <- survminer::ggsurvplot(fit,
        data = df,
        palette = c("dodgerblue4", "dodgerblue", "orange", "gold"),
        pval = F,
        conf.int = F,
        fun = "event",
        break.time.by = 365,
        risk.table = T,
        legend = "none",
        title = paste0("p: ", results$summary[feature, "p_value"])
      )

      # Feature value vs survival time - ELN + genetics
      df1 <- data[, c("relapse_status", "time_to_relapse", "transplant")] %>%
        tibble::rownames_to_column(var = "PatientID") %>%
        left_join(metadata[, c("PatientID", "ELN_risk_classification_at_diagnosis")], by = "PatientID") %>%
        left_join(subgroups_AML, by = "PatientID")
      cols <- c(
        "BCR-ABL1" = "#DB324D", "inv(16)" = "#68C3D4", "inv(3)" = "#F5D547",
        "NPM1_A" = "#E3C9F5", "NPM1_B" = "#BF96DC",
        "NPM1_C" = "#9C64C3", "NPM1_D" = "#7932AA", "NPM1_E" = "#560091",
        "mutated RUNX1" = "#41E454", "t(15;17)" = "#E08E45",
        "t(8;21)" = "deeppink", "other_A" = "#C7E0C5",
        "other_B" = "#86AB83", "other_C" = "#457641", "other_D" = "#054200"
      )
      df2 <- all_features_MCpctgs[, feature, drop = F] %>% tibble::rownames_to_column(var = "PatientID")
      plotdf2 <- full_join(df1, df2, by = "PatientID")
      colnames(plotdf2) <- c("PatientID", "relapse_status", "time", "transplant", "ELN", "WHO_genetic_abn", "feature")
      plotdf2$relapse_status <- factor(plotdf2$relapse_status, levels = c(0, 1), labels = c("no_relapse", "relapse"))
      plotdf2$WHO_genetic_abn <- factor(plotdf2$WHO_genetic_abn, levels = names(cols))
      med <- median(plotdf2$feature, na.rm = T)
      p1 <- ggplot(plotdf2) +
        geom_point(aes(x = time, y = log(feature), col = WHO_genetic_abn, shape = ELN), size = 4) +
        scale_color_manual(values = cols) +
        geom_hline(yintercept = log(med), col = "red") +
        theme_minimal() +
        guides(col = guide_legend(ncol = 2))
      p1_legend <- as_ggplot(get_legend(p1))
      p1 <- p1 + theme(legend.position = "none")

      # Feature value vs survival time
      plotdf2 <- data %>%
        select(all_of(feature)) %>%
        tibble::rownames_to_column("PatientID") %>%
        right_join(plotdf2, by = "PatientID")
      colnames(plotdf2) <- c(
        "PatientID", "feature_median", "relapse_status", "time",
        "transplant", "ELN", "subgroup", "feature"
      )
      cols2 <- c(
        "transplant.above_equal_median" = "orange",
        "no transplant.below_median" = "dodgerblue",
        "transplant.below_median" = "gold",
        "no transplant.above_equal_median" = "dodgerblue4"
      )
      plotdf2$`Feature/transplant` <- interaction(
        ifelse(plotdf2$transplant == "yes", "transplant", "no transplant"),
        plotdf2$feature_median
      )
      p2 <- ggplot(plotdf2) +
        geom_point(aes(
          x = time,
          y = log(feature),
          col = `Feature/transplant`,
          shape = relapse_status
        ), size = 4) +
        geom_hline(yintercept = log(med), col = "red") +
        scale_color_manual(values = cols2) +
        scale_shape_manual(values = c(1, 16), labels = c("no", "yes"), na.translate = F) +
        theme_minimal() +
        guides(col = "none")
      p2_legend <- as_ggplot(get_legend(p2))
      p2 <- p2 + theme(legend.position = "none")
      left_up <- ggarrange(p1, p2, p2_legend, ncol = 1, nrow = 3, heights = c(5, 5, 1))

      # 2DScatter
      max_channel <- max(fsom$data[, GetChannels(fsom, unlist(channelpairs))])
      p2Dscatters <- Plot2DScatters(fsom, channelpairs,
        metaclusters = MC, plotFile = NULL,
        xLim = c(0, max_channel), yLim = c(0, max_channel), centers = F
      )
      p3 <- ggarrange(plotlist = lapply(p2Dscatters, function(x) x + coord_fixed() + ggtitle("")), nrow = 2, ncol = 3)
      right_up <- ggarrange(ggarrange(p1_legend, NULL, ncol = 1, nrow = 2, heights = c(3, 1)), p3, ncol = 2, nrow = 1, widths = c(1, 4))

      # Arrange plot
      upper <- ggarrange(left_up, ggarrange(right_up, p$plot, ncol = 1, nrow = 2, heights = c(3, 2)), ncol = 2, nrow = 1, widths = c(1, 3))
      total_plot <- ggarrange(upper, p$table, ncol = 1, nrow = 2, heights = c(5, 1))
      total_plot <- ggpubr::annotate_figure(total_plot, top = paste0(feature, ": ", phenotype))
      pdf(paste0("./blasts/Survival_analysis/relapse/survival_time_relapse_cr_", feature, ".pdf"), width = 15, height = 12)
      print(total_plot)
      dev.off()
    }
  }
}


pdf("./blasts/Survival_analysis/relapse/Survival_time_relapse.pdf", width = 23, height = 8)
orders <- order(sapply(pl, function(x) x$order))
pls <- lapply(pl, function(x) x$plot)
print(pls[orders])
dev.off()

#----------------------Survival analysis w competing risks----------------------
all_features <- readRDS(paste0(input_RDS, "all_features_all.RDS"))
rownames(all_features) <- all_features$PatientID
all_features_MCpctgs <- all_features[, grepl("XMC", colnames(all_features))]
all_features_MCpctgs_median <- as.data.frame(apply(all_features_MCpctgs, 2, function(x) {
  med_x <- median(x, na.rm = T)
  return(ifelse(x > med_x, "above_equal_median", "below_median"))
}))
rownames(all_features_MCpctgs_median) <- all_features_MCpctgs_median$PatientID <- all_features$PatientID
relapse_HSCT <- list(
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
  )
)
relapse_df <- rbind(
  data.frame(PatientID = relapse_HSCT$relapse_vs_noRelapse_noTransplant$group1, relapse_HSCT = "relapse after chemo"),
  data.frame(PatientID = c(relapse_HSCT$relapse_vs_noRelapse_noTransplant$group2, relapse_HSCT$relapse_vs_noRelapse_transplant$group2), relapse_HSCT = "no relapse"),
  data.frame(PatientID = relapse_HSCT$relapse_vs_noRelapse_transplant$group1, relapse_HSCT = "relapse after transplant"),
  data.frame(PatientID = metadata %>% filter(relapse_.morphological. == "primary refractory") %>% pull(PatientID), relapse_HSCT = "primary refractory")
)


data <- metadata %>%
  select(
    PatientID, date_of_diagnosis, cause_of_death,
    date_of_death, survival_status_at_last_contact,
    Date_of_last_contact, transplant
  ) %>%
  full_join(relapse_df, by = "PatientID") %>%
  mutate(survival_time_days = as.numeric(date_of_death - date_of_diagnosis)) %>%
  mutate(survival_status = case_when(
    survival_status_at_last_contact == "alive" ~ 0,
    cause_of_death %in% c(
      "death from AML progressive disease",
      "death in aplasie",
      "death with (R/R) AML but from another cause"
    ) ~ 1,
    !cause_of_death %in% c(
      "death from AML progressive disease",
      "death in aplasie",
      "death with (R/R) AML but from another cause"
    ) ~ 2
  )) %>%
  mutate(survival_time_days = ifelse(is.na(survival_time_days) & survival_status == 0,
    Date_of_last_contact - date_of_diagnosis,
    survival_time_days
  )) %>%
  mutate(survival_status = factor(survival_status)) %>%
  select(PatientID, survival_status, survival_time_days, relapse_HSCT, transplant) %>%
  full_join(all_features_MCpctgs_median, by = "PatientID") %>%
  tibble::column_to_rownames("PatientID") %>%
  mutate(survival_status = factor(survival_status)) %>%
  filter(!is.na(survival_status))


results <- multiUnivariateCI(data %>% select(-relapse_HSCT),
  time = "survival_time_days",
  event = "survival_status"
)
sign_features <- results$significant
pl <- lapply(sign_features, function(feature) {
  f <- as.formula(paste("tidycmprsk::Surv(survivalTime, survivalStatus) ~ ", feature))
  df <- data.frame(
    feature = data[, feature], survivalTime = data[, "survival_time_days"],
    survivalStatus = data[, "survival_status"]
  )
  colnames(df)[1] <- feature
  cuminc(f, data = df) %>%
    ggcuminc() +
    labs(x = "Days") +
    add_confidence_interval() +
    add_risktable() +
    ggtitle(paste0(feature, " (p: ", results$summary[feature, "p_value1"], ")"))
})

pdf("./blasts/Survival_analysis/cause_of_death/Survival_time_cod_cr.pdf", width = 5 * 6, height = 4 * 4)
ggarrange(plotlist = pl)
dev.off()


for (t in seq(0, 6)) {
  tube_sign_features <- sign_features[grepl(paste0("Tube_", t), sign_features)]
  if (length(tube_sign_features) != 0) {
    fsom <- readRDS(paste0(input_RDS, "fsom_agg_toMapFiles_all_fullNames_tube_", t, ".RDS"))
    channelpairs <- plots[[as.character(t)]]
    for (f in tube_sign_features) {
      # Median feature vs survival time
      tube <- stringr::str_extract(f, "Tube_[0-9]")
      MC <- as.numeric(stringr::str_extract(f, "(?<=MC)[0-9]*"))
      phenotype <- gsub(".*: ", "", levels(fsom$metaclustering)[MC])
      phenotype <- gsub("FSC-A\\+, |SSC-A\\+, |, CD45\\+", "", phenotype)

      # Feature value vs survival time, subgroups/ELN risk
      df1 <- data[, c("survival_status", "survival_time_days", "relapse_HSCT", "transplant")] %>%
        tibble::rownames_to_column("PatientID") %>%
        left_join(metadata[, c("PatientID", "ELN_risk_classification_at_diagnosis")], by = "PatientID") %>%
        left_join(subgroups_AML, by = "PatientID")
      cols <- c(
        "BCR-ABL1" = "#DB324D", "inv(16)" = "#68C3D4", "inv(3)" = "#F5D547",
        "NPM1_A" = "#E3C9F5", "NPM1_B" = "#BF96DC",
        "NPM1_C" = "#9C64C3", "NPM1_D" = "#7932AA", "NPM1_E" = "#560091",
        "mutated RUNX1" = "#41E454", "t(15;17)" = "#E08E45",
        "t(8;21)" = "deeppink", "other_A" = "#C7E0C5",
        "other_B" = "#86AB83", "other_C" = "#457641", "other_D" = "#054200"
      )
      df2 <- all_features_MCpctgs[, f, drop = F] %>% tibble::rownames_to_column("PatientID")
      plotdf1 <- full_join(df1, df2, by = "PatientID")
      colnames(plotdf1) <- c("PatientID", "survival_status", "time", "relapse", "transplant", "ELN", "WHO_genetic_abn", "feature")
      plotdf1$WHO_genetic_abn <- factor(plotdf1$WHO_genetic_abn, levels = names(cols))
      med <- median(plotdf1$feature, na.rm = T)
      p1 <- ggplot(plotdf1) +
        geom_point(aes(x = time, y = log(feature), shape = ELN, col = WHO_genetic_abn), size = 3) +
        geom_hline(yintercept = log(med), col = "red") +
        scale_color_manual(values = cols) +
        guides(col = guide_legend(ncol = 2)) +
        theme_minimal()
      p1_legend <- as_ggplot(ggpubr::get_legend(p1))
      p1 <- p1 + theme(legend.position = "none")

      # Feature value vs survival time, relapse/transplant
      plotdf2 <- data %>%
        select(all_of(f)) %>%
        tibble::rownames_to_column("PatientID") %>%
        right_join(plotdf1, by = "PatientID")
      colnames(plotdf2) <- c(
        "PatientID", "feature_median", "survival_status", "time", "relapse_status", "transplant",
        "ELN", "WHO_genetic_abn", "feature"
      )
      cols2 <- c(
        "transplant.above_equal_median" = "orange",
        "no transplant.below_median" = "dodgerblue",
        "transplant.below_median" = "gold",
        "no transplant.above_equal_median" = "dodgerblue4"
      )
      plotdf2$`Feature/transplant` <- interaction(
        ifelse(plotdf2$transplant == "yes", "transplant", "no transplant"),
        plotdf2$feature_median
      )
      p2 <- ggplot(plotdf2) +
        geom_point(aes(
          x = time,
          y = log(feature),
          shape = survival_status,
          col = feature_median
        ), size = 4) +
        geom_hline(yintercept = log(med), col = "red") +
        scale_color_manual(values = c("#467388", "#EB7779")) +
        scale_shape_manual(values = c(1, 16, 13), labels = c("alive", "death_by_AML", "death_by_other_cause")) +
        theme_minimal() +
        guides(col = "none")
      p2_legend <- as_ggplot(ggpubr::get_legend(p2))
      p2 <- p2 + theme(legend.position = "none")

      upper_left <- ggarrange(p1, p2, p2_legend, ncol = 1, nrow = 3, heights = c(5, 5, 1))

      # 2DScatter
      max_channel <- max(fsom$data[, GetChannels(fsom, unlist(channelpairs))])
      p2Dscatters <- Plot2DScatters(fsom, channelpairs,
        metaclusters = MC, plotFile = NULL,
        xLim = c(0, max_channel), yLim = c(0, max_channel), centers = F
      )
      p3 <- ggarrange(plotlist = lapply(p2Dscatters, function(x) x + coord_fixed() + ggtitle("")), nrow = 2, ncol = 3)
      upper_right <- ggarrange(p1_legend, p3, nrow = 1, ncol = 2, widths = c(1, 4))

      # Cumulative incidence
      plotdf4 <- data.frame(
        feature = data[, f], survivalTime = data[, "survival_time_days"],
        survivalStatus = data[, "survival_status"]
      )
      colnames(plotdf4)[1] <- f
      formula <- as.formula(paste("tidycmprsk::Surv(survivalTime, survivalStatus) ~ ", f))
      p4 <- cuminc(formula, data = plotdf4) %>%
        ggcuminc(theme = ggplot2::theme_classic()) +
        add_censor_mark() +
        labs(x = "Days") +
        scale_color_manual(values = c("#467388", "#EB7779")) +
        add_confidence_interval() +
        scale_fill_manual(values = c("#467388", "#EB7779")) +
        add_risktable() +
        ggtitle(paste0("p: ", results$summary[f, "p_value1"]))

      # Parse figures
      total_plot <- ggarrange(upper_left, ggarrange(upper_right, p4, nrow = 2, ncol = 1), ncol = 2, nrow = 1, widths = c(1, 3))
      total_plot <- ggpubr::annotate_figure(total_plot, top = paste0(f, ": ", phenotype))
      pdf(paste0("./blasts/Survival_analysis/cause_of_death/survival_time_cod_cr_", f, ".pdf"), width = 15, height = 12)
      print(total_plot)
      dev.off()
    }
  }
}

#-----------------------------------Patient Overview----------------------------
# Parse metadata
metadata_lim_tmp <- metadata[, c(
  "PatientID", "gender",
  "ELN_risk_classification_at_diagnosis",
  "age_at_diagnosis"
)]

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

order_ELN <- metadata$PatientID[order(
  metadata$ELN_risk_classification_at_diagnosis,
  metadata$Date_of_last_contact - metadata$date_of_diagnosis
)]

# Prepare plot
metadata_plot <- metadata %>%
  select(
    PatientID, transplant, `relapse_.morphological.`, date_of_diagnosis,
    date_of_death, date_of_relapse, cause_of_death, Date_of_last_contact,
    survival_status_at_last_contact, ELN_risk_classification_at_diagnosis,
    date_of_transplantation, gender, age_at_diagnosis
  ) %>%
  mutate(lifetime = Date_of_last_contact - date_of_diagnosis) %>%
  mutate(relapse_time = date_of_relapse - date_of_diagnosis) %>%
  mutate(transplant_time = date_of_transplantation - date_of_diagnosis) %>%
  mutate(death_from_AML = cause_of_death %in% c(
    "death from AML progressive disease",
    "death in aplasie",
    "death with (R/R) AML but from another cause"
  )) %>%
  mutate(relapse = grepl("relapse", `relapse_.morphological.`)) %>%
  mutate(PatientID = factor(PatientID, levels = order_ELN))

sub_metadata_gender <- metadata[, c("PatientID", "gender")] %>% mutate(x = -60)
sub_metadata_age <- metadata[, c("PatientID", "age_at_diagnosis")] %>% mutate(x = -130)
sub_metadata_ELN <- metadata[, c("PatientID", "ELN_risk_classification_at_diagnosis")] %>% mutate(x = -200)


ticks <- c(
  -200, -130, -60,
  seq(0, as.numeric(max(metadata_plot$lifetime)), 365)
)
labels <- c(
  "ELN", "age", "gender",
  "0y", "1y", "2y", "3y", "4y", "5y"
)

metadata_plot$survival_status_AML <- "alive"
metadata_plot$survival_status_AML[metadata_plot$survival_status_at_last_contact == "died" & metadata_plot$death_from_AML] <- "death from AML"
metadata_plot$survival_status_AML[metadata_plot$survival_status_at_last_contact == "died" & !metadata_plot$death_from_AML] <- "death from other causes"
p <- ggplot(metadata_plot) +
  # Life time
  geom_segment(aes(x = 0, xend = lifetime, y = PatientID, yend = PatientID, col = survival_status_AML), lineend = "round", linewidth = 0.3) +
  scale_color_manual(values = c("alive" = "green4", "death from AML" = "orangered3", "death from other causes" = "grey70")) +
  guides(color = guide_legend(override.aes = list(linewidth = 1), title = "Survival status")) +
  geom_vline(xintercept = 730, color = "grey50", linewidth = 0.3, linetype = "dashed") +
  # Cause of dead
  ggnewscale::new_scale_color() +
  geom_point(
    data = metadata_plot %>% filter(grepl("death", survival_status_AML)),
    aes(x = lifetime, y = PatientID, shape = survival_status_AML, color = survival_status_AML), size = 0.6, stroke = 0.2
  ) +
  scale_shape_manual(values = c("death from other causes" = 13, "death from AML" = 16)) +
  scale_color_manual(values = c("death from other causes" = "grey70", "death from AML" = "orangered3")) +
  # Relapse
  geom_text(aes(x = relapse_time, y = PatientID), label = "R", size = 1.5, fontface = "bold") +
  # Transplant
  geom_text(aes(x = transplant_time, y = PatientID), label = "T", size = 1.5, fontface = "bold") +
  # Metadata
  ggnewscale::new_scale_color() +
  geom_point(data = sub_metadata_ELN, aes(x = x, y = PatientID, col = ELN_risk_classification_at_diagnosis), shape = 15, size = 1) +
  scale_color_manual(values = c("poor" = "firebrick", "intermediate" = "orange", "favorable" = "lightgreen")) +
  guides(color = guide_legend(override.aes = list(size = 2), title = "ELN risk")) +
  ggnewscale::new_scale_color() +
  geom_point(data = sub_metadata_gender, aes(x = x, y = PatientID, col = gender), shape = 15, size = 1) +
  scale_color_manual(values = c("male" = "#8372e0", "female" = "#FFF170")) +
  guides(color = guide_legend(override.aes = list(size = 2), title = "Gender")) +
  ggnewscale::new_scale_color() +
  geom_point(data = sub_metadata_age, aes(x = x, y = PatientID, col = age_at_diagnosis), shape = 15, size = 1) +
  guides(color = guide_colourbar(title = "Age at diagnosis")) +
  scale_color_viridis_c(direction = -1) +
  # Theme
  geom_vline(xintercept = 0, linewidth = 0.3) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
    axis.text.y = element_text(size = 4),
    axis.line.x = element_line(linewidth = 0.3),
    axis.line.y = element_blank(),
    axis.ticks = element_line(linewidth = 0.3)
  ) +
  scale_x_continuous(breaks = ticks, labels = labels) +
  xlab("Time from diagnosis")

p_legend <- as_ggplot(get_legend(p))

pdf(paste0("./blasts/overview_patients.pdf"), width = 3.34, height = 8)
print(p + theme(legend.position = "none"))
dev.off()
pdf(paste0("./blasts/overview_patients_legend.pdf"), width = 5, height = 7)
print(p_legend)
dev.off()
#------------------------------Balanced Accuracy plot---------------------------
order_plot <- c(
  "favorable vs poor", "inv16 vs noInv16", "NPM1 vs noNPM1", "relapse vs noRelapse transplant", "relapse vs noRelapse noTransplant",
  "A2 vs D2"
)
labels_plot <- c("ELN risk", "inv(16)", "NPM1", "relapse after transplant", "relapse after chemotherapy", "survival status 2y")
bestResults <- do.call(rbind, lapply(names(groupsList), function(x) {
  listOfResults <- readRDS(paste0("./blasts/XGBoost/RDS/QualityMetrics_", x, ".RDS")) %>%
    filter(name == "balanced_accuracies") %>%
    arrange(desc(value)) %>%
    dplyr::slice(c(1, which(featureSet == "only_metadata")))
  listOfResults$c_group <- x
  return(listOfResults)
})) %>%
  as.data.frame() %>%
  mutate(Features = factor(ifelse(featureSet == "only_metadata", "metadata", "metadata + phenotype"),
    levels = c("metadata", "metadata + phenotype")
  )) %>%
  mutate(c_group = gsub("_", " ", c_group)) %>%
  mutate(c_group = factor(c_group, levels = order_plot))
allResults <- do.call(rbind, lapply(names(groupsList), function(x) {
  listOfResults <- readRDS(paste0("./blasts/XGBoost/RDS/QualityMetrics_", x, ".RDS")) %>%
    filter(name == "balanced_accuracies") %>%
    filter(featureSet == "all")
  listOfResults$c_group <- x
  return(listOfResults)
}))
p <- ggplot(bestResults) +
  geom_col(aes(x = c_group, y = value, fill = Features, group = Features),
    position = position_dodge2()
  ) +
  scale_fill_manual(values = c("metadata" = "grey80", "metadata + phenotype" = "#cdeda6ff")) +
  scale_x_discrete(labels = labels_plot) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey40", linewidth = 0.3) +
  xlab("") +
  ylab("Balanced accuracy") +
  ylim(0, 1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, size = 4),
    axis.text.y = element_text(size = 4),
    axis.line.x = element_line(linewidth = 0.3),
    axis.line.y = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.title = element_text(size = 5)
  )

pdf(paste0("./blasts/XGBoost/BA.pdf"), width = 3.34, height = 2)
print(p)
dev.off()

#-----------------------------Patient mapping-----------------------------------
fsom <- readRDS(paste0(output_heatmaps, "fsom_agg_toMapFiles_all_fullNames_tube_0.RDS"))
ff <- read.FCS(samples_norm %>% filter(PatientID == "1627") %>% filter(Tube == 0) %>% pull(File))
nd <- NewData(fsom, ff)
p1 <- PlotStars(fsom, backgroundValues = fsom$metaclustering, list_insteadof_ggarrange = TRUE, maxNodeSize = 2)
p2 <- PlotStars(nd, backgroundValues = nd$metaclustering, list_insteadof_ggarrange = TRUE, maxNodeSize = 2)
a_p <- ggarrange(NULL, p1$starLegend, ggarrange(p1$tree, p2$tree, nrow = 1, ncol = 2), ncol = 1, nrow = 3, heights = c(0.5, 1, 9))

pdf("D:/PhD/Papers/PaperEuroflow/FlowSOM_mapping.pdf", width = 20, height = 15)
print(a_p)
dev.off()
