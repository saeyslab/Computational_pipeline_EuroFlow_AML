############################################################################### /
# LIBRARIES
############################################################################### /
library(flowCore)
library(FlowSOM)
library(Rtsne)
library(dplyr)
library(ggplot2)
library(flowDensity)
source("./3_helperFunctions.R")

############################################################################### /
# SCRIPT
############################################################################### /
#----------------------------------Define variables-----------------------------
seed <- 2022
set.seed(seed)
mc.cores <- 1

input_raw_FCS_files <- "./blasts/raw_FCS_files_blasts/"
input_compensation_matrices <- "./comp_matr/"
input_metadata <- "./Excel_metadata/"

output_unif_FCS_files <- "./blasts/FCS_files_blasts/"
output_prepro_FCS_files <- "./blasts/Preprocessed_fcs/"
output_norm_FCS_files <- "./blasts/Normalized_fcs/"
output_RDS <- "./blasts/RDS/"

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

#--------------------------------Parse raw data---------------------------------
fcs_files <- paste0("./blasts/raw_FCS_files_blasts/", list.files(input_raw_FCS_files, pattern = ".*fcs$"))
noSpaces_fcs_files <- gsub(" ", "_", fcs_files)
file.rename(fcs_files, noSpaces_fcs_files)
samples <- data.frame(
  File = noSpaces_fcs_files,
  Panel = stringr::str_match(noSpaces_fcs_files, "ALOT|AML")[, 1],
  PatientID = stringr::str_match(noSpaces_fcs_files, "([0-9]{4}-?[0-9]?)")[, 2],
  Tube = as.numeric(stringr::str_match(noSpaces_fcs_files, "(?<=Tube|AML_?)[1-7]")),
  stringsAsFactors = FALSE,
  PBOorNOT = grepl("PBO", noSpaces_fcs_files)
)

samples$Tube[is.na(samples$Tube)] <- 0

#-----------------------------Uniformize data-----------------------------------
channels_of_interest <- c(
  "FSC-A", "FSC-H", "SSC-A", "FITC-A", "PE-A",
  "PerCP-Cy5-5-A", "PE-Cy7-A", "APC-A", "APC-H7-A", "HV450-A",
  "HV500-A"
)

markers_of_interest_all <- list(
  "0" = c(
    "FSC-A", "FSC-H", "SSC-A", "cyMPO",
    "cyCD79a", "CD34", "CD19", "CD7", "CD3",
    "cyCD3", "CD45"
  ),
  "1" = c(
    "FSC-A", "FSC-H", "SSC-A", "CD16", "CD13",
    "CD34", "CD117", "CD11b", "CD10", "HLA-DR",
    "CD45"
  ),
  "2" = c(
    "FSC-A", "FSC-H", "SSC-A", "CD35", "CD64",
    "CD34", "CD117", "CD300e", "CD14", "HLA-DR",
    "CD45"
  ),
  "3" = c(
    "FSC-A", "FSC-H", "SSC-A", "CD36", "CD105",
    "CD34", "CD117", "CD33", "CD71", "HLA-DR",
    "CD45"
  ),
  "4" = c(
    "FSC-A", "FSC-H", "SSC-A", "nuTdT",
    "CD56", "CD34", "CD117", "CD7", "CD19",
    "HLA-DR", "CD45"
  ),
  "5" = c(
    "FSC-A", "FSC-H", "SSC-A", "CD15", "NG2",
    "CD34", "CD117", "CD22", "CD38", "HLA-DR",
    "CD45"
  ),
  "6" = c(
    "FSC-A", "FSC-H", "SSC-A", "CD42a+CD61",
    "CD203c", "CD34", "CD117", "CD123", "CD4",
    "HLA-DR", "CD45"
  ),
  "7" = c(
    "FSC-A", "FSC-H", "SSC-A", "CD41", "CD25",
    "CD34", "CD117", "CD42b", "CD9", "HLA-DR",
    "CD45"
  )
)
listTubes <- list(
  "0" = "_ALOT", "1" = "_AML1", "2" = "_AML2", "3" = "_AML3",
  "4" = "_AML4", "5" = "_AML5", "6" = "_AML6", "7" = "_AML7"
)

y <- 1
stopNow <- 0
wrongTubeFiles <- list()
for (tube in 0:7) {
  samples_tube <- samples %>% filter(Tube == tube)
  markers_o_i <- markers_of_interest_all[[as.character(tube)]]
  markersInTube <- tolower(sortMarkers(markers_o_i[4:11]))
  for (file in samples_tube$File) {
    tube_correct <- tube
    print(paste0("Tube ", tube, ": ", y, " of ", nrow(samples)))
    oldFile <- file
    ff <- read.FCS(file)
    tmp <- samples[samples$File == file, ]
    names_ff <- unname(ff@parameters@data[["name"]])
    desc_ff <- unname(ff@parameters@data[["desc"]])
    colnames(ff) <- gsub(" ", "-", colnames(ff))

    # PacB and PacO are changed to HV450 and HV500 respectively
    colnames(ff) <- gsub("PacB", "HV450", colnames(ff))
    colnames(ff) <- gsub("PacO", "HV500", colnames(ff))
    relapse_or_not <- grepl("_R", file)
    PBOorNOT <- grepl("PB", file)

    ff <- ff[, channels_of_interest]

    if (grepl("0020|0024|0025", file) & tube == 0) {
      ff@parameters@data[["desc"]][c(10, 11)] <- c("cyCD3", "CD45")
    }
    if (grepl("1608", file) & tube == 2) {
      ff@parameters@data[["desc"]][5] <- "CD64"
    }
    if (grepl("1604", file) & tube == 4) {
      ff@parameters@data[["desc"]][4] <- "nuTdT"
    }

    markersInFF <- ff@parameters@data[["desc"]]

    # Remove -, + and spaces in the marker names  and make marker names uniform
    markersInFF <- gsub(" |-|+", "", markersInFF)
    markersInFF <- gsub("HV500C|HV450", "", markersInFF) # 1747 ALOT
    markersInFF <- gsub("cyCD79a", "cyCD79a", markersInFF) # 1747 ALOT
    markersInFF <- gsub("IREM2", "CD300e", markersInFF)
    markersInFF <- gsub("Cy 79a|Cy79a", "cyCD79a", markersInFF)
    markersInFF <- gsub("CyCD3", "cyCD3", markersInFF)
    markersInFF <- gsub("HLADR", "HLA-DR", markersInFF)

    ff@parameters@data[["desc"]] <- markersInFF

    # Check if ff markers resemble markers from tube in title and if not change it
    markersInFF <- markersInFF[!is.na(markersInFF)]
    markersInFF <- tolower(sortMarkers(markersInFF))

    file <- paste0(
      tmp$PatientID, "_", tmp$Panel, ifelse(tube == 0, "", tmp$Tube),
      ifelse(PBOorNOT, "_PBO", ""), ifelse(relapse_or_not, "_R", ""), ".fcs"
    )

    if (!all(markersInTube %in% markersInFF)) {
      for (tube_markers in names(markers_of_interest_all)) {
        print(tube_markers)
        markers_o_i <- markers_of_interest_all[[tube_markers]]
        markersInTube_tmp <- tolower(sortMarkers(markers_o_i[4:11]))

        # Save wrong tubes in a file
        if (all(markersInFF %in% markersInTube_tmp)) {
          if (!file.exists(paste0(output_unif_FCS_files, "wrongTubes.txt"))) {
            write("Files in wrong tubes \n------------------------\n",
              file = paste0(output_unif_FCS_files, "wrongTubes.txt"), append = TRUE
            )
            write("Old filename\tNew Filename\tOld tube\tNew tube",
              file = paste0(output_unif_FCS_files, "wrongTubes.txt"), append = TRUE
            )
          }
          warning(paste0(file, ": Wrong tube! Should be tube ", tube_markers))
          wrongTubeFiles[[file]] <- tube_markers
          tube_correct <- tube_markers
          file <- gsub(listTubes[[as.character(tube)]], listTubes[[as.character(tube_markers)]], file)
          write(paste0(oldFile, "\t", file, "\t", tube, "\t", tube_markers, "\t"),
            file = paste0(output_unif_FCS_files, "wrongTubes.txt"), append = TRUE
          )
          y <- y + 1
          ff@parameters@data[["desc"]] <- markers_o_i
          break
        }
      }
      if (!all(markersInFF %in% markersInTube_tmp)) {
        warning(paste0("No correct tube found for ", oldFile, "!"))
        write(paste0("No correct tube found for ", oldFile, "!"),
          file = paste0(output_unif_FCS_files, "wrongTubes.txt"), append = TRUE
        )
        stopNow <- 1
      }
    } else {
      ff@parameters@data[["desc"]] <- markers_of_interest_all[[as.character(tube)]]
    }


    if (stopNow == 1) {
      stopNow <- 0
      next
    }

    if (is.null(ff@description$SPILL)) {
      file_tmp <- paste0(
        tmp$PatientID, "_", tmp$Panel, ifelse(tube == 0, "", tube_correct),
        ifelse(PBOorNOT, "_PBO", ""), ifelse(relapse_or_not, "_R", ""), ".fcs"
      )

      ff_orig <- read.FCS(paste0("./FCS_files/", basename(file_tmp)))
      ff@description$SPILL <- ff_orig@description$SPILL
    }

    colnames(ff@description$SPILL) <- gsub(" ", "-", colnames(ff@description$SPILL))
    colnames(ff@description$SPILL) <- gsub("PacB", "HV450", colnames(ff@description$SPILL))
    colnames(ff@description$SPILL) <- gsub("PacO", "HV500", colnames(ff@description$SPILL))
    ff@description$SPILL <- ff@description$SPILL[, channels_of_interest[4:11]]

    ff@description[["ORIGINALGUID"]] <- file
    ff@description[["$FIL"]] <- file
    ff@description[["GUID"]] <- file

    write.FCS(ff, filename = paste0(output_unif_FCS_files, file))
    y <- y + 1
  }
}


#-------------------------Parse uniformized files-------------------------------
unif_fcs_files <- paste0(
  output_unif_FCS_files,
  list.files(output_unif_FCS_files, pattern = ".*fcs$")
)

unif_samples <- data.frame(
  File = unif_fcs_files,
  Panel = stringr::str_match(unif_fcs_files, "ALOT|AML")[, 1],
  PatientID = stringr::str_match(unif_fcs_files, "([0-9]+-?[0-9])")[, 2],
  Tube = as.numeric(stringr::str_match(unif_fcs_files, "AML([1-7]*)")[, 2]),
  stringsAsFactors = FALSE
) %>%
  filter(!PatientID %in% patientsToRemove)

unif_samples$Tube[is.na(unif_samples$Tube)] <- 0

tmp <- as.factor(gsub("(^00|^15|^16|^17|^18|^19).*", "\\1", unif_samples$PatientID))
levels(tmp) <- c("2015", "2016", "2017", "2018", "2019")
unif_samples$Year <- tmp

#-------------------------Compensate, transform---------------------------------
manual_comp <- openxlsx::read.xlsx(paste0(input_compensation_matrices, "adapted_comp.xlsx"),
  colNames = FALSE
)
manual_comp[, 1] <- gsub("ALOT", "0", manual_comp[, 1])
manual_comp[, 1] <- gsub("tube", "", manual_comp[, 1])


for (tube in 0:6) {
  samples_tube_tmp <- unif_samples %>% filter(Tube == tube)
  samples_tube <- samples_tube_tmp %>% pull(File)
  manual_comp_tube <- manual_comp[manual_comp == tube, ]
  tfList <- readRDS(paste0(output_RDS, "tfList_tube_", tube, ".RDS"))

  for (samp in samples_tube) {
    title <- gsub("(.*/)((ALOT|AML)* *[0-9]{4}-?2?.*)\\.fcs$", "\\2", samp)
    print(title)
    ff <- read.FCS(samp)

    #----Compensation
    patientID_tmp <- samples_tube_tmp[samples_tube_tmp$File == samp, "PatientID"]
    if (patientID_tmp %in% manual_comp_tube[, 2]) {
      whichComp <- manual_comp_tube[manual_comp_tube[, 2] %in% patientID_tmp, 3]
      comp <- read.csv(paste0(input_compensation_matrices, whichComp, ".csv"))
      colnames(comp) <- gsub("\\.{4}.*", "", colnames(comp))
      comp$X <- NULL
      colnames(comp) <- gsub("\\.", "-", colnames(comp))
    } else {
      comp <- ff@description$SPILL
    }
    colnames(comp) <- gsub("PacB-A", "HV450-A", colnames(comp))
    colnames(comp) <- gsub("PacO-A", "HV500-A", colnames(comp))

    ff <- compensate(ff, comp)

    #----Transformation 
    ff <- transform(ff, tfList)

    #----Save results
    write.FCS(ff, paste0(output_prepro_FCS_files, title, "_prepro.fcs"))
  }
}

#----------------------Parse preprocessed files---------------------------------
fcs_files_prepro <- paste0(
  output_prepro_FCS_files,
  list.files(output_prepro_FCS_files, pattern = ".*\\.fcs$")
)
samples_prepro <- data.frame(
  File = fcs_files_prepro,
  Panel = stringr::str_match(fcs_files_prepro, "ALOT|AML")[, 1],
  PatientID = stringr::str_match(fcs_files_prepro, "([0-9]+-?[0-9])")[, 2],
  Tube = as.numeric(stringr::str_match(fcs_files_prepro, "AML([1-7]*)")[, 2]),
  stringsAsFactors = FALSE
) %>%
  filter(!PatientID %in% patientsToRemove)
samples_prepro$Tube[is.na(samples_prepro$Tube)] <- 0
tmp <- factor(substr(samples_prepro$PatientID, 1, 2))
levels(tmp) <- c("2015", "2016", "2017", "2018", "2019")
samples_prepro$Year <- tmp

#--------------------------------Normalisation----------------------------------
normalizationList <- list(
  "0" = c("CD45", "CD7"),
  "1" = c("CD45", "FSC-H"),
  "2" = c("CD45", "FSC-H"),
  "3" = c("CD45", "FSC-H"),
  "4" = c("CD45", "FSC-H"),
  "5" = c("CD45", "FSC-H"),
  "6" = c("CD45", "CD123", "FSC-H")
)
if (!dir.exists(output_norm_FCS_files)) dir.create(output_norm_FCS_files)
for (tube in seq(0, 6)) {
  message(paste0("Tube ", tube))
  samples_tube <- samples_prepro %>% dplyr::filter(Tube == tube)
  channelsToNorm <- GetChannels(
    read.FCS(samples_tube$File[2]),
    normalizationList[[as.character(tube)]]
  )
  ALOT_or_AML <- ifelse(tube == 0, "date_flow_acquisition_ALOT",
    "date_flow_acquisition_AML_tubes"
  )
  dateColMetadata <- metadata[, c("PatientID", ALOT_or_AML)] %>%
    dplyr::filter(PatientID %in% samples_tube$PatientID)
  dateColMetadata$markerChange <- factor(
    dateColMetadata[, ALOT_or_AML] > as.Date("2016-05-05")
  )
  samples_tube <- dplyr::inner_join(
    samples_tube, dateColMetadata,
    by = "PatientID"
  )

  normModel <- readRDS(paste0("./RDS/normModel_tube_markerchange_", tube, ".rds"))
  CytoNorm::QuantileNorm.normalize(
    model = normModel,
    files = samples_tube$File,
    labels = samples_tube$markerChange,
    transformList = NULL,
    transformList.reverse = NULL,
    outputDir = output_norm_FCS_files,
    prefix = "", truncate_max_range = FALSE
  )
  if (tube == 2) {
    fcs_files_norm <- paste0(
      output_norm_FCS_files,
      list.files(output_norm_FCS_files, pattern = ".*fcs$")
    )
    samples_tube <- data.frame(
      File = fcs_files_norm,
      Panel = stringr::str_match(fcs_files_norm, "ALOT|AML")[, 1],
      PatientID = stringr::str_match(fcs_files_norm, "([0-9]+-?[0-9])")[, 2],
      Tube = as.numeric(stringr::str_match(fcs_files_norm, "AML([1-7]*)")[, 2]),
      stringsAsFactors = FALSE
    ) %>% filter(Tube == 2)
    dateColMetadata <- metadata[, c("PatientID", "date_flow_acquisition_AML_tubes")] %>%
      dplyr::filter(PatientID %in% samples_tube$PatientID)
    dateColMetadata$markerChange <- factor(
      dateColMetadata[, ALOT_or_AML] > as.Date("2016-08-01")
    )
    samples_tube <- dplyr::inner_join(
      samples_tube, dateColMetadata,
      by = "PatientID"
    )

    channelsToNorm <- GetChannels(read.FCS(samples_tube$File[1]), c("CD64", "CD34", "CD117"))

    normModel <- readRDS(paste0("./RDS/normModel_tube_batch_tube_2_", tube, ".rds"))
    CytoNorm::QuantileNorm.normalize(
      model = normModel,
      files = samples_tube$File,
      labels = samples_tube$markerChange,
      transformList = NULL,
      transformList.reverse = NULL,
      outputDir = output_norm_FCS_files,
      prefix = "Norm_",
      removeOriginal = TRUE
    )
  }
  if (tube == 6) {
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
    ) %>% filter(Tube == 6)
    tmp <- factor(substr(samples_norm$PatientID, 1, 2))
    levels(tmp) <- c("2015", "2016", "2017", "2018", "2019")
    samples_norm$Year <- tmp == 2019
    channelsToNorm <- GetChannels(read.FCS(samples_norm$File[1]), c("CD123", "CD4"))

    normModel <- readRDS(paste0("./RDS/normModel_tube_batch_tube_6_", tube, ".rds"))
    CytoNorm::QuantileNorm.normalize(
      model = normModel,
      files = samples_norm$File,
      labels = samples_norm$Year,
      transformList = NULL,
      transformList.reverse = NULL,
      outputDir = output_norm_FCS_files,
      prefix = "Norm_",
      removeOriginal = TRUE
    )
  }
}
fcs_files_norm <- paste0(
  getwd(), output_norm_FCS_files,
  list.files(output_norm_FCS_files, pattern = ".*fcs$")
)
fcs_files_norm_newName <- gsub("prepro", "norm", fcs_files_norm)
fcs_files_norm_newName <- gsub("Norm_", "", fcs_files_norm_newName)
file.rename(fcs_files_norm, fcs_files_norm_newName)

#------------------------Parse normalized files---------------------------------
fcs_files_norm <- paste0(
  output_norm_FCS_files,
  list.files("./blasts/Normalized_fcs/", pattern = ".*fcs$")
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
  dplyr::filter(File != paste0(output_norm_FCS_files, "1747_ALOT_PBO_norm.fcs"))

#-----------------------------Plot filescatters---------------------------------
if (!dir.exists(paste0(output_norm_FCS_files, "FileScatters"))) {
  dir.create(paste0(output_norm_FCS_files, "FileScatters"))
  dir.create(paste0(output_norm_FCS_files, "FileScatters/Preprocessed"))
  dir.create(paste0(output_norm_FCS_files, "FileScatters/Normalized"))
  dir.create(paste0(output_norm_FCS_files, "FileScatters/Composed"))
}

normalizationList <- list(
  "0" = c("CD45", "CD7"),
  "1" = c("CD45"),
  "2" = c("CD45", "CD64", "CD34", "CD117"),
  "3" = c("CD45"),
  "4" = c("CD45"),
  "5" = c("CD45"),
  "6" = c("CD45", "CD123", "CD4")
)


for (tube in seq(0, 6)) {
  # Preprocessed
  message(paste0("Tube ", tube))
  cat("Plot preprocessed\n")
  samples_tube <- samples_prepro %>% filter(Tube == tube)
  ALOT_or_AML <- ifelse(tube == 0, "date_flow_acquisition_ALOT",
    "date_flow_acquisition_AML_tubes"
  )
  dateColMetadata <- metadata[, c("PatientID", ALOT_or_AML)] %>%
    filter(PatientID %in% samples_tube$PatientID)

  dateColMetadata$afterMarkerChange <- factor(dateColMetadata[, ALOT_or_AML] >
    as.Date("2016-05-05"))
  tmp <- inner_join(samples_tube, dateColMetadata, "PatientID")
  tmp <- tmp[order(tmp[, ALOT_or_AML]), ]

  ff <- read.FCS(tmp$File[1])
  norm_channels <- GetChannels(ff, normalizationList[[as.character(tube)]])
  channelnames <- colnames(ff@exprs)
  channelnames <- channelnames[which(!channelnames %in% c("Time", "FSC-H"))]
  PlotFileScatters(tmp$File,
    names = tmp$PatientID, groups = tmp$afterMarkerChange,
    channels = channelnames,
    plotFile = paste0(
      output_norm_FCS_files,
      "FileScatters/Preprocessed/FileScatters_preproc_markerChange",
      tube, ".png"
    ),
    ncol = 1, yLabel = c("marker", "channel"), silent = TRUE
  )

  PlotFileScatters(tmp$File,
    names = tmp$PatientID, groups = tmp$Year,
    channels = channelnames,
    plotFile = paste0(
      output_norm_FCS_files,
      "FileScatters/Preprocessed/FileScatters_preproc_year",
      tube, ".png"
    ),
    ncol = 1, yLabel = c("marker", "channel"),
    quantiles = c(0.5), silent = TRUE
  )

  p_pp <- PlotFileScatters(tmp$File,
    names = tmp$PatientID, groups = tmp$Year,
    channels = norm_channels,
    plotFile = NULL,
    yLabel = c("marker", "channel"), silent = TRUE
  )

  # Normalized
  if (nrow(samples_norm) == 1) next
  tube_samples <- samples_norm %>% filter(Tube == tube)
  ALOT_or_AML <- ifelse(tube == 0, "date_flow_acquisition_ALOT",
    "date_flow_acquisition_AML_tubes"
  )
  dateColMetadata <- metadata[, c("PatientID", ALOT_or_AML)] %>%
    filter(PatientID %in% tube_samples$PatientID)

  dateColMetadata$afterMarkerChange <- factor(dateColMetadata[, ALOT_or_AML] >
    as.Date("2016-05-05"))
  tmp <- inner_join(tube_samples, dateColMetadata, "PatientID")
  tmp <- tmp[order(tmp[, ALOT_or_AML]), ]

  ff <- read.FCS(tmp$File[1])
  channelnames <- colnames(ff@exprs)
  channelnames <- channelnames[which(!channelnames %in% c("Time", "FSC-H"))]
  cat("Plot normalized\n")
  PlotFileScatters(tmp$File,
    names = tmp$PatientID, groups = tmp$afterMarkerChange,
    channels = channelnames,
    plotFile = paste0(
      output_norm_FCS_files,
      "FileScatters/Normalized/FileScatters_norm_markerChange",
      tube, ".png"
    ),
    ncol = 1, yLabel = c("marker", "channel"), silent = TRUE
  )

  PlotFileScatters(tmp$File,
    names = tmp$PatientID, groups = tmp$Year,
    channels = channelnames,
    plotFile = paste0(
      output_norm_FCS_files,
      "FileScatters/Normalized/FileScatters_norm_year",
      tube, ".png"
    ),
    ncol = 1, yLabel = c("marker", "channel"),
    quantiles = c(0.5), silent = TRUE
  )

  p_norm <- PlotFileScatters(tmp$File,
    names = tmp$PatientID, groups = tmp$Year,
    channels = norm_channels,
    plotFile = NULL,
    yLabel = c("marker", "channel"), silent = TRUE
  )

  cat("Plot composed\n")
  png(paste0(output_norm_FCS_files, "FileScatters/Composed/FileScatters_comp_year_", tube, ".png"),
    width = length(norm_channels) * (60 + 15 * length(tmp$File)),
    height = 500
  )
  p <- ggpubr::ggarrange(
    plotlist = c(p_pp, p_norm),
    nrow = 2, ncol = length(norm_channels),
    common.legend = TRUE
  )
  print(p)
  dev.off()
}



#----------------------------Plot PCA of quantiles------------------------------
if (!dir.exists(paste0(output_norm_FCS_files, "PCA_on_quantiles"))) {
  dir.create(paste0(output_norm_FCS_files, "PCA_on_quantiles"))
  dir.create(paste0(output_norm_FCS_files, "PCA_on_quantiles/Preprocessed"))
  dir.create(paste0(output_norm_FCS_files, "PCA_on_quantiles/Normalized"))
}

for (tube in seq(0, 6)) {
  message(paste0("Tube ", tube))
  # Preprocessed
  samples_tube <- samples_prepro %>% filter(Tube == tube)
  quantiles <- matrix(NA,
    nrow = nrow(samples_tube),
    ncol = 0,
    dimnames = list(
      samples_tube$PatientID,
      NULL
    )
  )
  agg <- AggregateFlowFrames(samples_tube$File,
    cTotal = 3000000, keepOrder = TRUE,
    silent = TRUE
  )
  channels <- colnames(agg)[which(!colnames(agg) %in% c(
    "Time", "File", "File_scattered", "Original_ID", "FSC-H"
  ))]

  ALOT_or_AML <- ifelse(tube == 0, "date_flow_acquisition_ALOT",
    "date_flow_acquisition_AML_tubes"
  )
  dateColMetadata <- metadata[, c("PatientID", ALOT_or_AML)] %>%
    filter(PatientID %in% samples_tube$PatientID)

  dateColMetadata$afterMarkerChange <- factor(dateColMetadata[, ALOT_or_AML] >
    as.Date("2016-05-05"))
  tmp <- inner_join(samples_tube, dateColMetadata, "PatientID")
  tmp <- tmp[order(tmp[, ALOT_or_AML]), ]

  for (channel in channels) {
    x <- exprs(agg)[, channel]
    q <- tapply(
      x, agg@exprs[, "File"],
      quantile, c(0.25, 0.5, 0.75)
    )
    q <- do.call(rbind, q)
    colnames(q) <- paste0(c("q25_", "q50_", "q75_"), channel)
    quantiles <- cbind(quantiles, q)
  }

  quantiles <- quantiles[, apply(quantiles, 2, var) != 0]
  pca <- prcomp(scale(quantiles))
  plotdf <- data.frame(pca$x[, 1:2], "PatientID" = rownames(pca$x[, 1:2])) %>%
    dplyr::inner_join(tmp, by = "PatientID")
  plotdf$PBO <- grepl("PBO", plotdf$File)
  png(paste0(output_norm_FCS_files, "PCA_on_quantiles/Preprocessed/PCA_q_prepro_markerChange_", tube, ".png"),
    width = 500, height = 500
  )
  p <- ggplot(plotdf) +
    geom_point(aes(x = PC1, y = PC2, col = afterMarkerChange)) +
    theme_minimal() +
    ggtitle(paste0("Tube ", tube))
  p <- ggExtra::ggMarginal(p, groupFill = TRUE, alpha = 0.5)
  print(p)
  dev.off()

  png(paste0(output_norm_FCS_files, "PCA_on_quantiles/Preprocessed/PCA_q_prepro_year_", tube, ".png"),
    width = 500, height = 500
  )
  p <- ggplot(plotdf) +
    geom_point(aes(x = PC1, y = PC2, col = Year)) +
    theme_minimal() +
    ggtitle(paste0("Tube ", tube))
  p <- ggExtra::ggMarginal(p, groupFill = TRUE, alpha = 0.5)
  print(p)
  dev.off()

  png(paste0(output_norm_FCS_files, "PCA_on_quantiles/Preprocessed/PCA_prepro_var_", tube, ".png"),
    width = 500, height = 500
  )
  p <- factoextra::fviz_pca_var(pca,
    col.var = "contrib",
    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
    repel = TRUE
  )
  print(p)
  dev.off()

  # Normalized
  samples_tube <- samples_norm %>% filter(Tube == tube)
  quantiles <- matrix(NA,
    nrow = nrow(samples_tube),
    ncol = 0,
    dimnames = list(
      samples_tube$PatientID,
      NULL
    )
  )
  agg <- AggregateFlowFrames(samples_tube$File,
    cTotal = 3000000, keepOrder = TRUE,
    silent = TRUE
  )
  channels <- colnames(agg)[which(!colnames(agg) %in% c(
    "Time", "File", "File_scattered", "Original_ID", "FSC-H"
  ))]

  ALOT_or_AML <- ifelse(tube == 0, "date_flow_acquisition_ALOT",
    "date_flow_acquisition_AML_tubes"
  )
  dateColMetadata <- metadata[, c("PatientID", ALOT_or_AML)] %>%
    filter(PatientID %in% samples_tube$PatientID)

  dateColMetadata$afterMarkerChange <- factor(dateColMetadata[, ALOT_or_AML] >
    as.Date("2016-05-05"))
  tmp <- inner_join(samples_tube, dateColMetadata, "PatientID")
  tmp <- tmp[order(tmp[, ALOT_or_AML]), ]

  for (channel in channels) {
    x <- exprs(agg)[, channel]
    q <- tapply(
      x, agg@exprs[, "File"],
      quantile, c(0.25, 0.5, 0.75)
    )
    q <- do.call(rbind, q)
    colnames(q) <- paste0(c("q25_", "q50_", "q75_"), channel)
    quantiles <- cbind(quantiles, q)
  }

  quantiles <- quantiles[, apply(quantiles, 2, var) != 0]
  pca <- prcomp(scale(quantiles))
  plotdf <- data.frame(pca$x[, 1:2], "PatientID" = rownames(pca$x[, 1:2])) %>%
    dplyr::inner_join(tmp, by = "PatientID")
  plotdf$PBO <- grepl("PBO", plotdf$File)

  png(paste0(output_norm_FCS_files, "PCA_on_quantiles/Normalized/PCA_q_norm_markerChange_", tube, ".png"),
    width = 500, height = 500
  )
  p <- ggplot(plotdf) +
    geom_point(aes(x = PC1, y = PC2, col = afterMarkerChange)) +
    theme_minimal() +
    ggtitle(paste0("Tube ", tube))
  p <- ggExtra::ggMarginal(p, groupFill = TRUE, alpha = 0.5)
  print(p)
  dev.off()

  png(paste0(output_norm_FCS_files, "PCA_on_quantiles/Normalized/PCA_q_norm_year_", tube, ".png"),
    width = 500, height = 500
  )
  p <- ggplot(plotdf) +
    geom_point(aes(x = PC1, y = PC2, col = Year)) +
    theme_minimal() +
    ggtitle(paste0("Tube ", tube))
  p <- ggExtra::ggMarginal(p, groupFill = TRUE, alpha = 0.5)
  print(p)
  dev.off()

  png(paste0(output_norm_FCS_files, "PCA_on_quantiles/Normalized/PCA_norm_var_", tube, ".png"),
    width = 500, height = 500
  )
  p <- factoextra::fviz_pca_var(pca,
    col.var = "contrib",
    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
    repel = TRUE
  )
  print(p)
  dev.off()
}
