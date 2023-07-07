###############################################################################/
# LIBRARIES
###############################################################################/
library(flowCore)
library(FlowSOM)
library(Rtsne)
library(dplyr)   
library(ggplot2) 
library(flowDensity)
source("./3_helperFunctions.R")

###############################################################################/
# SCRIPT
###############################################################################/
#----------------------------------Define variables-----------------------------
seed <- 2022
set.seed(seed)

input_raw_FCS_files <- "./raw_FCS_files/"
input_compensation_matrices <- "./comp_matr/"
input_metadata <- "./Excel_metadata/" 

output_unif_FCS_files <- "./FCS_files/"
output_prepro_FCS_files <- "./Preprocessed_fcs/"
output_norm_FCS_files <- "./Normalized_fcs/"
output_RDS <- "./RDS/"

#--------------------------------Read metadata----------------------------------
metadata <- openxlsx::read.xlsx(paste0(input_metadata, "full_database_02082022.xlsx"), 
                                check.names = FALSE, 
                                sep.names = "_") 

for (x in colnames(metadata)[grepl("date|Date", colnames(metadata))]){
  y <- convDate(metadata[, x])
  metadata[, x] <- y
}

non_intensive_chemotherapy_patients <- metadata %>% 
  filter(treatment_modality != "intensive chemotherapy") %>% 
  pull(PatientID)
noAML_patients <- c("1861", "1870", "1873", "1901", "1903", "1914", "1919", 
                    "1922", "1923", "1934", "1939", "1943", "1947", "0011")
relapse_patients <- c("1925")
noInfo_patients <- c("1908", "1910", "1911", "1912", "1928", "1942", "1946", 
                     "1734")
tooYoung_patients <- c("1866", "1868", "1869", "1605", "1906")
secAML_or_otherDiagnosis_patients <- c("1616", "1723", "1812", "1904", "1740", 
                                       "1711", "1879")
strangeFlowData_patient <- c("1951")
patientsToRemove <- c(noAML_patients, noInfo_patients, tooYoung_patients, 
                      secAML_or_otherDiagnosis_patients, strangeFlowData_patient,
                      relapse_patients, non_intensive_chemotherapy_patients)

#--------------------------------File 1937 AML3---------------------------------
# File 1937 AML3 raised errors in read.FCS from flowCore, so we had to manually 
# adapt the function to read it in and save it
ff <- readFCS(paste0(input_raw_FCS_files, "Anonymous file 6_Export_File 1937 AML3.fcs"))
write.FCS(ff, paste0(input_raw_FCS_files, "Anonymous file 6_Export_File 1937 AML3.fcs"))

#--------------------------------Parse raw data---------------------------------
file.rename("1625 AML 5.fcs", "1625 AML5.fcs")
file.rename("1866 AML 5.fcs", "1866 AML5.fcs")
file.rename("1747 AML 5.fcs", "1747 AML5.fcs")

fcs_files <- list.files(input_raw_FCS_files, pattern = ".*fcs$")

samples <- data.frame(File = fcs_files,
                      Panel = stringr::str_match(fcs_files, "ALOT|AML")[,1],
                      PatientID = stringr::str_match(fcs_files, "([0-9]{4}-?[0-9]?)")[,2],
                      Tube = as.numeric(stringr::str_match(fcs_files, "(?<=Tube|AML)[1-7]")),
                      stringsAsFactors = FALSE,
                      PBOorNOT = grepl("PBO", fcs_files))

samples$Tube[is.na(samples$Tube)] <- 0

#-----------------------------Uniformize data-----------------------------------
channels_of_interest <- c("Time", "FSC-A","FSC-H","SSC-A", "FITC-A", "PE-A",
                          "PerCP-Cy5-5-A","PE-Cy7-A","APC-A","APC-H7-A","HV450-A",
                          "HV500-A")

markers_of_interest_all <- list("0" = c("FSC-A", "FSC-H", "SSC-A", "cyMPO", 
                                        "cyCD79a", "CD34", "CD19", "CD7", "CD3", 
                                        "cyCD3", "CD45"),
                                "1" = c("FSC-A","FSC-H", "SSC-A", "CD16", "CD13",
                                        "CD34", "CD117", "CD11b", "CD10", "HLA-DR",
                                        "CD45"),
                                "2" = c("FSC-A", "FSC-H", "SSC-A", "CD35", "CD64",
                                        "CD34", "CD117", "CD300e", "CD14", "HLA-DR",
                                        "CD45"),
                                "3" = c("FSC-A", "FSC-H", "SSC-A", "CD36", "CD105", 
                                        "CD34", "CD117", "CD33", "CD71", "HLA-DR",
                                        "CD45"),
                                "4" = c("FSC-A",  "FSC-H", "SSC-A", "nuTdT",
                                        "CD56", "CD34", "CD117", "CD7", "CD19",
                                        "HLA-DR", "CD45"),
                                "5" = c("FSC-A", "FSC-H", "SSC-A", "CD15", "NG2",
                                        "CD34", "CD117", "CD22", "CD38", "HLA-DR", 
                                        "CD45"),
                                "6" = c("FSC-A", "FSC-H", "SSC-A", "CD42a+CD61",
                                        "CD203c", "CD34", "CD117", "CD123", "CD4", 
                                        "HLA-DR", "CD45"),
                                "7" = c("FSC-A", "FSC-H", "SSC-A", "CD41", "CD25", 
                                        "CD34", "CD117","CD42b", "CD9", "HLA-DR",
                                        "CD45"))
listTubes <- list("0" = "_ALOT", "1" = "_AML1", "2" = "_AML2", "3" = "_AML3",
                  "4" = "_AML4", "5" = "_AML5", "6" = "_AML6", "7" = "_AML7")

y <- 1
stopNow <- 0
wrongTubeFiles <- list()
for (tube in 0:7){
  samples_tube <- samples %>% filter(Tube == tube)
  markers_o_i <- markers_of_interest_all[[as.character(tube)]]
  markersInTube <- tolower(sortMarkers(markers_o_i[4:11]))
  for(file in samples_tube$File){

    print(paste0("Tube ", tube, ": ", y, " of ,", nrow(samples)))
    oldFile <- file
    ff <- read.FCS(file)
    tmp <- samples[samples$File == file, ]
    names_ff <- unname(ff@parameters@data[["name"]])
    desc_ff <- unname(ff@parameters@data[["desc"]])
    colnames(ff) <- gsub(" ", "-", colnames(ff))
    
    # PacB and PacO are changed to HV450 and HV500 respectively
    colnames(ff) <- gsub("PacB", "HV450", colnames(ff))
    colnames(ff) <- gsub("PacO", "HV500", colnames(ff))
    colnames(ff@description$SPILL) <- gsub(" ", "-", colnames(ff@description$SPILL))
    colnames(ff@description$SPILL) <- gsub("PacB", "HV450", colnames(ff@description$SPILL))
    colnames(ff@description$SPILL) <- gsub("PacO", "HV500", colnames(ff@description$SPILL))
    ff <- ff[, channels_of_interest]
    markersInFF <- ff@parameters@data[["desc"]]
    
    # Remove -, + and spaces in the marker names  and make marker names uniform
    markersInFF <- gsub(" |-|+", "", markersInFF)
    markersInFF <- gsub("HV500C|HV450", "", markersInFF) #1747 ALOT
    markersInFF <- gsub("cyCD79a", "cyCD79a", markersInFF) #1747 ALOT
    markersInFF <- gsub("IREM2", "CD300e", markersInFF) 
    markersInFF <- gsub("Cy 79a|Cy79a", "cyCD79a", markersInFF) 
    markersInFF <- gsub("CyCD3", "cyCD3", markersInFF)  
    markersInFF <- gsub("HLADR", "HLA-DR", markersInFF)  
    ff@description$SPILL <- ff@description$SPILL[, channels_of_interest[5:12]]
    ff@parameters@data[["desc"]] <- markersInFF
    
    # Check if ff markers resemble markers from tube in title and if not change it
    markersInFF <- markersInFF[!is.na(markersInFF)]
    markersInFF <- tolower(sortMarkers(markersInFF))
    relapse_or_not <- grepl("_R", file)
    PBOorNOT <- grepl("PB", file)
    file <- paste0(tmp$PatientID, "_", tmp$Panel, ifelse(tube == 0, "", tmp$Tube), 
                   ifelse(PBOorNOT, "_PBO", ""), ifelse(relapse_or_not, "_R", ""), ".fcs")
    
    if (!all(markersInTube %in% markersInFF)){
      for (tube_markers in names(markers_of_interest_all)){
        markers_o_i <- markers_of_interest_all[[tube_markers]]
        markersInTube_tmp <- tolower(sortMarkers(markers_o_i[4:11]))
        
        # Save wrong tubes in a file
        if (all(markersInFF %in% markersInTube_tmp)){
          if (!file.exists(paste0(output_unif_FCS_files, "wrongTubes.txt"))){
            write("Files in wrong tubes \n------------------------\n", 
                  file = paste0(output_unif_FCS_files, "wrongTubes.txt"), append = TRUE)
            write("Old filename\tNew Filename\tOld tube\tNew tube",
                  file = paste0(output_unif_FCS_files, "wrongTubes.txt"), append = TRUE)
          }
          warning(paste0(file, ": Wrong tube! Should be tube ", tube_markers))
          wrongTubeFiles[[file]] <- tube_markers
          file <- gsub(listTubes[[as.character(tube)]], listTubes[[as.character(tube_markers)]], file)
          write(paste0(oldFile, "\t", file, "\t", tube, "\t", tube_markers, "\t"),
                file = paste0(output_unif_FCS_files, "wrongTubes.txt"), append = TRUE)
           y <- y + 1
           ff@parameters@data[["desc"]][2:12] <- markers_o_i
          break
        }
      } else {
        ff@parameters@data[["desc"]][2:12] <- markers_of_interest_all[[tube]]
      }
      if (!all(markersInFF %in% markersInTube_tmp)){
        warning(paste0("No correct tube found for ", oldFile, "!"))
        write(paste0("No correct tube found for ", oldFile, "!"),
              file = paste0(output_unif_FCS_files, "wrongTubes.txt"), append = TRUE)
        stopNow <- 1
      }
    }
    
    if (stopNow == 1){
      stopNow <- 0
      next
    }

    ff@description[["ORIGINALGUID"]] <- file
    ff@description[["$FIL"]] <- file
    ff@description[["GUID"]] <- file
    
    write.FCS(ff, filename = paste0(output_unif_FCS_files, file))
    y <- y + 1
  }
}

#-------------------------Parse uniformized files-------------------------------
unif_fcs_files <- paste0(output_unif_FCS_files, 
                    list.files(output_unif_FCS_files, pattern = ".*fcs$"))

unif_samples <- data.frame(File = unif_fcs_files,
                      Panel = stringr::str_match(unif_fcs_files, "ALOT|AML")[,1],
                      PatientID = stringr::str_match(unif_fcs_files, "([0-9]+-?[0-9])")[, 2],
                      Tube = as.numeric(stringr::str_match(unif_fcs_files, "AML([1-7]*)")[,2]),
                      stringsAsFactors = FALSE)  %>% 
  filter(!PatientID %in% patientsToRemove)  

unif_samples$Tube[is.na(unif_samples$Tube)] <- 0

tmp <- as.factor(gsub("(^00|^15|^16|^17|^18|^19).*","\\1",unif_samples$PatientID))
levels(tmp) <- c("2015", "2016", "2017", "2018", "2019")
unif_samples$Year <- tmp

#-------Compensate, transform, remove margin events, doublets and debris--------
manual_comp <-  openxlsx::read.xlsx(paste0(input_compensation_matrices, "adapted_comp.xlsx"), 
                                    colNames = FALSE)
manual_comp[, 1] <- gsub("ALOT", "0", manual_comp[, 1])
manual_comp[, 1] <- gsub("tube", "", manual_comp[, 1])

AMLs <- paste0("_", c("ALOT", paste0("AML", seq_len(6))))
more_strict_doublets <- c("0028_AML1", "0028_AML2", "0028_AML3", "0028_AML4", 
                          "0028_AML5", "0028_AML6", "0045_AML1", "0045_AML2",  
                          "0045_AML3", "0045_AML5", "0045_AML6", "1605_AML5", 
                          "1618_AML5", "1716_AML5", "1724_AML1", "1724_AML2", 
                          "1729_AML5", "1746_AML6", "1867_AML1", "1867_AML2", 
                          "1867_AML3", "1867_AML4", "1867_AML5", "1870_AML5", 
                          "1909_AML5", "1947_AML1", "1722_ALOT", "1809_AML5", 
                          "1706_AML6" , "0025_AML1", 
                          paste0(1807, AMLs), paste0(1810, AMLs), paste0("0022", AMLs),
                          paste0(1702, AMLs), paste0(1716, AMLs), paste0("0010", AMLs))
even_more_strict_doublets <- c("0026_AML1", "0026_AML2", "0026_AML3", "0026_AML4", 
                               "0026_AML5", "0026_AML6", "1910_AML5",
                               paste0("0021", AMLs), paste0("0029", AMLs))
less_strict_doublets <- c("1863_ALOT_PBO", "0032_ALOT", "0039_ALOT", "1717_ALOT_PBO", 
                          "1739_AML2", "1740_ALOT_PBO", "1742_ALOT", "1726_ALOT", 
                          "1726_AML1", "1726_AML2", "1726_AML3", "1726_AML4", 
                          "1726_AML5", "1726_AML6", "1607_ALOT", "1617_ALOT", 
                          "1622_ALOT", "1623_ALOT", "1627_ALOT", "1707_ALOT", 
                          "1703_ALOT", "1705_ALOT", "1707_ALOT", "1711_ALOT_PBO",
                          "1713_ALOT", "1719_ALOT_PBO", "1720_ALOT", "1727_ALOT",
                          "1744_ALOT", "1749_ALOT_PBO", "1751_ALOT", "1752_ALOT", 
                          "1927_ALOT_PBO", "1937_AML1", "1937_AML2", "1944_ALOT", 
                          "1946_ALOT_PBO", "1947_ALOT", "0007-2_ALOT", "0016_ALOT", 
                          "0019_ALOT", "0020_ALOT", "1601_ALOT", "1701_ALOT", "1817_ALOT",
                          paste0(1710, AMLs), paste0(1743, AMLs), paste0(1802, AMLs),
                          paste0(1875, AMLs), paste0(1753, AMLs), 
                          paste0(1935, paste0(AMLs, "_PBO")))
even_less_strict_doublets <- c("1735_ALOT", "0014_ALOT", "0043_ALOT",
                               "1716_ALOT_PBO", "1739_ALOT", "1739_AML1")
channels_of_interest <- c("Time", "FSC-A","FSC-H","SSC-A", "FITC-A", "PE-A",
                          "PerCP-Cy5-5-A","PE-Cy7-A","APC-A","APC-H7-A","HV450-A",
                          "HV500-A")
customDebris <- list("0009_AML3" = c(9, 18),
                     "0009_AML6" = c(80, 91),
                     "0014_AML2" = c(64, 73, 74),
                     "0014_AML3" = c(72, 77, 80, 81),
                     "0014_AML5" = c(18, 73:75),
                     "0014_AML6" = c(1, 10),
                     "0019_AML3" = c(61:62, 67:69, 72, 77, 81),
                     "0019_AML5" = c(52, 61, 69, 77, 78, 80, 81),
                     "0021_AML3" = c(9, 16:18, 25, 42, 43, 73),
                     "0021_AML6" = c(6, 8, 9, 45),
                     "0022_AML3" = c(12, 57),
                     "0023_AML5" = c(78, 79),
                     "0023_AML6" = c(9, 37, 38, 44, 45, 54),
                     "0025_AML1" = c(46, 47, 56, 73, 74),
                     "0025_AML2" = c(37, 38, 47, 57, 59, 65, 73),
                     "0025_AML3" = c(5, 8, 9, 14, 15, 24, 25, 26, 55),
                     "0025_AML5" = c(5:7, 16, 24:26, 33),
                     "0025_AML6" = c(6, 7, 14, 24:27, 35, 36),
                     "0032_AML1" = c(8, 27),
                     "0032_AML2" = c(8, 18, 27),
                     "0032_AML3" = c(2, 10),
                     "0032_AML4" = c(3, 4, 10, 12, 13, 19, 30, 31),
                     "0032_AML5" = c(53, 54, 63, 71),
                     "0032_AML6" = c(7, 16),
                     "0037_AML5" = c(1, 37, 38),
                     "0037_AML6" = c(1, 73, 74),
                     "1714_ALOT" = c(1:3, 6, 10, 11, 16, 19, 20, 29, 36:39, 46:47),
                     "1716_AML1" = c(19, 70:72, 78:81),
                     "1716_AML2" = c(8, 9, 17, 18, 27),
                     "1716_AML3" = c(51, 54, 63, 71, 72, 80, 81),
                     "1716_AML5" = c(8, 9, 17, 18, 27, 45),
                     "1716_AML6" = c(48, 73:76),
                     "1739_ALOT" = c(9, 55, 72, 80),
                     "1739_AML4" = c(9, 46),
                     "1739_AML5" = c(9, 74),
                     "1739_AML6" = c(2, 9),
                     "1802_AML4" = c(10),
                     "1864_AML1" = c(9, 18, 27, 51, 53),
                     "1864_AML2" = c(51, 63, 80),
                     "1864_AML3" = c(29),
                     "1864_AML4" = c(64),
                     "1864_AML5" = c(18, 43),
                     "1864_AML6" = c(10, 29, 43, 46),
                     "1885_AML1" = c(71, 72, 78, 79),
                     "1885_AML4" = c(19, 20),
                     "1885_AML6" = c(2, 11),
                     "1888_AML2" = c(18),
                     "1905_AML1" = c(9, 17, 46),
                     "1905_AML2" = c(69),
                     "1905_AML3" = c(13, 19),
                     "1905_AML4" = c(35, 36, 59),
                     "1905_AML5" = c(16, 24, 35, 43),
                     "1905_AML6" = c(23),
                     "1954_AML1" = c(1:4, 10, 11, 19:21, 23),
                     "1954_AML2" = c(8, 9, 15, 16),
                     "1954_AML3" = c(7:9, 16, 18, 20, 44),
                     "1954_AML5" = c(37, 39, 45, 47, 55:57, 66, 74, 75),
                     "1954_AML6" = c(9, 18, 24, 27, 35, 36, 74))

reprFile <- list("0" = "1627",
                 "1" = "1953",
                 "2" = "1953",
                 "3" = "1953",
                 "4" = "1953",
                 "5" = "1953", 
                 "6" = "1953")

resetRDS <- FALSE

for (tube in seq(0, 6)){
  message(paste("--------------------Tube", tube, "--------------------"))
  samples_tube_tmp <- samples_nc %>% filter(Tube == tube) 
  samples_tube <- samples_tube_tmp %>% pull(File)
  manual_comp_tube <- manual_comp[manual_comp == tube, ]
  query <- as.list(rep("low", length(channels_of_interest[-1])))
  names(query) <- channels_of_interest[-1]
  if (tube == 0){
    colsToUseDebris <- channels_of_interest[-c(1, 3)]
  } else {
    colsToUseDebris <-  c("FSC-A", "SSC-A", "HV450-A", 
                          "HV500-A", "PerCP-Cy5-5-A", "PE-Cy7-A") #Scatters and backbone markers
  }
  
  #----Determine transformlist
  ff <- read.FCS(samples_tube[grepl(reprFile[[as.character(tube)]],
                                    samples_tube)])
  ff <- PeacoQC::RemoveMargins(ff, c("FSC-A", "SSC-A", "FSC-H"))
  
  if (samples_tube[1] %in% manual_comp_tube[, 2]){
    whichComp <- manual_comp_tube[manual_comp_tube[, 2] %in% samples_tube[1], 3]
    comp <- read.csv(paste0("./comp_matr/", whichComp, ".csv"))
    colnames(comp) <- gsub("\\.{4}.*", "", colnames(comp))
    comp$X <- NULL
    colnames(comp) <- gsub("\\.", "-", colnames(comp))
  } else comp <- ff@description$SPILL
  
  ff <- compensate(ff, comp)
  tfList <- estimateLogicle(ff, channels = colnames(ff@description$SPILL))
  ff_t <- transform(ff, tfList)
  q01_goal <- median(apply(ff_t@exprs[, -c(1:4)], 2, quantile, 0.01))
  q99_goal <- median(apply(ff_t@exprs[, -c(1:4)], 2, quantile, 0.99))
  for (scatter in c("FSC-A", "FSC-H", "SSC-A")){
    
    q01 <- quantile(ff@exprs[, scatter], 0.01)
    q99 <- quantile(ff@exprs[, scatter], 0.99)
    lTransform <- linearTransform(a = (q99_goal - q01_goal) / (q99 - q01), 
                                  b = (q01_goal * q99 - q99_goal * q01) /
                                    (q99 - q01))
    lTransform<- transformList(scatter, lTransform)
    tfList <- c(tfList, lTransform)
  }
  saveRDS(tfList, paste0("../../../group/irc/personal/artuurc/Preprocessing/RDS/tfList_tube_", tube, ".RDS"))  
  
  # #----Determine cutoffs for debrisClusterfinder function on aggregate
  if (!file.exists(paste0(
    "../../../group/irc/personal/artuurc/Preprocessing/FCS/ff_agg_tube_",
    tube,
    ".fcs")) | resetRDS){
    ff_agg <- AggregateFlowFrames(samples_tube,
                                  cTotal = 50000 * length(samples_tube),
                                  keepOrder = TRUE,
                                  writeOutput = FALSE,
                                  silent = TRUE)
    # Remove margins
    ff_agg <- PeacoQC::RemoveMargins(ff_agg, c("FSC-A", "SSC-A", "FSC-H"))
    # Transformation
    ff_agg <- transform(ff_agg, tfList)
    # PeacoQC
    ff_agg <- suppressWarnings(suppressMessages(PeacoQC::PeacoQC(
      ff = ff_agg,
      channels = colnames(ff_agg)[-which(colnames(ff_agg) == "Time")],
      plot = FALSE,
      save_fcs = FALSE,
      report = FALSE,
      output_directory = NULL)$FinalFF
    ))
    write.FCS(ff_agg,
              paste0("../../../group/irc/personal/artuurc/Preprocessing/FCS/ff_agg_tube_",
                     tube,".fcs"))
  } else {
    ff_agg <- read.FCS(paste0(
      "../../../group/irc/personal/artuurc/Preprocessing/FCS/ff_agg_tube_",
      tube,".fcs"))
  }
  FSC_cutoff_agg <- min(deGate(ff_agg, "FSC-A", all.cuts = TRUE, upper = FALSE))
  SSC_cutoff_agg <- min(deGate(ff_agg, "SSC-A", all.cuts = TRUE, upper = TRUE))
  selection_highSSC_agg <- exprs(ff_agg)[,"SSC-A"] > SSC_cutoff_agg + 1
  ff_highSSC_agg <- ff_agg[selection_highSSC_agg,]
  FSC_cutoff_highSSC_agg <- max(deGate(ff_highSSC_agg, "FSC-A", all.cuts = TRUE,
                                       upper = TRUE))
  saveRDS(list("FSC_cutoff_agg" = FSC_cutoff_agg,
               "SSC_cutoff_agg" = SSC_cutoff_agg,
               "FSC_cutoff_highSSC_agg" = FSC_cutoff_highSSC_agg),
          paste0("../../../group/irc/personal/artuurc/Preprocessing/RDS/cutoffs_tube_", tube, ".RDS")) 
  #----Loop over every file
  parallel::mclapply(samples_tube, function(samp){
    title <- gsub(".*([0-9]{4}-?2?.*).fcs","\\1", samp)
    print(title)
    ff <- read.FCS(samp)
    #----Remove Margin events
    indices_margins <- suppressWarnings(PeacoQC::RemoveMargins(
      ff, c("FSC-A", "SSC-A", "FSC-H"), output = "full")$indices_margins)
    selection_inrange <- !seq_len(nrow(ff)) %in% indices_margins 
    selection <- selection_inrange
    perc_remove_margins <- round((sum(!selection_inrange) / 
                                    length(selection_inrange)) * 100, 2)
    cat(paste0("RemoveMargins removed ", perc_remove_margins, "% of cells.\n"))
    ff1 <- ff[selection_inrange, ]
    
    #----Compensation
    patientID_tmp <- samples_tube_tmp[samples_tube_tmp$File == samp, "PatientID"]
    if (patientID_tmp %in% manual_comp_tube[, 2]){
      whichComp <- manual_comp_tube[manual_comp_tube[, 2] %in% patientID_tmp, 3]
      comp <- read.csv(
        paste0("../../../group/irc/personal/artuurc/Compensation_matrices/", 
               whichComp, ".csv"))
      colnames(comp) <- gsub("\\.{4}.*", "", colnames(comp))
      comp$X <- NULL
      colnames(comp) <- gsub("\\.", "-", colnames(comp))
    } else comp <- ff1@description$SPILL
    colnames(comp) <- gsub("PacB-A", "HV450-A", colnames(comp))
    colnames(comp) <- gsub("PacO-A", "HV500-A", colnames(comp))
    
    ff1 <- compensate(ff1, comp)
    
    #----Transformation 
    ff1 <- transform(ff1, tfList)
    
    #----PeacoQC
    ff_qc <- suppressMessages(
      suppressWarnings(
        PeacoQC::PeacoQC(ff = ff1,
                         channels = colnames(ff)[-which(colnames(ff) == "Time")],
                         plot = TRUE,
                         save_fcs = FALSE,
                         report = TRUE,
                         output_directory = "../../../group/irc/personal/artuurc/Preprocessing")
      )
    )
    selection_quality <- ff_qc$GoodCells
    perc_selection_q <- round(ff_qc$PercentageRemoved, 2)
    ff2 <- ff1[selection_quality, ]
    cat(paste0("PeacoQC removed ", perc_selection_q, "% of cells.\n"))
    selection[selection] <- selection_quality
    
    #----No debris
    fsom <- FlowSOM(ff2,
                    colsToUse = colsToUseDebris,  
                    nClus = 5, 
                    xdim = 9, 
                    ydim = 9, 
                    seed = 2022,
                    scale = FALSE)
    
    debris <- debrisClusterfinder(fsom = fsom,
                                  SSC_cutoff_agg = SSC_cutoff_agg,
                                  FSC_cutoff_agg = FSC_cutoff_agg,
                                  FSC_cutoff_highSSC_agg = FSC_cutoff_highSSC_agg,
                                  title = title)
    debrisList <- if (title %in% names(customDebris)){
      Plot2DScatters(fsom, channelpairs = list(c("FSC-A", "SSC-A")), 
                     clusters = seq_len(NClusters(fsom)), plotFile = paste0("../../../group/irc/personal/artuurc/Preprocessing/2DScatters/", title, ".png"))
      customDebris[[title]]
    } else {
      debris$debrisClusters
    }
    cutoff <- debris$cutoff
    if (cutoff <= FSC_cutoff_agg / 2) cutoff <- FSC_cutoff_agg / 2
    selection_noDebriscl <- !GetClusters(fsom) %in% debrisList
    selection_rightofcutoff <- exprs(ff2)[, "FSC-A"] > cutoff
    selection_noDebris <- selection_noDebriscl & selection_rightofcutoff #TRUE = no debris, FALSE = debris
    perc_selection_debris <- round((sum(!selection_noDebris) / length(selection)) * 100, 2)
    cat(paste0("File has ", round(perc_selection_debris, 2), "% of debris.\n"))
    ff3 <- ff2[selection_noDebris, ]
    selection[selection] <- selection_noDebris
    
    #----No doublets
    factor <- 4
    if (title %in% more_strict_doublets) factor <- 3
    if (title %in% even_more_strict_doublets) factor <- 2
    if (title %in% less_strict_doublets) factor <- 6
    if (title %in% even_less_strict_doublets) factor <- 7
    selection_single <- !seq_len(nrow(ff3)) %in% PeacoQC::RemoveDoublets(
      ff3, output = "full", nmad = factor)$indices_doublets #TRUE = no doublets, FALSE = doublet
    perc_selection_doublets <- round((sum(!selection_single) / length(selection)) * 100, 2)
    cat(paste0("File has ", perc_selection_doublets, "% of doublets.\n"))
    selection[selection] <- selection_single
    ff4 <- ff3[selection_single, ]
    
    #----Total percentage removed
    perc_total <- round((sum(!selection) / length(selection)) * 100, 2)
    
    #----Plot densities
    if (!dir.exists("../../../group/irc/personal/artuurc/Preprocessing/Preprocess_densities")){
      dir.create("../../../group/irc/personal/artuurc/Preprocessing/Preprocess_densities")
    }
    png(paste0("../../../group/irc/personal/artuurc/Preprocessing/Preprocess_densities/",
               title, ".png"), width = 1800, height = 300)
    par(mfrow = c(1, 6))
    plotDens(ff, c("FSC-A", "SSC-A"), main = paste0(title,": All events (",
                                                    nrow(ff), ")" ))
    plotDens(ff, c("FSC-A", "SSC-A"), main = paste0("Margin events (",
                                                    perc_remove_margins,"%)"))
    points(ff@exprs[, c("FSC-A","SSC-A")][!selection_inrange,], col = "red", pch = ".")
    plotDens(ff1, c("Time", "FSC-A"), main = paste0("PeacoQC (",
                                                    perc_selection_q,"%)"))
    points(ff1@exprs[, c("FSC-A","SSC-A")][!selection_quality,], col = "red", pch = ".")
    plotDens(ff2, c("FSC-A", "SSC-A"), main = paste0("Debris (",
                                                     perc_selection_debris, "%)"))
    points(ff2@exprs[, c("FSC-A","SSC-A")][!selection_noDebris, ], col = "red", pch = ".")
    xmax <- max(exprs(ff4)[, "FSC-A"])
    plotDens(ff3, c("FSC-A", "FSC-H"), main = paste0("Doublets (",
                                                     perc_selection_doublets, "%)"),
             xlim = c(0, xmax))
    points(ff3@exprs[, c("FSC-A","FSC-H")][!selection_single,], col = "red", pch = ".")
    plotDens(ff4, c("FSC-A", "SSC-A"), main = paste0("Cleaned (", perc_total, "%)"),
             xlim = c(0, xmax))
    dev.off()
    
    
    #----Print QC_file
    if (!file.exists("../../../group/irc/personal/artuurc/Preprocessing/Quality_control.txt")){
      write("Quality control\n------------------------\n",
            "../../../group/irc/personal/artuurc/Preprocessing/Quality_control.txt",append = T)
      write("Filename\tMargins\tPeacoQC\tDebris\tDoublet\tTotal",
            file = "../../../group/irc/personal/artuurc/Preprocessing/Quality_control.txt",
            append = T)
    }
    write(paste0(ff4@description$FILENAME, "\t",
                 round(perc_remove_margins, 2), "\t",
                 round(perc_selection_q, 2), "\t",
                 round(perc_selection_debris, 2), "\t",
                 round(perc_selection_doublets, 2), "\t",
                 round(perc_total, 2)),
          file = "../../../group/irc/personal/artuurc/Preprocessing/Quality_control.txt",
          append = T)
    #----Save results
    write.FCS(ff4, paste0("../../../group/irc/personal/artuurc/Preprocessing/FCS/", title, "_prepro.fcs"))
    
  }, mc.cores = 24)
}

#----------------------Parse preprocessed files---------------------------------
fcs_files_prepro <- paste0(output_prepro_FCS_files,
                           list.files(output_prepro_FCS_files, pattern = ".*\\.fcs$"))
samples_prepro <- data.frame(File = fcs_files_prepro,
                             Panel = stringr::str_match(fcs_files_prepro, "ALOT|AML")[,1],
                             PatientID = stringr::str_match(fcs_files_prepro, "([0-9]+-?[0-9])")[, 2],
                             Tube = as.numeric(stringr::str_match(fcs_files_prepro, "AML([1-7]*)")[,2]),
                             stringsAsFactors = FALSE) %>% 
  filter(!PatientID %in% patientsToRemove) 
samples_prepro$Tube[is.na(samples_prepro$Tube)] <- 0
tmp <- factor(substr(samples_prepro$PatientID, 1, 2))
levels(tmp) <- c("2015", "2016", "2017", "2018", "2019")
samples_prepro$Year <- tmp

#--------------------------------Normalisation----------------------------------
normalizationList <- list("0" = c("CD45",  "CD7"), 
                          "1" = c("CD45", "FSC-H"), 
                          "2" = c("CD45", "FSC-H"),
                          "3" = c("CD45", "FSC-H"),
                          "4" = c("CD45", "FSC-H"),
                          "5" = c("CD45", "FSC-H"),
                          "6" = c("CD45", "CD123", "FSC-H"))
if (!dir.exists(output_norm_FCS_files)) dir.create(output_norm_FCS_files)
for (tube in seq(0, 6)){
  message(paste0("Tube ", tube))
  samples_tube <- samples_prepro %>% dplyr::filter(Tube == tube)
  channelsToNorm <- GetChannels(read.FCS(samples_tube$File[1]), 
                                normalizationList[[as.character(tube)]])
  ALOT_or_AML <- ifelse(tube == 0, "date_flow_acquisition_ALOT",
                        "date_flow_acquisition_AML_tubes")
  dateColMetadata <- metadata[, c("PatientID", ALOT_or_AML)] %>% 
    dplyr::filter(PatientID %in% samples_tube$PatientID)
  dateColMetadata$markerChange <- factor(
    dateColMetadata[, ALOT_or_AML] > as.Date("2016-05-05"))
  samples_tube <- dplyr::inner_join(
    samples_tube, dateColMetadata, by = "PatientID")
  normModel <- CytoNorm::QuantileNorm.train(files = samples_tube$File, 
                                            labels = samples_tube$markerChange,
                                            channels = channelsToNorm,
                                            transformList = NULL,
                                            nQ = 101,
                                            limit = NULL,
                                            goal = "mean")
  saveRDS(normModel, paste0(output_RDS, "normModel_tube_markerchange_", tube, ".rds"))
  CytoNorm::QuantileNorm.normalize(model = normModel,
                                   files = samples_tube$File,
                                   labels = samples_tube$markerChange,
                                   transformList = NULL,
                                   transformList.reverse = NULL,
                                   outputDir = output_norm_FCS_files,
                                   prefix = "")
  if (tube == 2){
    fcs_files_norm <- paste0(output_norm_FCS_files,
                             list.files(output_norm_FCS_files, pattern = ".*fcs$"))
    samples_tube <- data.frame(File = fcs_files_norm,
                               Panel = stringr::str_match(fcs_files_norm, "ALOT|AML")[,1],
                               PatientID = stringr::str_match(fcs_files_norm, "([0-9]+-?[0-9])")[, 2],
                               Tube = as.numeric(stringr::str_match(fcs_files_norm, "AML([1-7]*)")[,2]),
                               stringsAsFactors = FALSE) %>% filter(Tube == 2)
    dateColMetadata <- metadata[, c("PatientID", "date_flow_acquisition_AML_tubes")] %>% 
      dplyr::filter(PatientID %in% samples_tube$PatientID)
    dateColMetadata$markerChange <- factor(
      dateColMetadata[, ALOT_or_AML] > as.Date("2016-08-01"))
    samples_tube <- dplyr::inner_join(
      samples_tube, dateColMetadata, by = "PatientID")
    
    channelsToNorm <- GetChannels(read.FCS(samples_tube$File[1]), c( "CD64", "CD34", "CD117"))
    normModel <- CytoNorm::QuantileNorm.train(files = samples_tube$File, 
                                              labels = samples_tube$markerChange,
                                              channels = channelsToNorm,
                                              transformList = NULL)
    saveRDS(normModel, paste0(output_RDS, "normModel_tube_batch_tube_2_", tube, ".rds"))
    CytoNorm::QuantileNorm.normalize(model = normModel,
                                     files = samples_tube$File,
                                     labels = samples_tube$markerChange,
                                     transformList = NULL,
                                     transformList.reverse = NULL,
                                     outputDir = output_norm_FCS_files,
                                     prefix = "Norm_",
                                     removeOriginal = TRUE)
  }
  if (tube == 6){
    fcs_files_norm <- paste0(output_norm_FCS_files,
                             list.files(output_norm_FCS_files, pattern = ".*fcs$"))
    samples_norm <- data.frame(File = fcs_files_norm,
                               Panel = stringr::str_match(fcs_files_norm, "ALOT|AML")[,1],
                               PatientID = stringr::str_match(fcs_files_norm, "([0-9]+-?[0-9])")[, 2],
                               Tube = as.numeric(stringr::str_match(fcs_files_norm, "AML([1-7]*)")[,2]),
                               stringsAsFactors = FALSE) %>% filter(Tube == 6)
    tmp <- factor(substr(samples_norm$PatientID, 1, 2))
    levels(tmp) <- c("2015", "2016", "2017", "2018", "2019")
    samples_norm$Year <- tmp == 2019
    channelsToNorm <- GetChannels(read.FCS(samples_norm$File[1]), c("CD123", "CD4")) 
    
    normModel <- CytoNorm::QuantileNorm.train(files = samples_norm$File, 
                                              labels = samples_norm$Year,
                                              channels = channelsToNorm,
                                              transformList = NULL)
    saveRDS(normModel, paste0(output_RDS, "normModel_tube_batch_tube_6_", tube, ".rds"))
    CytoNorm::QuantileNorm.normalize(model = normModel,
                                     files = samples_norm$File,
                                     labels = samples_norm$Year,
                                     transformList = NULL,
                                     transformList.reverse = NULL,
                                     outputDir = output_norm_FCS_files,
                                     prefix = "Norm_",
                                     removeOriginal = TRUE)
  }
}
fcs_files_norm <- paste0(getwd(), output_norm_FCS_files,
                         list.files(output_norm_FCS_files, pattern = ".*fcs$"))
fcs_files_norm_newName <- gsub("prepro", "norm", fcs_files_norm)
fcs_files_norm_newName <- gsub("Norm_", "", fcs_files_norm_newName)
file.rename(fcs_files_norm, fcs_files_norm_newName)

#------------------------Parse normalized files---------------------------------
fcs_files_norm <- paste0(output_norm_FCS_files, 
                         list.files("./Normalized_fcs/", pattern = ".*fcs$"))
samples_norm <- data.frame(File = fcs_files_norm,
                           Panel = stringr::str_match(fcs_files_norm, "ALOT|AML")[,1],
                           PatientID = stringr::str_match(fcs_files_norm, "([0-9]+-?[0-9])")[, 2],
                           Tube = as.numeric(stringr::str_match(fcs_files_norm, "AML([1-7]*)")[,2]),
                           stringsAsFactors = FALSE)  %>% 
  filter(!PatientID %in% patientsToRemove) 
samples_norm$Tube[is.na(samples_norm$Tube)] <- 0
tmp <- factor(substr(samples_norm$PatientID, 1, 2))
levels(tmp) <- c("2015", "2016", "2017", "2018", "2019")
samples_norm$Year <- tmp
samples_norm <- samples_norm %>%  #filter out because double file
  dplyr::filter(File != paste0(output_norm_FCS_files, "1747_ALOT_PBO_norm.fcs"))

#-----------------------------Plot filescatters---------------------------------
if (!dir.exists(paste0(output_norm_FCS_files, "FileScatters"))) {
  dir.create(paste0(output_norm_FCS_files, "FileScatters"))
  dir.create(paste0(output_norm_FCS_files, "FileScatters/Preprocessed"))
  dir.create(paste0(output_norm_FCS_files, "FileScatters/Normalized"))
  dir.create(paste0(output_norm_FCS_files, "FileScatters/Composed"))
}

normalizationList <- list("0" = c("CD45",  "CD7"), 
                          "1" = c("CD45"), 
                          "2" = c("CD45", "CD64", "CD34", "CD117"),
                          "3" = c("CD45"),
                          "4" = c("CD45"),
                          "5" = c("CD45"),
                          "6" = c("CD45", "CD123", "CD4"))


for (tube in seq(0, 6)){
  # Preprocessed
  message(paste0("Tube ", tube))
  cat("Plot preprocessed\n")
  samples_tube <- samples_prepro %>% filter(Tube == tube)
  ALOT_or_AML <- ifelse(tube == 0, "date_flow_acquisition_ALOT",
                        "date_flow_acquisition_AML_tubes")
  dateColMetadata <- metadata[, c("PatientID", ALOT_or_AML)] %>% 
    filter(PatientID %in% samples_tube$PatientID)
  
  dateColMetadata$afterMarkerChange <- factor(dateColMetadata[ , ALOT_or_AML] > 
                                                as.Date("2016-05-05"))
  tmp <- inner_join(samples_tube, dateColMetadata, "PatientID")
  tmp <- tmp[order(tmp[,ALOT_or_AML]), ]

  ff <- read.FCS(tmp$File[1])
  norm_channels <- GetChannels(ff, normalizationList[[as.character(tube)]])
  channelnames <- colnames(ff@exprs)
  channelnames <- channelnames[which(!channelnames %in% c("Time", "FSC-H"))]
  PlotFileScatters(tmp$File, names = tmp$PatientID, groups = tmp$afterMarkerChange, 
                   channels = channelnames,
                   plotFile = paste0(
                     output_norm_FCS_files, 
                     "FileScatters/Preprocessed/FileScatters_preproc_markerChange", 
                     tube, ".png"
                   ),
                   ncol = 1, yLabel = c("marker", "channel"), silent = TRUE)
  
  PlotFileScatters(tmp$File, names = tmp$PatientID, groups = tmp$Year, 
                   channels = channelnames,
                   plotFile = paste0(
                     output_norm_FCS_files, 
                     "FileScatters/Preprocessed/FileScatters_preproc_year",
                     tube, ".png"
                   ),
                   ncol = 1, yLabel = c("marker", "channel"),
                   quantiles = c(0.5), silent = TRUE)
  
  p_pp <- PlotFileScatters(tmp$File, names = tmp$PatientID, groups = tmp$Year, 
                           channels = norm_channels,
                           plotFile = NULL,
                           yLabel = c("marker", "channel"), silent = TRUE)
  
  #Normalized
  if (nrow(samples_norm) == 1) next
  tube_samples <- samples_norm %>% filter(Tube == tube)
  ALOT_or_AML <- ifelse(tube == 0, "date_flow_acquisition_ALOT",
                        "date_flow_acquisition_AML_tubes")
  dateColMetadata <- metadata[, c("PatientID", ALOT_or_AML)] %>% 
    filter(PatientID %in% tube_samples$PatientID)
  
  dateColMetadata$afterMarkerChange <- factor(dateColMetadata[ ,ALOT_or_AML] > 
                                                as.Date("2016-05-05"))
  tmp <- inner_join(tube_samples, dateColMetadata, "PatientID")
  tmp <- tmp[order(tmp[,ALOT_or_AML]), ]

  ff <- read.FCS(tmp$File[1])
  channelnames <- colnames(ff@exprs)
  channelnames <- channelnames[which(!channelnames %in% c("Time", "FSC-H"))]
  cat("Plot normalized\n")
  PlotFileScatters(tmp$File, names = tmp$PatientID, groups = tmp$afterMarkerChange, 
                   channels = channelnames,
                   plotFile = paste0(
                     output_norm_FCS_files,
                     "FileScatters/Normalized/FileScatters_norm_markerChange", 
                     tube, ".png"
                   ),
                   ncol = 1, yLabel = c("marker", "channel"), silent = TRUE)
  
  PlotFileScatters(tmp$File, names = tmp$PatientID, groups = tmp$Year, 
                   channels = channelnames,
                   plotFile = paste0(
                     output_norm_FCS_files, 
                     "FileScatters/Normalized/FileScatters_norm_year", 
                     tube, ".png"
                   ),
                   ncol = 1, yLabel = c("marker", "channel"),
                   quantiles = c(0.5), silent = TRUE)
  
  p_norm <- PlotFileScatters(tmp$File, names = tmp$PatientID, groups = tmp$Year, 
                             channels =  norm_channels,
                             plotFile = NULL,
                             yLabel = c("marker", "channel"), silent = TRUE)

  cat("Plot composed\n")
  png(paste0(output_norm_FCS_files, "FileScatters/Composed/FileScatters_comp_year_", tube, ".png"), 
      width = length(norm_channels)  * (60 + 15 * length(tmp$File)), 
      height = 500)
  p <- ggpubr::ggarrange(plotlist = c(p_pp, p_norm),
                         nrow = 2, ncol = length(norm_channels), 
                         common.legend = TRUE)
  print(p)
  dev.off()
}

#----------------------------Plot PCA of quantiles------------------------------
if (!dir.exists(paste0(output_norm_FCS_files, "PCA_on_quantiles"))) {
  dir.create(paste0(output_norm_FCS_files, "PCA_on_quantiles"))
  dir.create(paste0(output_norm_FCS_files, "PCA_on_quantiles/Preprocessed"))
  dir.create(paste0(output_norm_FCS_files, "PCA_on_quantiles/Normalized"))
}

for (tube in seq(0, 6)){
  message(paste0("Tube ", tube))
  # Preprocessed
  samples_tube <- samples_prepro %>% filter(Tube == tube)
  quantiles <- matrix(NA,
                      nrow = nrow(samples_tube),
                      ncol = 0,
                      dimnames = list(samples_tube$PatientID,
                                      NULL))
  agg <- AggregateFlowFrames(samples_tube$File, cTotal = 3000000, keepOrder = TRUE,
                             silent = TRUE)
  channels <- colnames(agg)[which(!colnames(agg) %in% c(
    "Time", "File", "File_scattered", "Original_ID", "FSC-H"))]
  
  ALOT_or_AML <- ifelse(tube == 0, "date_flow_acquisition_ALOT",
                        "date_flow_acquisition_AML_tubes")
  dateColMetadata <- metadata[, c("PatientID", ALOT_or_AML)] %>% 
    filter(PatientID %in% samples_tube$PatientID)
  
  dateColMetadata$afterMarkerChange <- factor(dateColMetadata[ ,ALOT_or_AML] > 
                                                as.Date("2016-05-05"))
  tmp <- inner_join(samples_tube, dateColMetadata, "PatientID")
  tmp <- tmp[order(tmp[,ALOT_or_AML]), ]
  
  for(channel in channels){
    x <- exprs(agg)[, channel]
    q <- tapply(x, agg@exprs[, "File"],
                quantile, c(0.25, 0.5, 0.75))
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
      width = 500, height = 500)
  p <- ggplot(plotdf) + 
    geom_point(aes(x = PC1, y = PC2, col = afterMarkerChange)) +
    theme_minimal() +
    ggtitle(paste0("Tube ", tube))
  p <- ggExtra::ggMarginal(p, groupFill = TRUE,  alpha = 0.5)
  print(p)
  dev.off()
  
  png(paste0(output_norm_FCS_files, "PCA_on_quantiles/Preprocessed/PCA_q_prepro_year_", tube, ".png"), 
      width = 500, height = 500)
  p <- ggplot(plotdf) + 
    geom_point(aes(x = PC1, y = PC2, col = Year)) +
    theme_minimal() +
    ggtitle(paste0("Tube ", tube))
  p <- ggExtra::ggMarginal(p, groupFill = TRUE,  alpha = 0.5)
  print(p)
  dev.off()
  
  png(paste0(output_norm_FCS_files, "PCA_on_quantiles/Preprocessed/PCA_prepro_var_", tube, ".png"),
      width = 500, height = 500)
  p <- factoextra::fviz_pca_var(pca,
                                col.var = "contrib",
                                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                repel = TRUE)
  print(p)
  dev.off()
  
  # Normalized
  samples_tube <- samples_norm %>% filter(Tube == tube)
  quantiles <- matrix(NA,
                      nrow = nrow(samples_tube),
                      ncol = 0,
                      dimnames = list(samples_tube$PatientID,
                                      NULL))
  agg <- AggregateFlowFrames(samples_tube$File, cTotal = 3000000, keepOrder = TRUE,
                             silent = TRUE)
  channels <- colnames(agg)[which(!colnames(agg) %in% c(
    "Time", "File", "File_scattered", "Original_ID", "FSC-H"))]
  
  ALOT_or_AML <- ifelse(tube == 0, "date_flow_acquisition_ALOT",
                        "date_flow_acquisition_AML_tubes")
  dateColMetadata <- metadata[, c("PatientID", ALOT_or_AML)] %>% 
    filter(PatientID %in% samples_tube$PatientID)
  
  dateColMetadata$afterMarkerChange <- factor(dateColMetadata[ ,ALOT_or_AML] > 
                                                as.Date("2016-05-05"))
  tmp <- inner_join(samples_tube, dateColMetadata, "PatientID")
  tmp <- tmp[order(tmp[,ALOT_or_AML]), ]
  
  for(channel in channels){
    x <- exprs(agg)[, channel]
    q <- tapply(x, agg@exprs[, "File"],
                quantile, c(0.25, 0.5, 0.75))
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
      width = 500, height = 500)
  p <- ggplot(plotdf) + 
    geom_point(aes(x = PC1, y = PC2, col = afterMarkerChange)) +
    theme_minimal() +
    ggtitle(paste0("Tube ", tube))
  p <- ggExtra::ggMarginal(p, groupFill = TRUE,  alpha = 0.5)
  print(p)
  dev.off()
  
  png(paste0(output_norm_FCS_files, "PCA_on_quantiles/Normalized/PCA_q_norm_year_", tube, ".png"),
      width = 500, height = 500)
  p <- ggplot(plotdf) + 
    geom_point(aes(x = PC1, y = PC2, col = Year)) +
    theme_minimal() +
    ggtitle(paste0("Tube ", tube))
  p <- ggExtra::ggMarginal(p, groupFill = TRUE,  alpha = 0.5)
  print(p)
  dev.off()
  
  png(paste0(output_norm_FCS_files, "PCA_on_quantiles/Normalized/PCA_norm_var_", tube, ".png"),
      width = 500, height = 500)
  p <- factoextra::fviz_pca_var(pca,
                                col.var = "contrib",
                                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                repel = TRUE)
  print(p)
  dev.off()
}

he <- lapply(samples_tube$File, function(x){
  ff = read.FCS(x)
  return(GetMarkers(ff, colnames(exprs(ff))))
})
