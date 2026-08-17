## This file prepares for Rainbow analysis. Two files are required: an XLSX-file with the cumulative values per well, and a CSV-file with the clustered values for each partition, i.e. the combinations +++++, +--+-, etc

# Please install the MultiplexPCRAnalyser package v1.0.0 to use with this file. This can be found here: https://github.com/buchauer-lab/Q4ddPCR_Analysis
# This file is used to set plate and sample specific parameters for the analysis of Q4ddPCR

# version 1.0.0 (date 2026-08-14) 

# PART A: Creation of the XLSX-file
# PART  B: Creation of the CSV-file

#extras to accomodate Q4ddPCR script from Gaebler
# PART C: Functions that didn't work in the original script and have been updated

# run this script together with Q4ddPCR_Qiacuity_RainbowAnalysis.R


###########


library(dplyr)
library(tidyr)
library(readr)
library(writexl)
library(MultiplexPCRAnalyser)


## PART A:  Creating a Bio-rad format file from the Qiacuity analysis file (with composite data in wells)

createXLSXfile <- function(qiacuity_file, clusterRPP30_file, path = getwd(), xlsx_outfile=paste(substr(qiacuity_file, 1, nchar(qiacuity_file)-4),Sys.Date(),"_outfile.xlsx", sep=""), verbose_option = TRUE){
  
  #in_csv <- read_csv(qiacuity_file)
  in_csv <- read_csv(qiacuity_file, skip=1)  # analysis file has "," as separator
  in_clusterRPP30 <- read_csv(clusterRPP30_file,  skip=1)  # MO file has ";" as separator
  #in_clusterRPP30 <- read_csv(clusterRPP30_file)  # MO file has ";" as separator
  
  
  ## 1. Desired Bio-Rad column order
  excel_names <- c(
    "Well", "Sample description 1", "Sample description 2",
    "Sample description 3", "Sample description 4",
    "Target", "Conc(copies/µL)", "pg/µL",
    "Status", "Status Reason", "Experiment", "SampleType", "TargetType",
    "Supermix", "DyeName(s)", "Copies/20µLWell",
    "TotalConfMax", "TotalConfMin", "PoissonConfMax", "PoissonConfMin",
    "Accepted Droplets", "Positives", "Negatives",
    "Copies/uL linked molecules", "CNV", "TotalCNVMax", "TotalCNVMin",
    "PoissonCNVMax", "PoissonCNVMin", "ReferenceCopies", "UnknownCopies",
    "Threshold1", "Threshold2", "Threshold3",
    "ThresholdSigmaAbove", "ThresholdSigmaBelow",
    "ReferenceUsed", "Ratio", "TotalRatioMax", "TotalRatioMin",
    "PoissonRatioMax", "PoissonRatioMin",
    "Fractional Abundance",
    "TotalFractionalAbundanceMax", "TotalFractionalAbundanceMin",
    "PoissonFractionalAbundanceMax", "PoissonFractionalAbundanceMin",
    "MeanAmplitudeOfPositives", "MeanAmplitudeOfNegatives",
    "MeanAmplitudeTotal", "ExperimentComments", "MergedWells",
    "TotalConfidenceMax68", "TotalConfidenceMin68",
    "PoissonConfidenceMax68", "PoissonConfidenceMin68",
    "TotalCNVMax68", "TotalCNVMin68",
    "PoissonCNVMax68", "PoissonCNVMin68",
    "TotalRatioMax68", "TotalRatioMin68",
    "PoissonRatioMax68", "PoissonRatioMin68",
    "TotalFractionalAbundanceMax68", "TotalFractionalAbundanceMin68",
    "PoissonFractionalAbundanceMax68", "PoissonFractionalAbundanceMin68",
    "TiltCorrected",
    "Ch1+Ch2+", "Ch1+Ch2-", "Ch1-Ch2+", "Ch1-Ch2-",
    "Ch3+Ch4+", "Ch3+Ch4-", "Ch3-Ch4+", "Ch3-Ch4-",
    "Ch5+Ch6+", "Ch5+Ch6-", "Ch5-Ch6+", "Ch5-Ch6-"
  )
  
  ## 2. Mapping from target to dye name (adapt to your panel)
  dye_map <- c(
    PSI = "FAM",
    ENV = "HEX",
    POL = "TAMRA",
    RU5 = "ROX",
    GAG = "Cy5",
    RPP30 = "FAM",
    `RPP30 3p`= "HEX"
  )
  
  ## 3. Build a full Bio-Rad-style table from your current per-target data
  ## Replace 'per_target' and column names here with your real object
  
  biorad_like <- in_csv %>%
    mutate(
      `Sample description 1` = `Sample/NTC/Control`,
      `Sample description 2` = `Total volume [µl]`,
      `Sample description 3` = NA_character_,
      `Sample description 4` = NA_character_,
      `Target` = `Target (Name)`,
      `Conc(copies/µL)`      = `Conc. [cp/<U+00B5>L] (dPCR reaction)`,
      `pg/µL`                = NA_real_,
      Status                 = "Manual",
      `Status Reason`        = "Converted from Qiacuity",
      Experiment             = "Qiacuity_Rainbow",
      SampleType             = NA_real_,
      TargetType             = NA_real_,
      Supermix               = NA_real_,
      `DyeName(s)`           = unname(dye_map[`Target (Name)`]),
      `Copies/20µLWell`      = `Conc. [cp/<U+00B5>L] (undiluted sample)`,
      TotalConfMax           = NA_real_,
      TotalConfMin           = NA_real_,
      PoissonConfMax         = NA_real_,
      PoissonConfMin         = NA_real_,
      `Accepted Droplets`    =`Partitions (Valid)`,
      `Positives`              = `Partitions (Positive)`,
      `Negatives`             = `Partitions (Negative)`,
      `Copies/uL linked molecules` = NA_real_,
      CNV                    = NA_real_,
      TotalCNVMax            = NA_real_,
      TotalCNVMin            = NA_real_,
      PoissonCNVMax          = NA_real_,
      PoissonCNVMin          = NA_real_,
      ReferenceCopies        = NA_real_,
      UnknownCopies          = NA_real_,
      Threshold1             = NA_real_,
      Threshold2             = NA_real_,
      Threshold3             = NA_real_,
      ThresholdSigmaAbove    = NA_real_,
      ThresholdSigmaBelow    = NA_real_,
      ReferenceUsed          = NA_character_,
      Ratio                  = NA_real_,
      TotalRatioMax          = NA_real_,
      TotalRatioMin          = NA_real_,
      PoissonRatioMax        = NA_real_,
      PoissonRatioMin        = NA_real_,
      `Fractional Abundance`           = NA_real_,
      TotalFractionalAbundanceMax     = NA_real_,
      TotalFractionalAbundanceMin     = NA_real_,
      PoissonFractionalAbundanceMax   = NA_real_,
      PoissonFractionalAbundanceMin   = NA_real_,
      MeanAmplitudeOfPositives        = NA_real_,
      MeanAmplitudeOfNegatives        = NA_real_,
      MeanAmplitudeTotal              = NA_real_,
      ExperimentComments              = NA_character_,
      MergedWells                     = NA_character_,
      TotalConfidenceMax68            = NA_real_,
      TotalConfidenceMin68            = NA_real_,
      PoissonConfidenceMax68          = NA_real_,
      PoissonConfidenceMin68          = NA_real_,
      TotalCNVMax68                   = NA_real_,
      TotalCNVMin68                   = NA_real_,
      PoissonCNVMax68                 = NA_real_,
      PoissonCNVMin68                 = NA_real_,
      TotalRatioMax68                 = NA_real_,
      TotalRatioMin68                 = NA_real_,
      PoissonRatioMax68               = NA_real_,
      PoissonRatioMin68               = NA_real_,
      TotalFractionalAbundanceMax68   = NA_real_,
      TotalFractionalAbundanceMin68   = NA_real_,
      PoissonFractionalAbundanceMax68 = NA_real_,
      PoissonFractionalAbundanceMin68 = NA_real_,
      TiltCorrected                   = "No",
      `Ch1+Ch2+`                      = NA_real_,
      `Ch1+Ch2-`                      = NA_real_,
      `Ch1-Ch2+`                      = NA_real_,
      `Ch1-Ch2-`                      = NA_real_,
      `Ch3+Ch4+`                      = NA_real_,
      `Ch3+Ch4-`                      = NA_real_,
      `Ch3-Ch4+`                      = NA_real_,
      `Ch3-Ch4-`                      = NA_real_,
      `Ch5+Ch6+`                      = NA_real_,
      `Ch5+Ch6-`                      = NA_real_,
      `Ch5-Ch6+`                      = NA_real_,
      `Ch5-Ch6-`                      = NA_real_
      
      
    )%>%
    # Keep only needed columns and order them as in Bio-Rad export
    # Keep only needed columns and order them as in Bio-Rad export
    select(all_of(excel_names))
  
  
  
  #### Add data for RPP30 counts in Ch1 and Ch2 to account for shearing
  # The `Group` data (++,+-,...) for `Count categories` should go to the xlsx-sheet, in "Ch1+Ch2+", "Ch1+Ch2-", "Ch1-Ch2+", "Ch1-Ch2-"
  
  ## 4. Read Qiacuity RPP30 shearing control
  
  rpp30_counts <- in_clusterRPP30 %>%
    rename(
      well    = Well,
      grp     = Group,
      count   = `Count categories`
    ) %>%
    mutate(
      well  = as.character(well),
      grp   = as.character(grp),
      count = as.numeric(count)
    ) %>%
    filter(!is.na(grp)) %>%
    group_by(well, grp) %>%
    summarise(
      Count = sum(count, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  ## 5. Spread into quadrants: ++, +-, -+, --
  rpp30_wide <- rpp30_counts %>%
    tidyr::pivot_wider(
      id_cols     = well,
      names_from  = grp,
      values_from = Count,
      values_fill = list(Count = 0)
    ) %>%
    mutate(
      Ch1pCh2p = if ("++" %in% names(.)) `++` else 0,
      Ch1pCh2m = if ("+-" %in% names(.)) `+-` else 0,
      Ch1mCh2p = if ("-+" %in% names(.)) `-+` else 0,
      Ch1mCh2m = if ("--" %in% names(.)) `--` else 0
    ) %>%
    select(
      Well      = well,
      Ch1pCh2p,
      Ch1pCh2m,
      Ch1mCh2p,
      Ch1mCh2m
    )
  
  ## 6. Join into existing Bio-Rad table and overwrite the 4 gate columns
  
  biorad_like_with_rpp30 <- biorad_like %>%
    left_join(rpp30_wide, by = "Well") %>%
    mutate(
      `Ch1+Ch2+` = ifelse(!is.na(Ch1pCh2p), Ch1pCh2p, `Ch1+Ch2+`),
      `Ch1+Ch2-` = ifelse(!is.na(Ch1pCh2m), Ch1pCh2m, `Ch1+Ch2-`),
      `Ch1-Ch2+` = ifelse(!is.na(Ch1mCh2p), Ch1mCh2p, `Ch1-Ch2+`),
      `Ch1-Ch2-` = ifelse(!is.na(Ch1mCh2m), Ch1mCh2m, `Ch1-Ch2-`)
    ) %>%
    select(all_of(excel_names))
  
  
  writexl::write_xlsx(biorad_like_with_rpp30, xlsx_outfile)
  if (verbose_option) 
    cat("XLSX-file: ", xlsx_outfile, " written.\n")
}

# PART B: Create the clusterfiles: 

createCSVfile <- function(clusterRainbow_file, clusterRPP30_file, path = getwd(), csv_outfile=paste(substr(qiacuity_file, 1, nchar(qiacuity_file)-4),Sys.Date(),"_outfile.csv", sep=""), verbose_option = TRUE){
  
  in_clusterRPP30 <- read_csv(clusterRPP30_file,  skip=1)
  in_clusterRainbow <- read_csv(clusterRainbow_file, skip=1)
  
  hiv_raw<- bind_rows(
    in_clusterRPP30 %>% mutate(panel = "RPP30"),
    in_clusterRainbow %>% mutate(panel = "Rainbow")
  )
  
  hiv <- hiv_raw %>%
    rename(
      well    = Well,
      grp     = Group,                 # rename to a simple name
      count   = `Count categories`,
      targets = `Target names`
    ) %>%
    mutate(
      well  = as.character(well),
      grp   = as.character(grp),
      count = as.numeric(count)
    ) %>%
    filter(!is.na(grp))
  
  hiv_groups <- hiv %>%
    group_by(well, grp) %>%
    summarise(
      Count = sum(count, na.rm = TRUE),
      .groups = "drop"
    )
  
  decode_group <- function(g) {
    data.frame(
      Target1 = as.integer(substr(g, 1, 1) == "+"),
      Target2 = as.integer(substr(g, 2, 2) == "+"),
      Target3 = as.integer(substr(g, 3, 3) == "+"),
      Target4 = as.integer(substr(g, 4, 4) == "+"),
      Target5 = as.integer(substr(g, 5, 5) == "+")
    )
  }
  
  cluster_like <- cbind(
    hiv_groups,
    decode_group(hiv_groups$grp)
  ) %>%
    transmute(
      Well       = well,
      `Target 1` = Target1,
      `Target 2` = Target2,
      `Target 3` = Target3,
      `Target 4` = Target4,
      `Target 5` = Target5,
      Count      = Count
    )
  write_csv(data.frame(c("Target value of 0 = negative","Target value of 1 = positive", "Target value of u = unclassified(Advanced Classification Method)"," ")), csv_outfile, col_names=F) ## This just messed up for later read_csv in Gaebler script
  write_csv(cluster_like, csv_outfile, append=T, col_names=T)
  #write_csv(data.frame(c(" ", " ","Well,Cluster 1,Cluster 2,Angle,S Value"," "," ")), csv_outfile, append=T, col_names=F, quote="none") #extrapart neeted
  
  if (verbose_option) 
    cat("CSV-file: ", csv_outfile, " written.\n")
}

# PART C:
# This superseeds the create_confusion_matrix function as it is only for 4 colors.

create_confusion_matrix2 <- function (df1, data_table, ch_dye, target_channel) 
{
  target_dye_map <- setNames(ch_dye[target_channel], names(target_channel))
  dye_gene_map <- tibble::deframe(unique(data_table[, c("DyeName(s)", 
                                                        "Target")]))
  dye_gene_map <- unlist(lapply(dye_gene_map, function(x) {
    paste0(substr(x, nchar(x) - 2, nchar(x)), "+")
  }))
  idx <- grep("Target", names(df1))
  names(df1)[idx] <- dye_gene_map[target_dye_map[names(df1)[idx]]]
  df_long <- df1[, 1:(length(ch_dye)+2)] %>% rowwise() %>% mutate(name = {    #changed 2 to 1
    cols <- names(across(2:(length(ch_dye)+1)))[c_across(2:(length(ch_dye)+1)) == 1]
    if (length(cols) == 0) 
      "all_negative"
    else paste(cols, collapse = "")
  }) %>% ungroup()
  conf_mat <- df_long %>% tidyr::pivot_wider(id_cols = Well, 
                                             names_from = name, values_from = Count, values_fill = 0)
  return(conf_mat)
}

get_multi_pos2 <- function (df1, genes, tar_mio_factor) 
{
  if (is.null(genes)) {
    stop("Need to specify genes to compute multiple positives.")
  }
  multi_pos <- get_multiplet_count(df, genes) ## can't access
  conc_pos <- df$`Conc(copies/uL)` * as.vector(unlist(multi_pos))/df$Positives
  name <- paste0("Concentration ", paste(genes, collapse = "."), 
                 " positive for target (copies/ul)")
  if (length(conc_pos) != 0) {
    df[[name]] <- conc_pos
  }
  else {
    df[[name]] <- 0
  }
  rows_w_target <- Reduce("|", lapply(genes, function(x) (grepl(x, 
                                                                df$Target))))
  df[is.na(df[, name]), name] <- 0
  df[!(rows_w_target), name] <- NA
  df <- df %>% group_by(group_id, Target) %>% mutate(`:=`(!!sym(paste0("Mean ", 
                                                                       name)), mean(!!sym(name), na.rm = TRUE))) %>% mutate(`:=`(!!sym(paste0("SD ", 
                                                                                                                                              name)), sd(!!sym(name), na.rm = TRUE))) %>% ungroup()
  name2 <- paste0("Intact concentration ", paste(genes, collapse = "."), 
                  " (copies/ul)")
  df <- df %>% group_by(group_id) %>% mutate(`:=`(!!sym(name2), 
                                                  mean(!!sym(name), na.rm = TRUE))) %>% mutate(`:=`(!!sym(paste0("SD ", 
                                                                                                                 name2)), sd(!!sym(name), na.rm = TRUE))) %>% ungroup()
  name3 <- paste0("intact provirus/Mio cells ", paste(genes, 
                                                      collapse = "."))
  df <- df %>% mutate(`:=`(!!sym(name3), tar_mio_factor[`Sample description 1`] * 
                             10^6 * !!sym(name2)/(`Mean concentration RPP30 + RPP30Shear (copies/uL)`))) %>% 
    mutate(`:=`(!!sym(paste0("SD ", name3)), abs(tar_mio_factor[`Sample description 1`] * 
                                                   10^6/`Mean concentration RPP30 + RPP30Shear (copies/uL)`) * 
                  !!sym(paste0("SD ", name2))))
  df[!(rows_w_target), c(name2, name3)] <- NA
  df[[paste0(name3, ", corrected for shearing")]] <- df[[name3]]/df$`Mean unsheared`
  df[[paste0("SD ", name3, ", corrected for shearing")]] <- df[[paste0("SD ", 
                                                                       name3)]]/df$`Mean unsheared`
  return(df)
}

merge_tables2 <- function (data_table, conf_mat, shear_table) 
{
  tab <- merge(data_table, conf_mat, by = "Well", all = T)
  tab <- Filter(function(x) !all(is.na(x)), tab)
  if (any(water_name %in% unique(tab$`Sample description 1`))) {
    return(tab)
  }
  if (all(unique(tab$`Sample description 1`) %in% shear_table$`Sample description 1`)) {
    mean_rpp_conc <- tibble::deframe(unique(shear_table[, 
                                                        c("Sample description 1", "Mean concentration RPP30 + RPP30Shear (copies/uL)")]))
    mean_unsheared <- tibble::deframe(unique(shear_table[, 
                                                         c("Sample description 1", "Mean unsheared")]))
    tab$`Mean concentration RPP30 + RPP30Shear (copies/uL)` <- mean_rpp_conc[tab$`Sample description 1`]
    tab$`Mean unsheared` <- mean_unsheared[tab$`Sample description 1`]
  }
  else {
    stop("All sample descriptions in data_table must be present in shear_table.")
  }
  return(tab)
}

add_dilution_factor2 <- function (df1, dilution_factor) 
{
  if (any(!(df1$`Sample description 1` %in% names(dilution_factor)))) {
    warning("Assuming 1 as dilution factor if not specified otherwise.")
    if (any(!(names(dilution_factor) %in% df1$`Sample description 1`))) {
      warning("Not all specified dilution factor names appear in sample description.\n                Please, check if this is expected.")
    }
    dilution_factor[unique(df1$`Sample description 1`[!(df1$`Sample description 1` %in% 
                                                          names(dilution_factor))])] <- 1
  }
  df1$`Dilution Factor` <- dilution_factor[df1$`Sample description 1`]
  return(df1)
}

compute_total_HIV_envPsi2 <-function (df)  #not implemented. Don't see the need
{
  if (!("Mean Target/Mio cells" %in% names(df))) {
    stop("No Mean Target/Mio cells detected")
  }
  env <- unique(df[grepl("ENV", df$Target), "Mean Target/Mio cells"])
  psi <- unique(df[grepl("PSI", df$Target), "Mean Target/Mio cells"])
  name <- grep("^intact provirus/Mio cells(?=.*PSI)(?=.*ENV)(?!.*POL)(?!.*GAG)(?=.*shearing)", 
               names(df), value = TRUE, perl = TRUE)
  if (length(name) == 0) {
    warning("No Env+Psi+ detected. Will use 0.")
    envpsi <- 0
  }
  else {
    envpsi <- unique(na.omit(df[[name]]))
  }
  return(unlist(unname(env + psi - envpsi)))
}







