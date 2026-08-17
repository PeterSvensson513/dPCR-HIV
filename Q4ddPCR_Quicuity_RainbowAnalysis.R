# Modifications of the Gaebler script to execute dPCR analysis of intact/defective HIV, run together with the Q4ddPCR_convertQiacuityFiles.R for updated functions
# version 1.0.0 (date 2026-08-14) 


# Step 1: Set Parameters 
xlsx_file=paste(substr(qiacuity_file, 1, nchar(qiacuity_file)-4),Sys.Date(),"_outfile.xlsx", sep="")
csv_file=paste(substr(qiacuity_file, 1, nchar(qiacuity_file)-4),Sys.Date(),"_outfile.csv", sep="")

createXLSXfile(qiacuity_file, clusterRPP30_file)
createCSVfile(clusterRainbow_file, clusterRPP30_file)

# Step 2: Read files
information <- read_files(
  xlsx_file,
  csv_file,
  csv_skip,
  remove_channel,
  rm_zero_channel_wells
)


# Read input files
#in_csv <- information[[1]]
in_csv <- read_csv(csv_file, skip=csv_skip)
dtQC   <- information[[2]]

# Quality control and preprocessing



s_desc <- unique(dtQC$`Sample description 1`)
tar_mio_factor <- setNames(rep(rpp30c_per_cell, length(s_desc)), s_desc)
dilution_factor <- setNames(rep(dilution_factor_RPP30, length(s_desc)), s_desc)
#special solution for dilution experiments 2026
if (is_dilution_curve==TRUE)
  dilution_factor[1:5] <- c(455,112,28,7,1)
mean_cells_per_reac_factor <-setNames(rep(replicates_per_sample, length(s_desc)), s_desc)

dtQC <- sufficient_droplets(dtQC, threshold)
dtQC <- add_dilution_factor2(dtQC, dilution_factor)
if (!all(unique(dtQC$`Sample description 1`) %in% names(tar_mio_factor))) {
  warning("Set tar_mio_factor to 1")
  s_desc <- unique(dtQC$`Sample description 1`)
  tar_mio_factor <- setNames(rep(1, length(s_desc)), s_desc)
}

# Split wells and tables
shear_wells <- unique(dtQC[dtQC$Target %in% shear_name, "Well"])
water_wells <- unique(dtQC[dtQC$`Sample description 1` %in% water_name, "Well"]) # change 20260501 the description has the water name
data_wells  <- setdiff(unique(dtQC$Well), union(shear_wells, water_wells))

group_ids <- get_group_id(dtQC)
dtQC$group_id <- group_ids[dtQC$Well]
shear_table <- dtQC[dtQC$Well %in% shear_wells, ]
shear_table$`Sample description 1`<-unlist(strsplit(shear_table$`Sample description 1`, " RPP30"))
water_table <- dtQC[dtQC$Well %in% water_wells, ]
water_csv   <- in_csv[in_csv$Well %in% water_wells, ]

data_table  <- dtQC[dtQC$Well %in% data_wells, ]
data_csv    <- in_csv[in_csv$Well %in% data_wells, ]

data_table <- data_table[, grep("Ch", colnames(data_table), value = TRUE, invert = TRUE)]
data_table <- data_table |> filter(!Well %in% remove_channel)
shear_table <- shear_table |> filter(!Well %in% remove_channel)
data_table <- data_table |> filter(!Well %in% water_wells)
shear_table <- shear_table |> filter(!Well %in% water_wells)

# Shearing factor computation
shear_table <- compute_shearing_factor(
  shear_table,
  mean_copies_factor,
  mean_cells_per_reac_factor
)

#special solution for dilution experiments 2026
if (is_dilution_curve==TRUE){
  # if (sum(!unique(data_table[,2])%in%unique(unlist(shear_table[,2])))>0)
  #    missingdata<-unique(data_table[,2])[!unique(data_table[,2])%in%unique(unlist(shear_table[,2]))]
  
  shear_table2<-shear_table[c(1,2,1,2,1,2,1,2,1,2,3,4),]
  shear_table2[3:4,1]<-"I1"
  shear_table2[5:6,1]<-"I2"
  shear_table2[7:8,1]<-"I3"
  shear_table2[9:10,1]<-"I4"
  shear_table2[3:10,2]<-c("D2","D2","D3","D3", "D4","D4", "D5","D5")
  shear_table<-shear_table2
  
}

# Main table computations
conf_mat <- create_confusion_matrix2(
  data_csv,
  data_table,
  ch_dye,
  target_channel
)
tab <- merge_tables2(data_table, conf_mat, shear_table)
tab <- compute_target_means(tab)

# Multi-positive analysis
# this computes the possible combinations of target genes (doublets, triplets)
multi_positives <- get_multipos(compute_all_positives_for)
for (multip in multi_positives) {
  tab <- get_multi_pos(tab, multip, tar_mio_factor)
}

# Total HIV quantification
total_HIV_dict <- setNames(
  unlist(lapply(unique(tab$group_id), function(x) {
    compute_total_HIV(tab[tab$group_id == x, ])
  })),
  unique(tab$group_id)
)

tab[["total HIV DNA/Mio cells"]] <-
  total_HIV_dict[as.character(tab$group_id)]


tab_intact <- tab |> 
  group_by(group_id, `Sample description 1`) |>
  select(c(`Accepted Droplets`, Threshold, `Dilution Factor`,starts_with("intact provirus/Mio cells")&ends_with("shearing"))) |>
  summarize(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop_last")	

colnames(tab_intact)[-(1:5)]<- sapply(colnames(tab_intact)[-(1:5)], function(x) substr(x,27, nchar(x)-24))

output_file <- paste(substr(xlsx_file,1,nchar(xlsx_file)-12),"_Analysis.xlsx", sep="")
#openxlsx::write.xlsx(tab_intact, output_file)


# Water control analysis									
if (length(water_wells) > 0) {
  h2o_conf_mat <- create_confusion_matrix2(
    water_csv,
    water_table,
    ch_dye,
    target_channel
  )
  
  h2o_tab <- merge_tables(water_table, h2o_conf_mat, shear_table)
  h2o_tab <- h2o_tab[, 1:(ncol(h2o_tab) - 2)]
} else {
  h2o_tab <- NULL
}

# Output generation
out_tables <- lapply(unique(tab$group_id), function(x) {
  tab[
    tab$group_id == x,
    grep("group_id", names(tab), value = TRUE, invert = TRUE)
  ]
})
names(out_tables)<-unlist(tab_intact[,2])
out_tables <- append(out_tables,list(intact=tab_intact))
out_tables <- append(out_tables,list(h2o=h2o_tab))
out_tables <- append(out_tables,list(shear=shear_table))
out_tables <- append(out_tables,list(cmatrix=conf_mat))



## Add a summary sheet similar to results from Processing_Q4dPCR.R. 

tab_summary <- tab |> 
  group_by(group_id, `Sample description 1`) |>
  select(c(`intact provirus/Mio cells PSI.ENV.POL.RU5.GAG, corrected for shearing`,`total HIV DNA/Mio cells`, `intact provirus/Mio cells PSI.ENV, corrected for shearing`,`intact provirus/Mio cells PSI.ENV` ,`Mean unsheared`))  |>
             summarize(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop_last")	
           

tmp<-tab |> 
  group_by(group_id, `Sample description 1`, `Target`="RU5") |>
  select(c(`Mean Target/Mio cells`))  |>
  summarize(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop_last")
tab_summary$RU5 <- tmp$`Mean Target/Mio cells`


tmp<-tab |> 
  group_by(group_id, `Sample description 1`, `Target`="PSI") |>
  select(c(`Mean Target/Mio cells`))  |>
  summarize(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop_last")
tab_summary$PSI <- tmp$`Mean Target/Mio cells`

tmp<-tab |> 
  group_by(group_id, `Sample description 1`, `Target`="ENV") |>
  select(c(`Mean Target/Mio cells`))  |>
  summarize(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop_last")
tab_summary$ENV <- tmp$`Mean Target/Mio cells`


tab_summary$defective <- tab_summary$`total HIV DNA/Mio cells`-tab_summary$`intact provirus/Mio cells PSI.ENV.POL.RU5.GAG, corrected for shearing`
tab_summary$intact_fraction <- tab_summary$`intact provirus/Mio cells PSI.ENV.POL.RU5.GAG, corrected for shearing`/tab_summary$`total HIV DNA/Mio cells`*100
tab_summary$intact_fraction_RU5 <- tab_summary$`intact provirus/Mio cells PSI.ENV.POL.RU5.GAG, corrected for shearing`/tab_summary$RU5*100
tab_summary$IPDA_total <- tab_summary$PSI+tab_summary$ENV-tab_summary$`intact provirus/Mio cells PSI.ENV, corrected for shearing`
tab_summary$IPDA_defective <- tab_summary$IPDA_total-tab_summary$`intact provirus/Mio cells PSI.ENV, corrected for shearing`
tab_summary$IPDA_intactfraction <- tab_summary$`intact provirus/Mio cells PSI.ENV, corrected for shearing`/tab_summary$IPDA_total*100
tab_summary$defective_RU5 <- tab_summary$RU5-tab_summary$`intact provirus/Mio cells PSI.ENV.POL.RU5.GAG, corrected for shearing`

#tab_final<-tab_summary[,c(2,3,4,11,12,5,14,15,16,6,13),]
#colnames(tab_final) <- c("Participant","Intact (5-color) per E6 CD4","Total HIV DNA per E6 CD4",	"Defectives per E6 CD4",	"Intact fraction (%)","IPDA Intacts per E6 CD4","Total HIV DNA IPDA per E6 CD4","IPDA Defectives per E6 CD4","Intact fraction IPDA (%)", "shearing","Intact fraction 5 colors (%RU5)")

tab_final<-tab_summary[,c(2,3,8,17,13,5,14,15,16,7,4,12),]
colnames(tab_final) <- c("Sample","Intact (5-color) per E6 CD4","Total HIV DNA per E6 CD4 (RU5)",	"Defectives per E6 CD4",	"Intact fraction (%)","IPDA Intacts per E6 CD4","Total HIV DNA IPDA per E6 CD4","IPDA Defectives per E6 CD4","Intact fraction IPDA (%)", "shearing","Total HIV DNA per E6 CD4","Intact fraction 5 colors (%)")

out_tables2<-append(list(summary=tab_final),out_tables)
openxlsx::write.xlsx(out_tables2, output_file)
