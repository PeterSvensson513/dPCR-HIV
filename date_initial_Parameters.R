# set file names
# file path to the .csv file that you exported under 'Data table':
qiacuity_file = "analysis.csv"

# file names to the 2 .csv files with data from individual partitions that you exported by clicking 'Export Cluster Data' for Rainbow and RPP30:
clusterRPP30_file = "MO_RPP30.txt"       
clusterRainbow_file = "MO_Rainbow.csv" 

path=""
setwd(path)

Rscript_file <-"https://github.com/PeterSvensson513/dPCR-HIV/blob/main/Q4ddPCR_Quicuity_RainbowAnalysis.R"
Rfunctions_file <-"https://github.com/PeterSvensson513/dPCR-HIV/blob/main/Q4ddPCR_convertQiacuityFiles.R"

is_dilution_curve=TRUE

# add the dilution factor if you used a lower DNA concentration for the RPP30-wells than for the HIV-reaction wells
custom_dilution_factor <- TRUE
dilution_factor_RPP30=10

# reflects number of RPP30 copies per cell. 2 for primary cells, 4 for JLat10.6
rpp30c_per_cell=4

# volume per well 12 for 8.5K plate, 12 for 26K plate
mean_copies_factor <- 12 

# number of replicates 
replicates_per_sample=2

# minimum number of accepted droplets to continue with well
threshold <- 4000 

# names used in shear controls and water controls
shear_name <- c("RPP30", "RPP30 3p")
water_name <- c("WATER")

# add the wells you want to exclude from your analysis, add all the wells with positive, negative and no template controls. Using formula:
# remove_channel <- c("H1","H2","H3") 
remove_channel = NULL

# if you have used different dilution factors, number of replicates or RPP30 copies you can change that. But you also need to change the script then.
#dilution_factor <- c("Sample 1" = 10,"Sample 2" = 10) 




##### no need to change below


# remove wells that have concentration 0 for at least one channel
rm_zero_channel_wells <- FALSE 

csv_skip <- 4 # number of rows before the table starts

# define which dye was used in which channel
ch_dye <- c("Ch1" = "FAM",
            "Ch2" = "HEX",
            "Ch3" = "TAMRA",
            "Ch4" = "ROX",
            "Ch5" = "Cy5")

# define which Target column (in the csv file) corresponds to which channel
target_channel <- c("Target 1" = "Ch1",
                    "Target 2" = "Ch2",
                    "Target 3" = "Ch3",
                    "Target 4" = "Ch4",
                    "Target 5" = "Ch5")

# define the used target genes
compute_all_positives_for <- c("PSI", "ENV", "POL", "RU5", "GAG")


#import the adapted functions
source(Rfunctions_file)

# Step 1: Set Parameters 
xlsx_file=paste(substr(qiacuity_file, 1, nchar(qiacuity_file)-4),"_dPCRfile.xlsx", sep="")
csv_file=paste(substr(qiacuity_file, 1, nchar(qiacuity_file)-4),"_dPCRfile.csv", sep="")

# Step 2: Create input files for dPCR analysis
createXLSXfile(qiacuity_file, clusterRPP30_file)
createCSVfile(clusterRainbow_file, clusterRPP30_file)

# Step 2: Run the dPCR analysis
source(Rscript_file)
