# dPCR-HIV
Using Rainbow on Qiacuity

V1.0.0 (2026-08-14)

#Setting up run: 
1.	For each sample you need both HIV and RPP30 probes/primers
2.	Every biological sample should be given a unique name, reuse same for replicates with HIV primer/probes and RPP30 primer/probes (if reusing shearing controls, you need to make changes to the script).
3.	Water-only control should be labelled WATER. 
4.	Probes should be written: PSI = “FAM", ENV = "HEX",  POL = "TAMRA",  RU5 = "ROX", GAG = "Cy5", RPP30 = "FAM", RPP30 3p= "HEX"

#1. Export data, do initial quality check 

After run: Go to analysis in the Qiacuity software:
1.	Initial quality check: Check that the thresholds are set correctly. Select all wells, go to 1D Scatterplot. Adjust thresholds if needed (you cannot adjust this downstream). Generate report by checking box “add all to report”. Click Generate report and check all data is present as expected. If unexpected patterns are seen, you can add more PCR cycles and take another picture with a different exposure time. 
2.	Export data: Select all wells, then choose all targets (7) in dropdown menu.	
a.	Run analysis
b.	On top of table, click ”Export to CSV…”, choose Current results (becomes the “analysis” file)
3.	Select Rainbow wells (psi, env, …), choose all targets (5)
a.	Run analysis
b.	On top of table, click ”Export to CSV…”, choose Multiple occupancy (becomes the Rainbow “MO” file)
4.	Select RPP30 wells (rpp30, rpp30 3p), choose all targets (2)
a.	Run analysis
b.	On top of table, click ”Export to CSV…”, choose Multiple occupancy (becomes the RPP30 “MO” file)
5.	Save raw data for documentation: Export experiment (this is the raw data to be saved, will allow you to recreate all analysis in the Qiaquity software)
6.	Put the five files in a folder that starts with the date and your initials (e.g. 20260813PS)
7.	Put on a USB and transfer to your own computer (export file contains your raw data and goes to ELN), you can delete files on the computer when analysis in finished and experiment is properly documented.

#2. Analyze the intact/defective proviruses using R
Three files are needed: 
1.	Q4ddPCR_convertQiacuityFiles.R
2.	Q4ddPCR_Quicuity_RainbowAnalysis.R
3.	[date-initial]_parameters.R (template)

Open Rstudio (MultiplexPCRAnalyzer installed, devtools::install_github("buchauer-lab/Q4ddPCR_Analysis")
Change working directory to where your files are.

Open file “[date-initials]_parameters.R” in RStudio. Save with correct date and initials. Adjust the file names and anything else that is needed (e.g. well volume) and re-save. 
Run the parameters file: source("[PATH]/[date-initial]_parameters.R")
This should generate 3 files in the working directory. The one that ends with “_Analysis.xlsx” contains the output, the intact/defective proviruses can be found in the “Intact” tab. You also have the individual concentrations, shearing and original count data in this file. 

#3. Perform statistical analysis and plot data in Prism
Import the data into Prism GraphPad and perform statistical analysis and draw graphs. You can do some calculations in Excel, but the final figures and statistics are best made in GraphPad
Technical replicates (duplicate wells in the same plate) and Experimental replicates (different plates done on different days) should not be compared together. Technical replicates can be used to assess technical variation, experimental replicates to draw scientific conclusions. At least 3 experimental replicates are needed for statistical tests.
Compare experimental conditions use a standard test such as an unpaired two-tailed t-test. This can only be done once, otherwise you need to account for multiple testing.
When making figures, use standard settings for publications (e.g. Arial 6-8pt, minimal size but legible, line width 1)




