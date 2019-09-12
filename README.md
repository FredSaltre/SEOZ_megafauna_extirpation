SOUTH-EASTERN MEGAFAUNA EXTINCTIONS  
***********************************

Code required to run model for inferring timing of megafauna extinction and first human arrival in south-eastern Australia

Reference: Saltre et al. (in review) Climate-human interaction caused southeast Australian megafauna extinction

Main R script for inferring megafauna extinction timings (Megafauna_Australia_v2.r) and first human arival timings (HumanInference.r)
Main script return timing of extinction/arrival on geographic coordinate provided in GridCoordinate_v2.csv. 

Also included the validation procedure for these models as a function of 2 scenarios of colonisation/extinction:
(1) the first scenario describes a gradual colonisation (or gradual extirpation) across space: "MegafaunaValidation(gradient).r" and "MegafaunaValidation(2entrances).r"
(2) the second scenario describes two entrances of colonisation/extirpation in the landscape: "MegafaunaValidation(2entrances).r" and "HumanValidation(2entrances).r"
These scripts require the following associated function: "Genpop(gradient).r" for scenario 1 and "Genpop(2entrances).r" for scenario 2

Included an example of Matlab script to map the timing of extinction and arrival and calculate the bearing of each: "Megafauna_analyses.m"

Two examples of scripts are included to construct and compared generalized least-squares models 
to determine which predictor among the climate variables and initial human colonization best described spatial variation in the bearing of megafauna extirpations.
(1) script to analyse in CLimate only areas: ScriptGLS_Bearing(Climateonly).r
(2) script to analyse in Climate-Human areas: ScriptGLS_Bearing(HumanClimate).r

the excel file "SourceData_Figs3 -SI4-SI5(Saltre_etal).xlsx" contains the full distribution of values underlying all reported averages in graphs and charts of Fig. 3 and Supplementary Figs 4 & 5

F Saltre Flinders University frederik.saltre@flinders.edu.au January 2019
