# Southeast Australia Megafauna extinction  

Code required to run model for inferring timing of megafauna extinction and first human arrival in south-eastern Australia

Reference: Saltré, F., Chadoeuf, J., Peters, K.J., McDowell, M.C., Friedrich, T., Timmermann, A., Ulm, S. & Bradshaw, C.J.A. (2019) Climate-human interaction associated with southeast Australian megafauna extinction patterns. Nature Communications, 10, 53

## Abstract

The mechanisms leading to megafauna (>44 kg) extinctions in Late Pleistocene (126,000—12,000 years ago) Australia are highly contested because standard chronological analyses rely on scarce data of varying quality and ignore spatial complexity. Relevant archaeological and palaeontological records are most often also biased by differential preservation resulting in under-representated older events. Chronological analyses have attributed megafaunal extinctions to climate change, humans, or a combination of the two, but rarely consider spatial variation in extinction patterns, initial human appearance trajectories, and palaeoclimate change together. Here we develop a statistical approach to infer spatio-temporal trajectories of megafauna extirpations (local extinctions) and initial human appearance in south-eastern Australia. We identify a combined climate-human effect on regional extirpation patterns suggesting that small, mobile Aboriginal populations potentially needed access to drinkable water to survive arid ecosystems, but were simultaneously constrained by climate-dependent net landscape primary productivity. Thus, the co-drivers of megafauna extirpations were themselves constrained by the spatial distribution of climate-dependent water sources

## R-scripts

Main R script for inferring megafauna extinction timings (Megafauna_Australia_v2.r) and first human arival timings (HumanInference.r)
Main script return timing of extinction/arrival on geographic coordinate provided in GridCoordinate_v2.csv. 

Also included the validation procedure for these models as a function of 2 scenarios of colonisation/extinction:

(1) the first scenario describes a gradual colonisation (or gradual extirpation) across space: "MegafaunaValidation(gradient).r" and "MegafaunaValidation(2entrances).r"

(2) the second scenario describes two entrances of colonisation/extirpation in the landscape: "MegafaunaValidation(2entrances).r" and "HumanValidation(2entrances).r"

These scripts require the following associated function: "Genpop(gradient).r" for scenario 1 and "Genpop(2entrances).r" for scenario 2

## Matlab scripts

Included an example of Matlab script to map the timing of extinction and arrival and calculate the bearing of each: "Megafauna_analyses.m"

Two examples of scripts are included to construct and compared generalized least-squares model to determine which predictor among the climate variables and initial human colonization best described spatial variation in the bearing of megafauna extirpations.

(1) script to analyse in CLimate only areas: ScriptGLS_Bearing(Climateonly).r

(2) script to analyse in Climate-Human areas: ScriptGLS_Bearing(HumanClimate).r

## Excel files

the excel file "SourceData_Figs3 -SI4-SI5(Saltre_etal).xlsx" contains the full distribution of values underlying all reported averages in graphs and charts of Fig. 3 and Supplementary Figs 4 & 5




***************
Frederik Saltre, Flinders University, frederik.saltre@flinders.edu.au January 2019

![Screen Shot 2022-04-19 at 4 08 06 pm](https://user-images.githubusercontent.com/46954120/163941558-54b18035-5a74-44a2-984f-da05da11d048.png)
