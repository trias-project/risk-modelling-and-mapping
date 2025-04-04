# Alien species risk modelling and mapping

This repository contains the framework and R code for predicting the distribution of alien species throughout Belgium and greater Europe at 1 km<sup>2</sup> resolution as part of the TrIAS project. 
<br>

## Repo structure

```
├── README.md              : Description of this repository
├── LICENSE                : Repository license
├── risk-modelling-and-mapping.Rproj : RStudio project file
├── .gitignore             : Files and directories to be ignored by git
│
├── data
│   ├── external          : external files required to run the model. The majority of these files will be downloaded and stored in the right folders by running script 00_prepare_files_and_folders.R.
│   
│
└── src                    : R Code
```


## Requirements to run this workflow
1.  **RStudio** installed on your local computer.
2.   Have an active **GBIF account**. In script 01, you will need to enter your GBIF username, password, and email address to enable the species occurrence data download.
3. **Clone this repository** to your local computer.
4. Provide a **project name** and the **name(s) of the study species** in the respective scripts. 
<br>

## Executing the workflow
 
To execute this workflow, run the following scripts (stored in the `src` folder) in the designated order:

0. **Script 00_prepare_files_and_folders**: Sets up the folder structure and downloads the files (climate rasters, habitat predictors, spatial boundaries,...) necessary to run the workflow.
1. **Script 01_global_occurrence_download.R**: After specifying the name of your project and your species of interest, this script creates a project folder on your local computer and retrieves occurrence data for the respective species from the Global Biodiversity Information Facility (GBIF). To allow this data download, a pop-up will appear, requesting you to enter your GBIF username, password, and email address.
2. **Script 02_fit_global_model.R**: Builds a global-scale climate-only species distribution model (SDM) for each species specified in script 01. Ensure you use the same project name as in the former script.
3. **Script 03_fit_European_model.R**: Generates European-level SDMs for the specified species. Again, use the same project name as in script 01.
4. **Script 04_Make_country_level_predictions.R**: Predicts species distributions and generates confidence maps under different climate change scenarios (RCP 2.6, RCP 4.5, and RCP 8.5) for a country or region of interest. At the moment, the workflow is only operational for Belgium, but this will be adjusted soon to incorporate more countries.
<br>
 
## What does the Trias modeling workflow do?
1.	Automatically generates habitat suitability maps using machine learning. 
Our workflow requires only a species name and generates an ensemble of machine learning algorithms stacked together as a meta-model to produce the final habitat suitability map at 1 km<sup>2</sup> resolution. Maps are generated automatically for standard IPCC greenhouse gas emission scenarios (RCP's 2.6, 4.5, and 8.5).  
2.	Automatically generates confidence maps for each habitat suitability map. These illustrate confidence of each individual prediction across your study extent.
3.	Addresses geographic sampling bias
4.	Incorporates best practices for the placement of pseudo-absences: pseudo absences are placed in the same ecoregions where presences occur. We use the global model to restrict pseudo-absences to areas of low predicted suitability. We use the taxonomic occurrence grid (aka bias grid) to not place pseudoabsences in areas of low sampling effort. The taxonomic occurrence grid summarizes the sampling effort of the higher taxon the modelled species belongs to.
5.	Detects and removes highly correlated predictors. Highly correlated predictors can have undesirable effects and confuse the interpretation of variable importance
6.	Integrates multiple machine learning algorithms to predict habitat suitabilities. It has been consistently demonstrated that the choice of algorithm has the largest impact on predicted suitability.
7.	Assesses spatial autocorrelation in the residuals to assess the impacts of clustering. If high, we recommend the employment of thinning.
<br>

## Contributors

[List of contributors](https://github.com/trias-project/risk-modelling-and-mapping/contributors)
<br>
## License

[MIT License](https://github.com/trias-project/risk-modelling-and-mapping/blob/master/LICENSE)
