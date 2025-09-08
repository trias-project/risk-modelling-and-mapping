# Alien species risk modelling and mapping

This repository contains the framework and R code for predicting the distribution of alien species throughout Belgium and greater Europe at 5 km<sup>2</sup> (climate model) and 1 km<sup>2</sup> resolution (habitat model) as part of the TrIAS project. 
<br>

## Repo structure

```
├── README.md              : Description of this repository
├── LICENSE                : Repository license
├── risk-modelling-and-mapping.Rproj : RStudio project file
├── .gitignore             : Files and directories to be ignored by git
│
├── data
│   ├── external          : external files required to run the model. The majority of these files will be downloaded and stored in the right folders by running script 01_prepare_files_and_folders.R.
│   
│
└── src                    : R Code
```


## Requirements to run this workflow
1.  **RStudio** installed on your local computer.
2.   Have an active **GBIF account**. In script 02, you will need to enter your GBIF username, password, and email address to enable the species occurrence data download.
3. **Clone this repository** to your local computer.
4. Provide a **project name**, the **name(s) of the study species**, and **a country of interest** in the 00_configurations script. 
<br>

## Executing the workflow
 
To execute this workflow, run the following scripts (stored in the `src` folder) in the designated order:

0. **Script 00_configurations**: Specifies the workflow configurations, including project name, species to be modeled, and country of interest. IMPORTANT: users must actively set these fields themselves.
1. **Script 01_prepare_files_and_folders**: Sets up the folder structure and downloads the files (climate rasters, habitat predictors, spatial boundaries,...) necessary to run the workflow.
2. **Script 02_global_occurrence_download.R**: Retrieves occurrence data for the species of interest defined in script 00 from the Global Biodiversity Information Facility (GBIF). To allow the data download, a pop-up will appear, requesting you to enter your GBIF username, password, and email address.
3. **Script 03_fit_global_model.R**: Builds a global-scale climate-only species distribution model (SDM) for each species of interest, at a resolution of 5 km<sup>2</sup>. 
4. **Script 04_fit_European_model.R**: Generates European-level SDMs for the specified species at a resolution of 1 km<sup>2</sup>
5. **Script 05_Make_country_level_predictions.R**: Predicts species distributions and generates confidence maps under different climate change scenarios (RCP 2.6, RCP 7.0, and RCP 8.5) for a country or region of interest.
<br>
 
## What does the Trias modeling workflow do?
1.	Automatically generates habitat suitability maps using machine learning. 
Our workflow requires only a species name and generates an ensemble of machine learning algorithms stacked together as a meta-model to produce the final habitat suitability map at 1 km<sup>2</sup> resolution. Maps are generated automatically for standard IPCC greenhouse gas emission scenarios (RCP's 2.6, 7.0, and 8.5).  
2.	Automatically generates confidence maps for each habitat suitability map. These illustrate confidence of each individual prediction across your study extent.
3.	Addresses geographic sampling bias
4.	Incorporates best practices for the placement of pseudoabsences: <br>
* Global model: Pseudoabsences are sampled within the same biomes as species presences but excluded from presence grid cells. A taxonomic occurrence grid (bias grid) captures the sampling intensity of the higher taxon and is used to weight grid cells, assigning greater weight to well-sampled areas.
* European model: Pseudoabsences are sampled across all of Europe, excluding presence grid cells. In ecoregions with presences, grid cells are weighted using the bias grid, while outside these ecoregions all cells are assigned the minimum weight (1).
5.	Detects and removes highly correlated predictors. Highly correlated predictors can have undesirable effects and confuse the interpretation of variable importance
6.	Integrates multiple machine learning algorithms to predict habitat suitabilities. It has been consistently demonstrated that the choice of algorithm has the largest impact on predicted suitability.
7.	Assesses spatial autocorrelation in the residuals to assess the impacts of clustering. If high, we recommend the employment of thinning.
<br>

## Contributors

[List of contributors](https://github.com/trias-project/risk-modelling-and-mapping/contributors)
<br>
## License

[MIT License](https://github.com/trias-project/risk-modelling-and-mapping/blob/master/LICENSE)
