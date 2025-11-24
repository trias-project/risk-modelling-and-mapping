# Alien species risk modelling and mapping

This repository contains the framework and R code for predicting the distribution of alien species in Europe at 5 km<sup>2</sup> (climate model) and 1 km<sup>2</sup> resolution (habitat model) as part of the TrIAS project. 
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
2.   Have an active **GBIF account**. Before running the workflow for the first time, store your GBIF username, password, and email address in your `~/.Renviron` file as instructed in the script 00_configurations.R.
3. **Clone this repository** to your local computer.
4. Set your specific configurations for the workflow in the **00_configurations** script (see below). 
<br>

## Prior to running the workflow

Before you execute this workflow, specify the following configurations in the **00_configurations script**, which is stored in the `src` folder:

* **project**: The name of your project. A folder with this name will be created automatically in the `./data/projects` folder, and all workflow outputs will be stored there.

* **species_to_model**: Specify the species you want to model as a character string (e.g., `"Vespa velutina"`). Multiple species can be provided as a character vector (e.g., `c("Vespa velutina", "Aedes albopictus")`). Use Latin binomials only (no authorship or year).

* **occurrence_thinning_method**: When a species has more than 10,000 occurrences, records are thinned to 10,000 to match the number of pseudoabsences. Options:
  * `"random"`: randomly samples 10,000 occurrences.
  * `"kmeans_clustering"`: performs k-means clustering in environmental space and selects 10,000 cluster centroids, ensuring the thinned occurrences represent the broadest environmental variation.<br>
<br>
* **mtp_probabilities**: Defines the minimum training presence (MTP) thresholds used to convert continuous favorability predictions into binary presence/absence maps. For each value in `mtp_probabilities`, the workflow removes the lowest *x*% of occurrence probabilities and uses the next-lowest value as the threshold. Example: `mtp_probabilities = c(0.01, 0.05)` will produce binarized maps where the threshold corresponds to the lowest favorability that remains after the 1% and 5% lowest-favorability occurrences are removed.

* **country_of_interest**: Currently inactive, do not modify this parameter. The workflow currently produces predictions only for the whole of Europe. Future versions will allow masking outputs to a user-defined country.

## Executing the workflow

To execute this workflow, run the **06_run_wiSDM.R** script, stored in the `src` folder. This script automatically runs the following scripts  in the designated order:

0. **Script 00_configurations**: Specifies the workflow's configurations (e.g., project name, species to model, thinning method, MTP threshold settings). IMPORTANT: users must actively set these fields themselves prior to running the workflow.
1. **Script 01_prepare_files_and_folders**: Sets up the folder structure and downloads the files (climate rasters, habitat predictors, spatial boundaries,...) necessary to run the workflow.
2. **Script 02_global_occurrence_download.R**: Retrieves occurrence data for the species of interest defined in script 00 from the Global Biodiversity Information Facility (GBIF). 
3. **Script 03_fit_climate_model.R**: Builds a global-scale climate-only species distribution model (SDM) for each species of interest, at a resolution of 5 km<sup>2</sup>. The results of this model can be found in the folder `./data/projects/<your project>/species/Climate`.
4. **Script 04_fit_habitat_model.R**: Generates a European-scale habitat-only species distribution model for the specified species at a resolution of 1 km<sup>2</sup> and integrates these predictions with the 5 km² climate-only predictions from script 03. The two prediction layers are integrated by using the geometric mean to generate a final suitability map at 1 km<sup>2</sup> resolution that reflects both climate suitability and habitat suitability. Final predictions are generated for both current conditions and for two future periods, 2041-2070 and 2071-2100, under different climate change scenarios (SSP1-2.6, SSP3-7.0, and SSP5-8.5). The results of the habitat model can be found in the folder `./data/projects/<your project>/species/Habitat`, while the final predictions, combining both the habitat and the climate suitability, can be found in the `./data/projects/<your project>/species/Combined` folder.
<br>
 
## What does the TrIAS modeling workflow do?
1.	**Generates habitat suitability maps using machine learning.**
The workflow requires only a species name and then fits an ensemble of ten machine-learning algorithms to estimate both climate suitability (using global occurrences) and habitat suitability (using European occurrences). For each type (habitat and climate), the predictions from these ten algorithms are summarized using a PCA, and the median of the five models that explain the most variation along the first PCA axis is used to generate the final suitability map. These two maps are then combined using the geometric mean in order to produce a final overall suitability prediction. All maps are generated automatically for current conditions, and climate suitability maps, along with the final combined maps, are also produced for three standard Shared Socioeconomic Pathways (SSP1-2.6, SSP3-7.0, and SSP5-8.5) for the periods 2041–2070 and 2071–2100.
2.	**Addresses geographic sampling bias.**
3.	**Implements best practices for pseudoabsence placement:** <br>
* **Climate model**: Pseudoabsences are sampled within the same biomes as species presences but excluded from presence grid cells. A taxonomic occurrence grid (bias grid) captures the sampling intensity of the higher taxon and is used to weight grid cells, assigning greater weight to well-sampled areas.
* **Habitat model**: Pseudoabsences are sampled across all of Europe, excluding presence grid cells. In ecoregions with presences, grid cells are weighted using the bias grid, while outside these ecoregions all cells are assigned the minimum weight (1).
4.	**Detects and removes highly correlated predictors.**<br>
Highly correlated predictors can have undesirable effects and confuse the interpretation of variable importance.
<br>
<br>

## Future functionalities
1.	**Automatic generation of confidence maps** for each suitability map. These will illustrate prediction uncertainty across the study extent.
2.	**Evaluation of spatial autocorrelation** in model residuals to detect clustering effects. If autocorrelation is high, the workflow will recommend (but not automatically apply) occurrence thinning.
<br>

## Contributors

[List of contributors](https://github.com/trias-project/risk-modelling-and-mapping/contributors)
<br>

## License

[MIT License](https://github.com/trias-project/risk-modelling-and-mapping/blob/master/LICENSE)
