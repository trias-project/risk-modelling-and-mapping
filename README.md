# Alien species risk modelling and mapping

This repository contains the framework and R code for predicting the distribution of alien species in Europe at 1 km<sup>2</sup> (climate model) and 1 km<sup>2</sup> resolution (habitat/land cover model), initially as part of the TrIAS project. The workflows is further developed as part of the projects [GuardIAS](https://guardias.eu/) and [OneSTOP](https://onestop-project.eu/).
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

* **occurrence_thinning_method**: When a species has more than 10000 occurrences, records are thinned to 10000 to match the number of pseudoabsences. Options:<ul>
<li><code>"random"</code>: randomly samples 10000 occurrences.</li>
<li><code>"kmeans_clustering"</code>: performs k-means clustering in environmental space and selects 10000 cluster centroids, ensuring the thinned occurrences represent the broadest environmental variation.</li>
</ul>

* **mtp_probabilities**: Defines the minimum training presence (MTP) thresholds used to convert continuous favorability predictions into binary presence/absence maps. For each value in `mtp_probabilities`, the workflow removes the lowest *x*% of occurrence probabilities and uses the next-lowest value as the threshold. Example: `mtp_probabilities = c(0.01, 0.05)` will produce binarized maps where the threshold corresponds to the lowest favorability that remains after the 1% and 5% lowest-favorability occurrences are removed.

* **country_of_interest**: Either "Europe" or a European country. If a specific European country is provided, all output maps (PNG, PDF, and raster format) will be generated for that country only.

* **update_files**: Whether or not all needed raster input files should be downloaded again. Options are *"ask"*: the user will be prompted if a specific file should be downloaded again (any missing files are downloaded by default), *"yes"*: all files are downloaded again and *"no"*: only missing files are downloaded. When you select "ask", a pop-up window will appear for each layer asking whether it should be downloaded again. The workflow will pause and cannot continue automatically until you respond to each prompt. Note that these pop-ups may sometimes open behind the RStudio window, so you may need to minimize RStudio to locate them. 

## Executing the workflow

To execute this workflow, run the **06_run_wiSDM.R** script, stored in the `src` folder. This script automatically runs the following scripts in the designated order:

0. **Script 00_configurations**: Specifies the workflow's configurations (e.g., project name, species to model, thinning method, MTP threshold settings). IMPORTANT: users must actively set these fields themselves prior to running the workflow.
1. **Script 01_prepare_files_and_folders**: Sets up the folder structure and downloads the files (climate rasters, habitat predictors, spatial boundaries,...) necessary to run the workflow.
2. **Script 02_global_occurrence_download.R**: Retrieves occurrence data for the species of interest defined in script 00 from the Global Biodiversity Information Facility (GBIF). 
3. **Script 03_fit_climate_model.R**: Builds a global-scale climate-only species distribution model (SDM) for each species of interest, at a resolution of 1 km<sup>2</sup>. The results of this model are presented at the level of the country of interest that was specified in the configurations script and they can be found in the folder `./data/projects/<your project>/species/Climate`.
4. **Script 04_fit_habitat_model.R**: Generates a European-scale habitat-only species distribution model for the specified species at a resolution of 1 km<sup>2</sup> and integrates these predictions with the 1 km² climate-only predictions from script 03. The two prediction layers are integrated by using the geometric mean to generate a final suitability map at 1 km<sup>2</sup> resolution that reflects both climate suitability and habitat (land cover) suitability. Final predictions are generated for both current conditions and for two future periods, 2041-2070 and 2071-2100, under different climate change scenarios (SSP1-2.6, SSP3-7.0, and SSP5-8.5). The results of the habitat model can be found in the folder `./data/projects/<your project>/species/Habitat`, while the final predictions, combining both the habitat and the climate suitability, can be found in the `./data/projects/<your project>/species/Combined` folder.
<br>
 
## What does this modelling workflow do?
1.	**Generates habitat suitability maps using machine learning.**
The workflow requires only a species name and then fits an ensemble of ten machine-learning algorithms to estimate both climate suitability (using global occurrences) and habitat suitability (using European occurrences). For each type (habitat and climate), the favourability maps from all ten algorithms are first stacked and analysed using Principal Component Analysis (PCA). For each algorithm, the spatial variance of its predictions projected onto the first principal component (PC1) is calculated. The five algorithms with the highest variance along PC1, representing those that contribute most strongly to the dominant pattern of variation across the study area, are then used to generate the final suitability map. These two maps (one for the climate model and one for the habitat model) are then combined using the geometric mean in order to produce a final overall suitability prediction. All maps are generated automatically for current conditions. Climate suitability maps, along with the final combined maps, are also produced for three standard Shared Socioeconomic Pathways (SSP1-2.6, SSP3-7.0, and SSP5-8.5) for the periods 2041–2070 and 2071–2100.
2.	**Implements best practices for pseudoabsence placement:** 
    * **Climate model**: Pseudoabsences are sampled within the same biomes as species presences but excluded from presence grid cells. A taxonomic sampling effort grid (bias grid) captures the sampling intensity of the higher taxon and is used to weight grid cells, assigning greater weight to well-sampled areas.
    * **Habitat model**: Pseudoabsences are sampled across all of Europe, excluding presence grid cells. In ecoregions with presences, grid cells are weighted using the bias grid, while outside these ecoregions all cells are assigned the minimum weight (1).
3.	**Detects and removes highly correlated predictors.**<br>
Highly correlated predictors can have undesirable effects and confuse the interpretation of variable importance.
4.	**Automatic generation of confidence maps** for each suitability map. These maps illustrate prediction uncertainty across the study area by calculating the population standard deviation of the predictions produced by the five algorithms for both the climate and habitat models. The standard deviation of the combined prediction is then computed using the appropriate error-propagation formula for the geometric mean. Note that the standard deviation is always calculated from the mean of the five algorithms, whereas the final suitability predictions are based on their median. Consequently, the confidence maps provide a relative measure of model uncertainty but should not be interpreted as the uncertainty of the median prediction itself.
<br>

## Future functionalities
1.	**Evaluation of spatial autocorrelation** in model residuals to detect clustering effects. If autocorrelation is high, the workflow will recommend (but not automatically apply) occurrence thinning.
2.  **Cross-validated Boyce Index calculation** to provide a robust, presence-only model performance measure. This will be implemented in script 05_cross_validation.R, which is not functional yet.
<br>

## Contributors

[List of contributors](https://github.com/trias-project/risk-modelling-and-mapping/contributors)
<br>

## Funders

This workflow is being developed in the framework of the projects [GuardIAS](https://guardias.eu/) and [OneSTOP](https://onestop-project.eu/), both receiving funding from the European Union’s Horizon Europe Research and Innovation Programme:
- ID No 101181413 (GuardIAS)
- ID No 101180559 (OneSTOP)

## References
Davis AJS, Groom Q, Adriaens T, Vanderhoeven S, De Troch R, Oldoni D, Desmet P, Reyserhove L, Lens L and Strubbe D (2024) Reproducible WiSDM: a workflow for reproducible invasive alien species risk maps under climate change scenarios using standardized open data. Front. Ecol. Evol. 12:1148895. doi: [10.3389/fevo.2024.1148895](https://doi.org/10.3389/fevo.2024.1148895)

**Important note**: The wiSDM version on the current main branch (v2.0.0) represents a complete restructuring and substantial update of the workflow and does not exactly reproduce the implementation described in Davis et al. (2024). The published workflow is preserved in release v1.1.0.

## License

[MIT License](https://github.com/trias-project/risk-modelling-and-mapping/blob/master/LICENSE)
