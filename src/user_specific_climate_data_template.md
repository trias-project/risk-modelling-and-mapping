# User-specific climate data template

This file explains how to use [user_specific_climate_data_template.csv](user_specific_climate_data_template.csv) with the climate-model workflow.

## Purpose

Use the CSV template when you want to bypass the downloadable CHELSA climate inputs and provide your own climate rasters to `03_fit_climate_model.R`.

This feature is optional, experimental and only works with continuous raster layers/predictors.

- If `user_specific_climate_data <- NULL`, the original CHELSA workflow is used.
- If `user_specific_climate_data` points to your CSV file, the workflow uses your climate rasters instead.

## Required columns

The CSV manifest file must contain exactly these columns:

- `period`
- `scenario`
- `var_name`
- `file_path`

## Meaning of each column

- `period`: one of `current`, `2041-2070`, or `2071-2100` (currently these are the only options available)
- `scenario`: one of `current`, `ssp126`, `ssp370`, or `ssp585`  (currently these are the only options available)
- `var_name`: the predictor name that will be used inside the model
- `file_path`: full path or relative path to a single-layer raster file (multiband/multilayer are not supported)

Different period/scenario names are, at least for now, not supported.

## Required period and scenario combinations

The current implementation requires:

- `current,current`
- `2041-2070,ssp126`
- `2041-2070,ssp370`
- `2041-2070,ssp585`
- `2071-2100,ssp126`
- `2071-2100,ssp370`
- `2071-2100,ssp585`

## Important rules

- Rows with `period = current` and `scenario = current` define the master current global climate stack.
- That current stack is used for both model training and present-day prediction.
- The `var_name` values in `current,current` must exactly match the `var_name` values in every future combination.
- Each `file_path` must point to a single-layer raster.
- Within each period/scenario combination, all rasters must share the same extent, resolution, grid alignment, and CRS.
- Across all combinations, all rasters must share the same CRS.
- Future rasters may differ from current in extent and resolution, but not in variable names.

## How to use

1. Copy and edit [user_specific_climate_data_template.csv](user_specific_climate_data_template.csv).
2. Replace the example `file_path` values with your real raster paths.
3. Keep the same variable names across current and future combinations.
4. Set the configuration in [00_configurations.R](00_configurations.R:24):

```r
user_specific_climate_data <- "./src/user_specific_climate_data_template.csv"
```

You can also use an absolute path, for example:

```r
user_specific_climate_data <- "C:/my_project/custom_climate_manifest.csv"
```

## Notes for modelers

- The example CSV is only a template. Add as many predictors as you need.
- Every predictor listed under `current,current` must also be present in each future period/scenario block (to enable prediction/future projections).
- Habitat-data customization is not part of this file yet. This template covers climate inputs only.
