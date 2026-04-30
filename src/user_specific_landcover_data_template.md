# User-specific land-cover data template

This file explains how to use [user_specific_landcover_data_template.csv](user_specific_landcover_data_template.csv) with the habitat-model workflow.

## Purpose

Use the CSV template when you want to bypass the default downloaded land-cover predictors and provide your own land-cover rasters to `04_fit_habitat_model.R`.

Currently only works with continuous land cover predictors (e.g., % cover, or distance to a given land cover type or feature, such as rivers, roads, settlements)

This feature is optional and experimental.

- If `user_specific_landcover_data <- NULL`, the original wiSDM habitat workflow is used.
- If `user_specific_landcover_data` points to your CSV file, the workflow uses your land-cover rasters instead.

## Required columns

The CSV must contain exactly these columns:

- `period`
- `scenario`
- `var_name`
- `file_path`

## Meaning of each column

- `period`: one of `current`, `2041-2070`, or `2071-2100` (currently these are the only options available)
- `scenario`: one of `current`, `ssp126`, `ssp370`, or `ssp585` (currently available scenarios)
- `var_name`: the predictor name that will be used inside the habitat model
- `file_path`: full path or relative path to a single-layer raster file

Different period/scenario names are, at least for now, not supported.

## Allowed manifest formats

Two valid patterns are presently supported:

1. Current-only habitat model:
   - only `current,current` rows are provided
   - habitat projections stay static for future climate combinations

2. Dynamic habitat model:
   - `current,current` rows are provided
   - all 6 future combinations are also provided:
     - `2041-2070,ssp126`
     - `2041-2070,ssp370`
     - `2041-2070,ssp585`
     - `2071-2100,ssp126`
     - `2071-2100,ssp370`
     - `2071-2100,ssp585`

Partial future coverage is not allowed. If defined, all period/scenario combos must be supplied.

## Important rules

- Rows with `period = current` and `scenario = current` define the master current habitat stack.
- The `var_name` values in `current,current` must exactly match the `var_name` values in every future combination, if future combinations are provided.
- Each `file_path` must point to a single-layer raster. Multiband or multilayer files are not supported.
- Within each period/scenario combination, all rasters must share the same extent, resolution, grid alignment, and coordinate reference system (CRS).
- Across all land-cover combinations, all rasters must share the same CRS.
- User-specific land-cover rasters are assumed to already be transformed consistently.
- The workflow does not center or scale user-specific land-cover rasters internally.
- If `user_specific_climate_data` is also set, the user-specific climate and land-cover manifests must use the same CRS.

## How to use

1. Copy and edit [user_specific_landcover_data_template.csv](user_specific_landcover_data_template.csv).
2. Replace the example `file_path` values with your real raster paths.
3. Keep the same variable names across current and future combinations whenever future habitat inputs are included.
4. Set the configuration in [00_configurations.R](00_configurations.R:28):

```r
user_specific_landcover_data <- "./src/user_specific_landcover_data_template.csv"
```

You can also use an absolute path, for example:

```r
user_specific_landcover_data <- "C:/my_project/custom_landcover_manifest.csv"
```

## Notes for modelers

- The example CSV is only a template. Add as many predictors as you need.
- If you want the habitat model to stay static in future combined projections, only provide `current,current` rows. Ignore/remove the ones for future scenarios.
- If you want dynamic future habitat projections, provide all 6 future scenario combinations with no gaps.
