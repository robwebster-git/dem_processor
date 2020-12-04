# DEM Processor

## Introduction

This set of utilities was written by Rob Webster to enable cryosphere geodetic mass balance and DEM differencing operations, and to speed the whole process up.

It was originally written as part of an MSc dissertation project, investigating ice elevation change in South Georgia, a sub-Antarctic island in the Southern Ocean.

It is designed to work with groups of coregistered DEMs from two time periods or "years", and to produce elevation difference maps, vector files with zonal statistics, and mass balance and error estimates.

It is tested with DEMs produced from raw TanDEM-X data, processed using interferometric techniques by Gamma Remote Sensing software, and coregistered to a reference DEM using David Shean's fantastic open source [`dem_coreg`](https://github.com/dshean/demcoreg) software (**Shean et al. 2016**).  However, it should work with any groups of suitably coregistered DEMs.

The structure is based on a 12 step program, but bears no further resemblance to Alcoholics Anonymous.  

A configuration file is created for each process, from the template provided.  This configuration file sets all the relevant variables and settings for the process, including which steps to run.  They can all be run in a single long process, or subsets of steps can be run, as long as the input dependencies are satisfied from running the previous steps at an earlier time.  This means that you don't have to run the whole process again if some problems have been found and fixed at a later stage and you wish to recompute errors, for example.

Typically, you will have several DEMs over the course of a few months during one year, and another set of DEMs for the same area from a period of several months in another year, ideally several years apart.  You wish to merge these groups of DEMs to create a single elevation model for each "year", and to investigate the elevation differences and rates of change between the two.

In order to maintain accuracy for rates of change measurements, DEM Processor will retain source date information for each pixel of the merged DEMs, rather than using a median date for each year.

An important part of the process is the use of RGI (Randolph Glacier Inventory / GLIMS) polygons which delineate glacier basins.  Zonal statistics are then calculated for individual glacier basins.  Overall statistics for larger areas can then be calculated at the end of the process, by aggregating individual basins, if required.

Details on how the error estimates are implemented can be seen in the detailed Step 12 information below, with references to the relevant papers.  

If desired, the user can implement their own error calculation functions instead.

## Steps

0. Option, if not available already - create an elevation binned shapefile from reference DEM and RGI glacier basins shapefile
1. Take coregistered outputs and add "days" band to help with accurate dt computations
2. Merge DEMs with "days" bands to create a single DEM for each input year
3. Pad these merged DEMs with nodata values so that they have exactly the same extents
4. DEM Differencing - create dh/dt raster and "years" raster
5. Crop the above rasters to ice-only areas, using RGI glacier basins shapefile
6. Filter the ice-only areas dh/dt raster (filter by sigma, percentile, or fixed limits)
7. Calculate per-glacier basin and per-elevation bin dh/dt zonal statistics ie mean dh/dt
8. Produce a raster of mean dh/dt per elevation bin and glacier
9. Fill voids in dh/dt raster using the mean values raster produced in previous step
10. Create a radar penetration volume bias raster from reference DEM, ELA, upper elevation limit, and penetration estimates
11. Calculate stable ground area-weighted standard deviation using slope bins (default 10 degrees, 0-50)
12. Calculate overall errors and mass balance

## Usage

First, you will need to install a suitable conda environment (or satisfy the python packages required in some other way).  An `environment.yml` file is supplied here which can be used to create an environment which is known to work, though it does include some extra unneccesary dependencies.

```
conda env create -f environment.yml
conda activate dem_processor
```

The above will install and activate the new environment, called `dem-processor`.

To use the code, you will need to clone this repository in the usual way.  The main program is `dem_processor.py`.

Since almost all of the configuration options are defined in the configuration file, that is the only thing you need to pass to the code.

`dem_processor.py --config /full/path/to/configuration_file.ini`

## Setting Up Directory Structure & Input Data

The code is designed to work with certain standard directories and file naming conventions.  Following this list will make sure everything is created properly.

1. Create a base directory into which all outputs will go.  This is referred to in the configuration file as the `output_dir`.
2. Create two directories within this base directory - these are for the source DEMs for each of the time-periods.  I like to call these `dems_YYYY` where `YYYY` is the year.  Into these copy the DEMs you want to work with.
3. Create a directory called `reference_dem` and place into it a copy of your reference DEM, in GeoTiff format.
4. Create a directory called `shapefiles` and place into it your RGI glacier basins shapefile.  You may also wish to put a shapefile of the area to use for calculating stable i.e. non-glaciated ground statistics to feed into the error computations, and if you already have a shapefile of elevation bins intersected with glacier basins it can go in here too.  If not, you can create one using step 0 of the processing chain.
5. Inside the output dir, the program will add to these to create the following directory structure (this example is for DEMs in 2011 and 2013)

```
├── dems_2011
├── dems_2013
├── dems_with_days_2011
├── dems_with_days_2013
├── dhdt
├── dhdt_filtered
├── dhdt_filtered_and_filled
├── dhdt_mean_value_rasters
├── dhdt_zonal_stats
├── final_stats
├── merged_dems
├── merged_dems_padded
├── penetration_rasters
├── reference_dem
├── shapefiles
├── stable_ground_stats
└── years
```
The directory setup should be complete and you can now create a configuration file in order to run the program.  I usually place this configuration file in the base directory root, but you can place it anywhere as the full path must be passed to python, as you'll see below.

## Example Configuration File

A configuration file should be created for each pipeline, with the settings pertaining to that area and set of inputs.  An example for the above directory structure may be like the one below.  Note that the years here are only used to name files and to differentiate the two time periods that the DEMs fall into.  They are not used in any calculations.  The DEMs also do not have to be all in the same year - you may have a set of DEMs from an austral summer which cross from one year into another.  This is no problem - the dh/dt calculations use the datestamps in the filenames of the DEMs.

Save the following in a file with the `.ini` extension.

```
[main]
year_1 = 2011
year_2 = 2013

steps_to_run = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
output_dir = ~/dems_pipeline_2011_to_2013
dems_year_1 = ~/dems_pipeline_2011_to_2013/dems_2011
dems_year_2 = ~/dems_pipeline_2011_to_2013/dems_2013
refdem = ~/dems_pipeline_2011_to_2013/reference_dem/srtmgl1.tif
ice_shapefile = shapefiles/rgi6_4326.geojson
basin_bins_shapefile = shapefiles/rgi-basins-elevation-100m-bins.geojson

# If the flag below is set, the basin_bins_shapefile specified above is used in preference to any generated by step 0, even if it exists
use_basin_bins_from_config_file_in_preference_to_auto_generated = True
preferred_shapefile_crs = 3762

[elevation_bins]
# Elevation range for binning - units of reference DEM
min_elevation = 0
max_elevation = 2000
bin_size = 100
min_polygon_size = 10000 # m2

[filtering]
filter_by_sigma = True
sigma_order = 3
filter_by_percentile = False
lower_percentile = 2
upper_percentile = 98
filter_by_specific_values = False
lower_threshold = -35
upper_threshold = 35

[penetration]
ela = 300
ela_pen = 0
upper_limit = 2000
upper_limit_pen = 5

[stable_ground]
slope_binned_shapefile = ~/dems_pipeline_2011_to_2013/shapefiles/stable_ground.geojson

[errors]
p2a = 2.188
p2a_paulref = 5.03
uncertainty_in_density = 60
lag_distance = 340
```

## Steps in More Detail

`Step 0`

Optional initial step - if not available already, create an elevation binned shapefile the from reference DEM and RGI glacier basins shapefile

![South Georgia - RGI glacier basins and 100m elevation bins](images/sg-basins-elevation-bins-step0.png)

`Step 1`

Take coregistered outputs and add "days" band to help with accurate dt computations

`Step 2`

Merge DEMs with "days" bands to create a single DEM for each input year

`Step 3`

Pad these merged DEMs with nodata values so that they have exactly the same extents

`Step 4`

DEM Differencing - create dh/dt raster and "years" raster

![DEM differencing output over all areas](images/sg-dhdt-allareas.png)

`Step 5`

Crop the above rasters to ice-only areas, using RGI glacier basins shapefile

![DEM differencing output cropped to ice-only areas](images/sg-dhdt-iceareas.png)

`Step 6`

Filter the ice-only areas dh/dt raster

![Step 5 output filtered to remove outliers, in this case using a sigma filter](images/sg-sigmafiltered.png)

Filtering is useful to remove outliers in the elevation change dataset.  This can remove some errors resulting from phase unwrapping errors in the InSAR processing, for example.  There are a number of methods that can be employed here, and it is recommended that a few are tried and inspected to see how they perform.
They are set up in the configuration file.

```
[filtering]
filter_by_sigma = True
sigma_order = 3
filter_by_percentile = False
lower_percentile = 2
upper_percentile = 98
filter_by_specific_values = False
lower_threshold = -35
upper_threshold = 35
```

If the boolean values are set to `True` then the associated filtering is undertaken.  These produce separate rasters for each process, rather than being chained together.  Usually you would test each method and settle on a single mmethod to take into further processing.

`sigma` method is filtering based on mean and standard devaition of the dh/dt dataset.  A common approach is to use a 3-sigma filter, removing outliers 3 standard deviations beyond the mean.

`percentile` is used where a little more control is required, perhaps because the distribution is such that you wish to trim a little more off one end of the distribution that the other.  A common setting for this is a 1-99% percentile or perhaps 2-98% filter.

`specific_values` is used when you have inspected your data and found that there are clear thresholds that you wish to filter the data by.  For example, you may decide that elevation change values greater than -50m per year are unrealistic and arise from errors.  These can be filtered out by specifying min and max elevation change values and setting `filter_by_specific_values` to `True`. 

`Step 7`

Calculate per-glacier basin and per-elevation bin dh/dt zonal statistics ie mean dh/dt

`Step 8`

Produce a raster of mean dh/dt per elevation bin and glacier

![Rasterized zonal stats - mean values per elevation bin and glacier basin](images/sg-meanbinvalues-raster.png)

`Step 9`

Fill voids in dh/dt raster using the mean values raster produced in previous step

![dh/dt raster with voids filled with mean values raster from previous step](images/sg-dhdt-filledwithmeans.png)

`Step 10`

Create a radar penetration volume bias raster from reference DEM, ELA, upper elevation limit, and penetration estimates

![Estimated radar penetration between estimated Equilibrium Line Altitude and an upper elevation bound](images/sg-radar-penetration.png)

In the configuration file, the section on penetration is as follows:

```
[penetration]
ela = 300
ela_pen = 0
upper_limit = 2000
upper_limit_pen = 5
```

These are the settings used to generate the raster above, in conjuntion with the reference DEM.  The `ela` is an estimate of the equilibrium line altitude, below which radar penetration into bare ice can be treated as being zero.  The `upper_limit` is set to 2000m in this case, as only data below this is used in the calculation of mass balance as a result of the high slope and low accuracy of measurements above this level in the case of South Georgia.  An estimated penetration (`ela_pen`) of 0m at the ELA, rising linearly to a penetration at the upper elevation limit (`upper_limit_pen`) of 5m at 2000m is used to create the raster.  This is be integrated across the relevant basins in step 12, to give an estimated volume bias due to radar penetration for that basin.

`Step 11`

Calculate stable ground area-weighted standard deviation using slope bins (default 10 degrees, 0-50)

A calculation of an overall area-weighted standard deviation value for non-glaciated terrain is required to assess the accuracy of the dh/dt raster.  If the ground is stable and should therefore be the same elevation in both datasets, the mean of the actual elevation difference values should be close to zero and the standard deviation of the values gives an approximation to the errors.

In order to implement this, you must supply a shapefile delineating the "stable ground" areas.  Furthermore, it is expected that since accuracy varies with slope, this will be slope-binned into bins of 5 or 10 degrees.  The final error estimate will be an area-weighted average of these bins.

The path to this shapefile is specified in the configuration file.  Usually you would place this in the `shapefiles` directory within the base directory, during the setup of the project.

```
[stable_ground]
slope_binned_shapefile = ~/dems_pipeline_2011_to_2013/shapefiles/stable_ground.geojson
```

Once this value is calculated, it is written to the configuration file so that later steps which require it can be run without rerunning this processing intensive step.

`Step 12`

Calculate overall errors and mass balance

The final phase of the process calculates errors and mass balance estimates for each glacier basin and writes all these attributes, plus intermediate ones, to a shapefile in the subdirectory `final_stats`.

The main method used is to follow the uncertainty analysis of **Braun et al. (2019)**, **Malz et al. (2018)**, **Farias-Barahona et al. (2020)**, **Paul et al. (2013)**, **Rolstad et al. (2009)**, with earlier influence from **Nuth et al. (2007)**, see references below.

For more information on exactly how this works, please see the above papers, and inspect the Step 12 code within `dem_processor.py` and the function `calculate_error` within `dem_functions.py`.

A section of the configuration file lets you set some of the parameters.  

`p2a` is an estimated value of perimeter to area for the glacier basins in the study area.  For more information on this and where it comes from, see the papers referenced above.  

`p2a_paulref` is a reference value of perimeter to area ratio that was arrived at experimentally by Paul et al. and is reproduced here.  It is used to scale the estimates of error in glacier boundaries which varies with ratio of perimeter to area.

`uncertainty_in_density` is again lifted from the above papers and can be adjusted.  In the current setup, +/- 60kg/m3 is used for the two standard ice density values calculated for : 850 and 900 kg/m3

`lag_distance` is an average measure of the scale of spatial autocorrelation in elevation change in a glacier.  The number 340 is taken from the paper by Farias-Barahona et al. 2020, which was measured for South Georgia Island.

`area_weighted_std` is added to the configuration file by step 11, so that it is available for step 12 even if you wish to run step 12 again on it's own.

```
[errors]
p2a = 2.188
p2a_paulref = 5.03
uncertainty_in_density = 60
lag_distance = 340
area_weighted_std = 2.8991783665400344
```

The final zonal statistics shapefile, as well as all of the vector and raster layers for the intermediate steps, are now available under the base directory.

A CSV file with the same essential contents as the `final_stats_...` shapefile is also produced, to make it easier to use with software such as Excel.

![Final Stats Example - part of the attribute table](images/sg-finalstats-attributes.png)

These can of course be used in many ways, one example is to map ice mass balance by glacier basin, to see which basins are gaining mass and which are losing, as in the example below (in this case the scale is from red ~ -0.1755 to blue ~ 0.1662 Gigatonnes per year).

![Final Stats - Can be processed further and visualised as desired, with QGIS, python etc.](images/sg-finalstats-example.png)

![Final Stats - In this example the errors are mapped, showing that in general smaller glacier basins are mores susceptible to errors](images/sg-finalstats-example-errors.png)

The details of these examples are not important, it is just an illustration of the outputs and what can be done with them.

## Improvements Required

1. Add code to automatically generate a slope binned shapefile from the reference DEM and a shapefile with basic area outline polygons. Currently this is a bit fiddly to get correct (I could not get satisfactory results from `gdaldem slope`, or from the python package `richdem`, but instead had to use QGIS to generate this for my project.)

2. Adding a step that masks out areas of known errors for each year, such as phase unwrapping errors.  The creation of the mask is a manual step performed by drawing polygons around erroneous areas on the merged rasters.  However, there is currently no automated step to apply such a mask polygon file to a merged elevation raster.  This must be done manually between steps 3 and 4 if needed, modifying the outputs from step 3 but leaving their names and location the same, before continuing with the rest of the steps.

3. Adding different error calculation methods to extend flexibility

4. Testing with a variety of different input datasets.

# Conclusion

I hope this can be useful for academics, students, or enthusiasts to save time for a geodetic mass balance or similar project.

Rob Webster, December 2020

# REFERENCES

Braun et al. 2019, and Farias-Barahona et al. 2020, who in turn were
influenced by the work of Malz et al. 2018, Rolstad et al. 2009 and Nuth et al. 2007.

Braun, Matthias H, Philipp Malz, Christian Sommer, David Farias-Barahona, Tobias Sauter, Gino Casassa, Alvaro Soruco,
Pedro Skvarca, and Thorsten C Seehaus (2019). “Constraining glacier elevation and mass changes in South America”.
In: Nature Climate Change 9.2, pp. 130–136.

Farias-Barahona, David, Christian Sommer, Tobias Sauter, Daniel Bannister, Thorsten C Seehaus, Philipp Malz, Gino
Casassa, Paul A Mayewski, Jenny V Turton, and Matthias H Braun (2020). “Detailed quantification of glacier elevation
and mass changes in South Georgia”. In: Environmental Research Letters 15.3, p. 034036.

Malz, Philipp, Wolfgang Meier, Gino Casassa, Ricardo Ja˜na, Pedro Skvarca, and Matthias H Braun (2018). “Elevation
and mass changes of the Southern Patagonia Icefield derived from TanDEM-X and SRTM data”. In: Remote Sensing
10.2, p. 188.

Nuth, C, J Kohler, HF Aas, O Brandt, and JO Hagen (2007). “Glacier geometry and elevation changes on Svalbard
(1936–90): a baseline dataset”. In: Annals of Glaciology 46, pp. 106–116.

Paul, Frank, Nicholas E Barrand, S Baumann, Etienne Berthier, Tobias Bolch, K Casey, Holger Frey, SP Joshi, V
Konovalov, Raymond Le Bris, et al. (2013). “On the accuracy of glacier outlines derived from remote-sensing data”.
In: Annals of Glaciology 54.63, pp. 171–182.

Rolstad, C, T Haug, and B Denby (2009). “Spatially integrated geodetic glacier mass balance and its uncertainty based
on geostatistical analysis: application to the western Svartisen ice cap, Norway”. In: Journal of Glaciology 55.192,
pp. 666–680.

Shean, David E, Oleg Alexandrov, Zachary M Moratto, Benjamin E Smith, Ian R Joughin, Claire Porter, and Paul
Morin (2016). “An automated, open-source pipeline for mass production of digital elevation models (DEMs) from veryhigh-
resolution commercial stereo satellite imagery”. In: ISPRS Journal of Photogrammetry and Remote Sensing 116,
pp. 101–117.