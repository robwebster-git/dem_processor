import warnings; warnings.filterwarnings('ignore', 'GeoSeries.isna', UserWarning)

from dem_functions import *

from configparser import ConfigParser
import json
import click
import glob
import pandas as pd


@click.option('--config', type=click.Path(), help='Path to a configuration file for the processing')
@click.command()
def main(config):

    # Create ConfigParser() object to access the configuration file
    c = ConfigParser()
    c.read(config)

    # Populate the main variables from the configuration file
    year_1 = c.getint('main', 'year_1')
    year_2 = c.getint('main', 'year_2')
    output_dir = c.get('main', 'output_dir')
    dems_year_1 = c.get('main', 'dems_year_1')
    dems_year_2 = c.get('main', 'dems_year_2')
    refdem = c.get('main', 'refdem')
    ice_shapefile = c.get('main', 'ice_shapefile')
    basin_bins_shapefile = c.get('main', 'basin_bins_shapefile')
    preferred_shapefile_crs = c.getint('main', 'preferred_shapefile_crs')
    
    steps_to_run = json.loads(c.get('main', 'steps_to_run'))

    #  Change working directory to the output directory, and create the necessary directories if they don't exist
    build_structure(output_dir, year_1, year_2)

    print(f'\nPreparing to run step : {steps_to_run}\n')

    '''
    STEP 0

    Create a shapefile with elevation bins, from the reference DEM
    and a shapefile of RGI glacier basins.

    This is step 0 as it only needs to be done once

    '''


    if 0 in steps_to_run:
        
        min_elevation = c.getint('elevation_bins', 'min_elevation') 
        max_elevation = c.getint('elevation_bins', 'max_elevation')
        bin_size = c.getint('elevation_bins', 'bin_size')
        min_polygon_size = c.getint('elevation_bins', 'min_polygon_size')
    
        create_elevation_bin_shapefile(refdem, ice_shapefile, preferred_shapefile_crs, min_polygon_size, min_elevation, max_elevation, bin_size)


    '''     Step 1

    Create and write an extra band to each input DEM, with the number of days between
    the DEM and a reference date.

    Do this by creating a dictionary matching dates in YYYYMMDD format to total days elapsed
    since a reference date of 1st Jan 2000.  This is used to build a "days raster" 
    in order to improve dh/dt calculations.

    '''

    if 1 in steps_to_run:

        print('\nRunning step 1...')

        # Use glob to make lists of files from the two years
        dem_filepaths_1 = glob.glob(os.path.join(dems_year_1, '*.tif'))
        dem_filepaths_2 = glob.glob(os.path.join(dems_year_2, '*.tif'))

        # Populate a dictionary with the number of days between each dataset and a reference date
        days_since_refday = create_days_dictionary(dem_filepaths_1, dem_filepaths_2)

        # Loop across DEMs from year 1 and then year 2, and write a second band
        # to each with the number of days to reference date in each pixel
        for dem in dem_filepaths_1:
            write_dem_with_days_band(dem, days_since_refday, f'dems_with_days_{year_1}')
        for dem in dem_filepaths_2:
            write_dem_with_days_band(dem, days_since_refday, f'dems_with_days_{year_2}')


    '''     Step 2

    Merge DEMs with days band attached

    '''

    if 2 in steps_to_run:
        
        print('\nRunning Step 2...')

        for year in [year_1, year_2]:
            dem_filepaths = glob.glob(os.path.join(f'dems_with_days_{year}', '*.tif'))
            
            # Check if more than one dem, and merge them if so
            if len(dem_filepaths) > 1:

                print(f'\nMerging {len(dem_filepaths)} DEMs for year {year}\n')

                merge_dems_with_days_band(dem_filepaths, year, 'merged_dems')

            else:
                print(f'\n{len(dem_filepaths)} DEMs found for year {year} - no merging to do\n')


    '''     Step 3

    Pad the merged DEMs

    '''

    if 3 in steps_to_run:

        print('\nRunning Step 3...')

        # Create list of DEMs to pad
        dems_to_pad = glob.glob(os.path.join('merged_dems', 'merged_*_with_days_band.tif'))

        if len(dems_to_pad) > 0:
            # Pad the rasters so that they have exactly the same shape and extents
            pad_rasters(dems_to_pad, 'merged_dems_padded')

        else:
            print('\nNo rasters found to pad...')


    '''     Step 4

    DEM differencing

    '''
    if 4 in steps_to_run:

        print('\nRunning Step 4...')

        #  Paths to DEMs - the padded suffix is because the extents were slightly different so
        #  GDAL was used to pad the extents with nodata values so they are all the same

        dem_year_1 = os.path.join('merged_dems_padded', f'padded_merged_{year_1}_with_days_band.tif')
        dem_year_2 = os.path.join('merged_dems_padded', f'padded_merged_{year_2}_with_days_band.tif')

        dem_year_1_data = None
        dem_year_2_data = None
        dem_year_1_days = None
        dem_year_2_days = None

        profile = None

        # Set up the required output filenames
        dhdt_filename = os.path.join('dhdt', f'dhdt_{year_1}_to_{year_2}.tif') 
        years_filename = os.path.join('years', f'years_{os.path.basename(dhdt_filename)}')

        if len(dem_year_1) > 0:
            #  Read the data for both bands as masked arrays
            with rio.open(dem_year_1) as src:
                profile = src.profile.copy()
                dem_year_1_data = src.read(1, masked=True)
                dem_year_1_days = src.read(2, masked=True).astype(np.int16)

        else:
            print('\nNo padded DEMs found for year 1')

        if len(dem_year_2) > 0:
            with rio.open(dem_year_2) as src:
                profile = src.profile.copy()
                dem_year_2_data = src.read(1, masked=True)
                dem_year_2_days = src.read(2, masked=True).astype(np.int16)

        else:
            print('\nNo padded DEMs found for year 2')

        if all(v is not None for v in [dem_year_1_data, dem_year_1_days, dem_year_2_data, dem_year_2_days, profile]):
            try:
                # Create new arrays called days, in which each value is the number of days between the
                # data acquisitions in the two merged and padded input DEMs
                days = dem_year_2_days - dem_year_1_days

                # Create another array with the same info converted to years
                years = days / 365

                # DEM Differencing
                diff = dem_year_2_data - dem_year_1_data

                # Elevation difference divided by time elapsed, for each pixel with data in both datasets
                rate_of_elevation_change = diff / years

                # Update the profile since we are writing only a single band in the following rasters
                profile['count'] = 1

                # Write the elevation change dh/dt raster
                with rio.open(dhdt_filename, 'w', **profile) as dst:
                    dst.write(rate_of_elevation_change.astype(np.float32).filled(-9999), 1)

                # Write a raster where each pixel is the number of years between the two DEMs for that pixel
                with rio.open(years_filename, 'w', **profile) as dst:
                    dst.write(years.astype(np.float32).filled(profile['nodata']), 1)


            except ValueError as e:

                print(f'\nAn error occured trying to perform DEM differencing, probably as a result of missing input data : \n{e}')



    '''     Step 5

    Crop the DEM difference and years rasters with ice-areas shapefile

    '''

    output_dhdt_filename = None

    if 5 in steps_to_run:

        print('\nRunning Step 5...')

        # Set paths for DEM difference file etc
        dem_file = os.path.join('dhdt', f'dhdt_{year_1}_to_{year_2}.tif')
        output_dhdt_filename = os.path.join('dhdt', f'ice_areas_{os.path.basename(dem_file)}')

        # Crop the data and write to the output filename.  ice_shapefile is the RGI glacier basins shapefile
        write_cropped_with_shapefile(dem_file, ice_shapefile, output_dhdt_filename)

        # Do the same for the years file
        years_file = os.path.join('years', f'years_dhdt_{year_1}_to_{year_2}.tif')
        output_years_filename = os.path.join('years', f'ice_areas_{os.path.basename(years_file)}')

        # Write years file cropped to ice areas
        write_cropped_with_shapefile(years_file, ice_shapefile, output_years_filename)

    

    '''     Step 6

    Filter the dh/dt raster

    Options are set up in the configuration file

    '''

    if 6 in steps_to_run:

        print('\nRunning Step 6...')

        # Load filtering settings from configuration file

        filter_by_sigma = c.getboolean('filtering', 'filter_by_sigma')
        sigma_order = c.getint('filtering','sigma_order')
        filter_by_percentile = c.getboolean('filtering', 'filter_by_percentile')
        lower_percentile = c.getint('filtering', 'lower_percentile')
        upper_percentile = c.getint('filtering', 'upper_percentile')
        filter_by_specific_values = c.getboolean('filtering', 'filter_by_specific_values')
        lower_threshold = c.getint('filtering', 'lower_threshold')
        upper_threshold = c.getint('filtering', 'upper_threshold')


        # Input raster is the cropped DEM difference data from the previous step

        
        if output_dhdt_filename is not None:
            raster = output_dhdt_filename
        else:
            dem_file = os.path.join('dhdt', f'dhdt_{year_1}_to_{year_2}.tif')
            raster = os.path.join('dhdt', f'ice_areas_{os.path.basename(dem_file)}')

        #  Output options
        output_dir = 'dhdt_filtered'


        '''  Execute the filtering and write new rasters to the correct directory'''

        if filter_by_sigma == True:
            filter_raster(raster, output_dir, method='sigma', order=sigma_order)


        if filter_by_percentile == True:
            filter_raster(raster, output_dir, method='percentile', lower=lower_percentile, upper=upper_percentile)


        if filter_by_specific_values == True:
            filter_raster_by_specific_values(raster, output_dir, lower=lower_threshold, upper=upper_threshold)




    '''     Step 7

    Calculate zonal stats for each elevation bin and glacier basin

    '''

    if 7 in steps_to_run:

        print('\nRunning Step 7...')

        dhdts = glob.glob(os.path.join('dhdt_filtered', f'ice_areas_dhdt_*.tif'))

        # Use the shapefile newly created by STEP 0, if it exists.  If not, use the one specified in the configuration file
        bin_shapes = glob.glob(os.path.join('bins_elevation', 'basins_intersected_shapefile_*.geojson'))

        use_bins_shapefile_from_config = c.getboolean('main', 'use_basin_bins_from_config_file_in_preference_to_auto_generated')

        if len(bin_shapes) > 0 and not use_bins_shapefile_from_config:
            basin_bins_shapefile = bin_shapes[0]

        for dhdt in dhdts:
            print(f'\nProcessing file : {dhdt}')
            zonal_stats_by_basin(dhdt, basin_bins_shapefile, 'dhdt_zonal_stats', preferred_shapefile_crs)



    '''     Step 8

    Create Mean dh/dt Raster Per Elevation Bin and Basin

    '''

    if 8 in steps_to_run:

        print('\nRunning Step 8...')

        template_rasters = glob.glob(os.path.join('dhdt_filtered', 'ice_areas_dhdt_*.tif'))
        output_dir = 'dhdt_mean_value_rasters'
        stats_shapefiles = glob.glob(os.path.join('dhdt_zonal_stats', 'stats_*.geojson'))

        for stats_shapefile in stats_shapefiles:
            create_mean_value_raster(template_rasters[0], stats_shapefile, output_dir)


    '''     Step 9

    Fill a dhdt raster with a dhdt mean values raster

    '''

    if 9 in steps_to_run:
        
        print('\nRunning Step 9...')

        dhdt_files = glob.glob(os.path.join('dhdt_filtered', '*.tif'))
        means_file = glob.glob(os.path.join('dhdt_mean_value_rasters','*.tif'))[0]
        output_dir = 'dhdt_filtered_and_filled'

        for dhdt_file in dhdt_files:
            fill_dhdt_with_mean_values(dhdt_file, means_file, output_dir)


    '''     Step 10

    Create Radar Penetration Volume Bias Raster from DEM, ELA, & Upper Elevation Limit

    '''

    if 10 in steps_to_run:

        print('\nRunning Step 10...')

        # Load settings from configuration file
        ela = c.getint('penetration', 'ela')
        ela_pen = c.getint('penetration', 'ela_pen')
        upper_limit = c.getint('penetration', 'upper_limit')
        upper_limit_pen = c.getint('penetration', 'upper_limit_pen')

        output_dir = 'penetration_rasters'
        radar_penetration_raster(refdem, ice_shapefile, ela, ela_pen, upper_limit, upper_limit_pen, output_dir)


    '''     Step 11

    Calculate stable ground area weighted standard deviation

    '''
    if 11 in steps_to_run:

        print('\nRunning Step 11...')

        slope_binned_shapefile = c.get('stable_ground', 'slope_binned_shapefile')

        # Paths to DEMs - these should be dhdt rasters before being clipped to ice areas
        dem = os.path.join('dhdt', f'dhdt_{year_1}_to_{year_2}.tif')
        output_dir = 'stable_ground_stats'

        # Calculate the Area-wighted Standard Deviation of dh/dt values over stable ground
        area_weighted_std = calculate_stable_ground_area_weighted_std(dem, slope_binned_shapefile, output_dir, preferred_shapefile_crs)
        
        c.set('errors', 'area_weighted_std', str(area_weighted_std))
        
        with open(config, 'w') as f:
            c.write(f)

        print(f'Area weighted Standard Deviation is {area_weighted_std:.3f}')



    '''     Step 12

    Calculate Overall Errors and Mass Balance

    '''

    if 12 in steps_to_run:

        print('\nRunning Step 12...')

        # Load data
        stats_shp = glob.glob(os.path.join('dhdt_zonal_stats', 'stats_*.geojson'))[0]
        dhdt_file_before_filling = glob.glob(os.path.join('dhdt_filtered', 'ice_areas_dhdt*filtered.tif'))[0]
        dhdt_file = glob.glob(os.path.join('dhdt_filtered_and_filled', 'ice_areas*filtered_filled.tif'))[0]
        years_file = glob.glob(os.path.join('years', 'ice_areas_years_dhdt_*.tif'))[0]

        # Read shapefile and rename column to avoid overwriting
        gdf = gpd.read_file(stats_shp)
        gdf.rename(columns={'mean':'original_mean'}, inplace=True)

        # Remove the individual elevation bins by dissolving those together into glacier basins
        gdf.geometry = gdf.geometry.buffer(2)
        gdf = gdf.dissolve(by='RGIId').reset_index()

        # Get a rasterio profile from the dhdt file
        with rio.open(dhdt_file) as src:
            profile = src.profile.copy()
            pixel_area = src.res[0] * src.res[1]

        # Transform gdf to same CRS as the dhdt raster
        main_crs = profile['crs']
        gdf = gdf.to_crs(main_crs)

        # Calculate original valid data area by counting pixels in raster before filling
        gdf = gpd.GeoDataFrame.from_features(zonal_stats(gdf, dhdt_file_before_filling, geojson_out=True, stats='count', nodata=-9999), crs=gdf.crs)

        # Create attribute - area of each basin covered by original remotely sensed data
        gdf['rs_data_area'] = gdf['count'] * (pixel_area)

        # Remove glacier basins with no remotely sensed data
        gdf = gdf.loc[gdf['rs_data_area'] != 0]

        # Zonal stats on dhdt file
        gdf = gpd.GeoDataFrame.from_features(zonal_stats(gdf, dhdt_file, geojson_out=True, stats='count mean std', nodata=-9999), crs=gdf.crs)

        # A way of calculating basin area by multiplying count in filled rasters by pixel area
        gdf['total_basin_area'] = gdf['count'] * (pixel_area)
        del gdf['count']

        gdf.rename(columns={'mean':'dhdt_basin_mean', 'std':'dhdt_basin_std'}, inplace=True)


        ## Years Section

        # The idea here is to perform zonal statistics on each glacier basin using the years raster, 
        # in which each pixel represents the number of years between the two DEM dates.  A single number of years
        # is needed for the error calculation.  The method here uses the mean number of years for each basin, 
        # which is not perfect but will be better than using majority stats. 

        gdf = gpd.GeoDataFrame.from_features(zonal_stats(gdf, years_file, geojson_out=True, stats='mean', nodata=-9999), crs=gdf.crs)
        gdf.rename(columns={'mean': 'years'}, inplace=True)


        ## Radar Penetration Section

        penetration_raster = glob.glob(os.path.join('penetration_rasters', 'penetration-*.tif'))[0]

        with rio.open(penetration_raster) as src:
            profile_pen = src.profile.copy()

        gdf = gdf.to_crs(profile_pen['crs'])
        gdf = gpd.GeoDataFrame.from_features(zonal_stats(gdf, penetration_raster, geojson_out=True, stats='count mean', nodata=-9999), crs=gdf.crs)

        # Calculate penetration volume bias as the mean penetration per glacier basin x the number of pixels x pixel area
        gdf['Vpen_m3'] = gdf['count'] * gdf['mean'] * pixel_area


        ''' Now that all zonal statistics are done, set the GeoDataFrame to our desired output CRS '''
        preferred_shapefile_crs = c.getint('main', 'preferred_shapefile_crs')
        gdf = gdf.to_crs(epsg=preferred_shapefile_crs)          #  Note this should be based on metres and not degrees

        ''' This section covers the creation of new attributes and calculation of results and errors '''
        # Set list of ice densities to process
        densities = [850, 900]                          # kg per m3

        # Area and perimeter attributes
        gdf['basin_geom_area'] = gdf.geometry.area / (1_000_000) # convert square metres to square km
        gdf['basin_geom_perimeter'] = gdf.geometry.length / (1000) # convert metres to km

        # Create attribute of area to perimeter ratio for each basin
        gdf['basin_geom_p2a_ratio'] = gdf['basin_geom_perimeter'] / gdf['basin_geom_area']

        #  Create attribute with percentage of basin area covered by original remote sensing data
        gdf['pc_data_in_basin'] = (gdf['rs_data_area'] / gdf.geometry.area * 100).round(decimals=1)

        # Clean up old unused columns by keeping only the following columns:
        keep_cols = ['geometry', 'Name', 'RGIId',
            'dhdt_basin_mean', 'dhdt_basin_std', 'rs_data_area', 'total_basin_area', 'pc_data_in_basin',
            'basin_geom_area', 'basin_geom_perimeter', 'basin_geom_p2a_ratio', 'Vpen_m3', 'years']

        # Filter to keep only the columns specified above
        gdf = gdf[keep_cols]

        # If no glacier basin name, copy the RGI Id into the Name column
        gdf.loc[gdf['Name'].isnull(), 'Name'] = gdf['RGIId']

        # Clean out shapes with NULL values for Vpen_m3
        gdf = gdf.loc[~gdf['Vpen_m3'].isnull()]

        # Create attribute for ice volume change per basin  
        gdf['dV_m3'] = (gdf['total_basin_area'] * gdf['dhdt_basin_mean'])

        # Create Glacier Mass Balance attributes (in kg/year and in Gt/year) for each of the density values 
        for density in densities:
            gdf[f'GMB_kgyear_{density}'] = gdf['dV_m3'] * density
            gdf[f'GMB_gtyear_{density}'] = gdf[f'GMB_kgyear_{density}'] / 1_000_000_000_000
            
        # Round the numerical values in the dataframe to 4 decimal places    
        gdf = gdf.round(decimals=4)

        # Set all variable prior to calculating final mass balance errors per basin
        # NOTE - variables which are designated as strings are the keys for the correct data in the dataframe


            

        dhdt = 'dhdt_basin_mean'                    #  Mean elevation change         
        p2a = c.getfloat('errors', 'p2a')             # Perimeter (km) to area (km2) ratio 
        p2a_paulref = c.getfloat('errors', 'p2a_paulref')  # Perimeter to area ratio reported in Paul et al 2013 
        A = 'total_basin_area'                      # Basin total area in m2
        dA = 0.03 * (p2a / p2a_paulref)             # Approximate error in glacier area
        Sf = A                                      # Not really necessary but keeps the same terminology as the paper
        dp = c.getint('errors', 'uncertainty_in_density') / 1_000_000_000_000                 # Uncertainty in density
        
        Vpen = 'Vpen_m3'                            # Volume bias from radar signal penetration
        dt = 'years'                                # Time period in years
        sd_AW_dhdt = c.getfloat('errors', 'area_weighted_std')              # Area-weighted std deviation of elevation change over solid ground
        lag_distance = c.getint('errors', 'lag_distance')      # Taken from Farias-barahona, 2020

        # Feed the following variables into the error calculation function:
        # dMdt, dhdt, p2a, p2a_paulref, A, dp, p, Vpen, dt, sd_AW_dhdt, lag_distance
        total_errors = {}

        for density in densities:
            
            #define density dependent variables
            dMdt = f'GMB_gtyear_{density}'              # Gt / year mass balance
            p = density / 1_000_000_000_000             # Density

            # Calculate errors for each density value
            gdf[f'error_{density}'] = gdf.apply(lambda x: calculate_error(x[dMdt], x[dhdt], p2a, p2a_paulref, x[A], dA, dp, p, x[Vpen], x[dt], sd_AW_dhdt, lag_distance), axis=1)
            gdf[f'squared_error_{density}'] = gdf[f'error_{density}']**2
            total_errors[density] = math.sqrt(gdf[f'squared_error_{density}'].sum())


        '''  FINAL SECTION - '''
        # Setup final output shapefile path and filename
        final_stats_output = os.path.join('final_stats', f'final_stats_{year_1}_to_{year_2}.geojson')
        final_stats_output_csv = os.path.join('final_stats', f'final_stats_{year_1}_to_{year_2}.csv')

        # Write data to shapefile
        gdf.crs = CRS.from_epsg(preferred_shapefile_crs)
        gdf.to_file(final_stats_output, driver='GeoJSON')

        df = pd.DataFrame(gdf.drop(columns='geometry'))
        df.to_csv(final_stats_output_csv)

        print('\nStep 12 done...\n')


if __name__ == '__main__':
    main()
