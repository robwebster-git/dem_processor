import warnings; warnings.simplefilter(action='ignore', category=FutureWarning)
import geopandas as gpd
import numpy as np
import rasterio as rio
from rasterio.merge import merge
from rasterio.mask import mask as robmask
from rasterio.io import MemoryFile
from rasterio.features import shapes
from rasterstats import zonal_stats
from pyproj import CRS
import re
import subprocess
import os
import datetime
import math



def create_elevation_bin_shapefile(refdem, ice_shapefile, output_crs, min_polygon_size, min_elevation, max_elevation, bin_size):

    if not os.path.exists('bins_elevation'):
        os.mkdir('bins_elevation')

    # Read reference DEM data and copy the profile
    with rio.open(refdem) as src:
        ref = src.read(1, masked=True)
        profile = src.profile.copy()

    # Read input shapefile containing RGI glacier basins
    ice = gpd.read_file(ice_shapefile)

    '''Reclassify reference DEM elevation data into binned values'''

    # Output binned raster
    binned_raster_out = f'bins_elevation/binned_elevation_min_{min_elevation}_max_{max_elevation}_binsize_{bin_size}.tif'

    # Output binned shapefile
    binned_shapefile_out = os.path.join('bins_elevation', f'raw_shapefile_{os.path.basename(binned_raster_out)[:-4]}.geojson')

    # Output binned shapefile intersected with RGI glacier basins shapefile
    binned_intersected_out = os.path.join('bins_elevation', f'basins_intersected_shapefile_{os.path.basename(binned_raster_out)[:-4]}.geojson')

    # Create bins list
    bins = [(((x+1) * bin_size) + min_elevation) for x in range((max_elevation - min_elevation) // bin_size)]
    print(f'\nExtracting binned elevation using bins with breaks at : \n{bins}')

    binned = ((np.digitize(ref, bins) + 1) * bin_size) + min_elevation
    binned = binned.astype(np.int16)
    binned[ref <= min_elevation] = profile['nodata']
    binned[ref >= max_elevation] = profile['nodata']

    profile['dtype'] = np.int16

    if not os.path.exists(binned_raster_out):
        with rio.open(binned_raster_out, 'w', **profile) as dst:
            dst.write(binned.astype(np.int16), 1)
    else:
        print('\nRaster already exists, continuing to shape extraction...')


    features = list(shapes(binned, mask=(binned!=profile['nodata']), connectivity=4, transform=profile['transform']))

    feats = [{'geometry' : {'type' : feature[0]['type'], 'coordinates' : feature[0]['coordinates']}, 'properties':{'elevation_bin' : feature[1]}} for feature in features]

    gdf = gpd.GeoDataFrame.from_features(feats, crs=profile['crs'])
    gdf = gdf.to_crs(epsg=output_crs)
    ice = ice.to_crs(epsg=output_crs)

    # Clean empty and null geometries
    gdf = gdf.loc[~(gdf.geometry.is_empty | gdf.geometry.isna())]
    gdf = gdf.loc[gdf.geometry.area >= min_polygon_size]

    # Explode MultiPolygons to turn them into Polygons
    gdf_multipolygons = gdf.loc[gdf.geometry.geom_type == 'MultiPolygon'].explode().reset_index(drop=True)
    gdf_polygons = gdf.loc[gdf.geometry.geom_type == 'Polygon']

    # Glue polygons back together in a single dataframe
    gdf = gpd.pd.concat((gdf_polygons, gdf_multipolygons), axis=0)

    # Write out raw vectorised version of the elevation bins raster - not strictly needed but handy for checking
    gdf.to_file(binned_shapefile_out, driver='GeoJSON')

    # Perform an intersection between the glacier basins and the binned elevations to get both in one shapefile
    intersect = gpd.overlay(ice, gdf, how='intersection')
    intersect.to_file(binned_intersected_out, driver='GeoJSON')

    print(f'done creating elevation bin shapefile : {binned_intersected_out}')


def build_structure(base_directory, year_1, year_2):
    #  Change working directory to the output directory as set in the CONFIGURATION
    #  Paths will be relative to this  

    os.chdir(base_directory)

    # Create all the necessary directories for the full process

    ## DEMs with days bands added
    if not os.path.exists(os.path.join(f'dems_with_days_{year_1}')):
        os.mkdir(f'dems_with_days_{year_1}')

    if not os.path.exists(os.path.join(f'dems_with_days_{year_2}')):
        os.mkdir(f'dems_with_days_{year_2}')


    ## Merged & Padded DEMs
    if not os.path.exists(os.path.join('merged_dems')):
        os.mkdir('merged_dems')

    if not os.path.exists(os.path.join('merged_dems_padded')):
        os.mkdir('merged_dems_padded')


    ## Difference Maps
    if not os.path.exists(os.path.join('dhdt')):
        os.mkdir('dhdt')

    if not os.path.exists(os.path.join('dhdt_filtered')):
        os.mkdir('dhdt_filtered')

    if not os.path.exists(os.path.join('dhdt_mean_value_rasters')):
        os.mkdir('dhdt_mean_value_rasters')

    if not os.path.exists(os.path.join('dhdt_filtered_and_filled')):
        os.mkdir('dhdt_filtered_and_filled')


    ## Zonal Statistics
    if not os.path.exists(os.path.join('dhdt_zonal_stats')):
        os.mkdir('dhdt_zonal_stats')

    if not os.path.exists(os.path.join('stable_ground_stats')):
        os.mkdir('stable_ground_stats')

    if not os.path.exists(os.path.join('final_stats')):
        os.mkdir('final_stats')


    ## Auxiliary Data
    if not os.path.exists(os.path.join('penetration_rasters')):
        os.mkdir('penetration_rasters')

    if not os.path.exists(os.path.join('years')):
        os.mkdir('years')


def create_days_dictionary(dem_filepaths_1, dem_filepaths_2):

    # Initialize dictionary
    days_since_refday = {}

    # Initialize list for the dem filenames
    dem_files = []

    # Make a single list of all dems
    dem_filepaths = dem_filepaths_1 + dem_filepaths_2

    # Extract filenames from full paths
    dem_files = [os.path.basename(dem) for dem in dem_filepaths]

    # Use regular expressions to extract dates.  Date must be in the format YYYYMMDD
    dates = [re.search('\d{4}(0[1-9]|1[0-2])(0[1-9]|1[0-9]|2[0-9]|3[01])', dem_file).group(0) for dem_file in dem_files]

    # Create lists for days, months and years
    days = [int(d[-2:]) for d in dates]
    months = [int(d[-4:-2]) for d in dates]
    years = [int(d[0:4]) for d in dates] 

    # Create a datetime object representing the reference day, in this case 1st Jan 2000
    refday = datetime.date(2000, 1, 1)

    # Create a list of datetime objects from the day, month and year lists
    dtimes = [datetime.date(year, month, day) for (year, month, day) in zip(years, months, days)]

    # Use timedelta to find number of days between each DEM and the reference day
    days_since_refday_list = [(dtime - refday).days for dtime in dtimes]

    # Build the dictionary that will be used to create the date raster
    for i, d in enumerate(dtimes):
        days_since_refday[d.strftime('%Y%m%d')] = days_since_refday_list[i]


    return days_since_refday


def write_dem_with_days_band(dem, days_since_refday, output_dir):
    '''
    Takes in a DEM file (full or relative path), and a dictionary containing
    datestrings as keys and days between that date and a reference date as values.

    Writes a new raster with two bands, to the output directory specified.  

    Band 1 - the original raster DEM
    Band 2 - each pixel contains the number of days between reference date and the data for that pixel

    '''
    
    # Use regex to extract datestamp in form YYYYMMDD 
    datestamp = re.search('\d{4}(0[1-9]|1[0-2])(0[1-9]|1[0-9]|2[0-9]|3[01])', dem).group(0) 
    
    # Open original DEM
    with rio.open(dem) as src:
        data = src.read(1, masked=True)
        profile = src.profile.copy()
        days_raster = (data / data) * days_since_refday[str(datestamp)]
        profile['dtype'] = np.float32
        profile['nodata'] = -9999
        profile['count'] = 2
    
    with rio.open(os.path.join(output_dir, f'dem_with_days_{datestamp}.tif'), 'w', **profile) as dst:
        dst.write(data.astype(np.float32).filled(-9999), 1)
        dst.write(days_raster.astype(np.float32).filled(-9999), 2)
    
    return True




def merge_dems_with_days_band(dems_with_days_band, year, output_dir):
    '''
    Takes in a list of DEMs which have a days band, merges them,
    and writes a new merged raster with two bands:

    Band 1 - Elevation data
    Band 2 - Number of days between data and reference date

    '''

    open_dems_list = []

    for dem in dems_with_days_band:
        # Create a list of rasterio dataset reader objects i.e. opened files
        open_dems_list.append(rio.open(dem))

        # Create a profile from the first file in the list
        profile = open_dems_list[0].profile.copy()

        # Merge the files
        merge_data, merge_transform = merge(open_dems_list, indexes=[1, 2], method='last')

        # Update profile with new transform and size of the merged raster
        profile['transform'] = merge_transform
        profile['width'] = merge_data.shape[2]
        profile['height'] = merge_data.shape[1]

        # Tell the profile that there will be two bands in the output
        profile['count'] = 2

        # Write both bands of the newly merged raster
        with rio.open(os.path.join(output_dir, f'merged_{year}_with_days_band.tif'), 'w', **profile) as dst:
            dst.write(merge_data.astype(profile['dtype']))
            dst.set_band_description(1, 'elevation')
            dst.set_band_description(2, 'days since 01-01-2000')





def pad_rasters(rasters, output_dir):

    #  Only tested on projected co-ordinate system, not geographic
    #  Need to be careful of using integers if units are degrees lat/lon
    '''
    eg:

    gdalwarp -te -74000 -5589900 83350 -5481300 -tr 10 10 -dstnodata -9999 -co COMPRESS=DEFLATE -co TILED=YES -co
    BLOCKXSIZE=256 -co BLOCKYSIZE=256 merged_2011_with_descriptions.tif merged_2011_with_descriptions_padded.tif
    '''

    max_left = None
    max_bottom = None
    max_right = None
    max_top = None

    # Loop through rasters, finding the minimum bounding rectangle that encompasses all of them
    i=0
    for raster in rasters:
        with rio.open(raster) as src:
            if i==0:
                max_left = int(src.bounds[0])
                max_bottom = int(src.bounds[1])
                max_right = int(src.bounds[2])
                max_top = int(src.bounds[3])
                i+=1
            else:
                if src.bounds[0] < max_left:
                    max_left = int(src.bounds[0])
                if src.bounds[1] < max_bottom:
                    max_bottom = int(src.bounds[1])
                if src.bounds[2] > max_right:
                    max_right = int(src.bounds[2])
                if src.bounds[3] > max_top:
                    max_top = int(src.bounds[3])

    print(f'Bounding box to contain all rasters : Left:{max_left}, Bottom:{max_bottom}, Right:{max_right}, Top:{max_top}\n')

    # Run gdalwarp on each input raster, setting the appropriate bounding box and filling any new areas with nodata value

    for raster in rasters:
        raster_filename = os.path.basename(raster)
        output_filepath = os.path.join(output_dir, f'padded_{raster_filename}')

        if os.path.exists(output_filepath):
            print(f'Destination file {output_filepath} already exists, continuing...')
            continue

        # If padded file does not already exist, run gdalwarp, using the newly found
        # bounds above that will be the same for all of the output files.

        p = subprocess.call(f'gdalwarp -te {max_left} {max_bottom} {max_right} {max_top} -tr 10 10 -dstnodata -9999 -co COMPRESS=DEFLATE -co TILED=YES -co BLOCKXSIZE=256 -co BLOCKYSIZE=256 {raster} {output_filepath}', shell=True)
        
    return True



def write_cropped_with_shapefile(raster_file, ice_shapefile, output_file):
    '''
    Take in a raster and a shapefile

    Crop the raster with the shapefile and
    output the cropped raster

    '''

    with rio.open(raster_file) as src:
        gdf = gpd.read_file(ice_shapefile)
        gdf_reprojected = gdf.to_crs(src.crs)

        shapes = [s for s in gdf_reprojected['geometry']]
        
        out_image, out_transform = rio.mask.mask(src, shapes, crop=True)
        out_meta = src.meta
        out_meta.update({"driver": "GTiff",
                    "height": out_image.shape[1],
                    "width": out_image.shape[2],
                    "transform": out_transform})


    with rio.open(output_file, 'w', **out_meta) as dst:
        dst.write(out_image)

    return True



def filter_raster(raster, output_dir, method='sigma', order=3, lower=2, upper=98):
    with rio.open(raster) as src:

        output = None
        filtered_data = None

        # Read data as masked array
        data = src.read(1, masked=True)
        
        # Copy profile
        profile = src.profile.copy()

        # Filter outliers beyond a number of standard deviations
        if method == 'sigma':

            print(f'Filtering using {order}-sigma method')

            standard_deviation = data.std()
            sigma = standard_deviation * order

            filter_mask = ((data > (data.mean() + sigma)) | (data < (data.mean() - sigma)))
            filtered_data = np.ma.MaskedArray(data, mask=filter_mask)

            output = f'{os.path.basename(raster)[:-4]}_{order}sigma_filtered.tif'


        # Percentile method
        elif method == 'percentile':

                print(f'Filtering using {lower}-{upper} percentile method')

                upper_percentile = np.percentile(data, upper)
                lower_percentile = np.percentile(data, lower)

                filter_mask = (data > upper_percentile) | (data < lower_percentile)
                filtered_data = np.ma.MaskedArray(data, mask=filter_mask)

                output = f'{os.path.basename(raster)[:-4]}_{lower}-{upper}_percentile_filtered.tif'

    if (output is not None) and (filtered_data is not None):
        # Write resulting filtered raster
        with rio.open(os.path.join(output_dir, output), 'w', **profile) as dst:
            dst.write(filtered_data.filled(-9999), 1)

    else:
        return


def filter_raster_by_specific_values(raster, output_dir, lower=-40, upper=40):

    # Set output filename for filtered data
    output = f'{os.path.basename(raster)[:-4]}_{lower}-{upper}_filtered.tif'

    with rio.open(raster) as src:

        # Read data as masked array
        data = src.read(1, masked=True)
        
        # Copy profile
        profile = src.profile.copy()

        filter_mask = (data > upper) | (data < lower)
        filtered_data = np.ma.MaskedArray(data, mask=filter_mask)

    with rio.open(os.path.join(output_dir, output), 'w', **profile) as dst:
        dst.write(filtered_data.filled(-9999), 1)

    return True


def stats_mask2(geoms, datas, **mask_kw):
    mask_kw.setdefault('crop', True)
    mask_kw.setdefault('all_touched', False)
    mask_kw.setdefault('filled', False)
    result, mask_transform = robmask(dataset=datas, shapes=(geoms,), **mask_kw)
    return result


def create_mean_value_raster(template_raster, stats_shapefile, output_dir, column_name='mean'):

    output_filename = f'{os.path.basename(template_raster)[:-4]}_mean_per_elevation_bin.tif'

    # The template raster is only used to get a profile with
    # which to write the new raster
    with rio.open(template_raster) as src:
        profile = src.profile

    # Read the shapefile with the zonal stats in it
    gdf = gpd.read_file(stats_shapefile)

    # Reproject to same CRS as template raster
    gdf = gdf.to_crs(profile['crs'])

    # create a generator of geom, value pairs to use in rasterizing
    shapes = ((geom,value) for geom, value in zip(gdf.geometry, gdf[column_name]))

    # Write new output
    with rio.open(os.path.join(output_dir, output_filename), 'w+', **profile) as dst:
        data = dst.read(1)
        rasterized = rio.features.rasterize(shapes=shapes, out=data, fill=-9999, transform=dst.transform)
        dst.write(rasterized, 1)

    return True


def radar_penetration_raster(refdem, shapefile, ela, ela_pen, upper_limit, upper_limit_pen, output_dir):

    # Output filename
    penetration_raster_filename = os.path.join(output_dir, f'penetration-ela-{ela}m-ela-pen-{0}m-upper-limit-{upper_limit}-upper-limit-pen-{upper_limit_pen}m.tif')

    with rio.open(refdem) as src: 
        gdf = gpd.read_file(shapefile)
        gdf_reprojected = gdf.to_crs(src.crs)

        shapes = [s for s in gdf_reprojected['geometry']]
        
        data, out_transform = rio.mask.mask(src, shapes, crop=True)
        out_meta = src.meta
        out_meta.update({"driver": "GTiff",
                    "height": data.shape[1],
                    "width": data.shape[2],
                    "transform": out_transform})

    # keep only DEM data between the ELA and the upper limit set above
    mask = (data > upper_limit) | (data < ela) 
    dem = np.ma.MaskedArray(data, mask=mask, fill_value=out_meta['nodata'])

    # Subtract the ELA altitude from DEM values to remove areas where there is no penetration
    dem_salient = dem - ela     

    # Calculate penetration per pixel based on linear increase from 
    # "ela_pen" metres at ELA to "upper_limit_pen" at upper elevation limit
    penetration = (dem_salient / (upper_limit - ela)) * (upper_limit_pen - ela_pen) 

    # Write resulting raster back out
    with rio.open(penetration_raster_filename, 'w', **out_meta) as dst:
        dst.write(penetration.filled(-9999))
    
    return True


def fill_dhdt_with_mean_values(dhdt_file, means_file, output_dir):

    output_filename = f'{os.path.basename(dhdt_file)[:-4]}_filled.tif'

    with rio.open(dhdt_file) as src:
        dhdt = src.read(1, masked=True)
        profile = src.profile.copy()

    with rio.open(means_file) as src:
        means = src.read(1, masked=True)

    # Get the mask from the dhdt raster, which tells us where there are gaps
    holes_in_dhdt = np.ma.getmask(dhdt)

    # Make a copy of the dhdt data, which will be filled
    filled = dhdt.copy()

    # Set values in filled which are masked to the value of the means array at the same locations
    filled[holes_in_dhdt] = means[holes_in_dhdt]

    with rio.open(os.path.join(output_dir, output_filename), 'w', **profile) as dst:
        dst.write(filled, 1)

    return True


def zonal_stats_by_basin(raster, shapefile, output_dir, output_crs=3762):

    output_projection = CRS.from_epsg(output_crs)
    
    output_filename = f'stats_{os.path.basename(raster[:-4])}.geojson'

    with rio.open(raster) as src:

        # Read shapefile
        gdf_unprojected = gpd.read_file(shapefile)

        # Reproject to same CRS as input raster
        gdf = gdf_unprojected.to_crs(src.crs)

        # Clean empty and null geometries
        gdf = gdf.loc[~(gdf.geometry.is_empty | gdf.geometry.isna())]

        # Explode MultiPolygons to turn them into Polygons
        gdf = gdf.loc[gdf.geometry.geom_type == 'MultiPolygon'].explode().reset_index(drop=True)

        # Create new 'mean' attribute and apply stats_mask2 to populate this column
        #gdf['mean'] = gdf.geometry.apply(lambda x: stats_mask2(x, src)).apply(np.ma.mean)
        gdf = gpd.GeoDataFrame.from_features(zonal_stats(gdf, raster, geojson_out=True, stats='mean', nodata=-9999), crs=gdf.crs)

        # Reproject GeoDataFrame to output_crs
        gdf.to_crs(epsg=output_crs, inplace=True)

        # Write the new GeoJSON file to output_dir
        print(os.path.join(output_dir, output_filename))
        gdf.to_file(os.path.join(output_dir, output_filename), driver='GeoJSON')
    
    return True



def calculate_stable_ground_area_weighted_std(dem, slope_binned_shapefile, output_dir, output_crs):

    # Load slope binned shapefile, bins from 0 to 50.  0 = from 0 to 10 degrees etc
    slope_gdf = gpd.read_file(slope_binned_shapefile)

    # Add area column for info and to calculate area weighting of standard deviation
    slope_gdf['area_m2'] = slope_gdf.geometry.area

    with rio.open(dem) as src:
        data = src.read(1, masked=True)
        profile = src.profile.copy()

        # Filter input raster with 2-98% quantile filter
        mask = (data < np.percentile(data, 2)) | (data > np.percentile(data, 98))
        filtered = np.ma.MaskedArray(data.copy(), mask=mask)

        # Reproject slope shapefile to the same CRS as the input raster 
        slope_gdf = slope_gdf.to_crs(profile['crs'])

        # Sort by feature type eg Polygon, MultiPolygon in order to handle each differently
        features = slope_gdf.groupby(slope_gdf['geometry'].geom_type)

        with MemoryFile() as memfile:
            with memfile.open(**profile) as mem:
                mem.write(filtered, 1)

                # Process different feature types separately.  Currently only set up for MultiPolygons
                multipolygons = features.get_group('MultiPolygon')

                # Calculate and add stats            
                multipolygons['mean'] = multipolygons.geometry.apply(lambda x: stats_mask2(x, mem)).apply(np.ma.mean)
                multipolygons['std'] = multipolygons.geometry.apply(lambda x: stats_mask2(x, mem)).apply(np.ma.std)

                try:
                    multipolygons.to_crs(epsg=output_crs, inplace=True)
                    multipolygons.to_file(os.path.join(output_dir, f'stable_ground_stats_{os.path.basename(dem)[:-4]}.geojson'), driver='GeoJSON')
                except TypeError as e:
                    print(e)

    # Calculate the actual area weighted standard deviation for the whole stable ground area now
    total_area = multipolygons['area_m2'].sum()
    multipolygons['area_x_std'] = multipolygons['area_m2'] * multipolygons['std']
    area_weighted_std = multipolygons['area_x_std'].sum() / total_area

    return area_weighted_std


def calculate_error(dMdt, dhdt, p2a, p2a_paulref, A, dA, dp, p, Vpen, dt, sd_AW_dhdt, lag_distance):
    '''
    Calculates the overall mass balance error according to error formula in dissertation
    '''
    Sf = A
    A_km2 = A / 1_000_000   # convert m2 to km2 

    # Calculate Sc, which is the surface area of autocorrelation, from the lag_distance
    Sc = math.pi * (lag_distance ** 2)
    
    # Calculate the corrected error in dh/dt from Sc, plus glacier area Sf and the 
    # Area-Weighted standard deviations of elevation change on stable ground 
    error_dhdt = math.sqrt(Sc/(5 * Sf)) * sd_AW_dhdt
    
    # Calculate the mass balance error
    error_dMdt = math.sqrt((dMdt**2 * ((error_dhdt/dhdt)**2 + ((dA/A_km2)**2) + (dp/p)**2) + ((Vpen / dt) * p)))

    return error_dMdt