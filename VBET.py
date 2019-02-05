# imports
import geopandas as gpd
import rasterio
import rasterio.mask
from rasterio.features import shapes
from shapely.geometry import Point
from shapely.ops import cascaded_union
from rasterstats import zonal_stats
import numpy as np
import skimage.morphology as mo
from scipy.signal import convolve2d
import json


class VBET:
    """
    The Valley Bottom Extraction Tool (V-BET) extracts a valley bottom of floodplain from a DEM using a
    stream network.
    """
    def __init__(self, **kwargs):

        self.network = gpd.read_file(kwargs['network'])
        self.dem = kwargs['dem']
        self.out = kwargs['out']
        self.scratch = kwargs['scratch']
        self.lg_da = kwargs['lg_da']
        self.med_da = kwargs['med_da']
        self.lg_slope = kwargs['lg_slope']
        self.med_slope = kwargs['med_slope']
        self.sm_slope = kwargs['sm_slope']
        self.lg_buf = kwargs['lg_buf']
        self.med_buf = kwargs['med_buf']
        self.sm_buf = kwargs['sm_buf']
        self.min_buf = kwargs['min_buf']
        self.dr_area = kwargs['dr_area']
        self.lg_depth = kwargs['lg_depth']
        self.med_depth = kwargs['med_depth']

        self.crs_out = self.network.crs

        self.polygons = []

    def add_da(self):
        """
        Adds a drainage area attribute to each segment of the drainage network
        :return:
        """
        da_list = []

        for i in self.network.index:
            seg = self.network.loc[i]
            geom = seg['geometry']
            pos = int(len(geom.coords.xy[0])/2)
            mid_pt_x = geom.coords.xy[0][pos]
            mid_pt_y = geom.coords.xy[1][pos]

            pt = Point(mid_pt_x, mid_pt_y)
            buf = pt.buffer(50)

            zs = zonal_stats(buf, self.dem, stats='max')
            da_val = zs[0].get('max')

            da_list.append(da_val)

        self.network['Drain_Area'] = da_list

        # check for segments with lower DA value than upstream segment
        # maybe add this..? would have to add network topology

        return

    def add_elev(self):
        """
        Adds a median stream elevation value to each segment of the drainage network
        :return:
        """
        elev_list = []

        for i in self.network.index:
            segment = self.network.loc[i]

            seg_geom = segment['geometry']
            pos = int(len(seg_geom.coords.xy[0]) / 2)
            mid_pt_x = seg_geom.coords.xy[0][pos]
            mid_pt_y = seg_geom.coords.xy[1][pos]

            pt = Point(mid_pt_x, mid_pt_y)
            buf = pt.buffer(20)

            zs = zonal_stats(buf, self.dem, stats='min')
            elev_val = zs[0].get('min')

            elev_list.append(elev_val)

        self.network['Elev'] = elev_list

        return

    def slope(self, dem):
        """
        Finds the slope using partial derivative method
        :param dem: path to a digital elevation raster
        :return: a 2-D array with the values representing slope for the cell
        """
        with rasterio.open(dem, 'r') as src:
            arr = src.read()[0, :, :]

            xres = src.res[0]
            yres = src.res[1]

        x = np.array([[-1 / (8 * xres), 0, 1 / (8 * xres)],
                      [-2 / (8 * xres), 0, 2 / (8 * xres)],
                      [-1 / (8 * xres), 0, 1 / (8 * xres)]])
        y = np.array([[1 / (8 * yres), 2 / (8 * yres), 1 / (8 * yres)],
                      [0, 0, 0],
                      [-1 / (8 * yres), -2 / (8 * yres), -1 / (8 * yres)]])

        x_grad = convolve2d(arr, x, mode='same', boundary='fill', fillvalue=1)
        y_grad = convolve2d(arr, y, mode='same', boundary='fill', fillvalue=1)
        slope = np.arctan(np.sqrt(x_grad ** 2 + y_grad ** 2)) * (180. / np.pi)
        slope = slope.astype(src.dtypes[0])

        return slope

    def reclassify(self, array, ndval, thresh):
        """
        Splits an input array into two values: 1 and NODATA based on a threshold value
        :param array: a 2-D array
        :param ndval: NoData value
        :param thresh: The threshold value. Values < thresh are converted to 1
        and values > thresh are converted to NoData
        :return: a 2-D array of with values of 1 and NoData
        """
        rows, cols = array.shape

        out_array = np.full(array.shape, ndval)

        for j in range(0, rows - 1):
            for i in range(0, cols - 1):
                if array[j, i] == ndval:
                    out_array[j, i] = ndval
                elif array[j, i] > thresh:
                    out_array[j, i] = ndval
                else:
                    out_array[j, i] = 1

        return out_array

    def raster_overlap(self, array1, array2, ndval):
        """
        Finds the overlap between two orthogonal arrays (same dimensions)
        :param array1: first 2-D array
        :param array2: second 2-D array
        :param ndval: a no data value
        :return: 2-D array with a value of 1 where both input arrays have values and value of NoData where either of
        input arrays have NoData
        """
        if array1.shape != array2.shape:
            raise Exception('rasters are not same size')

        out_array = np.full(array1.shape, ndval)

        for j in range(0, array1.shape[0] - 1):
            for i in range(0, array1.shape[1] - 1):
                if array1[j, i] == 1. and array2[j, i] == 1.:
                    out_array[j, i] = 1.
                else:
                    out_array[j, i] = ndval

        return out_array

    def fill_raster_holes(self, array, thresh, ndval):
        """
        Fills in holes and gaps in an array of 1s and NoData
        :param array: 2-D array of 1s and NoData
        :param thresh: hole size (cells) below which should be filled
        :param ndval: NoData value
        :return: 2-D array like input array but with holes filled
        """
        binary = np.zeros_like(array, dtype=bool)
        for j in range(0, array.shape[0] - 1):
            for i in range(0, array.shape[1] - 1):
                if array[j, i] == 1:
                    binary[j, i] = 1

        b = mo.remove_small_holes(binary, thresh, 1)
        c = mo.binary_closing(b, selem=np.ones((7, 7)))
        d = mo.remove_small_holes(c, thresh, 1)

        out_array = np.full(d.shape, ndval, dtype=np.float32)
        for j in range(0, d.shape[0] - 1):
            for i in range(0, d.shape[1] - 1):
                if d[j, i] == True:
                    out_array[j, i] = 1.

        return out_array

    def array_to_raster(self, array, raster_like, raster_out):
        """
        Save an array as a raster dataset
        :param array: array to convert to raster
        :param raster_like: a raster from which to take metadata (e.g. spatial reference, nodata value etc.)
        :param raster_out: path to store output raster
        :return:
        """
        with rasterio.open(raster_like, 'r') as src:
            meta = src.profile
            dtype = src.dtypes[0]

        out_array = np.asarray(array, dtype)

        with rasterio.open(raster_out, 'w', **meta) as dst:
            dst.write(out_array, 1)

        return

    def raster_to_shp(self, array, raster_like):
        """
        Convert the 1 values in an array of 1s and NoData to a polygon
        :param array: 2-D array of 1s and NoData
        :param raster_like: a raster from which to take metadata (e.g. spatial reference)
        :param shp_out: path to store output shapefile
        :return:
        """
        with rasterio.open(raster_like) as src:
            transform = src.affine
            crs = src.crs

        results = (
            {'properties': {'raster_val': v}, 'geometry': s}
            for i, (s, v)
            in enumerate(
                shapes(array, mask=array == 1., transform=transform)))

        geoms = list(results)

        df = gpd.GeoDataFrame.from_features(geoms)
        df.crs = crs
        #df.to_file(shp_out)
        geom = df['geometry']

        for x in range(len(geom)):
            self.polygons.append(geom[x])

        return

    def getFeatures(self, gdf):
        """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""

        return [json.loads(gdf.to_json())['features'][0]['geometry']]

    def valley_bottom(self):
        """
        Run the VBET algorithm
        :return: saves a valley bottom shapefile
        """

        for i in self.network.index:
            seg = self.network.loc[i]
            elev = seg['Elev']
            da = seg['Drain_Area']
            geom = seg['geometry']

            print 'elev', elev

            if da >= self.lg_da:
                buf = geom.buffer(self.lg_buf)
            elif da < self.lg_da and da >= self.med_da:
                buf = geom.buffer(self.med_buf)
            else:
                buf = geom.buffer(self.sm_buf)

            bufds = gpd.GeoSeries(buf)
            coords = self.getFeatures(bufds)

            with rasterio.open(self.dem) as src:
                out_image, out_transform = rasterio.mask.mask(src, coords, crop=True)
                out_meta = src.meta.copy()

            out_meta.update({'driver': 'Gtiff',
                             'height': out_image.shape[1],
                             'width': out_image.shape[2],
                             'transform': out_transform})
            with rasterio.open(self.scratch + '/dem_sub.tif', 'w', **out_meta) as dest:
                dest.write(out_image)

            dem = self.scratch + "/dem_sub.tif"
            demsrc = rasterio.open(dem)
            demarray = demsrc.read()[0, :, :]
            ndval = demsrc.nodata

            slope = self.slope(dem)

            if da >= self.lg_da:
                slope_sub = self.reclassify(slope, ndval, self.lg_slope)
            elif da < self.lg_da and da >= self.med_da:
                slope_sub = self.reclassify(slope, ndval, self.med_slope)
            else:
                slope_sub = self.reclassify(slope, ndval, self.sm_slope)

            if da >= self.lg_da:
                depth = self.reclassify(demarray, ndval, elev + self.lg_depth)
            elif da < self.lg_da and da >= self.med_da:
                depth = self.reclassify(demarray, ndval, elev + self.med_depth)
            else:
                depth = None

            if depth is not None:
                overlap = self.raster_overlap(slope_sub, depth, ndval)
                filled = self.fill_raster_holes(overlap, 10000, ndval)
                self.raster_to_shp(filled, dem)
                #poly = gpd.read_file(self.scratch + '/poly' + str(i))
                #poly_geom = poly.loc[0, 'geometry']
                #polygons.append(geom)
            else:
                filled = self.fill_raster_holes(slope_sub, 10000, ndval)
                print filled[0:10, 0:10]
                self.raster_to_shp(filled, dem)
                #poly = gpd.read_file(self.scratch + '/poly' + str(i))
                #poly_geom = poly.loc[0, 'geometry']
                #polygons.append(geom)

        # merge all polygons in folder and dissolve
        vb = gpd.GeoSeries(cascaded_union(self.polygons))
        vb.crs = self.crs_out
        vb.to_file(self.out)

