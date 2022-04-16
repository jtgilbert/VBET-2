# imports
import geopandas as gpd
import rasterio
import rasterio.mask
from rasterio.features import shapes
from shapely.geometry import Point, LineString, Polygon, MultiPolygon
from shapely.ops import unary_union, cascaded_union
from rasterstats import zonal_stats
import numpy as np
import skimage.morphology as mo
from scipy.signal import convolve2d
from scipy.linalg import lstsq
import json
import os.path
from tqdm import tqdm


class VBET:
    """
    The Valley Bottom Extraction Tool (V-BET) extracts a valley bottom of floodplain from a DEM using a
    stream network.
    """
    def __init__(self, **kwargs):

        self.network = gpd.read_file(kwargs['network'])
        self.streams = kwargs['network']
        self.dem = kwargs['dem']
        self.out = kwargs['out']
        self.scratch = kwargs['scratch']  # make so if dir doesnt exist, script creates it
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
        self.da_field = kwargs['da_field']
        self.lg_depth = kwargs['lg_depth']
        self.med_depth = kwargs['med_depth']
        self.sm_depth = kwargs['sm_depth']

        # either use selected drainage area field, or pull drainage area from raster
        if self.da_field is not None:
            self.network['Drain_Area'] = self.network[self.da_field]

        # set crs for output
        self.crs_out = self.network.crs

        # check that scratch directory exists, make if not
        if os.path.exists(self.scratch):
            pass
        else:
            os.mkdir(self.scratch)

        # check that datasets are in projected coordinate system
        if self.network.crs is None:
            raise Exception('All geospatial inputs should be have the same projected coordinate reference system')
        if rasterio.open(self.dem).crs is None:
            raise Exception('All geospatial inputs should be have the same projected coordinate reference system')
        if self.dr_area is not None:
            if rasterio.open(self.dr_area).crs is None:
                raise Exception('All geospatial inputs should be have the same projected coordinate reference system')

        # check that all depth params have a value or are None
        if self.lg_depth is None and (self.med_depth is not None or self.sm_depth is not None):
            raise Exception('All depths (lg, med, sm) must have value or be None')
        elif self.med_depth is None and (self.lg_depth is not None or self.sm_depth is not None):
            raise Exception('All depths (lg, med, sm) must have value or be None')
        elif self.sm_depth is None and (self.lg_depth is not None or self.med_depth is not None):
            raise Exception('All depths (lg, med, sm) must have value or be None')

        # check that there are no segments with less than 5 vertices
        few_verts = []
        for i in self.network.index:
            if len(self.network.loc[i].geometry.xy[0]) <= 5:
                few_verts.append(i)
        if len(few_verts) > 0:
            raise Exception("Network segments with IDs ", few_verts, "don't have enough vertices for DEM detrending. "
                                                                     "Add vertices in GIS")

        # add container for individual valley bottom features and add the minimum buffer into it
        self.polygons = []

        network_geom = self.network['geometry']
        min_buf = network_geom.buffer(self.min_buf)

        for x in range(len(min_buf)):
            self.polygons.append(min_buf[x])

        # save total network length for use in later parameter
        self.seglengths = 0
        for x in self.network.index:
            self.seglengths += self.network.loc[x].geometry.length

    def clean_network(self):

        print('Cleaning up drainage network for VBET input')
        print('starting with {} network segments'.format(len(self.network)))
        # minimum length - remove short segments
        with rasterio.open(self.dem, 'r') as src:
            xres = src.res[0]
        self.network = self.network[self.network.geometry.length > 5*xres]

        # get rid of perfectly straight segments
        sin = []
        for i in self.network.index:
            seg_geom = self.network.loc[i, 'geometry']
            x_coord1 = seg_geom.boundary[0].xy[0][0]
            y_coord1 = seg_geom.boundary[0].xy[1][0]
            x_coord2 = seg_geom.boundary[1].xy[0][0]
            y_coord2 = seg_geom.boundary[1].xy[1][0]

            dist = ((x_coord1 - x_coord2) ** 2 + (y_coord1 - y_coord2) ** 2) ** 0.5
            sin_val = seg_geom.length / dist
            sin.append(sin_val)
        self.network['sinuos'] = sin

        self.network = self.network[self.network['sinuos'] >= 1.00001]

        print('cleaned to {} network segments'.format(len(self.network)))

    def add_da(self):
        """
        Adds a drainage area attribute to each segment of the drainage network
        :return:
        """
        print('Adding drainage area to network')
        da_list = []

        for i in self.network.index:
            seg = self.network.loc[i]
            geom = seg['geometry']
            pos = int(len(geom.coords.xy[0])/2)
            mid_pt_x = geom.coords.xy[0][pos]
            mid_pt_y = geom.coords.xy[1][pos]

            pt = Point(mid_pt_x, mid_pt_y)
            buf = pt.buffer(50)  # make buffer distance function of resolution (e.g. 5*res)

            zs = zonal_stats(buf, self.dr_area, stats='max')
            da_val = zs[0].get('max')

            da_list.append(da_val)

        self.network['Drain_Area'] = da_list

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

    def detrend(self, dem, seg_geom):
        with rasterio.open(dem) as src:
            meta = src.profile
            arr = src.read()[0, :, :]
            res_x = src.res[0]
            res_y = src.res[1]
            x_min = src.transform[2]
            y_max = src.transform[5]
            y_min = y_max - (src.height*res_y)

        # points along network in real coords
        _xs = seg_geom.xy[0][::2]
        _ys = seg_geom.xy[1][::2]

        zs = np.zeros_like(_xs)

        for i in range(len(_xs)):
            pt = Point(_xs[i], _ys[i])
            buf = pt.buffer(5)
            zonal = zonal_stats(buf, dem, stats='min')
            val = zonal[0].get('min')

            zs[i] = val

        # points in array coords
        xs = np.zeros_like(_xs)
        ys = np.zeros_like(_ys)

        for i in range(len(_xs)):
            xs[i] = int((_xs[i] - x_min) / res_x)  # column in array space
            ys[i] = int((y_max - _ys[i]) / res_y)  # row in array space

        xs = xs[np.isfinite(zs)]
        ys = ys[np.isfinite(zs)]
        zs = zs[np.isfinite(zs)]  # its currently possible to use only 2 points..?

        # do fit
        tmp_A = []
        tmp_b = []
        for i in range(len(xs)):
            tmp_A.append([xs[i], ys[i], 1])
            tmp_b.append(zs[i])
        b = np.array(tmp_b).T
        A = np.array(tmp_A)
        fit = lstsq(A, b)

        trend = np.full((src.height, src.width), src.nodata, dtype=src.dtypes[0])
        for j in range(trend.shape[0]):
            for i in range(trend.shape[1]):
                trend[j, i] = fit[0][0] * i + fit[0][1] * j + fit[0][2]

        out_arr = arr - trend

        return out_arr

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
                elif thresh >= array[j, i] > 0:
                    out_array[j, i] = 1
                else:
                    array[j, i] = ndval

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
            transform = src.transform
            crs = src.crs

        results = (
            {'properties': {'raster_val': v}, 'geometry': s}
            for i, (s, v)
            in enumerate(
                shapes(array, mask=array == 1., transform=transform)))

        geoms = list(results)
        if len(geoms) == 0:
            return 0
            #raise Exception('No raster cells to convert to shapefile, try changing model parameters')

        else:
            df = gpd.GeoDataFrame.from_features(geoms)
            df.crs = crs
            geom = df['geometry']

            area = []

            for x in range(len(geom)):
                area.append(geom[x].area)
                self.polygons.append(geom[x])

            return sum(area)

    def getFeatures(self, gdf):
        """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""

        return [json.loads(gdf.to_json())['features'][0]['geometry']]

    def chaikins_corner_cutting(self, coords, refinements=5):
        coords = np.array(coords)

        for _ in range(refinements):
            L = coords.repeat(2, axis=0)
            R = np.empty_like(L)
            R[0] = L[0]
            R[2::2] = L[1:-1:2]
            R[1:-1:2] = L[2::2]
            R[-1] = L[-1]
            coords = L * 0.75 + R * 0.25

        return coords

    def valley_bottom(self):
        """
        Run the VBET algorithm
        :return: saves a valley bottom shapefile
        """

        self.clean_network()

        print('Generating valley bottom for each network segment')
        for i in tqdm(self.network.index):
            seg = self.network.loc[i]
            da = seg['Drain_Area']
            geom = seg['geometry']

            # print('segment ', i+1, ' of ', len(self.network.index))

            ept1 = (geom.boundary[0].xy[0][0], geom.boundary[0].xy[1][0])
            ept2 = (geom.boundary[1].xy[0][0], geom.boundary[1].xy[1][0])
            line = LineString([ept1, ept2])

            if da >= self.lg_da:
                buf = line.buffer(self.lg_buf, cap_style=1)
            elif self.lg_da > da >= self.med_da:
                buf = line.buffer(self.med_buf, cap_style=1)
            else:
                buf = line.buffer(self.sm_buf, cap_style=1)

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
            elif self.lg_da > da >= self.med_da:
                slope_sub = self.reclassify(slope, ndval, self.med_slope)
            else:
                slope_sub = self.reclassify(slope, ndval, self.sm_slope)

            # set thresholds for hole filling
            avlen = int(self.seglengths / len(self.network))
            if da < self.med_da:
                thresh = avlen * self.sm_buf * 0.005
            elif self.med_da <= da < self.lg_da:
                thresh = avlen * self.med_buf * 0.005
            else:  # da >= self.lg_da:
                thresh = avlen * self.lg_buf * 0.005

            # only detrend if depth is being used
            if self.med_depth is not None:
                detr = self.detrend(dem, geom)  # might want to change this offset

                if da >= self.lg_da:
                    depth = self.reclassify(detr, ndval, self.lg_depth)
                elif self.lg_da > da >= self.med_da:
                    depth = self.reclassify(detr, ndval, self.med_depth)
                else:
                    depth = self.reclassify(detr, ndval, self.sm_depth)

                overlap = self.raster_overlap(slope_sub, depth, ndval)
                if 1 in overlap:
                    filled = self.fill_raster_holes(overlap, thresh, ndval)
                    a = self.raster_to_shp(filled, dem)
                    self.network.loc[i, 'fp_area'] = a
                else:
                    self.network.loc[i, 'fp_area'] = 0
            else:
                if 1 in slope_sub:
                    filled = self.fill_raster_holes(slope_sub, thresh, ndval)
                    a = self.raster_to_shp(filled, dem)
                    self.network.loc[i, 'fp_area'] = a
                else:
                    self.network.loc[i, 'fp_area'] = 0

            demsrc.close()

        self.network.to_file(self.streams)

        # merge all polygons in folder and dissolve
        print("Merging valley bottom segments")
        vb = gpd.GeoSeries(unary_union(self.polygons))  #
        vb.crs = self.crs_out
        vb.to_file(self.scratch + "/tempvb.shp")

        # simplify and smooth polygon
        print("Cleaning valley bottom")
        vbc = gpd.read_file(self.scratch + "/tempvb.shp")
        vbc = vbc.simplify(3, preserve_topology=True)  # make number a function of dem resolution
        vbc.to_file(self.scratch + "/tempvb.shp")

        # get rid of small unattached polygons
        self.network.to_file(self.scratch + "/dissnetwork.shp")
        network2 = gpd.read_file(self.scratch + "/dissnetwork.shp")
        network2['dissolve'] = 1
        network2 = network2.dissolve('dissolve')
        vb1 = gpd.read_file(self.scratch + "/tempvb.shp")
        vbm2s = vb1.explode()
        sub = []
        for i in vbm2s.index:
            segs = 0
            for j in network2.index:
                if network2.loc[j].geometry.intersects(vbm2s.loc[i].geometry):
                    segs += 1
            if segs > 0:
                sub.append(True)
            else:
                sub.append(False)


        #sub = vbm2s['geometry'].area >= 0.001*vbm2s['geometry'].area.max()
        vbcut = vbm2s[sub].reset_index(drop=True)
        vbcut.to_file(self.scratch + "/tempvb.shp")
        # do it a second time?
        #vb2 = gpd.read_file(self.scratch + "/tempvb.shp")
        #vbm2s2 = vb2.explode()
        #sub2 = vbm2s2['geometry'].area >= 0.05*vbm2s2['geometry'].area.max()
        #vbcut2 = vbm2s2[sub2].reset_index(drop=True)
        #vbcut2.to_file(self.scratch + "/tempvb.shp")

        # TRY THIS - have to dissolve before chaikins

        polys = []
        for i in vbcut.index:
            coords = list(vbcut.loc[i].geometry.exterior.coords)  # vbcut WAS vbc when using shapely simplify.
            new_coords = self.chaikins_corner_cutting(coords)
            polys.append(Polygon(new_coords))

        if len(polys) > 1:
            p = MultiPolygon(polys)
        else:
            p = polys[0]

        vbf = gpd.GeoDataFrame(index=[0], crs=self.crs_out, geometry=[p])
        vbf = vbf.explode()
        areas = []
        for i in vbf.index:
            areas.append(vbf.loc[i].geometry.area/1000000.)
        vbf['Area_km2'] = areas

        vbf.to_file(self.out)

        return

