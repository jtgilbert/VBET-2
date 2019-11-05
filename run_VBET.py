import VBET


params = {
    'network': '/path/to/stream/network.shp',
    'dem': '/path/to/dem.tif',
    'out': '/path/to/store/valleybottom/output.shp',
    'scratch': '/path/to/some/scratch/workspace/folder',
    'lg_da': 250,
    'med_da': 25,
    'lg_slope': 2,
    'med_slope': 4,
    'sm_slope': 6,
    'lg_buf': 500,
    'med_buf': 150,
    'sm_buf': 30,
    'min_buf': 8,
    'dr_area': '/path/to/drainage/area/raster.tif',
    'lg_depth': 2,
    'med_depth': 1.5,
    'sm_depth': 1.25
}

vb = VBET.VBET(**params)
vb.add_da()
vb.valley_bottom()

