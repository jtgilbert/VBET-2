import VBET
from datetime import datetime


params = {
    'network': '/path/to/stream/network.shp',
    'dem': '/path/to/DEM.tif',
    'out': '/path/to/valleybottom/output.shp',
    'scratch': '/path/for/scratch',
    'lg_da': 300,
    'med_da': 30,
    'lg_slope': 2,
    'med_slope': 3,
    'sm_slope': 4,
    'lg_buf': 1500,
    'med_buf': 400,
    'sm_buf': 100,
    'min_buf': 12,
    'dr_area': '/path/to/drainage/area/raster.tif',
    'lg_depth': 4,
    'med_depth': 2,
    'sm_depth': 1
}

print('started: ', datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
vb = VBET.VBET(**params)
vb.add_da()
vb.valley_bottom()
print('ended: ', datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
