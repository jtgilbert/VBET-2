import VBET


params = {
    'network': '/home/jordan/Documents/Geoscience/LH/NHD_LH_sub.shp',
    'dem': '/home/jordan/Documents/Geoscience/LH/output_be.tif',
    'out': '/home/jordan/Documents/Geoscience/LH/vb_lh.shp',
    'scratch': '/home/jordan/Documents/Geoscience/LH/scratch',
    'lg_da': 250,
    'med_da': 10,
    'lg_slope': 2,
    'med_slope': 4,
    'sm_slope': 6,
    'lg_buf': 250,
    'med_buf': 100,
    'sm_buf': 50,
    'min_buf': 8,
    'dr_area': '/home/jordan/Documents/Geoscience/LH/LH_DrArea.tif',
    'lg_depth': 2,
    'med_depth': 2,
    'sm_depth': 1.75
}

vb = VBET.VBET(**params)
vb.add_da()
vb.valley_bottom()

