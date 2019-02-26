import VBET


params = {
    'network': '/home/jordan/Documents/lost_horse/NHD_Lost_Horse.shp',
    'dem': '/home/jordan/Documents/lost_horse/output_be.tif',
    'out': '/home/jordan/Documents/lost_horse/vb_lh.shp',
    'scratch': '/home/jordan/Documents/lost_horse/scratch',
    'lg_da': 250,
    'med_da': 25,
    'lg_slope': 2,
    'med_slope': 4,
    'sm_slope': 6,
    'lg_buf': 250,
    'med_buf': 150,
    'sm_buf': 100,
    'min_buf': 8,
    'dr_area': '/home/jordan/Documents/lost_horse/topo/DrArea_Lost_Horse.tif',
    'lg_depth': 2,
    'med_depth': 1.5,
    'sm_depth': 1.
}

vb = VBET.VBET(**params)
vb.add_da()
vb.valley_bottom()

