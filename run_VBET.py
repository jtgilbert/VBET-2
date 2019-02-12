import VBET


params = {
    'network': '/home/jordan/Documents/RioGrande/Lower_RG.shp',
    'dem': '/home/jordan/Documents/RioGrande/4m_clipped.tif',
    'out': '/home/jordan/Documents/RioGrande/vb_lower_rg3.shp',
    'scratch': '/home/jordan/Documents/RioGrande/scratch2',
    'lg_da': 250,
    'med_da': 25,
    'lg_slope': 1.5,
    'med_slope': 3,
    'sm_slope': 10,
    'lg_buf': 1000,
    'med_buf': 150,
    'sm_buf': 10,
    'min_buf': 8,
    'dr_area': '/home/jordan/Documents/RioGrande/Rio_Grande_DrArea.tif',
    'lg_depth': 1.7,
    'med_depth': 1
}

vb = VBET.VBET(**params)
vb.add_da()
vb.valley_bottom()

# try running in current form with longer segments, ~2 km maybe